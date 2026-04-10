# -----------------------------
# 1) Chuẩn bị dữ liệu & hàm phụ
# -----------------------------
library(readxl); library(dplyr); library(tidyr); library(zoo)
library(xts); library(rugarch); library(forecast)
library(lmtest); library(sandwich); library(FinTS)
cpu <-read_excel("data/cpu.xlsx")
#cpu <- readxl::read_excel("C:/Users/trong/OneDrive/HUYNHTT/Year 2025/Thang 11/cpu.xlsx")
cpu <- cpu %>%
  mutate(ym = as.integer(ym),
         date = as.Date(as.yearmon(as.character(ym), "%Y%m")),
         RGOLD_RCPU = RGOLD * RCPU) %>%
  arrange(date)

dep_vars <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")

# Mapping FX theo thị trường
fx_map <- c(
  RSET  = "RTHB", # Thailand
  RVNI  = "RVND", # Vietnam (HOSE)
  RHNX  = "RVND", # Vietnam (HNX)
  RSTI  = "RSGD", # Singapore
  RKLSE = "RMYR", # Malaysia
  RJKSE = "RIDR", # Indonesia
  RPSE  = "RPHP"  # Philippines
)

other_exog <- c("RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")  # dùng chung

MAE  <- function(e) mean(abs(e), na.rm=TRUE)
RMSE <- function(e) sqrt(mean(e^2, na.rm=TRUE))

# -----------------------------
# 2) Hàm chạy từng thị trường (OLS + GARCH-X + OOS)
# -----------------------------
run_country_model <- function(yvar, data, roll_window=48) {
  message("== Running for ", yvar, " ==")
  
  fx_var    <- unname(fx_map[[yvar]])
  exog_vars <- c(fx_var, other_exog)
  
  df <- data %>%
    dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>%
    tidyr::drop_na()
  
  # ---- OLS in-sample ----
  f_ols <- as.formula(paste(yvar, "~",
                            paste(c(fx_var,"RWTI","RVIX","RGOLD","RCPU"), collapse=" + "),
                            "+ RGOLD:RCPU"))
  reg_ols <- lm(f_ols, data=df)
  bw <- sandwich::bwNeweyWest(reg_ols)
  se <- sandwich::NeweyWest(reg_ols, lag=bw, prewhite=TRUE, adjust=TRUE)
  ols_sum <- coeftest(reg_ols, vcov.=se)
  
  arch <- tryCatch(FinTS::ArchTest(residuals(reg_ols), lags=12)$p.value, error=function(e) NA_real_)
  bp   <- tryCatch(lmtest::bptest(reg_ols)$p.value,                     error=function(e) NA_real_)
  
  # ---- GARCH-X in-sample ----
  y_xts <- xts(df[[yvar]], order.by=df$date)
  X_mat <- as.matrix(df[, exog_vars])
  spec <- ugarchspec(
    mean.model     = list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X_mat),
    variance.model = list(model="sGARCH", garchOrder=c(1,1)),
    distribution.model="std"
  )
  fit <- tryCatch(
    ugarchfit(spec, data=y_xts, solver="hybrid", solver.control=list(trace=0)),
    error=function(e) NULL
  )
  
  # ---- Rolling OOS (48M) ----
  N <- nrow(df)
  oos <- vector("list", max(0, N-roll_window))
  k <- 0
  for (i in seq(roll_window, N-1)) {
    tr <- (i-roll_window+1):i
    te <- i+1
    y_tr <- df[[yvar]][tr]; y_te <- df[[yvar]][te]
    X_tr <- as.matrix(df[tr, exog_vars]); X_te <- as.matrix(df[te, exog_vars])
    
    # GARCH-X forecast
    sp <- ugarchspec(
      mean.model=list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X_tr),
      variance.model=list(model="sGARCH", garchOrder=c(1,1)),
      distribution.model="std"
    )
    fit_r <- tryCatch(ugarchfit(sp, data=y_tr, solver="hybrid",
                                solver.control=list(trace=0)), error=function(e) NULL)
    pred_gx <- if(!is.null(fit_r)) {
      fc <- tryCatch(ugarchforecast(fit_r, n.ahead=1,
                                    external.forecasts=list(mreg=matrix(X_te, nrow=1))),
                     error=function(e) NULL)
      if (is.null(fc)) NA_real_ else as.numeric(fitted(fc))
    } else NA_real_
    
    # RW và ARIMA
    pred_rw <- tail(y_tr,1)
    fit_a <- suppressWarnings(auto.arima(y_tr, ic="aic", seasonal=FALSE))
    pred_a <- as.numeric(forecast(fit_a,h=1)$mean)
    
    k <- k + 1
    oos[[k]] <- data.frame(date=df$date[te],
                           y_true=y_te,
                           yhat_garchx=pred_gx,
                           yhat_rw=pred_rw,
                           yhat_arima=pred_a)
  }
  oos_df <- bind_rows(oos) %>%
    mutate(e_gx=y_true-yhat_garchx, e_rw=y_true-yhat_rw, e_ar=y_true-yhat_arima)
  oos_mae  <- c(GARCHX=MAE(oos_df$e_gx), RW=MAE(oos_df$e_rw), ARIMA=MAE(oos_df$e_ar))
  oos_rmse <- c(GARCHX=RMSE(oos_df$e_gx), RW=RMSE(oos_df$e_rw), ARIMA=RMSE(oos_df$e_ar))
  
  list(
    yvar=yvar,
    FX=fx_var,
    OLS=ols_sum,
    pARCH=arch,
    pBP=bp,
    GARCHX_summary=if(!is.null(fit)) coef(fit) else NA,
    OOS_MAE=oos_mae,
    OOS_RMSE=oos_rmse
  )
}

results <- lapply(dep_vars, run_country_model, data=cpu, roll_window=48)

# Tóm tắt nhanh bảng MAE/RMSE
mae_tbl <- do.call(rbind, lapply(results, function(x)
  data.frame(Market=x$yvar, t(x$OOS_MAE))))
rmse_tbl <- do.call(rbind, lapply(results, function(x)
  data.frame(Market=x$yvar, t(x$OOS_RMSE))))

cat("\n== MAE across markets ==\n"); print(mae_tbl, row.names=FALSE)
cat("\n== RMSE across markets ==\n"); print(rmse_tbl, row.names=FALSE)


# ==================================================
# Robustness regression for ASEAN-7 stock indices
# (OLS HAC – bản rút gọn theo yêu cầu, dùng FX mapping)
# ==================================================


# 2) Markets & regressors
dep_vars  <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")
fx_map <- c(RSET="RTHB", RVNI="RVND", RHNX="RVND", RSTI="RSGD", RKLSE="RMYR", RJKSE="RIDR", RPSE="RPHP")
other_exog <- c("RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")

# 3) Storage
ols_results <- list()

# 4) Helper: tidy coeftest (position-based)
tidy_coeftest <- function(ct_mat) {
  stopifnot(is.matrix(ct_mat), ncol(ct_mat) >= 4)
  data.frame(
    term      = rownames(ct_mat),
    Estimate  = as.numeric(ct_mat[, 1]),
    Std.Error = as.numeric(ct_mat[, 2]),
    t.value   = as.numeric(ct_mat[, 3]),
    p.value   = as.numeric(ct_mat[, 4]),
    row.names = NULL, check.names = FALSE
  )
}

# 5) Loop across markets
for (yvar in dep_vars) {
  message("\n===========================\n== Running OLS for ", yvar, "\n===========================")
  
  fx_var    <- unname(fx_map[[yvar]])
  exog_vars <- c(fx_var, other_exog)
  
  df <- cpu %>% dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>% tidyr::drop_na()
  
  form <- as.formula(paste(
    yvar, "~",
    paste(c(fx_var, "RWTI","RVIX","RGOLD","RCPU"), collapse = " + "),
    "+ RGOLD:RCPU"
  ))
  reg  <- lm(form, data = df)
  
  # HAC vcov
  bw  <- sandwich::bwNeweyWest(reg)
  hac <- sandwich::NeweyWest(reg, lag = bw, prewhite = TRUE, adjust = TRUE)
  
  # Coefs (HAC)
  ct  <- lmtest::coeftest(reg, vcov. = hac)
  tbl <- tidy_coeftest(ct)
  
  # Diagnostics
  bp_p   <- tryCatch(lmtest::bptest(reg)$p.value, error = function(e) NA_real_)
  arch_p <- tryCatch(FinTS::ArchTest(residuals(reg), lags = 12)$p.value, error = function(e) NA_real_)
  
  # Print per-market
  cat("\n--- OLS (HAC) results for", yvar, "---\n")
  print(tbl, row.names = FALSE, digits = 6)
  cat("HAC bandwidth (bwNeweyWest):", bw, "\n")
  cat("Breusch–Pagan p-value:", formatC(bp_p, digits = 4, format = "fg"),
      "| ARCH-LM(12) p-value:", formatC(arch_p, digits = 4, format = "fg"), "\n")
  
  # Store
  ols_results[[yvar]] <- list(
    yvar = yvar, fx = fx_var, table = tbl, bw = bw, bp_p = bp_p, arch_p = arch_p
  )
}

# 6) Combined summary table
combine_one <- function(res) {
  tbl <- res$table
  pull_val <- function(term, col) {
    v <- tbl[[col]][match(term, tbl$term)]
    ifelse(is.na(v), NA_real_, as.numeric(v))
  }
  data.frame(
    Market = res$yvar,
    FX_var = res$fx,
    Intercept   = pull_val("(Intercept)", "Estimate"),
    Intercept_p = pull_val("(Intercept)", "p.value"),
    FX      = pull_val(res$fx, "Estimate"),
    FX_p    = pull_val(res$fx, "p.value"),
    RWTI    = pull_val("RWTI", "Estimate"),
    RWTI_p  = pull_val("RWTI", "p.value"),
    RVIX    = pull_val("RVIX", "Estimate"),
    RVIX_p  = pull_val("RVIX", "p.value"),
    RGOLD   = pull_val("RGOLD", "Estimate"),
    RGOLD_p = pull_val("RGOLD", "p.value"),
    RCPU    = pull_val("RCPU", "Estimate"),
    RCPU_p  = pull_val("RCPU", "p.value"),
    RGOLD_RCPU   = pull_val("RGOLD:RCPU", "Estimate"),
    RGOLD_RCPU_p = pull_val("RGOLD:RCPU", "p.value"),
    BP_p = res$bp_p, ARCHLM12_p = res$arch_p, HAC_bw = res$bw,
    check.names = FALSE
  )
}

summary_ols <- do.call(rbind, lapply(ols_results, combine_one)) %>% arrange(Market)

cat("\n==============================\n== OLS (HAC) – Combined summary\n==============================\n")
print(summary_ols, row.names = FALSE, digits = 6)

# 7) Pretty (stars) for paper
star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < .001, "***",
                ifelse(p < .01,  "**",
                       ifelse(p < .05,  "*",
                              ifelse(p < .10,  ".", "")))))
}

pretty_tbl <- summary_ols %>%
  mutate(
    `(Intercept)` = sprintf("%.3f%s", Intercept,   star(Intercept_p)),
    FX            = sprintf("%.3f%s", FX,          star(FX_p)),
    RWTI          = sprintf("%.3f%s", RWTI,        star(RWTI_p)),
    RVIX          = sprintf("%.3f%s", RVIX,        star(RVIX_p)),
    RGOLD         = sprintf("%.3f%s", RGOLD,       star(RGOLD_p)),
    RCPU          = sprintf("%.3f%s", RCPU,        star(RCPU_p)),
    `RGOLD:RCPU`  = sprintf("%.3f%s", RGOLD_RCPU,  star(RGOLD_RCPU_p)),
    `BP p`        = sprintf("%.4f", BP_p),
    `ARCH-LM(12) p` = sprintf("%.4f", ARCHLM12_p)
  ) %>%
  dplyr::select(Market, FX_var, `(Intercept)`, FX, RWTI, RVIX, RGOLD, RCPU, `RGOLD:RCPU`, `BP p`, `ARCH-LM(12) p`, HAC_bw)

cat("\n==============================\n== OLS (HAC) – Paper-ready table (coef with stars)\n==============================\n")
print(pretty_tbl, row.names = FALSE, right = FALSE)

# ==================================================
# GARCH-X for ASEAN-7 stock indices — Full script
# Mean: AR(1) với external regressors = (FX theo thị trường, RWTI, RVIX, RGOLD, RCPU, RGOLD×RCPU)
# Var : sGARCH(1,1), Student-t
# ==================================================

# Helpers
pstar <- function(p) ifelse(is.na(p), "",
                            ifelse(p<.001,"***", ifelse(p<.01,"**", ifelse(p<.05,"*", ifelse(p<.10,".","")))))
lb_pval <- function(x, lag=10) {
  out <- tryCatch(Box.test(x, lag=lag, type="Ljung-Box")$p.value, error=function(e) NA_real_); as.numeric(out)
}
tidy_garch <- function(fit) {
  est <- coef(fit)
  vc  <- tryCatch(vcov(fit, robust = FALSE), error=function(e) NULL)
  vcr <- tryCatch(vcov(fit, robust = TRUE ), error=function(e) NULL)
  
  se  <- if (!is.null(vc))  sqrt(diag(vc))  else rep(NA_real_, length(est))
  ser <- if (!is.null(vcr)) sqrt(diag(vcr)) else rep(NA_real_, length(est))
  
  z    <- est / se
  zr   <- est / ser
  p    <- 2 * pnorm(abs(z),  lower.tail = FALSE)
  pr   <- 2 * pnorm(abs(zr), lower.tail = FALSE)
  
  df <- data.frame(
    param    = names(est),
    Estimate = as.numeric(est),
    SE_std   = as.numeric(se),
    z_std    = as.numeric(z),
    p_std    = as.numeric(p),
    SE_rob   = as.numeric(ser),
    z_rob    = as.numeric(zr),
    p_rob    = as.numeric(pr),
    row.names = NULL,
    check.names = FALSE
  )
  
  a1  <- df$Estimate[df$param=="alpha1"]; b1 <- df$Estimate[df$param=="beta1"]
  apb <- if (length(a1) && length(b1)) as.numeric(a1 + b1) else NA_real_
  
  eps  <- residuals(fit, standardize = TRUE)
  lb10_p  <- lb_pval(eps, 10)
  lb10_p2 <- lb_pval(eps^2, 10)
  
  list(df = df,
       alpha1 = a1, beta1 = b1, alpha_plus_beta = apb,
       LB10_p = lb10_p, LB10_sq_p = lb10_p2)
}

# 3) Loop & estimate
garch_results <- list()
cat("\n================ GARCH-X estimation (AR(1)-sGARCH(1,1)-t) ================\n")

for (yvar in dep_vars) {
  message("\n== Running GARCH-X for ", yvar, " ==")
  
  fx_var    <- unname(fx_map[[yvar]])
  exog_vars <- c(fx_var, other_exog)
  
  df <- cpu %>% dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>% drop_na()
  y  <- xts(df[[yvar]], order.by = df$date)
  X  <- as.matrix(df[, exog_vars])
  
  spec <- ugarchspec(
    mean.model      = list(armaOrder = c(1,0), include.mean = TRUE,
                           external.regressors = X),
    variance.model  = list(model = "sGARCH", garchOrder = c(1,1)),
    distribution.model = "std"
  )
  
  fit <- tryCatch(
    ugarchfit(spec = spec, data = y, solver = "hybrid", solver.control = list(trace = 0)),
    error = function(e) { message("!! GARCH failed for ", yvar, " — ", e$message); NULL }
  )
  
  if (is.null(fit)) next
  
  tg <- tidy_garch(fit)
  
  garch_results[[yvar]] <- list(
    yvar = yvar,
    tidy = tg$df,
    alpha1 = tg$alpha1, beta1 = tg$beta1, alpha_plus_beta = tg$alpha_plus_beta,
    LB10_p = tg$LB10_p, LB10_sq_p = tg$LB10_sq_p
  )
  
  row_mx6 <- subset(tg$df, param == "mxreg6")
  if (nrow(row_mx6)) {
    cat(sprintf("RGOLD×RCPU (mxreg6): Estimate=%.5f | p_std=%.4g | p_rob=%.4g\n",
                row_mx6$Estimate, row_mx6$p_std, row_mx6$p_rob))
  } else cat("RGOLD×RCPU (mxreg6): <not found>\n")
  
  cat(sprintf("alpha=%.4f, beta=%.4f, alpha+beta=%.4f | LB-Q(10) p=%.3f; LB-Q^2(10) p=%.3f\n",
              as.numeric(tg$alpha1), as.numeric(tg$beta1), as.numeric(tg$alpha_plus_beta),
              as.numeric(tg$LB10_p), as.numeric(tg$LB10_sq_p)))
}

# 4) Combined summary table
pull_param <- function(res, pname, which = c("Estimate","p_std","p_rob")) {
  which <- match.arg(which)
  if (is.null(res) || is.null(res$tidy)) return(NA_real_)
  v <- res$tidy[[which]][match(pname, res$tidy$param)]
  ifelse(is.na(v), NA_real_, as.numeric(v))
}

summary_tbl <- do.call(rbind, lapply(dep_vars, function(mkt){
  res <- garch_results[[mkt]]
  data.frame(
    Market = mkt,
    mu      = pull_param(res, "mu",      "Estimate"),
    mu_p    = pull_param(res, "mu",      "p_rob"),
    ar1     = pull_param(res, "ar1",     "Estimate"),
    ar1_p   = pull_param(res, "ar1",     "p_rob"),
    
    # Mean regressors: mxreg1..mxreg6 = (FX, RWTI, RVIX, RGOLD, RCPU, RGOLD×RCPU)
    FX             = pull_param(res, "mxreg1",  "Estimate"),
    FX_p_rob       = pull_param(res, "mxreg1",  "p_rob"),
    RWTI           = pull_param(res, "mxreg2",  "Estimate"),
    RWTI_p_rob     = pull_param(res, "mxreg2",  "p_rob"),
    RVIX           = pull_param(res, "mxreg3",  "Estimate"),
    RVIX_p_rob     = pull_param(res, "mxreg3",  "p_rob"),
    RGOLD          = pull_param(res, "mxreg4",  "Estimate"),
    RGOLD_p_rob    = pull_param(res, "mxreg4",  "p_rob"),
    RCPU           = pull_param(res, "mxreg5",  "Estimate"),
    RCPU_p_rob     = pull_param(res, "mxreg5",  "p_rob"),
    RGOLD_RCPU     = pull_param(res, "mxreg6",  "Estimate"),
    RGOLD_RCPU_p_rob = pull_param(res, "mxreg6","p_rob"),
    
    alpha1        = if (!is.null(res)) as.numeric(res$alpha1) else NA_real_,
    beta1         = if (!is.null(res)) as.numeric(res$beta1)  else NA_real_,
    alpha_plus_beta = if (!is.null(res)) as.numeric(res$alpha_plus_beta) else NA_real_,
    LB_Q10_p      = if (!is.null(res)) as.numeric(res$LB10_p) else NA_real_,
    LB_Q10sq_p    = if (!is.null(res)) as.numeric(res$LB10_sq_p) else NA_real_,
    check.names = FALSE
  )
})) %>% arrange(Market)

cat("\n==============================\n== GARCH-X – Combined summary (robust p-values)\n==============================\n")
print(summary_tbl, row.names = FALSE, digits = 6)

# 5) Paper-ready
pretty_tbl <- summary_tbl %>%
  transmute(
    Market,
    mu        = sprintf("%.3f%s", mu,      pstar(mu_p)),
    ar1       = sprintf("%.3f%s", ar1,     pstar(ar1_p)),
    FX        = sprintf("%.3f%s", FX,      pstar(FX_p_rob)),
    RWTI      = sprintf("%.3f%s", RWTI,    pstar(RWTI_p_rob)),
    RVIX      = sprintf("%.3f%s", RVIX,    pstar(RVIX_p_rob)),
    RGOLD     = sprintf("%.3f%s", RGOLD,   pstar(RGOLD_p_rob)),
    RCPU      = sprintf("%.3f%s", RCPU,    pstar(RCPU_p_rob)),
    `RGOLD×RCPU` = sprintf("%.3f%s", RGOLD_RCPU, pstar(RGOLD_RCPU_p_rob)),
    `α`       = sprintf("%.3f", alpha1),
    `β`       = sprintf("%.3f", beta1),
    `α+β`     = sprintf("%.3f", alpha_plus_beta),
    `LB Q(10) p`   = sprintf("%.3f", LB_Q10_p),
    `LB Q²(10) p`  = sprintf("%.3f", LB_Q10sq_p)
  )

cat("\n==============================\n== GARCH-X – Paper-ready (mean coefs with stars; robust p)\n==============================\n")
print(pretty_tbl, row.names = FALSE, right = FALSE)

# 6) Save CSVs (bật nếu cần)
# write.csv(summary_ols, "ASEAN7_OLS_HAC_summary.csv", row.names = FALSE)
# write.csv(pretty_tbl,  "ASEAN7_OLS_HAC_pretty.csv",  row.names = FALSE)
# write.csv(summary_tbl, "ASEAN7_GARCHX_summary.csv",  row.names = FALSE)
# write.csv(pretty_tbl,  "ASEAN7_GARCHX_pretty.csv",   row.names = FALSE)


# =========================================================
# Granger causality cho TỪNG CẶP biến (X → Y)
# - Phương pháp: Toda–Yamamoto (mặc định) hoặc Classic
# - Yêu cầu cài: readxl, dplyr, tidyr, zoo, vars, urca, car, lmtest
# Dữ liệu: file cpu.xlsx (các cột đã có trong dự án của bạn)
# =========================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(zoo)
  library(vars); library(urca); library(car); library(lmtest)
})

# ---------- 0) Đọc & tiền xử lý ----------
cpu <- readxl::read_excel("C:/Users/trong/OneDrive/HUYNHTT/Year 2025/Thang 11/cpu.xlsx") %>%
  mutate(
    ym   = as.integer(ym),
    date = as.Date(as.yearmon(as.character(ym), "%Y%m"))
  ) %>%
  arrange(date)

# Map thị trường -> cặp (Stock, FX)
stock_fx_map <- list(
  RSET  = "RTHB", # Thailand
  RVNI  = "RVND", # Vietnam (HOSE)
  RHNX  = "RVND", # Vietnam (HNX)
  RSTI  = "RSGD", # Singapore
  RKLSE = "RMYR", # Malaysia
  RJKSE = "RIDR", # Indonesia
  RPSE  = "RPHP"  # Philippines
)

# Các biến vĩ mô chung
cpu_var <- "RCPU"
wti_var <- "RWTI"

# ---------- 1) Helper: chọn độ trễ theo chỉ số IC ----------
.pick_ic_lag <- function(sel, ic = c("AIC","BIC","HQ")) {
  ic <- match.arg(ic)
  # vars::VARselect trả list "selection" (vector tên) và các ma trận AIC(n), HQ(n), SC(n), FPE
  # An toàn: nếu thiếu, fallback về min theo cột tương ứng
  if (!is.null(sel$selection) && ic %in% names(sel$selection)) {
    return(as.integer(sel$selection[[ic]]))
  }
  # Fallback: tìm ma trận tương ứng
  mat_name <- switch(ic, AIC="AIC", BIC="SC", HQ="HQ")
  M <- sel[[mat_name]]
  if (is.null(M)) {
    # fallback cuối: lấy min AIC nếu có, nếu không thì lag=2
    if (!is.null(sel$AIC)) return(as.integer(which.min(sel$AIC)))
    return(2L)
  }
  as.integer(which.min(M))
}

# ---------- 2) Helper: xác định bậc tích phân bằng ADF (robust tên cval) ----------
adf_Iorder <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 30) return(NA_integer_)
  ur0 <- tryCatch(urca::ur.df(x, type = "drift", lags = 4, selectlags = "AIC"),
                  error = function(e) NULL)
  if (is.null(ur0)) return(NA_integer_)
  
  tst <- tryCatch(ur0@teststat, error = function(e) NULL)
  if (is.null(tst)) return(NA_integer_)
  tau_name <- names(tst)[grep("^tau", names(tst))]
  if (length(tau_name) == 0) return(NA_integer_)
  stat_tau <- as.numeric(tst[[tau_name[1]]])
  
  cval <- tryCatch(ur0@cval, error = function(e) NULL)
  if (is.null(cval) || !is.matrix(cval)) return(NA_integer_)
  
  col_5 <- which(colnames(cval) %in% c("5pct","5%","5 percent","5 pc"))
  if (length(col_5) == 0) col_5 <- ceiling(ncol(cval)/2)
  
  row_tau <- grep("^tau", rownames(cval))
  if (length(row_tau) == 0) row_tau <- 1
  
  crit_5 <- as.numeric(cval[row_tau[1], col_5[1]])
  if (!is.finite(stat_tau) || !is.finite(crit_5)) return(NA_integer_)
  
  # Bác bỏ H0 (đơn vị gốc) → I(0); ngược lại I(1)
  if (stat_tau < crit_5) 0L else 1L
}

# ---------- 3) Hàm kiểm định Granger cho TỪNG CẶP (X → Y) ----------
# df   : data.frame có cột date, Y, X
# y,x  : tên biến ký tự ("RSET", "RTHB", "RCPU", "RWTI", ...)
# max_lag: trễ tối đa cho VARselect
# ic     : "AIC"|"BIC"|"HQ" để chọn k
# method : "TY" (Toda–Yamamoto) hoặc "classic"
granger_pair <- function(df, y, x, max_lag = 12, ic = c("AIC","BIC","HQ"),
                         method = c("TY","classic")) {
  ic <- match.arg(ic); method <- match.arg(method)
  
  dsub <- df %>%
    dplyr::select(date, Y = dplyr::all_of(y), X = dplyr::all_of(x)) %>%
    tidyr::drop_na()
  
  # Nếu trùng tên hay không tồn tại, báo lỗi rõ ràng
  if (!all(c("Y","X") %in% names(dsub))) {
    stop("Không tìm thấy biến: ", y, " hoặc ", x)
  }
  if (nrow(dsub) < 40) stop("Chuỗi quá ngắn cho VAR/Granger.")
  
  # Chọn độ trễ k theo IC
  sel <- vars::VARselect(dsub[, c("Y","X")], lag.max = max_lag, type = "const")
  k <- .pick_ic_lag(sel, ic)
  k <- max(1L, min(k, max_lag))
  
  if (method == "classic") {
    # VAR(k) ở mức; dùng vars::causality
    fit <- vars::VAR(dsub[, c("Y","X")], p = k, type = "const")
    cxy <- vars::causality(fit, cause = "X")$Granger
    # Trích p-value và thống kê
    return(list(
      method = "Classic",
      lag_k  = k,
      stat   = unname(cxy$statistic),
      pvalue = unname(cxy$p.value),
      parameter = unname(cxy$parameter),
      decision = ifelse(cxy$p.value < 0.05, "Reject H0: X ↛ Y (Có nhân quả)", "Không bác bỏ H0")
    ))
  }
  
  # ----- Toda–Yamamoto -----
  # dmax = max bậc tích phân của (Y,X)
  dY <- adf_Iorder(dsub$Y); dX <- adf_Iorder(dsub$X)
  dmax <- max(c(0L, dY, dX), na.rm = TRUE)
  p_aug <- k + dmax
  
  fit <- vars::VAR(dsub[, c("Y","X")], p = p_aug, type = "const")
  
  # Wald test: trong MÔ HÌNH CỦA Y, ràng buộc các hệ số X.l1..X.lk = 0
  # Hệ số được đặt tên như "Y.l1","X.l1", ...
  eqY <- fit$varresult$Y
  # Xây hệ ràng buộc cho 1..k (không đụng đến phần dmax)
  hyp <- paste0("X.l", 1:k, " = 0")
  wald <- car::linearHypothesis(eqY, hyp, test = "Chisq") # theo TY ⇒ Chi-square
  
  pval <- as.numeric(wald[2, "Pr(>Chisq)"])
  stat <- as.numeric(wald[2, "Chisq"])
  df   <- as.numeric(wald[2, "Df"])
  
  list(
    method   = "Toda–Yamamoto",
    lag_k    = k,
    dmax     = dmax,
    p_aug    = p_aug,
    stat     = stat,
    df       = df,
    pvalue   = pval,
    decision = ifelse(is.finite(pval) && pval < 0.05, "Reject H0: X ↛ Y (Có nhân quả)", "Không bác bỏ H0"),
    note     = "Test trên phương trình Y; ràng buộc X.l1..X.lk = 0"
  )
}

# ---------- 4) Ví dụ gọi HAI-BIẾN cho TỪNG THỊ TRƯỜNG ----------
# Bạn chỉ việc chọn thị trường & cặp (X → Y) rồi gọi granger_pair(cpu, Y, X, ...)
# Ví dụ minh hoạ cho THÁI LAN (RSET ~ RTHB, RCPU, RWTI):

# Chọn thị trường
market <- "RSET"
fx_var <- stock_fx_map[[market]]  # "RTHB"

# (1) FX → Stock
gc_fx_to_stock_TY  <- granger_pair(cpu, y = market, x = fx_var,  max_lag = 12, ic = "AIC", method = "TY")
gc_fx_to_stock_cl  <- granger_pair(cpu, y = market, x = fx_var,  max_lag = 12, ic = "AIC", method = "classic")

# (2) Stock → FX
gc_stock_to_fx_TY  <- granger_pair(cpu, y = fx_var, x = market,  max_lag = 12, ic = "AIC", method = "TY")

# (3) CPU → Stock
gc_cpu_to_stock_TY <- granger_pair(cpu, y = market, x = cpu_var, max_lag = 12, ic = "AIC", method = "TY")

# (4) WTI → Stock
gc_wti_to_stock_TY <- granger_pair(cpu, y = market, x = wti_var, max_lag = 12, ic = "AIC", method = "TY")

# (5) CPU → FX
gc_cpu_to_fx_TY    <- granger_pair(cpu, y = fx_var, x = cpu_var, max_lag = 12, ic = "AIC", method = "TY")

# (6) WTI → FX
gc_wti_to_fx_TY    <- granger_pair(cpu, y = fx_var, x = wti_var, max_lag = 12, ic = "AIC", method = "TY")

# In thử một kết quả
print(gc_fx_to_stock_TY)
print(gc_cpu_to_stock_TY)

# ---------- 5) Gợi ý chạy nhanh cho 7 thị trường ----------
# Ví dụ chạy FX → Stock (TY) cho tất cả thị trường:
gc_fx2stock_all <- lapply(names(stock_fx_map), function(mkt) {
  fx <- stock_fx_map[[mkt]]
  out <- granger_pair(cpu, y = mkt, x = fx, max_lag = 12, ic = "AIC", method = "TY")
  data.frame(Market = mkt, Pair = paste0(fx," → ", mkt),
             lag_k = out$lag_k, dmax = ifelse(is.null(out$dmax), NA, out$dmax),
             stat = out$stat, df = out$df, pvalue = out$pvalue,
             decision = out$decision)
})
gc_fx2stock_all <- dplyr::bind_rows(gc_fx2stock_all)
print(gc_fx2stock_all)





# ==================================================
# ASEAN-7: OLS (HAC) + GARCH-X (AR(1)-sGARCH(1,1)-t)
# Biến dùng: RGEPU + RGOLD×RGEPU
# ==================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(zoo)
  library(lmtest); library(sandwich); library(FinTS)
  library(xts); library(rugarch)
})

# --- 1) Load & chuẩn bị dữ liệu ---
cpu <- read_excel("C:/Users/trong/OneDrive/HUYNHTT/Year 2025/Thang 11/cpu.xlsx") %>%
  mutate(
    date = as.Date(as.yearmon(as.character(ym <- as.integer(ym)), "%Y%m")),
    RGOLD_RGEPU = RGOLD * RGEPU
  ) %>%
  arrange(date)

dep_vars <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")
fx_map <- c(RSET="RTHB", RVNI="RVND", RHNX="RVND", RSTI="RSGD",
            RKLSE="RMYR", RJKSE="RIDR", RPSE="RPHP")

# thứ tự exog: FX, RWTI, RVIX, RGOLD, RGEPU, RGOLD_RGEPU
other_exog <- c("RWTI","RVIX","RGOLD","RGEPU","RGOLD_RGEPU")

star <- function(p) ifelse(is.na(p),"",
                           ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*",
                                                                        ifelse(p<.10,".","")))))

# --- 2) OLS (HAC) ---
cat("\n================ OLS (HAC) ================\n")
ols_out <- list()
for (yvar in dep_vars) {
  fx <- fx_map[[yvar]]
  df <- cpu %>% dplyr::select(dplyr::all_of(yvar), dplyr::all_of(fx), dplyr::all_of(other_exog)) %>% drop_na()
  f  <- as.formula(paste(
    yvar,"~",paste(c(fx,"RWTI","RVIX","RGOLD","RGEPU"), collapse="+"), "+ RGOLD:RGEPU"
  ))
  m   <- lm(f, data=df)
  hac <- NeweyWest(m, lag=bwNeweyWest(m), prewhite=TRUE, adjust=TRUE)
  ct  <- coeftest(m, vcov.=hac)
  bp  <- tryCatch(bptest(m)$p.value, error=function(e) NA)
  arch<- tryCatch(ArchTest(residuals(m), lags=12)$p.value, error=function(e) NA)
  
  cat("\n---", yvar, "---\n"); print(ct, digits=4)
  cat("BP p =", sprintf("%.4f", bp), "| ARCH-LM(12) p =", sprintf("%.4f", arch), "\n")
  
  ols_out[[yvar]] <- data.frame(term=rownames(ct), Estimate=ct[,1],
                                pval=ct[,4], Market=yvar, BP=bp, ARCH=arch)
}

# --- 3) GARCH-X (robust p-value) ---
cat("\n================ GARCH-X (robust p-values) ================\n")
garch_out <- list()
for (yvar in dep_vars) {
  fx <- fx_map[[yvar]]
  df <- cpu %>% dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(fx), dplyr::all_of(other_exog)) %>% drop_na()
  y  <- xts(df[[yvar]], order.by=df$date)
  X  <- as.matrix(df[, c(fx, other_exog)])
  
  spec <- ugarchspec(
    mean.model     = list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X),
    variance.model = list(model="sGARCH", garchOrder=c(1,1)),
    distribution.model="std"
  )
  fit <- tryCatch(ugarchfit(spec, data=y, solver="hybrid", solver.control=list(trace=0)),
                  error=function(e) NULL)
  if (is.null(fit)) next
  
  est  <- coef(fit)
  vcr  <- vcov(fit, robust=TRUE)
  se   <- sqrt(diag(vcr))
  pval <- 2*pnorm(abs(est/se), lower.tail=FALSE)
  alpha1 <- est["alpha1"]; beta1 <- est["beta1"]
  
  cat("\n---", yvar, "---\n")
  print(data.frame(param=names(est), Estimate=round(est,4), pval=round(pval,4)), row.names=FALSE)
  cat("alpha+beta =", round(alpha1+beta1,3), "\n")
  
  garch_out[[yvar]] <- data.frame(param=names(est), Estimate=est, pval=pval, Market=yvar)
}

# --- 4) Xem nhanh: hệ số RGEPU & RGOLD×RGEPU (mxreg5, mxreg6) ---
gpeu_tbl <- lapply(garch_out, function(x){
  subset(x, param %in% c("mxreg5","mxreg6")) %>%
    mutate(Sig = star(pval)) %>%
    transmute(Market, param, Estimate=round(Estimate,3), pval=round(pval,3), Sig)
}) %>% bind_rows()

cat("\n=== GARCH-X: RGEPU và RGOLD×RGEPU (với stars) ===\n")
print(gpeu_tbl)



# ==================================================
# ASEAN-7: OLS (HAC) + GARCH-X (AR(1)-sGARCH(1,1)-t)
# Biến dùng: RCPU + RGOLD×RCPU
# ==================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(zoo)
  library(lmtest); library(sandwich); library(FinTS)
  library(xts); library(rugarch)
})

# --- 1) Load & chuẩn bị dữ liệu ---
cpu <- read_excel("C:/Users/trong/OneDrive/HUYNHTT/Year 2025/Thang 11/cpu.xlsx") %>%
  mutate(
    date = as.Date(as.yearmon(as.character(ym <- as.integer(ym)), "%Y%m")),
    RGOLD_RCPU = RGOLD * RCPU
  ) %>%
  arrange(date)

dep_vars <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")
fx_map <- c(RSET="RTHB", RVNI="RVND", RHNX="RVND", RSTI="RSGD",
            RKLSE="RMYR", RJKSE="RIDR", RPSE="RPHP")

# thứ tự exog: FX, RWTI, RVIX, RGOLD, RCPU, RGOLD_RCPU
other_exog <- c("RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")

star <- function(p) ifelse(is.na(p),"",
                           ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*",
                                                                        ifelse(p<.10,".","")))))

# --- 2) OLS (HAC) ---
cat("\n================ OLS (HAC) ================\n")
ols_out <- list()
for (yvar in dep_vars) {
  fx <- fx_map[[yvar]]
  df <- cpu %>% dplyr::select(dplyr::all_of(yvar), dplyr::all_of(fx), dplyr::all_of(other_exog)) %>% tidyr::drop_na()
  f  <- as.formula(paste(
    yvar,"~",paste(c(fx,"RWTI","RVIX","RGOLD","RCPU"), collapse="+"), "+ RGOLD:RCPU"
  ))
  m   <- lm(f, data=df)
  hac <- NeweyWest(m, lag=bwNeweyWest(m), prewhite=TRUE, adjust=TRUE)
  ct  <- coeftest(m, vcov.=hac)
  bp  <- tryCatch(bptest(m)$p.value, error=function(e) NA)
  arch<- tryCatch(ArchTest(residuals(m), lags=12)$p.value, error=function(e) NA)
  
  cat("\n---", yvar, "---\n"); print(ct, digits=4)
  cat("BP p =", sprintf("%.4f", bp), "| ARCH-LM(12) p =", sprintf("%.4f", arch), "\n")
  
  ols_out[[yvar]] <- data.frame(term=rownames(ct), Estimate=ct[,1],
                                pval=ct[,4], Market=yvar, BP=bp, ARCH=arch)
}

# --- 3) GARCH-X (robust p-value) ---
cat("\n================ GARCH-X (robust p-values) ================\n")
garch_out <- list()
for (yvar in dep_vars) {
  fx <- fx_map[[yvar]]
  df <- cpu %>% dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(fx), dplyr::all_of(other_exog)) %>% tidyr::drop_na()
  y  <- xts(df[[yvar]], order.by=df$date)
  X  <- as.matrix(df[, c(fx, other_exog)])
  
  spec <- ugarchspec(
    mean.model     = list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X),
    variance.model = list(model="sGARCH", garchOrder=c(1,1)),
    distribution.model="std"
  )
  fit <- tryCatch(ugarchfit(spec, data=y, solver="hybrid", solver.control=list(trace=0)),
                  error=function(e) NULL)
  if (is.null(fit)) next
  
  est  <- coef(fit)
  vcr  <- vcov(fit, robust=TRUE)
  se   <- sqrt(diag(vcr))
  pval <- 2*pnorm(abs(est/se), lower.tail=FALSE)
  alpha1 <- est["alpha1"]; beta1 <- est["beta1"]
  
  cat("\n---", yvar, "---\n")
  print(data.frame(param=names(est), Estimate=round(est,4), pval=round(pval,4)), row.names=FALSE)
  cat("alpha+beta =", round(alpha1+beta1,3), "\n")
  
  garch_out[[yvar]] <- data.frame(param=names(est), Estimate=est, pval=pval, Market=yvar)
}

# --- 4) Xem nhanh: hệ số RCPU & RGOLD×RCPU (mxreg5, mxreg6) ---
rcpu_tbl <- lapply(garch_out, function(x){
  subset(x, param %in% c("mxreg5","mxreg6")) %>%
    mutate(Sig = star(pval)) %>%
    transmute(Market, param, Estimate=round(Estimate,3), pval=round(pval,3), Sig)
}) %>% bind_rows()

cat("\n=== GARCH-X: RCPU và RGOLD×RCPU (với stars) ===\n")
print(rcpu_tbl)











# ==================================================
#  THỐNG KÊ MÔ TẢ & MA TRẬN TƯƠNG QUAN
# ==================================================

library(psych)  # để mô tả gọn đẹp

# --- chọn các biến cần thống kê ---
vars_desc <- c("RSET","RVNI","RSTI","RKLSE","RJKSE","RPSE","RTHB","RVND","RSGD","RIDR","RMYR","RPHP","RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")

# --- lọc dữ liệu, bỏ NA ---
df_desc <- cpu %>% dplyr::select(any_of(vars_desc))

# --- 1) Thống kê mô tả ---
cat("\n==================== Thống kê mô tả ====================\n")
desc_tbl <- psych::describe(df_desc)[,c("n","mean","sd","min","max","skew","kurtosis")]
print(round(desc_tbl,6))

# --- 2) Ma trận tương quan Pearson ---
cat("\n==================== Ma trận tương quan ====================\n")
cor_mat <- cor(df_desc, use="pairwise.complete.obs", method="pearson")
print(round(cor_mat,3))


# ==================================================
# ASEAN-7: OLS (HAC) + GARCH-X (AR(1)-sGARCH(1,1)-t)
# So sánh OOS: GARCH-X vs OLS vs RW vs ARIMA
# Biến exog: FX, RWTI, RVIX, RGOLD, RCPU, RGOLD×RCPU
# ==================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(zoo)
  library(xts); library(rugarch); library(forecast)
  library(lmtest); library(sandwich); library(FinTS)
})

options(stringsAsFactors = FALSE)

# -----------------------------
# 1) Chuẩn bị dữ liệu & hàm phụ
# -----------------------------
cpu <- readxl::read_excel("C:/Users/trong/OneDrive/HUYNHTT/Year 2025/Thang 11/cpu.xlsx")
cpu <- cpu %>%
  mutate(
    ym   = as.integer(ym),
    date = as.Date(as.yearmon(as.character(ym), "%Y%m")),
    RGOLD_RCPU = RGOLD * RCPU
  ) %>%
  arrange(date)

# Danh sách thị trường phụ thuộc
dep_vars <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")

# Mapping FX theo thị trường
fx_map <- c(
  RSET  = "RTHB", # Thailand
  RVNI  = "RVND", # Vietnam (HOSE)
  RHNX  = "RVND", # Vietnam (HNX)
  RSTI  = "RSGD", # Singapore
  RKLSE = "RMYR", # Malaysia
  RJKSE = "RIDR", # Indonesia
  RPSE  = "RPHP"  # Philippines
)

# Biến ngoại sinh dùng chung
other_exog <- c("RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")

# Hàm lỗi
MAE  <- function(e) mean(abs(e), na.rm=TRUE)
RMSE <- function(e) sqrt(mean(e^2, na.rm=TRUE))

# -----------------------------
# 2) Hàm chạy từng thị trường
#    (OLS in-sample + GARCH-X in-sample + Rolling OOS 4 mô hình)
# -----------------------------
run_country_model <- function(yvar, data, roll_window=48) {
  message("== Running for ", yvar, " ==")
  
  # Xác định biến FX cho yvar
  fx_var    <- unname(fx_map[[yvar]])
  exog_vars <- c(fx_var, other_exog)
  
  # Chỉ giữ cột cần thiết & loại NA
  df <- data %>%
    dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>%
    drop_na()
  
  # -------- OLS in-sample (HAC) --------
  f_ols <- as.formula(paste(
    yvar, "~", paste(c(fx_var,"RWTI","RVIX","RGOLD","RCPU"), collapse=" + "),
    "+ RGOLD:RCPU"
  ))
  reg_ols <- lm(f_ols, data=df)
  # HAC (Newey-West) tự động chọn lag
  bw <- sandwich::bwNeweyWest(reg_ols)
  se <- sandwich::NeweyWest(reg_ols, lag=bw, prewhite=TRUE, adjust=TRUE)
  ols_sum <- coeftest(reg_ols, vcov.=se)
  
  # Một số test chẩn đoán
  arch_p <- tryCatch(FinTS::ArchTest(residuals(reg_ols), lags=12)$p.value, error=function(e) NA_real_)
  bp_p   <- tryCatch(lmtest::bptest(reg_ols)$p.value,                     error=function(e) NA_real_)
  
  # -------- GARCH-X in-sample --------
  y_xts <- xts(df[[yvar]], order.by=df$date)
  X_mat <- as.matrix(df[, exog_vars])
  spec <- ugarchspec(
    mean.model     = list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X_mat),
    variance.model = list(model="sGARCH", garchOrder=c(1,1)),
    distribution.model="std"
  )
  fit_garch_ins <- tryCatch(
    ugarchfit(spec, data=y_xts, solver="hybrid", solver.control=list(trace=0)),
    error=function(e) NULL
  )
  
  # -------- Rolling OOS (48M) --------
  N <- nrow(df)
  if (N <= roll_window) {
    warning("Chuỗi ", yvar, " quá ngắn cho rolling OOS.")
    return(list(
      yvar=yvar, FX=fx_var,
      OLS=ols_sum, pARCH=arch_p, pBP=bp_p,
      GARCHX_summary=if(!is.null(fit_garch_ins)) coef(fit_garch_ins) else NA,
      OOS_MAE=NA, OOS_RMSE=NA
    ))
  }
  
  oos <- vector("list", N - roll_window)
  k <- 0
  
  for (i in seq(roll_window, N-1)) {
    tr <- (i - roll_window + 1):i
    te <- i + 1
    
    y_tr <- df[[yvar]][tr]
    y_te <- df[[yvar]][te]
    
    # Ma trận exog theo spec GARCH (gồm cả tương tác đã có trong RGOLD_RCPU)
    X_tr <- as.matrix(df[tr, exog_vars])
    X_te <- as.matrix(df[te,  exog_vars])
    
    # ===== 1) GARCH-X forecast =====
    sp <- ugarchspec(
      mean.model=list(armaOrder=c(1,0), include.mean=TRUE, external.regressors=X_tr),
      variance.model=list(model="sGARCH", garchOrder=c(1,1)),
      distribution.model="std"
    )
    fit_r <- tryCatch(
      ugarchfit(sp, data=y_tr, solver="hybrid", solver.control=list(trace=0)),
      error=function(e) NULL
    )
    pred_gx <- if (!is.null(fit_r)) {
      fc <- tryCatch(
        ugarchforecast(fit_r, n.ahead=1,
                       external.forecasts=list(mreg=matrix(X_te, nrow=1))),
        error=function(e) NULL
      )
      if (is.null(fc)) NA_real_ else as.numeric(fitted(fc))
    } else NA_real_
    
    # ===== 2) OLS forecast (cùng biến ngoại sinh + tương tác) =====
    # Dùng khung dữ liệu train/test cho lm()
    df_tr_ols <- df[tr, c(yvar, fx_var,"RWTI","RVIX","RGOLD","RCPU")]
    df_te_ols <- df[te,  c(fx_var,"RWTI","RVIX","RGOLD","RCPU")]
    # Tạo biến tương tác RGOLD:RCPU cho OLS OOS (nếu muốn explicit)
    df_tr_ols$RGOLD_RCPU <- df_tr_ols$RGOLD * df_tr_ols$RCPU
    df_te_ols$RGOLD_RCPU <- df_te_ols$RGOLD * df_te_ols$RCPU
    
    f_ols_oos <- as.formula(paste(
      yvar, "~", paste(c(fx_var,"RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU"), collapse=" + ")
    ))
    fit_ols_oos <- tryCatch(lm(f_ols_oos, data=df_tr_ols), error=function(e) NULL)
    pred_ols <- if (!is.null(fit_ols_oos)) {
      as.numeric(predict(fit_ols_oos, newdata=df_te_ols))
    } else NA_real_
    
    # ===== 3) Random Walk =====
    pred_rw <- tail(y_tr, 1)
    
    # ===== 4) ARIMA =====
    fit_a <- tryCatch(
      suppressWarnings(auto.arima(y_tr, ic="aic", seasonal=FALSE)),
      error=function(e) NULL
    )
    pred_a <- if (!is.null(fit_a)) as.numeric(forecast(fit_a, h=1)$mean) else NA_real_
    
    k <- k + 1
    oos[[k]] <- data.frame(
      date = df$date[te],
      y_true = y_te,
      yhat_garchx = pred_gx,
      yhat_ols    = pred_ols,
      yhat_rw     = pred_rw,
      yhat_arima  = pred_a
    )
  }
  
  oos_df <- bind_rows(oos) %>%
    mutate(
      e_gx  = y_true - yhat_garchx,
      e_ols = y_true - yhat_ols,
      e_rw  = y_true - yhat_rw,
      e_ar  = y_true - yhat_arima
    )
  
  oos_mae  <- c(GARCHX = MAE(oos_df$e_gx),
                OLS    = MAE(oos_df$e_ols),
                RW     = MAE(oos_df$e_rw),
                ARIMA  = MAE(oos_df$e_ar))
  oos_rmse <- c(GARCHX = RMSE(oos_df$e_gx),
                OLS    = RMSE(oos_df$e_ols),
                RW     = RMSE(oos_df$e_rw),
                ARIMA  = RMSE(oos_df$e_ar))
  
  list(
    yvar = yvar,
    FX   = fx_var,
    OLS  = ols_sum,                    # Bảng hệ số OLS (HAC) in-sample
    pARCH= arch_p,                     # P-value ARCH LM
    pBP  = bp_p,                       # P-value Breusch-Pagan
    GARCHX_summary = if(!is.null(fit_garch_ins)) coef(fit_garch_ins) else NA,
    OOS_MAE  = oos_mae,
    OOS_RMSE = oos_rmse,
    OOS_raw  = oos_df                  # lưu kèm lỗi dự báo từng thời điểm (tiện báo cáo)
  )
}

# -----------------------------
# 3) Chạy tất cả thị trường
# -----------------------------
results <- lapply(dep_vars, run_country_model, data=cpu, roll_window=48)

# -----------------------------
# 4) Tóm tắt bảng MAE/RMSE (4 mô hình)
# -----------------------------
mae_tbl <- do.call(rbind, lapply(results, function(x)
  data.frame(Market=x$yvar, t(x$OOS_MAE), row.names=NULL)))
rmse_tbl <- do.call(rbind, lapply(results, function(x)
  data.frame(Market=x$yvar, t(x$OOS_RMSE), row.names=NULL)))

cat("\n== MAE across markets (GARCHX, OLS, RW, ARIMA) ==\n"); print(mae_tbl, row.names=FALSE)
cat("\n== RMSE across markets (GARCHX, OLS, RW, ARIMA) ==\n"); print(rmse_tbl, row.names=FALSE)

# (Tùy chọn) Lưu CSV kết quả tổng hợp
# write.csv(mae_tbl,  "mae_tbl.csv",  row.names=FALSE)
# write.csv(rmse_tbl, "rmse_tbl.csv", row.names=FALSE)

# (Tùy chọn) Lấy OOS raw cho 1 thị trường, ví dụ RVNI:
# oos_rvni <- results[[which(dep_vars=="RVNI")]]$OOS_raw
# head(oos_rvni)

############## Thêm code chỉnh sửa bản R0

# ==================================================
# ECONOMIC SIGNIFICANCE of Gold × CPU
# Ý tưởng:
# Economic significance = beta_(Gold×CPU) * SD(Gold×CPU)
# -> diễn giải: khi Gold×CPU tăng 1 độ lệch chuẩn, lợi suất thị trường thay đổi bao nhiêu
# Thêm cột chuẩn hoá theo SD(Y): (% of SD(Y))
# ==================================================

cat("\n================ ECONOMIC SIGNIFICANCE: Gold × CPU ================\n")

econ_sig_tbl <- lapply(dep_vars, function(yvar) {
  
  fx_var <- fx_map[[yvar]]
  exog_vars <- c(fx_var, "RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")
  
  df <- cpu %>%
    dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>%
    tidyr::drop_na()
  
  # --- SD của interaction và SD của biến phụ thuộc ---
  sd_interaction <- sd(df$RGOLD_RCPU, na.rm = TRUE)
  sd_y <- sd(df[[yvar]], na.rm = TRUE)
  
  # =========================
  # 1) OLS economic significance
  # =========================
  beta_ols <- summary_ols$RGOLD_RCPU[summary_ols$Market == yvar]
  p_ols    <- summary_ols$RGOLD_RCPU_p[summary_ols$Market == yvar]
  
  econ_ols_1sd <- beta_ols * sd_interaction
  econ_ols_pct_sdY <- ifelse(is.na(sd_y) || sd_y == 0, NA_real_,
                             100 * econ_ols_1sd / sd_y)
  
  # =========================
  # 2) GARCH-X economic significance
  # mxreg6 = RGOLD×RCPU theo đúng thứ tự exog hiện tại
  # =========================
  res_g <- garch_results[[yvar]]
  
  beta_garch <- if (!is.null(res_g)) {
    v <- res_g$tidy$Estimate[match("mxreg6", res_g$tidy$param)]
    ifelse(is.na(v), NA_real_, as.numeric(v))
  } else NA_real_
  
  p_garch <- if (!is.null(res_g)) {
    v <- res_g$tidy$p_rob[match("mxreg6", res_g$tidy$param)]
    ifelse(is.na(v), NA_real_, as.numeric(v))
  } else NA_real_
  
  econ_garch_1sd <- beta_garch * sd_interaction
  econ_garch_pct_sdY <- ifelse(is.na(sd_y) || sd_y == 0, NA_real_,
                               100 * econ_garch_1sd / sd_y)
  
  data.frame(
    Market = yvar,
    SD_RGOLD_RCPU = sd_interaction,
    SD_Y = sd_y,
    
    OLS_beta_RGOLD_RCPU = beta_ols,
    OLS_pvalue = p_ols,
    OLS_EconSig_1SD = econ_ols_1sd,
    OLS_EconSig_pct_SDY = econ_ols_pct_sdY,
    
    GARCH_beta_RGOLD_RCPU = beta_garch,
    GARCH_pvalue_rob = p_garch,
    GARCH_EconSig_1SD = econ_garch_1sd,
    GARCH_EconSig_pct_SDY = econ_garch_pct_sdY,
    
    check.names = FALSE
  )
}) %>% bind_rows()

print(econ_sig_tbl, row.names = FALSE, digits = 6)

# ---------- Bảng đẹp để đưa vào paper ----------
star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < .001, "***",
                ifelse(p < .01,  "**",
                       ifelse(p < .05,  "*",
                              ifelse(p < .10,  ".", "")))))
}

econ_sig_pretty <- econ_sig_tbl %>%
  mutate(
    OLS_beta = sprintf("%.4f%s", OLS_beta_RGOLD_RCPU, star(OLS_pvalue)),
    GARCH_beta = sprintf("%.4f%s", GARCH_beta_RGOLD_RCPU, star(GARCH_pvalue_rob)),
    OLS_EconSig_1SD = sprintf("%.4f", OLS_EconSig_1SD),
    GARCH_EconSig_1SD = sprintf("%.4f", GARCH_EconSig_1SD),
    OLS_EconSig_pct_SDY = sprintf("%.2f%%", OLS_EconSig_pct_SDY),
    GARCH_EconSig_pct_SDY = sprintf("%.2f%%", GARCH_EconSig_pct_SDY)
  ) %>%
  dplyr::select(
    Market,
    OLS_beta,
    OLS_EconSig_1SD,
    OLS_EconSig_pct_SDY,
    GARCH_beta,
    GARCH_EconSig_1SD,
    GARCH_EconSig_pct_SDY
  )

cat("\n================ Paper-ready: Economic Significance of Gold × CPU ================\n")
print(econ_sig_pretty, row.names = FALSE, right = FALSE)

# (Tuỳ chọn) lưu file
# write.csv(econ_sig_tbl, "ASEAN7_EconomicSignificance_GOLD_CPU.csv", row.names = FALSE)
# write.csv(econ_sig_pretty, "ASEAN7_EconomicSignificance_GOLD_CPU_pretty.csv", row.names = FALSE)



# ==================================================
# ROBUSTNESS BY SUBSAMPLES
# OLS (HAC) + GARCH-X
# Samples: Full, Post-2005, Post-2015
# Focus: sign consistency of RGOLD, RCPU, RGOLD×RCPU
# ==================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lmtest)
  library(sandwich)
  library(FinTS)
  library(xts)
  library(rugarch)
})

# --------------------------------------------------
# 0) Nếu chưa có biến interaction thì tạo lại cho chắc
# --------------------------------------------------
cpu <- cpu %>%
  mutate(
    RGOLD_RCPU = RGOLD * RCPU
  ) %>%
  arrange(date)

dep_vars <- c("RSET","RVNI","RHNX","RSTI","RKLSE","RJKSE","RPSE")

fx_map <- c(
  RSET  = "RTHB",
  RVNI  = "RVND",
  RHNX  = "RVND",
  RSTI  = "RSGD",
  RKLSE = "RMYR",
  RJKSE = "RIDR",
  RPSE  = "RPHP"
)

other_exog <- c("RWTI","RVIX","RGOLD","RCPU","RGOLD_RCPU")

star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < .001, "***",
                ifelse(p < .01,  "**",
                       ifelse(p < .05,  "*",
                              ifelse(p < .10, ".", "")))))
}

sign_label <- function(x) {
  ifelse(is.na(x), "NA",
         ifelse(x > 0, "Positive",
                ifelse(x < 0, "Negative", "Zero")))
}

# --------------------------------------------------
# 1) Định nghĩa mẫu
# --------------------------------------------------
sample_list <- list(
  Full      = cpu,
  Post_2005 = cpu %>% filter(date >= as.Date("2005-01-01")),
  Post_2015 = cpu %>% filter(date >= as.Date("2015-01-01"))
)

# --------------------------------------------------
# 2) Helper: tidy OLS HAC
# --------------------------------------------------
tidy_coeftest <- function(ct_mat) {
  data.frame(
    term      = rownames(ct_mat),
    Estimate  = as.numeric(ct_mat[, 1]),
    Std.Error = as.numeric(ct_mat[, 2]),
    t.value   = as.numeric(ct_mat[, 3]),
    p.value   = as.numeric(ct_mat[, 4]),
    row.names = NULL,
    check.names = FALSE
  )
}

# --------------------------------------------------
# 3) Helper: tidy GARCH
# --------------------------------------------------
tidy_garch <- function(fit) {
  est <- coef(fit)
  vcr <- tryCatch(vcov(fit, robust = TRUE), error = function(e) NULL)
  ser <- if (!is.null(vcr)) sqrt(diag(vcr)) else rep(NA_real_, length(est))
  zr  <- est / ser
  pr  <- 2 * pnorm(abs(zr), lower.tail = FALSE)
  
  data.frame(
    param    = names(est),
    Estimate = as.numeric(est),
    SE_rob   = as.numeric(ser),
    z_rob    = as.numeric(zr),
    p_rob    = as.numeric(pr),
    row.names = NULL,
    check.names = FALSE
  )
}

# --------------------------------------------------
# 4) Chạy OLS + GARCH-X theo sample × market
# --------------------------------------------------
robust_results <- list()

for (sname in names(sample_list)) {
  cat("\n====================================================\n")
  cat("Running sample:", sname, "\n")
  cat("====================================================\n")
  
  dat_s <- sample_list[[sname]]
  robust_results[[sname]] <- list()
  
  for (yvar in dep_vars) {
    cat("\n------------------------------\n")
    cat("Market:", yvar, "| Sample:", sname, "\n")
    cat("------------------------------\n")
    
    fx_var    <- unname(fx_map[[yvar]])
    exog_vars <- c(fx_var, other_exog)
    
    df <- dat_s %>%
      dplyr::select(date, dplyr::all_of(yvar), dplyr::all_of(exog_vars)) %>%
      tidyr::drop_na()
    
    # Bỏ qua nếu mẫu quá ngắn
    if (nrow(df) < 40) {
      cat("Sample too short, skipped.\n")
      robust_results[[sname]][[yvar]] <- list(
        OLS = NULL,
        GARCH = NULL,
        N = nrow(df)
      )
      next
    }
    
    # =========================
    # OLS (HAC)
    # =========================
    f_ols <- as.formula(paste(
      yvar, "~",
      paste(c(fx_var, "RWTI","RVIX","RGOLD","RCPU"), collapse = " + "),
      "+ RGOLD:RCPU"
    ))
    
    reg_ols <- lm(f_ols, data = df)
    bw <- sandwich::bwNeweyWest(reg_ols)
    hac <- sandwich::NeweyWest(reg_ols, lag = bw, prewhite = TRUE, adjust = TRUE)
    ct <- lmtest::coeftest(reg_ols, vcov. = hac)
    ols_tbl <- tidy_coeftest(ct)
    
    bp_p   <- tryCatch(lmtest::bptest(reg_ols)$p.value, error = function(e) NA_real_)
    arch_p <- tryCatch(FinTS::ArchTest(residuals(reg_ols), lags = 12)$p.value, error = function(e) NA_real_)
    
    # =========================
    # GARCH-X
    # mxreg1..mxreg6 = FX, RWTI, RVIX, RGOLD, RCPU, RGOLD_RCPU
    # =========================
    y_xts <- xts(df[[yvar]], order.by = df$date)
    X_mat <- as.matrix(df[, exog_vars])
    
    spec <- ugarchspec(
      mean.model = list(
        armaOrder = c(1,0),
        include.mean = TRUE,
        external.regressors = X_mat
      ),
      variance.model = list(
        model = "sGARCH",
        garchOrder = c(1,1)
      ),
      distribution.model = "std"
    )
    
    fit <- tryCatch(
      ugarchfit(spec, data = y_xts, solver = "hybrid", solver.control = list(trace = 0)),
      error = function(e) NULL
    )
    
    garch_tbl <- if (!is.null(fit)) tidy_garch(fit) else NULL
    
    robust_results[[sname]][[yvar]] <- list(
      OLS = ols_tbl,
      GARCH = garch_tbl,
      N = nrow(df),
      OLS_BP_p = bp_p,
      OLS_ARCH_p = arch_p
    )
  }
}

# --------------------------------------------------
# 5) Trích bảng dấu cho reviewer
# --------------------------------------------------
extract_sign_summary <- function(res_list, sample_name, market_name) {
  res <- res_list[[sample_name]][[market_name]]
  
  if (is.null(res) || is.null(res$OLS)) {
    return(data.frame(
      Sample = sample_name,
      Market = market_name,
      N = ifelse(is.null(res$N), NA, res$N),
      
      OLS_RGOLD = NA, OLS_RGOLD_p = NA, OLS_RGOLD_sign = NA,
      OLS_RCPU = NA, OLS_RCPU_p = NA, OLS_RCPU_sign = NA,
      OLS_RGOLD_RCPU = NA, OLS_RGOLD_RCPU_p = NA, OLS_RGOLD_RCPU_sign = NA,
      
      GARCH_RGOLD = NA, GARCH_RGOLD_p = NA, GARCH_RGOLD_sign = NA,
      GARCH_RCPU = NA, GARCH_RCPU_p = NA, GARCH_RCPU_sign = NA,
      GARCH_RGOLD_RCPU = NA, GARCH_RGOLD_RCPU_p = NA, GARCH_RGOLD_RCPU_sign = NA
    ))
  }
  
  ols_tbl <- res$OLS
  garch_tbl <- res$GARCH
  
  # OLS
  ols_rgold <- ols_tbl$Estimate[match("RGOLD", ols_tbl$term)]
  ols_rgold_p <- ols_tbl$p.value[match("RGOLD", ols_tbl$term)]
  
  ols_rcpu <- ols_tbl$Estimate[match("RCPU", ols_tbl$term)]
  ols_rcpu_p <- ols_tbl$p.value[match("RCPU", ols_tbl$term)]
  
  ols_int <- ols_tbl$Estimate[match("RGOLD:RCPU", ols_tbl$term)]
  ols_int_p <- ols_tbl$p.value[match("RGOLD:RCPU", ols_tbl$term)]
  
  # GARCH
  g_rgold <- if (!is.null(garch_tbl)) garch_tbl$Estimate[match("mxreg4", garch_tbl$param)] else NA_real_
  g_rgold_p <- if (!is.null(garch_tbl)) garch_tbl$p_rob[match("mxreg4", garch_tbl$param)] else NA_real_
  
  g_rcpu <- if (!is.null(garch_tbl)) garch_tbl$Estimate[match("mxreg5", garch_tbl$param)] else NA_real_
  g_rcpu_p <- if (!is.null(garch_tbl)) garch_tbl$p_rob[match("mxreg5", garch_tbl$param)] else NA_real_
  
  g_int <- if (!is.null(garch_tbl)) garch_tbl$Estimate[match("mxreg6", garch_tbl$param)] else NA_real_
  g_int_p <- if (!is.null(garch_tbl)) garch_tbl$p_rob[match("mxreg6", garch_tbl$param)] else NA_real_
  
  data.frame(
    Sample = sample_name,
    Market = market_name,
    N = res$N,
    
    OLS_RGOLD = ols_rgold,
    OLS_RGOLD_p = ols_rgold_p,
    OLS_RGOLD_sign = sign_label(ols_rgold),
    
    OLS_RCPU = ols_rcpu,
    OLS_RCPU_p = ols_rcpu_p,
    OLS_RCPU_sign = sign_label(ols_rcpu),
    
    OLS_RGOLD_RCPU = ols_int,
    OLS_RGOLD_RCPU_p = ols_int_p,
    OLS_RGOLD_RCPU_sign = sign_label(ols_int),
    
    GARCH_RGOLD = g_rgold,
    GARCH_RGOLD_p = g_rgold_p,
    GARCH_RGOLD_sign = sign_label(g_rgold),
    
    GARCH_RCPU = g_rcpu,
    GARCH_RCPU_p = g_rcpu_p,
    GARCH_RCPU_sign = sign_label(g_rcpu),
    
    GARCH_RGOLD_RCPU = g_int,
    GARCH_RGOLD_RCPU_p = g_int_p,
    GARCH_RGOLD_RCPU_sign = sign_label(g_int),
    
    check.names = FALSE
  )
}

sign_summary_tbl <- bind_rows(
  lapply(names(sample_list), function(sname) {
    bind_rows(lapply(dep_vars, function(yvar) {
      extract_sign_summary(robust_results, sname, yvar)
    }))
  })
)

cat("\n================ SIGN SUMMARY BY SAMPLE ================\n")
print(sign_summary_tbl, row.names = FALSE, digits = 6)

# --------------------------------------------------
# 6) Bảng đẹp cho báo cáo
# --------------------------------------------------
pretty_sign_tbl <- sign_summary_tbl %>%
  mutate(
    OLS_RGOLD = sprintf("%.3f%s (%s)", OLS_RGOLD, star(OLS_RGOLD_p), OLS_RGOLD_sign),
    OLS_RCPU = sprintf("%.3f%s (%s)", OLS_RCPU, star(OLS_RCPU_p), OLS_RCPU_sign),
    OLS_RGOLD_RCPU = sprintf("%.3f%s (%s)", OLS_RGOLD_RCPU, star(OLS_RGOLD_RCPU_p), OLS_RGOLD_RCPU_sign),
    
    GARCH_RGOLD = sprintf("%.3f%s (%s)", GARCH_RGOLD, star(GARCH_RGOLD_p), GARCH_RGOLD_sign),
    GARCH_RCPU = sprintf("%.3f%s (%s)", GARCH_RCPU, star(GARCH_RCPU_p), GARCH_RCPU_sign),
    GARCH_RGOLD_RCPU = sprintf("%.3f%s (%s)", GARCH_RGOLD_RCPU, star(GARCH_RGOLD_RCPU_p), GARCH_RGOLD_RCPU_sign)
  ) %>%
  dplyr::select(
    Sample, Market, N,
    OLS_RGOLD, OLS_RCPU, OLS_RGOLD_RCPU,
    GARCH_RGOLD, GARCH_RCPU, GARCH_RGOLD_RCPU
  )

cat("\n================ PAPER-READY SIGN TABLE ================\n")
print(pretty_sign_tbl, row.names = FALSE, right = FALSE)

# --------------------------------------------------
# 7) Tóm tắt consistency của dấu qua các mẫu
# --------------------------------------------------
consistency_tbl <- sign_summary_tbl %>%
  group_by(Market) %>%
  summarise(
    OLS_RGOLD_signs = paste(unique(na.omit(OLS_RGOLD_sign)), collapse = ", "),
    OLS_RCPU_signs = paste(unique(na.omit(OLS_RCPU_sign)), collapse = ", "),
    OLS_RGOLD_RCPU_signs = paste(unique(na.omit(OLS_RGOLD_RCPU_sign)), collapse = ", "),
    
    GARCH_RGOLD_signs = paste(unique(na.omit(GARCH_RGOLD_sign)), collapse = ", "),
    GARCH_RCPU_signs = paste(unique(na.omit(GARCH_RCPU_sign)), collapse = ", "),
    GARCH_RGOLD_RCPU_signs = paste(unique(na.omit(GARCH_RGOLD_RCPU_sign)), collapse = ", "),
    .groups = "drop"
  )

cat("\n================ SIGN CONSISTENCY ACROSS SUBSAMPLES ================\n")
print(consistency_tbl, row.names = FALSE)

# (Tuỳ chọn) lưu file
# write.csv(sign_summary_tbl, "robust_subsample_sign_summary.csv", row.names = FALSE)
# write.csv(pretty_sign_tbl, "robust_subsample_sign_pretty.csv", row.names = FALSE)
# write.csv(consistency_tbl, "robust_subsample_sign_consistency.csv", row.names = FALSE)















