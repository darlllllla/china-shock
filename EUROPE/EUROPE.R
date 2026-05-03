# ============================================================
#  China Shock Analysis — European NUTS-2 regions
#  Methodology: Autor, Dorn & Hanson (2013), adapted for Europe
#  Bartik (shift-share) IPW + IV-2SLS
#  Periods: 2000-2007, 2008-2011, 2012-2019
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(AER)
library(lmtest)
library(sandwich)
library(stargazer)
library(ggplot2)
library(ggrepel)

# Patch stargazer for compatibility with R >= 4.x.
# Fixes vectorised if() warnings in is.na(s) and s == "" checks.
local({
  ns <- asNamespace("stargazer")
  for (nm in ls(ns, all.names = TRUE)) {
    fn <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (!is.function(fn)) next
    src <- tryCatch(deparse(body(fn)), error = function(e) NULL)

    has_isna  <- any(grepl("if \\(is\\.na\\(s\\)\\)", src))
    has_empty <- any(grepl("if \\(s == \"\"\\)", src))
    if (!has_isna && !has_empty) next

    src_new <- src
    if (has_isna)
      src_new <- gsub("if \\(is\\.na\\(s\\)\\)", "if (all(is.na(s)))", src_new)
    if (has_empty)
      src_new <- gsub("if \\(s == \"\"\\)", "if (all(s == \"\"))", src_new)

    new_body <- tryCatch(
      parse(text = paste(src_new, collapse = "\n"))[[1]],
      error = function(e) NULL
    )
    if (is.null(new_body)) next

    tryCatch({
      unlockBinding(nm, ns)
      body(fn) <- new_body
      assign(nm, fn, envir = ns)
      lockBinding(nm, ns)
    }, error = function(e) NULL)
  }
})

# ============================================================
# 1. COUNTRY CONFIGURATION
# ============================================================

FILE         <- file.choose()   # select the country .xlsx file
COUNTRY      <- "COUNTRY"
COUNTRY_CODE <- "COUNTRY CODE"

# ============================================================
# 2. DATA LOADING
# ============================================================

employment <- read_excel(FILE, sheet = "employment") |>
  setNames(c("zone", "year", "industry", "emp")) |>
  mutate(year     = as.integer(year),
         emp      = as.numeric(emp),
         industry = trimws(industry),
         zone     = trimws(zone))

import_de <- read_excel(FILE, sheet = "import") |>
  setNames(c("year", "industry", "import")) |>
  mutate(year     = as.integer(year),
         import   = as.numeric(import),
         industry = trimws(industry))

industries_emp <- sort(unique(employment$industry))
industries_imp <- sort(unique(import_de$industry))
missing_in_imp <- setdiff(industries_emp, industries_imp)
missing_in_emp <- setdiff(industries_imp, industries_emp)
if (length(missing_in_imp) > 0)
  warning("Industries in employment but NOT in import: ", paste(missing_in_imp, collapse = ", "))
if (length(missing_in_emp) > 0)
  message("Industries in import but NOT in employment: ", paste(missing_in_emp, collapse = ", "))

import_oc_raw <- read_excel(FILE, sheet = "import OC",
                            col_names = FALSE, .name_repair = "minimal")

# Other-country imports used as IV instrument (Bartik-style exclusion restriction)
countries_oc <- c("Australia", "Finland", "Japan",
                  "Norway", "NewZealand", "Sweden")
col_starts   <- c(1L, 5L, 9L, 13L, 17L, 21L)

import_oc_list <- vector("list", length(countries_oc))
for (k in seq_along(countries_oc)) {
  cs  <- col_starts[k]
  blk <- import_oc_raw[3:nrow(import_oc_raw), cs:(cs + 2L)]
  colnames(blk) <- c("year", "industry", "import")
  blk <- blk |>
    mutate(year     = suppressWarnings(as.integer(year)),
           import   = suppressWarnings(as.numeric(import)),
           industry = trimws(as.character(industry)),
           country  = countries_oc[k]) |>
    filter(!is.na(year))
  import_oc_list[[k]] <- blk
}

import_oc_agg <- bind_rows(import_oc_list) |>
  group_by(year, industry) |>
  summarise(import_oc = sum(import, na.rm = TRUE), .groups = "drop")

# ============================================================
# 2b. MANUFACTURING INDUSTRY FILTER
# ============================================================

MFG_INDUSTRIES <- c(
  "c-e"
)

# ============================================================
# 3. BUILD REGRESSION DATASET
# ============================================================
# industry_filter:
#   NULL             -> full economy (default)
#   character vector -> selected industries only (e.g., MFG_INDUSTRIES)

build_regression_data <- function(t0, t1, industry_filter = NULL) {

  # Full employment data is always used in the L_r0 denominator for s_ri
  emp_t0_full <- filter(employment, year == t0)
  emp_t1_full <- filter(employment, year == t1)

  emp_t0 <- if (!is.null(industry_filter)) filter(emp_t0_full, industry %in% industry_filter) else emp_t0_full
  emp_t1 <- if (!is.null(industry_filter)) filter(emp_t1_full, industry %in% industry_filter) else emp_t1_full

  imp_data    <- if (!is.null(industry_filter)) filter(import_de,    industry %in% industry_filter) else import_de
  imp_oc_data <- if (!is.null(industry_filter)) filter(import_oc_agg, industry %in% industry_filter) else import_oc_agg

  if (!is.null(industry_filter) && nrow(emp_t0) == 0)
    stop("industry_filter matched no industries. ",
         "Check with: cat(sort(unique(employment$industry)), sep='\\n')")

  L_r0_filt <- emp_t0 |> group_by(zone) |>
    summarise(L_r0_filt = sum(emp, na.rm = TRUE), .groups = "drop")
  L_r1_filt <- emp_t1 |> group_by(zone) |>
    summarise(L_r1_filt = sum(emp, na.rm = TRUE), .groups = "drop")

  dep_var <- L_r0_filt |>
    inner_join(L_r1_filt, by = "zone") |>
    mutate(delta_lnL = log(L_r1_filt) - log(L_r0_filt))

  lost_t1 <- setdiff(L_r0_filt$zone, L_r1_filt$zone)
  lost_t0 <- setdiff(L_r1_filt$zone, L_r0_filt$zone)
  if (length(lost_t1) > 0)
    warning("Zones present in t0 but missing in t1: ", paste(lost_t1, collapse = ", "))
  if (length(lost_t0) > 0)
    warning("Zones present in t1 but missing in t0: ", paste(lost_t0, collapse = ", "))

  # s_ri denominator: total employment across all sectors (not filtered)
  L_r0_total <- emp_t0_full |> group_by(zone) |>
    summarise(L_r0 = sum(emp, na.rm = TRUE), .groups = "drop")

  s_ri <- emp_t0 |>
    inner_join(L_r0_total, by = "zone") |>
    mutate(s_ri = emp / L_r0) |>
    select(zone, industry, s_ri)

  L_i0 <- emp_t0 |> group_by(industry) |>
    summarise(L_i0 = sum(emp, na.rm = TRUE), .groups = "drop")

  delta_M_de <- imp_data |>
    filter(year %in% c(t0, t1)) |>
    pivot_wider(names_from = year, values_from = import, names_prefix = "M_") |>
    mutate(delta_M = .data[[paste0("M_", t1)]] - .data[[paste0("M_", t0)]]) |>
    select(industry, delta_M)

  shock_i <- L_i0 |>
    left_join(delta_M_de, by = "industry") |>
    mutate(shock = delta_M / L_i0) |>
    select(industry, shock)

  n_na_shock <- sum(is.na(shock_i$shock))
  if (n_na_shock > 0)
    warning(sprintf("Period %d-%d: %d industries have NA shock_i", t0, t1, n_na_shock))

  IPW_df <- s_ri |>
    left_join(shock_i, by = "industry") |>
    group_by(zone) |>
    summarise(IPW = sum(s_ri * shock, na.rm = TRUE) / 1000, .groups = "drop")

  delta_M_oc <- imp_oc_data |>
    filter(year %in% c(t0, t1)) |>
    pivot_wider(names_from = year, values_from = import_oc, names_prefix = "M_") |>
    mutate(delta_M_oc = .data[[paste0("M_", t1)]] - .data[[paste0("M_", t0)]]) |>
    select(industry, delta_M_oc)

  shock_oc <- L_i0 |>
    left_join(delta_M_oc, by = "industry") |>
    mutate(shock_oc = delta_M_oc / L_i0) |>
    select(industry, shock_oc)

  IPW_OC_df <- s_ri |>
    left_join(shock_oc, by = "industry") |>
    group_by(zone) |>
    summarise(IPW_OC = sum(s_ri * shock_oc, na.rm = TRUE) / 1000, .groups = "drop")

  out <- dep_var |>
    left_join(IPW_df,    by = "zone") |>
    left_join(IPW_OC_df, by = "zone") |>
    filter(!is.na(delta_lnL), !is.na(IPW), !is.na(IPW_OC),
           is.finite(delta_lnL), is.finite(IPW), is.finite(IPW_OC))

  all_zones <- unique(emp_t0$zone)
  dropped   <- setdiff(all_zones, out$zone)
  if (length(dropped) > 0)
    warning(sprintf("Period %d-%d: %d zone(s) dropped: %s",
                    t0, t1, length(dropped), paste(dropped, collapse = ", ")))

  return(out)
}

# ============================================================
# 4. BUILD DATASETS — FULL ECONOMY
# ============================================================

df_p1_all <- build_regression_data(t0 = 2000, t1 = 2007)
df_p2_all <- build_regression_data(t0 = 2008, t1 = 2011)
df_p3_all <- build_regression_data(t0 = 2012, t1 = 2019)

cat("\n=== FULL ECONOMY ===\n")
cat(COUNTRY, "- Period 1 (2000-2007): N =", nrow(df_p1_all), "zones\n")
cat(COUNTRY, "- Period 2 (2008-2011): N =", nrow(df_p2_all), "zones\n")
cat(COUNTRY, "- Period 3 (2012-2019): N =", nrow(df_p3_all), "zones\n\n")
cat("--- Descriptive statistics: Period 1 ---\n")
print(summary(df_p1_all[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 2 ---\n")
print(summary(df_p2_all[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 3 ---\n")
print(summary(df_p3_all[, c("delta_lnL", "IPW", "IPW_OC")]))

# ============================================================
# 4b. BUILD DATASETS — MANUFACTURING ONLY (C-E)
# ============================================================

df_p1_mfg <- build_regression_data(t0 = 2000, t1 = 2007, industry_filter = MFG_INDUSTRIES)
df_p2_mfg <- build_regression_data(t0 = 2008, t1 = 2011, industry_filter = MFG_INDUSTRIES)
df_p3_mfg <- build_regression_data(t0 = 2012, t1 = 2019, industry_filter = MFG_INDUSTRIES)

cat("\n=== MANUFACTURING ONLY (C-E) ===\n")
cat(COUNTRY, "- Period 1 (2000-2007): N =", nrow(df_p1_mfg), "zones\n")
cat(COUNTRY, "- Period 2 (2008-2011): N =", nrow(df_p2_mfg), "zones\n")
cat(COUNTRY, "- Period 3 (2012-2019): N =", nrow(df_p3_mfg), "zones\n\n")
cat("--- Descriptive statistics: Period 1 ---\n")
print(summary(df_p1_mfg[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 2 ---\n")
print(summary(df_p2_mfg[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 3 ---\n")
print(summary(df_p3_mfg[, c("delta_lnL", "IPW", "IPW_OC")]))

# ============================================================
# 5. REGRESSION
# ============================================================

run_china_shock <- function(df, period_label, run_iv = TRUE) {

  sep <- paste0(rep("=", 64), collapse = "")
  cat("\n", sep, "\n", sep = "")
  cat("  PERIOD:", period_label, "  N =", nrow(df), "zones\n")
  cat(sep, "\n", sep = "")

  n_clusters <- nrow(df)
  if (n_clusters < 50)
    warning(sprintf(
      "Period '%s': N = %d clusters < 50. Clustered SE may be unreliable. ",
      period_label, n_clusters,
      "Consider wild bootstrap (package 'boottest') or HC2 as a robustness check."))

  ols          <- lm(delta_lnL ~ IPW, data = df)
  ols_vcov     <- vcovCL(ols, cluster = ~zone)
  ols_vcov_hc2 <- vcovHC(ols, type = "HC2")

  cat("\n[OLS]  delta_lnL ~ IPW  [Clustered SE]\n")
  print(coeftest(ols, vcov = ols_vcov))
  cat("[OLS]  Robustness: HC2 SE\n")
  print(coeftest(ols, vcov = ols_vcov_hc2))

  ols_se <- sqrt(diag(ols_vcov))
  ols_p  <- coeftest(ols, vcov = ols_vcov)[, "Pr(>|t|)"]

  cooks <- cooks.distance(ols)
  influential <- which(cooks > 4 / nrow(df))
  if (length(influential) > 0)
    cat(sprintf("\n  [!] Influential observations (Cook's D > 4/N): %s\n",
                paste(df$zone[influential], collapse = ", ")))

  if (run_iv) {
    fs      <- lm(IPW ~ IPW_OC, data = df)
    fs_vcov <- vcovCL(fs, cluster = ~zone)
    fs_wald <- waldtest(fs, vcov = fs_vcov)
    f_stat  <- fs_wald$F[2]
    cat("\n[First Stage]  IPW ~ IPW_OC\n")
    print(coeftest(fs, vcov = fs_vcov))
    cat(sprintf("  Clustered Wald F = %.3f  [Staiger-Stock: F > 10]\n", f_stat))

    if (f_stat < 10)
      warning(sprintf(
        "Weak instrument: Wald F = %.2f < 10, period '%s'.",
        f_stat, period_label))

    iv      <- ivreg(delta_lnL ~ IPW | IPW_OC, data = df)
    iv_vcov <- vcovCL(iv, cluster = ~zone)
    cat("\n[IV-2SLS]  delta_lnL ~ IPW | IPW_OC\n")
    print(coeftest(iv, vcov = iv_vcov))

    iv_se <- sqrt(diag(iv_vcov))
    iv_p  <- coeftest(iv, vcov = iv_vcov)[, "Pr(>|t|)"]
  } else {
    # IV not identified: with a single filtered industry, IPW and IPW_OC are collinear
    cat("\n[IV-2SLS] — skipped: instrument not identified with a single industry.\n")
    iv <- ols; iv_se <- ols_se; iv_p <- ols_p; f_stat <- NA
  }

  if (run_iv) {
    hausman_stat <- tryCatch({
      b_ols <- coef(ols)["IPW"]
      b_iv  <- coef(iv)["IPW"]
      var_iv  <- iv_vcov["IPW", "IPW"]
      var_ols <- ols_vcov["IPW", "IPW"]

      H <- (b_iv - b_ols)^2 / (var_iv - var_ols)
      p <- pchisq(H, df = 1, lower.tail = FALSE)

      cat(sprintf("\n[Hausman test]  H = %.4f  (chi2(1)),  p-value = %.4f\n", H, p))
      if (p < 0.05)
        cat("  -> Reject H0: IPW is endogenous; IV preferred over OLS\n")
      else
        cat("  -> Fail to reject H0: IPW is exogenous; OLS is consistent\n")

      list(stat = H, p = p)
    }, error = function(e) {
      cat("\n[Hausman test]  Could not compute:", conditionMessage(e), "\n")
      NULL
    })
  }

  invisible(list(
    ols = ols, iv = iv, first_stage = if (run_iv) fs else NULL,
    ols_se = ols_se, iv_se = iv_se,
    ols_p = ols_p, iv_p = iv_p,
    f_stat = f_stat, period = period_label,
    cooks = cooks, run_iv = run_iv,
    n_clusters = n_clusters
  ))
}

cat("\n\n========== REGRESSIONS: FULL ECONOMY ==========\n")
res_p1_all <- run_china_shock(df_p1_all, "2000-2007 | All sectors",    run_iv = TRUE)
res_p2_all <- run_china_shock(df_p2_all, "2008-2011 | All sectors",    run_iv = TRUE)
res_p3_all <- run_china_shock(df_p3_all, "2012-2019 | All sectors",    run_iv = TRUE)

cat("\n\n========== REGRESSIONS: MANUFACTURING (C-E) ==========\n")
res_p1_mfg <- run_china_shock(df_p1_mfg, "2000-2007 | Manufacturing", run_iv = FALSE)
res_p2_mfg <- run_china_shock(df_p2_mfg, "2008-2011 | Manufacturing", run_iv = FALSE)
res_p3_mfg <- run_china_shock(df_p3_mfg, "2012-2019 | Manufacturing", run_iv = FALSE)

# ============================================================
# 6. SUMMARY TABLES
# ============================================================

cat("\n\n")

stargazer(
  res_p1_all$ols, res_p1_all$iv,
  res_p2_all$ols, res_p2_all$iv,
  res_p3_all$ols, res_p3_all$iv,
  type             = "text",
  title            = paste("China Shock - Full Economy:", COUNTRY),
  column.labels    = c("2000-2007", "2008-2011", "2012-2019"),
  column.separate  = c(2L, 2L, 2L),
  model.names      = FALSE,
  model.numbers    = TRUE,
  dep.var.labels   = "Delta ln(Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_all$ols_se, res_p1_all$iv_se,
            res_p2_all$ols_se, res_p2_all$iv_se,
            res_p3_all$ols_se, res_p3_all$iv_se),
  p  = list(res_p1_all$ols_p,  res_p1_all$iv_p,
            res_p2_all$ols_p,  res_p2_all$iv_p,
            res_p3_all$ols_p,  res_p3_all$iv_p),
  add.lines = list(
    c("Estimator",    "OLS", "2SLS",      "OLS", "2SLS",      "OLS", "2SLS"),
    c("Sample",       "All", "All",        "All", "All",        "All", "All"),
    c("Instrument",   "---", "OC imports", "---", "OC imports", "---", "OC imports"),
    c("Clustered SE", "Yes", "Yes",        "Yes", "Yes",        "Yes", "Yes"),
    c("1st-stage F",
      "---", formatC(res_p1_all$f_stat, digits = 2, format = "f"),
      "---", formatC(res_p2_all$f_stat, digits = 2, format = "f"),
      "---", formatC(res_p3_all$f_stat, digits = 2, format = "f"))
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

cat("\n\n")

stargazer(
  res_p1_mfg$ols, res_p2_mfg$ols, res_p3_mfg$ols,
  type             = "text",
  title            = paste("China Shock - Manufacturing C-E (OLS only):", COUNTRY),
  column.labels    = c("2000-2007", "2008-2011", "2012-2019"),
  model.names      = FALSE,
  dep.var.labels   = "Delta ln(Mfg Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_mfg$ols_se, res_p2_mfg$ols_se, res_p3_mfg$ols_se),
  p  = list(res_p1_mfg$ols_p,  res_p2_mfg$ols_p,  res_p3_mfg$ols_p),
  add.lines = list(
    c("Estimator",    "OLS", "OLS", "OLS"),
    c("Sample",       "C-E", "C-E", "C-E"),
    c("Clustered SE", "Yes", "Yes", "Yes"),
    c("Note", "IV not identified", "IV not identified", "IV not identified")
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

cat("\n\n")

stargazer(
  res_p1_all$iv, res_p1_mfg$ols,
  res_p2_all$iv, res_p2_mfg$ols,
  res_p3_all$iv, res_p3_mfg$ols,
  type             = "text",
  title            = paste("China Shock - IV (All) vs OLS (Mfg):", COUNTRY),
  column.labels    = c("2000-2007", "2008-2011", "2012-2019"),
  column.separate  = c(2L, 2L, 2L),
  model.names      = FALSE,
  model.numbers    = TRUE,
  dep.var.labels   = "Delta ln(Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_all$iv_se, res_p1_mfg$ols_se,
            res_p2_all$iv_se, res_p2_mfg$ols_se,
            res_p3_all$iv_se, res_p3_mfg$ols_se),
  p  = list(res_p1_all$iv_p,  res_p1_mfg$ols_p,
            res_p2_all$iv_p,  res_p2_mfg$ols_p,
            res_p3_all$iv_p,  res_p3_mfg$ols_p),
  add.lines = list(
    c("Estimator",    "2SLS", "OLS", "2SLS", "OLS", "2SLS", "OLS"),
    c("Sample",       "All",  "C-E", "All",  "C-E", "All",  "C-E"),
    c("Clustered SE", "Yes",  "Yes", "Yes",  "Yes", "Yes",  "Yes"),
    c("1st-stage F",
      formatC(res_p1_all$f_stat, digits = 2, format = "f"), "n/a",
      formatC(res_p2_all$f_stat, digits = 2, format = "f"), "n/a",
      formatC(res_p3_all$f_stat, digits = 2, format = "f"), "n/a")
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

# ============================================================
# 7. REGIONAL IMPACT RANKING TABLE
# ============================================================

make_impact_table_iv <- function(df, res_iv, period_label) {
  beta_iv <- coef(res_iv$iv)["IPW"]

  cat("\n")
  cat(paste0(rep("=", 64), collapse = ""), "\n")
  cat(sprintf("  REGIONS RANKED BY EXPOSURE TO CHINA SHOCK\n"))
  cat(sprintf("  Period: %s | Estimate: 2SLS  (b_IV = %.4f)\n",
              period_label, beta_iv))
  cat(sprintf("  Sorted by IPW (descending)\n"))
  cat(paste0(rep("=", 64), collapse = ""), "\n\n")

  tbl <- df |>
    select(zone, IPW, delta_lnL) |>
    mutate(
      fitted_iv = beta_iv * IPW
    ) |>
    arrange(desc(IPW)) |>
    mutate(
      rank       = row_number(),
      IPW_fmt    = formatC(IPW,       digits = 4, format = "f"),
      fitted_fmt = formatC(fitted_iv, digits = 4, format = "f"),
      dlnL_fmt   = formatC(delta_lnL, digits = 4, format = "f"),
      severity   = case_when(
        IPW > quantile(IPW, 0.75) ~ "HIGH",
        IPW > quantile(IPW, 0.25) ~ "MEDIUM",
        TRUE                       ~ "LOW"
      )
    )

  header <- sprintf("%-4s  %-30s  %10s  %14s  %12s  %8s",
                    "Rank", "Region",
                    "IPW (th.USD)", "b*IPW (2SLS)", "DlnL", "Severity")
  sep <- paste0(rep("-", nchar(header)), collapse = "")

  cat(header, "\n", sep, "\n", sep = "")

  for (i in seq_len(nrow(tbl))) {
    cat(sprintf("%-4d  %-30s  %10s  %14s  %12s  %8s\n",
                tbl$rank[i],
                tbl$zone[i],
                tbl$IPW_fmt[i],
                tbl$fitted_fmt[i],
                tbl$dlnL_fmt[i],
                tbl$severity[i]))
  }

  cat(sep, "\n")
  cat(sprintf("  Most exposed:  %s  (b*IPW = %s)\n",
              tbl$zone[1], tbl$fitted_fmt[1]))
  cat(sprintf("  Least exposed: %s  (b*IPW = %s)\n",
              tbl$zone[nrow(tbl)], tbl$fitted_fmt[nrow(tbl)]))
  cat(sprintf("  HIGH   : %d region(s)\n", sum(tbl$severity == "HIGH")))
  cat(sprintf("  MEDIUM : %d region(s)\n", sum(tbl$severity == "MEDIUM")))
  cat(sprintf("  LOW    : %d region(s)\n", sum(tbl$severity == "LOW")))
  cat("\n")

  invisible(tbl)
}

impact_p1_all <- make_impact_table_iv(df_p1_all, res_p1_all, "2000-2007")
impact_p2_all <- make_impact_table_iv(df_p2_all, res_p2_all, "2008-2011")
impact_p3_all <- make_impact_table_iv(df_p3_all, res_p3_all, "2012-2019")

# ============================================================
# 8. PLOTS: scatter (DlnL ~ IPW) + Cook's Distance
# ============================================================

theme_china <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(size = 10, face = "bold",    hjust = 0),
    plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0),
    axis.line     = element_line(linewidth = 0.4),
    axis.ticks    = element_line(linewidth = 0.3),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9, color = "black"),
    legend.position = "none",
    plot.margin   = margin(6, 8, 4, 6, "pt")
  )

COL_POINT  <- "grey20"
COL_FIT    <- "black"
COL_LABEL  <- "grey30"
COL_THRESH <- "grey40"
COL_HIGH   <- "grey15"

make_scatter <- function(df, res, period_label, sample_label = "Full economy") {
  b  <- coef(res$ols)
  se <- res$ols_se["IPW"]
  r2 <- summary(res$ols)$r.squared

  fit <- data.frame(IPW = seq(min(df$IPW), max(df$IPW), length.out = 200))
  fit$delta_lnL <- b[1] + b[2] * fit$IPW

  note <- sprintf("b = %.4f (SE = %.4f)   R2 = %.3f   N = %d",
                  b[2], se, r2, nrow(df))

  thresh_cook <- 4 / nrow(df)
  ipw_cutoff  <- quantile(abs(df$IPW - median(df$IPW)), 0.80)

  df <- df |>
    mutate(
      cook_d    = res$cooks,
      fitted    = b[1] + b[2] * IPW,
      nudge_y   = sign(delta_lnL - fitted) * diff(range(delta_lnL)) * 0.04,
      label_txt = ifelse(
        cook_d > thresh_cook | abs(IPW - median(IPW)) > ipw_cutoff,
        zone, ""
      )
    )

  ggplot(df, aes(IPW, delta_lnL)) +
    geom_point(shape = 21, fill = "white", color = COL_POINT,
               size = 2.5, stroke = 0.6) +
    geom_line(data = fit, color = COL_FIT, linewidth = 0.8) +
    geom_text_repel(
      aes(label = label_txt),
      nudge_y            = df$nudge_y,
      size               = 2.6,
      color              = COL_LABEL,
      segment.size       = 0.25,
      segment.color      = "grey60",
      min.segment.length = 0.2,
      box.padding        = 0.5,
      point.padding      = 0.4,
      force              = 3,
      force_pull         = 0.5,
      max.overlaps       = Inf,
      seed               = 42
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             label = note, size = 2.8, color = "grey35", fontface = "italic") +
    labs(
      title    = paste0(period_label, " - ", sample_label),
      subtitle = "OLS fit with region labels (Cook's D > 4/N or IPW outlier)",
      x = "Import exposure (thousand USD per worker)",
      y = expression(Delta * " ln(Employment)")
    ) +
    theme_china
}

make_cooks <- function(df, res, period_label, sample_label = "Full economy") {
  thresh   <- 4 / nrow(df)
  cooks_df <- data.frame(zone  = df$zone,
                         cooks = res$cooks) |>
    mutate(
      influential = cooks > thresh,
      bar_label   = ifelse(cooks > thresh, sprintf("%.3f", cooks), "")
    )

  n_inf <- sum(cooks_df$influential)
  note  <- sprintf("Threshold 4/N = %.3f  |  Influential: %d region(s)",
                   thresh, n_inf)

  if (n_inf > 0)
    message("Cook's D > 4/N [", period_label, "]: ",
            paste(cooks_df$zone[cooks_df$influential], collapse = ", "))

  ggplot(cooks_df, aes(reorder(zone, cooks), cooks)) +
    geom_col(aes(fill = influential), width = 0.65) +
    geom_hline(yintercept = thresh, linetype = "dashed",
               color = COL_THRESH, linewidth = 0.5) +
    geom_text(aes(label = bar_label), hjust = -0.15,
              size = 2.6, color = "grey20") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = COL_HIGH),
                      guide  = "none") +
    annotate("text",
             x = Inf, y = thresh, hjust = 1.05, vjust = -0.5,
             label = sprintf("4/N = %.3f", thresh),
             color = COL_THRESH, size = 2.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             label = note, size = 2.8, color = "grey35", fontface = "italic") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip() +
    labs(
      title    = paste0(period_label, " - ", sample_label),
      subtitle = "Cook's distance per region; dark bars exceed threshold 4/N",
      x = NULL,
      y = "Cook's distance"
    ) +
    theme_china +
    theme(axis.text.y = element_text(size = 7.5))
}

p_scatter_p1 <- make_scatter(df_p1_all, res_p1_all, "2000-2007", "Full economy")
p_scatter_p2 <- make_scatter(df_p2_all, res_p2_all, "2008-2011", "Full economy")
p_scatter_p3 <- make_scatter(df_p3_all, res_p3_all, "2012-2019", "Full economy")
p_cooks_p1   <- make_cooks  (df_p1_all, res_p1_all, "2000-2007", "Full economy")
p_cooks_p2   <- make_cooks  (df_p2_all, res_p2_all, "2008-2011", "Full economy")
p_cooks_p3   <- make_cooks  (df_p3_all, res_p3_all, "2012-2019", "Full economy")

invisible(lapply(list(
  p_scatter_p1, p_scatter_p2, p_scatter_p3,
  p_cooks_p1,   p_cooks_p2,   p_cooks_p3
), print))

dir.create("figures", showWarnings = FALSE)

plots_to_save <- list(
  list(p_scatter_p1, "scatter_p1"),
  list(p_scatter_p2, "scatter_p2"),
  list(p_scatter_p3, "scatter_p3"),
  list(p_cooks_p1,   "cooks_p1"),
  list(p_cooks_p2,   "cooks_p2"),
  list(p_cooks_p3,   "cooks_p3")
)

for (item in plots_to_save) {
  ggsave(
    sprintf("figures/%s_%s.pdf", COUNTRY_CODE, item[[2]]),
    plot = item[[1]], width = 5.5, height = 3.8, units = "in",
    dpi = 300, bg = "white"
  )
}
message("Plots saved to figures/")

# ============================================================
# 9. EXPORT DATASETS
# ============================================================

dir.create("data_out", showWarnings = FALSE)

write.csv(df_p1_all, paste0("data_out/", COUNTRY_CODE, "_df_p1_all.csv"), row.names = FALSE)
write.csv(df_p2_all, paste0("data_out/", COUNTRY_CODE, "_df_p2_all.csv"), row.names = FALSE)
write.csv(df_p3_all, paste0("data_out/", COUNTRY_CODE, "_df_p3_all.csv"), row.names = FALSE)
write.csv(df_p1_mfg, paste0("data_out/", COUNTRY_CODE, "_df_p1_mfg.csv"), row.names = FALSE)
write.csv(df_p2_mfg, paste0("data_out/", COUNTRY_CODE, "_df_p2_mfg.csv"), row.names = FALSE)
write.csv(df_p3_mfg, paste0("data_out/", COUNTRY_CODE, "_df_p3_mfg.csv"), row.names = FALSE)

message("Datasets saved to data_out/")
