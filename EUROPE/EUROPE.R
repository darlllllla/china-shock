# ============================================================
#  China Shock Analysis — European NUTS-2 regions
#  Methodology: Autor, Dorn & Hanson (2013), adapted for Europe
#  Bartik (shift-share) IPW + IV-2SLS
#  Periods: 1999-2011 and 2012-2019
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

# ============================================================
# 1. НАСТРОЙКИ СТРАНЫ
# ============================================================

FILE         <- file.choose()
COUNTRY      <- "Spain"
COUNTRY_CODE <- "sp"                   # code

# ============================================================
# 2. ЗАГРУЗКА ДАННЫХ
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
  warning("Отрасли в employment, но НЕ в import: ", paste(missing_in_imp, collapse = ", "))
if (length(missing_in_emp) > 0)
  message("Отрасли в import, но НЕ в employment: ", paste(missing_in_emp, collapse = ", "))

import_oc_raw <- read_excel(FILE, sheet = "import OC",
                            col_names = FALSE, .name_repair = "minimal")

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
# 2б. СПИСОК ОТРАСЛЕЙ ДЛЯ MANUFACTURING-ФИЛЬТРА
# ============================================================

MFG_INDUSTRIES <- c(
  "c-e"
)

# ============================================================
# 3. ПОСТРОЕНИЕ ДАТАСЕТА
# ============================================================
# Параметр industry_filter:
#   NULL             → полная экономика (поведение по умолчанию)
#   character vector → только указанные отрасли (например, MFG_INDUSTRIES)

build_regression_data <- function(t0, t1, industry_filter = NULL) {
  
  # --- ВСЕГДА берём полные данные занятости для знаменателя L_r0 ---
  emp_t0_full <- filter(employment, year == t0)
  emp_t1_full <- filter(employment, year == t1)
  
  # Для числителей delta_lnL и s_ri — фильтруем если нужно
  emp_t0 <- if (!is.null(industry_filter)) filter(emp_t0_full, industry %in% industry_filter) else emp_t0_full
  emp_t1 <- if (!is.null(industry_filter)) filter(emp_t1_full, industry %in% industry_filter) else emp_t1_full
  
  imp_data    <- if (!is.null(industry_filter)) filter(import_de,    industry %in% industry_filter) else import_de
  imp_oc_data <- if (!is.null(industry_filter)) filter(import_oc_agg, industry %in% industry_filter) else import_oc_agg
  
  if (!is.null(industry_filter) && nrow(emp_t0) == 0)
    stop("industry_filter не совпал ни с одной отраслью. ",
         "Проверь: cat(sort(unique(employment$industry)), sep='\\n')")
  
  # --- delta_lnL: рост занятости только в отфильтрованных отраслях ---
  L_r0_filt <- emp_t0 |> group_by(zone) |>
    summarise(L_r0_filt = sum(emp, na.rm = TRUE), .groups = "drop")
  L_r1_filt <- emp_t1 |> group_by(zone) |>
    summarise(L_r1_filt = sum(emp, na.rm = TRUE), .groups = "drop")
  
  dep_var <- L_r0_filt |>
    inner_join(L_r1_filt, by = "zone") |>
    mutate(delta_lnL = log(L_r1_filt) - log(L_r0_filt))
  
  # диагностика потерь из inner_join
  lost_t1 <- setdiff(L_r0_filt$zone, L_r1_filt$zone)
  lost_t0 <- setdiff(L_r1_filt$zone, L_r0_filt$zone)
  if (length(lost_t1) > 0)
    warning("Зоны есть в t0, но нет в t1: ", paste(lost_t1, collapse = ", "))
  if (length(lost_t0) > 0)
    warning("Зоны есть в t1, но нет в t0: ", paste(lost_t0, collapse = ", "))
  
  # --- L_r0 для знаменателя s_ri: ВСЕГДА вся экономика ---
  L_r0_total <- emp_t0_full |> group_by(zone) |>
    summarise(L_r0 = sum(emp, na.rm = TRUE), .groups = "drop")
  
  # --- s_ri: доля отфильтрованных отраслей в СУММАРНОЙ занятости региона ---
  s_ri <- emp_t0 |>
    inner_join(L_r0_total, by = "zone") |>
    mutate(s_ri = emp / L_r0) |>
    select(zone, industry, s_ri)
  
  # --- L_i0: национальная занятость в отфильтрованных отраслях ---
  L_i0 <- emp_t0 |> group_by(industry) |>
    summarise(L_i0 = sum(emp, na.rm = TRUE), .groups = "drop")
  
  # --- Всё остальное без изменений ---
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
    warning(sprintf("Период %d-%d: %d отраслей получили NA shock_i", t0, t1, n_na_shock))
  
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
    warning(sprintf("Период %d-%d: исключены %d зон(ы): %s",
                    t0, t1, length(dropped), paste(dropped, collapse = ", ")))
  
  return(out)
}

# ============================================================
# 4. ФОРМИРОВАНИЕ ДАТАСЕТОВ — ПОЛНАЯ ЭКОНОМИКА
# ============================================================

df_p1_all <- build_regression_data(t0 = 1999, t1 = 2011)
df_p2_all <- build_regression_data(t0 = 2012, t1 = 2019)

cat("\n=== ПОЛНАЯ ЭКОНОМИКА ===\n")
cat(COUNTRY, "— Period 1 (1999-2011): N =", nrow(df_p1_all), "zones\n")
cat(COUNTRY, "— Period 2 (2012-2019): N =", nrow(df_p2_all), "zones\n\n")
cat("--- Descriptive statistics: Period 1 ---\n")
print(summary(df_p1_all[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 2 ---\n")
print(summary(df_p2_all[, c("delta_lnL", "IPW", "IPW_OC")]))

# ============================================================
# 4б. ФОРМИРОВАНИЕ ДАТАСЕТОВ — ТОЛЬКО MANUFACTURING (C-E)
# ============================================================

df_p1_mfg <- build_regression_data(t0 = 1999, t1 = 2011, MFG_INDUSTRIES)
df_p2_mfg <- build_regression_data(t0 = 2012, t1 = 2019, MFG_INDUSTRIES)

cat("\n=== ТОЛЬКО MANUFACTURING (C-E) ===\n")
cat(COUNTRY, "— Period 1 (1999-2011): N =", nrow(df_p1_mfg), "zones\n")
cat(COUNTRY, "— Period 2 (2012-2019): N =", nrow(df_p2_mfg), "zones\n\n")
cat("--- Descriptive statistics: Period 1 ---\n")
print(summary(df_p1_mfg[, c("delta_lnL", "IPW", "IPW_OC")]))
cat("\n--- Descriptive statistics: Period 2 ---\n")
print(summary(df_p2_mfg[, c("delta_lnL", "IPW", "IPW_OC")]))

# ============================================================
# 5. ПОСТРОЕНИЕ РЕГРЕССИИ
# ============================================================

run_china_shock <- function(df, period_label, run_iv = TRUE) {
  
  sep <- paste0(rep("=", 64), collapse = "")
  cat("\n", sep, "\n", sep = "")
  cat("  PERIOD:", period_label, "  N =", nrow(df), "zones\n")
  cat(sep, "\n", sep = "")
  
  n_clusters <- nrow(df)
  if (n_clusters < 50)
    warning(sprintf(
      "Период '%s': N = %d кластеров < 50. Clustered SE может быть ненадёжным. ",
      period_label, n_clusters,
      "Рассмотри wild bootstrap (пакет 'boottest') или HC2 как robustness check."))
  
  ols      <- lm(delta_lnL ~ IPW, data = df)
  ols_vcov     <- vcovCL(ols, cluster = ~zone)    # clustered (основной)
  ols_vcov_hc2 <- vcovHC(ols, type = "HC2")       # HC2 robustness check
  
  cat("\n[OLS]  delta_lnL ~ IPW  [Clustered SE]\n")
  print(coeftest(ols, vcov = ols_vcov))
  cat("[OLS]  Robustness: HC2 SE\n")
  print(coeftest(ols, vcov = ols_vcov_hc2))
  
  ols_se <- sqrt(diag(ols_vcov))
  ols_p  <- coeftest(ols, vcov = ols_vcov)[, "Pr(>|t|)"]
  
  cooks <- cooks.distance(ols)
  influential <- which(cooks > 4 / nrow(df))
  if (length(influential) > 0)
    cat(sprintf("\n  [!] Влиятельные наблюдения (Cook's D > 4/N): %s\n",
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
        "Слабый инструмент: Wald F = %.2f < 10, период '%s'.",
        f_stat, period_label))
    
    iv      <- ivreg(delta_lnL ~ IPW | IPW_OC, data = df)
    iv_vcov <- vcovCL(iv, cluster = ~zone)
    cat("\n[IV-2SLS]  delta_lnL ~ IPW | IPW_OC\n")
    print(coeftest(iv, vcov = iv_vcov))
    
    iv_se <- sqrt(diag(iv_vcov))
    iv_p  <- coeftest(iv, vcov = iv_vcov)[, "Pr(>|t|)"]
  } else {
    cat("\n[IV-2SLS] — пропущен: при одной отрасли инструмент не идентифицирован.\n")
    iv <- ols; iv_se <- ols_se; iv_p <- ols_p; f_stat <- NA
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

# --- Полная экономика: OLS + IV (как обычно) ---
cat("\n\n========== РЕГРЕССИИ: ПОЛНАЯ ЭКОНОМИКА ==========\n")
res_p1_all <- run_china_shock(df_p1_all, "1999-2011 | All sectors",    run_iv = TRUE)
res_p2_all <- run_china_shock(df_p2_all, "2012-2019 | All sectors",    run_iv = TRUE)

# --- Manufacturing: только OLS (IV не идентифицирован при 1 отрасли) ---
cat("\n\n========== РЕГРЕССИИ: MANUFACTURING (C-E) ==========\n")
cat("ПРИМЕЧАНИЕ: IV пропущен — при одной отрасли IPW и IPW_OC коллинеарны.\n")
res_p1_mfg <- run_china_shock(df_p1_mfg, "1999-2011 | Manufacturing", run_iv = FALSE)
res_p2_mfg <- run_china_shock(df_p2_mfg, "2012-2019 | Manufacturing", run_iv = FALSE)

# ============================================================
# 6. СВОДНЫЕ ТАБЛИЦЫ
# ============================================================

cat("\n\n")

# --- Таблица 1: Полная экономика (как раньше) ---
stargazer(
  res_p1_all$ols, res_p1_all$iv,
  res_p2_all$ols, res_p2_all$iv,
  type             = "text",
  title            = paste("China Shock — Full Economy:", COUNTRY),
  column.labels    = c("1999-2011", "2012-2019"),
  column.separate  = c(2L, 2L),
  model.names      = FALSE,
  model.numbers    = TRUE,
  dep.var.labels   = "Delta ln(Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_all$ols_se, res_p1_all$iv_se,
            res_p2_all$ols_se, res_p2_all$iv_se),
  p  = list(res_p1_all$ols_p,  res_p1_all$iv_p,
            res_p2_all$ols_p,  res_p2_all$iv_p),
  add.lines = list(
    c("Estimator",    "OLS", "2SLS",       "OLS", "2SLS"),
    c("Sample",       "All", "All",         "All", "All"),
    c("Instrument",   "---", "OC imports",  "---", "OC imports"),
    c("Clustered SE", "Yes", "Yes",         "Yes", "Yes"),
    c("1st-stage F",
      "---", formatC(res_p1_all$f_stat, digits = 2, format = "f"),
      "---", formatC(res_p2_all$f_stat, digits = 2, format = "f"))
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

cat("\n\n")

# --- Таблица 2: Manufacturing — только OLS ---
stargazer(
  res_p1_mfg$ols, res_p2_mfg$ols,
  type             = "text",
  title            = paste("China Shock — Manufacturing C-E (OLS only):", COUNTRY),
  column.labels    = c("1999-2011", "2012-2019"),
  model.names      = FALSE,
  dep.var.labels   = "Delta ln(Mfg Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_mfg$ols_se, res_p2_mfg$ols_se),
  p  = list(res_p1_mfg$ols_p,  res_p2_mfg$ols_p),
  add.lines = list(
    c("Estimator",    "OLS",  "OLS"),
    c("Sample",       "C-E",  "C-E"),
    c("Clustered SE", "Yes",  "Yes"),
    c("Note", "IV not identified (single industry)", "IV not identified (single industry)")
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

cat("\n\n")

# --- Таблица 3: Сравнение All vs Manufacturing ---
stargazer(
  res_p1_all$iv, res_p1_mfg$ols,   # ← mfg: явно берём $ols, не $iv
  res_p2_all$iv, res_p2_mfg$ols,
  type             = "text",
  title            = paste("China Shock — IV (All) vs OLS (Mfg):", COUNTRY),
  column.labels    = c("1999-2011", "2012-2019"),
  column.separate  = c(2L, 2L),
  model.names      = FALSE,
  model.numbers    = TRUE,
  dep.var.labels   = "Delta ln(Employment)",
  dep.var.caption  = "",
  covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
  se = list(res_p1_all$iv_se, res_p1_mfg$ols_se,
            res_p2_all$iv_se, res_p2_mfg$ols_se),
  p  = list(res_p1_all$iv_p,  res_p1_mfg$ols_p,
            res_p2_all$iv_p,  res_p2_mfg$ols_p),
  add.lines = list(
    c("Estimator",    "2SLS",  "OLS",     "2SLS",  "OLS"),
    c("Sample",       "All",   "C-E",     "All",   "C-E"),
    c("Clustered SE", "Yes",   "Yes",     "Yes",   "Yes"),
    c("1st-stage F",
      formatC(res_p1_all$f_stat, digits = 2, format = "f"), "n/a",
      formatC(res_p2_all$f_stat, digits = 2, format = "f"), "n/a")
  ),
  omit.stat = c("f", "ser"),
  no.space  = TRUE
)

# ============================================================
# ГРАФИКИ: scatter (Δln L ~ IPW) + Cook's Distance
# ============================================================

theme_china <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(size = 10, face = "bold",   hjust = 0),
    plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0),
    axis.line     = element_line(linewidth = 0.4),
    axis.ticks    = element_line(linewidth = 0.3),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9, color = "black"),
    legend.position = "none",
    plot.margin   = margin(6, 8, 4, 6, "pt")
  )

COL_POINT <- "grey20"
COL_FIT   <- "black"
COL_LABEL <- "grey30"
COL_THRESH <- "grey40"
COL_HIGH  <- "grey15"

# ---- Scatter -----------------------------------------------
make_scatter <- function(df, res, period_label, sample_label = "Full economy") {
  b  <- coef(res$ols)
  se <- res$ols_se["IPW"]
  r2 <- summary(res$ols)$r.squared
  
  fit <- data.frame(IPW = seq(min(df$IPW), max(df$IPW), length.out = 200))
  fit$delta_lnL <- b[1] + b[2] * fit$IPW
  
  note <- sprintf("\u03b2\u2009=\u2009%.4f (SE\u2009=\u2009%.4f)   R\u00b2\u2009=\u2009%.3f   N\u2009=\u2009%d",
                  b[2], se, r2, nrow(df))
  
  thresh_cook <- 4 / nrow(df)
  ipw_cutoff  <- quantile(abs(df$IPW - median(df$IPW)), 0.80)
  
  df <- df |>
    mutate(
      cook_d    = res$cooks,
      fitted    = b[1] + b[2] * IPW,
      # Отталкиваем метку от линии регрессии по знаку остатка
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
      nudge_y            = df$nudge_y,   # начальный сдвиг от линии
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
      title    = paste0(period_label, " \u2014 ", sample_label),
      subtitle = "OLS fit with region labels (Cook\u2019s D > 4/N or IPW outlier)",
      x = "Import exposure (thousand USD per worker)",
      y = expression(Delta * " ln(Employment)")
    ) +
    theme_china
}

# ---- Cook's Distance ---------------------------------------
make_cooks <- function(df, res, period_label, sample_label = "Full economy") {
  thresh   <- 4 / nrow(df)
  cooks_df <- data.frame(zone  = df$zone,
                         cooks = res$cooks) |>
    mutate(
      influential = cooks > thresh,
      bar_label   = ifelse(cooks > thresh, sprintf("%.3f", cooks), "")
    )
  
  n_inf <- sum(cooks_df$influential)
  note  <- sprintf("Threshold 4/N\u2009=\u2009%.3f  |  Influential: %d region(s)",
                   thresh, n_inf)
  
  if (n_inf > 0)
    message("Cook's D > 4/N [", period_label, "]: ",
            paste(cooks_df$zone[cooks_df$influential], collapse = ", "))
  
  ggplot(cooks_df, aes(reorder(zone, cooks), cooks)) +
    geom_col(aes(fill = influential), width = 0.65) +
    geom_hline(yintercept = thresh, linetype = "dashed",
               color = COL_THRESH, linewidth = 0.5) +
    # Значение на столбце — только для влиятельных
    geom_text(aes(label = bar_label), hjust = -0.15,
              size = 2.6, color = "grey20") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = COL_HIGH),
                      guide  = "none") +
    # vjust = -0.4 поднимает подпись над пунктиром, не уезжает за ось
    annotate("text",
             x = Inf, y = thresh, hjust = 1.05, vjust = -0.5,
             label = sprintf("4/N\u2009=\u2009%.3f", thresh),
             color = COL_THRESH, size = 2.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             label = note, size = 2.8, color = "grey35", fontface = "italic") +
    # Запас справа чтобы подписи Cook's D на столбцах не обрезались
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip() +
    labs(
      title    = paste0(period_label, " \u2014 ", sample_label),
      subtitle = "Cook\u2019s distance per region; dark bars exceed threshold 4/N",
      x = NULL,
      y = "Cook's distance"
    ) +
    theme_china +
    theme(axis.text.y = element_text(size = 7.5))
}

# ---- Строим — только полная экономика ----------------------
p_scatter_p1 <- make_scatter(df_p1_all, res_p1_all, "1999–2011", "Full economy")
p_scatter_p2 <- make_scatter(df_p2_all, res_p2_all, "2012–2019", "Full economy")
p_cooks_p1   <- make_cooks  (df_p1_all, res_p1_all, "1999–2011", "Full economy")
p_cooks_p2   <- make_cooks  (df_p2_all, res_p2_all, "2012–2019", "Full economy")

# ---- Печатаем ----------------------------------------------
invisible(lapply(list(
  p_scatter_p1, p_scatter_p2,
  p_cooks_p1,   p_cooks_p2
), print))

# ---- Экспорт в PDF -----------------------------------------
dir.create("figures", showWarnings = FALSE)

plots_to_save <- list(
  list(p_scatter_p1, "scatter_p1"),
  list(p_scatter_p2, "scatter_p2"),
  list(p_cooks_p1,   "cooks_p1"),
  list(p_cooks_p2,   "cooks_p2")
)

for (item in plots_to_save) {
  ggsave(
    sprintf("figures/%s_%s.pdf", COUNTRY_CODE, item[[2]]),
    plot = item[[1]], width = 5.5, height = 3.8, units = "in",
    dpi = 300, bg = "white"
  )
}
message("Графики сохранены в figures/")


# ============================================================
# 10. ЭКСПОРТ В LATEX
# ============================================================

format_f <- function(x) {
  if (is.null(x) || !is.numeric(x) || is.na(x)) return(" ")
  formatC(x, digits = 2, format = "f")
}

# ---------------- FULL ECONOMY ----------------
stargazer(
  res_p1_all$ols, res_p1_all$iv,
  res_p2_all$ols, res_p2_all$iv,
  type = "latex",
  out  = paste0(COUNTRY_CODE, "_china_shock_all.tex"),
  
  title = paste("China Shock and Employment:", COUNTRY),
  label = paste0("tab:", COUNTRY_CODE, "_china_all"),
  
  column.labels   = c("1999--2011", "2012--2019"),
  column.separate = c(2, 2),
  
  model.names   = FALSE,
  model.numbers = TRUE,
  
  dep.var.labels  = "$\\Delta \\ln(\\mathrm{Employment})$",
  dep.var.caption = "",
  
  covariate.labels = c("Import exposure (IPW)", "Constant"),
  
  se = list(res_p1_all$ols_se, res_p1_all$iv_se,
            res_p2_all$ols_se, res_p2_all$iv_se),
  
  p  = list(res_p1_all$ols_p,  res_p1_all$iv_p,
            res_p2_all$ols_p,  res_p2_all$iv_p),
  
  add.lines = list(
    c("Estimator",    "OLS", "2SLS", "OLS", "2SLS"),
    c("Instrument",   "", "OC imports", "", "OC imports"),
    c("Clustered SE", "Yes", "Yes", "Yes", "Yes"),
    c("First-stage F",
      "", format_f(res_p1_all$f_stat),
      "", format_f(res_p2_all$f_stat))
  ),
  
  omit.stat = c("f", "ser"),
  no.space  = TRUE,
  float     = TRUE
)

cat("LaTeX (Full economy) сохранён\n")


# ---------------- MANUFACTURING ----------------
stargazer(
  res_p1_mfg$ols, res_p2_mfg$ols,
  type = "latex",
  out  = paste0(COUNTRY_CODE, "_china_shock_mfg.tex"),
  
  title = paste("China Shock and Employment (Manufacturing C--E):", COUNTRY),
  label = paste0("tab:", COUNTRY_CODE, "_china_mfg"),
  
  column.labels = c("1999--2011", "2012--2019"),
  
  model.names = FALSE,
  
  dep.var.labels = "$\\Delta \\ln(\\mathrm{Manufacturing\\ Employment})$",
  
  covariate.labels = c("Import exposure (IPW)", "Constant"),
  
  se = list(res_p1_mfg$ols_se, res_p2_mfg$ols_se),
  p  = list(res_p1_mfg$ols_p,  res_p2_mfg$ols_p),
  
  add.lines = list(
    c("Estimator", "OLS", "OLS"),
    c("Sample", "Manufacturing", "Manufacturing"),
    c("Clustered SE", "Yes", "Yes"),
    c("Note", "IV not identified", "IV not identified")
  ),
  
  omit.stat = c("f", "ser"),
  no.space  = TRUE,
  float     = TRUE
)

cat("LaTeX (Manufacturing) сохранён\n")


# ---------------- COMPARISON TABLE ----------------
stargazer(
  res_p1_all$iv, res_p1_mfg$ols,
  res_p2_all$iv, res_p2_mfg$ols,
  type = "latex",
  out  = paste0(COUNTRY_CODE, "_china_shock_compare.tex"),
  
  title = paste("China Shock: IV (Full) vs OLS (Manufacturing) —", COUNTRY),
  label = paste0("tab:", COUNTRY_CODE, "_china_compare"),
  
  column.labels   = c("1999--2011", "2012--2019"),
  column.separate = c(2, 2),
  
  model.names   = FALSE,
  model.numbers = TRUE,
  
  dep.var.labels = "$\\Delta \\ln(\\mathrm{Employment})$",
  
  covariate.labels = c("Import exposure (IPW)", "Constant"),
  
  se = list(res_p1_all$iv_se, res_p1_mfg$ols_se,
            res_p2_all$iv_se, res_p2_mfg$ols_se),
  
  p  = list(res_p1_all$iv_p,  res_p1_mfg$ols_p,
            res_p2_all$iv_p,  res_p2_mfg$ols_p),
  
  add.lines = list(
    c("Estimator", "2SLS", "OLS", "2SLS", "OLS"),
    c("Sample", "Full", "Mfg", "Full", "Mfg"),
    c("Clustered SE", "Yes", "Yes", "Yes", "Yes"),
    c("First-stage F",
      format_f(res_p1_all$f_stat), "",
      format_f(res_p2_all$f_stat), "")
  ),
  
  omit.stat = c("f", "ser"),
  no.space  = TRUE,
  float     = TRUE
)

cat("LaTeX (Comparison) сохранён\n")