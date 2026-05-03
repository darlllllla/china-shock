# ============================================================
#  China Shock Analysis — US Commuting Zones
#  Methodology: Autor, Dorn & Hanson (2013)
#  Bartik (shift-share) IPW + IV-2SLS
#  Periods: 1999-2011, 2012-2019
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(AER)
library(lmtest)
library(sandwich)
library(stargazer)
library(ggplot2)
library(ggrepel)

# ============================================================
# 1. CONFIGURATION
# ============================================================

PERIODS <- list(

  list(
    label            = "1999-2011",
    imp              = "./data/processed/delta_import_USA_1999_2011.csv",
    imp_iv           = "./data/processed/delta_import_IV_by_country_1999_2011.csv",
    cz_ind_mfg_t0    = "./data/raw/period_1999_2011/cz_industry_manuf_1999.csv",
    cz_ind_nonmfg_t0 = "./data/raw/period_1999_2011/cz_industry_non_manuf_1999.csv",
    cz_tot_t0        = "./data/raw/period_1999_2011/cz_total_manuf_1999.csv",
    cz_tot_t1        = "./data/raw/period_1999_2011/cz_total_manuf_2011.csv"
  ),

  list(
    label            = "2012-2019",
    imp              = "./data/processed/delta_import_USA_2012_2019.csv",
    imp_iv           = "./data/processed/delta_import_IV_by_country_2012_2019.csv",
    cz_ind_mfg_t0    = "./data/raw/period_2012_2019/cz_industry_manuf_2012.csv",
    cz_ind_nonmfg_t0 = "./data/raw/period_2012_2019/cz_industry_non_manuf_2012.csv",
    cz_tot_t0        = "./data/raw/period_2012_2019/cz_total_manuf_2012.csv",
    cz_tot_t1        = "./data/raw/period_2012_2019/cz_total_manuf_2019.csv"
  )

)

COUNTRY      <- "USA"
COUNTRY_CODE <- "usa"

# ============================================================
# 2. DATA LOADERS
# ============================================================

load_import <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  names(df) <- str_trim(names(df))
  naics_col <- names(df)[str_detect(tolower(names(df)), "naics")][1]
  dcol      <- names(df)[str_detect(tolower(names(df)), "import|delta")][1]
  if (is.na(naics_col)) stop("NAICS column not found in: ", path)
  if (is.na(dcol))      stop("d_import column not found in: ", path)
  df |>
    rename(NAICS4 = !!naics_col, d_import = !!dcol) |>
    mutate(
      NAICS4   = str_pad(str_sub(as.character(
        suppressWarnings(as.integer(as.numeric(NAICS4)))), 1, 4), 4, "left", "0"),
      d_import = suppressWarnings(as.numeric(d_import))
    ) |>
    filter(!is.na(NAICS4), !is.na(d_import), str_detect(NAICS4, "^\\d{4}$"))
}

load_import_iv <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  names(df) <- str_trim(names(df))
  naics_col   <- names(df)[str_detect(tolower(names(df)), "naics")][1]
  dimp_col    <- names(df)[str_detect(tolower(names(df)), "import|delta")][1]
  country_col <- names(df)[str_detect(tolower(names(df)), "country|reporter|nation")][1]
  if (is.na(naics_col)) stop("NAICS column not found in IV file: ", path)
  if (is.na(dimp_col))  stop("d_import column not found in IV file: ", path)
  df <- df |>
    rename(NAICS4 = !!naics_col, d_import_raw = !!dimp_col) |>
    mutate(
      NAICS4       = str_pad(str_sub(as.character(
        suppressWarnings(as.integer(as.numeric(NAICS4)))), 1, 4), 4, "left", "0"),
      d_import_raw = suppressWarnings(as.numeric(d_import_raw))
    ) |>
    filter(!is.na(NAICS4), !is.na(d_import_raw), str_detect(NAICS4, "^\\d{4}$"))
  if (!is.na(country_col))
    cat("    IV countries:", paste(unique(df[[country_col]]), collapse = ", "), "\n")
  df |>
    group_by(NAICS4) |>
    summarise(d_import_iv = sum(d_import_raw, na.rm = TRUE), .groups = "drop")
}

load_cz_industry <- function(path) {
  df <- read_delim(path, delim = ifelse(grepl(";", readLines(path, n = 1)), ";", ","),
                   show_col_types = FALSE)
  names(df) <- str_trim(names(df))
  cz_col    <- names(df)[str_detect(tolower(names(df)), "^cz$|^czone$")][1]
  naics_col <- names(df)[str_detect(tolower(names(df)), "naics")][1]
  emp_col   <- names(df)[str_detect(tolower(names(df)), "^emp$|employ")][1]
  if (is.na(cz_col))    stop("CZ column not found in: ", path)
  if (is.na(naics_col)) stop("NAICS column not found in: ", path)
  if (is.na(emp_col))   stop("emp column not found in: ", path)
  df |>
    rename(cz = !!cz_col, NAICS4 = !!naics_col, emp = !!emp_col) |>
    mutate(
      cz     = as.character(as.integer(as.numeric(cz))),
      NAICS4 = str_pad(str_sub(as.character(as.integer(as.numeric(NAICS4))), 1, 4), 4, "left", "0"),
      emp    = suppressWarnings(as.numeric(emp))
    ) |>
    filter(!is.na(cz), !is.na(NAICS4), !is.na(emp))
}

load_cz_total <- function(path) {
  df <- read_delim(path, delim = ifelse(grepl(";", readLines(path, n = 1)), ";", ","),
                   show_col_types = FALSE)
  names(df) <- str_trim(names(df))
  cz_col  <- names(df)[str_detect(tolower(names(df)), "^cz$|^czone$")][1]
  emp_col <- names(df)[str_detect(tolower(names(df)), "manuf|emp")][1]
  if (is.na(cz_col))  stop("CZ column not found in: ", path)
  if (is.na(emp_col)) stop("Employment column not found in: ", path)
  df |>
    rename(cz = !!cz_col, L_manuf = !!emp_col) |>
    mutate(
      cz      = as.character(as.integer(as.numeric(cz))),
      L_manuf = suppressWarnings(as.numeric(L_manuf))
    ) |>
    group_by(cz) |>
    summarise(L_manuf = first(L_manuf), .groups = "drop") |>
    filter(!is.na(cz), !is.na(L_manuf))
}

# ============================================================
# 3. BUILD IPW PANEL
# ============================================================
# sample = "manuf": manufacturing CZs only
# sample = "full":  manuf + non-manuf
# s_ri denominator: always manufacturing employment

build_ipw_panel <- function(cfg, imp_data, imp_col = "d_import",
                            sample = c("manuf", "full"), verbose = TRUE) {
  sample <- match.arg(sample)
  label  <- cfg$label

  if (verbose) {
    sep <- paste0(rep("=", 55), collapse = "")
    cat("\n", sep, "\n  Period:", label, "| Sample:", sample, "\n", sep, "\n", sep = "")
  }

  cz_ind_mfg <- load_cz_industry(cfg$cz_ind_mfg_t0)
  cz_ind <- if (sample == "full")
    bind_rows(cz_ind_mfg, load_cz_industry(cfg$cz_ind_nonmfg_t0))
  else
    cz_ind_mfg

  cz_tot_0 <- load_cz_total(cfg$cz_tot_t0)
  cz_tot_1 <- load_cz_total(cfg$cz_tot_t1)

  imp <- imp_data |> rename(d_import = !!imp_col)

  L_i0 <- cz_ind |>
    group_by(NAICS4) |>
    summarise(L_i0 = sum(emp, na.rm = TRUE), .groups = "drop")

  shock <- imp |>
    inner_join(L_i0, by = "NAICS4") |>
    filter(L_i0 > 0) |>
    mutate(shock_i = d_import / L_i0)
  if (verbose) cat("  Industries in shock:", nrow(shock), "\n")

  cz_shares <- cz_ind |>
    left_join(cz_tot_0, by = "cz") |>
    mutate(s_ri = emp / L_manuf) |>
    filter(!is.na(s_ri), !is.na(L_manuf))

  if (verbose) {
    share_sum <- cz_shares |> group_by(cz) |> summarise(s = sum(s_ri), .groups = "drop")
    cat("  Median share sum by CZ:", round(median(share_sum$s), 4), "\n")
  }

  IPW <- cz_shares |>
    inner_join(select(shock, NAICS4, shock_i), by = "NAICS4") |>
    mutate(component = s_ri * shock_i) |>
    group_by(cz) |>
    summarise(IPW = sum(component, na.rm = TRUE), .groups = "drop")

  panel <- cz_tot_0 |>
    inner_join(cz_tot_1, by = "cz", suffix = c("_t0", "_t1")) |>
    filter(L_manuf_t0 > 0, L_manuf_t1 > 0) |>
    mutate(d_ln_L = log(L_manuf_t1) - log(L_manuf_t0))

  panel <- panel |>
    left_join(IPW, by = "cz") |>
    mutate(IPW = replace_na(IPW, 0) / 1000, period = label, sample = sample)

  if (verbose) {
    cat("  CZ in panel:", nrow(panel), "\n")
    cat("  Median IPW:", round(median(panel$IPW), 2), "thou. USD / worker\n")
    cat("  Median Delta ln L:", round(median(panel$d_ln_L), 3), "\n")
  }
  panel
}

# ============================================================
# 4. RUN PIPELINE
# ============================================================

panels_mfg  <- list()
panels_full <- list()

for (cfg in PERIODS) {
  imp <- load_import(cfg$imp)
  lbl <- cfg$label
  panels_mfg[[lbl]]  <- build_ipw_panel(cfg, imp, "d_import", "manuf")
  panels_full[[lbl]] <- build_ipw_panel(cfg, imp, "d_import", "full")
}
cat("\nAll periods processed.\n")

# ============================================================
# 5. OLS + IV-2SLS
# ============================================================

run_china_shock <- function(panel, period_label) {
  sep <- paste0(rep("=", 64), collapse = "")
  cat("\n", sep, "\n", sep = "")
  cat("  PERIOD:", period_label, "  N =", nrow(panel), "CZ\n")
  cat(sep, "\n", sep = "")

  ols          <- lm(d_ln_L ~ IPW, data = panel)
  ols_vcov     <- vcovCL(ols, cluster = ~cz)
  ols_vcov_hc2 <- vcovHC(ols, type = "HC2")

  cat("\n[OLS]  d_ln_L ~ IPW  [Clustered SE]\n")
  print(coeftest(ols, vcov = ols_vcov))
  cat("[OLS]  Robustness: HC2 SE\n")
  print(coeftest(ols, vcov = ols_vcov_hc2))

  ols_se <- sqrt(diag(ols_vcov))
  ols_p  <- coeftest(ols, vcov = ols_vcov)[, "Pr(>|t|)"]
  cooks  <- cooks.distance(ols)
  n_inf  <- sum(cooks > 4 / nrow(panel))
  if (n_inf > 0) cat(sprintf("  [!] Influential CZ (Cook's D > 4/N): %d\n", n_inf))

  invisible(list(
    ols = ols, ols_vcov = ols_vcov, ols_se = ols_se, ols_p = ols_p,
    cooks = cooks, period = period_label, n_clusters = nrow(panel)
  ))
}

run_iv <- function(res_ols, panel_2sls, period_label) {
  fs      <- lm(IPW ~ IPW_iv, data = panel_2sls)
  fs_vcov <- vcovCL(fs, cluster = ~cz)
  f_stat  <- waldtest(fs, vcov = fs_vcov)$F[2]

  cat("\n[First Stage]  IPW ~ IPW_iv\n")
  print(coeftest(fs, vcov = fs_vcov))
  cat(sprintf("  Clustered Wald F = %.3f  [Staiger-Stock: F > 10]\n", f_stat))
  if (f_stat < 10)
    warning(sprintf("Weak instrument: Wald F = %.2f < 10, '%s'.", f_stat, period_label))

  iv      <- ivreg(d_ln_L ~ IPW | IPW_iv, data = panel_2sls)
  iv_vcov <- vcovCL(iv, cluster = ~cz)
  cat("\n[IV-2SLS]  d_ln_L ~ IPW | IPW_iv\n")
  print(coeftest(iv, vcov = iv_vcov))

  iv_se <- sqrt(diag(iv_vcov))
  iv_p  <- coeftest(iv, vcov = iv_vcov)[, "Pr(>|t|)"]

  c(res_ols, list(iv = iv, iv_vcov = iv_vcov, iv_se = iv_se, iv_p = iv_p,
                  f_stat = f_stat, first_stage = fs))
}

run_all <- function(panels_list, smp) {
  sample_label <- if (smp == "manuf") "Manufacturing" else "Full economy"
  cat("\n\n========== REGRESSIONS:", sample_label, "==========\n")
  results <- list()

  cat("\n--- OLS ---\n")
  for (lbl in names(panels_list))
    results[[lbl]] <- run_china_shock(panels_list[[lbl]], paste(lbl, "|", sample_label))

  cat("\n--- IV-2SLS ---\n")
  for (cfg in PERIODS) {
    lbl <- cfg$label
    cat("\n--- IV:", lbl, "---\n")
    imp_iv   <- load_import_iv(cfg$imp_iv)
    panel_iv <- build_ipw_panel(cfg, imp_iv, "d_import_iv", smp, verbose = FALSE)
    panel_iv <- rename(panel_iv, IPW_iv = IPW)
    p2sls    <- inner_join(panels_list[[lbl]], select(panel_iv, cz, IPW_iv), by = "cz")
    cat("  CZ in 2SLS panel:", nrow(p2sls), "\n")
    results[[lbl]] <- run_iv(results[[lbl]], p2sls, paste(lbl, "|", sample_label))
  }
  results
}

results_mfg  <- run_all(panels_mfg,  "manuf")
results_full <- run_all(panels_full, "full")

# ============================================================
# 6. SUMMARY TABLES
# ============================================================

format_f <- function(x) {
  if (is.null(x) || !is.numeric(x) || is.na(x)) return(" ")
  formatC(x, digits = 2, format = "f")
}

make_table <- function(results, title_sfx, dep_lbl = "Delta ln(Mfg Employment)") {
  lbl1 <- names(results)[1]; lbl2 <- names(results)[2]
  r1 <- results[[lbl1]];     r2 <- results[[lbl2]]
  stargazer(
    r1$ols, r1$iv, r2$ols, r2$iv,
    type             = "text",
    title            = paste("China Shock —", title_sfx, ":", COUNTRY),
    column.labels    = c(lbl1, lbl2),
    column.separate  = c(2L, 2L),
    model.names      = FALSE,
    model.numbers    = TRUE,
    dep.var.labels   = dep_lbl,
    dep.var.caption  = "",
    covariate.labels = c("IPW (thou. USD / worker)", "Constant"),
    se = list(r1$ols_se, r1$iv_se, r2$ols_se, r2$iv_se),
    p  = list(r1$ols_p,  r1$iv_p,  r2$ols_p,  r2$iv_p),
    add.lines = list(
      c("Estimator",    "OLS", "2SLS",      "OLS", "2SLS"),
      c("Instrument",   "---", "OC imports", "---", "OC imports"),
      c("Clustered SE", "Yes", "Yes",        "Yes", "Yes"),
      c("1st-stage F",
        "---", format_f(r1$f_stat), "---", format_f(r2$f_stat))
    ),
    omit.stat = c("f", "ser"),
    no.space  = TRUE
  )
}

cat("\n\n")
make_table(results_mfg,  "Manufacturing only")
cat("\n\n")
make_table(results_full, "Full economy", dep_lbl = "Delta ln(Employment)")

# ============================================================
# 7. PLOTS
# ============================================================

theme_china <- theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(size = 10, face = "bold",    hjust = 0),
    plot.subtitle   = element_text(size = 9,  color = "grey40", hjust = 0),
    axis.line       = element_line(linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9, color = "black"),
    legend.position = "none",
    plot.margin     = margin(6, 8, 4, 6, "pt")
  )

make_scatter <- function(df, res, period_label, sample_label, n_labels = 10) {
  b  <- coef(res$ols); se <- res$ols_se["IPW"]; r2 <- summary(res$ols)$r.squared
  note <- sprintf("b = %.4f (SE = %.4f)   R2 = %.3f   N = %d", b[2], se, r2, nrow(df))
  fit <- data.frame(IPW = seq(min(df$IPW), max(df$IPW), length.out = 200))
  fit$d_ln_L <- b[1] + b[2] * fit$IPW

  df2 <- df |>
    mutate(
      cook_d   = res$cooks,
      fitted_y = b[1] + b[2] * IPW,
      nudge_y  = sign(d_ln_L - fitted_y) * diff(range(d_ln_L)) * 0.035
    ) |>
    arrange(desc(cook_d)) |>
    mutate(label_txt = ifelse(row_number() <= n_labels, cz, ""))

  ggplot(df2, aes(IPW, d_ln_L)) +
    geom_hex(bins = 45, aes(fill = after_stat(log(count + 1))), color = NA) +
    scale_fill_gradient(low = "grey92", high = "grey10") +
    geom_line(data = fit, color = "black", linewidth = 0.9) +
    geom_text_repel(
      data = filter(df2, label_txt != ""),
      aes(label = label_txt),
      nudge_y            = filter(df2, label_txt != "")$nudge_y,
      size = 2.4, color = "grey15",
      segment.size = 0.25, segment.color = "grey55",
      min.segment.length = 0.1, box.padding = 0.4, point.padding = 0.3,
      force = 4, max.overlaps = Inf, seed = 42
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4,
             label = note, size = 2.7, color = "grey35", fontface = "italic") +
    labs(
      title    = paste0(period_label, " — ", sample_label),
      subtitle = paste0("Hex-binned density; top-", n_labels, " CZ by Cook's D labelled"),
      x = "Import exposure (thousand USD per worker)",
      y = expression(Delta * " ln(Mfg Employment)")
    ) +
    theme_china
}

make_cooks_top <- function(df, res, period_label, sample_label, top_n = 30) {
  thresh      <- 4 / nrow(df)
  n_total_inf <- sum(res$cooks > thresh)
  cooks_df <- data.frame(cz = df$cz, cooks = res$cooks) |>
    arrange(desc(cooks)) |>
    slice_head(n = top_n) |>
    mutate(influential = cooks > thresh,
           cz_f = factor(cz, levels = rev(cz)))

  note <- sprintf("4/N = %.4f  |  Influential: %d / %d CZ  |  Showing top %d",
                  thresh, n_total_inf, nrow(df), top_n)

  ggplot(cooks_df, aes(y = cz_f, x = cooks)) +
    geom_segment(aes(yend = cz_f, x = 0, xend = cooks, color = influential), linewidth = 0.6) +
    geom_point(aes(color = influential), size = 2.2) +
    geom_vline(xintercept = thresh, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey68", "TRUE" = "grey10"), guide = "none") +
    annotate("text", y = 0.5, x = thresh, hjust = -0.1, vjust = 0,
             label = sprintf("4/N=%.4f", thresh), color = "grey40", size = 2.5, angle = 90) +
    annotate("text", y = Inf, x = -Inf, hjust = -0.03, vjust = 1.4,
             label = note, size = 2.6, color = "grey35", fontface = "italic") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = paste0(period_label, " — ", sample_label),
      subtitle = paste0("Cook's distance, top-", top_n, " CZ; dark = exceed 4/N"),
      x = "Cook's distance", y = NULL
    ) +
    theme_china +
    theme(axis.text.y = element_text(size = 7))
}

make_hist <- function(df, period_label, sample_label) {
  med_ipw <- median(df$IPW); p90 <- quantile(df$IPW, 0.90)
  note <- sprintf("Median = %.2f  |  P90 = %.2f  |  N = %d CZ", med_ipw, p90, nrow(df))
  ggplot(df, aes(IPW)) +
    geom_histogram(bins = 50, fill = "grey30", color = "white", linewidth = 0.25) +
    geom_vline(xintercept = med_ipw, color = "black",  linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = p90,     color = "grey50", linetype = "dotted", linewidth = 0.7) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4,
             label = note, size = 2.7, color = "grey35", fontface = "italic") +
    labs(
      title    = paste0(period_label, " — ", sample_label),
      subtitle = "IPW distribution; dashed = median, dotted = P90",
      x = "Import exposure (thousand USD per worker)", y = "Number of CZ"
    ) +
    theme_china
}

make_all_plots <- function(panels, results, sample_label) {
  lbl1 <- names(panels)[1]; lbl2 <- names(panels)[2]
  list(
    scatter_p1 = make_scatter  (panels[[lbl1]], results[[lbl1]], lbl1, sample_label),
    scatter_p2 = make_scatter  (panels[[lbl2]], results[[lbl2]], lbl2, sample_label),
    cooks_p1   = make_cooks_top(panels[[lbl1]], results[[lbl1]], lbl1, sample_label),
    cooks_p2   = make_cooks_top(panels[[lbl2]], results[[lbl2]], lbl2, sample_label),
    hist_p1    = make_hist     (panels[[lbl1]], lbl1, sample_label),
    hist_p2    = make_hist     (panels[[lbl2]], lbl2, sample_label)
  )
}

plots_mfg  <- make_all_plots(panels_mfg,  results_mfg,  "Manufacturing")
plots_full <- make_all_plots(panels_full, results_full, "Full economy")

invisible(lapply(c(plots_mfg, plots_full), print))

dir.create("figures", showWarnings = FALSE)
save_plots <- function(plots, prefix) {
  for (nm in names(plots))
    ggsave(sprintf("figures/%s_%s_%s.pdf", COUNTRY_CODE, prefix, nm),
           plot = plots[[nm]], width = 5.5, height = 3.8, units = "in", dpi = 300, bg = "white")
}
save_plots(plots_mfg,  "mfg")
save_plots(plots_full, "full")
message("Plots saved to figures/")

# ============================================================
# 8. SAVE DATASETS
# ============================================================

dir.create("data_out", showWarnings = FALSE)
for (lbl in names(panels_mfg)) {
  safe <- str_replace_all(lbl, "-", "_")
  write.csv(panels_mfg[[lbl]],
            paste0("data_out/", COUNTRY_CODE, "_panel_mfg_",  safe, ".csv"), row.names = FALSE)
  write.csv(panels_full[[lbl]],
            paste0("data_out/", COUNTRY_CODE, "_panel_full_", safe, ".csv"), row.names = FALSE)
}
message("Datasets saved to data_out/")
