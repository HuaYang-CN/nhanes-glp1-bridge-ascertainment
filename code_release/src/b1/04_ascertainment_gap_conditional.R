
#!/usr/bin/env Rscript

# Configuration and guarded writers supplied by src/runtime.R.
OUT_DIR <- file.path(RESULT_ROOT, "ascertainment_gap")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
dat <- readRDS(file.path(DERIVED_DIR, "bridge_analysis_ready.rds"))

make_design <- function(d) {
  svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTPH2YR, nest = TRUE, data = d)
}
fmt_p <- function(p) dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p))
prop_est <- function(design, var) {
  x <- svyciprop(as.formula(paste0("~", var)), design = design, method = "logit", df = degf(design), na.rm = TRUE)
  ci <- confint(x)
  m <- svymean(as.formula(paste0("~", var)), design = design, na.rm = TRUE)
  est0 <- as.numeric(coef(m)[1]); se0 <- as.numeric(SE(m)[1]); tval <- est0 / se0
  pval <- 2 * pt(abs(tval), df = degf(design), lower.tail = FALSE)
  tibble(
    estimate = 100 * as.numeric(coef(x)[1]),
    lower = 100 * as.numeric(ci[1]),
    upper = 100 * as.numeric(ci[2]),
    p_value = pval,
    ci_label = sprintf("%.1f%% (%.1f%%–%.1f%%)", 100 * as.numeric(coef(x)[1]), 100 * as.numeric(ci[1]), 100 * as.numeric(ci[2])),
    p_label = fmt_p(pval)
  )
}

d <- dat %>% filter(base_eligible, !is.na(primary_observable), !is.na(diagnosis_proxy))
des <- make_design(d)

des_measurement_pos <- subset(des, primary_observable == 1)
p_key <- prop_est(des_measurement_pos, "measurement_only")
key_tab <- tibble(
  characteristic = "Measurement-only among measurement-positive adults",
  estimate = p_key$estimate,
  lower = p_key$lower,
  upper = p_key$upper,
  p_value = p_key$p_value,
  ci_label = p_key$ci_label,
  p_label = p_key$p_label
)
write_csv(key_tab, file.path(OUT_DIR, "KEY_measurement_only_among_measurement_positive.csv"))
