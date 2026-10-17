
#!/usr/bin/env Rscript

# Configuration and guarded writers supplied by src/runtime.R.
OUT_DIR <- file.path(RESULT_ROOT, "sensitivity_fpg")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
dat <- readRDS(file.path(DERIVED_DIR, "bridge_analysis_ready.rds"))

bmi_plus_binary <- function(bmi, threshold, binary01) {
  out <- rep(NA_integer_, length(bmi))
  out[!is.na(bmi) & bmi < threshold] <- 0L
  idx <- !is.na(bmi) & bmi >= threshold & !is.na(binary01)
  out[idx] <- as.integer(binary01[idx] == 1L)
  out
}
or_status <- function(...) {
  m <- cbind(...)
  apply(m, 1, function(z) {
    if (any(z == 1L, na.rm = TRUE)) return(1L)
    if (all(!is.na(z)) && all(z == 0L)) return(0L)
    NA_integer_
  })
}
make_design <- function(d) {
  svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAF2YR, nest = TRUE, data = d)
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

dat2 <- dat %>%
  mutate(
    a1c_or_fpg_prediabetes = ifelse(is.na(a1c_prediabetes) & is.na(fpg_prediabetes), NA_integer_,
                                    as.integer((!is.na(a1c_prediabetes) & a1c_prediabetes == 1) |
                                                 (!is.na(fpg_prediabetes) & fpg_prediabetes == 1))),
    crit_prediabetes_a1c_fasting = bmi_plus_binary(BMXBMI, 27, a1c_prediabetes),
    crit_prediabetes_a1c_or_fpg = bmi_plus_binary(BMXBMI, 27, a1c_or_fpg_prediabetes),
    primary_a1c_fasting = or_status(crit_bmi35, crit_egfr, crit_prediabetes_a1c_fasting, crit_mi, crit_stroke),
    primary_a1c_or_fpg = or_status(crit_bmi35, crit_egfr, crit_prediabetes_a1c_or_fpg, crit_mi, crit_stroke)
  )

fast_d <- dat2 %>%
  filter(RIDAGEYR >= 20, medicare == 1, no_reported_diabetes == 1, !is.na(BMXBMI),
         WTSAF2YR > 0, !is.na(LBXGH), !is.na(LBXGLU),
         !is.na(primary_a1c_fasting), !is.na(primary_a1c_or_fpg)) %>%
  mutate(
    gained_by_fpg = as.integer(primary_a1c_fasting == 0 & primary_a1c_or_fpg == 1)
  )

des <- make_design(fast_d)
labels <- c(
  primary_a1c_fasting = "HbA1c-only prediabetes definition",
  primary_a1c_or_fpg = "HbA1c or fasting-glucose prediabetes definition"
)

main_tab <- lapply(names(labels), function(v) {
  p <- prop_est(des, v)
  tibble(
    outcome = v,
    characteristic = labels[[v]],
    n = nrow(fast_d),
    estimate = p$estimate,
    lower = p$lower,
    upper = p$upper,
    p_value = p$p_value,
    ci_label = p$ci_label,
    p_label = p$p_label
  )
}) %>% bind_rows()
write_csv(main_tab, file.path(OUT_DIR, "fasting_glucose_sensitivity_main.csv"))

des_aug_pos <- subset(des, primary_a1c_or_fpg == 1)
p_cond <- prop_est(des_aug_pos, "gained_by_fpg")
conditional_tab <- tibble(
  characteristic = "Added by fasting glucose among HbA1c/FPG-positive adults",
  estimate = p_cond$estimate,
  lower = p_cond$lower,
  upper = p_cond$upper,
  p_value = p_cond$p_value,
  ci_label = p_cond$ci_label,
  p_label = p_cond$p_label
)
write_csv(conditional_tab, file.path(OUT_DIR, "fasting_glucose_conditional_gain.csv"))
