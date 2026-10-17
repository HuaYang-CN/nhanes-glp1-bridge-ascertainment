
#!/usr/bin/env Rscript

# Figure 1 v4: robust subgroup domain CIs + meaningful overall subgroup P values.
# IMPORTANT:
# P values in Panel D are Rao-Scott survey-adjusted overall tests for
# heterogeneity across levels of each subgroup variable.
#
# This version also preserves the original complex survey design when
# calculating subgroup estimates.

# Configuration and guarded writers supplied by src/runtime.R.

OUT_DIR <- file.path(RESULT_ROOT, "primary")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
  library(tibble)
})

options(survey.lonely.psu = "adjust")

dat <- readRDS(file.path(DERIVED_DIR, "bridge_analysis_ready.rds"))

make_design <- function(d, weight_var = "WTPH2YR") {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", weight_var)),
    nest = TRUE,
    data = d
  )
}

fmt_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

prop_est <- function(design, var, df_override = NULL) {
  f <- as.formula(paste0("~", var))
  df_use <- if (is.null(df_override)) degf(design) else df_override

  # Point estimate + SE from the Taylor-linearized survey mean.
  m <- svymean(f, design = design, na.rm = TRUE)
  est <- as.numeric(coef(m)[1])
  se  <- as.numeric(SE(m)[1])

  # Preferred CI: logit CI. If a sparse domain causes a non-finite result,
  # try beta CI, then a clipped t-based Wald CI.
  methods <- c("logit", "beta")
  for (meth in methods) {
    z <- try(
      svyciprop(
        f,
        design = design,
        method = meth,
        df = df_use,
        na.rm = TRUE
      ),
      silent = TRUE
    )

    if (!inherits(z, "try-error")) {
      ci <- try(confint(z), silent = TRUE)
      if (!inherits(ci, "try-error") &&
          length(ci) >= 2 &&
          all(is.finite(as.numeric(ci[1:2])))) {
        return(tibble(
          estimate = 100 * est,
          lower = 100 * as.numeric(ci[1]),
          upper = 100 * as.numeric(ci[2]),
          ci_method = meth
        ))
      }
    }
  }

  if (is.finite(est) && is.finite(se) && is.finite(df_use) && df_use > 0) {
    crit <- qt(0.975, df = df_use)
    lo <- max(0, est - crit * se)
    hi <- min(1, est + crit * se)
    return(tibble(
      estimate = 100 * est,
      lower = 100 * lo,
      upper = 100 * hi,
      ci_method = "wald_fallback"
    ))
  }

  tibble(
    estimate = 100 * est,
    lower = NA_real_,
    upper = NA_real_,
    ci_method = "not_estimable"
  )
}

# Domain proportion calculated as a survey-weighted ratio in the FULL design.
# This is preferable to rebuilding/subsetting the design for sparse subgroups:
# observations outside the domain remain in the design with domain indicator 0,
# preserving the full NHANES strata/PSU structure and design degrees of freedom.
domain_prop_est <- function(design, outcome_var, group_var, level) {
  vars <- design$variables

  in_domain <- !is.na(vars[[group_var]]) & vars[[group_var]] == level
  outcome <- vars[[outcome_var]]

  numerator <- ifelse(in_domain & !is.na(outcome), outcome, 0)
  denominator <- as.numeric(in_domain & !is.na(outcome))

  des2 <- design
  des2$variables$.domain_num <- numerator
  des2$variables$.domain_den <- denominator

  rr <- try(
    svyratio(
      ~.domain_num,
      ~.domain_den,
      design = des2,
      na.rm = TRUE
    ),
    silent = TRUE
  )

  n_unweighted <- sum(denominator == 1, na.rm = TRUE)

  if (inherits(rr, "try-error")) {
    return(tibble(
      n = n_unweighted,
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      ci_method = "not_estimable"
    ))
  }

  est <- as.numeric(coef(rr)[1])
  se <- as.numeric(SE(rr)[1])
  df_use <- degf(design)
  crit <- qt(0.975, df = df_use)

  # Preferred domain CI: logit-transformed Taylor CI.
  # This stays inside [0,1] and uses the FULL NHANES design df.
  if (is.finite(est) && is.finite(se) &&
      est > 0 && est < 1 &&
      is.finite(df_use) && df_use > 0) {
    logit_est <- qlogis(est)
    se_logit <- se / (est * (1 - est))

    if (is.finite(se_logit)) {
      lo <- plogis(logit_est - crit * se_logit)
      hi <- plogis(logit_est + crit * se_logit)

      if (all(is.finite(c(lo, hi)))) {
        return(tibble(
          n = n_unweighted,
          estimate = 100 * est,
          lower = 100 * lo,
          upper = 100 * hi,
          ci_method = "domain_ratio_logit"
        ))
      }
    }
  }

  # Fallback: clipped t-based Taylor CI.
  if (is.finite(est) && is.finite(se) &&
      is.finite(df_use) && df_use > 0) {
    lo <- max(0, est - crit * se)
    hi <- min(1, est + crit * se)

    return(tibble(
      n = n_unweighted,
      estimate = 100 * est,
      lower = 100 * lo,
      upper = 100 * hi,
      ci_method = "domain_ratio_wald"
    ))
  }

  tibble(
    n = n_unweighted,
    estimate = 100 * est,
    lower = NA_real_,
    upper = NA_real_,
    ci_method = "not_estimable"
  )
}

total_est <- function(design, var) {
  x <- svytotal(as.formula(paste0("~", var)), design = design, na.rm = TRUE)
  est <- as.numeric(coef(x)[1])
  se <- as.numeric(SE(x)[1])
  crit <- qt(0.975, df = degf(design))
  tibble(
    estimate_millions = est / 1e6,
    lower_millions = (est - crit * se) / 1e6,
    upper_millions = (est + crit * se) / 1e6
  )
}

estimate_row <- function(d, var, label, weight = "WTPH2YR") {
  des <- make_design(d, weight)
  p <- prop_est(des, var)
  t <- total_est(des, var)
  tibble(
    characteristic = label,
    n = nrow(d),
    estimate = p$estimate,
    lower = p$lower,
    upper = p$upper,
    ci_method = p$ci_method,
    ci_label = ifelse(
      is.finite(p$lower) & is.finite(p$upper),
      sprintf("%.1f (%.1f–%.1f)", p$estimate, p$lower, p$upper),
      sprintf("%.1f (CI not estimable)", p$estimate)
    ),
    population_millions = t$estimate_millions,
    population_lower_millions = t$lower_millions,
    population_upper_millions = t$upper_millions
  )
}

# ----------------------------------------------------------------------
# A. Primary + sensitivity analyses
# ----------------------------------------------------------------------
primary_d <- dat %>%
  filter(base_eligible, !is.na(primary_observable))

expanded_d <- dat %>%
  filter(base_eligible, !is.na(expanded_proxy))

a1c_exclusion_d <- primary_d %>%
  filter(is.na(LBXGH) | LBXGH < 6.5)

table_A <- bind_rows(
  estimate_row(primary_d, "primary_observable", "Primary observable criteria"),
  estimate_row(expanded_d, "expanded_proxy", "Expanded clinical proxies"),
  estimate_row(a1c_exclusion_d, "primary_observable", "Exclude HbA1c ≥ 6.5%")
)

write_csv(table_A, file.path(OUT_DIR, "Figure1A_primary_and_sensitivity.csv"))

# ----------------------------------------------------------------------
# B. Clinical qualification routes
# ----------------------------------------------------------------------
primary_design <- make_design(primary_d, "WTPH2YR")

pathways <- c(
  crit_bmi35 = "BMI ≥ 35",
  crit_egfr = "BMI ≥ 30 + eGFR <60 mL/min/1.73 m²",
  crit_a1c_prediabetes = "BMI ≥ 27 + prediabetes",
  crit_mi = "BMI ≥ 27 + prior MI",
  crit_stroke = "BMI ≥ 27 + prior stroke",
  crit_chf_proxy = "BMI ≥ 30 + heart failure proxy",
  crit_bp_proxy = "BMI ≥ 30 + uncontrolled BP proxy"
)

table_B <- lapply(names(pathways), function(v) {
  p <- prop_est(primary_design, v)
  tibble(
    characteristic = unname(pathways[[v]]),
    n = nrow(primary_d),
    estimate = p$estimate,
    lower = p$lower,
    upper = p$upper,
    ci_method = p$ci_method,
    ci_label = ifelse(
      is.finite(p$lower) & is.finite(p$upper),
      sprintf("%.1f (%.1f–%.1f)", p$estimate, p$lower, p$upper),
      sprintf("%.1f (CI not estimable)", p$estimate)
    )
  )
}) %>% bind_rows()

write_csv(table_B, file.path(OUT_DIR, "Figure1B_criterion_pathways.csv"))

# ----------------------------------------------------------------------
# C. Identification strategy
# ----------------------------------------------------------------------
id_d <- primary_d %>% filter(!is.na(diagnosis_proxy))
id_design <- make_design(id_d, "WTPH2YR")

id_vars <- c(
  primary_observable = "Measurement-based criteria",
  diagnosis_proxy = "Self-reported diagnosis proxy",
  both_identified = "Identified by both",
  measurement_only = "Measurement only",
  diagnosis_only = "Diagnosis proxy only"
)

table_C <- lapply(names(id_vars), function(v) {
  p <- prop_est(id_design, v)
  tibble(
    characteristic = unname(id_vars[[v]]),
    n = nrow(id_d),
    estimate = p$estimate,
    lower = p$lower,
    upper = p$upper,
    ci_method = p$ci_method,
    ci_label = ifelse(
      is.finite(p$lower) & is.finite(p$upper),
      sprintf("%.1f (%.1f–%.1f)", p$estimate, p$lower, p$upper),
      sprintf("%.1f (CI not estimable)", p$estimate)
    )
  )
}) %>% bind_rows()

write_csv(table_C, file.path(OUT_DIR, "Figure1C_identification_strategy.csv"))

# ----------------------------------------------------------------------
# D. Faithful descriptive subgroup reproduction (planning timing unverified)
#    IMPORTANT: create the survey design ONCE, then conduct domain analyses.
# ----------------------------------------------------------------------
subgroup_d <- primary_d %>%
  mutate(
    age_group = case_when(
      RIDAGEYR >= 20 & RIDAGEYR <= 64 ~ "20–64 y",
      RIDAGEYR >= 65 & RIDAGEYR <= 74 ~ "65–74 y",
      RIDAGEYR >= 75 ~ "≥75 y",
      TRUE ~ NA_character_
    ),
    sex = case_when(
      RIAGENDR == 1 ~ "Men",
      RIAGENDR == 2 ~ "Women",
      TRUE ~ NA_character_
    ),
    race_ethnicity = case_when(
      RIDRETH3 == 1 ~ "Mexican American",
      RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White",
      RIDRETH3 == 4 ~ "Non-Hispanic Black",
      RIDRETH3 == 6 ~ "Non-Hispanic Asian",
      RIDRETH3 == 7 ~ "Other or multiracial",
      TRUE ~ NA_character_
    ),
    income_group = case_when(
      !is.na(INDFMPIR) & INDFMPIR < 2 ~ "Income-to-poverty ratio < 2",
      !is.na(INDFMPIR) & INDFMPIR >= 2 ~ "Income-to-poverty ratio ≥ 2",
      TRUE ~ NA_character_
    ),
    coverage_group = case_when(
      medicare == 1 & medicaid == 1 ~ "Medicare and Medicaid",
      medicare == 1 ~ "Medicare only",
      TRUE ~ NA_character_
    )
  )

subgroup_design <- make_design(subgroup_d, "WTPH2YR")

subgroup_levels <- list(
  Age = c("20–64 y", "65–74 y", "≥75 y"),
  Sex = c("Men", "Women"),
  `Race and ethnicity` = c(
    "Mexican American", "Other Hispanic", "Non-Hispanic White",
    "Non-Hispanic Black", "Non-Hispanic Asian", "Other or multiracial"
  ),
  Income = c("Income-to-poverty ratio < 2", "Income-to-poverty ratio ≥ 2"),
  Coverage = c("Medicare and Medicaid", "Medicare only")
)

subgroup_var <- c(
  Age = "age_group",
  Sex = "sex",
  `Race and ethnicity` = "race_ethnicity",
  Income = "income_group",
  Coverage = "coverage_group"
)

domain_design <- function(design, var, value) {
  keep <- !is.na(design$variables[[var]]) & design$variables[[var]] == value
  design[keep, ]
}

overall_p <- function(design, var) {
  vars <- design$variables
  keep <- !is.na(vars[[var]]) & !is.na(vars$primary_observable)
  des <- design[keep, ]

  # Rao-Scott survey-adjusted overall association test.
  f <- as.formula(paste0("~factor(primary_observable) + ", var))
  tst <- try(svychisq(f, design = des, statistic = "F"), silent = TRUE)

  if (inherits(tst, "try-error")) return(NA_real_)
  as.numeric(tst$p.value)
}

table_D_parts <- list()
p_parts <- list()
k <- 1

for (section in names(subgroup_levels)) {
  var <- subgroup_var[[section]]

  pval <- overall_p(subgroup_design, var)
  p_parts[[section]] <- tibble(
    section = section,
    p_value = pval,
    p_label = fmt_p(pval)
  )

  for (lev in subgroup_levels[[section]]) {
    # Estimate the subgroup prevalence as a DOMAIN ratio in the full survey design.
    # This avoids loss of design df in sparse groups (e.g., Mexican American).
    p <- domain_prop_est(
      subgroup_design,
      outcome_var = "primary_observable",
      group_var = var,
      level = lev
    )

    if (p$n == 0) next

    table_D_parts[[k]] <- tibble(
      section = section,
      characteristic = lev,
      n = p$n,
      estimate = p$estimate,
      lower = p$lower,
      upper = p$upper,
      ci_method = p$ci_method,
      ci_label = ifelse(
        is.finite(p$lower) & is.finite(p$upper),
        sprintf("%.1f (%.1f–%.1f)", p$estimate, p$lower, p$upper),
        sprintf("%.1f (CI not estimable)", p$estimate)
      )
    )
    k <- k + 1
  }
}

table_D <- bind_rows(table_D_parts)
table_D_p <- bind_rows(p_parts)

write_csv(table_D, file.path(OUT_DIR, "Figure1D_subgroups.csv"))
write_csv(table_D_p, file.path(OUT_DIR, "Figure1D_subgroup_overall_P.csv"))

cat("\nPanel D overall subgroup P values:\n")
print(table_D_p)
cat("\nPanel D subgroup estimates:\n")
print(table_D)

nonfinite_ci <- table_D %>%
  filter(!is.finite(lower) | !is.finite(upper))

if (nrow(nonfinite_ci) > 0) {
  warning(
    "One or more Panel D confidence intervals remain non-estimable after ",
    "full-design domain estimation. They have NOT been fabricated. ",
    "See Figure1D_subgroups.csv and ci_method."
  )
  print(nonfinite_ci)
}

cat("\nCI methods used in Panel D:\n")
print(table(table_D$ci_method, useNA = "ifany"))
