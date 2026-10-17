#!/usr/bin/env Rscript

# Build analysis-ready dataset from NHANES August 2021-August 2023.

# Configuration and guarded writers supplied by src/runtime.R.

RESULT_DIR <- file.path(RESULT_ROOT, "audit")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
})

read_xpt_upper <- function(stem) {
  f <- file.path(RAW_DIR, paste0(stem, ".xpt"))
  if (!file.exists(f)) stop("Missing ", f, ". Obtain the eleven components manually as documented in README.md.")
  x <- haven::read_xpt(f)
  names(x) <- toupper(names(x))
  x
}

yes_no <- function(x) {
  # NHANES standard yes/no coding: 1 Yes, 2 No; 7/9/missing -> unknown
  out <- rep(NA_integer_, length(x))
  out[x == 1] <- 1L
  out[x == 2] <- 0L
  out
}

row_mean_any <- function(...) {
  m <- cbind(...)
  n_obs <- rowSums(!is.na(m))
  out <- rowMeans(m, na.rm = TRUE)
  out[n_obs == 0] <- NA_real_
  out
}

bmi_only_status <- function(bmi, threshold) {
  ifelse(is.na(bmi), NA_integer_, as.integer(bmi >= threshold))
}

bmi_plus_binary <- function(bmi, threshold, binary01) {
  out <- rep(NA_integer_, length(bmi))
  out[!is.na(bmi) & bmi < threshold] <- 0L
  idx <- !is.na(bmi) & bmi >= threshold & !is.na(binary01)
  out[idx] <- as.integer(binary01[idx] == 1L)
  out
}

bmi_plus_numeric_condition <- function(bmi, threshold, x, fun) {
  out <- rep(NA_integer_, length(bmi))
  out[!is.na(bmi) & bmi < threshold] <- 0L
  idx <- !is.na(bmi) & bmi >= threshold & !is.na(x)
  out[idx] <- as.integer(fun(x[idx]))
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

egfr_2021_ckdepi <- function(scr_mg_dl, age, sex_code) {
  # 2021 CKD-EPI creatinine equation (race-free)
  # sex_code: NHANES RIAGENDR, 1 = male, 2 = female
  kappa <- ifelse(sex_code == 2, 0.7,
                  ifelse(sex_code == 1, 0.9, NA_real_))
  alpha <- ifelse(sex_code == 2, -0.241,
                  ifelse(sex_code == 1, -0.302, NA_real_))
  female_factor <- ifelse(sex_code == 2, 1.012,
                          ifelse(sex_code == 1, 1, NA_real_))
  ratio <- scr_mg_dl / kappa

  ifelse(
    is.na(scr_mg_dl) | is.na(age) | is.na(kappa),
    NA_real_,
    142 *
      pmin(ratio, 1)^alpha *
      pmax(ratio, 1)^(-1.200) *
      0.9938^age *
      female_factor
  )
}

# ---- Read source datasets ----
demo <- read_xpt_upper("DEMO_L") %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, SDMVSTRA, SDMVPSU)

hiq <- read_xpt_upper("HIQ_L") %>%
  select(SEQN, HIQ032B, HIQ032D)

bmx <- read_xpt_upper("BMX_L") %>%
  select(SEQN, BMXBMI)

diq <- read_xpt_upper("DIQ_L") %>%
  select(SEQN, DIQ010, DIQ160)

mcq <- read_xpt_upper("MCQ_L") %>%
  select(SEQN, MCQ160B, MCQ160E, MCQ160F)

kiq <- read_xpt_upper("KIQ_U_L") %>%
  select(SEQN, KIQ022)

bpq <- read_xpt_upper("BPQ_L") %>%
  select(SEQN, BPQ150)

bpx <- read_xpt_upper("BPXO_L") %>%
  select(SEQN,
         BPXOSY1, BPXOSY2, BPXOSY3,
         BPXODI1, BPXODI2, BPXODI3)

ghb <- read_xpt_upper("GHB_L") %>%
  select(SEQN, WTPH2YR, LBXGH)

bio <- read_xpt_upper("BIOPRO_L") %>%
  select(SEQN, LBXSCR)

glu <- read_xpt_upper("GLU_L") %>%
  select(SEQN, WTSAF2YR, LBXGLU)

dat <- demo %>%
  left_join(hiq, by = "SEQN") %>%
  left_join(bmx, by = "SEQN") %>%
  left_join(diq, by = "SEQN") %>%
  left_join(mcq, by = "SEQN") %>%
  left_join(kiq, by = "SEQN") %>%
  left_join(bpq, by = "SEQN") %>%
  left_join(bpx, by = "SEQN") %>%
  left_join(ghb, by = "SEQN") %>%
  left_join(bio, by = "SEQN") %>%
  left_join(glu, by = "SEQN")

# ---- Derived variables ----
dat <- dat %>%
  mutate(
    medicare = as.integer(HIQ032B == 2),  # HIQ032B uses value 2 = covered by Medicare
    medicaid = as.integer(HIQ032D == 4),  # HIQ032D uses value 4 = covered by Medicaid

    diagnosed_diabetes = as.integer(DIQ010 == 1),
    diabetes_status_known = as.integer(DIQ010 %in% c(1, 2, 3)),
    no_reported_diabetes = as.integer(DIQ010 %in% c(2, 3)),

    dx_prediabetes = yes_no(DIQ160),
    dx_kidney = yes_no(KIQ022),
    dx_chf = yes_no(MCQ160B),
    dx_mi = yes_no(MCQ160E),
    dx_stroke = yes_no(MCQ160F),
    on_bp_med = yes_no(BPQ150),

    mean_sbp = row_mean_any(BPXOSY1, BPXOSY2, BPXOSY3),
    mean_dbp = row_mean_any(BPXODI1, BPXODI2, BPXODI3),
    bp_uncontrolled = ifelse(
      is.na(mean_sbp) & is.na(mean_dbp), NA_integer_,
      as.integer((!is.na(mean_sbp) & mean_sbp > 140) |
                   (!is.na(mean_dbp) & mean_dbp > 90))
    ),

    egfr = egfr_2021_ckdepi(LBXSCR, RIDAGEYR, RIAGENDR),

    a1c_prediabetes = ifelse(
      is.na(LBXGH), NA_integer_,
      as.integer(LBXGH >= 5.7 & LBXGH < 6.5)
    ),
    a1c_diabetes_range = ifelse(
      is.na(LBXGH), NA_integer_,
      as.integer(LBXGH >= 6.5)
    ),

    fpg_prediabetes = ifelse(
      is.na(LBXGLU), NA_integer_,
      as.integer(LBXGLU >= 100 & LBXGLU <= 125)
    ),
    fpg_diabetes_range = ifelse(
      is.na(LBXGLU), NA_integer_,
      as.integer(LBXGLU >= 126)
    ),

    # Primary survey-observable pathways
    crit_bmi35 = bmi_only_status(BMXBMI, 35),
    crit_egfr = bmi_plus_numeric_condition(
      BMXBMI, 30, egfr, function(x) x < 60
    ),
    crit_a1c_prediabetes = bmi_plus_binary(
      BMXBMI, 27, a1c_prediabetes
    ),
    crit_mi = bmi_plus_binary(BMXBMI, 27, dx_mi),
    crit_stroke = bmi_plus_binary(BMXBMI, 27, dx_stroke),

    primary_observable = or_status(
      crit_bmi35, crit_egfr, crit_a1c_prediabetes, crit_mi, crit_stroke
    ),

    # Diagnosis-proxy definition: replace measured prediabetes/kidney function
    # with self-reported prediabetes and weak/failing kidneys.
    crit_dx_kidney = bmi_plus_binary(BMXBMI, 30, dx_kidney),
    crit_dx_prediabetes = bmi_plus_binary(BMXBMI, 27, dx_prediabetes),

    diagnosis_proxy = or_status(
      crit_bmi35, crit_dx_kidney, crit_dx_prediabetes, crit_mi, crit_stroke
    ),

    # Expanded proxies used in the original sensitivity analysis.
    crit_chf_proxy = bmi_plus_binary(BMXBMI, 30, dx_chf),
    treated_uncontrolled_bp = ifelse(
      is.na(on_bp_med) | is.na(bp_uncontrolled), NA_integer_,
      as.integer(on_bp_med == 1 & bp_uncontrolled == 1)
    ),
    crit_bp_proxy = bmi_plus_binary(BMXBMI, 30, treated_uncontrolled_bp),

    expanded_proxy = or_status(
      primary_observable, crit_chf_proxy, crit_bp_proxy
    ),

    both_identified = ifelse(
      is.na(primary_observable) | is.na(diagnosis_proxy), NA_integer_,
      as.integer(primary_observable == 1 & diagnosis_proxy == 1)
    ),
    measurement_only = ifelse(
      is.na(primary_observable) | is.na(diagnosis_proxy), NA_integer_,
      as.integer(primary_observable == 1 & diagnosis_proxy == 0)
    ),
    diagnosis_only = ifelse(
      is.na(primary_observable) | is.na(diagnosis_proxy), NA_integer_,
      as.integer(primary_observable == 0 & diagnosis_proxy == 1)
    ),

    age_group = case_when(
      RIDAGEYR >= 20 & RIDAGEYR <= 64 ~ "20-64 y",
      RIDAGEYR >= 65 & RIDAGEYR <= 74 ~ "65-74 y",
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
    ),

    # Core study population used to reconstruct the primary analysis.
    # The original manuscript analyzed adults ≥20, reported Medicare coverage,
    # excluded self-reported diabetes, and used positive phlebotomy weights.
    base_eligible = (
      RIDAGEYR >= 20 &
      medicare == 1 &
      WTPH2YR > 0 &
      no_reported_diabetes == 1 &
      !is.na(BMXBMI)
    )
  )

saveRDS(dat, file.path(DERIVED_DIR, "bridge_analysis_ready.rds"))

# Small audit table for immediate inspection
audit <- tibble::tibble(
  step = c(
    "All NHANES records",
    "Age ≥20",
    "Medicare covered",
    "Positive phlebotomy weight",
    "No self-reported diabetes + BMI observed",
    "Primary status ascertainable",
    "Diagnosis proxy ascertainable",
    "Expanded proxy ascertainable",
    "Positive fasting-subsample weight"
  ),
  n = c(
    nrow(dat),
    sum(dat$RIDAGEYR >= 20, na.rm = TRUE),
    sum(dat$RIDAGEYR >= 20 & dat$medicare == 1, na.rm = TRUE),
    sum(dat$RIDAGEYR >= 20 & dat$medicare == 1 &
          dat$WTPH2YR > 0, na.rm = TRUE),
    sum(dat$base_eligible, na.rm = TRUE),
    sum(dat$base_eligible & !is.na(dat$primary_observable), na.rm = TRUE),
    sum(dat$base_eligible & !is.na(dat$diagnosis_proxy), na.rm = TRUE),
    sum(dat$base_eligible & !is.na(dat$expanded_proxy), na.rm = TRUE),
    sum(dat$RIDAGEYR >= 20 & dat$medicare == 1 &
          dat$no_reported_diabetes == 1 & !is.na(dat$BMXBMI) &
          dat$WTSAF2YR > 0, na.rm = TRUE)
  )
)
write_csv(audit, file.path(RESULT_DIR, "cohort_audit.csv"))

cat("Analysis-ready dataset written to:\n",
    file.path(DERIVED_DIR, "bridge_analysis_ready.rds"), "\n\n")
print(audit)
