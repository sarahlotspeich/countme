# Load packages
library(dplyr)
library(tidyr)

# Load in unvalidated ALI data + hospitalizations
unvalALI_hosp = read.csv("~/Documents/Allostatic_load_audits/deidentified_ali_hospitalizations.csv") |>
  group_by(PAT_MRN_ID) |>
  mutate(UNVAL_ALI = mean(c(A1C, ALB, BMI, CHOL, CRP, CREAT_C, HCST, TRIG, BP_DIASTOLIC, BP_SYSTOLIC), na.rm = TRUE)) |>
  select(PAT_MRN_ID, NUM_ADMIT, NUM_ED, AGE_AT_ENCOUNTER, SEX, RACE, ETHNICITY, UNVAL_ALI)

# Load patient ID crosswalk
patID_cw = read.csv("~/Documents/Allostatic_load_audits/deidentified_ali_hospitalizations_crosswalk.csv")

# Load in Pilot + Wave I validated ALI data
valALI = read.csv("~/Documents/Allostatic_load_audits/Grayson's Thesis/summary_data_audit1_grayson_thesis.csv") |>
  dplyr::rename(ORIG_PAT_MRN_ID = PAT_MRN_ID) |>
  dplyr::bind_rows(
    read.csv("~/Documents/Allostatic_load_audits/Grayson's Thesis/summary_data_audit2_grayson_thesis.csv") |>
      dplyr::rename(ORIG_PAT_MRN_ID = PAT_MRN_ID)
  ) |>
  dplyr::left_join(patID_cw)

## Stratification variables for the Seemen et al. ALI
strat_vals = read.csv("~/Documents/Allostatic_load_audits/Audit_Protocol/strat_vals.csv")

# Create a long version with one row per variable per patient ------------------
valALI_long = valALI |>
  gather(key = "COMP",  value = "VAL", -c(1, 12)) |>
  left_join(unvalALI_hosp |>
              select(PAT_MRN_ID, SEX)) |> # Merge in sex (needed for the stratification values)
  left_join(strat_vals) |> ## Merge in the stratification values
  mutate(
    FLAG = if_else(condition = !is.na(VAL) & VAL == -1,
                   true = TRUE,
                   false = FALSE), ## If component = -1 it was only found by auditors outside of study period
    VAL = if_else(condition = FLAG, ## Then replace flagged values with NA, so that
                  true = NA, ### They will be excluded from ALI calculation
                  false = VAL),
    POINT = case_when( ### Create 0/1 "points" for the component stressors
      INEQ == ">" ~ as.numeric(VAL > STRAT),
      INEQ == ">=" ~ as.numeric(VAL >= STRAT),
      INEQ == "<" ~ as.numeric(VAL < STRAT))
  )

valALI = valALI_long |>
  group_by(ORIG_PAT_MRN_ID, PAT_MRN_ID) |>
  summarize(VAL_ALI = mean(POINT, na.rm = TRUE)) |>
  select(ORIG_PAT_MRN_ID, PAT_MRN_ID, VAL_ALI)

# Load in unvalidated ALI data + hospitalizations
ALI_hosp = unvalALI_hosp |>
  left_join(valALI)
ALI_hosp |>
  select(-ORIG_PAT_MRN_ID) |>
  write.csv("~/Dropbox (Wake Forest University)/4 - PAPERS/Sellers-KMP-ME/hospital_analysis_data.csv",
            row.names = FALSE)
