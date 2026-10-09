# Load libraries
library(MASS) ## for glm.nb
library(dplyr) ## for mutate, if_else, rename, bind_cols
library(bizicount) ## for zic.reg

# Read in data
hosp_data = read.csv("~/Documents/countme/hospital_analysis/hospital_analysis_data.csv")
## Rescale age to be in 10-year increments (rather than 1-year)
hosp_data = hosp_data |>
  mutate(AGE_AT_ENCOUNTER = (AGE_AT_ENCOUNTER - 18) / 10,
         UNVAL_ALI = UNVAL_ALI / 0.1,
         VAL_ALI = VAL_ALI / 0.1,
         VALIDATED = as.numeric(!is.na(VAL_ALI)),
         ANY_ADMIT = as.numeric(NUM_ADMIT > 0))
## Collapse smaller race categories (American Indian or Alaska Native --> Other)
table(hosp_data$RACE)
hosp_data = hosp_data |>
  mutate(RACE = if_else(condition = RACE == "American Indian or Alaska Native",
                        true = "Other",
                        false = RACE),
         RACE = factor(x = RACE,
                       levels = c("White or Caucasian", "Black or African American", "Asian Indian", "Other")))
table(hosp_data$RACE)

# Fit models of hospitalizations ~ ALI + covariates (age, sex, race) -----------
## Fit naive model (using unvalidated ALI on all 1000 patients ) ---------------
### Negative binomial ----------------------------------------------------------
n_fit = glm.nb(formula = NUM_ADMIT ~ UNVAL_ALI + SEX,
               data = hosp_data)
n_fit$coefficients
# (Intercept)   UNVAL_ALI     SEXMale
# -1.6836643   0.1710420   0.4340504
# > n_fit$theta
# [1] 0.2515643

### Zero-inflated negative binomial --------------------------------------------
n_fit_zi = zic.reg(
  fmla = NUM_ADMIT ~ UNVAL_ALI + SEX | SEX,
  data = hosp_data,
  dist = "nbinom",
  optimizer = "optim"
)
n_fit_zi$coef
# ct_(Intercept)  ct_UNVAL_ALI   ct_SEXMale     zi_(Intercept) zi_SEXMale     Theta
# -1.7361879      0.1860665      0.8064816     -7.3877147      6.5799973      0.3259274
