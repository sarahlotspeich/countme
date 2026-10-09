# Load libraries
library(MASS) ## for glm.nb
library(countme) ## for smle_nb and smle_zi_nb
library(splines) ## for bs
library(dplyr) ## for mutate, if_else, rename, bind_cols
library(bizicount) ## for zic.reg
library(ggplot2) ## for plots

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
n_fit = glm.nb(formula = NUM_ADMIT ~ UNVAL_ALI + AGE_AT_ENCOUNTER + SEX,
               data = hosp_data)
### Save estimates for forest plot ---------------------------------------------
n_coeff = n_fit |>
  summary() |>
  coef() |>
  data.frame() |>
  bind_cols(
    confint(n_fit) |>
      data.frame()
  ) |>
  rename(SE = Std..Error,
         p.value = Pr...z..,
         ci.lb = X2.5..,
         ci.ub = X97.5..)  |>
  mutate(term = toupper(names(n_fit$coefficients))) |>
  bind_rows(
    data.frame(Estimate = n_fit$theta,
               SE = n_fit$SE.theta,
               z.value = n_fit$theta / n_fit$SE.theta,
               p.value = 2 * pnorm(q = abs(n_fit$theta / n_fit$SE.theta), lower.tail = FALSE),
               row.names = "DISPERSION") |>
      mutate(ci.lb = Estimate - 1.96 * SE,
             ci.ub = Estimate + 1.96 * SE,
             term = "DISPERSION")
  )
### Zero-inflated negative binomial --------------------------------------------
n_fit_zi = zic.reg(
  fmla = NUM_ADMIT ~ UNVAL_ALI + AGE_AT_ENCOUNTER + SEX | AGE_AT_ENCOUNTER,
  data = hosp_data,
  dist = "nbinom",
  optimizer = "optim"
)
n_zi_coeff = data.frame(
  Estimate = summary(n_fit_zi)$coef,
  SE = summary(n_fit_zi)$se,
  z.value = summary(n_fit_zi)$coef / summary(n_fit_zi)$se,
  p.value = 2 * pnorm(q = abs(summary(n_fit_zi)$coef / summary(n_fit_zi)$se),
                      lower.tail = FALSE)) |>
  mutate(ci.lb = Estimate - 1.96 * SE,
         ci.ub = Estimate + 1.96 * SE,
         term = toupper(names(summary(n_fit_zi)$coef)))
n_zi_coeff

## Fit complete-case models
### Weighting model for validation sampling
val_wt_fit = glm(formula = VALIDATED ~ ANY_ADMIT + UNVAL_ALI + AGE_AT_ENCOUNTER,
                 family = "binomial",
                 data = hosp_data)
hosp_data$val_wt = 1 / predict(object = val_wt_fit,
                               type = "response")
### Negative binomial
cc_hosp_data = hosp_data[!is.na(hosp_data$VAL_ALI), ]
cc_fit = glm.nb(formula = NUM_ADMIT ~ VAL_ALI + AGE_AT_ENCOUNTER + SEX,
                data = cc_hosp_data,
                weights = cc_hosp_data$val_wt)
#### Save estimates for forest plot
cc_coeff = cc_fit |>
  summary() |>
  coef() |>
  data.frame() |>
  bind_cols(
    confint(cc_fit) |>
      data.frame()
  ) |>
  rename(SE = Std..Error,
         p.value = Pr...z..,
         ci.lb = X2.5..,
         ci.ub = X97.5..) |>
  mutate(term = toupper(names(cc_fit$coefficients))) |>
  bind_rows(
    data.frame(Estimate = cc_fit$theta,
               SE = cc_fit$SE.theta,
               z.value = cc_fit$theta / cc_fit$SE.theta,
               p.value = 2 * pnorm(q = abs(cc_fit$theta / cc_fit$SE.theta), lower.tail = FALSE),
               row.names = "DISPERSION"
    ) |>
      mutate(ci.lb = Estimate - 1.96 * SE,
             ci.ub = Estimate + 1.96 * SE,
             term = "DISPERSION")
  )
### Zero-inflated negative binomial
cc_fit_zi = zic.reg(
  fmla = NUM_ADMIT ~ VAL_ALI + AGE_AT_ENCOUNTER + SEX | AGE_AT_ENCOUNTER,
  data = cc_hosp_data,
  weights = cc_hosp_data$val_wt,
  dist = "nbinom",
  optimizer = "optim" #### Didn't return SEs with "nlm"
)

cc_zi_coeff = data.frame(
  Estimate = summary(cc_fit_zi)$coef,
  SE = summary(cc_fit_zi)$se,
  z.value = summary(cc_fit_zi)$coef / summary(cc_fit_zi)$se,
  p.value = 2 * pnorm(q = abs(summary(cc_fit_zi)$coef / summary(cc_fit_zi)$se),
                      lower.tail = FALSE)) |>
  mutate(ci.lb = Estimate - 1.96 * SE,
         ci.ub = Estimate + 1.96 * SE,
         term = toupper(names(summary(cc_fit_zi)$coef)))
cc_zi_coeff

## Fit SMLE model
### Setup B-splines on UNVAL_ALI (X*)
sn = 14
B = bs(x = hosp_data$UNVAL_ALI, ## Error-prone ALI (from EHR)
       df = sn, ## number of sieves (decided based on Phase I/Phase II sample sizes)
       Boundary.knots = range(hosp_data$UNVAL_ALI),
       intercept = TRUE)
colnames(B) = paste0("bs", seq(1, sn))
colSums(B[!is.na(hosp_data$VAL_ALI), ]) #### Ensure that there are no missing sieves in validated data

# Column bind B-splines to data
hosp_data = cbind(hosp_data, B)
### Negative binomial
smle_fit = smle_nb(analysis_formula = NUM_ADMIT ~ VAL_ALI + AGE_AT_ENCOUNTER + SEX,
                   error_formula = paste("VAL_ALI ~", paste(paste0("bs", 1:sn), collapse = "+")),
                   data = hosp_data,
                   no_se = FALSE)
smle_coeff = smle_fit$coefficients |>
  data.frame() |>
  rename(SE = Std..Error,
         p.value = Pr...z..) |>
  mutate(ci.lb = Estimate - 1.96 * SE,
         ci.ub = Estimate + 1.96 * SE,
         term = toupper(rownames(smle_fit$coefficients)))
### Zero-inflated negative binomial
smle_fit_zi = smle_zi_nb(analysis_formula = NUM_ADMIT ~ VAL_ALI + AGE_AT_ENCOUNTER + SEX | AGE_AT_ENCOUNTER,
                         error_formula = paste("VAL_ALI ~", paste(paste0("bs", 1:sn), collapse = "+")),
                         data = hosp_data,
                         no_se = FALSE,
                         optimizer = "optim")
smle_zi_coeff = smle_fit_zi$coefficients |>
  data.frame() |>
  rename(SE = Std..Error,
         p.value = Pr...z..) |>
  mutate(ci.lb = Estimate - 1.96 * SE,
         ci.ub = Estimate + 1.96 * SE,
         term = toupper(rownames(smle_fit_zi$coefficients)))

# Make a forest plot of estimates (95% confidence intervals)
n_coeff |>
  mutate(Model = "Negative Binomial", Analysis = "Naive") |>
  bind_rows(
    n_zi_coeff |>
      mutate(Model = "Zero-Inflated Negative Binomial", Analysis = "Naive")
  ) |>
  bind_rows(
    cc_coeff |>
      mutate(Model = "Negative Binomial", Analysis = "IPW")
  )  |>
  bind_rows(
    cc_zi_coeff |>
      mutate(Model = "Zero-Inflated Negative Binomial", Analysis = "IPW")
  ) |>
  bind_rows(
    smle_coeff |>
      mutate(Model = "Negative Binomial", Analysis = "SMLE")
  )  |>
  bind_rows(
    smle_zi_coeff |>
      mutate(Model = "Zero-Inflated Negative Binomial", Analysis = "SMLE")
  ) |>
  mutate(term = sub("\\(", "", term),
         term = sub("\\)", "", term),
         term = sub("CT_|COUNT:", "", term),
         term = sub("ZI_", "ZERO:", term),
         term = ifelse(term == "ZERO:INT", "ZERO:INTERCEPT", term),
         term = ifelse(term == "INT", "INTERCEPT", term),
         term = ifelse(term == "THETA", "DISPERSION", term),
         term = sub("UNVAL_", "", term),
         term = sub("VAL_", "", term)) |>
  filter(!grepl(pattern = "ZERO:", term)) |>
  mutate(term = factor(x = term,
                       levels = c("INTERCEPT", "ALI", "AGE_AT_ENCOUNTER", "SEXMALE", "DISPERSION"),
                       labels = c("Intercept (male, 18 years old, ALI of 0)", "ALI (per component)", "Age (per 10 years)", "Sex (male versus female)", "Dispersion"))
  ) |>
  ggplot(aes(x = Model, y = Estimate, color = Analysis,
             shape = Model, group = interaction(Analysis, Model))) +
  geom_hline(aes(yintercept = as.numeric(term == "Dispersion", 1, 0)),
             linetype = "dashed", color = "gray") +
  geom_point(position = position_dodge(1)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub),
                position = position_dodge(1)) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top",
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white", face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggthemes::scale_color_colorblind() +
  facet_wrap(~term, scales = "free")
ggsave(filename = "~/Documents/countme/hospital_analysis/forest_plot.png",
       device = "png", width = 10, height = 7, units = "in")
