# =============================================================
# 08_subgroups.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Pre-specified subgroup analysis and formal tests of
#          effect modification.
# Input:   data/processed/analysis_transformed.rds
#
# PRE-SPECIFIED SUBGROUPS (from the literature, not the results)
#   Sex             - Bohdanowicz-Pawlak (2023): BsmI and TaqI
#                     effects in postmenopausal women only
#   Vitamin D       - Xie (2020): FokI effect only in vitamin D
#                     deficient children; Salabat: ApaI-diabetes
#                     association only below 30 ng/ml
#   Obesity         - stated objective in the proposal
#   Glycaemic status- Salabat found ApaI associated with
#                     preDM/DM in this same cohort
#
# METHOD
#   The formal test of effect modification is the INTERACTION
#   term, assessed by nested model comparison. Stratified
#   estimates are reported for display only. A difference in
#   significance between strata is not a significant difference.
#
# CONTAMINATION RULE
#   A stratifying variable must not be a component of the
#   outcome. log_VAI is therefore excluded from the obesity
#   analysis, since BMI and waist circumference are inside the
#   VAI formula.
# =============================================================

source("R/00_setup.R")
library(broom)

data <- readRDS(file.path(path_processed, "analysis_transformed.rds"))


# ---- 1. Define the stratifying variables ----
# vitd_status and glycaemic already exist from script 03.
# Obesity is collapsed to binary: the underweight group (n = 22)
# is too small to stand alone.

data <- data %>%
  mutate(
    obesity = if_else(BMI >= 25, "overweight_obese", "normal"),
    obesity = factor(obesity, levels = c("normal", "overweight_obese")),
    # Collapse glycaemic to binary, as Salabat did for regression
    glyc2 = if_else(glycaemic == "healthy", "healthy", "preDM_DM"),
    glyc2 = factor(glyc2, levels = c("healthy", "preDM_DM"))
  )

strata_vars <- c("sex", "vitd_status", "obesity", "glyc2")

for (s in strata_vars) {
  cat("\n---", s, "---\n"); print(table(data[[s]], useNA = "ifany"))
}

# ---- 2. Check for sparse cells BEFORE modelling ----
# A genotype group with fewer than about 10 people in a stratum
# gives an unstable, uninterpretable coefficient.

genotype_vars <- c("apa_geno", "taq_geno", "ecorv_geno",
                   "fok_geno", "bsm_geno")

sparse <- do.call(rbind, lapply(strata_vars, function(s) {
  do.call(rbind, lapply(genotype_vars, function(g) {
    cells <- table(data[[g]], data[[s]])
    data.frame(Stratifier = s, Genotype = g,
               Smallest_cell = min(cells),
               Sparse = min(cells) < 10)
  }))
}))

sparse %>% filter(Sparse)

# ---- 3. Covariate sets ----
# Within a stratified model the stratifying variable no longer
# varies, so it is dropped from the covariates. The interaction
# model keeps it, because the product term requires its main
# effect to be present.

covariates_for <- function(outcome, stratifier) {
  
  covs <- c("age", "sex", "BMI", "VitD_Imputed")
  
  # BMI is inside the VAI formula
  if (outcome == "log_VAI") covs <- setdiff(covs, "BMI")
  
  # Drop whichever covariate duplicates the stratifier
  if (stratifier == "sex")         covs <- setdiff(covs, "sex")
  if (stratifier == "vitd_status") covs <- setdiff(covs, "VitD_Imputed")
  if (stratifier == "obesity")     covs <- setdiff(covs, "BMI")
  
  covs
}

# Which outcomes are valid for a given stratifier?
# log_VAI contains BMI and waist circumference, so it cannot be
# analysed when the stratifying variable is BMI-based.
valid_outcomes <- function(stratifier) {
  all_out <- c("AIP", "nonHDL_C", "log_LCI", "log_VAI")
  if (stratifier == "obesity") setdiff(all_out, "log_VAI") else all_out
}

# ---- 4. Interaction test ----
# Two nested models, differing ONLY by the interaction term:
#   A: outcome ~ genotype + stratifier + covariates
#   B: outcome ~ genotype * stratifier + covariates
#
# In R, "a * b" expands to "a + b + a:b" - both main effects
# plus their product. So model B contains everything in A, plus
# the product term. anova() then asks whether that addition
# explains significantly more variation.
#
# This gives ONE p-value for the whole interaction, which matters
# for a three-level genotype where the interaction spans two
# coefficients. Reading individual interaction coefficients would
# not test the question properly.

test_interaction <- function(outcome, predictor, stratifier, data) {
  
  covs <- covariates_for(outcome, stratifier)
  covs_text <- paste(covs, collapse = " + ")
  
  f_main <- as.formula(paste(outcome, "~", predictor, "+",
                             stratifier, "+", covs_text))
  f_int  <- as.formula(paste(outcome, "~", predictor, "*",
                             stratifier, "+", covs_text))
  
  m_main <- lm(f_main, data = data)
  m_int  <- lm(f_int,  data = data)
  
  comparison <- anova(m_main, m_int)
  
  data.frame(
    Stratifier    = stratifier,
    Outcome       = outcome,
    SNP           = predictor,
    n             = nobs(m_int),
    df_int        = comparison$Df[2],
    F_stat        = round(comparison$F[2], 3),
    p_interaction = comparison$`Pr(>F)`[2]
  )
}

# Test on one combination first
test_interaction("nonHDL_C", "ecorv_geno", "vitd_status", data)

# ---- 5. All interaction tests ----

predictors <- c(
  "apa_geno",   "apa_dom",   "apa_rec",
  "taq_geno",   "taq_dom",   "taq_rec",
  "ecorv_geno", "ecorv_dom", "ecorv_rec",
  "fok_geno",   "fok_dom",   "fok_rec",
  "bsm_geno",   "bsm_dom",   "bsm_rec"
)

interaction_results <- do.call(rbind, lapply(strata_vars, function(s) {
  combos <- expand.grid(outcome   = valid_outcomes(s),
                        predictor = predictors,
                        stringsAsFactors = FALSE)
  do.call(rbind,
          mapply(function(o, p) test_interaction(o, p, s, data),
                 combos$outcome, combos$predictor, SIMPLIFY = FALSE))
}))

rownames(interaction_results) <- NULL

# FDR within each stratifier - each is its own family of tests
interaction_results <- interaction_results %>%
  group_by(Stratifier) %>%
  mutate(p_FDR = p.adjust(p_interaction, method = "BH")) %>%
  ungroup()

# Overview
interaction_results %>%
  group_by(Stratifier) %>%
  summarise(tests        = n(),
            nominal_p05  = sum(p_interaction < 0.05),
            expected     = round(0.05 * n(), 1),
            survived_FDR = sum(p_FDR < 0.05),
            .groups = "drop")

# Anything nominally significant
interaction_results %>%
  filter(p_interaction < 0.05) %>%
  arrange(p_interaction)

# ---- 6. Stratified estimates (for display) ----
# These do NOT test effect modification - the interaction above
# does that. These show what the effect looks like in each group.

run_stratified <- function(outcome, predictor, stratifier,
                           stratum_label, data) {
  
  d <- data[!is.na(data[[stratifier]]) &
              data[[stratifier]] == stratum_label, ]
  
  covs <- covariates_for(outcome, stratifier)
  m <- lm(as.formula(paste(outcome, "~",
                           paste(c(predictor, covs), collapse = " + "))),
          data = d)
  
  tidy(m, conf.int = TRUE) %>%
    filter(startsWith(term, predictor)) %>%
    mutate(Stratifier = stratifier, Stratum = stratum_label,
           Outcome = outcome, SNP = predictor, n = nobs(m))
}

stratified <- do.call(rbind, lapply(strata_vars, function(s) {
  combos <- expand.grid(outcome   = valid_outcomes(s),
                        predictor = predictors,
                        stratum   = levels(data[[s]]),
                        stringsAsFactors = FALSE)
  do.call(rbind,
          mapply(function(o, p, st) run_stratified(o, p, s, st, data),
                 combos$outcome, combos$predictor, combos$stratum,
                 SIMPLIFY = FALSE))
}))

rownames(stratified) <- NULL

write.csv(interaction_results,
          file.path(path_tables, "interaction_tests.csv"), row.names = FALSE)
write.csv(stratified,
          file.path(path_tables, "stratified_estimates.csv"), row.names = FALSE)
#######

stratified %>%
  filter(SNP %in% c("ecorv_geno", "ecorv_rec"),
         Stratifier == "vitd_status", Outcome == "nonHDL_C") %>%
  select(Stratum, term, n, estimate, std.error, conf.low, conf.high, p.value)

stratified %>%
  filter(SNP %in% c("fok_geno", "fok_rec"),
         Stratifier == "sex", Outcome == "nonHDL_C") %>%
  select(Stratum, term, n, estimate, std.error, conf.low, conf.high, p.value)

#######
# Does the FokI x sex interaction persist after accounting for
# vitamin D status, or is it explained by it?
summary(lm(nonHDL_C ~ fok_rec * sex + fok_rec * vitd_status +
             age + BMI, data = data))

# Three-way: is the FokI effect specific to vitamin D deficient men?
summary(lm(nonHDL_C ~ fok_rec * sex * vitd_status + age + BMI, data = data))