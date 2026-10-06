# =============================================================
# 07_models.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Association between VDR polymorphisms and lipid
#          indices, using sequential adjustment.
# Input:   data/processed/analysis_transformed.rds
#
# ADJUSTMENT SEQUENCE (following Salabat et al.)
#   Model 1: crude
#   Model 2: + age, sex
#   Model 3: + BMI
#   Model 4: + vitamin D
# =============================================================

source("R/00_setup.R")
library(broom)

data <- readRDS(file.path(path_processed, "analysis_transformed.rds"))


# ---- 1. Model 1: crude ----
m1 <- lm(AIP ~ apa_geno, data = data)
summary(m1)

m2 <- lm(AIP ~ apa_geno + age + sex, data = data)
m3 <- lm(AIP ~ apa_geno + age + sex + BMI, data = data)
m4 <- lm(AIP ~ apa_geno + age + sex + BMI + VitD_Imputed, data = data)

# ---- 2. Models 2 to 4: sequential adjustment ----

m2 <- lm(AIP ~ apa_geno + age + sex, data = data)
m3 <- lm(AIP ~ apa_geno + age + sex + BMI, data = data)
m4 <- lm(AIP ~ apa_geno + age + sex + BMI + VitD_Imputed, data = data)

summary(m2)
summary(m3)
summary(m4)

# ---- 3. Tidy output ----
# tidy() turns the coefficient table into a data frame, which is
# what makes automation possible in the next session.

tidy(m4, conf.int = TRUE)
glance(m4)

nobs(m1)
nobs(m4)

# Does genotype as a whole explain anything, beyond covariates?

# Define the complete-case sample explicitly, then fit both
# models to it. Nested model comparison requires identical data.

vars_needed <- c("AIP", "apa_geno", "age", "sex", "BMI", "VitD_Imputed")
data_cc <- data[complete.cases(data[, vars_needed]), ]

nrow(data_cc)

m4_full   <- lm(AIP ~ apa_geno + age + sex + BMI + VitD_Imputed, data = data_cc)
m4_nogeno <- lm(AIP ~ age + sex + BMI + VitD_Imputed, data = data_cc)

anova(m4_nogeno, m4_full)