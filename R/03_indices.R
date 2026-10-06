# =============================================================
# 03_indices.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Compute the five lipid/atherogenic indices named in
#          the proposal.
# Input:   data/processed/geno.rds
# Output:  data/processed/analysis.rds
#
# UNITS - READ BEFORE EDITING
#   Source lipids are mg/dL. The formulas do NOT share units:
#     AIP, VAI       -> mmol/L
#     TG/HDL         -> unit-free (both in the same unit)
#     non-HDL-C, LCI -> mg/dL
#   We therefore keep both versions. This is deliberate, not
#   redundant. Using mg/dL in the VAI formula inflates the
#   result about 30-fold.
# =============================================================

source("R/00_setup.R")

data <- readRDS(file.path(path_processed, "geno.rds"))


# ---- 1. Convert to mmol/L ----
# Divisors differ because the molecular weights differ.

data <- data %>%
  mutate(
    TG_mmol  = TG  / 88.57,
    HDL_mmol = HDL / 38.67
  )

summary(data$TG)
summary(data$TG_mmol)

# ---- 2. Was LDL measured or calculated? ----
# Many labs do not measure LDL directly - they calculate it with
# the Friedewald equation from TC, HDL and TG.
#
#   Friedewald (mg/dL): LDL = TC - HDL - (TG / 5)
#
# This matters for LCI. If LDL was calculated, it contains no
# information beyond the other three lipids, and LCI becomes
# largely a rearrangement of them rather than a new measurement.
# Either way you must state which in your Methods.

data <- data %>%
  mutate(LDL_friedewald = chol - HDL - (TG / 5))

# If LDL was CALCULATED, the difference is ~0 for everyone.
# If MEASURED, the differences scatter.
summary(data$LDL - data$LDL_friedewald)

# ---- 3. AIP - Atherogenic Index of Plasma ----
# AIP = log10( TG / HDL-C ), both in mmol/L. Unitless.
# Published risk bands assume the mmol/L basis:
#   < 0.11 low | 0.11 to 0.21 intermediate | > 0.21 high

data <- data %>%
  mutate(AIP = log10(TG_mmol / HDL_mmol))


# ---- 4. TG/HDL-C ratio ----
# Computed in mg/dL, the basis for the usual clinical cut-points.
#
# NOTE: this is mathematically the same variable as AIP.
#   AIP = log10(TG_HDL_ratio) - 0.360
# A log transform plus a constant, so the two rank every
# participant identically. We report both because the proposal
# names both, but they are not independent findings.

data <- data %>%
  mutate(TG_HDL_ratio = TG / HDL)


# ---- 5. Non-HDL cholesterol ----
# non-HDL-C = TC - HDL-C, in mg/dL.
# All the cholesterol carried by atherogenic (apoB-containing)
# lipoproteins, combined.

data <- data %>%
  mutate(nonHDL_C = chol - HDL)


# ---- 6. LCI - Lipoprotein Combine Index ----
# LCI = (TC x TG x LDL-C) / HDL-C, all in mg/dL.
# Three concentrations multiplied together, so values are very
# large (order 10^5) and strongly right-skewed. Expect to need
# a log transform before regression - we check in session 5.

data <- data %>%
  mutate(LCI = (chol * TG * LDL) / HDL)

# ---- 7a. Convert gender to a factor ----
# Stored as labelled: 0 = male, 1 = female.
# Male is set as reference - an arbitrary but conventional
# choice. State it in your Methods.

data <- data %>%
  mutate(sex = relevel(as_factor(gender), ref = "male"))

table(data$sex)

# ---- 7. VAI - Visceral Adiposity Index ----
# Sex-specific. WC in cm, BMI in kg/m2, TG and HDL in mmol/L.
# Separate constants for men and women because the index was
# derived separately in each sex. A value near 1 is the reference.

data <- data %>%
  mutate(
    VAI = if_else(
      sex == "male",
      (Waist_circumference / (39.68 + 1.88 * BMI)) * (TG_mmol / 1.03) * (1.31 / HDL_mmol),
      (Waist_circumference / (36.58 + 1.89 * BMI)) * (TG_mmol / 0.81) * (1.52 / HDL_mmol)
    )
  )

# ---- 8. Checks ----

# (a) Confirm the AIP / TG-HDL relationship holds.
# This verifies BOTH formulas at once: if either had an error,
# the identity would break. Should print TRUE.
all(abs(data$AIP - (log10(data$TG_HDL_ratio) - 0.3599)) < 0.001, na.rm = TRUE)

# (b) Distributions of all five
data %>%
  select(AIP, TG_HDL_ratio, nonHDL_C, LCI, VAI) %>%
  summary()

# (c) Hand-check one VAI value
data %>%
  select(sex, Waist_circumference, BMI, TG_mmol, HDL_mmol, VAI) %>%
  head(3)

# ---- 9. Save ----
saveRDS(data, file.path(path_processed, "analysis.rds"))
