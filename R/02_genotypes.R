# =============================================================
# 02_genotypes.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Convert Stata's labelled genotype columns into R
#          factors with a sensible reference level, and build
#          carrier-state variables.
# Input:   data/processed/clean.rds
# Output:  data/processed/geno.rds
#
# SOURCE COLUMNS
#   We build from the *_ordinal columns, NOT the raw Apa/Taq/etc.
#   columns. The raw columns code the heterozygote LAST
#   (1 = AA, 2 = CC, 3 = AC), which is nominal, not ordered.
#   The *_ordinal columns are ordered correctly:
#   1 = homozygous major, 2 = heterozygous, 3 = homozygous minor.
# =============================================================

source("R/00_setup.R")

data <- readRDS(file.path(path_processed, "clean.rds"))

# ---- 1. See what a labelled column actually contains ----
head(data$apa_oerdinal)          # the raw numbers plus attached labels
attr(data$apa_oerdinal, "labels")   # the dictionary itself

# ---- 2. A helper function ----

set_reference_most_common <- function(f) {
  counts      <- table(f)                                   # count each level
  most_common <- names(sort(counts, decreasing = TRUE))[1]  # name of the biggest
  relevel(f, ref = most_common)                             # move it to first
}

# ---- 3. Build the genotype factors ----
# as_factor() (from haven) converts a labelled number into a
# factor using its Stata value labels, keeping their order.
# So levels come out as: homozygous major, heterozygous,
# homozygous minor - the biologically sensible order.

data <- data %>%
  mutate(
    apa_geno   = set_reference_most_common(as_factor(apa_oerdinal)),
    taq_geno   = set_reference_most_common(as_factor(Taq_ordinal)),
    ecorv_geno = set_reference_most_common(as_factor(EcoRV_ordinal)),
    fok_geno   = set_reference_most_common(as_factor(Fok_ordinal)),
    bsm_geno   = set_reference_most_common(as_factor(Bsm_ordinal))
  )

# ---- 4. Check the conversion ----
# The first level listed is the reference.

table(data$apa_geno)
table(data$taq_geno)
table(data$ecorv_geno)
table(data$fok_geno)
table(data$bsm_geno)

# ---- 5. Cross-check against the original ----
# Every row should map to exactly one column. Any off-diagonal
# value means the conversion lost or scrambled something.

table(as_factor(data$apa_oerdinal), data$apa_geno)

# ---- 6. Allele counts (0, 1, 2) ----
# as.numeric() on a labelled column returns 1, 2, 3.
# Subtracting 1 gives the number of copies of the counted allele:
#   0 = homozygous major, 1 = heterozygous, 2 = homozygous minor
#
# This is an intermediate variable used to build the carrier
# states below. It can also serve as an additive (per-allele)
# predictor if you later decide to add that model.

data <- data %>%
  mutate(
    apa_n   = as.numeric(apa_oerdinal)  - 1,
    taq_n   = as.numeric(Taq_ordinal)   - 1,
    ecorv_n = as.numeric(EcoRV_ordinal) - 1,
    fok_n   = as.numeric(Fok_ordinal)   - 1,
    bsm_n   = as.numeric(Bsm_ordinal)   - 1
  )

# Verify against the variable Stata already built for ApaI.
# If TRUE, our recoding logic is independently confirmed - and
# the same logic applies to the other four.
all(data$apa_n == data$apa2_ordinal, na.rm = TRUE)

# ---- 7. Carrier states ----
# Two contrasts per SNP, following Salabat et al.:
#
#   _dom : carries at least one minor allele (1 or 2 copies)
#          vs none. Tests "does carrying the allele matter at all?"
#
#   _rec : carries two minor alleles vs everyone else.
#          Tests "does the effect need both copies?"
#
# These are the two carrier states in the published paper, and
# they correspond to the dominant and recessive genetic models.
#
# We derive them from the allele counts rather than using the
# existing Stata columns (Apa_cat2, Apa_cat3, ...), because
# those are coded in inconsistent directions across SNPs -
# ApaI runs opposite to the other four. Deriving them ourselves
# guarantees all five are coded identically.

make_dominant <- function(n_copies) {
  f <- factor(if_else(n_copies >= 1, "carrier", "non-carrier"))
  relevel(f, ref = "non-carrier")
}

make_recessive <- function(n_copies) {
  f <- factor(if_else(n_copies == 2, "homozygous", "other"))
  relevel(f, ref = "other")
}

data <- data %>%
  mutate(
    apa_dom   = make_dominant(apa_n),   apa_rec   = make_recessive(apa_n),
    taq_dom   = make_dominant(taq_n),   taq_rec   = make_recessive(taq_n),
    ecorv_dom = make_dominant(ecorv_n), ecorv_rec = make_recessive(ecorv_n),
    fok_dom   = make_dominant(fok_n),   fok_rec   = make_recessive(fok_n),
    bsm_dom   = make_dominant(bsm_n),   bsm_rec   = make_recessive(bsm_n)
  )

# ---- 7. Carrier states ----
# Two contrasts per SNP, following Salabat et al.:
#
#   _dom : carries at least one minor allele (1 or 2 copies)
#          vs none. Tests "does carrying the allele matter at all?"
#
#   _rec : carries two minor alleles vs everyone else.
#          Tests "does the effect need both copies?"
#
# These are the two carrier states in the published paper, and
# they correspond to the dominant and recessive genetic models.
#
# We derive them from the allele counts rather than using the
# existing Stata columns (Apa_cat2, Apa_cat3, ...), because
# those are coded in inconsistent directions across SNPs -
# ApaI runs opposite to the other four. Deriving them ourselves
# guarantees all five are coded identically.

make_dominant <- function(n_copies) {
  f <- factor(if_else(n_copies >= 1, "carrier", "non-carrier"))
  relevel(f, ref = "non-carrier")
}

make_recessive <- function(n_copies) {
  f <- factor(if_else(n_copies == 2, "homozygous", "other"))
  relevel(f, ref = "other")
}

data <- data %>%
  mutate(
    apa_dom   = make_dominant(apa_n),   apa_rec   = make_recessive(apa_n),
    taq_dom   = make_dominant(taq_n),   taq_rec   = make_recessive(taq_n),
    ecorv_dom = make_dominant(ecorv_n), ecorv_rec = make_recessive(ecorv_n),
    fok_dom   = make_dominant(fok_n),   fok_rec   = make_recessive(fok_n),
    bsm_dom   = make_dominant(bsm_n),   bsm_rec   = make_recessive(bsm_n)
  )

# ---- 8. Check ----
# Each genotype should fall into exactly one carrier category.
table(data$apa_geno, data$apa_dom)
table(data$apa_geno, data$apa_rec)

# Counts for all five
sapply(data[c("apa_dom","taq_dom","ecorv_dom","fok_dom","bsm_dom")], table)
sapply(data[c("apa_rec","taq_rec","ecorv_rec","fok_rec","bsm_rec")], table)

# ---- 9. Save ----
saveRDS(data, file.path(path_processed, "geno.rds"))