# =============================================================
# 04_describe.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Table 1 - participant characteristics for the whole
#          cohort and by sex.
# Input:   data/processed/analysis.rds
# Output:  output/tables/table1_continuous.csv
#          output/tables/table1_categorical.csv
#
# REPORTING CONVENTIONS (following Salabat et al.)
#   Normally distributed continuous -> mean (SD)
#   Skewed continuous               -> median (IQR)
#   Categorical                     -> n (%)
#   Two-group comparison by sex     -> t-test / Mann-Whitney /
#                                      chi-square as appropriate
# =============================================================

source("R/00_setup.R")

data <- readRDS(file.path(path_processed, "analysis.rds"))

# ---- 1. Which variables are skewed? ----
# Decision rule, stated in the Methods so it is reproducible:
#   mean / median > 1.2  ->  report median (IQR)
# Shapiro-Wilk is reported separately for the record, but at
# n ~ 950 it rejects trivially small departures from normality
# and is not used to make this decision.

vars_continuous <- c(
  "age", "Height", "Weight", "BMI",
  "Waist_circumference", "Hip_circumference",
  "Sys_BP", "Dia_BP", "glu", "HbA1c", "VitD_Imputed",
  "TG", "chol", "HDL", "LDL",
  "AIP", "TG_HDL_ratio", "nonHDL_C", "LCI", "VAI"
)

skewness <- function(x) {
  mean((x - mean(x))^3) / sd(x)^3
}

skew_check <- do.call(rbind, lapply(vars_continuous, function(v) {
  x <- data[[v]][!is.na(data[[v]])]
  data.frame(
    Variable   = v,
    Mean       = round(mean(x), 2),
    Median     = round(median(x), 2),
    Skew_ratio = round(mean(x) / median(x), 3),
    Skewness   = round(skewness(x), 2),
    Shapiro_p  = signif(shapiro.test(x)$p.value, 3)
  )
}))

skew_check

# ---- 2. Identify skewed variables from the data ----
# Rule (stated in Methods): skewness coefficient > 1 is reported
# as median (IQR); otherwise mean (SD).
#
# The skewness coefficient is used in preference to the
# mean/median ratio, because the ratio is unreliable for
# variables centred near zero - AIP has a ratio of 1.17 despite
# being near-symmetric, simply because its mean is close to 0.

skew_threshold <- 1

vars_skewed <- skew_check$Variable[skew_check$Skewness > skew_threshold]
vars_normal <- setdiff(vars_continuous, vars_skewed)

cat("mean (SD)    :", paste(vars_normal, collapse = ", "), "\n")
cat("median (IQR) :", paste(vars_skewed, collapse = ", "), "\n")

# ---- 3. Helper: summarise one continuous variable ----
# Returns the appropriate summary for the whole cohort.
# The 'skewed' argument makes the formatting choice EXPLICIT
# rather than hidden inside the function - so your Methods can
# state exactly which variables were reported which way.

describe_continuous <- function(x, var_name, skewed = FALSE) {
  
  x_complete <- x[!is.na(x)]
  n_missing  <- sum(is.na(x))
  
  if (skewed) {
    q <- quantile(x_complete, c(0.25, 0.75))
    summary_text <- paste0(round(median(x_complete), 2),
                           " (", round(q[1], 2),
                           ", ", round(q[2], 2), ")")
    statistic <- "median (IQR)"
  } else {
    summary_text <- paste0(round(mean(x_complete), 2),
                           " (", round(sd(x_complete), 2), ")")
    statistic <- "mean (SD)"
  }
  
  data.frame(
    Variable  = var_name,
    Statistic = statistic,
    Total     = summary_text,
    n         = length(x_complete),
    Missing   = n_missing
  )
}

# ---- 4. Apply to every continuous variable ----

table1_total <- do.call(rbind, lapply(vars_continuous, function(v) {
  describe_continuous(data[[v]], v, skewed = v %in% vars_skewed)
}))

table1_total

# ---- 5. Helper: one complete Table 1 row ----
# Returns the total column, both sex columns, and a p-value.
# The 'skewed' argument routes BOTH the summary format and the
# choice of test, so the two can never disagree.

describe_row <- function(var_name, data, skewed = FALSE) {
  
  x      <- data[[var_name]]
  male   <- x[data$sex == "male"   & !is.na(x)]
  female <- x[data$sex == "female" & !is.na(x)]
  all_x  <- x[!is.na(x)]
  
  # Format one group's summary
  fmt <- function(v) {
    if (skewed) {
      q <- quantile(v, c(0.25, 0.75))
      paste0(round(median(v), 2), " (",
             round(q[1], 2), ", ", round(q[2], 2), ")")
    } else {
      paste0(round(mean(v), 2), " (", round(sd(v), 2), ")")
    }
  }
  
  # Choose the test to match the summary
  if (skewed) {
    p         <- wilcox.test(male, female)$p.value
    test_used <- "Mann-Whitney"
  } else {
    p         <- t.test(male, female)$p.value
    test_used <- "t-test"
  }
  
  data.frame(
    Variable  = var_name,
    Statistic = if (skewed) "median (IQR)" else "mean (SD)",
    Total     = fmt(all_x),
    Men       = fmt(male),
    Women     = fmt(female),
    P_value   = if (p < 0.001) "<0.001" else sprintf("%.3f", p),
    Test      = test_used,
    n         = length(all_x),
    Missing   = sum(is.na(x))
  )
}

# ---- 6. Build the continuous part of Table 1 ----

table1_continuous <- do.call(rbind, lapply(vars_continuous, function(v) {
  describe_row(v, data, skewed = v %in% vars_skewed)
}))

table1_continuous

# ---- 7. Helper: categorical variable ----
# Reports n (%) per level for the whole cohort and by sex,
# with a chi-square test of association with sex.
# One row per level, so the variable name appears once and
# levels are listed beneath it.

describe_categorical <- function(var_name, data) {
  
  x <- data[[var_name]]
  keep <- !is.na(x) & !is.na(data$sex)
  x    <- x[keep]
  sex  <- data$sex[keep]
  
  tab <- table(x, sex)
  
  # chisq.test on a 2-way table tests association between
  # the variable and sex.
  p <- chisq.test(tab)$p.value
  
  data.frame(
    Variable  = var_name,
    Level     = rownames(tab),
    Total     = paste0(rowSums(tab), " (",
                       round(100 * rowSums(tab) / sum(tab), 1), ")"),
    Men       = paste0(tab[, "male"], " (",
                       round(100 * tab[, "male"] / sum(tab[, "male"]), 1), ")"),
    Women     = paste0(tab[, "female"], " (",
                       round(100 * tab[, "female"] / sum(tab[, "female"]), 1), ")"),
    P_value   = if (p < 0.001) "<0.001" else sprintf("%.3f", p),
    Test      = "chi-square",
    row.names = NULL
  )
}

# ---- 8. Build the categorical part ----

vars_categorical <- c("bmi_cat", "vitd_status", "vitd_cat3",
                      "glycaemic", "AIP_band", "nonHDL_band")

table1_categorical <- do.call(rbind, lapply(vars_categorical, function(v) {
  describe_categorical(v, data)
}))

table1_categorical

# Sex itself, reported separately since it is the stratifier
table(data$sex)
round(100 * prop.table(table(data$sex)), 1)

# ---- 9. Save ----
write.csv(table1_continuous,
          file.path(path_tables, "table1_continuous.csv"), row.names = FALSE)
write.csv(table1_categorical,
          file.path(path_tables, "table1_categorical.csv"), row.names = FALSE)
write.csv(skew_check,
          file.path(path_tables, "normality_check.csv"), row.names = FALSE)