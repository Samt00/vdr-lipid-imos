# =============================================================
# 05_distributions.R
# Project: VDR polymorphisms and lipid indices (IMOS)
# Purpose: Visualise the distribution of the five lipid indices,
#          decide on transformations, and show each index across
#          genotypes.
# Input:   data/processed/analysis.rds
# Output:  output/figures/fig1_distributions.png
#          output/figures/fig2_by_genotype.png
#          data/processed/analysis_transformed.rds
# =============================================================

source("R/00_setup.R")
library(ggplot2)
library(tidyr)

data <- readRDS(file.path(path_processed, "analysis.rds"))


# ---- 1. Reshape the five indices to long format ----
# ggplot needs one row per observation-variable pair in order
# to facet. pivot_longer() stacks the five columns into two:
# one holding the index name, one holding the value.

indices <- c("AIP", "TG_HDL_ratio", "nonHDL_C", "LCI", "VAI")

long_data <- data %>%
  select(all_of(indices)) %>%
  pivot_longer(cols = everything(),
               names_to  = "index",
               values_to = "value") %>%
  filter(!is.na(value))

head(long_data)
dim(long_data)

# ---- 2. Histograms with density overlay ----
# facet_wrap() draws one panel per index.
# scales = "free" lets each panel set its own axes - essential
# here, since AIP spans about 2 units and LCI spans 700,000.

fig_hist <- ggplot(long_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey70", colour = "white",
                 linewidth = 0.2) +
  geom_density(colour = "steelblue", linewidth = 0.7) +
  facet_wrap(~ index, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Density") +
  theme_minimal(base_size = 11)

fig_hist

# ---- 3. Q-Q plots ----
# A Q-Q plot compares your data's quantiles against those of a
# perfect normal distribution. Points on the diagonal = normal.
# An upward curve at the right end = right skew (the largest
# values are larger than a normal distribution would produce).

fig_qq <- ggplot(long_data, aes(sample = value)) +
  stat_qq(size = 0.5, alpha = 0.4) +
  stat_qq_line(colour = "firebrick", linewidth = 0.6) +
  facet_wrap(~ index, scales = "free", ncol = 3) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal(base_size = 11)

fig_qq

# ---- 4. Boxplots ----
fig_box <- ggplot(long_data, aes(y = value)) +
  geom_boxplot(fill = "grey85", outlier.size = 0.6,
               outlier.alpha = 0.4) +
  facet_wrap(~ index, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank())

fig_box
# The bare minimum
# Set bins and colours
# Add a density curve for shape
ggplot(data, aes(x = AIP)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey70", colour = "white") +
  geom_density(colour = "steelblue", linewidth = 0.8) +
  labs(title = "AIP", x = NULL, y = "Density") +
  theme_minimal()

# ---- 5. Log transformation ----
# LCI, VAI and the TG/HDL ratio show marked right skew on both
# the skewness coefficient (2.60, 2.63, 3.12) and the Q-Q plots.
# All three are strictly positive products or ratios, so a
# natural log is appropriate.
#
# NOT transformed:
#   AIP      - skewness 0.24; it is already a log10 transform
#   nonHDL_C - skewness 0.71, approximately symmetric
#
# NOTE ON TG/HDL: log(TG_HDL_ratio) is algebraically identical
# to AIP up to a constant (AIP = log10(ratio) - 0.360), so the
# logged ratio and AIP are the same variable. We transform it
# for the figure, but AIP is what enters the models.

data <- data %>%
  mutate(
    log_LCI = log(LCI),
    log_VAI = log(VAI)
  )

# Check the transformation worked
skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3

c(LCI      = skewness(data$LCI[!is.na(data$LCI)]),
  log_LCI  = skewness(data$log_LCI[!is.na(data$log_LCI)]),
  VAI      = skewness(data$VAI[!is.na(data$VAI)]),
  log_VAI  = skewness(data$log_VAI[!is.na(data$log_VAI)]))

# ---- 6. Q-Q plots after transformation ----

long_transformed <- data %>%
  select(AIP, nonHDL_C, log_LCI, log_VAI) %>%
  pivot_longer(cols = everything(),
               names_to = "index", values_to = "value") %>%
  filter(!is.na(value))

fig_qq_after <- ggplot(long_transformed, aes(sample = value)) +
  stat_qq(size = 0.5, alpha = 0.4) +
  stat_qq_line(colour = "firebrick", linewidth = 0.6) +
  facet_wrap(~ index, scales = "free", ncol = 2) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal(base_size = 11)

fig_qq_after

# ---- 7. Save ----
saveRDS(data, file.path(path_processed, "analysis_transformed.rds"))

# ---- 8. Reshape for genotype boxplots ----
# We need a long format with THREE columns: index name, value,
# and genotype. pivot_longer is applied twice: once to stack
# the outcomes, once to stack the five SNPs.

outcomes <- c("AIP", "nonHDL_C", "log_LCI", "log_VAI")
genotypes <- c("apa_geno", "taq_geno", "ecorv_geno",
               "fok_geno", "bsm_geno")

geno_long <- data %>%
  select(all_of(c(outcomes, genotypes))) %>%
  pivot_longer(cols = all_of(outcomes),
               names_to = "index", values_to = "value") %>%
  pivot_longer(cols = all_of(genotypes),
               names_to = "snp", values_to = "genotype") %>%
  filter(!is.na(value), !is.na(genotype))

head(geno_long)
dim(geno_long)

# ---- 9. Clean up the SNP labels ----
# Replace variable names with the enzyme names a reader expects.

geno_long <- geno_long %>%
  mutate(
    snp = recode(snp,
                 apa_geno   = "ApaI",
                 taq_geno   = "TaqI",
                 ecorv_geno = "EcoRV",
                 fok_geno   = "FokI",
                 bsm_geno   = "BsmI"),
    snp = factor(snp, levels = c("ApaI", "TaqI", "EcoRV",
                                 "FokI", "BsmI"))
  )

# ---- 10. Boxplots for one SNP across all four outcomes ----

plot_snp <- function(snp_name) {
  geno_long %>%
    filter(snp == snp_name) %>%
    ggplot(aes(x = genotype, y = value)) +
    geom_boxplot(fill = "grey85", outlier.size = 0.5,
                 outlier.alpha = 0.3, width = 0.6) +
    facet_wrap(~ index, scales = "free_y", nrow = 1) +
    labs(title = snp_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major.x = element_blank())
}

plot_snp("ApaI")
plot_snp("TaqI")
plot_snp("EcoRV")
plot_snp("FokI")
plot_snp("BsmI")

# ---- 11. Save figures ----

ggsave(file.path(path_figures, "fig1_histograms.png"),
       fig_hist, width = 10, height = 6, dpi = 300)
ggsave(file.path(path_figures, "fig1_qq.png"),
       fig_qq, width = 10, height = 6, dpi = 300)
ggsave(file.path(path_figures, "fig1_qq_transformed.png"),
       fig_qq_after, width = 8, height = 7, dpi = 300)
ggsave(file.path(path_figures, "fig1_boxplots.png"),
       fig_box, width = 10, height = 4, dpi = 300)

for (s in levels(geno_long$snp)) {
  ggsave(file.path(path_figures, paste0("fig2_", s, ".png")),
         plot_snp(s), width = 10, height = 3.2, dpi = 300)
}

# ---- 12. Save the transformed dataset ----
saveRDS(data, file.path(path_processed, "analysis_transformed.rds"))