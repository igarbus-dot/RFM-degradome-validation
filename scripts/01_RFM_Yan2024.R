# =============================================================================
# Script 01: RFM Classification — Yan et al. (2024)
# Arabidopsis thaliana Col-0 floral tissue
# BioProject PRJNA1092576 | SRR28470273-75 (PARE) | SRR28470277-78 (RNA-seq)
# =============================================================================
# Authors: [YOUR NAME]
# Date: June 2026
# R version: 4.3.1
# =============================================================================

library(ggplot2)

WORK <- "/media/ingrid/D056C45556C43DCA/Bibliotecas/Documents/aPAPERS_IG/Review_degradome/RFM_validation"

# --- 1. Load data ---
rfm <- read.csv(file.path(WORK, "quantification/RFM_classified.csv"))
cat("Pairs loaded:", nrow(rfm), "\n")

# --- 2. RFM classification ---
# Cleavage confidence score: CS = 4 - DegradomeCategory (range 0-4)
# Transcript abundance: T = log2(RPKM + 1)
# Linear regression baseline: CS ~ T
# Residuals > +1 SD = PTD; < -1 SD = TD; intermediate = CR

fit <- lm(cleavage_score ~ log_RPKM, data = rfm)
rfm$residual <- residuals(fit)
sd_r <- sd(rfm$residual)

rfm$RFM <- ifelse(rfm$residual > sd_r,  "PTD",
           ifelse(rfm$residual < -sd_r, "TD", "CR"))

cat("R² of baseline regression:", round(summary(fit)$r.squared, 3), "\n")
cat("\nRFM classification:\n")
print(table(rfm$RFM))

# --- 3. Statistical validation ---
# Fisher's exact test: PTD vs non-PTD x Category 0-1 vs 2-4
rfm$high_conf <- rfm$DegradomeCategory <= 1
rfm$is_PTD    <- rfm$RFM == "PTD"

cont_table <- table(PTD = rfm$is_PTD, HighConf = rfm$high_conf)
fisher_result <- fisher.test(cont_table)
cat("\nFisher's exact test:\n")
cat("  OR =", round(fisher_result$estimate, 3), "\n")
cat("  p  =", fisher_result$p.value, "\n")

# --- 4. Permutation test ---
# Are validated targets enriched in PTD zone?
validated_genes <- c("AT3G51140", "AT1G27360", "AT2G33810",
                     "AT5G60910", "AT1G69170")
set.seed(42)
n_perm <- 10000

rfm$is_validated <- rfm$GeneID %in% validated_genes
obs_ptd_frac <- mean(rfm$RFM[rfm$is_validated == TRUE] == "PTD", na.rm = TRUE)

perm_fracs <- replicate(n_perm, {
  idx <- sample(nrow(rfm), sum(rfm$is_validated, na.rm = TRUE))
  mean(rfm$RFM[idx] == "PTD")
})

perm_p <- mean(perm_fracs >= obs_ptd_frac)
cat("\nPermutation test (n =", n_perm, "):\n")
cat("  Observed PTD fraction in validated targets:", round(obs_ptd_frac, 3), "\n")
cat("  Permutation p-value:", perm_p, "\n")

# --- 5. Save results ---
write.csv(rfm, file.path(WORK, "RFM_repository/results/RFM_Yan2024_classified.csv"),
          row.names = FALSE)
cat("\nResults saved.\n")
