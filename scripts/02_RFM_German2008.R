# =============================================================================
# Script 02: RFM Classification — German et al. (2008)
# Arabidopsis thaliana — foundational PARE dataset
# SRA: SRP000036 | SRR000222 + SRR000223 (454 sequencing)
# =============================================================================
# Authors: [YOUR NAME]
# Date: June 2026
# R version: 4.3.1
# =============================================================================

WORK <- "/media/ingrid/D056C45556C43DCA/Bibliotecas/Documents/aPAPERS_IG/Review_degradome/RFM_validation"

# --- 1. Load CleaveLand4 results ---
s222 <- read.table(file.path(WORK, "cleaveland/german2008/SRR000222_results.txt"),
                   header = FALSE, sep = "\t", comment.char = "#",
                   fill = TRUE, quote = "")
s223 <- read.table(file.path(WORK, "cleaveland/german2008/SRR000223_results.txt"),
                   header = FALSE, sep = "\t", comment.char = "#",
                   fill = TRUE, quote = "")

# Extract relevant columns: miRNA, transcript, category, pval
s222_clean <- s222[s222[,1] != "SiteID" & ncol(s222) >= 16, c(2,3,15,16)]
s223_clean <- s223[s223[,1] != "SiteID" & ncol(s223) >= 16, c(2,3,15,16)]
names(s222_clean) <- names(s223_clean) <- c("miRNA","transcript","category","pval")

cat("=== German et al. (2008) — CleaveLand4 results ===\n")
cat("SRR000222 interactions:", nrow(s222_clean), "\n")
cat("SRR000223 interactions:", nrow(s223_clean), "\n")

# --- 2. Combine and deduplicate ---
all_pairs <- unique(rbind(s222_clean, s223_clean))
all_pairs$GeneID_short <- gsub("\\..*", "", all_pairs$transcript)
cat("Total unique pairs:", nrow(all_pairs), "\n")

# Cross-replicate reproducibility
s222_pairs <- paste(s222_clean$miRNA, s222_clean$transcript, sep = "__")
s223_pairs <- paste(s223_clean$miRNA, s223_clean$transcript, sep = "__")
shared <- intersect(s222_pairs, s223_pairs)
cat("Pairs reproduced in both libraries:", length(shared), "\n")
cat("Reproduced pairs:", shared, "\n")

# --- 3. Load RPKM from Yan 2024 RNA-seq (TAIR10 reference) ---
rpkm_all <- read.csv(file.path(WORK, "quantification/rpkm.csv"))

# Match genes
matches <- merge(all_pairs, rpkm_all,
                 by.x = "GeneID_short", by.y = "Geneid",
                 all.x = TRUE)
cat("\nGenes with RPKM available:", sum(!is.na(matches$RPKM_mean)), "\n")

# --- 4. RFM classification ---
# All pairs are Category 4 -> cleavage_score = 0
matches$cleavage_score <- 4 - as.numeric(as.character(matches$category))
matches$log_RPKM <- log2(matches$RPKM_mean + 1)

# Apply regression parameters from Yan 2024 baseline
rfm_yan <- read.csv(file.path(WORK, "quantification/RFM_classified.csv"))
fit_yan <- lm(cleavage_score ~ log_RPKM, data = rfm_yan)
sd_yan  <- sd(residuals(fit_yan))

matches$residual <- matches$cleavage_score -
                    predict(fit_yan, newdata = data.frame(log_RPKM = matches$log_RPKM))
matches$RFM <- ifelse(matches$residual > sd_yan,  "PTD",
               ifelse(matches$residual < -sd_yan, "TD", "CR"))

matches$miRNA_family <- gsub("ath-(miR[0-9]+).*", "\\1", matches$miRNA)

cat("\n=== RFM Classification German 2008 ===\n")
print(table(matches$RFM))
cat("\nNote: All pairs Category 4 (cleavage_score = 0); all classified CR.\n")
cat("Reflects insufficient library depth of 454 sequencing (<500K reads).\n")
cat("miR156 module recovered: miR156h-AT1G73830 (SRR222); miR156a-AT5G16000 (SRR223)\n")

# --- 5. Combined summary ---
cat("\n=== Combined Yan 2024 + German 2008 ===\n")
cat("Yan 2024 CR fraction: 33/40 = 82.5%\n")
cat("German 2008 CR fraction: 13/13 = 100%\n")
cat("Combined CR fraction: 46/53 = 86.8%\n")

# --- 6. Save ---
write.csv(matches, file.path(WORK, "RFM_repository/results/RFM_German2008_classified.csv"),
          row.names = FALSE)
cat("\nResults saved.\n")
