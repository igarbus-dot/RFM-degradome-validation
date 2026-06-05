# =============================================================================
# Script 03: Predictive RFM — Eragrostis curvula
# Apomictic (Tanganyika) vs Sexual (OTA-S) genotypes
# sRNA-seq: BioProject PRJNA378998 | RNA-seq: Garbus et al. 2017
# =============================================================================
# Authors: [YOUR NAME]
# Date: June 2026
# R version: 4.3.1
# Dependencies: edgeR v4.10.1
# =============================================================================

library(edgeR)

WORK <- "/media/ingrid/D056C45556C43DCA/Bibliotecas/Documents/aPAPERS_IG/Review_degradome/RFM_validation"
EC   <- file.path(WORK, "ecurvula")

# --- 1. Load miRNA counts and RPKM ---
counts  <- read.csv(file.path(EC, "miRNA_abundance_by_family.csv"), row.names = 1)
rpkm    <- read.csv(file.path(EC, "transcript_RPKM.csv"))

cat("miRNA families loaded:", nrow(counts), "\n")
cat("Transcripts with RPKM:", nrow(rpkm), "\n")

# --- 2. EdgeR differential expression ---
# Genotypes: T3P1, T3P2 = Tanganyika (apomictic)
#            O2P1, O2P2 = OTA-S (sexual)
group <- factor(c("Tanganyika","Tanganyika","OTA-S","OTA-S"))
dge   <- DGEList(counts = counts, group = group)
dge   <- calcNormFactors(dge)
dge   <- estimateCommonDisp(dge)
et    <- exactTest(dge, pair = c("OTA-S","Tanganyika"))
res   <- topTags(et, n = Inf)$table
res$miRNA_family <- rownames(res)

cat("\n=== Differential expression (edgeR) ===\n")
cat("Criteria: |logFC| > 2, p < 0.01\n")
sig <- res[abs(res$logFC) > 2 & res$PValue < 0.01, ]
print(sig[, c("logFC","PValue","FDR")])

# Confirm miR156 and miR8175 not differential (concordance with Garbus 2019)
cat("\nmiR156  — logFC:", round(res["miR156","logFC"],2),
    "FDR:", round(res["miR156","FDR"],3), "\n")
cat("miR8175 — logFC:", round(res["miR8175","logFC"],2),
    "FDR:", round(res["miR8175","FDR"],3), "\n")

# Save edgeR results
write.csv(res, file.path(WORK, "RFM_repository/results/edgeR_miRNA_Ecurvula.csv"))

# --- 3. Predictive RFM classification ---
# Proxy: miRNA abundance (reads) as cleavage signal
# Target abundance: RPKM from 454 RNA-seq
# Thresholds: tertile-based (33rd/67th percentiles)

pairs <- data.frame(
  miRNA      = c("miR156","miR156","miR156","miR8175"),
  isotig     = c("isotig37981","isotig18002","isotig40670","isotig39554"),
  target_ann = c("SPL protein","SPL protein","SPL protein","TE retrotrans_gag")
)

# Extract miRNA counts per genotype
get_mirna_mean <- function(fam, geno_cols) {
  if (fam %in% rownames(counts))
    mean(as.numeric(counts[fam, geno_cols]))
  else NA
}

tan_cols <- c("T3P1","T3P2")
ota_cols <- c("O2P1","O2P2")

pairs$miRNA_Tanganyika <- sapply(pairs$miRNA, get_mirna_mean, tan_cols)
pairs$miRNA_OTA_S      <- sapply(pairs$miRNA, get_mirna_mean, ota_cols)

# Extract transcript RPKM per genotype
get_rpkm <- function(isotig, geno) {
  row <- rpkm[rpkm$transcript == isotig, ]
  if (nrow(row) == 0) return(NA)
  if (geno == "Tanganyika") return(row$RPKM_Tanganyika)
  return(row$RPKM_OTA_S)
}

pairs$RPKM_Tanganyika <- sapply(pairs$isotig, get_rpkm, "Tanganyika")
pairs$RPKM_OTA_S      <- sapply(pairs$isotig, get_rpkm, "OTA-S")

# Compute RFM ratio: log2(miRNA+1) / log2(RPKM+0.1)
rfm_ratio <- function(mirna_all, rpkm_all) {
  ratios <- log2(mirna_all + 1) / log2(rpkm_all + 0.1)
  t33 <- quantile(ratios, 0.33, na.rm = TRUE)
  t67 <- quantile(ratios, 0.67, na.rm = TRUE)
  ifelse(ratios > t67, "PTD",
  ifelse(ratios < t33, "TD", "CR"))
}

all_mirna <- c(pairs$miRNA_Tanganyika, pairs$miRNA_OTA_S)
all_rpkm  <- c(pairs$RPKM_Tanganyika, pairs$RPKM_OTA_S)

states <- rfm_ratio(all_mirna, all_rpkm)
pairs$RFM_Tanganyika <- states[1:nrow(pairs)]
pairs$RFM_OTA_S      <- states[(nrow(pairs)+1):(2*nrow(pairs))]

cat("\n=== Predicted RFM states ===\n")
print(pairs[, c("miRNA","isotig","target_ann",
                "miRNA_OTA_S","RPKM_OTA_S","RFM_OTA_S",
                "miRNA_Tanganyika","RPKM_Tanganyika","RFM_Tanganyika")])

# --- 4. Save ---
write.csv(pairs, file.path(WORK, "RFM_repository/results/RFM_predicted_Ecurvula.csv"),
          row.names = FALSE)
cat("\nResults saved.\n")
