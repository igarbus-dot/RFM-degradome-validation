# =============================================================================
# Script 04: Figure generation — RFM validation and E. curvula case study
# Produces Figures 7 and 8 of the manuscript
# =============================================================================
# Authors: [YOUR NAME]
# Date: June 2026
# R version: 4.3.1
# Dependencies: ggplot2
# =============================================================================

library(ggplot2)

WORK <- "/media/ingrid/D056C45556C43DCA/Bibliotecas/Documents/aPAPERS_IG/Review_degradome/RFM_validation"

# ---- FIGURE 7A: RFM scatter plot — Yan et al. 2024 ----
rfm <- read.csv(file.path(WORK, "RFM_repository/results/RFM_Yan2024_classified.csv"))

rfm$RFM <- factor(rfm$RFM, levels = c("PTD","CR","TD"))
colors   <- c("PTD" = "#1B6CA8", "CR" = "#27AE60", "TD" = "#E67E22")

p7a <- ggplot(rfm, aes(x = log_RPKM, y = cleavage_score, color = RFM)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40",
              linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = colors) +
  annotate("text", x = 8.5, y = 4.2,
           label = "miR156d-AT3G51140\n(Cat. 0, 3 replicates)",
           size = 3, hjust = 0, color = "#1B6CA8") +
  labs(x = expression(log[2](RPKM)),
       y = "Cleavage confidence score",
       color = "RFM zone",
       title = "A. RFM Classification — Yan et al. (2024)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

# ---- FIGURE 7B: Category 0-1 enrichment per RFM zone ----
rfm$high_conf <- rfm$DegradomeCategory <= 1
zone_pct <- aggregate(high_conf ~ RFM, data = rfm,
                      FUN = function(x) round(mean(x) * 100, 1))
names(zone_pct)[2] <- "pct_high_conf"

p7b <- ggplot(zone_pct, aes(x = RFM, y = pct_high_conf, fill = RFM)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = colors) +
  annotate("text", x = "PTD", y = 105, label = "***", size = 6) +
  labs(x = "RFM zone",
       y = "% Category 0-1 interactions",
       title = "B. High-confidence cleavage by RFM zone") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

# Save Figure 7
png(file.path(WORK, "RFM_repository/figures/Figure7_RFM_Arabidopsis.png"),
    width = 2400, height = 1200, res = 300)
gridExtra::grid.arrange(p7a, p7b, ncol = 2)
dev.off()
cat("Figure 7 saved.\n")

# ---- FIGURE 8: E. curvula RFM state transitions ----
ec <- read.csv(file.path(WORK, "RFM_repository/results/RFM_predicted_Ecurvula.csv"))

# Panel A: miR156 -> isotig37981 (SPL)
spl <- ec[ec$isotig == "isotig37981", ]
df_spl <- data.frame(
  Genotype = c("OTA-S\n(sexual)", "Tanganyika\n(apomictic)"),
  RPKM     = c(spl$RPKM_OTA_S, spl$RPKM_Tanganyika),
  RFM      = c(spl$RFM_OTA_S, spl$RFM_Tanganyika),
  miRNA    = c(spl$miRNA_OTA_S, spl$miRNA_Tanganyika)
)
df_spl$fill_col <- ifelse(df_spl$RFM == "PTD", "#1B6CA8", "#E67E22")

p8a <- ggplot(df_spl, aes(x = Genotype, y = RPKM, fill = RFM)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = paste0("miR156\n", round(miRNA,1), " reads")),
            vjust = 0.5, color = "white", size = 3.2, fontface = "bold",
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("PTD" = "#1B6CA8", "TD" = "#E67E22")) +
  annotate("segment", x = 1.3, xend = 1.7, y = 780, yend = 780,
           arrow = arrow(length = unit(0.2,"cm"), ends = "both")) +
  annotate("text", x = 1.5, y = 820, label = "TD \u2192 PTD",
           size = 3.5, fontface = "bold") +
  labs(x = "", y = "SPL transcript (RPKM)",
       fill = "RFM state",
       title = "A. miR156 \u2192 isotig37981 (SPL)") +
  theme_classic(base_size = 12)

# Panel B: miR8175 -> isotig39554 (TE)
te <- ec[ec$isotig == "isotig39554", ]
df_te <- data.frame(
  Genotype = c("OTA-S\n(sexual)", "Tanganyika\n(apomictic)"),
  RPKM     = c(te$RPKM_OTA_S, te$RPKM_Tanganyika),
  RFM      = c(te$RFM_OTA_S, te$RFM_Tanganyika),
  miRNA    = c(te$miRNA_OTA_S, te$miRNA_Tanganyika)
)

p8b <- ggplot(df_te, aes(x = Genotype, y = RPKM, fill = RFM)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = paste0("miR8175\n", round(miRNA,1), " reads")),
            vjust = 0.5, color = "white", size = 3.2, fontface = "bold",
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("PTD" = "#1B6CA8", "TD" = "#E67E22")) +
  annotate("segment", x = 1.3, xend = 1.7, y = 720, yend = 720,
           arrow = arrow(length = unit(0.2,"cm"), ends = "both")) +
  annotate("text", x = 1.5, y = 760, label = "PTD \u2192 TD",
           size = 3.5, fontface = "bold") +
  labs(x = "", y = "TE transcript (RPKM)",
       fill = "RFM state",
       title = "B. miR8175 \u2192 isotig39554 (TE)") +
  theme_classic(base_size = 12)

# Save Figure 8
png(file.path(WORK, "RFM_repository/figures/Figure8_RFM_Ecurvula.png"),
    width = 2400, height = 1200, res = 300)
gridExtra::grid.arrange(p8a, p8b, ncol = 2)
dev.off()
cat("Figure 8 saved.\n")
