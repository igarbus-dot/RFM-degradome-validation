# Regulatory Flux Model (RFM) — Validation Scripts

R scripts for computational validation of the Regulatory Flux Model introduced in:

> Garbus Ingrid (2026). Degradome Sequencing in Plants: From Cleavage Mapping to Regulatory Flux Modeling. *[JOURNAL]*. DOI: [INSERT UPON PUBLICATION]

## Description

The Regulatory Flux Model (RFM) integrates degradome cleavage signal intensity with transcript abundance to classify miRNA–target interactions into three mechanistically distinct regulatory states:
- **PTD** (Post-transcriptional Dominance): high cleavage signal relative to transcript abundance
- **CR** (Compensatory Regulation): proportional cleavage and transcript levels
- **TD** (Transcriptional Dominance): high transcript abundance with low cleavage signal

## Repository structure
(base) ingrid@ingrid-Z270X-Gaming-5:~$ cat > "$WORK/RFM_repository/README.md" << 'EOF'
> # Regulatory Flux Model (RFM) — Validation Scripts
> 
> R scripts for computational validation of the Regulatory Flux Model introduced in:
> 
> > Garbus Ingrid  (2026). Degradome Sequencing in Plants: From Cleavage Mapping to Regulatory Flux Modeling. *[JOURNAL]*. DOI: [INSERT UPON PUBLICATION]
> 
> ## Description
> 
> The Regulatory Flux Model (RFM) integrates degradome cleavage signal intensity with transcript abundance to classify miRNA–target interactions into three mechanistically distinct regulatory states:
> - **PTD** (Post-transcriptional Dominance): high cleavage signal relative to transcript abundance
> - **CR** (Compensatory Regulation): proportional cleavage and transcript levels
> - **TD** (Transcriptional Dominance): high transcript abundance with low cleavage signal
> 
> ## Repository structure
# Session info
R --vanilla --quiet << 'EOF'
cat("=== R Session Info ===\n")
sessionInfo()
