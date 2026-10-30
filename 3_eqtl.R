## Tumor cis-eQTL of the mediating SNVs
##
## expr   gene x patient matrix of TPM (tumor RNA-seq)
## d      data.frame, one row per patient: X (genotype), sex, age, TMB, karyotype,
##        PEER1-PEER15
## PEER factors: computed with the GTEx v8 eQTL pipeline (genes with TPM >= 0.1 in >= 20% of
## samples, 15 factors.

covars <- c("sex", "age", "TMB", "karyotype", paste0("PEER", 1:15))

## Target gene(s) of an SNV: its GTEx v8 whole-blood eGenes that are expressed in our RNA-seq; if none, the nearest
## protein-coding gene by TSS (GENCODE) within 1 Mb, otherwise the closest one.
## gtex   data.frame(chr, pos, gene_id) of significant whole-blood variant-gene pairs (GRCh38)
## genes  data.frame(gene_id, chr, tss) of protein-coding genes expressed in our RNA-seq
target_genes <- function(chr, pos, gtex, genes) {
  hit <- intersect(gtex$gene_id[gtex$chr == chr & gtex$pos == pos], genes$gene_id)
  if (length(hit)) return(hit)
  g <- genes[genes$chr == chr, ]; dist <- abs(g$tss - pos)
  if (any(dist <= 1e6)) { g <- g[dist <= 1e6, ]; dist <- dist[dist <= 1e6] }
  g$gene_id[which.min(dist)]
}

## eQTL test: log2(TPM + 1) ~ genotype + covariates; significant: P < 0.05
eqtl <- function(gene, d) {
  d$E <- log2(as.numeric(expr[gene, ]) + 1)
  fit <- lm(as.formula(paste("E ~ X +", paste(covars, collapse = " + "))), data = d)
  coef(summary(fit))["X", c("Estimate", "Std. Error", "Pr(>|t|)")]
}

