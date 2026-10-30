## Tumor cis-eQTL of the mediating SNVs
##
## expr   gene x patient matrix of TPM (tumor RNA-seq)
## d      data.frame, one row per patient: X (genotype), sex, age, TMB, karyotype,
##        PEER1-PEER15
## PEER factors: computed with the GTEx v8 eQTL pipeline (genes with TPM >= 0.1 in >= 20% of
## samples, 15 factors.

MIN_N  <- 5                                                   # patients required in every coded genotype group
covars <- c("sex", "age", "TMB", "karyotype", paste0("PEER", 1:15))
code   <- function(g, model) switch(model, add = g, dom = as.numeric(g >= 1), rec = as.numeric(g == 2))
enough <- function(g, model) {
  n <- c(sum(g == 0, na.rm = TRUE), sum(g == 1, na.rm = TRUE), sum(g == 2, na.rm = TRUE))
  switch(model, add = min(n) >= MIN_N, dom = min(n[1], n[2] + n[3]) >= MIN_N, rec = min(n[1] + n[2], n[3]) >= MIN_N)
}

## Target gene(s) of an SNV: its GTEx v8 whole-blood eGenes that are in our RNA-seq (rownames of expr); if none, the
## nearest protein-coding gene by TSS (GENCODE) within 1 Mb, otherwise the closest one.
## gtex   data.frame(chr, pos, gene_id) of significant whole-blood variant-gene pairs (GRCh38)
## genes  data.frame(gene_id, chr, tss) of protein-coding genes in our RNA-seq
target_genes <- function(chr, pos, gtex, genes) {
  hit <- intersect(gtex$gene_id[gtex$chr == chr & gtex$pos == pos], rownames(expr))
  if (length(hit)) return(hit)
  g <- genes[genes$chr == chr, ]; dist <- abs(g$tss - pos)
  if (any(dist <= 1e6)) { g <- g[dist <= 1e6, ]; dist <- dist[dist <= 1e6] }
  g$gene_id[which.min(dist)]
}

## eQTL test: log2(TPM + 1) ~ genotype + covariates, genotype g (0/1/2) coded in the SNV's GWAS model (tested only if
## every coded group has >= 5 patients); significant: P < 0.05
eqtl <- function(gene, g, model, d) {
  if (!enough(g, model)) return(c(NA, NA, NA))
  d$X <- code(g, model)
  d$E <- log2(as.numeric(expr[gene, ]) + 1)
  fit <- lm(as.formula(paste("E ~ X +", paste(covars, collapse = " + "))), data = d)
  coef(summary(fit))["X", c("Estimate", "Std. Error", "Pr(>|t|)")]
}

