## Genome-wide association of tumor evolutionary parameters with germline SNVs
##
## geno   SNV x patient matrix of genotype dosages (0, 1, 2) 
## snv    data.frame(id, chr, pos) for the rows of geno 
## pheno  data.frame, one row per patient  mu, tau, sex, age, TMB, karyotype (hyperdiploid = 1),
##        PC1-PC10 (germline principal components)

library(data.table)
library(qvalue)

MIN_N  <- 5                                           
covars <- c("sex", "age", "TMB", "karyotype", paste0("PC", 1:10))

## 1. Linear regression per SNV under additive, dominant and recessive coding
code   <- function(g, model) switch(model, add = g, dom = as.numeric(g >= 1), rec = as.numeric(g == 2))
enough <- function(g, model) {
  n <- c(sum(g == 0), sum(g == 1), sum(g == 2))
  switch(model, add = min(n) >= MIN_N, dom = min(n[1], n[2] + n[3]) >= MIN_N, rec = min(n[1] + n[2], n[3]) >= MIN_N)
}
test_snv <- function(g, y, model) {
  ok <- !is.na(g) & !is.na(y)
  if (!enough(g[ok], model)) return(c(NA, NA))
  d <- data.frame(y = y, G = code(g, model), pheno[, covars])[ok, ]
  coef(summary(lm(y ~ ., data = d)))["G", c("Estimate", "Pr(>|t|)")]
}
gwas <- rbindlist(lapply(c("mu", "tau"), function(param) rbindlist(lapply(c("add", "dom", "rec"), function(m) {
  r <- t(apply(geno, 1, test_snv, y = pheno[[param]], model = m))
  data.table(snv, parameter = param, model = m, beta = r[, 1], P = r[, 2])
}))))
lambda <- gwas[!is.na(P), .(lambda_GC = median(qchisq(P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)), by = .(parameter, model)]

## 2. LD blocks with TreeMap (Liu et al., Bioinformatics 2021)
source("TreeMap/R/my.functions.pre.processing.R")           
unit_of <- rbindlist(lapply(unique(snv$chr), function(ch) {
  i   <- which(snv$chr == ch); i <- i[order(snv$pos[i])]
  r2  <- suppressWarnings(cor(t(geno[i, ]), use = "pairwise.complete.obs"))^2; r2[is.na(r2)] <- 0
  blk <- find.blocks.graph(r2, cutoffs = 0.9)               
  grp <- integer(length(i)); for (k in seq_len(nrow(blk))) grp[blk$start[k]:blk$end[k]] <- k
  data.table(id = snv$id[i], unit = ifelse(grp > 0, paste0(ch, "_block", grp), snv$id[i]))
}))

## 3.  significant: q < 0.1
units <- merge(gwas[!is.na(P)], unit_of, by = "id")[,
  .(chr = chr[1], start = min(pos), end = max(pos), lead = id[which.min(P)], P = min(P)), by = .(parameter, model, unit)]
units[, q := qvalue(P)$qvalues, by = .(parameter, model)]
sig <- units[q < 0.1, .(chr = chr[1], start = min(start), end = max(end)), by = .(parameter, unit)]

## 4. Loci: significant units on the same chromosome merged when < 250 kb apart (chained; FUMA rule),
##    mu and tau pooled (GCE loci)
merge_loci <- function(u, gap = 2.5e5) {
  u <- u[order(chr, start)]; id <- integer(nrow(u)); k <- 0L; edge <- -Inf; prev <- ""
  for (i in seq_len(nrow(u))) {
    if (u$chr[i] != prev || u$start[i] - edge >= gap) { k <- k + 1L; edge <- -Inf }
    id[i] <- k; edge <- max(edge, u$end[i]); prev <- u$chr[i]
  }
  u[, locus := id][]
}
loci <- merge_loci(sig)
