# Analysis code

Code for the main analyses of the study.

| Script | Analysis |
|---|---|
| `1_genome_wide_association.R` | Genome-wide association of the evolutionary parameters (mu, tau) under additive, dominant and recessive models|
| `2_mediation.R` | Causal mediation: germline SNV -> evolutionary parameter -> overall survival |
| `3_eqtl.R` | Tumor cis-eQTL of the mediating SNVs |

Each script defines the models used; input data are described at the top of each script.

**Software:** R (4.4), packages `data.table`, `qvalue`, `survival`; PLINK 1.9 for genotype QC; TreeMap (Liu et al., Bioinformatics 2021) for LD blocks.

**mediator:** `2_mediation.R` loads the R package `mediator` (GerkeLab, MIT licence) from (https://github.com/GerkeLab/mediator).

**Data:** MMRF CoMMpass data.
