## Causal mediation: germline (X) -> evolutionary parameter (M) -> overall survival (Y)
##
## d   data.frame, one row per patient: OS_time, OS_event, X (genotype), the mediator, sex, age, TMB, karyotype
##
## Outcome model   Cox  Surv(OS_time, OS_event) ~ X + M + covariates      b = effect of M, c' = effect of X
## Mediator model  lm   M ~ X + covariates                                a = effect of X
##   NIE = exp(a b), NDE = exp(c');
## Estimated with the 'mediator' R package (GerkeLab, MIT licence; delta method)

library(survival)
library(dplyr)
library(tibble)
for (f in list.files("mediator/R", pattern = "\\.R$", full.names = TRUE)) source(f)   # mediator()

covars <- c("sex", "age", "TMB", "karyotype")

mediation <- function(d, M) {
  cox <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~ X +", M, "+", paste(covars, collapse = " + "))), data = d)
  med <- lm(as.formula(paste(M, "~ X +", paste(covars, collapse = " + "))), data = d)
  r <- as.data.frame(mediator(data = d, out.model = cox, med.model = med, treat = "X", mediator_name = M))
  colnames(r) <- c("Effect", "Estimate", "low", "high", "P")
  e <- function(k, col) r[r$Effect == k, col]
  data.frame(NIE = e("NIE", "Estimate"), NIE_low = e("NIE", "low"), NIE_high = e("NIE", "high"), NIE_P = e("NIE", "P"),
             NDE = e("NDE", "Estimate"), NDE_P = e("NDE", "P"),
             TE = e("Total Effect", "Estimate"), TE_P = e("Total Effect", "P"),
             PM = e("Proportion Mediated", "Estimate"))
}

mediating <- function(r) with(r, NIE_P < 0.05 & (NIE_low > 1 | NIE_high < 1) &
                                 ((NIE > 1 & TE > 1) | (NIE < 1 & TE < 1)) & PM > 0 & PM <= 1)
