# evalue.R -- E-values for the key HRs (VanderWeele-Ding). Quantifies how strong an
# unmeasured confounder (e.g. residential self-selection) would have to be, on the
# risk-ratio scale, to explain away each association. ~5.5% incidence -> rare=TRUE ok.
# install.packages("EValue") if needed
library(EValue)
show <- function(lab, est, lo, hi){
  e <- evalues.HR(est, lo, hi, rare = TRUE)
  ci <- e["E-values", c("lower","upper")]
  ci <- ci[!is.na(ci)][1]                     # for protective HR the CI E-value is the non-NA bound
  cat(sprintf("%-28s HR %.3f (%.3f-%.3f)  E-value point %.2f | CI %.2f\n",
              lab, est, lo, hi, e["E-values","point"], ci))
}
cat("=== E-values (point estimate | CI limit nearest null) ===\n")
show("Primary: cohesion->dep",        0.882, 0.857, 0.908)  # main cohort
show("Steps->incident depression",    0.753, 0.702, 0.807)  # device behaviour
show("Wellbeing(happy)->depression",  0.611, 0.581, 0.642)  # self-report affect
show("Cohesion->dep (SWB cohort)",    0.861, 0.813, 0.912)
