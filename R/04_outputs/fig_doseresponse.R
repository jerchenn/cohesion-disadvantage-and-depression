# fig_doseresponse.R -- Dose-response of incident depression across the cohesion distribution,
# modeled with a restricted cubic (natural) spline; hazard ratio relative to median cohesion, with
# a 95% CI band. Shows the graded protective effect (and lets readers see any nonlinearity).
# Aggregate curve coordinates only.
#
# INPUT: cox_dat.rds       OUTPUT: fig_doseresponse_data.csv, figure_doseresponse.{png,pdf}
# RUN AFTER build_cox.R.   Needs: survival, splines, ggplot2.

library(survival); library(splines); library(ggplot2)

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
d <- collapse(readRDS("cox_dat3.rds")); d <- d[!is.na(d$z_cohesion), ]

## natural spline (3 df) on cohesion; keep the basis object to rebuild it at the grid + reference
bas <- ns(d$z_cohesion, df=3)
fit <- coxph(Surv(time_days,event) ~ bas + age + sex_c + race_c + ethn_c +
               income_m + educ_m + log_util, d)
sp <- grep("^bas", names(coef(fit)))          # spline coefficient positions
b  <- coef(fit)[sp]; V <- vcov(fit)[sp, sp]

zg   <- seq(quantile(d$z_cohesion,.02), quantile(d$z_cohesion,.98), length.out=120)
ref  <- median(d$z_cohesion)
D    <- sweep(predict(bas, zg), 2, as.numeric(predict(bas, ref)))   # basis contrast vs median
lp   <- as.vector(D %*% b)
se   <- sqrt(rowSums((D %*% V) * D))
out  <- data.frame(z_cohesion=zg, hr=exp(lp), lo=exp(lp-1.96*se), hi=exp(lp+1.96*se))
write.csv(out, "fig_doseresponse_data.csv", row.names=FALSE)

p <- ggplot(out, aes(z_cohesion, hr)) +
  geom_hline(yintercept=1, linetype=2, colour="grey55") +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15) +
  geom_line(linewidth=0.8) +
  scale_y_log10() +
  labs(x="Perceived neighborhood cohesion (SD from mean)",
       y="Hazard ratio for incident depression\n(relative to median cohesion)",
       title="Dose-response: cohesion and incident depression") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title.position="plot",
        plot.title=element_text(face="bold", size=11))
ggsave("figure_doseresponse.png", p, width=7.0, height=4.6, dpi=300)
ggsave("figure_doseresponse.pdf", p, width=7.0, height=4.6)
cat("Dose-response spline written (HR vs median cohesion). Nonlinearity is testable via\n",
    "anova(coxph(~z_cohesion+covs), coxph(~ns(z_cohesion,3)+covs)); see hetero_cox.R.\n", sep="")
