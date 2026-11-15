# r2_multilevel.R -- Reviewer 2: multilevel structure (ZIP3 random effects) and contextual vs individual cohesion.
#  (a) primary on the ZIP3-available sample: naive, cluster-robust, and ZIP3 gaussian frailty (random intercept)
#  (b) area-mean cohesion (leave-one-out mean of the respondents' scores in the ZIP3) as a contextual term
#  (c) within/between decomposition: individual deviation from the area mean + area mean (Mundlak form)
#  (d) (c) plus area deprivation, to see whether the contextual term is distinct from deprivation
# No coxme on the Workbench image, so the random intercept uses survival::frailty(). Aggregate output only.
source("R/07_revision/_common.R", local = TRUE)

d <- cox[!is.na(cox$zip3_as_string) & !is.na(cox$z_cohesion), ]
d$clus <- factor(d$zip3_as_string)
cs <- as.numeric(table(d$clus))
cat(sprintf("ZIP3-available participants with a cohesion score: %s in %s ZIP3 areas; participants per area median %.0f (IQR %.0f-%.0f); areas with fewer than 20 participants: %s\n",
            mask(nrow(d)), mask(nlevels(d$clus)), median(cs), quantile(cs, .25), quantile(cs, .75), mask(sum(cs < 20))))

## leave-one-out area mean so that a person's own score does not enter the contextual term
S <- tapply(d$z_cohesion, d$clus, sum); N <- tapply(d$z_cohesion, d$clus, length)
d$area_mean <- as.numeric((S[d$clus] - d$z_cohesion) / pmax(N[d$clus] - 1, 1))
d$area_mean[N[d$clus] == 1] <- NA
d$within <- d$z_cohesion - d$area_mean
cat(sprintf("area-mean cohesion: SD across participants %.3f (individual score SD %.3f); share of total variance that is between areas: %.1f%%\n",
            sd(d$area_mean, na.rm = TRUE), sd(d$z_cohesion), 100 * var(ave(d$z_cohesion, d$clus), na.rm = TRUE) / var(d$z_cohesion)))

hdr("(a) Primary on the ZIP3-available sample")
hr(fit(RHS, d), lab = "naive SE")
hr(fit(paste(RHS, "+ cluster(clus)"), d), lab = "ZIP3 cluster-robust SE")
mf <- fit(paste(RHS, '+ frailty(clus, distribution = "gaussian")'), d)
hr(mf, lab = "ZIP3 gaussian frailty (random intercept)")
th <- tryCatch(mf$history[[1]]$theta, error = function(e) NA)
cat(sprintf("  random-intercept variance %.4f (SD %.3f)\n", th, sqrt(th)))

hdr("(b) Individual score plus area-mean (contextual) cohesion, ZIP3 frailty")
dd <- d[!is.na(d$area_mean), ]
mb <- fit(paste(RHS, '+ area_mean + frailty(clus, distribution = "gaussian")'), dd)
hr(mb, "z_cohesion", "individual cohesion, holding area mean fixed")
hr(mb, "area_mean", "area-mean cohesion (per 1 SD of individual scale), holding individual fixed")

hdr("(c) Within/between decomposition (Mundlak), ZIP3 frailty")
mc <- fit(paste(sub("z_cohesion", "within + area_mean", RHS), '+ frailty(clus, distribution = "gaussian")'), dd)
hr(mc, "within",    "within-area (individual deviation from area mean)")
hr(mc, "area_mean", "between-area (area mean); difference from within = contextual effect")
b <- coef(mc)[c("within", "area_mean")]; V <- vcov(mc)[c("within", "area_mean"), c("within", "area_mean")]
dlt <- b[["area_mean"]] - b[["within"]]; se <- sqrt(V[1,1] + V[2,2] - 2*V[1,2])
cat(sprintf("  between minus within (log HR): %.3f (95%% CI %.3f to %.3f), P=%.2g\n", dlt, dlt - 1.96*se, dlt + 1.96*se, 2*pnorm(-abs(dlt/se))))

hdr("(d) Decomposition plus area deprivation")
dd$dep_z <- as.numeric(scale(dd$deprivation_index))
ddd <- dd[!is.na(dd$dep_z), ]
md <- fit(paste(sub("z_cohesion", "within + area_mean + dep_z", RHS), '+ frailty(clus, distribution = "gaussian")'), ddd)
hr(md, "within",    "within-area cohesion, adjusted for area deprivation")
hr(md, "area_mean", "between-area cohesion, adjusted for area deprivation")
hr(md, "dep_z",     "area deprivation index (per SD)")
cat("\nNote for the text: cohesion is measured at the individual level and aggregated only as a contextual covariate; ZIP3 is the finest geography released, so the area mean averages over many neighbourhoods.\n")
