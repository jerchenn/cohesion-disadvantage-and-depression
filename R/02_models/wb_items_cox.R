# wb_items_cox.R -- item-level affect -> incident clinical depression.
# Each affect item's HR per SD (oriented higher = better mental state => protective HR<1),
# tagged by construct-overlap tier; full + 180-day washout. The novel, non-circular
# signal lives in the LOW-overlap tiers (wellbeing, social) that still forecast Dx.
library(survival)
wb <- readRDS("wb_item_dat.rds")
wb_names <- attr(wb, "wb_names"); tier <- attr(wb, "tier")
wb <- wb[!is.na(wb$event), ]                                   # at-risk cohort only

## collapsed covariate block (same as Study A / SWB models)
wb$sex_c  <- factor(ifelse(wb$sex %in% c("Female","Male"), wb$sex, "Other"))
wb$race_c <- factor(ifelse(wb$race %in% c("White","Black or African American"), wb$race, "Other"))
wb$ethn_c <- factor(ifelse(wb$ethnicity=="Hispanic or Latino","Hispanic",
                    ifelse(wb$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
wb$income_m <- ifelse(is.na(wb$income_n), median(wb$income_n, na.rm=TRUE), wb$income_n)
wb$educ_m   <- ifelse(is.na(wb$educ_n),   median(wb$educ_n,   na.rm=TRUE), wb$educ_n)
covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"

one <- function(v, dat){
  m <- tryCatch(coxph(as.formula(paste("Surv(time_days, event) ~", v, "+", covs)), dat),
                error=function(e) NULL)
  if (is.null(m)) return(c(hr=NA, lo=NA, hi=NA, p=NA, nev=NA))
  ci <- confint(m)
  c(hr=exp(coef(m)[[v]]), lo=exp(ci[v,1]), hi=exp(ci[v,2]),
    p=summary(m)$coefficients[v,"Pr(>|z|)"], nev=m$nevent)
}

res <- data.frame()
for (nm in wb_names){
  v <- paste0("z_", nm)
  f  <- one(v, wb)
  w  <- one(v, wb[wb$time_days > 180, ])
  res <- rbind(res, data.frame(item=nm, tier=tier[nm],
               HR=f["hr"], lo=f["lo"], hi=f["hi"], p=f["p"], events=f["nev"],
               HR_wash=w["hr"], lo_w=w["lo"], hi_w=w["hi"]))
}
res <- res[order(factor(res$tier, levels=c("wellbeing","social","gad","phq")), res$HR), ]
fmt <- function(h,l,u) sprintf("%.3f (%.3f-%.3f)", h,l,u)
out <- data.frame(item=res$item, tier=res$tier,
                  HR_per_SD=fmt(res$HR,res$lo,res$hi),
                  HR_180d_washout=fmt(res$HR_wash,res$lo_w,res$hi_w),
                  p=signif(res$p,2), events=res$events)
cat("=== Affect item -> incident depression (HR per SD, higher=better => HR<1 protective) ===\n")
cat("Ordered within tier. LOW-overlap tiers (wellbeing, social) are the non-circular evidence.\n\n")
print(out, row.names=FALSE)

## joint model: the three LOW/MOD-overlap channels together (independent contributions)
low <- c("z_happy","z_meaning","z_cutoff")
mj <- coxph(as.formula(paste("Surv(time_days, event) ~", paste(low, collapse=" + "), "+", covs)), wb)
cat("\n=== Joint low-overlap channels (happy + meaning + cutoff), mutually adjusted ===\n")
for (v in low) cat(sprintf("  %-10s HR %.3f (%.3f-%.3f)  p=%.2g\n", v,
    exp(coef(mj)[[v]]), exp(confint(mj)[v,1]), exp(confint(mj)[v,2]),
    summary(mj)$coefficients[v,"Pr(>|z|)"]))
cat("  model events:", mj$nevent, " N:", mj$n, "\n")

saveRDS(res, "wb_items_hr.rds")   # for the bridge test in env_wb_linkage.R
