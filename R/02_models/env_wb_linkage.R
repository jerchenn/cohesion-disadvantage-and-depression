# env_wb_linkage.R -- linkage between INDIVIDUAL environment questions and affect items,
# then the BRIDGE test: are the affect channels most coupled to cohesion the same ones
# that most forecast clinical depression?
library(sandwich); library(lmtest)
wb <- readRDS("wb_item_dat.rds")
wb_names <- attr(wb, "wb_names"); tier <- attr(wb, "tier")

env_items <- c("help","getalong","trust","values","watchout",           # cohesion
               "safe","crime","unsafe_night","unsafe_day",               # safety/danger
               "clean","upkeep","graffiti","vandalism","abandoned",      # decay
               "shops","transit","sidewalks","recreation")               # walkability
# z-score environment items (already oriented higher=better)
for (nm in env_items) wb[[paste0("ze_",nm)]] <- scale(wb[[nm]])[,1]
wb$income_z <- scale(wb$income_n)[,1]

## ---- 1. linkage matrix: standardized beta of each affect item on each env item ----
##        (bivariate, adjusted for age + income; HC3 SE)
B <- matrix(NA, length(env_items), length(wb_names),
            dimnames=list(env_items, wb_names))
for (e in env_items) for (w in wb_names){
  f <- as.formula(sprintf("z_%s ~ ze_%s + age + income_z", w, e))
  m <- try(lm(f, wb), silent=TRUE)
  if (!inherits(m,"try-error")) B[e,w] <- coef(m)[[paste0("ze_",e)]]
}
cat("=== Environment item x affect item: standardized beta (higher env better -> higher wellbeing) ===\n")
cat("(rows = individual neighbourhood questions; cols = affect items; adj age+income)\n\n")
print(round(B, 3))

## strongest couplings overall
lng <- data.frame(env=rep(rownames(B), ncol(B)),
                  wb =rep(colnames(B), each=nrow(B)),
                  beta=as.vector(B))
lng$tier <- tier[lng$wb]
cat("\n--- Top 15 environment->affect couplings ---\n")
print(head(lng[order(-lng$beta),], 15), row.names=FALSE)

## ---- 2. BRIDGE: cohesion-coupling of each affect item vs its depression HR ----
## cohesion-coupling = beta of affect item on cohesion COMPOSITE (adj age+income)
coh_beta <- sapply(wb_names, function(w){
  m <- lm(as.formula(sprintf("z_%s ~ z_cohesion + age + income_z", w)), wb)
  coef(m)[["z_cohesion"]]
})
hr <- readRDS("wb_items_hr.rds")                 # from wb_items_cox.R
hr <- hr[match(wb_names, hr$item), ]
protect <- -log(hr$HR)                            # >0 = better mental state lowers hazard
bridge <- data.frame(item=wb_names, tier=tier[wb_names],
                     cohesion_beta=round(coh_beta,3),
                     dep_protect_logHR=round(protect,3))
cat("\n=== BRIDGE: cohesion coupling vs depression-protection, by affect item ===\n")
print(bridge[order(-bridge$cohesion_beta),], row.names=FALSE)
ok <- is.finite(coh_beta) & is.finite(protect)
cat(sprintf("\nRank concordance (Spearman) cohesion-coupling vs depression-protection: rho = %.3f  (n=%d items)\n",
    cor(coh_beta[ok], protect[ok], method="spearman"), sum(ok)))
cat("Among LOW/MOD-overlap items only (wellbeing+social+gad):\n")
lowmod <- ok & tier[wb_names] %in% c("wellbeing","social","gad")
cat(sprintf("  rho = %.3f  (n=%d)\n", cor(coh_beta[lowmod], protect[lowmod], method="spearman"), sum(lowmod)))
