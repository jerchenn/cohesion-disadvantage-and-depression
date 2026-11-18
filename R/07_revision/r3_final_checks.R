# r3_final_checks.R -- last items for the R1 manuscript: software versions, Table 1 gender rows by cohesion
# tertile, and one definition of each subcohort size. Aggregate output only; counts under 20 masked.
source("R/07_revision/_common.R", local = TRUE)

hdr("Software")
cat(R.version.string, "| survival", as.character(packageVersion("survival")), "\n")

hdr("Table 1 gender rows (same tertiles as table1.R)")
cox$coh_t <- cut(cox$z_cohesion, quantile(cox$z_cohesion, 0:3/3, na.rm = TRUE), include.lowest = TRUE, labels = c("Low", "Middle", "High"))
g <- ifelse(cox$sex == "Female", "Women", ifelse(cox$sex == "Male", "Men", "Other or not reported"))
cell <- function(k, idx) { n <- sum(g[idx] == k); sprintf("%s (%.1f)", mask(n), 100 * n / length(idx)) }
cols <- list(Overall = seq_len(nrow(cox)), Low = which(cox$coh_t == "Low"), Middle = which(cox$coh_t == "Middle"), High = which(cox$coh_t == "High"))
print(t(sapply(c("Women", "Men", "Other or not reported"), function(k) sapply(cols, function(i) cell(k, i)))), quote = FALSE)
cat("sex values present:", paste(names(table(cox$sex)), collapse = " | "), "\n")

hdr("Subcohort sizes under each definition")
wb <- readRDS("wb_item_dat.rds"); fb <- readRDS("fitbit_dat.rds")
wb_ids <- as.character(wb$person_id[!is.na(wb$event)]); fb_ids <- as.character(fb$person_id[!is.na(fb$steps) & !is.na(fb$event)])
cat("SWB: rows in wb_item_dat.rds with an event status (eTable 1 definition):", mask(length(wb_ids)), "\n")
cat("SWB: of these, also in the main cohort file:", mask(sum(cox$person_id %in% wb_ids)), "\n")
m <- fit(RHS, cox[cox$person_id %in% wb_ids, ]); cat("SWB: used in the within-subcohort Cox model (complete cases):", mask(m$n), "events", mask(m$nevent), "\n")
cat("Fitbit: rows with steps and event status:", mask(length(fb_ids)), "| in main cohort:", mask(sum(cox$person_id %in% fb_ids)), "\n")
m <- fit(RHS, cox[cox$person_id %in% fb_ids, ]); cat("Fitbit: used in the within-subcohort Cox model:", mask(m$n), "events", mask(m$nevent), "\n")

hdr("Primary model, 4 decimals")
hr(fit(RHS, cox), lab = "primary (fully adjusted)")
