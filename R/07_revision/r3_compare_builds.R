# r3_compare_builds.R -- why the fresh rebuild (~/replication_check) differs from the cached cohort file used for the
# manuscript (~/cohesion-disadvantage-and-depression) in the third decimal. Aggregate output only: for each
# variable, the number of participants whose value differs between the two files (counts under 20 masked).
mask <- function(x) ifelse(x > 0 & x < 20, "<20", format(x, big.mark = ","))
A <- readRDS("~/cohesion-disadvantage-and-depression/cox_dat3.rds"); B <- readRDS("~/replication_check/cox_dat3.rds")
cat("cached file modified:", format(file.mtime("~/cohesion-disadvantage-and-depression/cox_dat3.rds")), "| fresh:", format(file.mtime("~/replication_check/cox_dat3.rds")), "\n")
cat("rows", mask(nrow(A)), mask(nrow(B)), "| same person_ids:", setequal(A$person_id, B$person_id), "\n")
cat("columns only in cached:", paste(setdiff(names(A), names(B)), collapse = ", "), "| only in fresh:", paste(setdiff(names(B), names(A)), collapse = ", "), "\n")
B <- B[match(A$person_id, B$person_id), ]
for (v in intersect(names(A), names(B))) {
  a <- A[[v]]; b <- B[[v]]
  if (is.numeric(a) && is.numeric(b)) d <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9))
  else d <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b)))
  n <- sum(d, na.rm = TRUE)
  if (n > 0) {
    cat(sprintf("%-22s differs for %s participants", v, mask(n)))
    if (is.numeric(a) && is.numeric(b)) { ok <- !is.na(a) & !is.na(b)
      cat(sprintf(" | correlation %.4f, mean abs difference %.4f, missing in cached %s vs fresh %s", cor(a[ok], b[ok]), mean(abs(a[ok] - b[ok])), mask(sum(is.na(a))), mask(sum(is.na(b))))) }
    else { tb <- table(cached = addNA(factor(a))[d], fresh = addNA(factor(b))[d]); tb[tb > 0 & tb < 20] <- NA; cat("\n"); print(tb) }
    cat("\n") }
}
cat("variables not listed are identical in both files\n")
cat("\nlast commits touching build_cox3.R:\n"); system("cd ~/cohesion-disadvantage-and-depression && git log -3 --format='%h %ad %s' --date=short -- R/01_build/build_cox3.R")
