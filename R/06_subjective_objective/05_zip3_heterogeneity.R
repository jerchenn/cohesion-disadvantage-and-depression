# R/06_subjective_objective/05_zip3_heterogeneity.R
# How much heterogeneity is WITHIN vs BETWEEN ZIP3? Decides whether ZIP3 is a usable reference group for
# the A-probe. ZIP3 is coarse (a 3-digit prefix spans a metro/region), so within-ZIP3 likely swamps between.
# Compute the ICC (between-ZIP3 variance / total) and the spread of ZIP3 MEANS vs individual spread, for
# BMI and self-rated physical health; plus between-ZIP3 spread of the area SES measures. Aggregate only.
# Low ICC (<~0.05) => ZIP3 means are a weak/noisy reference; A-probe underpowered; lean on B.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

geo <- q("SELECT person_id, zip3_as_string, median_income, deprivation_index, fraction_poverty
          FROM `__CDR__.ds_zip_code_socioeconomic`") |> distinct(person_id, .keep_all=TRUE)
bmi <- q("SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
          WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id")
lik5 <- function(a) dplyr::case_when(grepl('excellent|excllent',a)~5, grepl('very good',a)~4,
                                     grepl('\\bgood\\b',a)~3, grepl('fair',a)~2, grepl('poor',a)~1, TRUE~NA_real_)
srh <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
          WHERE question LIKE 'Overall Health: General Physical Health'") |>
  transmute(person_id, sr_phys=lik5(ans)) |> filter(!is.na(sr_phys)) |>
  group_by(person_id) |> summarise(sr_phys=mean(sr_phys), .groups="drop")

d <- geo |> inner_join(bmi, by="person_id") |> left_join(srh, by="person_id") |> filter(!is.na(zip3_as_string))

## variance decomposition: ICC = between-ZIP3 var / (between + within)
icc <- function(df, y){
  df <- df[!is.na(df[[y]]), ]
  g  <- df |> group_by(zip3_as_string) |> summarise(m=mean(.data[[y]]), v=var(.data[[y]]), n=n(), .groups="drop") |> filter(n>=20)
  gm <- weighted.mean(g$m, g$n)
  between <- sum(g$n*(g$m-gm)^2)/sum(g$n)          # weighted between-ZIP3 variance of means
  within  <- weighted.mean(g$v, g$n, na.rm=TRUE)   # avg within-ZIP3 variance
  total   <- var(df[[y]])
  cat(sprintf("  %-10s: ICC=%.4f | between-SD=%.3f within-SD=%.3f total-SD=%.3f | ZIP3-mean p10-p90=[%.2f, %.2f] | indiv p10-p90=[%.2f, %.2f]\n",
      y, between/(between+within), sqrt(between), sqrt(within), sqrt(total),
      quantile(g$m,.1), quantile(g$m,.9), quantile(df[[y]],.1,na.rm=T), quantile(df[[y]],.9,na.rm=T)))
}
cat("== within- vs between-ZIP3 variance decomposition (ZIP3s with n>=20) ==\n")
icc(d, "bmi")
icc(d, "sr_phys")

## between-ZIP3 spread of the AREA SES measures (these are ZIP3 constants -> pure between variation)
cat("\n== between-ZIP3 spread of area SES (one value per ZIP3) ==\n")
z <- d |> group_by(zip3_as_string) |> summarise(n=n(), med_inc=first(median_income),
        dep=first(deprivation_index), pov=first(fraction_poverty), .groups="drop") |> filter(n>=20)
z$med_inc <- as.numeric(z$med_inc); z$dep <- as.numeric(z$dep); z$pov <- as.numeric(z$pov)
cat(sprintf("  median_income across ZIP3s: p10=%.0f median=%.0f p90=%.0f\n", quantile(z$med_inc,.1,na.rm=T), median(z$med_inc,na.rm=T), quantile(z$med_inc,.9,na.rm=T)))
cat(sprintf("  deprivation_index across ZIP3s: p10=%.3f median=%.3f p90=%.3f\n", quantile(z$dep,.1,na.rm=T), median(z$dep,na.rm=T), quantile(z$dep,.9,na.rm=T)))
cat(sprintf("  fraction_poverty across ZIP3s: p10=%.3f median=%.3f p90=%.3f\n", quantile(z$pov,.1,na.rm=T), median(z$pov,na.rm=T), quantile(z$pov,.9,na.rm=T)))
cat("  n ZIP3s (>=20):", nrow(z), "\n")

cat("\nREAD: low BMI/sr_phys ICC => ZIP3 lumps heterogeneous neighborhoods; neighbors'-mean reference is noisy\n",
    "=> A-probe weak, lean on B. High between-ZIP3 SES spread but low health ICC = SES differs by area but health doesn't cluster.\n")
cat("\nDONE\n")
