# R/05_recognition/killtests.R -- Paper 2 (recognition study) kill tests: K1 (stratified power) + K2 (positive control).
# Same repo/workflow as Paper 1: run from the repo working dir. Reads cox_dat2.rds, cox_dat3.rds, and
# zip3_mh_supply.csv from "." (already in the working dir, like the other pipeline data). Aggregate output only (<20 suppressed).
# Encodings verified against the live ds_survey answer distributions (2026-09).

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR  <- Sys.getenv("WORKSPACE_CDR")
q    <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")
supp <- function(x) ifelse(x < 20, NA, x)                       # DUCC: suppress cells < 20

## ---- cohort (re-anchored; baseline = MAX(SDOH, EHHWB)) + area MH-supply tertile ----
cx <- readRDS("cox_dat2.rds")
if (!"zip3_as_string" %in% names(cx))                          # cox_dat2 already carries it; don't re-join
  cx <- left_join(cx, readRDS("cox_dat3.rds")[, c("person_id", "zip3_as_string")], by = "person_id")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)      # AoU format "010**" -> "010"
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))
cat("cohort n:", nrow(cx), " events:", sum(cx$event),
    " | supply linked:", sum(!is.na(cx$mh_prov_per_100k)), "/", nrow(cx), "\n")

## ---- exposure: experienced healthcare discrimination (frequency scale: any = 1, none = 0) ----
disc <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%Delayed Or No Care%'")
disc <- disc |> transmute(person_id, discrim = case_when(
    grepl('always|most of the time|some of the time', ans) ~ 1L,
    grepl('none of the time', ans)                         ~ 0L,
    TRUE ~ NA_integer_))                                       # pmi: skip / dont know -> NA

## ---- distress: PTSD 'felt distant or cut off' (moderate+ = 1) ----
dsr <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'")
dsr <- dsr |> transmute(person_id, distress = case_when(
    grepl('moderate|quite a bit|extremely', ans) ~ 1L,
    grepl('not at all|a little', ans)            ~ 0L,
    TRUE ~ NA_integer_))                                       # pmi: skip -> NA

## ---- assemble + kill tests ----
d <- cx |> select(person_id, event, sup_t) |>
  left_join(disc, by = "person_id") |> left_join(dsr, by = "person_id") |>
  filter(!is.na(sup_t), !is.na(discrim), !is.na(distress))
cat("analytic n:", nrow(d), " events:", sum(d$event), "\n")
cat("discrim:", paste(names(table(d$discrim)), table(d$discrim), sep = "="),
    " | distress:", paste(names(table(d$distress)), table(d$distress), sep = "="), "\n")

cat("\n== K2 positive control: recognition rate by distress ==\n")
print(as.data.frame(d |> group_by(distress) |> summarise(n = n(), events = supp(sum(event)), rate = round(mean(event), 4))))

cat("\n== discrim x distress (collapse supply) ==\n")
print(as.data.frame(d |> group_by(discrim, distress) |> summarise(n = n(), events = supp(sum(event)), .groups = "drop")))

cat("\n== sup_t x distress (collapse discrim) ==\n")
print(as.data.frame(d |> group_by(sup_t, distress) |> summarise(n = n(), events = supp(sum(event)), .groups = "drop")))

cat("\n== K1 full 3-way: sup_t x discrim x distress (<20 suppressed) ==\n")
print(as.data.frame(d |> group_by(sup_t, discrim, distress) |> summarise(n = n(), events = supp(sum(event)), .groups = "drop")))
cat("\nDONE\n")
