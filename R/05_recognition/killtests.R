# R/05_recognition/killtests.R
# Paper 2 (recognition study) -- lean kill-test build: K1 (stratified power) + K2 (positive control).
# Same repo/workflow as Paper 1: run from the repo working dir; reads cox_dat2.rds, cox_dat3.rds, and
# zip3_mh_supply.csv from "." (the working dir, as with the other pipeline scripts). Aggregate output only.

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)),
  bigint = "character")

## cohort: re-anchored (baseline = MAX(SDOH, EHHWB)); has event/time_days
cx <- readRDS("cox_dat2.rds")
cat("cohort n:", nrow(cx), " events:", sum(cx$event), "\n cols:", paste(names(cx), collapse = ", "), "\n")

## area supply tertile (zip3 sourced from cox_dat3 to be safe; supply lookup uploaded to the working dir)
z3  <- readRDS("cox_dat3.rds")[, c("person_id", "zip3_as_string")]
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, z3, by = "person_id")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)                 # "010**" -> "010"
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))

## DISCOVERY: EHHWB items (to finalize the distress measure)
cat("\n== EHHWB items ==\n")
print(q("SELECT question, COUNT(DISTINCT person_id) n FROM `__CDR__.ds_survey`
         WHERE survey='Emotional Health History and Well-Being' GROUP BY question ORDER BY question"), n = 100)

## access exposure: experienced healthcare discrimination (binary)
disc <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
           WHERE question LIKE '%Delayed Or No Care%'")
cat("\n== discrimination-item answers ==\n"); print(table(disc$ans))
disc <- disc |> transmute(person_id, discrim = as.integer(grepl('yes', ans)))

## distress proxy: EHHWB 'distant or cut off' (binary moderate+); REFINE after seeing items
dsr <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
          WHERE survey='Emotional Health History and Well-Being'
            AND (LOWER(question) LIKE '%distant%' OR LOWER(question) LIKE '%cut off%')")
cat("\n== distress-item answers ==\n"); print(table(dsr$ans))
dsr <- dsr |> transmute(person_id, distress = as.integer(!grepl('not at all|a little|never|skip|prefer not', ans)))

## assemble + run K2, K1
d <- cx |> select(person_id, event, sup_t) |>
  left_join(disc, by = "person_id") |> left_join(dsr, by = "person_id") |>
  filter(!is.na(sup_t), !is.na(discrim), !is.na(distress))
cat("\nanalytic n:", nrow(d), " events:", sum(d$event), "\n")

cat("\n== K2 (positive control): recognition rate by distress ==\n")
print(d |> group_by(distress) |> summarise(n = n(), events = sum(event), rate = round(mean(event), 4)))

cat("\n== K1 (stratified power): EVENTS per supply x discrim x distress cell (<20 suppressed) ==\n")
print(d |> group_by(sup_t, discrim, distress) |>
        summarise(events = sum(event), n = n(), .groups = "drop") |>
        mutate(events = ifelse(events < 20, NA, events), n = ifelse(n < 20, NA, n)), n = 100)
