# R/05_recognition/modality_check.R
# THE decisive alternative for "higher EHR recognition of distress in low-supply/rural areas":
# is it a MODALITY-of-care effect? In high-supply areas the distressed reach therapists/psychologists/
# counselors whose activity often does NOT generate an OMOP depression code; in low-supply/rural areas the
# only route is primary care, which DOES code. So EHR "recognition" may just track codeable encounters.
# Direct test: self-reported MH-professional CONTACT (access survey #44/#32) by supply tertile, among the
# distressed. If MH-professional contact is HIGHER in high-supply while EHR codes are LOWER -> modality
# explanation -> reframe from "recognition gap" to "care modality differs by area, EHR follows modality".
# Aggregate only. Self-contained; run from repo root.

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cx <- readRDS("cox_dat2.rds")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))

dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all', ans) ~ 0L, grepl('a little', ans) ~ 1L, grepl('moderate', ans) ~ 2L,
    grepl('quite a bit', ans) ~ 3L, grepl('extremely', ans) ~ 4L, TRUE ~ NA_integer_))

## self-reported MH-professional contact (ever spoken to / any visits)
spoke <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
            WHERE question LIKE '%Spoken To Mental Health Professional%'") |>
  transmute(person_id, spoke_mh = case_when(grepl('yes', ans) ~ 1L, grepl('no', ans) ~ 0L, TRUE ~ NA_integer_))
visit <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
            WHERE question LIKE '%Mental Health Professional Visits%'") |>
  transmute(person_id, mh_visit_ans = ans)

d <- cx |> left_join(dord, by="person_id") |> left_join(spoke, by="person_id") |> left_join(visit, by="person_id") |>
  filter(!is.na(sup_t), !is.na(distress_ord), distress_ord >= 2)          # among the distressed
cat("distressed n:", nrow(d), "\n")

cat("\n== self-reported 'spoken to MH professional' rate by supply tertile (among distressed) ==\n")
print(as.data.frame(d |> filter(!is.na(spoke_mh)) |> group_by(sup_t) |>
  summarise(n = n(), spoke_mh_rate = round(mean(spoke_mh), 4))))

cat("\n== 'MH professional visits' answer distribution by supply tertile (among distressed) ==\n")
print(as.data.frame(d |> filter(!is.na(mh_visit_ans)) |> count(sup_t, mh_visit_ans) |>
  group_by(sup_t) |> mutate(pct = round(100*n/sum(n),1)) |> ungroup() |> filter(n >= 20)))

cat("\nREAD: if high-supply has HIGHER self-reported MH-professional contact but LOWER EHR depression coding,\n",
    "the 'recognition' reversal is a care-MODALITY effect (rural care is more codeable), not a recognition gap.\n")
cat("\nDONE\n")
