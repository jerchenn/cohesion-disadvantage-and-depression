# R/05_recognition/modality_confirm.R
# The two LAST tests that can kill the "EHR depression ascertainment is modality-dependent" finding.
# (After these: estimation + writing, no more diagnostics.)
# TEST 1 (within-person crossover): among distressed, 2x2 of self-reported MH-professional contact
#   (intensity: >=16 visits) x subsequent EHR depression code, by supply tertile. Modality substitution
#   predicts "contact-but-no-code" skews HIGH-supply and "code-but-no-contact" skews LOW-supply.
# TEST 2 (mechanism -- is specialty care actually less codeable?): among self-reported MH-contact users,
#   what fraction have ANY MH-related EHR footprint, and does it THIN in high-supply areas?
# Aggregate only (<20 suppressed). Self-contained; run from repo root.

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")
supp <- function(x) ifelse(x < 20, NA, x)

cx <- readRDS("cox_dat2.rds")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))

dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all',ans)~0L, grepl('a little',ans)~1L, grepl('moderate',ans)~2L,
    grepl('quite a bit',ans)~3L, grepl('extremely',ans)~4L, TRUE~NA_integer_))

## intensity of specialty MH care: >=16 visits (the stronger exposure)
vis <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
          WHERE question LIKE '%Mental Health Professional Visits%'") |>
  transmute(person_id, mh_hi = case_when(
    grepl('16 or more|10 to 12|8 to 9|6 to 7', ans) ~ 1L,          # frequent specialty contact
    grepl('visits:', ans) ~ 0L, TRUE ~ NA_integer_))

d <- cx |> left_join(dord, by="person_id") |> left_join(vis, by="person_id") |>
  filter(!is.na(sup_t), !is.na(distress_ord), distress_ord >= 2, !is.na(mh_hi))
cat("distressed w/ visit data n:", nrow(d), "\n")

## TEST 1: 2x2 (frequent MH contact) x (EHR depression code = event) by supply tertile
cat("\n== TEST 1: MH-contact(hi) x EHR-code(event), counts by supply tertile ==\n")
t1 <- d |> group_by(sup_t, mh_hi, event) |> summarise(n = supp(n()), .groups="drop")
print(as.data.frame(t1))
cat("\n-- key cells as row%: within each tertile, share of frequent-MH-contact people who got an EHR code --\n")
print(as.data.frame(d |> filter(mh_hi==1) |> group_by(sup_t) |>
  summarise(n_contact = n(), got_code_pct = round(100*mean(event),1))))
cat("-- and share of EHR-coded people who had NO frequent MH contact --\n")
print(as.data.frame(d |> filter(event==1) |> group_by(sup_t) |>
  summarise(n_coded = n(), no_contact_pct = round(100*mean(mh_hi==0),1))))

## TEST 2: mechanism -- among frequent-MH-contact people, any MH-related EHR condition code? by tertile
## MH-related concept set: descendants of Mental disorder (SNOMED 432586 via concept_ancestor) is broad;
## here use the study's own depression phenotype presence in EHR ever, as the codeability proxy.
ehrdep <- q("SELECT DISTINCT person_id FROM `__CDR__.condition_occurrence` co
             JOIN `__CDR__.concept_ancestor` ca ON ca.descendant_concept_id = co.condition_concept_id
             WHERE ca.ancestor_concept_id = 440383")   # Depressive disorder (SNOMED) -- confirm concept id
ehrdep$any_dep_code <- 1L
m <- d |> filter(mh_hi==1) |> left_join(ehrdep, by="person_id")
m$any_dep_code[is.na(m$any_dep_code)] <- 0L
cat("\n== TEST 2: among frequent-MH-contact distressed, share with ANY depression EHR code, by supply ==\n")
print(as.data.frame(m |> group_by(sup_t) |> summarise(n = n(), any_dep_code_pct = round(100*mean(any_dep_code),1))))
cat("READ: if high-supply MH-contact users have LOWER any-code% -> specialty care is less codeable (mechanism shown).\n")

cat("\nDONE\n")
