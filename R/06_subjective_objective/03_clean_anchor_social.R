# R/06_subjective_objective/03_clean_anchor_social.R
# After 02: age discordance was largely observation/coding artifact; BMI (accumulation-free) showed NO age
# shift. So (a) use PHYSIOLOGICAL anchors not EHR counts, (b) test discordance across NON-AGE social groups.
# QUESTION: does the self-rated-physical-health -> BMI mapping differ by RACE, INCOME, EDUCATION? If the same
# objective BMI maps to different self-ratings across social groups, that's socially-patterned reporting
# heterogeneity on a CLEAN anchor -- the SSM thesis. Print social-var distributions first (lock mappings).
# Aggregate only. Individual-level.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

lik5 <- function(a) case_when(grepl('excellent|excllent',a)~5L, grepl('very good',a)~4L,
                              grepl('\\bgood\\b',a)~3L, grepl('fair',a)~2L, grepl('poor',a)~1L, TRUE~NA_integer_)
phys <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
           WHERE question LIKE 'Overall Health: General Physical Health'") |>
  transmute(person_id, v=lik5(ans)) |> filter(!is.na(v)) |> group_by(person_id) |> summarise(sr_phys=mean(v),.groups="drop")
bmi <- q("SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
          WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id")

## social moderators
race <- q("SELECT p.person_id, c.concept_name race FROM `__CDR__.person` p
           LEFT JOIN `__CDR__.concept` c ON c.concept_id=p.race_concept_id")
inc <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE LOWER(question) LIKE '%annual income%'")
edu <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE LOWER(question) LIKE '%highest grade%' OR LOWER(question) LIKE '%education%level%'")
cat("== income answer distribution ==\n"); print(sort(table(inc$ans),decreasing=TRUE))
cat("\n== education answer distribution ==\n"); print(sort(table(edu$ans),decreasing=TRUE))
cat("\n== race distribution ==\n"); print(sort(table(race$race),decreasing=TRUE))

inc <- inc |> mutate(income = case_when(
    grepl('less than 10k|10k.*25k|25k.*35k', ans) ~ "<35k",
    grepl('35k.*50k|50k.*75k', ans) ~ "35-75k",
    grepl('75k.*100k|100k.*150k|150k.*200k|more than 200k|200k', ans) ~ "75k+",
    TRUE ~ NA_character_)) |> filter(!is.na(income)) |> distinct(person_id,.keep_all=TRUE) |> select(person_id,income)
edu <- edu |> mutate(educ = case_when(
    grepl('never|grade|one through four|five through eight|nine through eleven|twelve.*no diploma', ans) ~ "<HS",
    grepl('high school|ged|twelve.*diploma', ans) ~ "HS",
    grepl('college.*no degree|some college|advanced degree.*no|one through three', ans) ~ "SomeCol",
    grepl('college graduate|advanced degree|bachelor|master|doctor', ans) ~ "College+",
    TRUE ~ NA_character_)) |> filter(!is.na(educ)) |> distinct(person_id,.keep_all=TRUE) |> select(person_id,educ)
race <- race |> mutate(race=case_when(grepl('White',race)~"White", grepl('Black|African',race)~"Black",
    grepl('Asian',race)~"Asian", grepl('More than one|Multi',race)~"Multi", TRUE~"Other/NA")) |>
  distinct(person_id,.keep_all=TRUE)

d <- phys |> inner_join(bmi,by="person_id") |> left_join(race,by="person_id") |>
  left_join(inc,by="person_id") |> left_join(edu,by="person_id")
d$bmi_g <- cut(d$bmi, c(0,25,30,35,100), labels=c("<25","25-30","30-35","35+"))
cat("\nanalytic n (self-rated phys + BMI):", nrow(d), "\n")

## CLEAN-ANCHOR DISCORDANCE: mean self-rated physical health at FIXED BMI category, across social groups.
## If sr_phys at the SAME BMI differs across groups -> socially-patterned reporting heterogeneity.
cat("\n== mean self-rated physical health at fixed BMI, by RACE ==\n")
print(as.data.frame(d |> filter(!is.na(bmi_g),race!="Other/NA") |> group_by(bmi_g,race) |>
  summarise(n=n(), sr_phys_mean=round(mean(sr_phys),2), .groups="drop") |> filter(n>=20)))
cat("\n== ... by INCOME ==\n")
print(as.data.frame(d |> filter(!is.na(bmi_g),!is.na(income)) |> group_by(bmi_g,income) |>
  summarise(n=n(), sr_phys_mean=round(mean(sr_phys),2), .groups="drop") |> filter(n>=20)))
cat("\n== ... by EDUCATION ==\n")
print(as.data.frame(d |> filter(!is.na(bmi_g),!is.na(educ)) |> group_by(bmi_g,educ) |>
  summarise(n=n(), sr_phys_mean=round(mean(sr_phys),2), .groups="drop") |> filter(n>=20)))

cat("\nDONE\n")
