# broad_discover.R -- BROAD inventory across the two source surveys (SDOH + EHHWB) plus the selected
# Basics items (disability battery, tenure, employment). Prints, untruncated:
#   (A) full question catalog: concept_id | survey | coverage n | full question text
#   (B) answer options per item (one compact line per concept_id)
# and writes both to CSV for a durable reference. This maps every candidate predictor/outcome to its
# concept_id and response scale so we can build recodes for all families. Aggregate only.
#
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

surveys <- c("Social Determinants of Health", "Emotional Health History and Well-Being")
basics  <- c(903573,903574,903575,903576,903577,903578,   # disability battery
             1585370,                                       # home ownership / tenure
             1585952)                                       # employment status
filt <- sprintf("(survey IN (%s) OR question_concept_id IN (%s))",
                paste(sprintf("'%s'", surveys), collapse=","), paste(basics, collapse=","))

## ---- (A) question catalog ----
q <- run_sql(sprintf("SELECT survey, question_concept_id, question, COUNT(DISTINCT person_id) AS n
                      FROM `%s.ds_survey` WHERE %s
                      GROUP BY survey, question_concept_id, question", cdr, filt))
q$n <- as.numeric(q$n); q$question <- gsub("\\s+", " ", q$question)
q <- q[order(q$survey, -q$n), ]
write.csv(q, "broad_catalog.csv", row.names=FALSE)
for (s in unique(q$survey)){
  sub <- q[q$survey==s, ]
  cat(sprintf("\n########## %s : %d questions ##########\n", s, nrow(sub)))
  for (i in seq_len(nrow(sub)))
    cat(sprintf("%s | n=%6d | %s\n", sub$question_concept_id[i], sub$n[i], sub$question[i]))
}

## ---- (B) answer options per item (compact, ordered by frequency) ----
a <- run_sql(sprintf("SELECT question_concept_id, answer, COUNT(*) AS n
                      FROM `%s.ds_survey` WHERE %s
                      GROUP BY question_concept_id, answer", cdr, filt))
a$n <- as.numeric(a$n)
write.csv(a, "broad_answers.csv", row.names=FALSE)
cat("\n\n########## answer options per item (concept_id : options high->low freq) ##########\n")
for (id in unique(a$question_concept_id)){
  ai <- a[a$question_concept_id==id, ]; ai <- ai[order(-ai$n), ]
  cat(sprintf("%s : %s\n", id, paste(ai$answer, collapse="  |  ")))
}
cat(sprintf("\nCatalog: %d questions across %d surveys+basics. CSVs: broad_catalog.csv, broad_answers.csv\n",
            nrow(q), length(unique(q$survey))))
