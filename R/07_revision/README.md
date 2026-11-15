# Revision scripts (JAMA Network Open JNO26-18095, decision of 9 Nov 2026)

Every script is an additional sensitivity analysis. The primary specification in cox_model3.R is not changed. All output is aggregate; counts below 20 print as "<20"; no minima or maxima of participant values are printed.

Run after the usual startup block, in this order (each is independent; the BigQuery ones take longer):

```r
run("R/07_revision/r2_reri_subsamples_vif.R")      # no BigQuery, fast
run("R/07_revision/r2_multilevel.R")               # no BigQuery; frailty fits take a few minutes
run("R/07_revision/r2_ipcw.R")                     # no BigQuery
run("R/07_revision/r2_mi_primary.R")               # no BigQuery; 6 x 20 imputations, ~10 min
run("R/07_revision/r1_exclusions_lags_ehrlength.R") # BigQuery: psychiatric concept sets
run("R/07_revision/r2_ascertainment.R")            # BigQuery: contact counts, visits, negative controls
run("R/07_revision/r1_responders.R")               # BigQuery: all EHR-linked participants (large pull, stays in memory)
```

| Script | Letter item | What it prints |
|---|---|---|
| r1_exclusions_lags_ehrlength.R | R1 psychiatric exclusions; R1 EHR length; R2 longer lags | primary after excluding prior bipolar, psychotic, substance use (and anxiety); prior EHR years as covariate; 365- and 730-day lags |
| r1_responders.R | R1 selection of SDOH responders | responders vs all EHR-linked participants; response-propensity weighted primary |
| r2_multilevel.R | R2 multilevel, contextual vs individual | ZIP3 frailty model; area-mean cohesion; within/between decomposition; with deprivation |
| r2_ascertainment.R | R2 informative contact; primary-care visits; negative controls | post-baseline contact by cohesion, income, deprivation; outpatient visits as covariate; cataract, hernia, fracture as negative controls |
| r2_ipcw.R | R2 informative censoring | censoring model; stabilised IPCW; weighted primary |
| r2_mi_primary.R | R2 MI as primary, auxiliaries, MNAR | MI (m = 20) with the imputation model stated; delta-adjusted MNAR scenarios |
| r2_reri_subsamples_vif.R | R2 additive interaction and multiplicity; subsample replication; item collinearity | RERI, AP, synergy index; Bonferroni; cohesion HR within SWB and Fitbit subcohorts; VIF and condition number |

Text-only items from the letter (no script): entry-criterion justification, deprivation index description, ZIP3 coarseness in abstract and conclusions, "clinically ascertained" wording, causal language and "predict", race and ethnicity and gender reporting, Key Points and Conclusions phrasing, Table 1 male row, forest plots on log scale, code list in the supplement (from the concept-set definition).
