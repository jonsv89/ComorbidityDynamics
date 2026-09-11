## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 06_future_defined_controls_v2.R \
#        --disease_a M51 --data_dir ./CaseControlStudy/ --out_dir ./future_controls/

library(data.table)

args      <- commandArgs(trailingOnly = TRUE)
disease_a <- args[which(args == "--disease_a") + 1]
data_dir  <- args[which(args == "--data_dir")  + 1]
out_dir   <- args[which(args == "--out_dir")   + 1]

FOLLOWUP_YEARS <- 5

cat(sprintf("\n=== Future-defined controls (v2, corregido): disease_a = %s ===\n", disease_a))

cohort <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort)
valid_diseases_any <- readRDS(file.path(data_dir, "valid_diseases_any.rds"))

patients_with_A <- cohort[cod == disease_a, unique(idp)]
cat(sprintf("Pacientes que tienen %s en algun momento: %d\n", disease_a, length(patients_with_A)))

all_first <- cohort[, {
  idx <- which.min(dat)
  .(index_date = dat[idx], index_cod = cod[idx],
    followup_end = followup_end[idx], followup_years = followup_years[idx],
    sexe = sexe[idx], rangos = rangos[idx])
}, by = idp]

candidates <- all_first[index_cod != disease_a &
                        index_cod %in% valid_diseases_any &
                        followup_years >= FOLLOWUP_YEARS]

candidates[, group := fifelse(idp %in% patients_with_A, "FUTURE_A", "NEVER_A")]

cat(sprintf("\nCandidates control-elegible (valid index, >=5y follow up): %d\n", nrow(candidates)))
cat(sprintf("  NEVER_A  (population used as control): %d (%.1f%%)\n",
            sum(candidates$group == "NEVER_A"), 100*mean(candidates$group == "NEVER_A")))
cat(sprintf("  FUTURE_A (excluded under the original approach): %d (%.1f%%)\n",
            sum(candidates$group == "FUTURE_A"), 100*mean(candidates$group == "FUTURE_A")))

fa_dates <- cohort[idp %in% candidates[group == "FUTURE_A", idp] & cod == disease_a,
                   .(dat_A = min(dat)), by = idp]
candidates <- merge(candidates, fa_dates, by = "idp", all.x = TRUE)
candidates[, years_to_A := as.numeric(dat_A - index_date) / 365.25]

cat("\nTime till developing A among the FUTURE_A (years since index):\n")
print(summary(candidates[group == "FUTURE_A", years_to_A]))

cat("\nCalculating comorbidity load before index date...\n")

cohort_sub <- cohort[idp %in% candidates$idp & cod != disease_a,
                     .(first_dat = min(dat)), by = .(idp, cod)]
rm(cohort); gc()

lut <- candidates[, .(idp, threshold = index_date)]
m <- merge(cohort_sub, lut, by = "idp", allow.cartesian = TRUE)
m <- m[first_dat <= threshold]
burden <- m[, .(n_comorbid = uniqueN(cod)), by = idp]

candidates <- merge(candidates, burden, by = "idp", all.x = TRUE)
candidates[is.na(n_comorbid), n_comorbid := 0L]

cat("\n=== Comparison: population used as control (NEVER_A) vs. the excluded one (FUTURE_A) ===\n")

summary_tbl <- candidates[, .(
  n = .N,
  mean_comorbid_previa = mean(n_comorbid),
  median_comorbid_previa = median(as.numeric(n_comorbid))
), by = group]
print(summary_tbl)

test_comorbid <- wilcox.test(n_comorbid ~ group, data = candidates)
cat(sprintf("\nWilcoxon (comorbilidad previa, NEVER_A vs FUTURE_A): p = %s\n",
            format.pval(test_comorbid$p.value, digits = 3)))

print(table(candidates$rangos, candidates$group))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("future_controls_v2_%s.csv", disease_a))
fwrite(candidates, out_file)
