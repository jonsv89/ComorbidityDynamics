## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 08_standardized_exclusion_ratio.R \
#        --data_dir ./CaseControlStudy/ --out_dir ./diagnostics_scale/
library(data.table)

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- args[which(args == "--data_dir") + 1]
out_dir  <- args[which(args == "--out_dir")  + 1]

FOLLOWUP_YEARS <- 5

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))

t0 <- Sys.time()
cat("=== Standardised Exclusion Ratio: cases vs. controls by age range ===\n")

cohort <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort)
valid_diseases_A   <- readRDS(file.path(data_dir, "valid_diseases_A.rds"))
valid_diseases_any <- readRDS(file.path(data_dir, "valid_diseases_any.rds"))
mem_log(sprintf("cohort.rds loaded: %d rows, %d patients", nrow(cohort), uniqueN(cohort$idp)))

setorder(cohort, idp, dat)
first_idx <- cohort[, .I[1], by = idp]$V1
all_first <- cohort[first_idx, .(idp, index_date_first = dat, index_cod_first = cod,
                                 followup_years_first = followup_years, sexe, rangos)]
rm(first_idx); gc()
mem_log(sprintf("all_first: %d pacientes", nrow(all_first)))

ref_pool <- all_first[index_cod_first %in% valid_diseases_any]
ref_pool[, excluded := followup_years_first < FOLLOWUP_YEARS]

ref_rates <- ref_pool[, .(
  n_total = .N,
  n_excluded = sum(excluded),
  exclusion_rate = mean(excluded)
), by = rangos]

print(ref_rates[order(rangos)])
fwrite(ref_rates, file.path(out_dir, "reference_exclusion_rates_by_age.csv"))

setkey(cohort, cod)
ref_lookup <- ref_rates[, .(rangos, exclusion_rate)]

results <- vector("list", length(valid_diseases_A))

for (i in seq_along(valid_diseases_A)) {
  disease_a <- valid_diseases_A[i]

  a_rows <- cohort[.(disease_a), .(idp, dat, followup_years, rangos), on = "cod"]
  if (nrow(a_rows) == 0 || is.na(a_rows$idp[1])) next

  a_events <- a_rows[, .SD[which.min(dat)], by = idp]
  if (nrow(a_events) < 20) next

  a_events[, excluded := followup_years < FOLLOWUP_YEARS]
  a_events <- merge(a_events, ref_lookup, by = "rangos", all.x = TRUE)

  n_obs <- sum(a_events$excluded)
  n_exp <- sum(a_events$exclusion_rate, na.rm = TRUE)   # suma de probabilidades = esperado
  n_total <- nrow(a_events)

  if (n_exp > 0) {
    ser <- n_obs / n_exp
    # Test de Poisson (igual que en SMR): observado ~ Poisson(esperado)
    p_val <- poisson.test(n_obs, T = n_exp)$p.value
  } else {
    ser <- NA_real_; p_val <- NA_real_
  }

  results[[i]] <- data.table(
    disease_a = disease_a,
    n_total = n_total,
    n_observed_excluded = n_obs,
    n_expected_excluded = round(n_exp, 1),
    observed_rate = n_obs / n_total,
    expected_rate = n_exp / n_total,
    SER = ser,
    p_value = p_val
  )

  if (i %% 100 == 0) {
    cat(sprintf("  ...%d/%d (%.1f min)\n", i, length(valid_diseases_A),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

ser_table <- rbindlist(results)
ser_table[, fdr := p.adjust(p_value, method = "BH")]
ser_table <- ser_table[order(-SER)]

fwrite(ser_table, file.path(out_dir, "standardized_exclusion_ratio_by_disease.csv"))

cat(sprintf("Disease with calculated SER: %d\n", nrow(ser_table)))
cat(sprintf("medium SER: %.2f | SER median: %.2f\n",
            mean(ser_table$SER, na.rm = TRUE), median(ser_table$SER, na.rm = TRUE)))
cat(sprintf("Diseases with SER > 1.5 and significant (FDR<0.05): %d (%.1f%%)\n",
            sum(ser_table$SER > 1.5 & ser_table$fdr < 0.05, na.rm = TRUE),
            100 * mean(ser_table$SER > 1.5 & ser_table$fdr < 0.05, na.rm = TRUE)))
cat(sprintf("Diseases with SER significatively < 1 (loose less than expected, FDR<0.05): %d (%.1f%%)\n",
            sum(ser_table$SER < 1 & ser_table$fdr < 0.05, na.rm = TRUE),
            100 * mean(ser_table$SER < 1 & ser_table$fdr < 0.05, na.rm = TRUE)))

cat("\nTop 15 with larger exclusion excess (higher SER):\n")
print(ser_table[1:15, .(disease_a, n_total, observed_rate = round(observed_rate, 3),
                         expected_rate = round(expected_rate, 3), SER = round(SER, 2), fdr = round(fdr, 4))])

cat("\nTop 15 with lower relative exclusion (lower SER):\n")
print(tail(ser_table, 15)[, .(disease_a, n_total, observed_rate = round(observed_rate, 3),
                               expected_rate = round(expected_rate, 3), SER = round(SER, 2), fdr = round(fdr, 4))])