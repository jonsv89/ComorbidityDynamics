## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 07_survivor_bias_5year_conditioning.R \
#        --disease_a M51 --data_dir ./CaseControlStudy/ --out_dir ./survivor_bias/

library(data.table)

args      <- commandArgs(trailingOnly = TRUE)
disease_a <- args[which(args == "--disease_a") + 1]
data_dir  <- args[which(args == "--data_dir")  + 1]
out_dir   <- args[which(args == "--out_dir")   + 1]

FOLLOWUP_YEARS <- 5

cat(sprintf("\n=== Survivor bias / 5-year conditioning: disease_a = %s ===\n", disease_a))

cohort <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort)

a_diag <- cohort[cod == disease_a, .(index_date = min(dat)), by = idp]
a_diag <- merge(a_diag, unique(cohort[cod == disease_a, .(idp, followup_end, followup_years,
                                                            sexe, rangos)]),
                by = "idp")
cat(sprintf("Patients diagnosed with %s (at any moment): %d\n", disease_a, nrow(a_diag)))

a_diag[, eligible := followup_years >= FOLLOWUP_YEARS]

n_elig <- sum(a_diag$eligible)
n_excl <- sum(!a_diag$eligible)
cat(sprintf("  Elegible (>=5 years of follow up):    %d (%.1f%%)\n",
            n_elig, 100 * n_elig / nrow(a_diag)))
cat(sprintf("  Excluidos (<5 years of follow up):     %d (%.1f%%)\n",
            n_excl, 100 * n_excl / nrow(a_diag)))

cat("\nCalculating comorbidity load before A diagnosis...\n")

cohort_sub <- cohort[idp %in% a_diag$idp & cod != disease_a,
                     .(first_dat = min(dat)), by = .(idp, cod)]
rm(cohort); gc()

lut <- a_diag[, .(idp, threshold = index_date)]
m <- merge(cohort_sub, lut, by = "idp", allow.cartesian = TRUE)
m <- m[first_dat <= threshold]
burden <- m[, .(n_comorbid = uniqueN(cod)), by = idp]

a_diag <- merge(a_diag, burden, by = "idp", all.x = TRUE)
a_diag[is.na(n_comorbid), n_comorbid := 0L]

summary_tbl <- a_diag[, .(
  n                    = .N,
  pct_mujer            = 100 * mean(sexe == "D", na.rm = TRUE),   # ajustar codigo de sexo si difiere
  mean_comorbid_previa = mean(n_comorbid),
  median_comorbid_previa = median(as.numeric(n_comorbid))
), by = eligible]

print(summary_tbl)

test_comorbid <- wilcox.test(n_comorbid ~ eligible, data = a_diag)
cat(sprintf("\nWilcoxon (crior comorbidity, elegible vs excluded): p = %s\n",
            format.pval(test_comorbid$p.value, digits = 3)))

print(table(a_diag$rangos, a_diag$eligible))

excl <- a_diag[eligible == FALSE]
print(summary(excl$followup_years))

bins <- c(-0.01, 0.5, 1, 2, 3, 4, 5)
excl[, followup_bin := cut(followup_years, breaks = bins,
                            labels = c("<0.5y","0.5-1y","1-2y","2-3y","3-4y","4-5y"))]
print(excl[, .N, by = followup_bin][order(followup_bin)])

a_diag[, rangos_f := factor(rangos)]
a_diag[, sexe_f    := factor(sexe)]

ipcw_model <- glm(eligible ~ rangos_f + sexe_f + n_comorbid, data = a_diag, family = binomial)
a_diag[, p_eligible := predict(ipcw_model, newdata = a_diag, type = "response")]
a_diag[, ipcw_weight := ifelse(eligible, 1 / p_eligible, NA_real_)]

print(summary(a_diag$p_eligible))

print(summary(a_diag[eligible == TRUE, ipcw_weight]))

cat(sprintf("\nExtreme IPCW weights (>10): %d de %d elegible (%.1f%%)\n",
            sum(a_diag[eligible == TRUE, ipcw_weight] > 10, na.rm = TRUE),
            n_elig,
            100 * sum(a_diag[eligible == TRUE, ipcw_weight] > 10, na.rm = TRUE) / n_elig))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("survivor_bias_%s.csv", disease_a))
fwrite(a_diag[, .(idp, eligible, followup_years, n_comorbid, sexe, rangos, p_eligible, ipcw_weight)],
       out_file)
cat(sprintf("\Saved: %s\n", out_file))
