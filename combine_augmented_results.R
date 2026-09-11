## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript combine_augmented_results.R --aug_dir ./augmented_rr/ --out_dir ./augmented_summary/
library(data.table)

args    <- commandArgs(trailingOnly = TRUE)
aug_dir <- args[which(args == "--aug_dir") + 1]
out_dir <- args[which(args == "--out_dir") + 1]

files <- list.files(aug_dir, pattern = "^augmented_rr_pairs_batch_.*\\.csv$", full.names = TRUE)

all_data <- rbindlist(lapply(files, fread), fill = TRUE)

all_data[, sig_orig := RR_original > 1.01]
all_data[, sig_aug   := RR_augmented > 1.01]
all_data[, flips     := sig_orig != sig_aug]

if ("c_orig" %in% names(all_data)) {
  for (thr in c(0, 10, 20, 50, 100)) {
    sub <- all_data[c_orig >= thr]
    cat(sprintf("  c_orig >= %-4d: n=%6d (%.1f%% del total) | correlacion=%.4f | mediana |cambio|=%.2f%% | cambio>=20%%: %.1f%%\n",
                thr, nrow(sub), 100*nrow(sub)/nrow(all_data),
                cor(sub$RR_original, sub$RR_augmented, use = "complete.obs"),
                median(abs(sub$pct_change), na.rm = TRUE),
                100*mean(abs(sub$pct_change) >= 20, na.rm = TRUE)))
  }
}

## ==================================================================== ##
## Top movers                                                            ##
## ==================================================================== ##
cat("\nTop 15 pairs with higher RR relative increase\n")
print(all_data[order(-pct_change)][1:15,
      .(disease_a, disease_b, RR_original = round(RR_original, 3),
        RR_augmented = round(RR_augmented, 3), pct_change = round(pct_change, 1),
        n_future_a_added)])

cat("\nTop 15 pairs with higher RR relative decrease\n")
print(all_data[order(pct_change)][1:15,
      .(disease_a, disease_b, RR_original = round(RR_original, 3),
        RR_augmented = round(RR_augmented, 3), pct_change = round(pct_change, 1),
        n_future_a_added)])

by_disease <- all_data[, .(
  n_pairs = .N,
  mean_pct_change = mean(pct_change, na.rm = TRUE),
  median_pct_change = median(as.numeric(pct_change)),
  mean_n_future_a = mean(n_future_a_added, na.rm = TRUE)
), by = disease_a][order(-abs(mean_pct_change))]

print(by_disease[1:15])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(all_data, file.path(out_dir, "augmented_all_pairs_combined.csv"))
fwrite(by_disease, file.path(out_dir, "augmented_by_disease_summary.csv"))