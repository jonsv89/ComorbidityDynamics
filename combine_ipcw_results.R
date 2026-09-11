## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript combine_ipcw_results.R --ipcw_dir ./ipcw_rr/ --out_dir ./ipcw_summary/
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
ipcw_dir <- args[which(args == "--ipcw_dir") + 1]
out_dir <- args[which(args == "--out_dir")  + 1]

files <- list.files(ipcw_dir, pattern = "^ipcw_rr_pairs_batch_.*\\.csv$", full.names = TRUE)

all_data <- rbindlist(lapply(files, fread), fill = TRUE)
cat(sprintf("Pares con RR ponderado calculado: %d\n", nrow(all_data)))

cat("\nTop 15 pairs with the higher relative increase of the RR\n")
print(all_data[order(-pct_change)][1:15,
      .(disease_a, disease_b, RR_naive = round(RR_naive, 3),
        RR_weighted = round(RR_weighted, 3), pct_change = round(pct_change, 1))])

cat("\nTop 15 pairs with the higher relative decrease of the RR\n")
print(all_data[order(pct_change)][1:15,
      .(disease_a, disease_b, RR_naive = round(RR_naive, 3),
        RR_weighted = round(RR_weighted, 3), pct_change = round(pct_change, 1))])

all_data[, sig_naive := RR_naive > 1.01]
all_data[, sig_weighted := RR_weighted > 1.01]
all_data[, flips := sig_naive != sig_weighted]

by_disease <- all_data[, .(
  n_pairs = .N,
  mean_pct_change = mean(pct_change, na.rm = TRUE),
  median_pct_change = median(as.numeric(pct_change))
), by = disease_a][order(-abs(mean_pct_change))]

print(by_disease[1:15])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(all_data, file.path(out_dir, "ipcw_all_pairs_combined.csv"))
fwrite(by_disease, file.path(out_dir, "ipcw_by_disease_summary.csv"))