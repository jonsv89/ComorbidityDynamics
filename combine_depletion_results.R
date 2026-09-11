## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript combine_depletion_results.R --depletion_dir ./depletion/ --out_dir ./depletion_summary/
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
depletion_dir <- args[which(args == "--depletion_dir") + 1]
out_dir <- args[which(args == "--out_dir")        + 1]

files <- list.files(depletion_dir, pattern = "^depletion_pairs_batch_.*\\.csv$", full.names = TRUE)
all_data <- rbindlist(lapply(files, fread), fill = TRUE)

global_by_window <- all_data[, .(
  n_pairs = .N,
  mean_excess = mean(excess, na.rm = TRUE),
  median_excess = median(as.numeric(excess))
), by = .(window, t_start)][order(t_start)]

print(global_by_window)

pair_trends <- all_data[, {
  if (.N >= 3 && !any(is.na(excess))) {
    fit <- tryCatch(lm(excess ~ t_start), error = function(e) NULL)
    if (!is.null(fit) && !is.na(coef(fit)["t_start"])) {
      s <- summary(fit)$coefficients
      .(slope = s["t_start", "Estimate"], p = s["t_start", "Pr(>|t|)"], n_windows = .N)
    } else {
      .(slope = NA_real_, p = NA_real_, n_windows = .N)
    }
  } else {
    .(slope = NA_real_, p = NA_real_, n_windows = .N)
  }
}, by = .(disease_a, disease_b)]

pair_trends <- pair_trends[!is.na(slope)]
pair_trends[, fdr := p.adjust(p, method = "BH")]

pair_trends[, direction := fcase(
  fdr < 0.05 & slope < 0, "Narrows (supports depletion)",
  fdr < 0.05 & slope > 0, "Widens (against depletion)",
  rep(TRUE, .N), "Without significant tendency"
)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(global_by_window, file.path(out_dir, "depletion_global_by_window.csv"))
fwrite(pair_trends, file.path(out_dir, "depletion_pair_trends.csv"))