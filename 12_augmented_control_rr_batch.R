## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 12_augmented_control_rr_batch.R \
#        --pairs_file ./batches/pairs_batch_0001.txt \
#        --data_dir ./CaseControlStudy/ --matched_dir ./matched/ \
#        --future_a_dir ./future_a_candidates/ --out_dir ./augmented_rr/
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pairs_file <- args[which(args == "--pairs_file")   + 1]
data_dir <- args[which(args == "--data_dir")     + 1]
matched_dir <- args[which(args == "--matched_dir")  + 1]
future_a_dir <- args[which(args == "--future_a_dir") + 1]
out_dir <- args[which(args == "--out_dir")      + 1]

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))
t0 <- Sys.time()

cat(sprintf("\n=== Augmented-control RR BATCH: %s ===\n", pairs_file))
mem_log("inicio")

pairs_tbl <- fread(pairs_file, header = FALSE, col.names = c("disease_a", "disease_b"))

cohort <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort)
setkey(cohort, cod)
mem_log(sprintf("cohort.rds cargado: %d filas", nrow(cohort)))

future_a_all <- readRDS(file.path(future_a_dir, "future_a_candidates.rds"))
setDT(future_a_all)
mem_log(sprintf("future_a_candidates saved: %d rows", nrow(future_a_all)))

set.seed(42) 

MAX_FUTURE_A_MULTIPLE <- 1

compute_rr_simple <- function(n1, a, n0, c) {
  if (a == 0 || c == 0 || n1 == 0 || n0 == 0) return(NA_real_)
  (a / n1) / (c / n0)
}

results <- vector("list", nrow(pairs_tbl))
n_ok <- 0L; n_skip <- 0L

for (i in seq_len(nrow(pairs_tbl))) {
  da_val <- pairs_tbl$disease_a[i]
  db_val <- pairs_tbl$disease_b[i]

  in_file <- file.path(matched_dir, sprintf("matched_complete_%s_%s.rds", da_val, db_val))
  if (!file.exists(in_file)) { n_skip <- n_skip + 1L; next }

  res <- tryCatch({
    obj     <- readRDS(in_file)
    matched <- obj$matched
    setDT(matched)

    case_has_B <- !is.na(matched$case_dat_B)
    ctrl_has_B <- !is.na(matched$ctrl_dat_B)

    n1 <- uniqueN(matched$case_idp)
    a  <- uniqueN(matched[case_has_B, case_idp])
    n0_orig <- uniqueN(matched$ctrl_idp)
    c_orig  <- uniqueN(matched[ctrl_has_B, ctrl_idp])

    rr_original <- compute_rr_simple(n1, a, n0_orig, c_orig)

    fa_sub <- future_a_all[disease_a == da_val]

    if (nrow(fa_sub) < 5) {
      list(rr_original = rr_original, rr_augmented = NA_real_, n_future_a = nrow(fa_sub),
           n0_orig = n0_orig, c_orig = c_orig)
    } else {
      cap <- n0_orig * MAX_FUTURE_A_MULTIPLE
      n_available <- nrow(fa_sub)
      if (n_available > cap) {
        fa_sub <- fa_sub[sample(.N, cap)]
      }

      fa_sub[, window_end := index_date_first + round(5 * 365.25)]

      b_rows <- cohort[.(db_val), .(idp, dat), on = "cod"]
      if (nrow(b_rows) == 0 || is.na(b_rows$idp[1])) {
        fa_sub[, dat_B := as.Date(NA)]
      } else {
        b_first <- b_rows[, .(dat_B = min(dat)), by = idp]
        fa_sub <- merge(fa_sub, b_first, by = "idp", all.x = TRUE)
      }

      fa_sub[, has_B := !is.na(dat_B) & dat_B >= index_date_first & dat_B <= window_end]

      n0_aug <- n0_orig + uniqueN(fa_sub$idp)
      c_aug  <- c_orig + sum(fa_sub$has_B)

      rr_augmented <- compute_rr_simple(n1, a, n0_aug, c_aug)

      list(rr_original = rr_original, rr_augmented = rr_augmented, n_future_a = nrow(fa_sub),
           n_future_a_available = n_available,
           n0_orig = n0_orig, c_orig = c_orig, n0_aug = n0_aug, c_aug = c_aug)
    }
  }, error = function(e) NULL)

  if (is.null(res) || is.na(res$rr_original) || is.na(res$rr_augmented)) {
    n_skip <- n_skip + 1L
    next
  }

  results[[i]] <- data.table(
    disease_a = da_val, disease_b = db_val,
    RR_original = res$rr_original, RR_augmented = res$rr_augmented,
    n_future_a_added = res$n_future_a,
    n_future_a_available = res$n_future_a_available,
    n0_orig = res$n0_orig, n0_aug = res$n0_aug,
    c_orig = res$c_orig, c_aug = res$c_aug
  )
  n_ok <- n_ok + 1L

  if (i %% 25 == 0) {
    mem_log(sprintf("processed %d/%d (ok=%d, skip=%d)", i, nrow(pairs_tbl), n_ok, n_skip))
  }
}

combined <- rbindlist(results)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
batch_id <- gsub("\\.txt$", "", basename(pairs_file))
out_file <- file.path(out_dir, sprintf("augmented_rr_%s.csv", batch_id))

if (nrow(combined) > 0) {
  combined[, pct_change := 100 * (RR_augmented - RR_original) / RR_original]
  fwrite(combined, out_file)
} else {
  cat("None pair produced valid results in this batch\n")
}

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("\nBatch completed: %d ok, %d without data/candidates (%.1f min)\n", n_ok, n_skip, elapsed))