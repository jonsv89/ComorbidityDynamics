## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 05_depletion_footprint_batch.R \
#        --pairs_file ./batches/pairs_batch_0001.txt \
#        --data_dir ./CaseControlStudy/ --matched_dir ./matched/ \
#        --out_dir ./depletion/
library(data.table)

args        <- commandArgs(trailingOnly = TRUE)
pairs_file  <- args[which(args == "--pairs_file")  + 1]
data_dir    <- args[which(args == "--data_dir")    + 1]
matched_dir <- args[which(args == "--matched_dir") + 1]
out_dir     <- args[which(args == "--out_dir")     + 1]

mem_log <- function(msg) {
  mem_mb <- sum(gc()[, 2])
  cat(sprintf("[MEM %6.0f MB] %s\n", mem_mb, msg))
}

t_start_job <- Sys.time()
cat(sprintf("\n=== Depletion footprint BATCH: %s ===\n", pairs_file))
mem_log("inicio")

pairs_tbl <- fread(pairs_file, header = FALSE, col.names = c("disease_a", "disease_b"))
cat(sprintf("Pairs in this batch: %d\n", nrow(pairs_tbl)))

cohort_full <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort_full)
setkey(cohort_full, idp)
mem_log(sprintf("cohort.rds loaded: %d rows, %d patients",
                nrow(cohort_full), uniqueN(cohort_full$idp)))

windows_cont <- list(
  "1_2" = c(1, 2), "2_3" = c(2, 3), "3_4" = c(3, 4), "4_5" = c(4, 5)
)

at_risk_start <- function(index_date, dat_B, followup_end, t_start, t_end) {
  days_start   <- t_start * 365.25
  days_end     <- t_end   * 365.25
  had_B_before <- !is.na(dat_B) & as.numeric(dat_B - index_date) <= days_start
  has_followup <- as.numeric(followup_end - index_date) >= days_end
  !had_B_before & has_followup
}

comorbidity_burden <- function(cohort_sub, idps, threshold_dates) {
  lut <- data.table(idp = idps, threshold = threshold_dates)
  m <- merge(cohort_sub, lut, by = "idp", allow.cartesian = TRUE)
  m <- m[first_dat <= threshold]
  burden <- m[, .(n_comorbid = uniqueN(cod)), by = idp]
  out <- merge(data.table(idp = idps), burden, by = "idp", all.x = TRUE)
  out[is.na(n_comorbid), n_comorbid := 0L]
  out$n_comorbid
}

all_results <- list()
n_ok <- 0L; n_skip <- 0L; n_err <- 0L

for (i in seq_len(nrow(pairs_tbl))) {

  disease_a <- pairs_tbl$disease_a[i]
  disease_b <- pairs_tbl$disease_b[i]

  in_file <- file.path(matched_dir,
                       sprintf("matched_complete_%s_%s.rds", disease_a, disease_b))
  if (!file.exists(in_file)) { n_skip <- n_skip + 1L; next }

  res_pair <- tryCatch({

    obj     <- readRDS(in_file)
    matched <- obj$matched

    relevant_idp <- unique(c(matched$case_idp, matched$ctrl_idp))
    cohort_sub <- cohort_full[.(relevant_idp), on = "idp", nomatch = 0]
    cohort_sub <- cohort_sub[!(cod %in% c(disease_a, disease_b)),
                             .(first_dat = min(dat)), by = .(idp, cod)]

    pair_results <- vector("list", length(windows_cont))

    for (w in seq_along(windows_cont)) {
      wname <- names(windows_cont)[w]
      ts    <- windows_cont[[wname]][1]
      te    <- windows_cont[[wname]][2]

      case_dt <- unique(matched[, .(idp = case_idp, index_date = case_index_date,
                                     dat_B = case_dat_B, followup_end = case_followup_end)])
      case_dt[, in_risk := at_risk_start(index_date, dat_B, followup_end, ts, te)]
      case_risk <- case_dt[in_risk == TRUE]

      ctrl_dt <- matched[, .(idp = ctrl_idp, index_date = ctrl_index_date,
                              dat_B = ctrl_dat_B, followup_end = ctrl_followup_end)]
      ctrl_dt[, in_risk := at_risk_start(index_date, dat_B, followup_end, ts, te)]
      ctrl_risk <- ctrl_dt[in_risk == TRUE]

      if (nrow(case_risk) == 0 || nrow(ctrl_risk) == 0) next

      case_risk[, threshold := index_date + round(ts * 365.25)]
      ctrl_risk[, threshold := index_date + round(ts * 365.25)]

      case_risk[, n_comorbid := comorbidity_burden(cohort_sub, idp, threshold)]
      ctrl_risk[, n_comorbid := comorbidity_burden(cohort_sub, idp, threshold)]

      pair_results[[w]] <- data.table(
        disease_a = disease_a, disease_b = disease_b,
        window = wname, t_start = ts,
        n_at_risk_case = nrow(case_risk), n_at_risk_ctrl = nrow(ctrl_risk),
        mean_comorbidity_case = mean(case_risk$n_comorbid),
        mean_comorbidity_ctrl = mean(ctrl_risk$n_comorbid)
      )
    }

    rbindlist(pair_results)

  }, error = function(e) {
    cat(sprintf("  ERROR in %s_%s: %s\n", disease_a, disease_b, conditionMessage(e)))
    NULL
  })

  if (is.null(res_pair) || nrow(res_pair) == 0) {
    n_err <- n_err + 1L
  } else {
    all_results[[length(all_results) + 1]] <- res_pair
    n_ok <- n_ok + 1L
  }

  if (i %% 25 == 0) {
    mem_log(sprintf("%d/%d pairs processed (ok=%d, skip=%d, err=%d)",
                     i, nrow(pairs_tbl), n_ok, n_skip, n_err))
  }
}

combined <- rbindlist(all_results)
combined[, excess := mean_comorbidity_case - mean_comorbidity_ctrl]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
batch_id <- gsub("\\.txt$", "", basename(pairs_file))
out_file <- file.path(out_dir, sprintf("depletion_%s.csv", batch_id))
fwrite(combined, out_file)

elapsed <- round(as.numeric(difftime(Sys.time(), t_start_job, units = "mins")), 1)
cat(sprintf("\n=== Complete batch: %d ok, %d without matched_complete, %d with error (%.1f min) ===\n",
            n_ok, n_skip, n_err, elapsed))
cat(sprintf("Saved: %s\n", out_file))
