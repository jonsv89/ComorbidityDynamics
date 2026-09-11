## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 05_depletion_footprint.R --disease_a E11 --disease_b I10 \
#              --data_dir ./CaseControlStudy/ --matched_dir ./matched/ \
#              --out_dir ./depletion/
library(data.table)

args        <- commandArgs(trailingOnly = TRUE)
disease_a   <- args[which(args == "--disease_a") + 1]
disease_b   <- args[which(args == "--disease_b") + 1]
data_dir    <- args[which(args == "--data_dir") + 1]
matched_dir <- args[which(args == "--matched_dir") + 1]
out_dir     <- args[which(args == "--out_dir") + 1]

mem_log <- function(msg) {
  mem_mb <- sum(gc()[, 2])
  cat(sprintf("[MEM %6.0f MB] %s\n", mem_mb, msg))
}

cat(sprintf("\n=== Depletion footprint %s -> %s ===\n", disease_a, disease_b))
mem_log("inicio")

in_file <- file.path(matched_dir,
                     sprintf("matched_complete_%s_%s.rds", disease_a, disease_b))
if (!file.exists(in_file)) {
  cat(sprintf("SKIP: no existe %s\n", in_file))
  quit(status = 0)
}

obj     <- readRDS(in_file)
matched <- obj$matched
mem_log(sprintf("pares cargados: %d filas", nrow(matched)))

cohort <- readRDS(file.path(data_dir, "cohort.rds"))

relevant_idp <- unique(c(matched$case_idp, matched$ctrl_idp))
cohort_sub <- cohort[idp %in% relevant_idp & !(cod %in% c(disease_a, disease_b)),
                     .(first_dat = min(dat)), by = .(idp, cod)]
rm(cohort); gc()
mem_log(sprintf("historial diagnostico (excl. A y B): %d filas, %d pacientes",
                nrow(cohort_sub), uniqueN(cohort_sub$idp)))

windows_cont <- list("1_2" = c(1, 2), "2_3" = c(2, 3), "3_4" = c(3, 4), "4_5" = c(4, 5))

at_risk_start <- function(index_date, dat_B, followup_end, t_start, t_end) {
  days_start   <- t_start * 365.25
  days_end     <- t_end   * 365.25
  had_B_before <- !is.na(dat_B) & as.numeric(dat_B - index_date) <= days_start
  has_followup <- as.numeric(followup_end - index_date) >= days_end
  !had_B_before & has_followup
}

comorbidity_burden <- function(idps, threshold_dates) {
  lut <- data.table(idp = idps, threshold = threshold_dates)
  m <- merge(cohort_sub, lut, by = "idp", allow.cartesian = TRUE)
  m <- m[first_dat <= threshold]
  burden <- m[, .(n_comorbid = uniqueN(cod)), by = idp]
  out <- merge(data.table(idp = idps), burden, by = "idp", all.x = TRUE)
  out[is.na(n_comorbid), n_comorbid := 0L]
  out$n_comorbid
}

results <- vector("list", length(windows_cont))

for (i in seq_along(windows_cont)) {
  wname <- names(windows_cont)[i]
  ts    <- windows_cont[[wname]][1]
  te    <- windows_cont[[wname]][2]

  case_dt <- unique(matched[, .(idp = case_idp, index_date = case_index_date,
                                 dat_B = case_dat_B,
                                 followup_end = case_followup_end)])
  case_dt[, in_risk := at_risk_start(index_date, dat_B, followup_end, ts, te)]
  case_risk <- case_dt[in_risk == TRUE]

  ctrl_dt <- matched[, .(idp = ctrl_idp, index_date = ctrl_index_date,
                          dat_B = ctrl_dat_B, followup_end = ctrl_followup_end)]
  ctrl_dt[, in_risk := at_risk_start(index_date, dat_B, followup_end, ts, te)]
  ctrl_risk <- ctrl_dt[in_risk == TRUE]

  if (nrow(case_risk) == 0 || nrow(ctrl_risk) == 0) {
    cat(sprintf("  %s: without patients at risk in any arm, jumping to the next\n", wname))
    next
  }

  case_risk[, threshold := index_date + round(ts * 365.25)]
  ctrl_risk[, threshold := index_date + round(ts * 365.25)]

  case_risk[, n_comorbid := comorbidity_burden(idp, threshold)]
  ctrl_risk[, n_comorbid := comorbidity_burden(idp, threshold)]

  results[[i]] <- data.table(
    window            = wname,
    t_start           = ts,
    arm               = c("case", "control"),
    n_at_risk         = c(nrow(case_risk), nrow(ctrl_risk)),
    mean_comorbidity  = c(mean(case_risk$n_comorbid), mean(ctrl_risk$n_comorbid)),
    median_comorbidity = c(median(case_risk$n_comorbid), median(ctrl_risk$n_comorbid))
  )

  cat(sprintf("  %s: n_casos=%d (media=%.2f) | n_ctrl=%d (media=%.2f) | exceso=%.2f\n",
              wname, nrow(case_risk), mean(case_risk$n_comorbid),
              nrow(ctrl_risk), mean(ctrl_risk$n_comorbid),
              mean(case_risk$n_comorbid) - mean(ctrl_risk$n_comorbid)))
}

footprint <- rbindlist(results)

cat("\nSummary: evolution of the prior comorbidity load\n")
print(dcast(footprint, window + t_start ~ arm, value.var = "mean_comorbidity"))

wide <- dcast(footprint, window + t_start ~ arm, value.var = "mean_comorbidity")
wide[, excess := case - control]
cat("\nComorbidity excess (cases - controls) by window:\n")
print(wide[, .(window, t_start, case, control, excess)])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("depletion_footprint_%s_%s.csv", disease_a, disease_b))
fwrite(footprint, out_file)
cat(sprintf("\nGuardado: %s\n", out_file))
