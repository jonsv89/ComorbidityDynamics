## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript matching_success_and_balance.R \
#        --matched_dir ./matched/ --weights_dir ./ipcw_weights/ \
#        --out_dir ./matching_balance/
library(data.table)

args        <- commandArgs(trailingOnly = TRUE)
matched_dir <- args[which(args == "--matched_dir") + 1]
weights_dir <- args[which(args == "--weights_dir") + 1]
out_dir     <- args[which(args == "--out_dir")     + 1]

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))
t0 <- Sys.time()

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

case_weights <- readRDS(file.path(weights_dir, "case_weights.rds"))
setDT(case_weights)
eligible_counts <- case_weights[, .(n_eligible = uniqueN(idp)), by = disease_a]
mem_log(sprintf("Elegible cases by loaded disease: %d diseases", nrow(eligible_counts)))

valid_A <- readRDS("CaseControlStudy/valid_diseases_A.rds")

all_files <- list.files(matched_dir, pattern = "^matched_complete_", full.names = FALSE)

results <- vector("list", length(valid_A))
temporal_gaps_list <- vector("list", length(valid_A))

for (i in seq_along(valid_A)) {
  da_val <- valid_A[i]
  candidate_files <- all_files[startsWith(all_files, paste0("matched_complete_", da_val, "_"))]
  if (length(candidate_files) == 0) next

  obj <- tryCatch(readRDS(file.path(matched_dir, candidate_files[1])), error = function(e) NULL)
  if (is.null(obj)) next
  m <- obj$matched
  setDT(m)

  n_cases_matched <- uniqueN(m$case_idp)
  controls_per_case <- m[, .N, by = case_idp]

  results[[i]] <- data.table(
    disease_a = da_val,
    n_cases_matched = n_cases_matched,
    median_controls_per_case = median(controls_per_case$N),
    pct_cases_5plus_controls = 100 * mean(controls_per_case$N >= 5),
    pct_cases_lt5_controls   = 100 * mean(controls_per_case$N < 5),
    min_controls = min(controls_per_case$N),
    max_controls = max(controls_per_case$N)
  )

  gap_days <- as.numeric(m$ctrl_index_date - m$case_index_date)
  temporal_gaps_list[[i]] <- data.table(disease_a = da_val, gap_days = gap_days)

  if (i %% 200 == 0) {
    cat(sprintf("  ...%d/%d (%.1f min)\n", i, length(valid_A),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

matching_summary <- rbindlist(results)
matching_summary <- merge(matching_summary, eligible_counts, by = "disease_a", all.x = TRUE)
matching_summary[, matching_yield_pct := 100 * n_cases_matched / n_eligible]

temporal_gaps <- rbindlist(temporal_gaps_list)
fwrite(matching_summary, file.path(out_dir, "matching_success_by_disease.csv"))
fwrite(temporal_gaps[, .(median_gap = median(gap_days), n = .N), by = disease_a], file.path(out_dir, "temporal_balance_by_disease.csv"))