## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 11_future_a_candidates_all.R \
#        --data_dir ./CaseControlStudy/ --out_dir ./future_a_candidates/
library(data.table)

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- args[which(args == "--data_dir") + 1]
out_dir  <- args[which(args == "--out_dir")  + 1]

FOLLOWUP_YEARS <- 5

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))
t0 <- Sys.time()

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort <- readRDS(file.path(data_dir, "cohort.rds"))
setDT(cohort)
valid_diseases_A   <- readRDS(file.path(data_dir, "valid_diseases_A.rds"))
valid_diseases_any <- readRDS(file.path(data_dir, "valid_diseases_any.rds"))
mem_log(sprintf("cohort.rds loaded: %d rows, %d patients", nrow(cohort), uniqueN(cohort$idp)))

setorder(cohort, idp, dat)
first_idx <- cohort[, .I[1], by = idp]$V1
all_first <- cohort[first_idx, .(idp, index_date_first = dat, index_cod_first = cod,
                                 followup_end, followup_years_first = followup_years)]
rm(first_idx); gc()
mem_log(sprintf("all_first: %d patients", nrow(all_first)))

candidates <- all_first[index_cod_first %in% valid_diseases_any &
                        followup_years_first >= FOLLOWUP_YEARS]
rm(all_first); gc()

setkey(cohort, cod)

future_a_list <- vector("list", length(valid_diseases_A))

for (i in seq_along(valid_diseases_A)) {
  da_val <- valid_diseases_A[i]

  a_rows <- cohort[.(da_val), .(idp, dat), on = "cod"]
  if (nrow(a_rows) == 0 || is.na(a_rows$idp[1])) next

  a_first_date <- a_rows[, .(dat_A = min(dat)), by = idp]

  fa <- merge(candidates[index_cod_first != da_val, .(idp, index_date_first, followup_end)],
             a_first_date, by = "idp")
  fa <- fa[dat_A > index_date_first]

  if (nrow(fa) < 10) next

  fa[, disease_a := da_val]
  future_a_list[[i]] <- fa[, .(idp, disease_a, index_date_first, followup_end, dat_A)]

  if (i %% 100 == 0) {
    cat(sprintf("  ...%d/%d (%.1f min)\n", i, length(valid_diseases_A),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

future_a_candidates <- rbindlist(future_a_list)
saveRDS(future_a_candidates, file.path(out_dir, "future_a_candidates.rds"))

cat(sprintf("\Saved: %s\n", file.path(out_dir, "future_a_candidates.rds")))