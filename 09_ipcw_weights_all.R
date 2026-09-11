## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 09_ipcw_weights_all.R \
#        --data_dir ./CaseControlStudy/ --out_dir ./ipcw_weights/
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
mem_log(sprintf("cohort.rds cargado: %d filas, %d pacientes", nrow(cohort), uniqueN(cohort$idp)))

setorder(cohort, idp, dat)
first_idx <- cohort[, .I[1], by = idp]$V1
all_first <- cohort[first_idx, .(idp, index_date_first = dat, index_cod_first = cod,
                                 followup_years_first = followup_years, sexe, rangos)]
rm(first_idx); gc()
mem_log(sprintf("all_first: %d pacientes", nrow(all_first)))

ctrl_pop <- all_first[index_cod_first %in% valid_diseases_any]
ctrl_pop[, eligible := followup_years_first >= FOLLOWUP_YEARS]
ctrl_pop[, rangos_f := factor(rangos)]
ctrl_pop[, sexe_f    := factor(sexe)]

ctrl_model <- glm(eligible ~ rangos_f + sexe_f, data = ctrl_pop, family = binomial)
ctrl_pop[, p_eligible := predict(ctrl_model, newdata = ctrl_pop, type = "response")]
ctrl_pop[, control_weight := ifelse(eligible, 1 / p_eligible, NA_real_)]

print(summary(ctrl_pop[eligible == TRUE, control_weight]))

control_weights <- ctrl_pop[eligible == TRUE, .(idp, control_weight)]
saveRDS(control_weights, file.path(out_dir, "control_weights.rds"))
cat(sprintf("Guardado: control_weights.rds (%d pacientes)\n", nrow(control_weights)))

rm(ctrl_pop); gc()
mem_log("controls weights saved")

setkey(cohort, cod)
case_weights_list <- vector("list", length(valid_diseases_A))

for (i in seq_along(valid_diseases_A)) {
  disease_a <- valid_diseases_A[i]

  a_rows <- cohort[.(disease_a), .(idp, dat, followup_years, rangos, sexe), on = "cod"]
  if (nrow(a_rows) == 0 || is.na(a_rows$idp[1])) next

  a_events <- a_rows[, .SD[which.min(dat)], by = idp]
  if (nrow(a_events) < 30) next

  a_events[, eligible := followup_years >= FOLLOWUP_YEARS]
  a_events[, rangos_f := factor(rangos)]
  a_events[, sexe_f    := factor(sexe)]

  if (uniqueN(a_events$eligible) < 2 || min(table(a_events$eligible)) < 10) next

  model <- tryCatch(glm(eligible ~ rangos_f + sexe_f, data = a_events, family = binomial),
                    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(model)) next

  a_events[, p_eligible := predict(model, newdata = a_events, type = "response")]
  a_events[, case_weight := ifelse(eligible, 1 / p_eligible, NA_real_)]

  case_weights_list[[i]] <- a_events[eligible == TRUE, .(idp, disease_a = disease_a, case_weight)]

  if (i %% 100 == 0) {
    cat(sprintf("  ...%d/%d (%.1f min)\n", i, length(valid_diseases_A),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

case_weights <- rbindlist(case_weights_list)
saveRDS(case_weights, file.path(out_dir, "case_weights.rds"))
cat("\nGlobal summary of cases' weights:\n")
print(summary(case_weights$case_weight))
