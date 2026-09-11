## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript cohort_attrition.R \
#        --raw_file ./Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt \
#        --data_dir ./CaseControlStudy/ --weights_dir ./ipcw_weights/ \
#        --out_dir ./attrition/
library(data.table)

args        <- commandArgs(trailingOnly = TRUE)
raw_file    <- args[which(args == "--raw_file")    + 1]
data_dir    <- args[which(args == "--data_dir")    + 1]
weights_dir <- args[which(args == "--weights_dir") + 1]
out_dir     <- args[which(args == "--out_dir")     + 1]

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))
t0 <- Sys.time()

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

STUDY_START  <- as.Date("2008-01-01")
STUDY_END    <- as.Date("2018-12-31")
MIN_PATIENTS <- 100L
FOLLOWUP_YEARS <- 5

attrition <- data.table(step = character(), n_patients = integer(), n_rows = integer())
add_step <- function(step, n_patients, n_rows) {
  attrition <<- rbind(attrition, data.table(step = step, n_patients = n_patients, n_rows = n_rows))
  cat(sprintf("  [%s] pacientes=%d, filas=%d\n", step, n_patients, n_rows))
}
dt <- fread(raw_file)
mem_log(sprintf("bruto cargado: %d filas", nrow(dt)))
add_step("0_bruto_sidiap", dt[, uniqueN(idp)], nrow(dt))

int_to_date <- function(x) as.Date(as.character(x), format = "%Y%m%d")
dt[, birth   := int_to_date(birth)]
dt[, sortida := int_to_date(sortida)]
dt[, dat     := int_to_date(dat)]

dt <- dt[dat >= STUDY_START & dat <= STUDY_END]
add_step("1_filtro_fechas_estudio", dt[, uniqueN(idp)], nrow(dt))

dt[, followup_end   := pmin(sortida, STUDY_END)]
dt[, followup_years := as.numeric(followup_end - dat) / 365.25]

N_COHORT <- dt[, uniqueN(idp)]
disease_counts <- dt[, .(N = .N), by = cod]
disease_counts[, prevalence := N / N_COHORT]
valid_diseases_any <- disease_counts[N >= MIN_PATIENTS, cod]

dt <- dt[cod %in% valid_diseases_any]
add_step("2_filtro_enfermedades_validas_N100", dt[, uniqueN(idp)], nrow(dt))

mem_log("cohorte tras pasos 0-2 (deberia coincidir con cohort.rds)")
rm(dt); gc()

control_weights <- readRDS(file.path(weights_dir, "control_weights.rds"))
n_control_eligible <- uniqueN(control_weights$idp)
add_step("3_elegible_como_control_5y_followup", n_control_eligible, NA_integer_)
rm(control_weights); gc()

case_weights <- readRDS(file.path(weights_dir, "case_weights.rds"))
n_case_eligible <- uniqueN(case_weights$idp)
add_step("4_elegible_como_caso_alguna_enfermedad_indice", n_case_eligible, NA_integer_)

cat(sprintf("\n=== Completed (%.1f min) ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

attrition[, pct_of_step0 := round(100 * n_patients / attrition$n_patients[1], 2)]
attrition[, pct_lost_this_step := round(100 * (1 - n_patients / shift(n_patients)), 2)]

print(attrition)

fwrite(attrition, file.path(out_dir, "cohort_attrition_funnel.csv"))
cat("\nGuardado: cohort_attrition_funnel.csv\n")