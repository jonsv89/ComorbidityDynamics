## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
# Use: Rscript 10_ipcw_reweighted_rr_batch.R \
#        --pairs_file ./batches/pairs_batch_0001.txt \
#        --matched_dir ./matched/ --weights_dir ./ipcw_weights/ \
#        --out_dir ./ipcw_rr/
library(data.table)

args         <- commandArgs(trailingOnly = TRUE)
pairs_file   <- args[which(args == "--pairs_file")   + 1]
matched_dir  <- args[which(args == "--matched_dir")  + 1]
weights_dir  <- args[which(args == "--weights_dir")  + 1]
out_dir      <- args[which(args == "--out_dir")      + 1]

mem_log <- function(msg) cat(sprintf("[MEM %6.0f MB] %s\n", sum(gc()[, 2]), msg))
t0 <- Sys.time()

cat(sprintf("\n=== IPCW-reweighted RR BATCH: %s ===\n", pairs_file))
mem_log("inicio")

pairs_tbl <- fread(pairs_file, header = FALSE, col.names = c("disease_a", "disease_b"))

control_weights <- readRDS(file.path(weights_dir, "control_weights.rds"))
setkey(control_weights, idp)

case_weights_all <- readRDS(file.path(weights_dir, "case_weights.rds"))
setkey(case_weights_all, disease_a, idp)
mem_log(sprintf("Pesos cargados: %d control, %d caso-enfermedad",
                nrow(control_weights), nrow(case_weights_all)))

compute_rr <- function(case_idp, case_has_B, ctrl_idp, ctrl_has_B,
                       case_w = NULL, ctrl_w = NULL) {
  if (is.null(case_w)) case_w <- rep(1, length(case_idp))
  if (is.null(ctrl_w)) ctrl_w <- rep(1, length(ctrl_idp))

  n1 <- sum(case_w)
  a  <- sum(case_w[case_has_B])
  n0 <- sum(ctrl_w)
  c  <- sum(ctrl_w[ctrl_has_B])

  if (a == 0 || c == 0 || n1 == 0 || n0 == 0) return(c(RR = NA_real_, a = a, n1 = n1, c = c, n0 = n0))

  rr <- (a / n1) / (c / n0)
  c(RR = rr, a = a, n1 = n1, c = c, n0 = n0)
}

results <- vector("list", nrow(pairs_tbl))
n_ok <- 0L; n_skip <- 0L
n_diag_printed <- 0L
MAX_DIAG <- 3L

for (i in seq_len(nrow(pairs_tbl))) {
  da_val <- pairs_tbl$disease_a[i]
  db_val <- pairs_tbl$disease_b[i]

  in_file <- file.path(matched_dir, sprintf("matched_complete_%s_%s.rds", da_val, db_val))
  if (!file.exists(in_file)) { n_skip <- n_skip + 1L; next }

  show_diag <- n_diag_printed < MAX_DIAG

  res <- tryCatch({
    obj     <- readRDS(in_file)
    matched <- obj$matched
    setDT(matched)

    case_has_B <- !is.na(matched$case_dat_B)
    ctrl_has_B <- !is.na(matched$ctrl_dat_B)

    naive <- compute_rr(matched$case_idp, case_has_B, matched$ctrl_idp, ctrl_has_B)

    case_w_sub <- case_weights_all[disease_a == da_val]
    setkey(case_w_sub, idp)
    cw <- case_w_sub[.(matched$case_idp), case_weight]
    ow <- control_weights[.(matched$ctrl_idp), control_weight, on = "idp"]

    if (show_diag) {
      n_diag_printed <<- n_diag_printed + 1L
    }

    if (any(is.na(cw)) || any(is.na(ow))) {
      list(naive = naive, weighted = NULL, missing_weights = TRUE)
    } else {
      weighted <- compute_rr(matched$case_idp, case_has_B, matched$ctrl_idp, ctrl_has_B,
                             case_w = cw, ctrl_w = ow)
      list(naive = naive, weighted = weighted, missing_weights = FALSE)
    }
  }, error = function(e) {
    if (show_diag) cat(sprintf("\n--- Real error in pair %s -> %s: %s ---\n",
                                da_val, db_val, conditionMessage(e)))
    NULL
  })

  if (is.null(res) || res$missing_weights || is.null(res$weighted) ||
      is.na(res$naive["RR"]) || is.na(res$weighted["RR"])) {
    n_skip <- n_skip + 1L
    next
  }

  results[[i]] <- data.table(
    disease_a = da_val, disease_b = db_val,
    RR_naive = res$naive["RR"], RR_weighted = res$weighted["RR"],
    n1_naive = res$naive["n1"], n1_weighted_eff = res$weighted["n1"],
    n0_naive = res$naive["n0"], n0_weighted_eff = res$weighted["n0"]
  )
  n_ok <- n_ok + 1L
}

combined <- rbindlist(results)

if (nrow(combined) == 0 || !("RR_weighted" %in% names(combined))) {
  quit(status = 1)
}

combined[, pct_change := 100 * (RR_weighted - RR_naive) / RR_naive]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
batch_id <- gsub("\\.txt$", "", basename(pairs_file))
out_file <- file.path(out_dir, sprintf("ipcw_rr_%s.csv", batch_id))
fwrite(combined, out_file)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("\n=== Completed: %d ok, %d without data/weights (%.1f min) ===\n", n_ok, n_skip, elapsed))