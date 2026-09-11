## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
args<-"Power_Anchored_60pct_v3"
if(args[1] == "Power_Anchored_60pct_v3"){
  library(data.table)
  set.seed(42)
  cat("\nSTEP A: Loading cumulative windows\n")

  windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")

  win_data <- lapply(windows_cumul, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt <- fread(f)
    dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk) & is.finite(SE_RR_shrunk)]
    dt[, pair := paste(disease_a, disease_b, sep = "_")]
    dt[, sig  := lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, .(pair, RR_shrunk, log_RR_shrunk, SE_RR_shrunk, cases_event, ctrl_events, sig)]
  })
  names(win_data) <- windows_cumul

  cat("STEP B/C: Building taxonomy (exact replica of Temporal_robustness_of_comorbidities)\n")

  windows_9   <- c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  window_type <- c(rep("cumulative",5), rep("conditional",4))

  all_9 <- rbindlist(lapply(seq_along(windows_9), function(i) {
    w  <- windows_9[i]
    wt <- window_type[i]
    f  <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt <- fread(f)
    dt <- dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, `:=`(window = w, window_type = wt, pair = paste(disease_a, disease_b, sep = "_"))]
    dt
  }))

  pair_windows <- all_9[, .(
    in_0_1 = any(window == "0_1"), in_0_2 = any(window == "0_2"), in_0_3 = any(window == "0_3"),
    in_0_4 = any(window == "0_4"), in_0_5 = any(window == "0_5"),
    in_1_2 = any(window == "1_2"), in_2_3 = any(window == "2_3"),
    in_3_4 = any(window == "3_4"), in_4_5 = any(window == "4_5"),
    n_cumul = sum(window_type == "cumulative"),
    n_cond  = sum(window_type == "conditional")
  ), by = pair]

  cat(sprintf("  Total significant in >=1 of 9 windows: %d  [manuscript: 146,591]\n", nrow(pair_windows)))

  pair_windows[, bio_class := fcase(
    n_cumul == 5 & n_cond == 4, "Omnipresent (all 9 windows)",
    n_cumul == 5 & n_cond > 0,  "Persistent cumulative + partial conditional",
    n_cumul == 5 & n_cond == 0, "Persistent cumulative only (no conditional risk)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond == 0, "Long-window only (late cumulative, no early signal)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond >  0, "Late cumulative + conditional risk",
    n_cumul == 0 & n_cond > 0,  "Conditional only (no cumulative signal)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond == 0, "Early transient (short window only)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond >  0, "Early cumulative + late conditional",
    rep(TRUE, .N), "Mixed/partial"
  )]

  cat("\n--- SANITY CHECK: reconstructed taxonomy vs manuscript ---\n")
  cat("Manuscript: Omnipresent=45,495 | Persistent cumul+partial cond=7,908 |\n")
  cat("            Persistent cumul only=2,445 | Long-window only=72,495 |\n")
  cat("            Late cumul+conditional=15,645 | Total sig(>=1 window)=146,591\n\n")
  print(pair_windows[, .N, by = bio_class][order(-N)])

  pair_windows[, late_emerging := bio_class %in% c(
    "Long-window only (late cumulative, no early signal)",
    "Late cumulative + conditional risk"
  )]

  n_late <- sum(pair_windows$late_emerging)
  cat(sprintf("\nLate-emerging (reconstructed): %d (%.1f%% of %d)  [manuscript: 88,140 / 60.1%%]\n",
              n_late, 100*n_late/nrow(pair_windows), nrow(pair_windows)))

  cat("\n*** STOP: compare every number above against the manuscript before reading below ***\n\n")

  fwrite(pair_windows, "ManuscriptFiles/Results/Power_anchored_taxonomy_v3.txt", sep = "\t", quote = FALSE)

  late_pairs <- pair_windows[late_emerging == TRUE, pair]

  cat("PART 1: Why is each late-emerging pair not significant at 0-1?\n")

  d01 <- win_data[["0_1"]][pair %in% late_pairs, .(pair, cases_event_01 = cases_event)]
  part1 <- merge(data.table(pair = late_pairs), d01, by = "pair", all.x = TRUE)
  part1[is.na(cases_event_01), cases_event_01 := 0L]
  part1[, reason := fifelse(cases_event_01 < 100, "Insufficient events (n<100)", "Enough events, effect/precision criterion not met")]

  p1 <- part1[, .N, by = reason][order(-N)]
  p1[, pct := round(100 * N / sum(N), 1)]
  print(p1)

  fwrite(p1, "ManuscriptFiles/Results/Power_anchored_part1_reasons_v2.txt", sep = "\t", quote = FALSE)

  cat("\nPART 2: Would the FINAL (0-5y) effect size survive at n=100 events?\n")

  d05 <- win_data[["0_5"]][pair %in% late_pairs,
                            .(pair, RR_05 = RR_shrunk, a_05 = cases_event, c_05 = ctrl_events)]
  d05 <- d05[!is.na(a_05) & !is.na(c_05) & a_05 >= 100 & c_05 > 0]
  cat(sprintf("  Pairs entering Part 2: %d of %d\n", nrow(d05), length(late_pairs)))

  d05[, thin_frac := pmin(1, 100 / a_05)]

  N_REP   <- 500
  n_pairs <- nrow(d05)
  times_sig <- integer(n_pairs)
  a05v <- d05$a_05; c05v <- d05$c_05; RRv <- d05$RR_05; fv <- d05$thin_frac

  for (r in seq_len(N_REP)) {
    a_thin <- rbinom(n_pairs, a05v, fv)
    c_thin <- rbinom(n_pairs, c05v, fv)
    valid  <- a_thin >= 1 & c_thin >= 1

    RR_thin <- rep(NA_real_, n_pairs)
    RR_thin[valid] <- RRv[valid] * (a_thin[valid] * c05v[valid]) / (c_thin[valid] * a05v[valid])

    SE_thin <- rep(NA_real_, n_pairs)
    SE_thin[valid] <- sqrt(1/a_thin[valid] + 1/c_thin[valid])

    CI_low_thin <- rep(NA_real_, n_pairs)
    CI_low_thin[valid] <- exp(log(RR_thin[valid]) - 1.96 * SE_thin[valid])

    sig_thin <- valid & a_thin >= 100 & !is.na(CI_low_thin) & CI_low_thin >= 1.01
    times_sig <- times_sig + as.integer(sig_thin)

    if (r %% 100 == 0) cat(sprintf("    ...replicate %d/%d\n", r, N_REP))
  }

  d05[, pct_still_significant := 100 * times_sig / N_REP]

  cat(sprintf("\nAcross %d late-emerging pairs (thinned to ~100 events):\n", n_pairs))
  cat(sprintf("  Median %% of replicates still significant: %.1f%%\n",
              median(d05$pct_still_significant)))
  cat(sprintf("  Pairs still significant in >=80%% of replicates: %d (%.1f%%)\n",
              sum(d05$pct_still_significant >= 80), 100*mean(d05$pct_still_significant >= 80)))
  cat(sprintf("  Pairs still significant in <20%% of replicates:  %d (%.1f%%)\n",
              sum(d05$pct_still_significant < 20), 100*mean(d05$pct_still_significant < 20)))

  cat("\nDistribution summary:\n")
  print(summary(d05$pct_still_significant))

  fwrite(d05, "ManuscriptFiles/Results/Power_anchored_part2_thinning_v2.txt", sep = "\t", quote = FALSE)
  cat("Saved: Power_anchored_part2_thinning_v2.txt\n")

  cat("\nPART 3: Trajectory z-test (log RR at 0-1 vs 0-5)\n")

  d01_full <- win_data[["0_1"]][pair %in% late_pairs,
                                 .(pair, logRR_01 = log_RR_shrunk, SE_01 = SE_RR_shrunk,
                                   a_01 = cases_event)]
  d05_full <- win_data[["0_5"]][pair %in% late_pairs,
                                 .(pair, logRR_05 = log_RR_shrunk, SE_05 = SE_RR_shrunk)]

  traj <- merge(d01_full, d05_full, by = "pair")
  traj <- traj[a_01 >= 5]

  traj[, z   := (logRR_05 - logRR_01) / sqrt(SE_01^2 + SE_05^2)]
  traj[, p   := 2 * pnorm(-abs(z))]
  traj[, fdr := p.adjust(p, method = "BH")]
  traj[, grew_significantly := !is.na(fdr) & fdr < 0.05 & logRR_05 > logRR_01]

  cat(sprintf("  Late-emerging pairs with a usable 0-1 estimate (a_01>=5): %d of %d\n",
              nrow(traj), length(late_pairs)))
  cat(sprintf("  ...log(RR) grew significantly from 0-1 to 0-5 (FDR<0.05): %d (%.1f%%)\n",
              sum(traj$grew_significantly), 100*mean(traj$grew_significantly)))
  cat(sprintf("  ...statistically indistinguishable from the 0-1 estimate: %d (%.1f%%)\n",
              sum(!traj$grew_significantly), 100*mean(!traj$grew_significantly)))

  fwrite(traj, "ManuscriptFiles/Results/Power_anchored_part3_trajectory_v2.txt", sep = "\t", quote = FALSE)
  cat("Saved: Power_anchored_part3_trajectory_v2.txt\n")
}