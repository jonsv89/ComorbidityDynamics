## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
args<-"Clock_Axes_vs_Prevalence"
if(args[1] == "Clock_Axes_vs_Prevalence"){

  library(data.table)

  code <- c(paste0("A0",0:9), paste0("A",10:99),
            paste0("B0",0:9), paste0("B",10:99),
            paste0("C0",0:9), paste0("C",10:99),
            paste0("D0",0:9), paste0("D",10:48), paste0("D",50:89),
            paste0("E0",0:9), paste0("E",10:99),
            paste0("F0",0:9), paste0("F",10:99),
            paste0("G0",0:9), paste0("G",10:99),
            paste0("H0",0:9), paste0("H",10:59), paste0("H",60:95),
            paste0("I0",0:9), paste0("I",10:99),
            paste0("J0",0:9), paste0("J",10:99),
            paste0("K0",0:9), paste0("K",10:93),
            paste0("L0",0:9), paste0("L",10:99),
            paste0("M0",0:9), paste0("M",10:99),
            paste0("N0",0:9), paste0("N",10:99),
            paste0("O0",0:9), paste0("O",10:99),
            paste0("P0",0:9), paste0("P",10:96),
            paste0("Q0",0:9), paste0("Q",10:99),
            paste0("R0",0:9), paste0("R",10:99),
            paste0("S0",0:9), paste0("S",10:99),
            paste0("T0",0:9), paste0("T",10:98))

  cate <- c(rep("I",200), rep("II",149), rep("III",40),
            rep("IV",100), rep("V",100), rep("VI",100),
            rep("VII",60), rep("VIII",36), rep("IX",100),
            rep("X",100), rep("XI",94), rep("XII",100),
            rep("XIII",100), rep("XIV",100), rep("XV",100),
            rep("XVI",97), rep("XVII",100), rep("XVIII",100),
            rep("XIX",199))

  catname <- c("Infectious","Neoplasms","Blood/Immune",
               "Endocrine/Metabolic","Mental","Nervous",
               "Eye","Ear","Circulatory","Respiratory",
               "Digestive","Skin","Musculoskeletal",
               "Genitourinary","Pregnancy","Perinatal",
               "Congenital","Symptoms","Injury")
  names(catname) <- unique(cate)
  names(cate)    <- code

  SIG_FILTER <- quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)

  cat("\nSTEP A: Rebuilding `slopes` and 0-5y event counts\n")

  windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")

  all_cumul <- rbindlist(lapply(windows_cumul, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt <- fread(f)
    dt <- dt[eval(SIG_FILTER)]
    dt[, `:=`(window     = w,
              window_num = as.integer(gsub("0_","",w)),
              pair       = paste(disease_a, disease_b, sep = "_"))]
    dt
  }))

  pairs_per_window <- all_cumul[, .(n_windows = uniqueN(window)), by = pair]
  persistent_pairs <- pairs_per_window[n_windows == 5, pair]

  rr_traj <- all_cumul[pair %in% persistent_pairs,
                       .(pair, window_num, RR_shrunk, SE_RR_shrunk, disease_a, disease_b)]
  rr_traj[, catA := catname[cate[disease_a]]]
  rr_traj[, catB := catname[cate[disease_b]]]

  slopes <- rr_traj[, {
    x <- window_num; y <- RR_shrunk; se <- SE_RR_shrunk
    if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
      .(slope = NA_real_)
    } else {
      fit <- lm(y ~ x, weights = 1/se^2)
      .(slope = coef(fit)["x"])
    }
  }, by = .(pair, catA, catB)]

  f05 <- "Results/shrinkage_rr_events/RR_contingency_both_0_5_shrunk.txt"
  d05 <- fread(f05)
  d05 <- d05[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
  d05[, pair := paste(disease_a, disease_b, sep = "_")]
  d05 <- d05[, .(pair, cases_event_05 = cases_event, ctrl_events_05 = ctrl_events)]

  slopes <- merge(slopes, d05, by = "pair")
  cat(sprintf("  Persistent pairs with slope + 0-5y proxies: %d\n", nrow(slopes)))

  cat("\nSTEP B: Rebuilding `cont_wide` and joining 0-5y event counts\n")

  windows_cont <- c("1_2","2_3","3_4","4_5")

  all_cont <- rbindlist(lapply(windows_cont, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt <- fread(f)
    dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
    dt[, `:=`(window     = w,
              window_num = as.integer(gsub("_.+","",w)),
              pair       = paste(disease_a, disease_b, sep = "_"),
              sig        = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
    dt[, catA := catname[cate[disease_a]]]
    dt
  }))

  pairs_cont4 <- all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]

  cont_wide <- data.table::dcast(
    all_cont[pair %in% pairs_cont4, .(pair, window, sig, catA)],
    pair + catA ~ window, value.var = "sig"
  )
  setnames(cont_wide, old = windows_cont, new = paste0("sig_", windows_cont))

  cont_wide[, pattern := fcase(
    (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == FALSE & sig_4_5 == FALSE), "Early only",
    (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == TRUE  | sig_4_5 == TRUE),  "Persistent conditional risk",
    (sig_1_2 == FALSE & sig_2_3 == FALSE) & (sig_3_4 == TRUE | sig_4_5 == TRUE), "Late emerging risk",
    rep(TRUE, .N), "Not significant"
  )]

  cont_wide[, is_persistent := pattern == "Persistent conditional risk"]

  cont_wide <- merge(cont_wide, d05, by = "pair")
  cat(sprintf("  Pairs in cont_wide (all 4 conditional windows) with 0-5y proxies: %d\n", nrow(cont_wide)))

  clock_2d <- merge(
    slopes[!is.na(slope) & !is.na(catA), .(med_slope = median(slope)), by = catA],
    cont_wide[!is.na(catA), .(pct_persist = mean(is_persistent) * 100), by = catA],
    by = "catA"
  )
  
  cat("\nSTEP D: Attenuation axis vs. background incidence\n")

  sd <- slopes[!is.na(slope) & !is.na(catA) & ctrl_events_05 > 0]
  sd[, log_incidence_proxy := log(ctrl_events_05)]

  cat_slope_inc <- sd[, .(med_slope = median(slope),
                           med_log_incidence = median(log_incidence_proxy)), by = catA]

  ct <- cor.test(cat_slope_inc$med_slope, cat_slope_inc$med_log_incidence, method = "pearson")
  cat(sprintf("Category-level: slope vs. log(background incidence proxy)\n"))
  cat(sprintf("  Pearson r = %.3f, 95%% CI [%.3f, %.3f], p = %.3f (n=%d categories)\n",
              ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value, nrow(cat_slope_inc)))

  m_reduced <- lm(slope ~ log_incidence_proxy, data = sd)
  m_full    <- lm(slope ~ catA + log_incidence_proxy, data = sd)
  a <- anova(m_reduced, m_full)
  print(a)

  r2_reduced <- summary(m_reduced)$r.squared
  r2_full    <- summary(m_full)$r.squared
  cat(sprintf("\n  R^2 (incidence only):            %.4f\n", r2_reduced))
  cat(sprintf("  R^2 (category + incidence):      %.4f\n", r2_full))
  cat(sprintf("  Incremental R^2 from category:   %.4f\n", r2_full - r2_reduced))
  cat(sprintf("  F-test for category effect after adjustment: F=%.2f, p=%s\n",
              a$F[2], format.pval(a[["Pr(>F)"]][2], digits = 3)))

  fwrite(cat_slope_inc, "ManuscriptFiles/Results/Clock_slope_vs_incidence_by_category.txt",
         sep = "\t", quote = FALSE)

  cat("\nSTEP E: Persistent-risk axis vs. sample size and background incidence\n")

  cw <- cont_wide[!is.na(catA) & cases_event_05 > 0 & ctrl_events_05 > 0]
  cw[, log_exposed_n     := log(cases_event_05)]
  cw[, log_incidence_proxy := log(ctrl_events_05)]

  cat_persist_conf <- cw[, .(
    pct_persist             = mean(is_persistent) * 100,
    median_exposed_n        = median(as.numeric(cases_event_05)),
    median_event_count      = median(as.numeric(cases_event_05)),
    median_incidence_proxy  = median(as.numeric(ctrl_events_05))
  ), by = catA]

  cat("\nCategory-level correlations (pct_persist vs. confounders):\n")
  for (v in c("median_exposed_n", "median_incidence_proxy")) {
    ctv <- cor.test(cat_persist_conf$pct_persist, log(cat_persist_conf[[v]]), method = "pearson")
    cat(sprintf("  pct_persist vs. log(%s): r = %.3f, 95%% CI [%.3f, %.3f], p = %.3f\n",
                v, ctv$estimate, ctv$conf.int[1], ctv$conf.int[2], ctv$p.value))
  }

  g_reduced <- glm(is_persistent ~ log_exposed_n + log_incidence_proxy,
                    data = cw, family = binomial)
  g_full    <- glm(is_persistent ~ catA + log_exposed_n + log_incidence_proxy,
                    data = cw, family = binomial)
  lrt <- anova(g_reduced, g_full, test = "Chisq")
  print(lrt)

  cat(sprintf("\n  Residual deviance (confounders only): %.0f\n", g_reduced$deviance))
  cat(sprintf("  Residual deviance (category + confounders): %.0f\n", g_full$deviance))
  cat(sprintf("  Likelihood-ratio test for category effect: Chisq=%.1f, df=%d, p=%s\n",
              lrt$Deviance[2], lrt$Df[2], format.pval(lrt[["Pr(>Chi)"]][2], digits = 3)))

  fwrite(cat_persist_conf, "ManuscriptFiles/Results/Clock_persist_vs_confounders_by_category.txt",
         sep = "\t", quote = FALSE)

  cat("\nSTEP F: Adjusted (residual) persistent-risk ranking by category\n")

  cw[, predicted_prob := predict(g_reduced, newdata = cw, type = "response")]

  cat_adjusted <- cw[, .(
    n_pairs      = .N,
    observed_pct = mean(is_persistent) * 100,
    expected_pct = mean(predicted_prob) * 100
  ), by = catA]

  cat_adjusted[, residual_pct := observed_pct - expected_pct]

  cat_adjusted[, se_pct := sqrt(expected_pct/100 * (1 - expected_pct/100) / n_pairs) * 100]
  cat_adjusted[, z   := residual_pct / se_pct]
  cat_adjusted[, p   := 2 * pnorm(-abs(z))]
  cat_adjusted[, fdr := p.adjust(p, method = "BH")]

  cat_adjusted <- cat_adjusted[order(-residual_pct)]
  print(cat_adjusted[, .(catA, n_pairs,
                          observed_pct = round(observed_pct, 1),
                          expected_pct = round(expected_pct, 1),
                          residual_pct = round(residual_pct, 1),
                          fdr = round(fdr, 4))])

  cat(sprintf("\nAbove expectation (FDR<0.05): %s\n",
              paste(cat_adjusted[residual_pct > 0 & fdr < 0.05, catA], collapse = ", ")))
  cat(sprintf("Below expectation (FDR<0.05): %s\n",
              paste(cat_adjusted[residual_pct < 0 & fdr < 0.05, catA], collapse = ", ")))

  fwrite(cat_adjusted, "ManuscriptFiles/Results/Clock_persist_adjusted_residuals.txt",
         sep = "\t", quote = FALSE)
  cat("Saved: Clock_persist_adjusted_residuals.txt\n")

  cat("\nSTEP G: Adjusted (residual) slope ranking by category\n")

  sd[, predicted_slope := predict(m_reduced, newdata = sd)]
  sd[, residual_slope  := slope - predicted_slope]

  cat_slope_adjusted <- sd[, .(
    n_pairs             = .N,
    observed_mean_slope = mean(slope),
    expected_mean_slope = mean(predicted_slope),
    mean_residual        = mean(residual_slope),
    se_residual          = stats::sd(residual_slope) / sqrt(.N)
  ), by = catA]

  cat_slope_adjusted[, t_stat := mean_residual / se_residual]
  cat_slope_adjusted[, p      := 2 * pt(-abs(t_stat), df = pmax(n_pairs - 1, 1))]
  cat_slope_adjusted[, fdr    := p.adjust(p, method = "BH")]

  cat_slope_adjusted <- cat_slope_adjusted[order(mean_residual)]
  print(cat_slope_adjusted[, .(catA, n_pairs,
                                observed_mean_slope = round(observed_mean_slope, 4),
                                expected_mean_slope = round(expected_mean_slope, 4),
                                mean_residual        = round(mean_residual, 4),
                                fdr = round(fdr, 4))])

  cat(sprintf("\nFaster attenuation than expected (more negative residual, FDR<0.05): %s\n",
              paste(cat_slope_adjusted[mean_residual < 0 & fdr < 0.05, catA], collapse = ", ")))
  cat(sprintf("Slower attenuation than expected (less negative/positive residual, FDR<0.05): %s\n",
              paste(cat_slope_adjusted[mean_residual > 0 & fdr < 0.05, catA], collapse = ", ")))

  fwrite(cat_slope_adjusted, "ManuscriptFiles/Results/Clock_slope_adjusted_residuals.txt",
         sep = "\t", quote = FALSE)
  cat("Saved: Clock_slope_adjusted_residuals.txt\n")
}