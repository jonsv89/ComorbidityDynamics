## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
args<-"Cluster_Validity_Biological_Clock"
if(args[1] == "Cluster_Validity_Biological_Clock"){
  library(data.table)
  library(cluster)
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

  cat("\nSTEP A: Rebuilding `slopes` from cumulative windows\n")

  windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")

  all_cumul <- rbindlist(lapply(windows_cumul, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) { cat("  Missing:", f, "\n"); return(NULL) }
    dt <- fread(f)
    dt <- dt[eval(SIG_FILTER)]
    dt[, `:=`(window = w,
              window_num = as.integer(gsub("0_","",w)),
              pair = paste(disease_a, disease_b, sep = "_"))]
    dt
  }))

  pairs_per_window <- all_cumul[, .(n_windows = uniqueN(window)), by = pair]
  persistent_pairs <- pairs_per_window[n_windows == 5, pair]

  rr_traj <- all_cumul[pair %in% persistent_pairs,
                       .(pair, window_num, RR_shrunk, SE_RR_shrunk,
                         disease_a, disease_b)]
  rr_traj[, catA := catname[cate[disease_a]]]
  rr_traj[, catB := catname[cate[disease_b]]]

  slopes <- rr_traj[, {
    x  <- window_num
    y  <- RR_shrunk
    se <- SE_RR_shrunk
    if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
      .(slope = NA_real_)
    } else {
      w_fit <- 1 / se^2
      fit   <- lm(y ~ x, weights = w_fit)
      .(slope = coef(fit)["x"])
    }
  }, by = .(pair, catA, catB)]

  traj_by_cat <- slopes[!is.na(catA) & !is.na(slope), .(
    n_total    = .N,
    med_slope  = median(slope)
  ), by = catA][order(med_slope)]

  cat("\nSTEP B: Rebuilding `cont_wide` from conditional windows\n")

  windows_cont <- c("1_2","2_3","3_4","4_5")

  all_cont <- rbindlist(lapply(windows_cont, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
    dt[, `:=`(window     = w,
              window_num = as.integer(gsub("_.+","",w)),
              pair       = paste(disease_a, disease_b, sep = "_"),
              sig        = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
    dt[, catA := catname[cate[disease_a]]]
    dt[, catB := catname[cate[disease_b]]]
    dt
  }))

  pairs_cont4 <- all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]

  cont_wide <- data.table::dcast(
    all_cont[pair %in% pairs_cont4,
             .(pair, window, sig, RR_shrunk, catA, catB)],
    pair + catA + catB ~ window,
    value.var = c("sig", "RR_shrunk")
  )

  cont_wide[, pattern := fcase(
    (sig_1_2 == TRUE | sig_2_3 == TRUE) &
      (sig_3_4 == FALSE & sig_4_5 == FALSE),
    "Early only (risk fades)",
    (sig_1_2 == TRUE | sig_2_3 == TRUE) &
      (sig_3_4 == TRUE  | sig_4_5 == TRUE),
    "Persistent conditional risk",
    (sig_1_2 == FALSE & sig_2_3 == FALSE) &
      (sig_3_4 == TRUE  | sig_4_5 == TRUE),
    "Late emerging risk",
    rep(TRUE, .N),
    "Not significant"
  )]

  pat_by_cat <- cont_wide[!is.na(catA) & !is.na(pattern), .(
    n_total     = .N,
    pct_persist = mean(pattern == "Persistent conditional risk") * 100
  ), by = catA]

  cat("\nSTEP C: Rebuilding clock_2d\n")

  clock_2d <- merge(
    traj_by_cat[, .(catA, med_slope, n_total)],
    pat_by_cat[,  .(catA, pct_persist)],
    by = "catA"
  )

  x_mid <- median(clock_2d$med_slope)
  y_mid <- median(clock_2d$pct_persist)

  clock_2d[, quadrant := fcase(
    med_slope <= x_mid & pct_persist >= y_mid, "Biphasic",
    med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive",
    med_slope <= x_mid & pct_persist <  y_mid, "Purely episodic",
    rep(TRUE, .N),                              "Chronic stable"
  )]

  cat("\n--- clock_2d (should match your validated STEP C output) ---\n")
  print(clock_2d[order(-pct_persist), .(catA, med_slope = round(med_slope,4),
                                          pct_persist = round(pct_persist,1),
                                          quadrant)])

  cat("\nSTEP E1: Correlation between clock axes (n = 18 categories)\n")
  ct_pearson  <- cor.test(clock_2d$med_slope, clock_2d$pct_persist, method = "pearson")
  ct_spearman <- suppressWarnings(
    cor.test(clock_2d$med_slope, clock_2d$pct_persist, method = "spearman")
  )

  cat(sprintf("Pearson  r = %.3f, 95%% CI [%.3f, %.3f], p = %.3f\n",
              ct_pearson$estimate, ct_pearson$conf.int[1], ct_pearson$conf.int[2],
              ct_pearson$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.3f (n=18 — CI omitted, not reliable at this n)\n",
              ct_spearman$estimate, ct_spearman$p.value))

  cor_summary <- data.table(
    method   = c("Pearson","Spearman"),
    estimate = c(ct_pearson$estimate, ct_spearman$estimate),
    ci_low   = c(ct_pearson$conf.int[1], NA),
    ci_high  = c(ct_pearson$conf.int[2], NA),
    p_value  = c(ct_pearson$p.value, ct_spearman$p.value)
  )
  fwrite(cor_summary, "ManuscriptFiles/Results/Biological_clock_axis_correlation.txt",
         sep = "\t", quote = FALSE)
  cat("Saved: Biological_clock_axis_correlation.txt\n")

  cat("\nSTEP E2: Silhouette width across k (k-means, scaled coordinates)\n")

  X <- clock_2d[, .(med_slope, pct_persist)]
  X_scaled <- scale(X)
  rownames(X_scaled) <- clock_2d$catA
  d_scaled <- dist(X_scaled)

  set.seed(42)
  sil_by_k <- rbindlist(lapply(2:8, function(k) {
    km <- kmeans(X_scaled, centers = k, nstart = 50, iter.max = 100)
    sil <- silhouette(km$cluster, d_scaled)
    data.table(k = k, mean_silhouette_width = mean(sil[, "sil_width"]))
  }))

  print(sil_by_k)
  cat(sprintf("\nBest k by silhouette width: k = %d (width = %.3f)\n",
              sil_by_k$k[which.max(sil_by_k$mean_silhouette_width)],
              max(sil_by_k$mean_silhouette_width)))
  cat("k=4 silhouette width:", sil_by_k[k == 4, mean_silhouette_width], "\n")
  cat("(Rule of thumb: <0.25 = no substantial cluster structure, ",
      "0.25-0.5 = weak/artificial, 0.5-0.7 = reasonable, >0.7 = strong.)\n")

  fwrite(sil_by_k, "ManuscriptFiles/Results/Biological_clock_silhouette_by_k.txt",
         sep = "\t", quote = FALSE)

  cat("\nSTEP E3: Gap statistic (k-means, B=500 bootstrap references)\n")

  set.seed(42)
  gap_km <- clusGap(X_scaled,
                     FUNcluster = function(x, k) kmeans(x, centers = k, nstart = 50),
                     K.max = 8, B = 500)

  print(gap_km, method = "Tibs2001SEmax")

  best_k_gap <- maxSE(gap_km$Tab[, "gap"], gap_km$Tab[, "SE.sim"], method = "Tibs2001SEmax")
  cat(sprintf("\nOptimal k by gap statistic (Tibshirani SE rule): k = %d\n", best_k_gap))
  if (best_k_gap == 1) {
    cat("NOTE: gap statistic selects k=1 — does NOT support discrete clusters;\n",
        "consistent with a continuous distribution rather than 4 separable classes.\n")
  } else if (best_k_gap == 4) {
    cat("NOTE: gap statistic supports k=4, consistent with the proposed quadrants.\n")
  } else {
    cat(sprintf("NOTE: gap statistic supports k=%d, NOT k=4 as originally proposed.\n", best_k_gap))
  }

  gap_table <- as.data.table(gap_km$Tab)
  gap_table[, k := 1:.N]
  fwrite(gap_table, "ManuscriptFiles/Results/Biological_clock_gap_statistic.txt",
         sep = "\t", quote = FALSE)
  cat("Saved: Biological_clock_gap_statistic.txt\n")

  cat("\nSTEP E4: k=4 stability under alternative specifications\n")

  adj_rand_index <- function(x, y) {
    tab <- table(x, y)
    a       <- sum(choose(tab, 2))
    bsum    <- sum(choose(rowSums(tab), 2))
    csum    <- sum(choose(colSums(tab), 2))
    n       <- sum(tab)
    expected  <- bsum * csum / choose(n, 2)
    maxindex  <- (bsum + csum) / 2
    if (maxindex == expected) return(1)
    (a - expected) / (maxindex - expected)
  }

  X_raw <- as.matrix(X)
  rownames(X_raw) <- clock_2d$catA

  set.seed(42)
  cl_original      <- clock_2d$quadrant 
  cl_kmeans_scaled <- kmeans(X_scaled, centers = 4, nstart = 50)$cluster
  cl_kmeans_raw    <- kmeans(X_raw,    centers = 4, nstart = 50)$cluster
  cl_hclust_ward   <- cutree(hclust(dist(X_scaled), method = "ward.D2"), k = 4)
  cl_pam_manhattan <- pam(X_scaled, k = 4, metric = "manhattan")$clustering

  cluster_specs <- list(
    "Median split (published)" = cl_original,
    "K-means, scaled, Euclidean" = cl_kmeans_scaled,
    "K-means, raw, Euclidean" = cl_kmeans_raw,
    "Hierarchical (Ward), scaled" = cl_hclust_ward,
    "PAM, scaled, Manhattan" = cl_pam_manhattan
  )

  spec_names <- names(cluster_specs)
  ari_matrix <- matrix(NA_real_, length(spec_names), length(spec_names),
                        dimnames = list(spec_names, spec_names))
  for (i in seq_along(spec_names)) {
    for (j in seq_along(spec_names)) {
      ari_matrix[i, j] <- adj_rand_index(cluster_specs[[i]], cluster_specs[[j]])
    }
  }

  cat("\nAdjusted Rand Index between clustering specifications (1 = perfect agreement, 0 = chance):\n")
  print(round(ari_matrix, 2))

  cat("\nAgreement with published median-split assignment specifically:\n")
  ari_vs_original <- data.table(
    specification = spec_names,
    ARI_vs_published = round(ari_matrix[, "Median split (published)"], 2)
  )
  print(ari_vs_original)

  fwrite(as.data.table(ari_matrix, keep.rownames = "specification"),
         "ManuscriptFiles/Results/Biological_clock_ARI_matrix.txt",
         sep = "\t", quote = FALSE)
  cat("Saved: Biological_clock_ARI_matrix.txt\n")

  cat("\n=== Cluster validity analysis complete ===\n")
}