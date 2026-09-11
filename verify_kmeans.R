## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({ library(data.table) })

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

windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")
all_cumul <- rbindlist(lapply(windows_cumul, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  dt <- fread(f); dt <- dt[eval(SIG_FILTER)]
  dt[, `:=`(window_num = as.integer(gsub("0_","",w)), pair = paste(disease_a, disease_b, sep = "_"))]
  dt
}))
persistent_pairs <- all_cumul[, .(n_windows = uniqueN(window_num)), by = pair][n_windows == 5, pair]
rr_traj <- all_cumul[pair %in% persistent_pairs, .(pair, window_num, RR_shrunk, SE_RR_shrunk, disease_a)]
rr_traj[, catA := catname[cate[disease_a]]]
slopes <- rr_traj[, {
  x <- window_num; y <- RR_shrunk; se <- SE_RR_shrunk
  if (any(is.na(se)) || any(se <= 0) || length(x) < 3) .(slope = NA_real_)
  else { fit <- lm(y ~ x, weights = 1/se^2); .(slope = coef(fit)["x"]) }
}, by = .(pair, catA)]
traj_by_cat <- slopes[!is.na(catA) & !is.na(slope), .(n_total = .N, med_slope = median(slope)), by = catA]

windows_cont <- c("1_2","2_3","3_4","4_5")
all_cont <- rbindlist(lapply(windows_cont, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  dt <- fread(f); dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
  dt[, `:=`(window = w, pair = paste(disease_a, disease_b, sep = "_"),
            sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  dt
}))
pairs_cont4 <- all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]
cont_wide <- dcast(all_cont[pair %in% pairs_cont4, .(pair, window, sig)], pair ~ window, value.var = "sig")
cont_wide[, catA := catname[cate[gsub("_.+","",pair)]]]
cont_wide[, pattern := fcase(
  (`1_2`==TRUE|`2_3`==TRUE)&(`3_4`==FALSE&`4_5`==FALSE), "Early only",
  (`1_2`==TRUE|`2_3`==TRUE)&(`3_4`==TRUE|`4_5`==TRUE),   "Persistent conditional risk",
  (`1_2`==FALSE&`2_3`==FALSE)&(`3_4`==TRUE|`4_5`==TRUE), "Late emerging risk",
  rep(TRUE,.N), "Not significant"
)]
pat_by_cat <- cont_wide[!is.na(catA)&!is.na(pattern), .(pct_persist = mean(pattern=="Persistent conditional risk")*100), by=catA]

clock_2d <- merge(traj_by_cat, pat_by_cat, by = "catA")
x_mid <- median(clock_2d$med_slope); y_mid <- median(clock_2d$pct_persist)
clock_2d[, quadrant := fcase(
  med_slope <= x_mid & pct_persist >= y_mid, "Late-emerging",
  med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive",
  med_slope <= x_mid & pct_persist <  y_mid, "Episodic",
  rep(TRUE,.N), "Chronic stable"
)]

cat("=== clock_2d (verify identical to previous scripts before trusting below) ===\n")
print(clock_2d[order(catA), .(catA, med_slope=round(med_slope,4), pct_persist=round(pct_persist,1), quadrant)])

## ---- K-means, k=4, exact same specification as originally reported ----
clock_norm <- copy(clock_2d)
clock_norm[, slope_z := scale(med_slope)[,1]]
clock_norm[, persist_z := scale(pct_persist)[,1]]

set.seed(42)
km <- kmeans(clock_norm[, .(slope_z, persist_z)], centers = 4, nstart = 50, iter.max = 100)
clock_norm[, cluster_kmeans := km$cluster]

centers <- as.data.table(km$centers)
centers[, cluster := 1:.N]
centers[, quadrant_km := fcase(
  slope_z <= 0 & persist_z >= 0, "Late-emerging",
  slope_z >  0 & persist_z >= 0, "Chronic progressive",
  slope_z <= 0 & persist_z <  0, "Episodic",
  rep(TRUE, .N), "Chronic stable"
)]
clock_norm <- merge(clock_norm, centers[, .(cluster, quadrant_km)], by.x = "cluster_kmeans", by.y = "cluster")

concordance <- merge(clock_norm[, .(catA, quadrant_km)], clock_2d[, .(catA, quadrant_orig = quadrant)], by = "catA")
concordance[, same := quadrant_km == quadrant_orig]

cat("\n=== K-MEANS VERIFICATION (settles the Circulatory question) ===\n")
print(concordance[order(catA)])
cat(sprintf("\nOverall concordance: %d/%d categories (%.1f%%)\n",
            sum(concordance$same), nrow(concordance), mean(concordance$same)*100))

cat("\n=== Circulatory specifically ===\n")
print(concordance[catA == "Circulatory"])