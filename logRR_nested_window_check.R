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

cat("STEP 1: Loading cumulative windows\n")

windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")
all_cumul <- rbindlist(lapply(windows_cumul, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  dt <- fread(f)
  dt <- dt[eval(SIG_FILTER)]
  dt[, `:=`(window_num = as.integer(gsub("0_","",w)), pair = paste(disease_a, disease_b, sep = "_"))]
  dt
}))

persistent_pairs <- all_cumul[, .(n_windows = uniqueN(window_num)), by = pair][n_windows == 5, pair]
cat(sprintf("  Persistent pairs: %d\n", length(persistent_pairs)))

rr_traj <- all_cumul[pair %in% persistent_pairs,
                     .(pair, window_num, RR_shrunk, SE_RR_shrunk, log_RR_shrunk, disease_a)]
rr_traj[, catA := catname[cate[disease_a]]]
## log-scale SE via delta method: SE(log RR) ~ SE(RR) / RR
rr_traj[, SE_logRR := SE_RR_shrunk / RR_shrunk]

cat("\nSTEP 2: Original slope (raw RR, 5-point OLS)\n")

slopes_raw <- rr_traj[, {
  x <- window_num; y <- RR_shrunk; se <- SE_RR_shrunk
  if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
    .(slope = NA_real_)
  } else {
    fit <- lm(y ~ x, weights = 1/se^2)
    .(slope = coef(fit)["x"])
  }
}, by = .(pair, catA)]

cat("STEP 3: Log-RR slope (5-point OLS on log scale)\n")

slopes_log <- rr_traj[, {
  x <- window_num; y <- log_RR_shrunk; se <- SE_logRR
  if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
    .(slope_log = NA_real_)
  } else {
    fit <- lm(y ~ x, weights = 1/se^2)
    .(slope_log = coef(fit)["x"])
  }
}, by = .(pair, catA)]

cat("STEP 4: Two-point log-attenuation (0-1 vs 0-5, nesting-free)\n")

wide_traj <- dcast(rr_traj, pair + catA ~ window_num,
                   value.var = c("log_RR_shrunk", "SE_logRR"))
setnames(wide_traj,
        c("log_RR_shrunk_1","log_RR_shrunk_5","SE_logRR_1","SE_logRR_5"),
        c("logRR_w1","logRR_w5","SE_w1","SE_w5"))

wide_traj[, twopoint_atten := logRR_w5 - logRR_w1]

cat("\nSTEP 5: Category-level attenuation, three methods\n")

cat_raw <- slopes_raw[!is.na(catA) & !is.na(slope), .(med_slope_raw = median(slope)), by = catA]
cat_log <- slopes_log[!is.na(catA) & !is.na(slope_log), .(med_slope_log = median(slope_log)), by = catA]
cat_2pt <- wide_traj[!is.na(catA) & !is.na(twopoint_atten), .(med_2pt = median(twopoint_atten)), by = catA]

comp <- merge(cat_raw, cat_log, by = "catA")
comp <- merge(comp, cat_2pt, by = "catA")

comp[, rank_raw := rank(med_slope_raw)]
comp[, rank_log := rank(med_slope_log)]
comp[, rank_2pt := rank(med_2pt)]

setorder(comp, rank_log)
print(comp[, .(catA, med_slope_raw = round(med_slope_raw,4),
               med_slope_log = round(med_slope_log,4),
               med_2pt = round(med_2pt,4),
               rank_raw, rank_log, rank_2pt)])

cat("\n=== Spearman correlation between category-level rankings ===\n")
cat(sprintf("  Raw-RR slope vs. log-RR slope:      rho = %.3f\n",
            cor(comp$med_slope_raw, comp$med_slope_log, method = "spearman")))
cat(sprintf("  Log-RR 5-point slope vs. two-point: rho = %.3f\n",
            cor(comp$med_slope_log, comp$med_2pt, method = "spearman")))
cat(sprintf("  Raw-RR slope vs. two-point measure: rho = %.3f\n",
            cor(comp$med_slope_raw, comp$med_2pt, method = "spearman")))

cat("\nSTEP 6: Quadrant assignment -- original vs. log-RR two-point\n")

windows_cont <- c("1_2","2_3","3_4","4_5")
all_cont <- rbindlist(lapply(windows_cont, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  dt <- fread(f)
  dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
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
pat_by_cat <- cont_wide[!is.na(catA)&!is.na(pattern),
                        .(pct_persist = mean(pattern=="Persistent conditional risk")*100), by=catA]

clock_orig <- merge(comp[, .(catA, med_slope_raw)], pat_by_cat, by = "catA")
x_mid_orig <- median(clock_orig$med_slope_raw)
y_mid <- median(clock_orig$pct_persist)
clock_orig[, quadrant_orig := fcase(
  med_slope_raw <= x_mid_orig & pct_persist >= y_mid, "Late-emerging",
  med_slope_raw >  x_mid_orig & pct_persist >= y_mid, "Chronic progressive",
  med_slope_raw <= x_mid_orig & pct_persist <  y_mid, "Episodic",
  rep(TRUE,.N), "Chronic stable"
)]

clock_2pt <- merge(comp[, .(catA, med_2pt)], pat_by_cat, by = "catA")
x_mid_2pt <- median(clock_2pt$med_2pt)
clock_2pt[, quadrant_2pt := fcase(
  med_2pt <= x_mid_2pt & pct_persist >= y_mid, "Late-emerging",
  med_2pt >  x_mid_2pt & pct_persist >= y_mid, "Chronic progressive",
  med_2pt <= x_mid_2pt & pct_persist <  y_mid, "Episodic",
  rep(TRUE,.N), "Chronic stable"
)]

quad_compare <- merge(clock_orig[, .(catA, quadrant_orig)], clock_2pt[, .(catA, quadrant_2pt)], by = "catA")
quad_compare[, same := quadrant_orig == quadrant_2pt]
print(quad_compare[order(catA)])
cat(sprintf("\nCategories with unchanged quadrant assignment: %d / %d (%.1f%%)\n",
            sum(quad_compare$same), nrow(quad_compare), 100*mean(quad_compare$same)))

fwrite(comp, "ManuscriptFiles/Results/LogRR_nested_window_comparison.txt", sep = "\t", quote = FALSE)
fwrite(quad_compare, "ManuscriptFiles/Results/LogRR_nested_window_quadrants.txt", sep = "\t", quote = FALSE)