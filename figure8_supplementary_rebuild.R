## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

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

colcod <- c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
            "#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54",
            "#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
names(colcod) <- unique(cate)
catcol <- as.character(colcod)
names(catcol) <- catname

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
x_mid <- median(clock_2d$med_slope)
y_mid <- median(clock_2d$pct_persist)
clock_2d[, catA_short := fcase(
  catA == "Endocrine/Metabolic", "Endocrine/\nMetabolic",
  catA == "Musculoskeletal",     "Musculo-\nskeletal",
  catA == "Blood/Immune",        "Blood/\nImmune",
  catA == "Genitourinary",       "Genito-\nurinary",
  rep(TRUE, .N),                  catA
)]
cat(sprintf("clock_2d built. x_mid=%.5f, y_mid=%.3f (verify against Table 2 before trusting below)\n", x_mid, y_mid))
print(clock_2d[order(catA), .(catA, med_slope=round(med_slope,4), pct_persist=round(pct_persist,1))])

set.seed(42)
N_BOOT <- 1000
slopes_src  <- slopes[!is.na(slope) & !is.na(catA), .(catA, slope)]
persist_src <- cont_wide[!is.na(catA), .(catA, pattern)]
cats <- sort(unique(clock_2d$catA))

boot_one <- function() {
  rbindlist(lapply(cats, function(cat_i) {
    s <- slopes_src[catA == cat_i, slope]
    p <- persist_src[catA == cat_i, pattern]
    data.table(
      catA = cat_i,
      med_slope = if (length(s)>=1) median(sample(s, length(s), replace=TRUE), na.rm=TRUE) else NA_real_,
      pct_persist = if (length(p)>=1) mean(sample(p, length(p), replace=TRUE)=="Persistent conditional risk", na.rm=TRUE)*100 else NA_real_
    )
  }))
}
boot_results <- rbindlist(lapply(seq_len(N_BOOT), function(i) { res <- boot_one(); res[, replicate := i]; res }))
cat(sprintf("Bootstrap complete: %d replicates x %d categories\n", N_BOOT, uniqueN(boot_results$catA)))

clock_2d[, quadrant_orig := fcase(
  med_slope<=x_mid & pct_persist>=y_mid, "Late-emerging",
  med_slope>x_mid & pct_persist>=y_mid, "Chronic progressive",
  med_slope<=x_mid & pct_persist<y_mid, "Episodic",
  rep(TRUE,.N), "Chronic stable"
)]
boot_results[, quadrant_boot := fcase(
  med_slope<=x_mid & pct_persist>=y_mid, "Late-emerging",
  med_slope>x_mid & pct_persist>=y_mid, "Chronic progressive",
  med_slope<=x_mid & pct_persist<y_mid, "Episodic",
  rep(TRUE,.N), "Chronic stable"
)]
boot_summary <- merge(boot_results, clock_2d[, .(catA, quadrant_orig)], by = "catA")
boot_summary <- boot_summary[, .(retention = mean(quadrant_boot == quadrant_orig)*100), by = catA]
cat("\nQuadrant retention (sanity check -- compare against Supplementary Table 2):\n")
print(boot_summary[order(-retention)])

clock_2d <- merge(clock_2d, boot_summary, by = "catA")
clock_2d[, label_full := sprintf("%s\n(%.0f%%)", catA_short, retention)]

boot_plot <- merge(boot_results, clock_2d[, .(catA, catA_short)], by = "catA")

p8 <- ggplot(clock_2d, aes(x = med_slope, y = pct_persist, color = catA)) +
  geom_vline(xintercept = x_mid, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = y_mid, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(data = boot_plot[replicate <= 200],
             aes(x = med_slope, y = pct_persist, color = catA),
             alpha = 0.06, size = 0.8, inherit.aes = FALSE) +
  stat_ellipse(data = boot_results, aes(x = med_slope, y = pct_persist, color = catA),
               level = 0.95, linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(size = 3.5, color = "black") +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label_full), size = 3.6, color = "grey20",
                  max.overlaps = 20, box.padding = 0.35, show.legend = FALSE) +
  scale_color_manual(values = catcol, guide = "none") +
  labs(
    x = "Rate of RR attenuation (median slope, cumulative windows)",
    y = "% of comorbidities with persistent conditional risk"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("ManuscriptFiles/Plots/Supplementary_Figure8_bootstrap_stability.pdf",
       p8, width = 12, height = 9, useDingbats = FALSE)
cat("\nSaved: Supplementary_Figure8_bootstrap_stability.pdf\n")

fwrite(boot_summary, "ManuscriptFiles/Results/Figure8_bootstrap_retention_check.txt", sep = "\t", quote = FALSE)