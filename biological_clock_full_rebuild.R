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

cat("STEP A: Rebuilding `slopes` from cumulative windows\n")

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
cat(sprintf("  Persistent pairs (all 5 cumulative windows): %d\n", length(persistent_pairs)))

rr_traj <- all_cumul[pair %in% persistent_pairs,
                     .(pair, window_num, RR_shrunk, SE_RR_shrunk, disease_a, disease_b)]
rr_traj[, catA := catname[cate[disease_a]]]

slopes <- rr_traj[, {
  x  <- window_num; y <- RR_shrunk; se <- SE_RR_shrunk
  if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
    .(slope = NA_real_)
  } else {
    fit <- lm(y ~ x, weights = 1/se^2)
    .(slope = coef(fit)["x"])
  }
}, by = .(pair, catA)]

traj_by_cat <- slopes[!is.na(catA) & !is.na(slope),
                      .(n_total = .N, med_slope = median(slope)), by = catA]

cat("\nSTEP B: Rebuilding `cont_wide` from conditional windows\n")
windows_cont <- c("1_2","2_3","3_4","4_5")

all_cont <- rbindlist(lapply(windows_cont, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
  dt[, `:=`(window = w,
            pair = paste(disease_a, disease_b, sep = "_"),
            sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  dt
}))

pairs_cont4 <- all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]
cat(sprintf("  Pairs present in all 4 conditional windows: %d\n", length(pairs_cont4)))

cont_wide <- dcast(
  all_cont[pair %in% pairs_cont4, .(pair, window, sig)],
  pair ~ window, value.var = "sig"
)
cont_wide[, catA := catname[cate[gsub("_.+","",pair)]]]

cont_wide[, pattern := fcase(
  (`1_2` == TRUE | `2_3` == TRUE) & (`3_4` == FALSE & `4_5` == FALSE),
  "Early only (risk fades)",
  (`1_2` == TRUE | `2_3` == TRUE) & (`3_4` == TRUE  | `4_5` == TRUE),
  "Persistent conditional risk",
  (`1_2` == FALSE & `2_3` == FALSE) & (`3_4` == TRUE | `4_5` == TRUE),
  "Late emerging risk",
  rep(TRUE, .N),
  "Not significant"
)]

pat_by_cat <- cont_wide[!is.na(catA) & !is.na(pattern),
                        .(pct_persist = mean(pattern == "Persistent conditional risk") * 100,
                          pct_late    = mean(pattern == "Late emerging risk") * 100),
                        by = catA]
cat("  pct_persist range (should be single digits to ~20%, NOT 55-95%):\n")
print(summary(pat_by_cat$pct_persist))

cat("\nSTEP C: Building clock_2d\n")

clock_2d <- merge(traj_by_cat, pat_by_cat[, .(catA, pct_persist, pct_late)], by = "catA")

x_mid <- median(clock_2d$med_slope)
y_mid <- median(clock_2d$pct_persist)
cat(sprintf("  x_mid = %.5f, y_mid = %.3f\n", x_mid, y_mid))

clock_2d[, quadrant := fcase(
  med_slope <= x_mid & pct_persist >= y_mid, "Late-emerging",
  med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive",
  med_slope <= x_mid & pct_persist <  y_mid, "Episodic",
  rep(TRUE, .N),                              "Chronic stable"
)]

clock_2d[, catA_short := fcase(
  catA == "Endocrine/Metabolic", "Endocrine/\nMetabolic",
  catA == "Musculoskeletal",     "Musculo-\nskeletal",
  catA == "Blood/Immune",        "Blood/\nImmune",
  catA == "Genitourinary",       "Genito-\nurinary",
  rep(TRUE, .N),                  catA
)]

print(clock_2d[order(quadrant, -med_slope), .(catA, med_slope, pct_persist, quadrant)])

fwrite(clock_2d, "ManuscriptFiles/Results/Biological_clock_2D_verified.txt", sep = "\t", quote = FALSE)

cat("STEP D: Pair-level bootstrap (1,000 replicates, axes independent)\n")

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
      med_slope = if (length(s) >= 1) median(sample(s, length(s), replace = TRUE), na.rm = TRUE) else NA_real_,
      pct_persist = if (length(p) >= 1) mean(sample(p, length(p), replace = TRUE) ==
                                                "Persistent conditional risk", na.rm = TRUE) * 100 else NA_real_
    )
  }))
}

boot_results <- rbindlist(lapply(seq_len(N_BOOT), function(i) {
  res <- boot_one()
  res[, replicate := i]
  res
}))
cat(sprintf("  Bootstrap complete: %d replicates x %d categories\n", N_BOOT, uniqueN(boot_results$catA)))

boot_summary <- boot_results[, {
  quad <- fcase(
    med_slope <= x_mid & pct_persist >= y_mid, "Late-emerging",
    med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive",
    med_slope <= x_mid & pct_persist <  y_mid, "Episodic",
    rep(TRUE, .N),                              "Chronic stable"
  )
  orig_quad <- clock_2d[catA == .BY$catA, quadrant]
  if (length(orig_quad) == 0) orig_quad <- NA_character_
  .(slope_lo = quantile(med_slope, 0.025, na.rm=TRUE), slope_hi = quantile(med_slope, 0.975, na.rm=TRUE),
    persist_lo = quantile(pct_persist, 0.025, na.rm=TRUE), persist_hi = quantile(pct_persist, 0.975, na.rm=TRUE),
    persist_mean = mean(pct_persist, na.rm=TRUE),
    quadrant_retention = mean(quad == orig_quad, na.rm = TRUE) * 100,
    original_quadrant = orig_quad[1])
}, by = catA]

cat("\n=== BOOTSTRAP SUMMARY (verified-correct methodology) ===\n")
print(boot_summary[order(-quadrant_retention),
      .(catA, original_quadrant, quadrant_retention = round(quadrant_retention, 1))])

cat("\nSanity check -- persist_mean should resemble clock_2d's pct_persist values above\n")
cat("(single digits to ~20%), NOT 55-95%:\n")
print(boot_summary[, .(catA, persist_mean = round(persist_mean,1))][order(-persist_mean)])

fwrite(boot_summary, "ManuscriptFiles/Results/Biological_clock_bootstrap_verified.txt", sep = "\t", quote = FALSE)
cat("\nSaved: Biological_clock_bootstrap_verified.txt\n")

cat("\nSTEP E: Building Figure 2a\n")

p2d <- ggplot(clock_2d, aes(x = med_slope, y = pct_persist,
                            color = catA, size = n_total)) +
  geom_vline(xintercept = x_mid, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = y_mid, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  stat_ellipse(data = boot_results,
               aes(x = med_slope, y = pct_persist, color = catA),
               level = 0.95, linewidth = 0.5, alpha = 0.6, inherit.aes = FALSE) +
  geom_point(alpha = 0.9) +
  geom_text_repel(aes(label = catA_short), size = 6, color = "grey20",
                  max.overlaps = 20, box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = catcol, guide = "none") +
  scale_size_continuous(range = c(4, 12), name = "N comorbidities") +
  annotate("text", x = x_mid - 0.003, y = Inf, label = "Faster attenuation",
           hjust = 1, vjust = 1.5, color = "grey40", size = 5, fontface = "italic") +
  annotate("text", x = x_mid + 0.001, y = Inf, label = "Slower attenuation",
           hjust = 0, vjust = 1.5, color = "grey40", size = 5, fontface = "italic") +
  labs(
    x = "Rate of RR attenuation (median slope, cumulative windows)",
    y = "% of comorbidities with persistent conditional risk",
    title = "The two-dimensional biological clock of multimorbidity",
    subtitle = paste0(
      "X-axis computed over persistent pairs (all 5 cumulative windows); ",
      "Y-axis computed over the full conditional-window population.\n",
      "Each point = one ICD-10 category. Size = number of comorbidity pairs. ",
      "Ellipses = 95% bootstrap confidence regions (1,000 replicates, axes resampled independently)."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey40", size = 9))

ggsave("ManuscriptFiles/Plots/Biological_clock_2D_continuous_verified.pdf",
       p2d, width = 13, height = 9, useDingbats = FALSE)
cat("Saved: Biological_clock_2D_continuous_verified.pdf\n")