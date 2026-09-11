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

build_taxonomy <- function(sig_filter_expr, label) {

  cat(sprintf("\n--- Building taxonomy: %s ---\n", label))

  windows_9   <- c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  window_type <- c(rep("cumulative",5), rep("conditional",4))

  all_9 <- rbindlist(lapply(seq_along(windows_9), function(i) {
    w  <- windows_9[i]
    wt <- window_type[i]
    f  <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt <- fread(f)
    dt <- dt[eval(sig_filter_expr)]
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

  cat(sprintf("  Total significant in >=1 of 9 windows: %d\n", nrow(pair_windows)))

  net_size <- sapply(c("0_1","0_2","0_3","0_4","0_5"), function(w) {
    pair_windows[get(paste0("in_", w)) == TRUE, .N]
  })

  list(pair_windows = pair_windows,
       class_counts = pair_windows[, .N, by = bio_class][order(-N)],
       net_size = net_size)
}

SIG_ORIGINAL <- quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)
SIG_RELAXED  <- quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01)   ## no n>=100 floor

res_original <- build_taxonomy(SIG_ORIGINAL, "ORIGINAL (n>=100)")
res_relaxed  <- build_taxonomy(SIG_RELAXED,  "RELAXED (no n>=100)")

cat("\n=== Sanity check: original criterion vs manuscript ===\n")
cat("Manuscript: Total=146,591 | Omnipresent=45,495 | Long-window only=72,495 |\n")
cat("            Early transient=157 | Late cumul+cond=15,645\n")
print(res_original$class_counts)

cat("\n=== STEP: Side-by-side comparison, all 9 categories ===\n")

comp <- merge(res_original$class_counts, res_relaxed$class_counts,
              by = "bio_class", all = TRUE, suffixes = c("_original", "_relaxed"))
comp[is.na(N_original), N_original := 0]
comp[is.na(N_relaxed),  N_relaxed  := 0]
comp[, pct_original := round(N_original / sum(N_original) * 100, 2)]
comp[, pct_relaxed  := round(N_relaxed  / sum(N_relaxed)  * 100, 2)]
comp[, fold_change  := round(N_relaxed / pmax(N_original, 1), 2)]
setorder(comp, -N_relaxed)
print(comp)

cat("\n=== Specifically: Early transient (the category this question is about) ===\n")
print(comp[bio_class == "Early transient (short window only)"])

cat("\n=== STEP: Network size (significant pairs) at each cumulative window ===\n")

net_comp <- data.table(
  window = c("0-1","0-2","0-3","0-4","0-5"),
  n_original = res_original$net_size,
  n_relaxed  = res_relaxed$net_size
)
net_comp[, fold_increase_original := round(n_original / n_original[1], 2)]
net_comp[, fold_increase_relaxed  := round(n_relaxed  / n_relaxed[1], 2)]
print(net_comp)

for (res in list(list(name="ORIGINAL", d=res_original), list(name="RELAXED", d=res_relaxed))) {
  pw <- res$d$pair_windows
  late <- pw[, .(late_emerging = bio_class %in% c(
    "Long-window only (late cumulative, no early signal)",
    "Late cumulative + conditional risk"))]
  n_late <- sum(late$late_emerging)
  cat(sprintf("\n[%s] Late-emerging: %d (%.1f%% of %d)\n",
              res$name, n_late, 100*n_late/nrow(pw), nrow(pw)))
}

fwrite(comp, "ManuscriptFiles/Results/Taxonomy_original_vs_relaxed_comparison.txt", sep = "\t", quote = FALSE)
fwrite(net_comp, "ManuscriptFiles/Results/Taxonomy_network_size_original_vs_relaxed.txt", sep = "\t", quote = FALSE)