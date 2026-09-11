## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

if (!dir.exists("ManuscriptFiles/Plots/Subcategory"))
  dir.create("ManuscriptFiles/Plots/Subcategory", recursive = TRUE)
if (!dir.exists("ManuscriptFiles/Results/Subcategory"))
  dir.create("ManuscriptFiles/Results/Subcategory", recursive = TRUE)

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
windows_cont  <- c("1_2","2_3","3_4","4_5")

colcod <- c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
            "#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54",
            "#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
names(colcod) <- unique(cate)

subcat_map <- fread("./icd10_3digitos_categoria_subcategoria_who2016.csv")
subcat_map <- subcat_map[, .(
  icd10_3             = as.character(icd10_3),
  subcategoria_rango  = as.character(subcategoria_rango),
  subcategoria_nombre = as.character(subcategoria_nombre),
  categoria_nombre    = as.character(categoria_nombre)
)]
subcat_rango  <- setNames(subcat_map$subcategoria_rango,  subcat_map$icd10_3)
subcat_nombre <- setNames(subcat_map$subcategoria_nombre, subcat_map$icd10_3)

all_cumul <- rbindlist(lapply(windows_cumul, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f); dt <- dt[eval(SIG_FILTER)]
  dt[, `:=`(window = w, window_num = as.integer(gsub("0_","",w)),
            pair = paste(disease_a, disease_b, sep = "_"))]
  dt[, subcatA := subcat_rango[disease_a]]
  dt[, subcatA_nombre := subcat_nombre[disease_a]]
  dt[, catA := catname[cate[disease_a]]]
  dt
}))
persistent_pairs <- all_cumul[, .(n_windows = uniqueN(window)), by = pair][n_windows == 5, pair]
cat(sprintf("  Persistent pairs: %d\n", length(persistent_pairs)))

rr_traj <- all_cumul[pair %in% persistent_pairs,
                     .(pair, window_num, RR_shrunk, SE_RR_shrunk, subcatA, subcatA_nombre, catA)]
slopes <- rr_traj[, {
  x <- window_num; y <- RR_shrunk; se <- SE_RR_shrunk
  if (any(is.na(se)) || any(se <= 0) || length(x) < 3) .(slope = NA_real_)
  else { fit <- lm(y ~ x, weights = 1/se^2); .(slope = coef(fit)["x"]) }
}, by = .(pair, subcatA, subcatA_nombre, catA)]

traj_by_subcat <- slopes[!is.na(subcatA) & !is.na(slope), .(
  n_total = .N, med_slope = median(slope), catA = first(catA), subcatA_nombre = first(subcatA_nombre)
), by = subcatA][order(med_slope)]

MIN_PAIRS_SUBCAT <- 10
traj_by_subcat_filt <- traj_by_subcat[n_total >= MIN_PAIRS_SUBCAT]
cat(sprintf("  Subcategories after filtering (>=%d pairs): %d\n", MIN_PAIRS_SUBCAT, nrow(traj_by_subcat_filt)))

all_cont <- rbindlist(lapply(windows_cont, function(w) {
  f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f); dt <- dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
  dt[, `:=`(window = w, pair = paste(disease_a, disease_b, sep = "_"),
            sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  dt[, subcatA := subcat_rango[disease_a]]
  dt[, subcatA_nombre := subcat_nombre[disease_a]]
  dt[, catA := catname[cate[disease_a]]]
  dt
}))
pairs_cont4 <- all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]

all_cont_filt <- all_cont[pair %in% pairs_cont4 & !is.na(subcatA),
                          .(pair, window, sig, subcatA, subcatA_nombre, catA)]

cont_wide <- dcast(all_cont_filt, pair + subcatA + subcatA_nombre + catA ~ window, value.var = "sig")
cont_wide[, pattern := fcase(
  (`1_2`==TRUE|`2_3`==TRUE)&(`3_4`==FALSE&`4_5`==FALSE), "Early only",
  (`1_2`==TRUE|`2_3`==TRUE)&(`3_4`==TRUE|`4_5`==TRUE),   "Persistent conditional risk",
  (`1_2`==FALSE&`2_3`==FALSE)&(`3_4`==TRUE|`4_5`==TRUE), "Late emerging risk",
  rep(TRUE,.N), "Not significant"
)]

pat_by_subcat <- cont_wide[!is.na(subcatA) & !is.na(pattern), .(
  pct_persist = mean(pattern == "Persistent conditional risk") * 100
), by = subcatA]

clock_2d_subcat <- merge(
  traj_by_subcat_filt[, .(subcatA, subcatA_nombre, catA, med_slope, n_total)],
  pat_by_subcat[, .(subcatA, pct_persist)],
  by = "subcatA"
)

x_mid_sub <- median(clock_2d_subcat$med_slope)
y_mid_sub <- median(clock_2d_subcat$pct_persist)
cat(sprintf("  x_mid_sub = %.5f, y_mid_sub = %.3f\n", x_mid_sub, y_mid_sub))
cat(sprintf("  Subcategories in final figure: %d\n", nrow(clock_2d_subcat)))

cat_colors_sub <- setNames(
  as.character(colcod[names(colcod) %in% names(catname)[catname %in% unique(clock_2d_subcat$catA)]]),
  catname[catname %in% unique(clock_2d_subcat$catA)]
)
unmatched <- setdiff(unique(clock_2d_subcat$catA), names(cat_colors_sub))
if (length(unmatched) > 0) {
  cat_colors_sub <- c(cat_colors_sub, setNames(rep("grey50", length(unmatched)), unmatched))
}

dat_a <- clock_2d_subcat[catA != "Pregnancy"]
dat_b <- clock_2d_subcat[catA == "Pregnancy"]

build_panel <- function(dat, x_mid, y_mid, show_legend, title_txt = NULL) {
  ggplot(dat, aes(x = med_slope, y = pct_persist, color = catA, size = n_total)) +
    geom_vline(xintercept = x_mid, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = y_mid, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(alpha = 0.8) +
    geom_text_repel(aes(label = subcatA), size = 2.3, color = "grey20",
                    max.overlaps = 40, box.padding = 0.3, segment.size = 0.3,
                    min.segment.length = 0.2, show.legend = FALSE) +
    scale_color_manual(values = cat_colors_sub, name = "ICD-10 category") +
    scale_size_continuous(range = c(2,10), name = "N persistent pairs") +
    labs(x = "Rate of RR attenuation (median slope, cumulative windows)",
        y = "% comorbidities with persistent conditional risk",
        title = title_txt) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11)
    )
}

p_a <- build_panel(dat_a, x_mid_sub, y_mid_sub, show_legend = TRUE,
                   "a) All ICD-10 categories except Pregnancy")
x_mid_preg <- median(dat_b$med_slope)
y_mid_preg <- median(dat_b$pct_persist)
p_b <- build_panel(dat_b, x_mid_preg, y_mid_preg, show_legend = FALSE,
                   "b) Pregnancy subcategories (note different x-axis scale)")

p2d_sub <- p_a + p_b + patchwork::plot_layout(widths = c(2.2, 1))

ggsave("ManuscriptFiles/Plots/Subcategory/Clock_2D_subcategory.pdf",
       p2d_sub, width = 18, height = 9, useDingbats = FALSE)
cat("\nSaved: Clock_2D_subcategory.pdf\n")

fwrite(clock_2d_subcat[order(catA, med_slope)], "ManuscriptFiles/Results/Subcategory/Clock_2D_subcategory.txt", sep = "\t", quote = FALSE)