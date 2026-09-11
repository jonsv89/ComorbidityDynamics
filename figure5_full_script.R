## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

## ------------------------------------------------------------------ ##
## ICD-10 category mapping (identical to rest of pipeline)             ##
## ------------------------------------------------------------------ ##
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

abbrev <- c(
  "Infectious"          = "Infect.",
  "Neoplasms"           = "Neopl.",
  "Blood/Immune"        = "Blood/Imm.",
  "Endocrine/Metabolic" = "Endocr.",
  "Mental"              = "Mental",
  "Nervous"             = "Nervous",
  "Eye"                 = "Eye",
  "Ear"                 = "Ear",
  "Circulatory"         = "Circ.",
  "Respiratory"         = "Resp.",
  "Digestive"           = "Digest.",
  "Skin"                = "Skin",
  "Musculoskeletal"     = "Muscul.",
  "Genitourinary"       = "Genito.",
  "Pregnancy"           = "Preg.",
  "Perinatal"           = "Perinat.",
  "Congenital"          = "Congenit.",
  "Symptoms"            = "Sympt.",
  "Injury"              = "Injury"
)

net_w_05 <- fread("Results/networks/RR_net_women_0_5.txt")[
  lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
net_m_05 <- fread("Results/networks/RR_net_men_0_5.txt")[
  lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]

net_w_05[, pair := paste(disease_a, disease_b, sep = "_")]
net_m_05[, pair := paste(disease_a, disease_b, sep = "_")]
net_w_05[, catA := catname[cate[disease_a]]]
net_w_05[, catB := catname[cate[disease_b]]]
net_m_05[, catA := catname[cate[disease_a]]]
net_m_05[, catB := catname[cate[disease_b]]]

pairs_w <- net_w_05$pair
pairs_m <- net_m_05$pair
shared  <- intersect(pairs_w, pairs_m)

net_w_05[, status := ifelse(pair %in% shared, "shared", "only_women")]
net_m_05[, status := ifelse(pair %in% shared, "shared", "only_men")]

net_sex <- rbindlist(list(
  net_w_05[, .(pair, catA, catB, status, sex = "women")],
  net_m_05[!pair %in% shared, .(pair, catA, catB, status, sex = "men")]
))

cat(sprintf("  Shared: %d | Only women: %d | Only men: %d\n",
            length(shared), sum(net_sex$status == "only_women"), sum(net_sex$status == "only_men")))

universe_AB <- net_sex[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
shared_AB   <- net_sex[status == "shared" & !is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
setnames(shared_AB, "N", "n_shared")

enrich_shared <- merge(universe_AB, shared_AB, by = c("catA","catB"), all.x = TRUE)
enrich_shared[is.na(n_shared), n_shared := 0L]

total_shared   <- length(shared)
total_universe <- nrow(net_sex[!is.na(catA) & !is.na(catB)])

enrich_shared <- enrich_shared[, {
  a1 <- n_shared; b1 <- total_shared - a1
  a0 <- N - a1;    b0 <- (total_universe - total_shared) - a0
  if (any(c(a1,b1,a0,b0) < 0)) {
    .(odds = NA_real_, p = NA_real_)
  } else {
    ft <- fisher.test(matrix(c(a1,b1,a0,b0), nrow = 2, byrow = TRUE))
    .(odds = as.numeric(ft$estimate), p = ft$p.value)
  }
}, by = .(catA, catB, N, n_shared)]

enrich_shared[, p_adj := p.adjust(p, "BH")]
enrich_shared[, log_odds := log(odds)]
enrich_shared[, sig := !is.na(p_adj) & p_adj <= 0.05]

women_AB <- net_w_05[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
men_AB   <- net_m_05[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
setnames(women_AB, "N", "n_women")
setnames(men_AB,   "N", "n_men")

enrich_sex <- merge(women_AB, men_AB, by = c("catA","catB"), all = TRUE)
enrich_sex[is.na(n_women), n_women := 0L]
enrich_sex[is.na(n_men),   n_men   := 0L]

total_women <- nrow(net_w_05[!is.na(catA) & !is.na(catB)])
total_men   <- nrow(net_m_05[!is.na(catA) & !is.na(catB)])

enrich_sex <- enrich_sex[, {
  a1 <- n_women; b1 <- total_women - a1
  a0 <- n_men;   b0 <- total_men   - a0
  if (any(c(a1,b1,a0,b0) < 0)) {
    .(odds = NA_real_, p = NA_real_)
  } else {
    ft <- fisher.test(matrix(c(a1,b1,a0,b0), nrow = 2, byrow = TRUE))
    .(odds = as.numeric(ft$estimate), p = ft$p.value)
  }
}, by = .(catA, catB, n_women, n_men)]

enrich_sex[, p_adj := p.adjust(p, "BH")]
enrich_sex[, log_odds := log(odds)]
enrich_sex[, sig := !is.na(p_adj) & p_adj <= 0.05]

plot_enrich_heatmap <- function(dt, low_col = "#C53030", high_col = "#2B6CB0",
                                exclude_cats = NULL) {
  dt2 <- copy(dt)
  if (!is.null(exclude_cats)) {
    dt2 <- dt2[!catA %in% exclude_cats & !catB %in% exclude_cats]
  }
  dt2[, catA_s := abbrev[catA]]
  dt2[, catB_s := abbrev[catB]]

  catA_levels <- sort(unique(dt2$catA_s))
  catB_levels <- sort(unique(dt2$catB_s))

  grid <- CJ(catA_s = catA_levels, catB_s = catB_levels)
  grid <- merge(grid, dt2[, .(catA_s, catB_s, log_odds, sig)],
                by = c("catA_s","catB_s"), all.x = TRUE)
  grid[, cell_status := fcase(
    is.na(log_odds), "No data",
    sig == TRUE,      "Significant",
    rep(TRUE, .N),     "Not significant"
  )]

  ggplot(grid, aes(x = catB_s, y = catA_s)) +
    geom_tile(data = grid[cell_status == "No data"],
              fill = "grey80", color = "white", linewidth = 0.3) +
    geom_tile(data = grid[cell_status != "No data"],
              aes(fill = log_odds), color = "white", linewidth = 0.3) +
    geom_text(data = grid[cell_status == "Significant"],
              label = "*", color = "black", fontface = "bold",
              size = 4, nudge_y = -0.08) +
    scale_fill_gradient2(
      midpoint = 0, low = low_col, mid = "white", high = high_col,
      na.value = "grey80",
      name   = "log(OR)\n(all tested\ncells; * = FDR\u22640.05)",
      limits = c(-3, 3), oob = scales::squish
    ) +
    labs(x = "Secondary disease category (B)", y = "Index disease category (A)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y     = element_text(size = 10),
      axis.title      = element_text(size = 14),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 10)
    )
}

p_shared <- plot_enrich_heatmap(enrich_shared, exclude_cats = c("Pregnancy", "Perinatal"))
p_sex    <- plot_enrich_heatmap(enrich_sex,    exclude_cats = c("Pregnancy", "Perinatal"))

cat("Combining panels (side by side)...\n")
p_combined <- wrap_plots(p_shared, p_sex, ncol = 2) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("ManuscriptFiles/Plots/Sex_network_enrichment_heatmaps.pdf", p_combined, width = 20, height = 11, useDingbats = FALSE)