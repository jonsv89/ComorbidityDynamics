## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
library(data.table)
library(ggplot2)
library(patchwork)

all_windows_9 <- c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
window_type   <- c(rep("cumulative",5), rep("conditional",4))

all_9 <- rbindlist(lapply(seq_along(all_windows_9), function(i) {
  w  <- all_windows_9[i]
  wt <- window_type[i]
  f  <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt <- dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  dt[, `:=`(window = w, window_type = wt, pair = paste(disease_a, disease_b, sep = "_"))]
  dt
}))

pair_windows <- all_9[, .(
  in_0_1 = any(window == "0_1"), in_0_2 = any(window == "0_2"),
  in_0_3 = any(window == "0_3"), in_0_4 = any(window == "0_4"),
  in_0_5 = any(window == "0_5"), in_1_2 = any(window == "1_2"),
  in_2_3 = any(window == "2_3"), in_3_4 = any(window == "3_4"),
  in_4_5 = any(window == "4_5"),
  n_cumul = sum(window_type == "cumulative"),
  n_cond  = sum(window_type == "conditional")
), by = pair]

pair_windows[, bio_class := fcase(
  n_cumul == 5 & n_cond == 4, "Omnipresent (all 9 windows)",
  n_cumul == 5 & n_cond > 0,  "Persistent cumulative + partial conditional",
  n_cumul == 5 & n_cond == 0, "Persistent cumulative only (no conditional risk)",
  in_0_5 == TRUE & in_0_1 == FALSE & n_cond == 0, "Long-window only (late cumulative, no early signal)",
  in_0_5 == TRUE & in_0_1 == FALSE & n_cond > 0,  "Late cumulative + conditional risk",
  n_cumul == 0 & n_cond > 0, "Conditional only (no cumulative signal)",
  in_0_1 == TRUE & in_0_5 == FALSE & n_cond == 0, "Early transient (short window only)",
  in_0_1 == TRUE & in_0_5 == FALSE & n_cond > 0,  "Early cumulative + late conditional",
  rep(TRUE, .N), "Mixed/partial"
)]

class_counts <- pair_windows[, .N, by = bio_class][, pct := round(100 * N / sum(N), 1)][order(-N)]
cat("=== Nine-class breakdown (compare against manuscript) ===\n")
print(class_counts)

pair_windows[, narrative_group := fcase(
  bio_class %in% c("Omnipresent (all 9 windows)", "Persistent cumulative + partial conditional",
                    "Persistent cumulative only (no conditional risk)"), "Robust",
  bio_class %in% c("Long-window only (late cumulative, no early signal)",
                    "Late cumulative + conditional risk"), "Late-emerging",
  bio_class %in% c("Conditional only (no cumulative signal)", "Early transient (short window only)",
                    "Early cumulative + late conditional"), "Transient/conditional",
  rep(TRUE, .N), "Mixed/partial"
)]

four_group_summary <- pair_windows[, .N, by = narrative_group][, pct := round(100 * N / sum(N), 1)][order(-N)]
cat("\n=== Four-group summary (compare against Results paragraph) ===\n")
print(four_group_summary)
cat("Expected: Late-emerging 60.1% (88,140); Robust 38.1% (55,848);\n")
cat("          Transient/conditional 0.8% (1,100); Mixed/partial 1.0% (1,503)\n")

fwrite(class_counts, "nine_class_summary.csv")
fwrite(four_group_summary, "four_group_summary.csv")

group_of_class <- c(
  "Omnipresent (all 9 windows)" = "Robust",
  "Persistent cumulative + partial conditional" = "Robust",
  "Persistent cumulative only (no conditional risk)" = "Robust",
  "Long-window only (late cumulative, no early signal)" = "Late-emerging",
  "Late cumulative + conditional risk" = "Late-emerging",
  "Conditional only (no cumulative signal)" = "Transient/conditional",
  "Early transient (short window only)" = "Transient/conditional",
  "Early cumulative + late conditional" = "Transient/conditional",
  "Mixed/partial" = "Mixed/partial"
)
group_colors <- c("Late-emerging" = "#7F77DD", "Robust" = "#639922",
                   "Transient/conditional" = "#D4537E", "Mixed/partial" = "#B4B2A9")

windows <- c("0-1","0-2","0-3","0-4","0-5","1-2","2-3","3-4","4-5")
S <- "Significant"; N0 <- "Not significant"
P <- "Partial (>=1 of 4)"; U <- "Unconstrained"

grid_def <- list(
  "Omnipresent (all 9 windows)" = c(S,S,S,S,S, S,S,S,S),
  "Persistent cumulative + partial conditional" = c(S,S,S,S,S, P,P,P,P),
  "Persistent cumulative only (no conditional risk)" = c(S,S,S,S,S, N0,N0,N0,N0),
  "Long-window only (late cumulative, no early signal)" = c(N0,U,U,U,S, N0,N0,N0,N0),
  "Late cumulative + conditional risk" = c(N0,U,U,U,S, P,P,P,P),
  "Conditional only (no cumulative signal)" = c(N0,N0,N0,N0,N0, P,P,P,P),
  "Early transient (short window only)" = c(S,U,U,U,N0, N0,N0,N0,N0),
  "Early cumulative + late conditional" = c(S,U,U,U,N0, P,P,P,P),
  "Mixed/partial" = rep(U, 9)
)

class_order <- names(grid_def)
mat <- rbindlist(lapply(class_order, function(cl) {
  data.table(class = cl, window = windows, state = grid_def[[cl]])
}))
mat[, window := factor(window, levels = windows)]
mat[, class  := factor(class, levels = rev(class_order))]
mat[, state  := factor(state, levels = c(S, N0, P, U))]

label_dt <- copy(class_counts)
setnames(label_dt, "bio_class", "class")
label_dt[, class := factor(class, levels = rev(class_order))]
label_dt[, label := sprintf("%s (%s%%)", format(N, big.mark = ",", trim = TRUE), pct)]

strip_dt <- data.table(class_name = class_order, group = group_of_class[class_order])
strip_dt[, class := factor(class_name, levels = rev(class_order))]

panel_a <- ggplot(mat, aes(x = window, y = class, fill = state)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.8) +
  geom_tile(data = strip_dt, aes(x = -0.6, y = class, fill = NULL),
            inherit.aes = FALSE, width = 0.5, height = 0.8,
            fill = group_colors[strip_dt$group]) +
  geom_text(data = label_dt, aes(x = 9.9, y = class, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3.2) +
  scale_fill_manual(values = setNames(c("#1D9E75", "#E24B4A", "#F2A623", "#B4B2A9"),
                                       c(S, N0, P, U)), name = NULL) +
  scale_x_discrete(position = "top") +
  coord_cartesian(xlim = c(-1.1, 13), clip = "off") +
  annotate("text", x = 3, y = length(class_order) + 1.1,
           label = "Cumulative windows", size = 3.3, fontface = "bold") +
  annotate("text", x = 7.5, y = length(class_order) + 1.1,
           label = "Conditional windows", size = 3.3, fontface = "bold") +
  labs(x = NULL, y = NULL, title = "a) Nine temporal detection classes") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1, size = 8.5),
    axis.text.x.top = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(20, 5, 5, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

window_totals <- all_9[, .N, by = window][
  , window := factor(window, levels = all_windows_9)][order(window)]
window_totals[, window_lab := factor(gsub("_", "-", as.character(window)), levels = windows)]
window_totals[, type := rep(c("Cumulative","Conditional"), c(5,4))]

cat("\n=== Total significant associations per window (Figure Panel B) ===\n")
print(window_totals[, .(window, N)])
cat("Sanity check -- 0_1 should be 56,084 and 0_5 should be 144,030\n")
fwrite(window_totals, "window_totals_for_figure.csv")

panel_b <- ggplot(window_totals, aes(x = window_lab, y = N, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Cumulative" = "#378ADD", "Conditional" = "#D85A30"), name = NULL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Significant associations",
       title = "b) Total significant associations by window") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11)
  )

four_group_summary[, narrative_group := factor(narrative_group,
  levels = rev(c("Late-emerging","Robust","Transient/conditional","Mixed/partial")))]
four_group_summary[, label := sprintf("%s%%  (%s pairs)",
  pct, format(N, big.mark = ",", trim = TRUE))]

panel_c <- ggplot(four_group_summary, aes(x = pct, y = narrative_group, fill = narrative_group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.4) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "% of significant associations", y = NULL,
       title = "c) Summary: four narrative groups") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 11)
  )

combined <- panel_a / panel_b / panel_c + plot_layout(heights = c(1.3, 1, 0.6))

ggsave("Supplementary_Figure_nine_temporal_classes.pdf", combined, width = 9, height = 13, useDingbats = FALSE)