## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

windows_cumul <- c("0_1","0_2","0_3","0_4","0_5")
windows_cont  <- c("1_2","2_3","3_4","4_5")

plot_comorbidity_clock <- function(pair_id, title_label, show_legend = TRUE) {
  da <- gsub("_.+","",pair_id)
  db <- gsub(".+_","",pair_id)

  cumul_ex <- rbindlist(lapply(windows_cumul, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt <- dt[disease_a == da & disease_b == db]
    if (nrow(dt) == 0) return(NULL)
    dt[, `:=`(window_num = as.integer(gsub("0_","",w)),
              type = "Cumulative (0 to k years)")]
    dt
  }))

  cond_ex <- rbindlist(lapply(windows_cont, function(w) {
    f <- sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt <- dt[disease_a == da & disease_b == db]
    if (nrow(dt) == 0) return(NULL)
    # Conditional windows positioned at their endpoint:
    # 1_2 -> 2, 2_3 -> 3, 3_4 -> 4, 4_5 -> 5
    dt[, `:=`(window_num = as.integer(gsub("_.+","",w)) + 1L,
              type = "Conditional (year k-1 to k)")]
    dt
  }))

  plot_dt <- rbindlist(list(cumul_ex, cond_ex), fill = TRUE)
  if (nrow(plot_dt) == 0) {
    stop(sprintf("No data found for pair %s -- check disease_a/disease_b codes and file paths.", pair_id))
  }

  plot_dt[, sig := lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]

  ggplot(plot_dt, aes(x = window_num, y = RR_shrunk,
                      color = type, fill = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey60") +
    geom_ribbon(aes(ymin = CI_low_RR_shrunk, ymax = CI_high_RR_shrunk),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(aes(linetype = type), linewidth = 1.2, na.rm = TRUE) +
    geom_point(aes(shape = sig), size = 3, na.rm = TRUE) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name   = "",
      drop   = FALSE,
      guide  = if (show_legend) "legend" else "none"
    ) +
    scale_linetype_manual(
      values = c("Cumulative (0 to k years)"   = "solid",
                 "Conditional (year k-1 to k)" = "dashed"),
      guide  = "none"
    ) +
    scale_color_manual(
      values = c("Cumulative (0 to k years)"   = "#1B6CA8",
                 "Conditional (year k-1 to k)" = "#B30000"),
      name   = "Window type",
      guide  = if (show_legend) "legend" else "none"
    ) +
    scale_fill_manual(
      values = c("Cumulative (0 to k years)"   = "#1B6CA8",
                 "Conditional (year k-1 to k)" = "#B30000"),
      guide  = "none"
    ) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("Year 1", "Year 2", "Year 3", "Year 4", "Year 5")
    ) +
    labs(
      x        = "Follow-up endpoint (years since index diagnosis)",
      y        = "Relative Risk (shrunk)",
      title    = title_label,
      subtitle = sprintf("%s \u2192 %s", da, db)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position  = if (show_legend) "bottom" else "none",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5)
    )
}

cat("Building panel 1/4: A63 -> A64 (Episodic)...\n")
p_ex1 <- plot_comorbidity_clock("A63_A64",
                                "Episodic\n(A63: Other predominantly sexually transmitted diseases, NEC\n\u2192 A64: Unspecified sexually transmitted disease)",
                                show_legend = FALSE)

cat("Building panel 2/4: F10 -> K70 (Chronic stable)...\n")
p_ex2 <- plot_comorbidity_clock("F10_K70",
                                "Chronic stable\n(F10: Mental and behavioural disorders due to use of alcohol\n\u2192 K70: Alcoholic liver disease)",
                                show_legend = FALSE)

cat("Building panel 3/4: J01 -> J32 (Chronic progressive)...\n")
p_ex3 <- plot_comorbidity_clock("J01_J32",
                                "Chronic progressive\n(J01: Acute sinusitis \u2192 J32: Chronic sinusitis)",
                                show_legend = FALSE)

cat("Building panel 4/4: R06 -> J96 (Late-emerging)...\n")
p_ex4 <- plot_comorbidity_clock("R06_J96",
                                "Late-emerging\n(R06: Abnormalities of breathing\n\u2192 J96: Respiratory failure, NEC)",
                                show_legend = TRUE)

cat("\nCombining panels...\n")
p_examples <- wrap_plots(p_ex1, p_ex2, p_ex3, p_ex4, ncol = 2) +
  plot_layout(axes = "collect", guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 13),
    plot.subtitle   = element_text(size = 11)
  )

ggsave("ManuscriptFiles/Plots/Biological_clock_examples_B.pdf",
       p_examples, width = 16, height = 11, useDingbats = FALSE)
cat("\nSaved: ManuscriptFiles/Plots/Biological_clock_examples_B.pdf\n")

examples_corrected <- c(
  "Purely episodic"     = "A63\u2192A64 (Other predominantly sexually transmitted diseases, NEC\u2192Unspecified sexually transmitted disease)",
  "Chronic stable"      = "F10\u2192K70 (Mental and behavioural disorders due to use of alcohol\u2192Alcoholic liver disease)",
  "Chronic progressive" = "J01\u2192J32 (Acute sinusitis\u2192Chronic sinusitis)",
  "Biphasic"            = "R06\u2192J96 (Abnormalities of breathing\u2192Respiratory failure, not elsewhere classified)"
)
cat("\nCorrected 'examples' labels (for the summary table elsewhere in the pipeline):\n")
print(examples_corrected)