## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dt <- fread("ManuscriptFiles/Results/Power_anchored_part2_thinning_v2.txt")
cat(sprintf("Loaded %d pairs\n", nrow(dt)))

med <- median(dt$pct_still_significant, na.rm = TRUE)
pct_reliable <- mean(dt$pct_still_significant >= 80, na.rm = TRUE) * 100

p <- ggplot(dt, aes(x = pct_still_significant)) +
  geom_histogram(binwidth = 2, boundary = 0, fill = "#1B6CA8", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = med, linetype = "dashed", color = "#B30000", linewidth = 0.7) +
  annotate("text", x = med + 3, y = Inf, label = sprintf("Median: %.1f%%", med),
           hjust = 0, vjust = 1.5, color = "#B30000", size = 4, fontface = "bold") +
  annotate("rect", xmin = 80, xmax = 100, ymin = 0, ymax = Inf,
           fill = "#2C6E3F", alpha = 0.15) +
  annotate("text", x = 90, y = Inf, label = sprintf("Reliably detectable\n(>=80%%): %.1f%%", pct_reliable),
           hjust = 0.5, vjust = 1.3, color = "#2C6E3F", size = 3.5, fontface = "bold") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "% of binomial-thinning replicates remaining significant (of 500)",
    y = "Number of late-emerging pairs",
    title = "Detectability of late-emerging associations at the minimum event threshold",
    subtitle = paste0(
      "Each pair's final (0-5-year) effect size, thinned to 100 events (the pipeline's own\n",
      "minimum detection threshold) and re-tested for significance across 500 replicates.\n",
      "N = ", format(nrow(dt), big.mark = ","), " late-emerging pairs."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

ggsave("ManuscriptFiles/Plots/Power_anchored_thinning_histogram.pdf", p, width = 10, height = 7, useDingbats = FALSE)
cat("\nSaved: Power_anchored_thinning_histogram.pdf\n")