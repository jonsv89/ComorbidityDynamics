## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
## Reuse the same ICD-10 category colour palette as Figures 1, 4 and 5
colcod <- c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
            "#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54",
            "#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
names(colcod) <- unique(cate)
catcol <- as.character(colcod)
names(catcol) <- catname

x_mid <- median(clock_2d$med_slope)
y_mid <- median(clock_2d$pct_persist)

p2d <- ggplot(clock_2d, aes(x = med_slope, y = pct_persist,
                            color = catA, size = n_total)) +
  geom_vline(xintercept = x_mid, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = y_mid, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  stat_ellipse(data = boot_results_v2,
               aes(x = med_slope, y = pct_persist, color = catA),
               level = 0.95, linewidth = 0.5, alpha = 0.6,
               inherit.aes = FALSE) +
  geom_point(alpha = 0.9) +
  geom_text_repel(aes(label = catA_short), size = 4,
                  color = "grey20", max.overlaps = 20,
                  box.padding = 0.3, show.legend = FALSE) +
  scale_color_manual(values = catcol, guide = "none") +
  scale_size_continuous(range = c(4, 12), name = "N comorbidities") +
  annotate("text", x = x_mid - 0.003, y = Inf,
           label = "Faster attenuation", hjust = 1, vjust = 1.5,
           color = "grey40", size = 3.5, fontface = "italic") +
  annotate("text", x = x_mid + 0.001, y = Inf,
           label = "Slower attenuation", hjust = 0, vjust = 1.5,
           color = "grey40", size = 3.5, fontface = "italic") +
  annotate("text", x = -Inf, y = y_mid + 0.3,
           label = "High persistent\nrisk", hjust = -0.1, vjust = 0,
           color = "grey40", size = 3.5, fontface = "italic") +
  annotate("text", x = -Inf, y = y_mid - 0.3,
           label = "Low persistent\nrisk", hjust = -0.1, vjust = 1,
           color = "grey40", size = 3.5, fontface = "italic") +
  labs(
    x = "Rate of RR attenuation (median slope, cumulative windows)",
    y = "% of comorbidities with persistent conditional risk (1-2 to 4-5 yr)",
    title = "The two-dimensional biological clock of multimorbidity",
    subtitle = paste0(
      "X-axis: how fast comorbidity strength attenuates over cumulative follow-up.\n",
      "Y-axis: proportion of disease pairs with persistent risk ",
      "even if not detected in early windows.\n",
      "Each point = one ICD-10 category (coloured as in Figs. 1, 4-5). ",
      "Size = number of comorbidity pairs.\n",
      "Ellipses = 95% bootstrap confidence regions per category (1,000 pair-level ",
      "resampling replicates);\n",
      "overlap across the reference lines indicates a continuous space ",
      "rather than discrete classes."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey40", size = 9)
  )

ggsave("ManuscriptFiles/Plots/Biological_clock_2D_continuous.pdf", p2d, width = 13, height = 9, useDingbats = FALSE)