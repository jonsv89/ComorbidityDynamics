## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

cat("=== Betweenness centrality sensitivity check ===\n\n")
SIG_FILTER <- quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)
HUBS <- c("M54", "J00", "T14")

cat("STEP 1: Loading 0-5-year significant network\n")

net <- fread("Results/shrinkage_rr_events/RR_contingency_both_0_5_shrunk.txt")
net <- net[eval(SIG_FILTER)]
cat(sprintf("  Significant pairs: %d\n", nrow(net)))

cat("\nSTEP 2: Hub removal and out-dangling-node pruning\n")

net_pruned <- net[!disease_a %in% HUBS & !disease_b %in% HUBS]
cat(sprintf("  After removing 3 hubs: %d pairs\n", nrow(net_pruned)))

repeat {
  outgoing_diseases <- unique(net_pruned$disease_a)
  before_n <- nrow(net_pruned)
  net_pruned <- net_pruned[disease_b %in% outgoing_diseases]
  if (nrow(net_pruned) == before_n) break
}

n_final <- length(unique(c(net_pruned$disease_a, net_pruned$disease_b)))
cat(sprintf("  Final network: %d diseases, %d pairs\n", n_final, nrow(net_pruned)))

cat("\nSTEP 3: Computing centrality measures\n")

g <- graph_from_data_frame(net_pruned[, .(disease_a, disease_b, weight = RR_shrunk)], directed = TRUE)

weighted_outdegree <- strength(g, mode = "out", weights = E(g)$weight)
pagerank_score <- page_rank(g, directed = TRUE, weights = E(g)$weight)$vector
betweenness_score <- betweenness(g, directed = TRUE, weights = 1 / E(g)$weight, normalized = TRUE)

centrality_dt <- data.table(
  disease = names(weighted_outdegree),
  outdegree = as.numeric(weighted_outdegree),
  pagerank = pagerank_score[names(weighted_outdegree)],
  betweenness = betweenness_score[names(weighted_outdegree)]
)

centrality_dt[, rank_outdegree := rank(-outdegree)]
centrality_dt[, rank_pagerank := rank(-pagerank)]
centrality_dt[, rank_betweenness := rank(-betweenness)]
centrality_dt[, discrepancy := rank_pagerank - rank_outdegree]  # published initiator/sink score

cat("\nSTEP 4: Does betweenness agree with the published axis?\n")

ct <- cor.test(centrality_dt$discrepancy, centrality_dt$betweenness, method = "spearman")
cat(sprintf("Spearman rho (discrepancy score vs. betweenness): %.3f (p=%.4g)\n",
            ct$estimate, ct$p.value))

ct2 <- cor.test(centrality_dt$rank_outdegree, centrality_dt$rank_betweenness, method = "spearman")
cat(sprintf("Spearman rho (outdegree rank vs. betweenness rank): %.3f (p=%.4g)\n",
            ct2$estimate, ct2$p.value))

cat("\nSTEP 5: Published initiators/sinks vs. betweenness\n")

N_TOP <- 19
N_BOT <- 21

published_initiators <- centrality_dt[order(-discrepancy)][1:N_TOP, disease]
published_sinks <- centrality_dt[order(discrepancy)][1:N_BOT, disease]

top_betweenness_initiators <- centrality_dt[disease %in% published_initiators]
top_betweenness_sinks <- centrality_dt[disease %in% published_sinks]

cat(sprintf("Median betweenness percentile of published initiators: %.1f\n",
            median(rank(centrality_dt$betweenness)[centrality_dt$disease %in% published_initiators]) /
              nrow(centrality_dt) * 100))
cat(sprintf("Median betweenness percentile of published sinks:      %.1f\n",
            median(rank(centrality_dt$betweenness)[centrality_dt$disease %in% published_sinks]) /
              nrow(centrality_dt) * 100))
cat("(50 = no different from a random disease in this network; further from 50 = distinctly high/low betweenness)\n")

## Top betweenness diseases NOT already flagged as initiator/sink
top_betweenness_novel <- centrality_dt[order(-betweenness)][1:20][
  !disease %in% c(published_initiators, published_sinks)]
cat(sprintf("\n  Of the top 20 diseases by betweenness, %d are NOT already published initiators/sinks:\n",
            nrow(top_betweenness_novel)))
print(top_betweenness_novel[, .(disease, betweenness = round(betweenness, 4), discrepancy)])

fwrite(centrality_dt, "ManuscriptFiles/Results/Centrality_betweenness_comparison.txt", sep = "\t", quote = FALSE)
cat("\Saved: Centrality_betweenness_comparison.txt\n")

