## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "reconcile_universes") {
  library(data.table)
  valid_any <- readRDS("CaseControlStudy/valid_diseases_any.rds")
  valid_A   <- readRDS("CaseControlStudy/valid_diseases_A.rds")

  cat(sprintf("valid_diseases_any.rds:  %d enfermedades  (esperado: 1,498)\n", length(valid_any)))
  cat(sprintf("valid_diseases_A.rds:    %d enfermedades  (esperado: 1,315)\n", length(valid_A)))

  is_subset <- all(valid_A %in% valid_any)
  cat(sprintf("valid_diseases_A es subconjunto de valid_diseases_any: %s\n\n", is_subset))

  soren_raw <- fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",
                     stringsAsFactors = FALSE, sep = "\t")
  dis_den_universe <- unique(c(soren_raw$A, soren_raw$B))
  cat(sprintf("Universo Dinamarca (todas las enfermedades en el fichero Westergaard): %d\n",
              length(dis_den_universe)))

  dis_common <- intersect(valid_any, dis_den_universe)
  cat(sprintf("Interseccion valid_diseases_any x Dinamarca: %d  (esperado: 1,204)\n\n",
              length(dis_common)))

  f <- "Results/networks/RR_net_both_0_5.txt"
  if (!file.exists(f)) f <- "Results/shrinkage_rr_events/RR_contingency_both_0_5_shrunk.txt"
  net <- fread(f)
  net_sig <- net[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]

  all_nodes <- unique(c(net_sig$disease_a, net_sig$disease_b))
  cat(sprintf("Nodos unicos en la red significativa 0-5y (both): %d\n", length(all_nodes)))

  exclude <- c("M54","J00","T14")
  present <- exclude[exclude %in% all_nodes]
  cat(sprintf("De M54/J00/T14, cuales aparecen como nodo: %s\n",
              if (length(present)) paste(present, collapse=", ") else "ninguna"))

  nodes_after_exclusion <- setdiff(all_nodes, exclude)
  cat(sprintf("Nodos tras excluir M54/J00/T14: %d  (esperado: 845)\n\n", length(nodes_after_exclusion)))

  nodes_as_index <- intersect(all_nodes, valid_A)
  cat(sprintf("Nodos de la red 0-5y (both) que ADEMAS son validos como indice (valid_diseases_A): %d\n",
              length(nodes_as_index)))
  nodes_as_index_no_excl <- setdiff(nodes_as_index, exclude)
  cat(sprintf("...tras excluir M54/J00/T14 de esos: %d  (esperado: 845)\n\n", length(nodes_as_index_no_excl)))

  nodes_with_outgoing <- unique(net_sig$disease_a)
  cat(sprintf("Nodos con >=1 arista saliente significativa (outdegree real > 0): %d\n",
              length(nodes_with_outgoing)))
  nodes_outgoing_no_excl <- setdiff(nodes_with_outgoing, exclude)
  cat(sprintf("...tras excluir M54/J00/T14 de esos: %d  (esperado: 845)\n\n", length(nodes_outgoing_no_excl)))

  node_partners <- function(node) {
    as_a <- net_sig[disease_a == node, disease_b]
    as_b <- net_sig[disease_b == node, disease_a]
    unique(c(as_a, as_b))
  }
  isolated_after_removal <- vapply(nodes_after_exclusion, function(n) {
    partners <- node_partners(n)
    length(partners) > 0 && all(partners %in% exclude)
  }, logical(1))
  n_isolated <- sum(isolated_after_removal)
  final_after_isolation <- length(nodes_after_exclusion) - n_isolated
}