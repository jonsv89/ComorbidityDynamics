## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
library(data.table)

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

SIG <- function(dt) dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]

net_w <- fread("Results/networks/RR_net_women_0_5.txt")
net_m <- fread("Results/networks/RR_net_men_0_5.txt")

sig_w <- SIG(net_w)
sig_m <- SIG(net_m)

sig_w[, pair     := paste(disease_a, disease_b, sep = "_")]
sig_w[, pair_rev := paste(disease_b, disease_a, sep = "_")]
sig_m[, pair     := paste(disease_a, disease_b, sep = "_")]

rev_w_to_m <- merge(
  sig_w[preferred_direction == TRUE,
        .(disease_a, disease_b, pair_rev,
          RR_w = RR_shrunk, theta_w = theta, n_w = cases_event)],
  sig_m[preferred_direction == TRUE,
        .(pair, RR_m = RR_shrunk, theta_m = theta, n_m = cases_event)],
  by.x = "pair_rev", by.y = "pair"
)
rev_w_to_m[, direction := "A->B en mujeres, B->A en hombres"]

sig_m[, pair_rev := paste(disease_b, disease_a, sep = "_")]
rev_m_to_w <- merge(
  sig_m[preferred_direction == TRUE,
        .(disease_a, disease_b, pair_rev,
          RR_m = RR_shrunk, theta_m = theta, n_m = cases_event)],
  sig_w[preferred_direction == TRUE,
        .(pair, RR_w = RR_shrunk, theta_w = theta, n_w = cases_event)],
  by.x = "pair_rev", by.y = "pair"
)
rev_m_to_w[, direction := "A->B en hombres, B->A en mujeres"]

all_rev <- rbindlist(list(
  rev_w_to_m[, .(disease_a, disease_b, RR_w, theta_w, n_w, RR_m, theta_m, n_m, direction)],
  rev_m_to_w[, .(disease_a, disease_b, RR_w, theta_w, n_w, RR_m, theta_m, n_m, direction)]
))

all_rev[, catA := catname[cate[disease_a]]]
all_rev[, catB := catname[cate[disease_b]]]
all_rev[, min_n := pmin(n_w, n_m)]
all_rev[, min_RR := pmin(RR_w, RR_m)]  # RR minimo entre ambos sexos, en su propia direccion preferida

print(all_rev[order(-min_n)][1:15,
      .(disease_a, catA, disease_b, catB, direction,
        RR_w = round(RR_w,2), n_w, RR_m = round(RR_m,2), n_m)])

print(all_rev[catA %in% c("Blood/Immune","Skin","Symptoms") | catB %in% c("Blood/Immune","Skin","Symptoms")][
      order(-min_n)][1:15,
      .(disease_a, catA, disease_b, catB, direction,
        RR_w = round(RR_w,2), n_w, RR_m = round(RR_m,2), n_m)])

fwrite(all_rev[order(-min_n)], "sex_directional_reversals_verified.csv")