# random_mixup_test.R
# Do the discordant WGS APOE genotypes look like RANDOM reassignment (sample /
# manifest mix-up) rather than structured single-SNP miscalls?
#
# Model under test ("random mix-up"): a wrong WGS genotype is an independent draw
# from the population genotype distribution pi, conditioned on != the true (Sanger)
# genotype. Expected off-diagonal confusion E[h->g] = D_h * pi_g/(1-pi_h).
# We test goodness-of-fit with a Monte Carlo permutation null and also compare a
# genotype-distance statistic (how many errors are 2-allele-substitutions away,
# which structured single-SNP miscalling cannot easily produce).

suppressMessages(library(tidyverse))
set.seed(1)

fc <- read.csv("../../Metadata_Merge/output_files/filtered/filtered_combined_post_QC.csv")

norm_g <- function(g) ifelse(is.na(g) | g == "", NA, gsub("[/ ]", "", trimws(g)))
d <- fc %>%
  mutate(sanger = norm_g(APOE.geno), wgs = norm_g(APOE_WGS)) %>%
  filter(!is.na(sanger), !is.na(wgs))

# Affected pool: Hispanic dual-typed (where ~all discordance lives = the plausibly
# scrambled set). pi = TRUE genotype distribution (Sanger) of that pool.
pool <- d %>% filter(Ethnicity == "HI")
genos <- c("23", "24", "33", "34", "44")
pool <- pool %>% filter(sanger %in% genos, wgs %in% genos)

pi <- prop.table(table(factor(pool$sanger, levels = genos)))
cat("Population genotype frequencies (Sanger truth, Hispanic pool, n =", nrow(pool), "):\n")
print(round(pi, 3))

disc <- pool %>% filter(sanger != wgs)
D <- nrow(disc)
cat("\nDiscordant in pool:", D, " (", round(100 * D / nrow(pool), 1), "% )\n", sep = "")

# Implied scramble fraction if the whole pool were partially permuted
r_max <- 1 - sum(pi^2)                       # discordance under FULL random shuffle
f_hat <- (D / nrow(pool)) / r_max
cat(sprintf("\nFull-shuffle discordance ceiling 1-sum(pi^2) = %.3f\n", r_max))
cat(sprintf("Implied fraction of pool scrambled to explain observed rate: f = %.2f\n", f_hat))

# ── Observed vs expected confusion (off-diagonal) under random-draw model ──────
O <- table(Sanger = factor(disc$sanger, genos), WGS = factor(disc$wgs, genos))
Dh <- rowSums(O)
E <- matrix(0, 5, 5, dimnames = list(genos, genos))
for (h in genos) {
  if (Dh[h] == 0) next
  pr <- pi; pr[h] <- 0; pr <- pr / sum(pr)   # draw != true genotype
  E[h, ] <- Dh[h] * pr
}
cat("\n=== Observed confusion (Sanger -> WGS), discordant ===\n"); print(O)
cat("\n=== Expected under random reassignment ===\n"); print(round(E, 1))

# G-test statistic on off-diagonal cells
offdiag <- which(row(O) != col(O))
Ov <- as.vector(O)[offdiag]; Ev <- as.vector(E)[offdiag]
keep <- Ev > 0
G_obs <- 2 * sum(ifelse(Ov[keep] > 0, Ov[keep] * log(Ov[keep] / Ev[keep]), 0))
cat(sprintf("\nObserved G (fit to random-draw model): %.2f\n", G_obs))

# ── Monte Carlo null: simulate the random-mix-up model, recompute G ───────────
B <- 20000
sim_G <- numeric(B)
sim_d2 <- numeric(B)

# genotype allele-substitution distance (min substitutions between allele multisets)
gdist <- function(a, b) {
  aa <- as.integer(strsplit(a, "")[[1]]); bb <- as.integer(strsplit(b, "")[[1]])
  # match best pairing
  min(sum(sort(aa) != sort(bb)),
      sum(sort(aa) != sort(rev(bb))))
}
dist_mat <- outer(genos, genos, Vectorize(gdist))
dimnames(dist_mat) <- list(genos, genos)

# observed share of errors that are distance-2 (two allele substitutions)
obs_d2 <- sum(O * (dist_mat == 2)) / D

draw_ne <- function(h, n) {            # n draws from pi excluding genotype h
  pr <- pi; pr[h] <- 0; pr <- pr / sum(pr)
  sample(genos, n, replace = TRUE, prob = pr)
}
true_h <- disc$sanger
for (b in 1:B) {
  sim_wgs <- character(D)
  for (h in unique(true_h)) {
    idx <- which(true_h == h)
    sim_wgs[idx] <- draw_ne(h, length(idx))
  }
  Os <- table(Sanger = factor(true_h, genos), WGS = factor(sim_wgs, genos))
  Ovs <- as.vector(Os)[offdiag]
  sim_G[b] <- 2 * sum(ifelse(Ovs[keep] > 0, Ovs[keep] * log(Ovs[keep] / Ev[keep]), 0))
  sim_d2[b] <- sum(Os * (dist_mat == 2)) / D
}

p_G <- mean(sim_G >= G_obs)
cat(sprintf("\nMonte Carlo p-value (G_obs vs random-mix-up null, B=%d): %.3f\n", B, p_G))
cat("  -> large p = data CONSISTENT with random reassignment; small p = structured (poor fit)\n")

cat(sprintf("\nDistance-2 (two-allele-substitution) error share:\n  observed = %.3f\n  random-mixup null = %.3f [%.3f, %.3f]\n",
            obs_d2, mean(sim_d2), quantile(sim_d2, .025), quantile(sim_d2, .975)))
cat("  (single-SNP MISCALLING would push observed well BELOW the null band; mix-up sits inside it)\n")

cat("\nDONE\n")
