#!/usr/bin/env Rscript
# Synthetic illustrations of textbook statistics, NOT biological observations.
f <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(f)), "_common.R"))
out <- output_dir("statistics")
set.seed(1234)

pois <- rpois(10000, lambda = 20)
nb <- rnbinom(10000, mu = 20, size = 1 / 0.1)
write.csv(data.frame(distribution = c("Poisson", "NB"), mean = c(mean(pois), mean(nb)),
                     variance = c(var(pois), var(nb))),
          file.path(out, "count_distributions.csv"), row.names = FALSE)
png(file.path(out, "count_distributions.png"), width = 1100, height = 450)
par(mfrow = c(1, 2))
hist(pois, breaks = 40, main = "SIMULATION: Poisson", xlab = "count", col = "steelblue")
hist(nb, breaks = 40, main = "SIMULATION: negative binomial", xlab = "count", col = "coral")
dev.off()

# Null simulation: cells share a donor effect, with no true group effect.
sim <- replicate(300, {
  donor <- rep(seq_len(12), each = 100)
  group_donor <- rep(c("control", "treated"), each = 6)
  y <- rnorm(12, sd = 1)[donor] + rnorm(1200, sd = 1)
  donor_means <- tapply(y, donor, mean)
  c(naive_cell = t.test(y ~ group_donor[donor])$p.value,
    donor_mean = t.test(donor_means ~ group_donor)$p.value)
})
write.csv(t(sim), file.path(out, "pseudoreplication.csv"), row.names = FALSE)
write.csv(data.frame(method = rownames(sim), null_rejection_rate = rowMeans(sim < 0.05)),
          file.path(out, "pseudoreplication_summary.csv"), row.names = FALSE)
png(file.path(out, "pseudoreplication.png"), width = 1000, height = 450)
par(mfrow = c(1, 2))
for (i in 1:2) hist(sim[i, ], breaks = seq(0, 1, 0.05), main = rownames(sim)[i], xlab = "Null p-value")
dev.off()

p <- c(0.001, 0.02, 0.03, 0.2, 0.8)
write.csv(data.frame(p = p, BH = p.adjust(p, "BH")), file.path(out, "bh_example.csv"), row.names = FALSE)
condition <- factor(rep(c("control", "treated"), each = 4))
confounded <- model.matrix(~ factor(condition) + condition)
balanced <- model.matrix(~ factor(rep(c("A", "B"), 4)) + condition)
write.csv(data.frame(design = c("confounded", "balanced"),
                     rank = c(qr(confounded)$rank, qr(balanced)$rank),
                     columns = c(ncol(confounded), ncol(balanced))),
          file.path(out, "design_rank.csv"), row.names = FALSE)

# Left-censoring exercise. Same simulated abundance distribution in both groups,
# different detection limits create an apparent difference.
group <- rep(c("control", "treated"), each = 30)
truth <- rlnorm(60, meanlog = 3, sdlog = 0.6)
observed <- truth
observed[truth < ifelse(group == "treated", 18, 10)] <- NA_real_
half_min <- observed
half_min[is.na(half_min)] <- min(observed, na.rm = TRUE) / 2
calc <- function(x) c(log2FC = mean(log2(x[group == "treated"]), na.rm = TRUE) -
                       mean(log2(x[group == "control"]), na.rm = TRUE),
                     p = t.test(log2(x) ~ group)$p.value)
sens <- rbind(complete_truth = calc(truth), observed_complete_case = calc(observed),
              half_minimum = calc(half_min))
write.csv(sens, file.path(out, "missingness_sensitivity.csv"))
record_session(out, "SIMULATED DATA ONLY; seed=1234; 300 null experiments, 12 donors x 100 cells.")
cat("Finished synthetic statistics lab:", out, "\n")
