.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressPackageStartupMessages({library(BayesFactor); library(bayestestR)})
options(width = 200)

set.seed(4)
d <- data.frame(y = c(rnorm(12, 10), rnorm(12, 11.2), rnorm(12, 10.1), rnorm(12, 10)),
                g = factor(rep(c("L1","L2","L3","L4"), each = 12)))
lev <- levels(d$g)
cat("Group means:\n"); print(tapply(d$y, d$g, mean))

# 1. Partitions: restricted growth strings of length 4
rgs <- list(); rec <- function(v) {
  if (length(v) == 4) { rgs[[length(rgs) + 1]] <<- v; return(invisible()) }
  for (k in 1:(max(v) + 1)) rec(c(v, k))
}
rec(1L)
stopifnot(length(rgs) == 15)
lab <- function(v) paste(sapply(split(lev, v), paste, collapse = ""), collapse = "|")

# 2. log BF
logBF <- sapply(rgs, function(v) {
  if (max(v) == 1) return(0)
  h <- factor(v[as.integer(d$g)])
  dd <- data.frame(y = d$y, h = h)
  b <- lmBF(y ~ h, data = dd, rscaleFixed = 0.5)
  as.vector(b@bayesFactor$bf)
})
errBF <- sapply(rgs, function(v) {
  if (max(v) == 1) return(0)
  dd <- data.frame(y = d$y, h = factor(v[as.integer(d$g)]))
  lmBF(y ~ h, data = dd, rscaleFixed = 0.5)@bayesFactor$error
})

# 3. prior
alpha <- 1
ew <- sapply(rgs, function(v) { s <- tabulate(v); alpha^length(s) * prod(factorial(s - 1)) })
K <- sapply(rgs, max)
prior <- ifelse(K == 1, 0.5, 0.5 * ew / sum(ew[K > 1]))

# 4. posterior
lp <- log(prior) + logBF
post <- exp(lp - max(lp)); post <- post / sum(post)

tab <- data.frame(partition = sapply(rgs, lab), K = K, ewens_w = ew,
                  logBF = round(logBF, 4), err = signif(errBF, 2),
                  prior = round(prior, 5), posterior = round(post, 5))
tab <- tab[order(-tab$posterior), ]
cat("\nPartitions:\n"); print(tab, row.names = FALSE)

pairs <- combn(4, 2)
pdiff <- apply(pairs, 2, function(p) sum(post[sapply(rgs, function(v) v[p[1]] != v[p[2]])]))
names(pdiff) <- apply(pairs, 2, function(p) paste(lev[p], collapse = "-"))

# 5. effect sizes
set.seed(1)
ch <- posterior(lmBF(y ~ g, data = d, rscaleFixed = 0.5), iterations = 20000, progress = FALSE)
ch <- as.matrix(ch)
cat("\nChain columns:", colnames(ch), "\n")
eff <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
  a <- lev[pairs[1, j]]; b <- lev[pairs[2, j]]
  dr <- ch[, paste0("g-", a)] - ch[, paste0("g-", b)]
  de <- dr / sqrt(ch[, "sig2"])
  h1 <- hdi(dr, ci = 0.95); h2 <- hdi(de, ci = 0.95)
  data.frame(pair = paste(a, b, sep = "-"), P_diff = round(pdiff[j], 4),
             diff_med = round(median(dr), 3), diff_lo = round(h1$CI_low, 3), diff_hi = round(h1$CI_high, 3),
             delta_med = round(median(de), 3), delta_lo = round(h2$CI_low, 3), delta_hi = round(h2$CI_high, 3))
}))
cat("\nPairs:\n"); print(eff, row.names = FALSE)

# 6. checks
cat("\nSum of prior:", format(sum(prior), digits = 15), "\n")
cat("Sum of posterior:", format(sum(post), digits = 15), "\n")
P <- matrix(0, 4, 4)
for (j in seq_len(ncol(pairs))) { P[pairs[1, j], pairs[2, j]] <- pdiff[j]; P[pairs[2, j], pairs[1, j]] <- pdiff[j] }
viol <- 0; minslack <- Inf
for (i in 1:4) for (j in 1:4) for (k in 1:4) if (length(unique(c(i, j, k))) == 3) {
  s <- P[i, j] + P[j, k] - P[i, k]; minslack <- min(minslack, s); if (s < -1e-12) viol <- viol + 1
}
cat("Triangle checks (24 ordered triplets): violations =", viol, "; min slack =", round(minslack, 5), "\n")
cat("P(at least one difference) =", round(1 - post[K == 1], 5), "\n")
write.csv(tab, "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/relecture/agent_53/partitions.csv", row.names = FALSE)
write.csv(eff, "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/relecture/agent_53/paires.csv", row.names = FALSE)
