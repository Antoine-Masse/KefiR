# agent_21 : heteroscedasticite sous H0, n=(5,5,40), sd=(4,4,1) et (1,1,4).
# Eclaire l'ecart 16 % (agent_12, N=300, graine 4242) / 25,5 % (r3c, N=200, graine 43).
# Usage : Rscript agent_21_heterosc_grand_N.R <N> <graine>
# Chaque jeu est regenerable : set.seed(graine_jeu) ; rnorm(50, 0, rep(sd, n)) avec graine_jeu dans le .rds.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress = FALSE)
args <- commandArgs(TRUE); N <- as.integer(args[1]); seed <- as.integer(args[2])
sortie <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats"
n <- c(5, 5, 40); g <- factor(rep(1:3, n))
un <- function(y) {
	b <- anovaBF(y ~ g, data = data.frame(y, g), progress = FALSE)@bayesFactor
	P <- plogis(b$bf)
	c(log_bf = b$bf, err = b$error, P = P, p_aov = anova(lm(y ~ g))[1, 5],
		p_welch = oneway.test(y ~ g)$p.value, p_kw = kruskal.test(y, g)$p.value)
}
sim <- function(s, base) {
	seeds <- base + seq_len(N)
	m <- t(vapply(seeds, function(sd0) { set.seed(sd0); y <- rnorm(sum(n), 0, rep(s, n)); un(y) }, numeric(6)))
	data.frame(jeu = seq_len(N), graine_jeu = seeds, m)
}
ic <- function(k, m) { b <- binom.test(k, m)$conf.int; sprintf("%.3f [%.3f ; %.3f] (%d/%d)", k / m, b[1], b[2], k, m) }

t0 <- Sys.time()
A <- sim(c(4, 4, 1), seed); ta <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
B <- sim(c(1, 1, 4), seed + 1000000L); tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Reproductions exactes des deux chaines d'origine
t1 <- Sys.time()
# r3c.R l.31, 40-41 : memes appels, memes graines
Pa <- function(y, f) { b <- extractBF(anovaBF(y ~ f, data.frame(y, f), progress = FALSE))$bf; c(P = b / (1 + b), p = summary(aov(y ~ f))[[1]][1, 5], kw = kruskal.test(y, f)$p.value) }
f2 <- factor(rep(1:3, c(5, 5, 40))); set.seed(43)
r3c <- t(replicate(200, { y <- rnorm(50, 0, c(4, 4, 1)[f2]); Pa(y, f2) }))
# agent_12 l.2-4 : memes appels, meme graine
set.seed(4242); s <- c(4, 4, 1)
a12 <- t(replicate(300, { y <- rnorm(sum(n), 0, rep(s, n)); bf <- exp(anovaBF(y ~ g, data = data.frame(y, g), progress = FALSE)@bayesFactor$bf)
	c(P = bf / (1 + bf), pa = anova(lm(y ~ g))[1, 5], pw = oneway.test(y ~ g)$p.value, pk = kruskal.test(y, g)$p.value) }))

# Meme jeu, les deux codes : difference de calcul ?
set.seed(777)
cmp <- t(replicate(200, { y <- rnorm(50, 0, c(4, 4, 1)[f2]); a <- Pa(y, f2); bf <- exp(anovaBF(y ~ g, data = data.frame(y, g), progress = FALSE)@bayesFactor$bf)
	c(a["P"] - bf / (1 + bf), a["p"] - anova(lm(y ~ g))[1, 5], a["kw"] - kruskal.test(y, g)$p.value) }))
tc <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

saveRDS(list(plan = list(n = n, sd_A = c(4, 4, 1), sd_B = c(1, 1, 4)), graine = seed, regeneration = "set.seed(graine_jeu); rnorm(50, 0, rep(sd, n))",
	A = A, B = B, r3c = r3c, a12 = a12, R = R.version.string, BayesFactor = as.character(packageVersion("BayesFactor"))),
	file.path(sortie, "agent_21_heterosc.rds"))

zz <- file(file.path(sortie, "agent_21_heterosc.txt"), "wt", encoding = "UTF-8"); sink(zz)
cat(R.version.string, "; BayesFactor", as.character(packageVersion("BayesFactor")), "; graine", seed, "; N", N, "\n")
cat("Duree : sd(4,4,1)", round(ta), "s ; total deux plans", round(tb), "s ; reproductions", round(tc), "s\n")
for (nm in c("A", "B")) {
	d <- get(nm); ok <- !is.na(d$P)
	cat("\nPlan", if (nm == "A") "sd=(4,4,1)" else "sd=(1,1,4)", "; NA :", sum(!ok), "; erreur mediane", signif(median(d$err, na.rm = TRUE), 3), "max", signif(max(d$err, na.rm = TRUE), 3), "\n")
	cat(" P>0,5   :", ic(sum(d$P[ok] > .5), sum(ok)), "\n P>0,95  :", ic(sum(d$P[ok] > .95), sum(ok)), "\n")
	cat(" P mediane :", round(median(d$P, na.rm = TRUE), 3), "\n")
	cat(" aov<.05 :", ic(sum(d$p_aov < .05), N), "\n Welch<.05 :", ic(sum(d$p_welch < .05), N), "\n Kruskal<.05 :", ic(sum(d$p_kw < .05), N), "\n")
}
cat("\nReproduction r3c (graine 43, N=200) : P>0,5", mean(r3c[, 1] > .5), "; P>0,95", mean(r3c[, 1] > .95), "; aov", mean(r3c[, 2] < .05), "; KW", mean(r3c[, 3] < .05), "; medP", round(median(r3c[, 1]), 3), "\n")
cat("Reproduction agent_12 (graine 4242, N=300) : P>0,5", mean(a12[, 1] > .5), "; P>0,95", mean(a12[, 1] > .95), "; aov", mean(a12[, 2] < .05), "; Welch", mean(a12[, 3] < .05), "; KW", mean(a12[, 4] < .05), "; medP", round(median(a12[, 1]), 3), "\n")
cat("\nMeme jeu, deux codes (200 jeux), |difference| maximale : P", signif(max(abs(cmp[, 1])), 3), "; p aov", signif(max(abs(cmp[, 2])), 3), "; p KW", signif(max(abs(cmp[, 3])), 3), "\n")
pt <- function(lbl, k1, k2) { r <- prop.test(c(k1, k2), c(300, 200), correct = TRUE); f <- fisher.test(matrix(c(k1, 300 - k1, k2, 200 - k2), 2)); cat(sprintf(" %-10s %d/300 contre %d/200 : prop.test p=%.4f ; Fisher p=%.4f\n", lbl, k1, k2, r$p.value, f$p.value)) }
cat("\nDeux proportions, agent_12 contre r3c :\n")
pt("P>0,95", sum(a12[, 1] > .95), sum(r3c[, 1] > .95)); pt("P>0,5", sum(a12[, 1] > .5), sum(r3c[, 1] > .5)); pt("aov", sum(a12[, 2] < .05), sum(r3c[, 2] < .05)); pt("KW", sum(a12[, 4] < .05), sum(r3c[, 3] < .05))
cat("\nChacun contre le grand N (binomial exact, taux de reference = estimation sur A) :\n")
pr <- mean(A$P > .95)
cat(sprintf(" P>0,95 : ref %.4f ; agent_12 p=%.4f ; r3c p=%.4f\n", pr, binom.test(sum(a12[, 1] > .95), 300, pr)$p.value, binom.test(sum(r3c[, 1] > .95), 200, pr)$p.value))
pa <- mean(A$p_aov < .05)
cat(sprintf(" aov    : ref %.4f ; agent_12 p=%.4f ; r3c p=%.4f\n", pa, binom.test(sum(a12[, 2] < .05), 300, pa)$p.value, binom.test(sum(r3c[, 2] < .05), 200, pa)$p.value))
cat("\nDispersion par lots du grand N, P>0,95 : lots de 300 et de 200\n")
for (m in c(300, 200)) { l <- split(A$P > .95, ceiling(seq_len(N) / m)); l <- l[lengths(l) == m]; v <- vapply(l, mean, 0); cat(sprintf(" lots de %d (%d lots) : min %.3f, max %.3f, ET %.3f\n", m, length(v), min(v), max(v), sd(v))) }
sink(); close(zz)
cat("FIN\n")
