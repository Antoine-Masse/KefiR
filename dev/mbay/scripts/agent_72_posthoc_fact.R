# agent_72 : calibration du post-hoc a plusieurs facteurs (D20). Plan 3 x 2 (A a 3 niveaux,
# B a 2), n = 10 par case, bruit normal. Scenarios :
#   H0  : aucun effet ; part des jeux avec une paire dite « differente » (P >= .75, >= .95).
#   A   : effet de A seul (a3 = +1) ; famille choisie = niveaux de A ? a1/a3 et a2/a3 dits
#         differents, a1/a2 non ?
#   AB  : interaction (a3 x b2 = +1,5, rien ailleurs) ; famille choisie = cases A x B ?
#         paires impliquant a3 x b2 dites differentes, les autres non ?
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 60L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(720000L + i + c(H0 = 0L, A = 1000L, AB = 2000L)[[sc]])
	d <- expand.grid(rep = 1:10, A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")))
	mu <- switch(sc, H0 = 0, A = 1 * (d$A == "a3"), AB = 1.5 * (d$A == "a3" & d$B == "b2"))
	d$y <- mu + stats::rnorm(nrow(d))
	rr <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	ph <- rr$posthoc
	fam <- if (length(ph$families)) vapply(ph$families, `[[`, character(1), "term") else character(0)
	pr <- ph$pairs
	any75 <- !is.null(pr) && any(pr$p_differ >= 0.75)
	any95 <- !is.null(pr) && any(pr$p_differ >= 0.95)
	if (sc == "A") {
		ok_fam <- identical(fam, "A")
		pa <- if (ok_fam) ph$families[[1L]]$pairs else NULL
		hit <- if (ok_fam) mean(pa$p_differ[pa$group1 != "a3" & pa$group2 == "a3" | pa$group1 == "a3"] >= 0.75) else NA
		fa <- if (ok_fam) pa$p_differ[pa$group1 == "a1" & pa$group2 == "a2"] >= 0.75 else NA
		return(c(ok_fam, hit, fa, any95, rr$effects$p_h1[1L] >= 0.75))
	}
	if (sc == "AB") {
		ok_fam <- "A:B" %in% fam
		pa <- if (ok_fam) ph$families[[which(fam == "A:B")]]$pairs else NULL
		inv <- if (ok_fam) grepl("a3 x b2", pa$group1, fixed = TRUE) | grepl("a3 x b2", pa$group2, fixed = TRUE) else NULL
		hit <- if (ok_fam) mean(pa$p_differ[inv] >= 0.75) else NA
		fa <- if (ok_fam) mean(pa$p_differ[!inv] >= 0.75) else NA
		return(c(ok_fam, hit, fa, any95, rr$effects$p_h1[3L] >= 0.75))
	}
	c(length(fam) > 0L, any75, any95, NA, NA)
}, numeric(5)))
lab <- switch(sc,
	H0 = "une famille detaillee : %.3f | une paire P>=.75 : %.3f | P>=.95 : %.3f",
	A = "famille = niveaux de A : %.3f | paires avec a3 dites differentes : %.3f | a1/a2 dite differente : %.3f | une paire P>=.95 : %.3f | A P>=.75 : %.3f",
	AB = "famille = cases A x B : %.3f | paires avec a3 x b2 dites differentes : %.3f | autres paires dites differentes : %.3f | une paire P>=.95 : %.3f | A:B P>=.75 : %.3f")
vals <- colMeans(r, na.rm = TRUE)
out <- c(R.version.string, sprintf("agent_72 : post-hoc a plusieurs facteurs, 3 x 2, n = 10, scenario %s, N = %d", sc, N),
	do.call(sprintf, c(list(lab), as.list(vals[seq_len(lengths(regmatches(lab, gregexpr("%", lab))))]))))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_72_%s.txt", sc))); cat(out, sep = "\n")
