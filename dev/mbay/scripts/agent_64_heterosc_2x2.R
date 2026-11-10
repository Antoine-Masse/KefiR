# agent_64 : la mise en defaut de D6 rejouee a 2 x 2. H0 vraie (aucun effet), effectifs et
# ecarts-types par case inegaux : n = (5, 5, 5, 40), sd = (4, 4, 4, 1). Part des jeux ou un
# terme (A, B ou A:B) est dit "oui" a P > .95 par le moteur a plusieurs facteurs (normale a
# variance commune), a cote d'aov (p < .05) et d'un test de Welch sur les cases (oneway.test).
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
n <- c(5, 5, 5, 40); sdv <- c(4, 4, 4, 1)
r <- t(vapply(seq_len(N), function(i) {
	set.seed(640000L + i)
	cell <- rep(1:4, n)
	d <- data.frame(A = factor(c("a1", "a2", "a1", "a2")[cell]), B = factor(c("b1", "b1", "b2", "b2")[cell]))
	d$y <- rnorm(length(cell), 0, sdv[cell])
	rr <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	av <- summary(aov(y ~ A * B, data = d))[[1]][1:3, "Pr(>F)"]
	c(any95 = any(rr$effects$p_h1 > 0.95), any75 = any(rr$effects$p_h1 > 0.75), aov = any(av < 0.05), strong = any(rr$effects$p_h1 > 0.95) && !isTRUE(rr$het_cells),
		welch = oneway.test(y ~ interaction(A, B), data = d)$p.value < 0.05)
}, numeric(5)))
out <- c(R.version.string, sprintf("agent_64 : 2 x 2, H0, n = (5,5,5,40), sd = (4,4,4,1), N = %d", N),
	sprintf("m.bay un terme P>.95 : %.3f | P>.75 : %.3f | aov un terme p<.05 : %.3f | Welch sur les cases p<.05 : %.3f | oui fort affiche (P>.95 sans plafond) : %.3f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 5]), mean(r[, 4])))
writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_64.txt")); cat(out, sep = "\n")
