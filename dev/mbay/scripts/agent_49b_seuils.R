# agent_49b : relecture des cas 2^3 d'agent_49 (aucun nouveau calcul de BF).
# 1) Faux termes a P > 0.5, > 0.75 (au-dela de "faible") et > 0.95, par loi.
# 2) Cas A:B:C seule : sous la question "seul ou en interaction", A, B, C agissent (via
#    A:B:C) ; seuls A:B, A:C, B:C sont faux. On les recompte ainsi.
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_49_cas")
src <- readLines(file.path(BASE, "dev/mbay/scripts/agent_49_scenario_additif.R"))
fn <- function(name) {
	i <- grep(paste0("^", name, " <- function"), src)
	j <- i + which(src[(i + 1):length(src)] == "}")[1]
	eval(parse(text = src[i:j]), envir = globalenv())
}
fn("espace"); fn("inclusion")
LOIS <- c("half", "beta", "beta_ordre", "hop_ch", "taille", "unif")
sp <- espace(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
nm <- length(sp$mods)
cas <- list(c("n10_A-B-C_f50", "A B C"), c("n5_A-B-C-AB_f50", "A B C A:B"), c("n10_A-B-C-AB_f50", "A B C A:B"),
	c("n10_ABC_f50", "A B C A:B:C"))
out <- c(R.version.string, "Faux terme : P>.5 / P>.75 / P>.95 (au moins un dans le jeu), par loi. A:B:C seule : A, B, C comptes comme vrais (ils agissent via A:B:C).")
for (cc in cas) {
	f <- list.files(CAS, full.names = TRUE); f <- f[startsWith(basename(f), paste0("2x2x2_", cc[1], "_"))]
	r <- do.call(rbind, lapply(f, readRDS))
	vrai <- strsplit(cc[2], " ")[[1]]
	faux <- which(!sp$termes %in% vrai)
	li <- sprintf("%s (N=%d, faux : %s) |", cc[1], nrow(r), paste(sp$termes[faux], collapse = ","))
	for (lo in LOIS) {
		inc <- t(apply(r[, seq_len(nm), drop = FALSE], 1L, inclusion, sp = sp, lo = lo))
		mf <- apply(inc[, faux, drop = FALSE], 1L, max)
		li <- paste0(li, sprintf(" %s %.3f %.3f %.3f |", lo, mean(mf > .5), mean(mf > .75), mean(mf > .95)))
	}
	out <- c(out, li)
}
writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_49b.txt"))
cat(out, sep = "\n")
