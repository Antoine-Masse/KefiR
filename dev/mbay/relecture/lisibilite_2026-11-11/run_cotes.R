.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths())); suppressMessages(devtools::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE))
options(width = 80)
# Audit de lisibilite du 11/11/2026 : m.test puis m.bay sur 13 plans de TP.
# Usage : Rscript run_cotes.R [plan1 plan2 ...]  (sans argument : tous les plans)
dossier <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/relecture/lisibilite_2026-11-11"

plans <- list()

plans$g2 <- function() {
	set.seed(101)
	d <- data.frame(g = factor(rep(c("temoin", "traite"), each = 10)))
	d$y <- round(c(rnorm(10, 12, 2), rnorm(10, 14.5, 2)), 1)
	list(d = d, mt = quote(m.test(y ~ g, data = d)), mb = quote(m.bay(y ~ g, data = d)))
}
plans$g4 <- function() {
	set.seed(102)
	d <- data.frame(souche = factor(rep(c("A", "B", "C", "D"), each = 8)))
	d$y <- round(rnorm(32, c(20, 20, 23, 25)[as.integer(d$souche)], 2.5), 1)
	list(d = d, mt = quote(m.test(y ~ souche, data = d)), mb = quote(m.bay(y ~ souche, data = d)))
}
plans$f2x3 <- function() {
	set.seed(103)
	d <- expand.grid(rep = 1:5, milieu = c("M1", "M2"), temp = c("20C", "30C", "37C"))
	d$milieu <- factor(d$milieu); d$temp <- factor(d$temp)
	d$y <- round(rnorm(nrow(d), c(5, 7, 8.5)[as.integer(d$temp)], 1.2), 2)
	list(d = d, mt = quote(m.test(y ~ milieu * temp, data = d)), mb = quote(m.bay(y ~ milieu * temp, data = d)))
}
plans$f2x2int <- function() {
	set.seed(104)
	d <- expand.grid(rep = 1:8, engrais = c("sans", "avec"), lumiere = c("faible", "forte"))
	d$engrais <- factor(d$engrais); d$lumiere <- factor(d$lumiere)
	m <- with(d, 10 + 0 * (engrais == "avec") + 1 * (lumiere == "forte") + 4 * (engrais == "avec" & lumiere == "forte"))
	d$y <- round(rnorm(nrow(d), m, 1.8), 1)
	list(d = d, mt = quote(m.test(y ~ engrais * lumiere, data = d)), mb = quote(m.bay(y ~ engrais * lumiere, data = d)))
}
plans$ancova <- function() {
	set.seed(105)
	d <- data.frame(regime = factor(rep(c("R1", "R2", "R3"), each = 10)))
	d$poids0 <- round(rnorm(30, 250, 20))
	d$poids <- round(d$poids0 * 1.1 + c(0, 8, 15)[as.integer(d$regime)] + rnorm(30, 0, 8))
	list(d = d, mt = quote(m.test(poids ~ regime + poids0, data = d)), mb = quote(m.bay(poids ~ regime + poids0, data = d)))
}
plans$apparie <- function() {
	set.seed(106)
	n <- 12; base <- rnorm(n, 80, 8)
	d <- data.frame(id = factor(rep(1:n, 2)), t = factor(rep(c("avant", "apres"), each = n), levels = c("avant", "apres")))
	d$y <- round(c(base + rnorm(n, 0, 2), base - 3 + rnorm(n, 0, 2)), 1)
	list(d = d, mt = quote(m.test(y ~ t | id, data = d)), mb = quote(m.bay(y ~ t | id, data = d)))
}
plans$repet3 <- function() {
	set.seed(107)
	n <- 10; base <- rnorm(n, 50, 6)
	d <- data.frame(id = factor(rep(1:n, 3)), t = factor(rep(c("J0", "J7", "J14"), each = n), levels = c("J0", "J7", "J14")))
	d$y <- round(rep(base, 3) + rep(c(0, 2, 4), each = n) + rnorm(3 * n, 0, 2.5), 1)
	list(d = d, mt = quote(m.test(y ~ t | id, data = d)), mb = quote(m.bay(y ~ t | id, data = d)))
}
plans$tranches <- function() {
	set.seed(108)
	n <- 8
	d <- expand.grid(temps = c("T1", "T2", "T3"), id = 1:(2 * n))
	d$groupe <- factor(ifelse(d$id <= n, "placebo", "actif")); d$id <- factor(d$id)
	d$temps <- factor(d$temps)
	base <- rnorm(2 * n, 30, 4)
	d$y <- round(base[as.integer(d$id)] + (d$groupe == "actif") * c(0, 2, 4)[as.integer(d$temps)] +
		rnorm(nrow(d), 0, 1.5), 1)
	list(d = d, mt = quote(m.test(y ~ groupe * temps + Error(id/temps), data = d)),
		mb = quote(m.bay(y ~ groupe * temps + Error(id/temps), data = d)))
}
plans$ouinon1 <- function() {
	set.seed(109)
	d <- data.frame(lot = factor(rep(c("L1", "L2", "L3"), each = 30)))
	d$germe <- factor(ifelse(runif(90) < c(0.5, 0.6, 0.8)[as.integer(d$lot)], "oui", "non"))
	list(d = d, mt = quote(m.test(germe ~ lot, data = d)), mb = quote(m.bay(germe ~ lot, data = d)))
}
plans$comptages <- function() {
	set.seed(110)
	d <- data.frame(boite = factor(rep(c("B1", "B2", "B3"), each = 8)))
	d$colonies <- rpois(24, c(12, 15, 22)[as.integer(d$boite)])
	list(d = d, mt = quote(m.test(colonies ~ boite, data = d)), mb = quote(m.bay(colonies ~ boite, data = d)))
}
plans$notes15 <- function() {
	set.seed(111)
	d <- data.frame(recette = factor(rep(c("R1", "R2", "R3"), each = 15)))
	p <- list(c(.1, .25, .35, .2, .1), c(.1, .2, .3, .25, .15), c(.02, .08, .2, .4, .3))
	d$note <- unlist(lapply(1:3, function(k) sample(1:5, 15, TRUE, prob = p[[k]])))
	d$note <- factor(d$note, levels = 1:5, ordered = TRUE)
	list(d = d, mt = quote(m.test(note ~ recette, data = d)), mb = quote(m.bay(note ~ recette, data = d)))
}
plans$ouinon2x2 <- function() {
	set.seed(112)
	d <- expand.grid(rep = 1:25, souche = c("S1", "S2"), milieu = c("riche", "pauvre"))
	d$souche <- factor(d$souche); d$milieu <- factor(d$milieu)
	pr <- with(d, 0.3 + 0.3 * (souche == "S2") + 0.1 * (milieu == "riche"))
	d$pousse <- factor(ifelse(runif(nrow(d)) < pr, "oui", "non"))
	list(d = d, mt = quote(m.test(pousse ~ souche * milieu, data = d)), mb = quote(m.bay(pousse ~ souche * milieu, data = d)))
}
plans$notes15num <- function() {
	# Variante : les notes saisies comme nombres entiers 1 a 5 (cas le plus courant en TP).
	p <- plans$notes15()
	p$d$note <- as.integer(as.character(p$d$note))
	p
}
plans$manova <- function() {
	set.seed(113)
	d <- data.frame(g = factor(rep(c("A", "B", "C"), each = 10)))
	d$y1 <- round(rnorm(30, c(10, 11, 13)[as.integer(d$g)], 1.5), 2)
	d$y2 <- round(0.5 * d$y1 + rnorm(30, c(5, 5, 6)[as.integer(d$g)], 1), 2)
	list(d = d, mt = quote(m.test(cbind(y1, y2) ~ g, data = d)), mb = quote(m.bay(cbind(y1, y2) ~ g, data = d)))
}

lancer <- function(nom) {
	p <- plans[[nom]]()
	d <- p$d
	f <- file.path(dossier, paste0("sorties_", nom, ".txt"))
	con <- file(f, open = "wt", encoding = "UTF-8")
	sink(con); sink(con, type = "message")
	on.exit({ sink(type = "message"); sink(); close(con) })
	cat("############ PLAN :", nom, "############\n\n")
	cat("== Donnees (head) ==\n"); print(head(d)); cat("n =", nrow(d), "\n\n")
	cat("============ m.test ============\n")
	cat(deparse(p$mt), "\n\n")
	t0 <- Sys.time()
	r <- tryCatch({ invisible(eval(p$mt)); "OK" }, error = function(e) paste("ERREUR m.test :", conditionMessage(e)))
	cat("\n[m.test :", r, "-", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s]\n\n")
	cat("============ m.bay ============\n")
	cat(deparse(p$mb), "\n\n")
	t0 <- Sys.time()
	r <- tryCatch({ invisible(eval(p$mb)); "OK" }, error = function(e) paste("ERREUR m.bay :", conditionMessage(e)))
	dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
	cat("\n[m.bay :", r, "-", round(dt, 1), "s]\n")
	if (dt > 60) {
		cat("\n[> 60 s : nouveau lancement avec sensitivity = FALSE]\n\n")
		mb2 <- p$mb; mb2$sensitivity <- FALSE
		t0 <- Sys.time()
		r <- tryCatch({ invisible(eval(mb2)); "OK" }, error = function(e) paste("ERREUR m.bay :", conditionMessage(e)))
		cat("\n[m.bay sensitivity = FALSE :", r, "-", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s]\n")
	}
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- names(plans)
for (nom in args) {
	message("Plan ", nom, " ...")
	lancer(nom)
	message("  fini.")
}
