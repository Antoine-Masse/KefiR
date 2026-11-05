# Audit systematique de m.test() : parcourt les combinaisons de plans et de donnees et
# signale les plantages et les incoherences internes du bilan.
#
# Execution, depuis la racine du package :
#   Rscript dev/test/agent_audit_mtest.R [fichier_sortie.csv]
# Colonnes du CSV : cas, statut (OK / ERREUR / INCOHERENT), detail.

racine <- getwd()
options(warn = -1)
pdf(NULL)
suppressMessages(suppressWarnings(capture.output(source("dev/load_all_kefir.R", encoding = "UTF-8"))))
setwd(racine)
suppressMessages({library(emmeans); library(multcomp)})
sortie <- commandArgs(TRUE)[1]

resultats <- list()
noter <- function(cas, statut, detail = "") {
  resultats[[length(resultats) + 1]] <<- data.frame(cas = cas, statut = statut, detail = substr(detail, 1, 200))
  if (!is.na(sortie)) write.csv(do.call(rbind, resultats), sortie, row.names = FALSE)
}

# Controles d'un resultat renvoye par m.test()
inspecter <- function(cas, r, o) {
  pb <- character(0)
  if (is.null(r)) pb <- c(pb, "resultat NULL")
  if (is.list(r)) {
    if (is.null(r$groups) || nrow(r$groups) == 0) pb <- c(pb, "$groups vide")
    if (is.null(r$method)) pb <- c(pb, "$method absent")
    p <- r$global_pvalue
    if (is.null(p) || length(p) != 1 || is.na(p) || p < 0 || p > 1) pb <- c(pb, paste("$global_pvalue anormale:", paste(format(p), collapse = ",")))
    if (!is.null(r$groups) && any(grepl("..", colnames(r$groups), fixed = TRUE))) pb <- c(pb, "nom de colonne altere")
    lettres <- unlist(r$groups[, -1, drop = FALSE])
    if (length(lettres) && all(lettres %in% c("?", NA))) pb <- c(pb, "lettres inconnues")
  }
  # Traces de debogage, ou erreurs internes signalees dans le bilan
  if (any(grepl("^\\[1\\] \"|av bootstrap|catego de base", o))) pb <- c(pb, "print de debogage")
  err <- grep("^.{0,6}(Erreur|Error)[ :]|chou[eé]|introuvable", o, value = TRUE)
  if (length(err)) pb <- c(pb, paste("dans le bilan:", err[1]))
  # Conclusion et test global doivent concorder (hors p tres proche du seuil)
  # (Non applicable si plusieurs facteurs : le bilan conclut sur chacun d'eux)
  if (is.list(r) && is.null(r$groups_by_factor) && is.null(r$effects) &&
      !is.null(r$global_pvalue) && !is.na(r$global_pvalue[1])) {
    p <- r$global_pvalue[1]
    if (p < 0.01 && any(grepl("Aucune diff.rence (de moyennes )?significative|Aucun effet significatif", o))) pb <- c(pb, paste("p =", signif(p, 3), "mais conclusion negative"))
    if (p > 0.2 && any(grepl("Diff.rences (de moyennes )?significatives entre", o))) pb <- c(pb, paste("p =", signif(p, 3), "mais conclusion positive"))
  }
  if (length(pb)) noter(cas, "INCOHERENT", paste(pb, collapse = " ; ")) else noter(cas, "OK")
}

essai <- function(cas, appel) {
  e <- substitute(appel)
  env <- parent.frame()
  r <- NULL
  o <- tryCatch(capture.output(r <- suppressMessages(suppressWarnings(eval(e, env)))),
                error = function(err) { noter(cas, "ERREUR", conditionMessage(err)); NULL })
  if (!is.null(o)) inspecter(cas, r, o)
  cat(".")
}

# Generateurs de reponses : normale, heteroscedastique, asymetrique, discrete, outliers
genere <- function(type, g, effet) {
  n <- length(g); m <- effet * (as.integer(factor(g)) - 1)
  switch(type,
    normal = rnorm(n, 10 + m, 2),
    hetero = rnorm(n, 10 + m, 1 + 2 * (as.integer(factor(g)) - 1)),
    asym = rexp(n, 1 / (2 + m)),
    discret = rpois(n, 3 + m),
    outlier = { y <- rnorm(n, 10 + m, 2); y[seq(1, n, by = 12)] <- 60; y })
}

set.seed(2026)
for (type in c("normal", "hetero", "asym", "discret", "outlier")) {
  for (effet in c(0, 3)) {
    for (ng in c(2, 3, 5)) for (ni in c(12, 40)) {
      g <- rep(letters[1:ng], each = ni)
      d <- data.frame(y = genere(type, g, effet), g = g, stringsAsFactors = FALSE)
      lab <- paste(type, "effet", effet, ng, "groupes n", ni)
      essai(paste("1F", lab), m.test(y ~ g, data = d, plot = FALSE, boot = FALSE, verbose = TRUE))
    }
    # Deux facteurs, equilibre et desequilibre
    d2 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), rep = 1:12, stringsAsFactors = FALSE)
    d2$y <- genere(type, d2$A, effet)
    essai(paste("2F equil", type, effet), m.test(y ~ A * B, data = d2, plot = FALSE, boot = FALSE, verbose = TRUE))
    essai(paste("2F additif", type, effet), m.test(y ~ A + B, data = d2, plot = FALSE, boot = FALSE, verbose = TRUE))
    d3 <- d2[-c(1:5, 20:21), ]
    essai(paste("2F desequil", type, effet), m.test(y ~ A * B, data = d3, plot = FALSE, boot = FALSE, verbose = TRUE))
    # ANCOVA
    d4 <- data.frame(F = rep(c("x", "y", "z"), each = 20), cov = runif(60, 50, 60))
    d4$y <- genere(type, d4$F, effet) + 0.3 * d4$cov
    essai(paste("ANCOVA +", type, effet), m.test(y ~ F + cov, data = d4, plot = FALSE, verbose = TRUE))
    essai(paste("ANCOVA *", type, effet), m.test(y ~ F * cov, data = d4, plot = FALSE, verbose = TRUE))
    d4$G <- rep(c("u", "v"), 30)
    essai(paste("ANCOVA 2F", type, effet), m.test(y ~ F + G + cov, data = d4, plot = FALSE, verbose = TRUE))
    # Apparie, 2 et 3 conditions
    d5 <- data.frame(id = rep(paste0("s", 1:20), times = 2), t = rep(c("t1", "t2"), each = 20))
    d5$y <- genere(type, d5$t, effet) + rep(rnorm(20), times = 2)
    essai(paste("apparie 2", type, effet), m.test(y ~ t, data = d5, paired = TRUE, id = "id", plot = FALSE, boot = FALSE, verbose = TRUE))
    d6 <- data.frame(id = rep(paste0("s", 1:20), times = 3), t = rep(c("t1", "t2", "t3"), each = 20))
    d6$y <- genere(type, d6$t, effet) + rep(rnorm(20), times = 3)
    essai(paste("apparie 3", type, effet), m.test(y ~ t, data = d6, paired = TRUE, id = "id", plot = FALSE, boot = FALSE, verbose = TRUE))
    # MANOVA
    d7 <- data.frame(g = rep(c("a", "b", "c"), each = 25)); d7$y1 <- genere(type, d7$g, effet); d7$y2 <- d7$y1 + rnorm(75, 0, 3)
    essai(paste("MANOVA", type, effet), m.test(cbind(y1, y2) ~ g, data = d7, plot = FALSE, verbose = TRUE))
    # Sorties alternatives
    p <- tryCatch(suppressMessages(suppressWarnings(m.test(y ~ g, data = d, return = FALSE))), error = function(e) e)
    if (inherits(p, "error") || !is.numeric(p) || length(p) != 1 || is.na(p)) noter(paste("return=FALSE 1F", type, effet), "INCOHERENT", paste(format(p), collapse = " "))
    p <- tryCatch(suppressMessages(suppressWarnings(m.test(y ~ A * B, data = d2, return = FALSE))), error = function(e) e)
    if (inherits(p, "error") || !is.numeric(p) || length(p) != 1 || is.na(p)) noter(paste("return=FALSE 2F", type, effet), "INCOHERENT", paste(format(p), collapse = " "))
    essai(paste("code=TRUE 1F", type, effet), m.test(y ~ g, data = d, plot = FALSE, boot = FALSE, code = TRUE))
    essai(paste("code=TRUE ANCOVA", type, effet), m.test(y ~ F + cov, data = d4, plot = FALSE, code = TRUE))
    essai(paste("plot 1F", type, effet), m.test(y ~ g, data = d, plot = TRUE, boot = FALSE))
    essai(paste("plot ANCOVA", type, effet), m.test(y ~ F + cov, data = d4, plot = TRUE))
  }
}

# Cas limites
set.seed(7)
n <- 60
b <- data.frame(y = rnorm(n, 10, 2), g = rep(c("a", "b", "c"), 20), f = factor(rep(c("a", "b", "c"), 20)), num = rep(1:3, 20), c2 = runif(n))
na <- b; na$y[c(3, 9)] <- NA; na$g[5] <- NA
essai("NA reponse et groupe", m.test(y ~ g, data = na, plot = FALSE, boot = FALSE))
essai("x,g NA", m.test(na$y, na$g, plot = FALSE, boot = FALSE))
essai("facteur type factor", m.test(y ~ f, data = b, plot = FALSE, boot = FALSE))
essai("groupe numerique", m.test(y ~ num, data = b, plot = FALSE, boot = FALSE))
essai("factor(num)", m.test(y ~ factor(num), data = b, plot = FALSE, boot = FALSE))
essai("data$", m.test(b$y ~ b$g, plot = FALSE, boot = FALSE))
essai("control", m.test(y ~ g, data = b, control = "a", plot = FALSE, boot = FALSE))
essai("control ANCOVA", m.test(y ~ g + c2, data = b, control = "a", plot = FALSE))
essai("variable absente", m.test(y ~ zzz, data = b, plot = FALSE))
z0 <- b; z0$y[z0$g == "a"] <- 0
essai("groupe de variance nulle", m.test(y ~ g, data = z0, plot = FALSE, boot = FALSE))
essai("tres petit groupe", m.test(c(1, 2, 3, 4, 5, 9), c("a", "a", "a", "b", "b", "b"), plot = FALSE, boot = FALSE))
essai("groupe singleton", m.test(c(1, 2, 3, 4, 5, 9, 7), c("a", "a", "a", "b", "b", "b", "c"), plot = FALSE, boot = FALSE))
essai("mixte Error", m.test(y ~ g + Error(id), data = transform(b, id = rep(paste0("s", 1:20), each = 3)), plot = FALSE))
essai("transforme 1F", m.test(log(y) ~ g, data = b, plot = FALSE, boot = FALSE))
essai("transforme 2F", m.test(sqrt(y) ~ g * f, data = transform(b, f = rep(c("u", "v"), 30)), plot = FALSE, boot = FALSE))
essai("anglais data", m.test(mpg ~ factor(cyl) + wt, data = mtcars, plot = FALSE))

tot <- do.call(rbind, resultats)
cat("\n\n"); print(table(tot$statut))
print(tot[tot$statut != "OK", ], right = FALSE)
