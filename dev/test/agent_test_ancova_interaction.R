# Test de non-regression : ANCOVA avec interaction facteur x covariable dans m.test()
#
# Bug corrige (2026-11-05) : m.test(y ~ F * cov) annoncait "interaction supprimee"
# mais gardait l'interaction dans le modele, et testait le facteur a cov = 0.
#
# Execution, depuis la racine du package (source le code de R/ par dev/load_all_kefir.R) :
#   Rscript dev/test/agent_test_ancova_interaction.R
# Chaque controle affiche OK ou ECHEC ; le script s'arrete en erreur si un controle echoue.

racine <- getwd()
suppressMessages(suppressWarnings(capture.output(source("dev/load_all_kefir.R", encoding = "UTF-8"))))
setwd(racine)
suppressMessages(library(emmeans)); suppressMessages(library(multcomp))

echecs <- 0
controle <- function(nom, condition) {
  cat(if (isTRUE(condition)) "OK    " else "ECHEC ", nom, "\n", sep = "")
  if (!isTRUE(condition)) echecs <<- echecs + 1
}
bilan <- function(...) {
  sortie <- capture.output(r <- suppressWarnings(m.test(..., plot = FALSE, verbose = TRUE)))
  list(r = r, sortie = sortie)
}
proche <- function(a, b, tol = 1e-6) isTRUE(abs(a - b) < tol)
type3 <- function(formule, d) {
  ct <- list(F = "contr.sum")
  car::Anova(lm(formule, data = d, contrasts = ct), type = "III")
}

# Donnees : covariable loin de 0, deux groupes separes, plan equilibre
set.seed(2)  # graine retenue : voie parametrique (assomptions respectees)
n <- 40
d <- data.frame(F = rep(c("A", "B"), each = n / 2), x = runif(n, 200, 220))
d$y <- 10 + 2 * (d$F == "B") + 0.05 * d$x + rnorm(n, 0, 1)

# ---------------------------------------------------------------------------
cat("\n1) Interaction NON significative : elle quitte reellement le modele\n")
b <- bilan(y ~ F * x, data = d)
ref <- type3(y ~ F + x, d)["F", "Pr(>F)"]
controle("la formule ajustee n'a plus d'interaction", any(grepl("Formule : y ~ F \\+ x", b$sortie)))
controle("aucune ligne F:x dans le tableau", !any(grepl("^F:x", b$sortie)))
controle("pas de 'Pentes heterogenes autorisees'", !any(grepl("Pentes h.t.rog.nes autoris", b$sortie)))
controle("5 assomptions annoncees", any(grepl("ASSOMPTION 1/5", b$sortie)))
controle("p globale = Type III du modele additif", proche(b$r$global_pvalue, ref))
controle("post-hoc coherent avec le tableau", proche(b$r$p.value["A", "B"], ref))
controle("return=FALSE = meme p", proche(m.test(y ~ F * x, data = d, return = FALSE), ref))
controle("meme resultat que la formule additive saisie",
         proche(bilan(y ~ F + x, data = d)$r$global_pvalue, ref))

# ---------------------------------------------------------------------------
cat("\n2) Interaction SIGNIFICATIVE : gardee, covariable centree\n")
set.seed(1)
s <- data.frame(F = rep(c("A", "B", "C"), each = 30), x = runif(90, 100, 120))
s$y <- 50 + c(0, 4, 8)[factor(s$F)] + c(1, 3, -2)[factor(s$F)] * (s$x - 110) + rnorm(90, 0, 2)
b <- bilan(y ~ F * x, data = s)
sc <- s; sc$x <- s$x - mean(s$x)
tab <- type3(y ~ F * x, sc)
controle("le modele garde l'interaction", any(grepl("^F:x", b$sortie)))
controle("le centrage est annonce", any(grepl("centr.e\\(s\\) sur leur moyenne", b$sortie)))
controle("p du facteur = Type III centre", proche(b$r$global_pvalue, tab["F", "Pr(>F)"]))
controle("la conclusion cite l'interaction", any(grepl("Interactions : F:x", b$sortie)))
brut <- emmeans(lm(y ~ F * x, data = s), "F")
controle("moyennes ajustees identiques avec ou sans centrage",
         proche(sum(abs(sort(b$r$descriptive_stats$adjusted_mean) - sort(summary(brut)$emmean))), 0))
controle("note post-hoc : moyenne d'origine de la covariable",
         grepl(paste0("x = ", round(mean(s$x), 2)), b$r$details$note, fixed = TRUE))

# ---------------------------------------------------------------------------
cat("\n3) La p globale designe le facteur, pas la premiere ligne\n")
ref <- anova(lm(y ~ x + F, data = d))["F", "Pr(>F)"]
controle("covariable en premier, return=FALSE", proche(m.test(y ~ x + F, data = d, return = FALSE), ref))
set.seed(1)
t3 <- data.frame(F = rep(c("A", "B", "C"), each = 15), x = runif(45, 20, 30))
t3$y <- c(0, 1, 2)[factor(t3$F)] + 0.2 * t3$x + rnorm(45, 0, 1)
b <- bilan(y ~ x + F, data = t3)
controle("3 groupes : $global_pvalue = test F du facteur",
         proche(b$r$global_pvalue, anova(lm(y ~ x + F, data = t3))["F", "Pr(>F)"]))

# ---------------------------------------------------------------------------
cat("\n4) Cas annexes\n")
m <- mtcars; m$am <- factor(m$am)
b <- bilan(mpg ~ am + wt, data = m)
pmat <- b$r$p.value; diag(pmat) <- NA
controle("niveaux numeriques : matrice de p remplie", min(pmat, na.rm = TRUE) < 1)
b <- bilan(log(mpg) ~ am * wt, data = m)
controle("reponse transformee log(y) : pas d'erreur", !is.null(b$r))
des <- d[-(1:6), ]
b <- bilan(y ~ F + x, data = des)
controle("plan desequilibre : pas de formule suggeree", !any(grepl("Facultatif", b$sortie)))
b <- bilan(y ~ F + x, data = d)
controle("plan equilibre, covariable apres : Type I propose", any(grepl("Facultatif.*y ~ x \\+ F", b$sortie)))
controle("message anglais 'tukey was changed' supprime", !any(grepl("was changed to", b$sortie)))
sortie_code <- capture.output(suppressWarnings(m.test(y ~ F * x, data = s, plot = FALSE, code = TRUE)))
controle("code=TRUE reproduit le centrage", any(grepl("dt\\$x <- dt\\$x - mean", sortie_code)))

# ---------------------------------------------------------------------------
cat("\n5) Anomalies annexes corrigees le 2026-11-05\n")
# a) $method a 2 groupes : donnees a queue lourde, voie Wilcoxon
set.seed(5)
w <- data.frame(x = c(rexp(40, 1), rexp(40, 0.2)), g = rep(c("A", "B"), each = 40))
invisible(capture.output(r <- suppressWarnings(m.test(w$x, w$g, plot = FALSE, boot = FALSE))))
controle("$method a 2 groupes annonce Wilcoxon", grepl("Wilcoxon", r$method))
# b) dunnTest ne doit plus rien ecrire par message()
msg <- character(0)
invisible(withCallingHandlers(capture.output(pairwise(iris[, 1], iris[, 5], type = "dunn")),
                    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") }))
if (length(msg) > 0) cat("  messages recus :", substr(msg, 1, 80), sep = "\n  ")
controle("pairwise(type='dunn') silencieux", !any(grepl("Dunn|Kruskal", msg)))
# c) valreg indique comment depasser la limite
o <- capture.output(v <- valreg(lm(Sepal.Length ~ Species * Sepal.Width, data = iris), verbose = TRUE))
controle("valreg annonce nvar", isFALSE(v) && any(grepl("nvar = 6", o)))
# d) medpb2 : deux groupes non differents partagent une lettre
set.seed(6)
ech <- data.frame(x = c(rnorm(60, 0), rnorm(60, 20), rnorm(60, 0.01), rnorm(60, 40)),
                  g = rep(c("P", "Q", "R", "S"), each = 60))
l <- .pairwise_medpb2(ech$x, ech$g, nboot = 300)$groups
controle("medpb2 : P et R meme lettre, Q et S a part",
         l$medpb2_Holm[l$categories == "P"] == l$medpb2_Holm[l$categories == "R"] &&
           length(unique(l$medpb2_Holm)) == 3)
# f) noms de colonnes stables dans $groups
invisible(capture.output(r <- suppressWarnings(m.test(iris[, 1], iris[, 5], plot = FALSE, iter = 30))))
controle("colonnes non alterees par data.frame()", !any(grepl("..", colnames(r$groups), fixed = TRUE)))

cat("\n", if (echecs == 0) "TOUS LES CONTROLES PASSENT" else paste(echecs, "ECHEC(S)"), "\n", sep = "")
if (echecs > 0) quit(status = 1)
