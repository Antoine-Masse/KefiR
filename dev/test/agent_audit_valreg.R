# Audit de valreg() : le verdict est compare aux tests pris un a un.
# Execution, depuis la racine du package :
#   Rscript dev/test/agent_audit_valreg.R

racine <- getwd()
options(warn = -1)
pdf(NULL)
suppressMessages(suppressWarnings(capture.output(source("dev/load_all_kefir.R", encoding = "UTF-8"))))
setwd(racine)

verdict <- function(reg, ...) {
  o <- character(0)
  r <- tryCatch({ o <- capture.output(v <- suppressMessages(valreg(reg, boot = FALSE, ...))); v }, error = function(e) paste("ERREUR:", conditionMessage(e)))
  list(v = r, alerte = paste(trimws(o[grepl("Attention|Warning|Plus de|Pas assez", o)]), collapse = " | "), o = o)
}
ligne <- function(nom, reg, attendu, ...) {
  x <- verdict(reg, ...)
  cat(if (identical(x$v, attendu)) "OK    " else "ECART ", nom, " : valreg = ", format(x$v), ", attendu = ", attendu,
      if (nzchar(x$alerte)) paste0("\n        alerte : ", substr(x$alerte, 1, 150)), "\n", sep = "")
}

set.seed(11); n <- 120
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), f7 = factor(rep(letters[1:8], 15)), f3 = factor(rep(c("u", "v", "w"), 40)))
d$y <- 2 + 1.5 * d$x1 - d$x2 + rnorm(n, 0, 0.8)
d$y7 <- d$y + as.integer(d$f7) * 0.8

ligne("modele sain", lm(y ~ x1 + x2, d), TRUE)
ligne("modele sain, facteur 8 modalites : pas de fausse colinearite", lm(y7 ~ x1 + f7, d), TRUE)
ligne("facteur sans effet : terme non significatif", lm(y ~ x1 + x2 + f3, d), FALSE)
d$yc <- d$y - 2 - mean(d$y - 2)
ligne("constante non significative : pas un motif de rejet", lm(yc ~ x1 + x2, d), TRUE)
d$x3 <- d$x1 + rnorm(n, 0, 0.01)
ligne("colinearite vraie", lm(y ~ x1 + x3 + x2, d), FALSE)
d$yb <- rnorm(n)
ligne("modele non significatif", lm(yb ~ x1, d), FALSE)
d$yh <- d$x1 + rnorm(n, 0, abs(d$x1) * 2 + 0.1)
ligne("heteroscedasticite", lm(yh ~ x1, d), FALSE)
d$yq <- d$x1^2 * 3 + rnorm(n, 0, 0.3)
ligne("non-linearite", lm(yq ~ x1, d), FALSE)
d$ye <- d$x1 + rexp(n, 0.3)^2
ligne("residus non normaux", lm(ye ~ x1, d), FALSE)
ligne("limite explicite depassee", lm(y ~ x1 + x2 + f3 + f7, d), FALSE, nvar = 5)
ligne("trop peu de valeurs", lm(y ~ x1, d[1:4, ]), FALSE)
petit <- d[1:18, ]
x <- verdict(lm(y ~ x1 + x2 + f3, petit)); cat("info  petit echantillon, 5 coefficients :", format(x$v), substr(x$alerte, 1, 120), "\n")
cat("\n--- modele sans predicteur\n"); x <- verdict(lm(y ~ 1, d)); cat(format(x$v), "\n")
cat("--- interaction\n"); x <- verdict(lm(y ~ x1 * f3, d)); cat(format(x$v), substr(x$alerte, 1, 200), "\n")
cat("--- k imbrique\n"); x <- tryCatch(capture.output(v <- valreg(lm(y ~ x1, d), k = 3, boot = FALSE)), error = function(e) conditionMessage(e)); print(v)
cat("--- bootstrap\n"); x <- capture.output(v <- valreg(lm(y ~ x1 + x2, d), boot = TRUE)); print(v)
cat("--- mixte\n")
suppressMessages(library(lme4)); dm <- transform(d, s = factor(rep(1:20, 6))); dm$ym <- dm$y + as.integer(dm$s) * 0.3
x <- tryCatch(capture.output(v <- suppressMessages(valreg(lmer(ym ~ x1 + (1 | s), dm), boot = FALSE))), error = function(e) { v <<- paste("ERREUR", conditionMessage(e)); character(0) }); print(v)
cat("\n--- corrigraph(mu=TRUE) (todolist)\n")
source("R/corrigraph.R", encoding = "UTF-8")
dc <- data.frame(q = c(rnorm(30), rep(NA, 30)), cat = rep(c("a", "b"), each = 30), q2 = rnorm(60), q3 = rnorm(60))
r <- tryCatch({ capture.output(corrigraph(dc, mu = TRUE, prop = TRUE)); "sans erreur" }, error = function(e) paste("ERREUR:", conditionMessage(e))); print(r)
