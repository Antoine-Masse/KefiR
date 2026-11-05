# Compare .bay_inclusion a bayestestR::bayesfactor_inclusion, et elucide l'ecart
# A:B apparie 0,4547 (a la main) contre 0,4415 (bayestestR) de agent_02.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress = FALSE)
suppressMessages(library(bayestestR))
pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE)

out <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats/agent_10_inclusion_vs_bayestestR.txt"
sink(out, split = TRUE)
on.exit(sink())
cat("BayesFactor", as.character(packageVersion("BayesFactor")),
    "; bayestestR", as.character(packageVersion("bayestestR")), "\n")

# Tableau de modeles a partir d'un objet BFBayesFactor, sans aucun recalcul :
# modele nul (log BF = 0), plus extractBF.
table_of <- function(ba) {
  x <- extractBF(ba, logbf = TRUE)
  list(log_bf = c(0, x$bf), labels = c("1", rownames(x)), error = c(0, x$error))
}
tab_bt <- function(ba, match) {
  r <- as.data.frame(bayesfactor_inclusion(ba, match_models = match))
  r$term <- rownames(r)
  r
}
compare <- function(ba, type, match) {
  t0 <- table_of(ba)
  a <- .bay_inclusion(t0$log_bf, t0$labels, type)
  b <- tab_bt(ba, match)
  b <- b[match(a$term, b$term), ]
  print(data.frame(term = a$term, log_bf_kefir = a$log_bf, log_bf_bayestestR = b$log_BF,
                   ecart = a$log_bf - b$log_BF, prior_kefir = a$prior_h1, prior_bt = b$p_prior,
                   p_h1_kefir = a$p_h1, p_post_bt = b$p_posterior), digits = 6, row.names = FALSE)
  invisible(a)
}

cat("\n=== 1. Source de bayestestR : la definition appariee ===\n")
cat("La regle est lue dans bayesfactor_inclusion.bayesfactor_models et .includes_interaction :\n")
print(getAnywhere(".includes_interaction")$objs[[1]])
print(getAnywhere("bayesfactor_inclusion.BFBayesFactor")$objs[[1]])
fn <- deparse(getAnywhere("bayesfactor_inclusion.bayesfactor_models")$objs[[1]])
cat(fn[grep("isTRUE\\(match_models\\)", fn)[2]:grep("mwith <- ", fn)], sep = "\n")
cat("\nextractBF lit des valeurs deja stockees ; le source de bayesfactor_models.BFBayesFactor\n")
cat("ne contient aucun appel d'echantillonnage (recompute est absent) :\n")
src <- deparse(getAnywhere("bayesfactor_models.BFBayesFactor")$objs[[1]])
cat("  appel a recompute dans le source :", any(grepl("recompute", src)), "\n")
cat("  extractBF dans le source :", any(grepl("extractBF", src)), "\n")

cat("\n=== 2. Le calcul de agent_02 : graine 221, whichModels = 'all' ===\n")
set.seed(221)
d <- expand.grid(A = gl(2, 1, labels = c("a1", "a2")), B = gl(3, 1, labels = c("b1", "b2", "b3")), r = 1:10)
d$y <- rnorm(nrow(d)) + .7 * (d$A == "a2") + .6 * (d$A == "a2" & d$B == "b3")
ba_all <- anovaBF(y ~ A * B, d, whichModels = "all")
x <- extractBF(ba_all); x$modele <- rownames(x)
print(x[, c("modele", "bf", "error")], row.names = FALSE)
b_get <- function(m) x$bf[x$modele == m]
before <- extractBF(ba_all)$bf
bt_all <- tab_bt(ba_all, TRUE)
after <- extractBF(ba_all)$bf
cat("Objet inchange apres l'appel a bayestestR :", identical(before, after), "\n")
bf_bt <- exp(bt_all["A:B", "log_BF"])
main <- b_get("A + B + A:B") / b_get("A + B")
viol <- mean(c(b_get("A:B"), b_get("A + A:B"), b_get("B + A:B"), b_get("A + B + A:B"))) / b_get("A + B")
cat(sprintf("  bayestestR apparie, A:B                                    : %.6f\n", bf_bt))
cat(sprintf("  A la main (A+B+A:B)/(A+B)                                  : %.6f\n", main))
cat(sprintf("  mean(A:B, A+A:B, B+A:B, A+B+A:B)/(A+B), MEME objet         : %.6f\n", viol))
cat(sprintf("  ecart bayestestR - formule 4 modeles                       : %.2e\n", bf_bt - viol))
cat(sprintf("  p_prior bayestestR de A:B = %.3f, soit 4/8 modeles dans M+ (non renormalise)\n",
            bt_all["A:B", "p_prior"]))
cat("  .bay_inclusion sur cet espace, 'matched' :",
    tryCatch({ .bay_inclusion(c(0, x$bf * 0 + log(x$bf)), c("1", x$modele), "matched"); "accepte" },
             error = function(e) paste("refus :", conditionMessage(e))), "\n")

cat("\n=== 3. Variabilite Monte Carlo de ce rapport, a titre de comparaison ===\n")
cat("Meme jeu de donnees, nouveaux anovaBF. L'ecart agent_02 a log = ",
    sprintf("%.4f", log(main / bf_bt)), "\n", sep = "")
set.seed(9001)
rep_main <- rep_viol <- numeric(8)
for (i in seq_len(8)) {
  xi <- extractBF(anovaBF(y ~ A * B, d, whichModels = "all"))
  gi <- function(m) xi$bf[rownames(xi) == m]
  rep_main[i] <- gi("A + B + A:B") / gi("A + B")
  rep_viol[i] <- mean(c(gi("A:B"), gi("A + A:B"), gi("B + A:B"), gi("A + B + A:B"))) / gi("A + B")
}
cat("  (A+B+A:B)/(A+B)   :", round(rep_main, 4), "\n")
cat("  formule 4 modeles :", round(rep_viol, 4), "\n")
cat(sprintf("  ET de log(rapport) : %.4f ; ecart moyen des deux definitions (log, meme objet) : %.4f\n",
            sd(log(rep_main)), mean(log(rep_main / rep_viol))))

cat("\n=== 4. Espace marginal : anovaBF whichModels = 'withmain', graine 221 ===\n")
set.seed(221)
ba <- anovaBF(y ~ A * B, d, whichModels = "withmain")
print(extractBF(ba)[, c("bf", "error")])
cat("\n-- match_models = TRUE contre type = 'matched' (meme objet)\n")
compare(ba, "matched", TRUE)
cat("\n-- match_models = FALSE contre type = 'all' (meme objet)\n")
compare(ba, "all", FALSE)

cat("\n=== 5. Second jeu, desequilibre, graine fixe ===\n")
set.seed(1010)
d2 <- data.frame(A = factor(sample(c("a1", "a2", "a3"), 90, TRUE, prob = c(.5, .3, .2))),
                 B = factor(sample(c("b1", "b2"), 90, TRUE)))
d2$y <- rnorm(90) + .5 * (d2$A == "a3") + .4 * (d2$B == "b2") + .5 * (d2$A == "a2" & d2$B == "b2")
ba2 <- anovaBF(y ~ A * B, d2, whichModels = "withmain")
print(extractBF(ba2)[, c("bf", "error")])
cat("\n-- match_models = TRUE contre type = 'matched'\n")
compare(ba2, "matched", TRUE)
cat("\n-- match_models = FALSE contre type = 'all'\n")
compare(ba2, "all", FALSE)
sink()
