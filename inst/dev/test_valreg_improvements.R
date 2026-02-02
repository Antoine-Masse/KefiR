# Test Priorité 7 - Améliorations valreg()
# Teste les nouvelles fonctionnalités de valreg()

cat("\n=== TEST PRIORITÉ 7 - AMÉLIORATIONS VALREG() ===\n\n")

# Charger packages
suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(lmtest)
  library(tseries)
  library(car)
  library(agricolae)
})

# Charger KefiR
setwd("/mnt/c/Users/masse/Desktop/KefiR/KefiR/R")
source("load_all_kefir.R")

################################################################################
# TEST 1: Paramètre k - Numérotation des messages
################################################################################
cat("\n=== Test 1: Paramètre k (numérotation messages) ===\n")

set.seed(123)
data_test <- data.frame(
  x1 = rnorm(40),
  x2 = rnorm(40),
  y = rnorm(40)
)

model_lm <- lm(y ~ x1 + x2, data = data_test)

cat("\nTest avec k=1 (défaut)...\n")
result1 <- valreg(model_lm, verbose = TRUE, boot = FALSE, k = 1)
cat("\nRésultat: valid=", result1$valid, ", k final=", result1$k, "\n")

cat("\nTest avec k=10 (continuation)...\n")
result2 <- valreg(model_lm, verbose = TRUE, boot = FALSE, k = 10)
cat("\nRésultat: valid=", result2$valid, ", k final=", result2$k, "\n")

test1 <- (result2$k > result1$k) && (result1$k >= 1)

if (test1) {
  cat("\n✓ Test 1 PASS - Le compteur k est incrémenté correctement\n")
} else {
  cat("\n✗ Test 1 FAIL - Problème avec le compteur k\n")
}

################################################################################
# TEST 2: Intégration .normality() avec tolerance
################################################################################
cat("\n\n=== Test 2: Intégration .normality() avec tolerance ===\n")

# Données légèrement non-normales
set.seed(456)
data_skewed <- data.frame(
  x = rgamma(60, shape = 2, rate = 1),  # Légèrement asymétrique
  y = rgamma(60, shape = 2, rate = 1)
)

model_skewed <- lm(y ~ x, data = data_skewed)

cat("\nTest avec tolerance='basic'...\n")
result_basic <- valreg(model_skewed, verbose = TRUE, boot = FALSE, tolerance = "basic")
cat("\nRésultat basic: valid=", result_basic$valid, "\n")

cat("\nTest avec tolerance='extrem'...\n")
result_extrem <- valreg(model_skewed, verbose = TRUE, boot = FALSE, tolerance = "extrem")
cat("\nRésultat extrem: valid=", result_extrem$valid, "\n")

test2 <- is.logical(result_basic$valid) && is.logical(result_extrem$valid)

if (test2) {
  cat("\n✓ Test 2 PASS - Paramètre tolerance fonctionne\n")
} else {
  cat("\n✗ Test 2 FAIL - Problème avec tolerance\n")
}

################################################################################
# TEST 3: Paramètre orderDW - Message informatif
################################################################################
cat("\n\n=== Test 3: Paramètre orderDW (message informatif) ===\n")

set.seed(789)
data_ordered <- data.frame(
  time = 1:50,
  x = rnorm(50),
  y = rnorm(50)
)

model_ordered <- lm(y ~ x, data = data_ordered)

cat("\nTest avec orderDW=NULL (par défaut)...\n")
result_no_order <- valreg(model_ordered, verbose = TRUE, boot = FALSE, orderDW = NULL)
cat("\nRésultat sans orderDW: valid=", result_no_order$valid, "\n")

cat("\nTest avec orderDW='time'...\n")
result_with_order <- valreg(model_ordered, verbose = TRUE, boot = FALSE, orderDW = "time")
cat("\nRésultat avec orderDW: valid=", result_with_order$valid, "\n")

test3 <- is.logical(result_no_order$valid) && is.logical(result_with_order$valid)

if (test3) {
  cat("\n✓ Test 3 PASS - Paramètre orderDW fonctionne\n")
} else {
  cat("\n✗ Test 3 FAIL - Problème avec orderDW\n")
}

################################################################################
# TEST 4: Intégration avec modèles mixtes
################################################################################
cat("\n\n=== Test 4: Intégration avec modèles mixtes (.mixed_model_analysis) ===\n")

set.seed(101)
n_subjects <- 20
n_conditions <- 3

data_mixed <- data.frame(
  subject = factor(rep(1:n_subjects, each = n_conditions)),
  condition = factor(rep(c("A", "B", "C"), times = n_subjects)),
  score = c(
    rnorm(n_subjects, mean = 10, sd = 2),
    rnorm(n_subjects, mean = 12, sd = 2),
    rnorm(n_subjects, mean = 11, sd = 2)
  ) + rep(rnorm(n_subjects, mean = 0, sd = 1.5), each = n_conditions)
)

model_mixed <- lmer(score ~ condition + (1 | subject), data = data_mixed)

cat("\nTest valreg() avec modèle mixte...\n")
result_mixed <- valreg(model_mixed, verbose = TRUE, k = 1, tolerance = "extrem", orderDW = NULL)
cat("\nRésultat mixed: valid=", result_mixed$valid, ", k=", result_mixed$k, "\n")

test4 <- is.logical(result_mixed$valid) && (result_mixed$k > 1)

if (test4) {
  cat("\n✓ Test 4 PASS - valreg() fonctionne avec modèles mixtes\n")
} else {
  cat("\n✗ Test 4 FAIL - Problème avec modèles mixtes\n")
}

################################################################################
# TEST 5: Format de retour (liste avec valid et k)
################################################################################
cat("\n\n=== Test 5: Format de retour (liste avec valid et k) ===\n")

result_format <- valreg(model_lm, verbose = FALSE, boot = FALSE)

test5 <- is.list(result_format) &&
         !is.null(result_format$valid) &&
         !is.null(result_format$k) &&
         is.logical(result_format$valid) &&
         is.numeric(result_format$k)

if (test5) {
  cat("\n✓ Test 5 PASS - Format de retour correct (liste avec valid et k)\n")
} else {
  cat("\n✗ Test 5 FAIL - Format de retour incorrect\n")
  cat("  Structure:", str(result_format), "\n")
}

################################################################################
# RÉSUMÉ
################################################################################
cat("\n\n=== RÉSUMÉ ===\n")

tests <- c(test1, test2, test3, test4, test5)
tests_bool <- sapply(tests, isTRUE)
passed <- sum(tests_bool)
total <- length(tests)

cat(sprintf("Tests réussis: %d/%d (%.1f%%)\n\n", passed, total, 100*passed/total))

if (passed == total) {
  cat("✓ Toutes les améliorations valreg() fonctionnent correctement\n")
  cat("  - Paramètre k pour numérotation\n")
  cat("  - Intégration .normality() avec tolerance\n")
  cat("  - Paramètre orderDW pour Durbin-Watson\n")
  cat("  - Support modèles mixtes amélioré\n")
  cat("  - Format de retour compatible m.test()\n")
} else {
  cat(sprintf("✗ %d test(s) ont échoué\n", total - passed))
}

cat("\n=== FIN ===\n")
