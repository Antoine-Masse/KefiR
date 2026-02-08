# Script de diagnostic pour localiser l'erreur NA exacte
# Exécuter APRÈS source("dev/load_all_kefir.R")

# Activer le traceback détaillé
options(error = function() {
  calls <- sys.calls()
  cat("\n========== TRACEBACK DÉTAILLÉ ==========\n")
  for (i in seq_along(calls)) {
    cat(sprintf("[%d] %s\n", i, paste(deparse(calls[[i]]), collapse = "\n     ")))
  }
  cat("=========================================\n")
})

# Test avec les vraies données
cat("\n=== Exécution du test ===\n")
tryCatch({
  result <- m.test(imc ~ Sport * Conso_alcool * Sexe, data = data, verbose = TRUE)
  cat("Succès!\n")
}, error = function(e) {
  cat("\nERREUR CAPTURÉE:", e$message, "\n")
})
