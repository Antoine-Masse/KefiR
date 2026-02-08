# Script pour tracer l'erreur exacte dans .posthoc()
# Exécuter avec: source("dev/debug_posthoc.R")

# Charger les fichiers sources
for (f in list.files('R', pattern='\\.R$', full.names=TRUE)) {
  tryCatch(source(f), error=function(e) message('Skip: ', f))
}

set.seed(123)
n <- 200
data <- data.frame(
  imc = rnorm(n, 25, 5),
  Sport = factor(sample(c('Oui', 'Non', 'Parfois'), n, replace=TRUE)),
  Conso_alcool = factor(sample(c('Jamais', 'Occasionnel', 'Regulier'), n, replace=TRUE)),
  Sexe = factor(sample(c('H', 'F'), n, replace=TRUE))
)

# Appeler avec traceback activé
options(error = function() {
  cat("\n=== TRACEBACK ===\n")
  traceback(2)
})

result <- m.test(imc ~ Sport * Conso_alcool * Sexe, data = data, verbose = TRUE, debug = TRUE)
