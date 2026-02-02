# Audit dépendances - Priorité 9

cat("\n=== AUDIT DÉPENDANCES - PRIORITÉ 9 ===\n\n")

# Packages utilisés dans le code
packages_code <- c(
  # Core stats
  "stats", "graphics", "grDevices",
  
  # Mixed models
  "lme4", "lmerTest",
  
  # Tests et ANOVA
  "car", "agricolae", "emmeans",
  
  # MANOVA
  "MVN", "biotools", "rrcov", "MASS",
  
  # Robustes
  "WRS2",
  
  # Mesures répétées
  "ez",
  
  # Autres
  "lmtest", "tseries"
)

cat("Packages identifiés dans le code:", length(packages_code), "\n\n")
for (p in sort(packages_code)) {
  cat(sprintf("  - %s\n", p))
}

cat("\n=== VÉRIFICATION INSTALLATION ===\n\n")

installed <- sapply(packages_code, requireNamespace, quietly = TRUE)

cat(sprintf("Installés: %d/%d (%.1f%%)\n\n", 
            sum(installed), length(installed),
            100*sum(installed)/length(installed)))

if (sum(!installed) > 0) {
  cat("⚠ Packages manquants:\n")
  for (p in packages_code[!installed]) {
    cat(sprintf("  - %s\n", p))
  }
}

cat("\n=== RECOMMANDATIONS DESCRIPTION ===\n\n")

cat("Imports (requis):\n")
cat("  lme4, lmerTest, car, agricolae, emmeans, lmtest\n\n")

cat("Suggests (optionnels pour fonctionnalités avancées):\n")
cat("  MVN (MANOVA), biotools (MANOVA), rrcov (MANOVA robuste),\n")
cat("  WRS2 (analyses robustes), ez (mesures répétées),\n")
cat("  tseries (tests série temporelle)\n\n")

cat("Depends:\n")
cat("  R (>= 4.0.0)\n\n")

cat("=== FIN AUDIT ===\n")
