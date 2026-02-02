# Script d'audit pour Priorité 9
# Identifie les fonctions utilisées par m.test() et ses sous-fonctions

cat("\n=== AUDIT PACKAGE KefiR - PRIORITÉ 9 ===\n\n")

# Note: This audit script is for development only
# It should be run manually from the R/ directory, not loaded with the package

# Liste tous les fichiers .R
all_files <- list.files(pattern = "\\.R$", all.files = TRUE)
hidden_files <- grep("^\\.", all_files, value = TRUE)
public_files <- grep("^[^\\.]", all_files, value = TRUE)

cat("=== INVENTAIRE FICHIERS ===\n")
cat(sprintf("Fichiers cachés (.*.R):  %d\n", length(hidden_files)))
cat(sprintf("Fichiers publics (*.R):  %d\n", length(public_files)))
cat(sprintf("Total fichiers .R:       %d\n\n", length(all_files)))

# Fichiers chargés par load_all_kefir.R
fichiers_kefir <- c(
  ".msg.R",
  ".dbg.R",
  ".vbse.R",
  ".format_pval.R",
  ".drop_error_term.R",
  ".strip_data_dollar_safe.R",
  ".normalize_formula_dollar.R",
  ".normalize_from_formula.R",
  ".formulator.R",
  ".formulator_safe.R",
  ".detect_model_type.R",
  ".exit.R",
  ".align_pairs.R",
  ".auto_preprocess_g.R",
  ".control_independence.R",
  ".normality.R",
  ".variance.R",
  ".boots.R",
  "00_jb.norm.test.R",
  "discret.test.R",
  "valreg.R",
  ".one_factor_analysis.R",
  ".manova_analysis.R",
  ".multi_factor_analysis02.R",
  ".mixed_model_analysis.R",
  ".posthoc.R",
  ".posthoc_MANOVA.R",
  ".posthoc_ANCOVA.R",
  ".plot_with_letters.R",
  "m.test_temp18.R"
)

cat("=== FICHIERS CHARGÉS PAR load_all_kefir.R ===\n")
cat(sprintf("Nombre: %d\n\n", length(fichiers_kefir)))
for (f in fichiers_kefir) {
  cat(sprintf("  ✓ %s\n", f))
}

# Identifier fichiers NON chargés
cat("\n=== FICHIERS NON CHARGÉS (potentiellement obsolètes) ===\n\n")

fichiers_non_charges <- setdiff(all_files, c(fichiers_kefir, "load_all_kefir.R", "audit_package.R"))

# Catégoriser
tests <- grep("^test_", fichiers_non_charges, value = TRUE)
testeurs <- grep("testeur", fichiers_non_charges, value = TRUE)
autres <- setdiff(fichiers_non_charges, c(tests, testeurs))

cat(sprintf("Tests (test_*.R): %d fichiers\n", length(tests)))
for (f in tests) {
  cat(sprintf("  - %s\n", f))
}

cat(sprintf("\nTesteurs (.testeur*.R): %d fichiers\n", length(testeurs)))
for (f in testeurs) {
  cat(sprintf("  - %s\n", f))
}

cat(sprintf("\nAutres fonctions: %d fichiers\n", length(autres)))
for (f in autres) {
  cat(sprintf("  - %s\n", f))
}

# Logs et fichiers texte
logs <- list.files(pattern = "\\.(log|txt|md)$")
cat(sprintf("\n=== LOGS ET FICHIERS TEXTE ===\n"))
cat(sprintf("Nombre: %d\n\n", length(logs)))
for (l in logs) {
  cat(sprintf("  - %s\n", l))
}

# Recommandations
cat("\n=== RECOMMANDATIONS ===\n\n")

cat("À DÉPLACER DANS R/poubelle/ :\n\n")
cat("1. Fichiers de test (sauf test_code_mode.R et test_valreg_improvements.R) :\n")
for (f in setdiff(tests, c("test_code_mode.R", "test_valreg_improvements.R"))) {
  cat(sprintf("   - %s\n", f))
}

cat("\n2. Anciens testeurs :\n")
for (f in testeurs) {
  cat(sprintf("   - %s\n", f))
}

cat("\n3. Anciens logs (garder kefir.log et bp.log) :\n")
logs_a_deplacer <- setdiff(logs, c("kefir.log", "bp.log"))
for (l in logs_a_deplacer) {
  cat(sprintf("   - %s\n", l))
}

cat("\n4. Fonctions potentiellement obsolètes (à vérifier manuellement) :\n")
fonctions_a_verifier <- c("m.test_temp17.R", "load_all_deps.R", "run_testeur_complet.R",
                          "run_tests.R", "run_tests_simple.R")
for (f in intersect(fonctions_a_verifier, autres)) {
  cat(sprintf("   - %s\n", f))
}

cat("\n5. Autres fonctions non chargées (potentiellement utiles hors m.test) :\n")
autres_restants <- setdiff(autres, fonctions_a_verifier)
for (f in autres_restants) {
  cat(sprintf("   - %s\n", f))
}

# Statistiques finales
cat("\n=== STATISTIQUES FINALES ===\n")
total_a_deplacer <- length(setdiff(tests, c("test_code_mode.R", "test_valreg_improvements.R"))) +
                    length(testeurs) +
                    length(logs_a_deplacer) +
                    length(intersect(fonctions_a_verifier, autres))
cat(sprintf("Fichiers à déplacer (recommandé): %d\n", total_a_deplacer))
cat(sprintf("Fichiers actifs (load_all_kefir): %d\n", length(fichiers_kefir)))
cat(sprintf("Autres fonctions à analyser: %d\n", length(autres_restants)))

cat("\n=== FIN AUDIT ===\n")
