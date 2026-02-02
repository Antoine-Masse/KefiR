# Audit Roxygen et en-têtes - Priorité 9
# Vérifie documentation des 30 fonctions actives

cat("\n=== AUDIT ROXYGEN & EN-TÊTES - PRIORITÉ 9 ===\n\n")

setwd("/mnt/c/Users/masse/Desktop/KefiR/KefiR/R")

# Fichiers actifs (load_all_kefir.R)
fichiers_actifs <- c(
  ".msg.R", ".dbg.R", ".vbse.R", ".format_pval.R", ".drop_error_term.R",
  ".strip_data_dollar_safe.R", ".normalize_formula_dollar.R",
  ".normalize_from_formula.R", ".formulator.R", ".formulator_safe.R",
  ".detect_model_type.R", ".exit.R", ".align_pairs.R", ".auto_preprocess_g.R",
  ".control_independence.R", ".normality.R", ".variance.R", ".boots.R",
  "00_jb.norm.test.R", "discret.test.R", "valreg.R", ".one_factor_analysis.R",
  ".manova_analysis.R", ".multi_factor_analysis02.R", ".mixed_model_analysis.R",
  ".posthoc.R", ".posthoc_MANOVA.R", ".posthoc_ANCOVA.R", ".plot_with_letters.R",
  "m.test_temp18.R"
)

cat("Fichiers à auditer:", length(fichiers_actifs), "\n\n")

results <- data.frame(
  fichier = character(),
  roxygen = character(),
  description = character(),
  params = character(),
  return_doc = character(),
  export = character(),
  stringsAsFactors = FALSE
)

for (f in fichiers_actifs) {
  if (!file.exists(f)) {
    results <- rbind(results, data.frame(
      fichier = f,
      roxygen = "MANQUANT",
      description = "N/A",
      params = "N/A",
      return_doc = "N/A",
      export = "N/A"
    ))
    next
  }

  content <- readLines(f, warn = FALSE, n = 100)

  # Vérifier présence #'
  has_roxygen <- any(grepl("^#'", content))

  # Vérifier @description ou description
  has_desc <- any(grepl("@description|^#' [A-Z]", content[1:30]))

  # Vérifier @param
  has_params <- any(grepl("@param", content))
  n_params <- sum(grepl("@param", content))

  # Vérifier @return
  has_return <- any(grepl("@return", content))

  # Vérifier @export (pour fonctions publiques)
  has_export <- any(grepl("@export", content))
  is_hidden <- grepl("^\\.", f)

  results <- rbind(results, data.frame(
    fichier = f,
    roxygen = ifelse(has_roxygen, "OUI", "NON"),
    description = ifelse(has_desc, "OUI", "NON"),
    params = ifelse(has_params, paste0("OUI (", n_params, ")"), "NON"),
    return_doc = ifelse(has_return, "OUI", "NON"),
    export = ifelse(is_hidden, "HIDDEN", ifelse(has_export, "OUI", "NON"))
  ))
}

cat("=== RÉSULTATS AUDIT ===\n\n")

# Fonctions SANS Roxygen
sans_roxygen <- results[results$roxygen == "NON", ]
if (nrow(sans_roxygen) > 0) {
  cat("⚠ FONCTIONS SANS DOCUMENTATION ROXYGEN:", nrow(sans_roxygen), "\n")
  for (i in 1:nrow(sans_roxygen)) {
    cat(sprintf("  - %s\n", sans_roxygen$fichier[i]))
  }
  cat("\n")
}

# Fonctions SANS description
sans_desc <- results[results$roxygen == "OUI" & results$description == "NON", ]
if (nrow(sans_desc) > 0) {
  cat("⚠ FONCTIONS SANS @description:", nrow(sans_desc), "\n")
  for (i in 1:nrow(sans_desc)) {
    cat(sprintf("  - %s\n", sans_desc$fichier[i]))
  }
  cat("\n")
}

# Fonctions SANS @param
sans_params <- results[results$roxygen == "OUI" & results$params == "NON", ]
if (nrow(sans_params) > 0) {
  cat("⚠ FONCTIONS SANS @param:", nrow(sans_params), "\n")
  for (i in 1:nrow(sans_params)) {
    cat(sprintf("  - %s\n", sans_params$fichier[i]))
  }
  cat("\n")
}

# Fonctions SANS @return
sans_return <- results[results$roxygen == "OUI" & results$return_doc == "NON", ]
if (nrow(sans_return) > 0) {
  cat("⚠ FONCTIONS SANS @return:", nrow(sans_return), "\n")
  for (i in 1:nrow(sans_return)) {
    cat(sprintf("  - %s\n", sans_return$fichier[i]))
  }
  cat("\n")
}

# Fonctions publiques SANS @export
publiques_sans_export <- results[results$export == "NON", ]
if (nrow(publiques_sans_export) > 0) {
  cat("⚠ FONCTIONS PUBLIQUES SANS @export:", nrow(publiques_sans_export), "\n")
  for (i in 1:nrow(publiques_sans_export)) {
    cat(sprintf("  - %s (à exporter ?)\n", publiques_sans_export$fichier[i]))
  }
  cat("\n")
}

# Statistiques globales
cat("=== STATISTIQUES GLOBALES ===\n\n")
cat(sprintf("Total fichiers audités:        %d\n", nrow(results)))
cat(sprintf("Avec Roxygen:                  %d (%.1f%%)\n",
            sum(results$roxygen == "OUI"),
            100*sum(results$roxygen == "OUI")/nrow(results)))
cat(sprintf("Avec @description:             %d (%.1f%%)\n",
            sum(results$description == "OUI"),
            100*sum(results$description == "OUI")/nrow(results)))
cat(sprintf("Avec @param:                   %d (%.1f%%)\n",
            sum(results$params != "NON"),
            100*sum(results$params != "NON")/nrow(results)))
cat(sprintf("Avec @return:                  %d (%.1f%%)\n",
            sum(results$return_doc == "OUI"),
            100*sum(results$return_doc == "OUI")/nrow(results)))
cat(sprintf("Fonctions cachées (hidden):    %d\n", sum(results$export == "HIDDEN")))
cat(sprintf("Fonctions publiques exportées: %d\n", sum(results$export == "OUI")))

cat("\n=== FIN AUDIT ===\n")
