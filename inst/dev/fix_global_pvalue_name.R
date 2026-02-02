# Script pour corriger le bug du nom "c" dans global_pvalue

cat("=== CORRECTION DU BUG global_pvalue ===\n\n")

# Fichier .one_factor_analysis.R
file_path <- "R/.one_factor_analysis.R"
content <- readLines(file_path, warn = FALSE)

# Trouver la ligne avec le return problématique
idx <- grep("return.*check.*<-.*list.*global_pvalue.*pvals", content)

if (length(idx) > 0) {
  cat("Ligne trouvée:", idx, "\n")
  cat("Ancienne:", content[idx], "\n\n")

  # Nouvelle version avec noms explicites
  content[idx] <- "\t  return(check <- list(x = x, g = g, check_normality = check_normality, check_variance_equal = check_variance_equal, k = k, global_pvalue = pvals))"

  cat("Nouvelle:", content[idx], "\n")

  writeLines(content, file_path)
  cat("\n✓ .one_factor_analysis.R corrigé\n\n")
} else {
  cat("⚠ Ligne non trouvée dans .one_factor_analysis.R\n\n")
}

# Aussi corriger m.test.R si pas déjà fait
cat("=== VÉRIFICATION m.test.R ===\n")
file_path2 <- "R/m.test.R"
content2 <- readLines(file_path2, warn = FALSE)

idx2 <- grep("^[[:space:]]*return\\(global_pvalue\\)$", content2)

if (length(idx2) > 0) {
  cat("Ligne trouvée:", idx2, "\n")
  cat("Ancienne:", content2[idx2], "\n")

  content2[idx2] <- "    return(as.numeric(global_pvalue))"

  cat("Nouvelle:", content2[idx2], "\n")

  writeLines(content2, file_path2)
  cat("\n✓ m.test.R corrigé\n\n")
} else {
  cat("✓ m.test.R déjà corrigé (ou ligne non trouvée)\n\n")
}

cat("=== FIN DES CORRECTIONS ===\n")
