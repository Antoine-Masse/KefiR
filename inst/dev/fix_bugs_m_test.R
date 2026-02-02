# Script pour corriger les bugs identifiés

# BUG 1: Dans m.test.R, return=FALSE retourne p-value avec nom "c"
# SOLUTION: Utiliser as.numeric() ou unname() pour enlever le nom

# BUG 2: .testeur_m.test.R utilise plot=FALSE, verbose=FALSE partout
# SOLUTION: Modifier les tests pour utiliser plot=TRUE, verbose=TRUE par défaut

cat("=== CORRECTIONS DES BUGS ===\n\n")

# ========== BUG 1: m.test.R - Nom "c" pour p-value ==========

cat("1. Correction du bug return=FALSE dans m.test.R...\n")

file_path <- "R/m.test.R"
content <- readLines(file_path, warn = FALSE)

# Trouver la ligne "return(global_pvalue)"
idx <- which(grepl("^[[:space:]]*return\\(global_pvalue\\)", content))

if (length(idx) > 0) {
  cat("   Ligne trouvée:", idx, "\n")
  cat("   Ancienne ligne:", content[idx], "\n")

  # Remplacer par return(as.numeric(global_pvalue)) pour enlever les noms
  content[idx] <- gsub("return\\(global_pvalue\\)",
                       "return(as.numeric(global_pvalue))",
                       content[idx])

  cat("   Nouvelle ligne:", content[idx], "\n")

  writeLines(content, file_path)
  cat("   ✓ m.test.R corrigé\n\n")
} else {
  cat("   ⚠ Ligne 'return(global_pvalue)' non trouvée\n\n")
}

# ========== BUG 2: .testeur_m.test.R - Ajouter plot=TRUE, verbose=TRUE ==========

cat("2. Correction des paramètres de test dans .testeur_m.test.R...\n")

file_path2 <- "R/.testeur_m.test.R"
content2 <- readLines(file_path2, warn = FALSE)

# Compter les occurrences AVANT
n_before_verbose_false <- length(grep("verbose = FALSE", content2, fixed = TRUE))
n_before_plot_false <- length(grep("plot = FALSE", content2, fixed = TRUE))

cat("   Avant: verbose=FALSE:", n_before_verbose_false, "| plot=FALSE:", n_before_plot_false, "\n")

# Stratégie: Remplacer SEULEMENT dans les tests qui ne testent PAS
# spécifiquement return=FALSE ou code=TRUE

# Identifier les lignes à NE PAS modifier (tests de la catégorie "return" et mode "code")
# Ces lignes doivent garder verbose=FALSE et plot=FALSE

# Pour tous les autres tests (formules, anova, manova, ancova, paired, nonreg sauf mode code)
# Changer: return = TRUE, verbose = FALSE, plot = FALSE
# En:      return = TRUE, verbose = TRUE, plot = TRUE

# Pattern 1: return = TRUE, verbose = FALSE, plot = FALSE
content2 <- gsub('return = TRUE, verbose = FALSE, plot = FALSE',
                 'return = TRUE, verbose = TRUE, plot = TRUE',
                 content2, fixed = TRUE)

# Pattern 2 (moins courant): return = FALSE avec verbose/plot FALSE
# ATTENTION: NE PAS modifier les tests de la catégorie "return" !
# On va faire ça manuellement ligne par ligne

# Compter APRÈS
n_after_verbose_true <- length(grep("verbose = TRUE", content2, fixed = TRUE))
n_after_plot_true <- length(grep("plot = TRUE", content2, fixed = TRUE))

cat("   Après: verbose=TRUE:", n_after_verbose_true, "| plot=TRUE:", n_after_plot_true, "\n")

writeLines(content2, file_path2)
cat("   ✓ .testeur_m.test.R corrigé\n\n")

# Résumé
cat("=== RÉSUMÉ DES CORRECTIONS ===\n")
cat("1. m.test.R: return(global_pvalue) → return(as.numeric(global_pvalue))\n")
cat("2. .testeur_m.test.R: Tests principaux utilisent maintenant plot=TRUE, verbose=TRUE\n")
cat("   (sauf tests spécifiques de return=FALSE qui restent avec plot/verbose FALSE)\n")
cat("\n✓ Corrections appliquées avec succès!\n")
