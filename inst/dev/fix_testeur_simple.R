# Script simple pour corriger .testeur_m.test.R

cat("=== CORRECTION SIMPLE .testeur_m.test.R ===\n\n")

file_path <- "R/.testeur_m.test.R"
content <- readLines(file_path, warn = FALSE)

# Compter AVANT
n_false_before <- length(grep("return = FALSE, verbose = FALSE, plot = FALSE", content, fixed = TRUE))
cat("AVANT: 'return = FALSE, verbose = FALSE, plot = FALSE' :", n_false_before, "\n\n")

# Ne PAS modifier les tests de return=FALSE (catégorie "return" et "formules")
# Modifier SEULEMENT ceux qui ont return = TRUE, verbose = FALSE, plot = FALSE

# Identifier les sections à NE PAS modifier
# Section formules: lignes 244-269
# Section return: lignes 275-306

# Simple stratégie: modifier ligne par ligne en vérifiant le numéro de ligne
formules_start <- 244
formules_end <- 269
return_start <- 275
return_end <- 306

for (i in 1:length(content)) {
  # Skip sections formules et return
  if ((i >= formules_start && i <= formules_end) ||
      (i >= return_start && i <= return_end)) {
    next
  }

  # Remplacer return = FALSE, verbose = FALSE, plot = FALSE
  # par return = FALSE, verbose = TRUE, plot = TRUE
  # (pour les tests hors formules/return qui utilisent return=FALSE)
  if (grepl("return = FALSE, verbose = FALSE, plot = FALSE", content[i], fixed = TRUE)) {
    content[i] <- gsub("return = FALSE, verbose = FALSE, plot = FALSE",
                       "return = FALSE, verbose = TRUE, plot = TRUE",
                       content[i], fixed = TRUE)
  }

  # Plus important: remplacer return = TRUE, verbose = FALSE, plot = FALSE
  if (grepl("return = TRUE, verbose = FALSE, plot = FALSE", content[i], fixed = TRUE)) {
    content[i] <- gsub("return = TRUE, verbose = FALSE, plot = FALSE",
                       "return = TRUE, verbose = TRUE, plot = TRUE",
                       content[i], fixed = TRUE)
  }
}

# Compter APRÈS
n_verbose_true <- length(grep("verbose = TRUE", content, fixed = TRUE))
n_plot_true <- length(grep("plot = TRUE", content, fixed = TRUE))

cat("APRÈS:\n")
cat("  verbose = TRUE:", n_verbose_true, "\n")
cat("  plot = TRUE:", n_plot_true, "\n\n")

# Sauvegarder
writeLines(content, file_path)

cat("✓ .testeur_m.test.R corrigé\n")
