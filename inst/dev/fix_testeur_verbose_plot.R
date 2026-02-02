# Script pour corriger .testeur_m.test.R : mettre verbose=TRUE, plot=TRUE
# SAUF pour les tests de la catégorie "return" et les tests de "code=TRUE"

cat("=== CORRECTION .testeur_m.test.R ===\n\n")

file_path <- "R/.testeur_m.test.R"
content <- readLines(file_path, warn = FALSE)

# Compter AVANT
n_verbose_false_before <- length(grep("verbose = FALSE", content, fixed = TRUE))
n_plot_false_before <- length(grep("plot = FALSE", content, fixed = TRUE))
cat("AVANT:\n")
cat("  verbose = FALSE:", n_verbose_false_before, "\n")
cat("  plot = FALSE:", n_plot_false_before, "\n\n")

# Stratégie: Remplacer TOUS les ", verbose = FALSE, plot = FALSE)"
# par ", verbose = TRUE, plot = TRUE)"
# SAUF dans la section "return" (lignes 275-306)

# Identifier les lignes de la section "return" à NE PAS modifier
return_section_start <- grep("CATÉGORIE 2 : return=FALSE", content, fixed = TRUE)
return_section_end <- grep("CATÉGORIE 3 : ANOVA", content, fixed = TRUE)

cat("Section 'return' : lignes", return_section_start, "à", return_section_end, "\n\n")

# Parcourir toutes les lignes et remplacer sauf dans section "return"
for (i in 1:length(content)) {
  # Ne pas modifier les lignes de la section "return"
  if (i >= return_section_start && i <= return_section_end) {
    next
  }

  # Remplacer verbose = FALSE par verbose = TRUE
  if (grepl("verbose = FALSE", content[i], fixed = TRUE)) {
    content[i] <- gsub("verbose = FALSE", "verbose = TRUE", content[i], fixed = TRUE)
  }

  # Remplacer plot = FALSE par plot = TRUE
  if (grepl("plot = FALSE", content[i], fixed = TRUE)) {
    content[i] <- gsub("plot = FALSE", "plot = TRUE", content[i], fixed = TRUE)
  }
}

# Compter APRÈS
n_verbose_false_after <- length(grep("verbose = FALSE", content, fixed = TRUE))
n_plot_false_after <- length(grep("plot = FALSE", content, fixed = TRUE))
n_verbose_true_after <- length(grep("verbose = TRUE", content, fixed = TRUE))
n_plot_true_after <- length(grep("plot = TRUE", content, fixed = TRUE))

cat("APRÈS:\n")
cat("  verbose = FALSE:", n_verbose_false_after, "\n")
cat("  plot = FALSE:", n_plot_false_after, "\n")
cat("  verbose = TRUE:", n_verbose_true_after, "\n")
cat("  plot = TRUE:", n_plot_true_after, "\n\n")

# Sauvegarder
writeLines(content, file_path)

cat("✓ .testeur_m.test.R corrigé\n")
cat("✓", n_verbose_false_before - n_verbose_false_after, "lignes verbose modifiées\n")
cat("✓", n_plot_false_before - n_plot_false_after, "lignes plot modifiées\n")
