# Script pour corriger la double incrémentation du compteur

cat("=== CORRECTION DOUBLE INCRÉMENTATION COMPTEUR ===\n\n")

file_path <- "R/valreg.R"
content <- readLines(file_path, warn = FALSE)

# Trouver les lignes avec "counter <- counter + 1" qui sont suivies (dans les 10 lignes)
# par un .vbse(..., cpt="on")

# Stratégie: supprimer les lignes qui ont "counter <- counter + 1" et qui sont
# suivies d'un .vbse avec cpt="on"

lines_to_remove <- c()

for (i in 1:length(content)) {
  if (grepl("^[[:space:]]*counter <- counter \\+ 1[[:space:]]*$", content[i])) {
    # Vérifier les 10 lignes suivantes pour un .vbse avec cpt="on"
    found_vbse_on <- FALSE
    for (j in (i+1):min(i+10, length(content))) {
      if (grepl('cpt="on"', content[j], fixed=TRUE)) {
        found_vbse_on <- TRUE
        break
      }
    }

    if (found_vbse_on) {
      cat("Ligne redondante trouvée:", i, "->", content[i], "\n")
      lines_to_remove <- c(lines_to_remove, i)
    }
  }
}

if (length(lines_to_remove) > 0) {
  cat("\nSuppression de", length(lines_to_remove), "lignes redondantes\n")
  content <- content[-lines_to_remove]

  writeLines(content, file_path)
  cat("✓ Corrections appliquées\n")
} else {
  cat("⚠ Aucune ligne redondante trouvée\n")
}
