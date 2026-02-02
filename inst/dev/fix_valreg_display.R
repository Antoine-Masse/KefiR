# Script pour corriger les problèmes d'affichage dans valreg.R

cat("=== CORRECTIONS AFFICHAGE valreg.R ===\n\n")

file_path <- "R/valreg.R"
content <- readLines(file_path, warn = FALSE)

# Problème 1: Message VIF en anglais de car::vif()
# Solution: Ajouter suppressMessages() autour de vif(reg)
cat("1. Suppression du message VIF en anglais...\n")

idx_vif <- grep("vif\\(reg\\) -> vif_reg", content, fixed = FALSE)
if (length(idx_vif) > 0) {
  cat("   Ligne trouvée:", idx_vif, "\n")
  cat("   Ancienne:", content[idx_vif], "\n")

  content[idx_vif] <- "      vif_reg <- suppressMessages(vif(reg))"

  cat("   Nouvelle:", content[idx_vif], "\n\n")
} else {
  cat("   ⚠ Ligne non trouvée\n\n")
}

# Problème 2: P-value DW non arrondie
# Solution: Utiliser .format_pval() ou round()
cat("2. Formatage de la p-value Durbin-Watson...\n")

# Trouver les lignes avec pvalt non formaté dans DW
idx_dw1 <- grep('paste0\\("Durbin-Watson test.*p.value.*informative only.*", pvalt\\)', content)
idx_dw2 <- grep('paste0\\("Test de Durbin-Watson.*informatif seulement.*", pvalt\\)', content)

if (length(idx_dw1) > 0) {
  cat("   Ligne anglaise trouvée:", idx_dw1, "\n")
  cat("   Ancienne:", content[idx_dw1], "\n")

  content[idx_dw1] <- '        .vbse(ang=paste0("Durbin-Watson test (dwtest()) - p.value (informative only): ", .format_pval(pvalt)),'

  cat("   Nouvelle:", content[idx_dw1], "\n")
}

if (length(idx_dw2) > 0) {
  cat("   Ligne française trouvée:", idx_dw2, "\n")
  cat("   Ancienne:", content[idx_dw2], "\n")

  content[idx_dw2] <- '              fr=paste0("Test de Durbin-Watson (dwtest()) - p.value (informatif seulement) : ", .format_pval(pvalt)),'

  cat("   Nouvelle:", content[idx_dw2], "\n\n")
}

# Sauvegarder
writeLines(content, file_path)

cat("✓ Corrections appliquées\n")
cat("\nRésumé des modifications:\n")
cat("  1. vif(reg) → suppressMessages(vif(reg))\n")
cat("  2. pvalt → .format_pval(pvalt) dans message DW\n")
