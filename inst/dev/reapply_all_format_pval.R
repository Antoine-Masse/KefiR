# Script pour réappliquer .format_pval() à TOUTES les p-values dans valreg.R

file_path <- "/mnt/c/Users/masse/Desktop/KefiR/KefiR/R/valreg.R"
content <- readLines(file_path, warn = FALSE)

# 1. Dans les .vbse() - Pattern: p-value: ", pval_mdl) → p-value: ", .format_pval(pval_mdl))
content <- gsub('p-value: ", pval_mdl\\)', 'p-value: ", .format_pval(pval_mdl))', content)
content <- gsub('p-value : ", pval_mdl\\)', 'p-value : ", .format_pval(pval_mdl))', content)

# 2. Dans les .vbse() - Pattern: p.value: ", pvalt) → p.value: ", .format_pval(pvalt))
content <- gsub('p\\.value: ", pvalt\\)', 'p.value: ", .format_pval(pvalt))', content)
content <- gsub('p\\.value : ", pvalt\\)', 'p.value : ", .format_pval(pvalt))', content)

# 3. Dans les .vbse() - Pattern: max(p.value): ", max(pval_coeff)) → max(p.value): ", .format_pval(max(pval_coeff)))
content <- gsub('max\\(p\\.value\\): ", max\\(pval_coeff\\)\\)', 'max(p.value): ", .format_pval(max(pval_coeff)))', content)
content <- gsub('max\\(p\\.value\\) : ", max\\(pval_coeff\\)\\)', 'max(p.value) : ", .format_pval(max(pval_coeff)))', content)
content <- gsub('max\\(pval_coeff\\): ", max\\(pval_coeff\\)\\)', 'max(pval_coeff): ", .format_pval(max(pval_coeff)))', content)
content <- gsub('max\\(pval_coeff\\) : ", max\\(pval_coeff\\)\\)', 'max(pval_coeff) : ", .format_pval(max(pval_coeff)))', content)

# 4. Dans les if(verbose) cat() - Pattern: ", .format_pval(pvalt), "n") → ", .format_pval(pvalt), "\n")
content <- gsub(', \\.format_pval\\(pvalt\\), "n"\\)', ', .format_pval(pvalt), "\\n")', content)

# 5. VIF et Cook's distance - arrondir les valeurs
content <- gsub('max\\(vif_reg, na.rm = TRUE\\), "\\\\n"\\)', 'round(max(vif_reg, na.rm = TRUE), 2), "\\n")', content)
content <- gsub('max\\(cooksd, na.rm = TRUE\\), "\\\\n"\\)', 'round(max(cooksd, na.rm = TRUE), 4), "\\n")', content)

# Sauvegarder
writeLines(content, file_path)

cat("✅ .format_pval() réappliqué à TOUTES les p-values\n")
cat("✅ round() appliqué aux VIF et Cook's distance\n")
cat("✅ Bug de syntaxe corrigé (, \"n\") → , \"\\n\"))\n")
cat("Le fichier valreg.R a été mis à jour.\n")
