# Script pour appliquer .format_pval() aux p-values dans les messages if(verbose) cat() restants

file_path <- "/mnt/c/Users/masse/Desktop/KefiR/KefiR/R/valreg.R"
content <- readLines(file_path, warn = FALSE)

# Pattern 1: p.value:"), pvalt, "\n") → p.value: ", .format_pval(pvalt), "\n")
content <- gsub('p\\.value:"), pvalt, "\\\\n"\\)', 'p.value: ", .format_pval(pvalt), "\\n")', content)

# Pattern 2: p.value :"), pvalt, "\n") → p.value : ", .format_pval(pvalt), "\n")
content <- gsub('p\\.value :"), pvalt, "\\\\n"\\)', 'p.value : ", .format_pval(pvalt), "\\n")', content)

# Pattern 3: max(cooksd, na.rm = TRUE), "\n") → round(max(cooksd, na.rm = TRUE), 4), "\n")
content <- gsub('max\\(cooksd, na.rm = TRUE\\), "\\\\n"\\)', 'round(max(cooksd, na.rm = TRUE), 4), "\\n")', content)

# Pattern 4: max(vif_reg, na.rm = TRUE), "\n") → round(max(vif_reg, na.rm = TRUE), 2), "\n")
content <- gsub('max\\(vif_reg, na.rm = TRUE\\), "\\\\n"\\)', 'round(max(vif_reg, na.rm = TRUE), 2), "\\n")', content)

# Sauvegarder
writeLines(content, file_path)

cat("✅ .format_pval() et round() appliqués aux messages if(verbose) cat() restants\n")
cat("Le fichier valreg.R a été mis à jour.\n")
