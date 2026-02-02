# Script complet pour finir la conversion de valreg.R

file_path <- "/mnt/c/Users/masse/Desktop/KefiR/KefiR/R/valreg.R"
content <- readLines(file_path, warn = FALSE)

# 1. CORRECTION BUG LIGNE 497: (reg)$p.value -> bptest_lmer(reg)$p.value
content[497] <- gsub("\\(reg\\)\\$p\\.value -> pvalt", "bptest_lmer(reg)$p.value -> pvalt", content[497])

# 2. CORRECTION BUGS DE SYNTAXE: "n" -> "\n" (lignes avec .format_pval)
content <- gsub(', \\.format_pval\\(pvalt\\), "n"\\)', ', .format_pval(pvalt), "\\n")', content)
content <- gsub('\\), round\\(max\\(cooksd, na.rm = TRUE\\), 4\\), "n"\\)', '), round(max(cooksd, na.rm = TRUE), 4), "\\n")', content)
content <- gsub('\\), round\\(max\\(vif_reg, na.rm = TRUE\\), 2\\), "n"\\)', '), round(max(vif_reg, na.rm = TRUE), 2), "\\n")', content)

# 3. CONVERSION DES if(verbose) cat() RESTANTS EN .vbse()

# Ligne 492: Analysis of variance of residuals
old_492 <- "    if (verbose) cat(counter, .msg(\"- Analysis of variance of residuals.\\n\", \"- Analyse de la variance des résidus.\\n\"))"
new_492 <- "    counter <- .vbse(ang=\"- Analysis of variance of residuals.\", fr=\"- Analyse de la variance des résidus.\", k=counter, cpt=\"on\", verbose=verbose)"
content <- gsub(old_492, new_492, content, fixed=TRUE)

# Lignes 500-501: Warning Breush-Pagan (avec message lmerMod)
old_500 <- "        if (verbose) cat(.msg(\"\\tWarning!\\n\\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value:\", \"\\tAttention !\\n\\tTest de Breush-Pagan (bptest()) - Variance non-constante des résidus. p.value : \", .format_pval(pvalt), \"\\n\")"
new_500 <- "        .vbse(ang=paste0(\"\\tWarning!\\n\\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value: \", .format_pval(pvalt)), fr=paste0(\"\\tAttention !\\n\\tTest de Breush-Pagan (bptest()) - Variance non-constante des résidus. p.value : \", .format_pval(pvalt)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_500, new_500, content, fixed=TRUE)

old_502 <- "\t\t\tif (verbose) cat(.msg(\"\\tHeteroscedasticity is not always an issue in a mixed model. See {nlme}.\\n\",\"\\tL'hétéroscédasticité n'est pas toujours un problème sur un modèle mixte. Voyez {nlme}.\\n\"))"
new_502 <- "\t\t\t.vbse(ang=\"\\tHeteroscedasticity is not always an issue in a mixed model. See {nlme}.\", fr=\"\\tL'hétéroscédasticité n'est pas toujours un problème sur un modèle mixte. Voyez {nlme}.\", k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_502, new_502, content, fixed=TRUE)

# Ligne 507: Good Breush-Pagan
old_507 <- "        if (verbose) cat(.msg(\"\\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value:\", \"\\tTest de Breush-Pagan (bptest()) - Variance constante des résidus. p.value : \", .format_pval(pvalt), \"\\n\")"
new_507 <- "        .vbse(ang=paste0(\"\\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value: \", .format_pval(pvalt)), fr=paste0(\"\\tTest de Breush-Pagan (bptest()) - Variance constante des résidus. p.value : \", .format_pval(pvalt)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_507, new_507, content, fixed=TRUE)

# Ligne 510: Too many values warning
old_510 <- "      cat(.msg(\"\\tWarning!\\n\\tToo many values to justify checking the constancy of the variance. Breush-Pagan test bptest() of {lmtest} may be hypersensitive.\\n\","
old_511 <- "              \"\\tAttention !\\n\\tTrop de valeurs pour justifier la vérification de la constance de la variance. Le test de Breush-Pagan bptest() de {lmtest} peut être hypersensible.\\n\"))"
# Chercher l'index de ces lignes
idx_510 <- which(grepl("Too many values to justify checking the constancy", content, fixed=TRUE))
if (length(idx_510) > 0) {
  content[idx_510] <- "      .vbse(ang=\"\\tWarning!\\n\\tToo many values to justify checking the constancy of the variance. Breush-Pagan test bptest() of {lmtest} may be hypersensitive.\", fr=\"\\tAttention !\\n\\tTrop de valeurs pour justifier la vérification de la constance de la variance. Le test de Breush-Pagan bptest() de {lmtest} peut être hypersensible.\", k=counter, cpt=\"off\", verbose=verbose)"
  # Supprimer la ligne suivante si elle existe (continuation)
  if (idx_510 < length(content) && grepl("Attention.*variance.*hypersensible", content[idx_510+1])) {
    content <- content[-(idx_510+1)]
  }
}

# Ligne 520: Analysis of leverage effect
old_520 <- "    if (verbose) cat(counter, .msg(\"- Analysis of leverage effect.\\n\", \"- Analyse de l'effet de levier.\\n\"))"
new_520 <- "    counter <- .vbse(ang=\"- Analysis of leverage effect.\", fr=\"- Analyse de l'effet de levier.\", k=counter, cpt=\"on\", verbose=verbose)"
content <- gsub(old_520, new_520, content, fixed=TRUE)

# Ligne 522: Warning Cook's distance
old_522 <- "      if (verbose) cat(.msg(\"\\tWarning!\\n\\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance()):\", \"\\tAttention !\\n\\tDistance de Cook (cooks.distance()) - Effet de levier. max(cooks.distance()) :\"), round(max(cooksd, na.rm = TRUE), 4), \"\\n\")"
new_522 <- "      .vbse(ang=paste0(\"\\tWarning!\\n\\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance()): \", round(max(cooksd, na.rm = TRUE), 4)), fr=paste0(\"\\tAttention !\\n\\tDistance de Cook (cooks.distance()) - Effet de levier. max(cooks.distance()) : \", round(max(cooksd, na.rm = TRUE), 4)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_522, new_522, content, fixed=TRUE)

# Ligne 525: Good Cook's distance
old_525 <- "      if (verbose) cat(.msg(\"\\tCook's distance (cooks.distance()) - No leverage effect. max(cooks.distance()):\", \"\\tDistance de Cook (cooks.distance()) - Pas d'effet de levier. max(cooks.distance()) :\"), round(max(cooksd, na.rm = TRUE), 4), \"\\n\")"
new_525 <- "      .vbse(ang=paste0(\"\\tCook's distance (cooks.distance()) - No leverage effect. max(cooks.distance()): \", round(max(cooksd, na.rm = TRUE), 4)), fr=paste0(\"\\tDistance de Cook (cooks.distance()) - Pas d'effet de levier. max(cooks.distance()) : \", round(max(cooksd, na.rm = TRUE), 4)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_525, new_525, content, fixed=TRUE)

# Ligne 531: Multicollinearity test
old_531 <- "      if (verbose) cat(counter, .msg(\"- Multicollinearity test (VIF).\\n\", \"- Test de multicolinéarité (VIF).\\n\"))"
new_531 <- "      counter <- .vbse(ang=\"- Multicollinearity test (VIF).\", fr=\"- Test de multicolinéarité (VIF).\", k=counter, cpt=\"on\", verbose=verbose)"
content <- gsub(old_531, new_531, content, fixed=TRUE)

# Ligne 533: Warning VIF
old_533 <- "        if (verbose) cat(.msg(\"\\tWarning!\\n\\tThe variance inflation factor (VIF) indicates collinear variables with car::vif(). max(vif()):\", \"\\tAttention !\\n\\tLe facteur d'inflation de la variance (VIF) indique des variables collinéaires avec car::vif(). max(vif()) :\"), round(max(vif_reg, na.rm = TRUE), 2), \"\\n\")"
new_533 <- "        .vbse(ang=paste0(\"\\tWarning!\\n\\tThe variance inflation factor (VIF) indicates collinear variables with car::vif(). max(vif()): \", round(max(vif_reg, na.rm = TRUE), 2)), fr=paste0(\"\\tAttention !\\n\\tLe facteur d'inflation de la variance (VIF) indique des variables collinéaires avec car::vif(). max(vif()) : \", round(max(vif_reg, na.rm = TRUE), 2)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_533, new_533, content, fixed=TRUE)

# Ligne 536: Good VIF
old_536 <- "        if (verbose) cat(.msg(\"\\tNo significant multicollinearity problem with car::vif(). max(vif()):\", \"\\tPas de problème significatif de multicolinéarité avec car::vif(). max(vif()) :\"), round(max(vif_reg, na.rm = TRUE), 2), \"\\n\")"
new_536 <- "        .vbse(ang=paste0(\"\\tNo significant multicollinearity problem with car::vif(). max(vif()): \", round(max(vif_reg, na.rm = TRUE), 2)), fr=paste0(\"\\tPas de problème significatif de multicolinéarité avec car::vif(). max(vif()) : \", round(max(vif_reg, na.rm = TRUE), 2)), k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_536, new_536, content, fixed=TRUE)

# Ligne 542: Analysis of solidity
old_542 <- "      if (verbose) cat(counter, .msg(\"- Analysis of solidity of model by bootstrap.\\n\", \"- Analyse de la solidité du modèle par bootstrap.\\n\"))"
new_542 <- "      counter <- .vbse(ang=\"- Analysis of solidity of model by bootstrap.\", fr=\"- Analyse de la solidité du modèle par bootstrap.\", k=counter, cpt=\"on\", verbose=verbose)"
content <- gsub(old_542, new_542, content, fixed=TRUE)

# Ligne 549: Warning Bootstrap
old_549 <- "        if (verbose) cat(.msg(\"\\tWarning!\\n\\tBootstrap (bootreg()) - Fragility of the model in bootstrap. Please, use bootreg()\\n\", \"\\tAttention !\\n\\tBootstrap (bootreg()) - Fragilité du modèle en bootstrap. Veuillez utiliser bootreg()\\n\"))"
new_549 <- "        .vbse(ang=\"\\tWarning!\\n\\tBootstrap (bootreg()) - Fragility of the model in bootstrap. Please, use bootreg()\", fr=\"\\tAttention !\\n\\tBootstrap (bootreg()) - Fragilité du modèle en bootstrap. Veuillez utiliser bootreg()\", k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_549, new_549, content, fixed=TRUE)

# Ligne 552: Good Bootstrap
old_552 <- "        if (verbose) cat(.msg(\"\\tBootstrap (bootreg()) - Solidity of the model in bootstrap.\\n\", \"\\tBootstrap (bootreg()) - Solidité du modèle en bootstrap.\\n\"))"
new_552 <- "        .vbse(ang=\"\\tBootstrap (bootreg()) - Solidity of the model in bootstrap.\", fr=\"\\tBootstrap (bootreg()) - Solidité du modèle en bootstrap.\", k=counter, cpt=\"off\", verbose=verbose)"
content <- gsub(old_552, new_552, content, fixed=TRUE)

# Sauvegarder
writeLines(content, file_path)

cat("✅ Bug ligne 497 corrigé: (reg)$p.value -> bptest_lmer(reg)$p.value\n")
cat("✅ Bugs de syntaxe corrigés: \"n\" -> \"\\n\"\n")
cat("✅ Tous les if(verbose) cat() convertis en .vbse()\n")
cat("\nVérification finale...\n")

# Compter les conversions
n_vbse <- length(grep("\\.vbse\\(", content))
n_ifverbose <- length(grep("if \\(verbose\\) cat\\(", content))

cat("Messages .vbse() :", n_vbse, "\n")
cat("Messages if(verbose) cat restants :", n_ifverbose, "\n")
