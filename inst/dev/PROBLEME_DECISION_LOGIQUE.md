# Problème de logique de décision dans `.multi_factor_analysis()`

## Date: 2025-11-05

## 🎯 Constat de l'utilisateur (CORRECT)

> "Je ne dis pas qu'il faudrait Friedman ici... Je trouve que le bilan n'est pas clair...
> À priori, j'aurais tendance à penser que c'est le modèle mixte car H n'est pas renseigné pour chaque F,G....
> Mais en 9 cela parlait encore de Friedman robuste : commentaire type pot commun."

---

## 🐛 Problème identifié

### Votre cas : `m.test(A~F*G, id="H", data=data)`

**Structure réelle des données :**
- Multi-facteurs : F × G (interaction)
- Mesures répétées : id = H
- **Déséquilibré** : "5 sujet(s) n'ont pas le nombre attendu de mesures"
- **Doublons** : "Certaines combinaisons id × facteur within ont plusieurs observations"
- H n'est **pas croisé** avec toutes les combinaisons F×G

**Test approprié** : **Modèle mixte (lmer)** ✅

**Ce que fait le code** : Affiche un message générique "pot-pourri" :
```
Approche robuste suggérée pour mesures répétées :
- Test de Friedman (friedman.test) pour k >= 3 conditions  ❌ PAS ADAPTÉ
- Mesures répétées robustes du package {WRS2}              ❌ PAS ADAPTÉ
- Modèle mixte avec estimation robuste (lmer avec REML)    ✅ ADAPTÉ
Note : L'implémentation de ces tests n'est PAS ENCORE complètement automatique
```

---

## 🔍 Cause dans le code

**Fichier** : `R/.multi_factor_analysis.R` lignes 2119-2140

```r
} else {
  # Cas sans id ou multi-facteurs : message générique  ⬅️ PROBLÈME ICI
  k <- .vbse(
    "Suggested robust approach for repeated measures:
     - Friedman test (friedman.test) for k >= 3 conditions
     - Robust repeated measures from {WRS2} package
     - Mixed model with robust estimation (lmer with REML)
     Note: Implementation of these tests is NOT YET fully automatic"
  )

  robust_results$method <- "Repeated_Measures_Suggested"  ⬅️ AUCUN TEST EXÉCUTÉ
}
```

**Problème** : Le code n'analyse PAS la structure des données pour choisir LA bonne méthode.

---

## 📊 Décision intelligente attendue

### Critères de décision :

| Situation | Test approprié | Critères |
|-----------|----------------|----------|
| **1 facteur, équilibré, complet** | Friedman | - 1 facteur within<br>- Chaque sujet a toutes les conditions<br>- Pas de données manquantes |
| **1 facteur, déséquilibré** | lmer | - 1 facteur within<br>- Données manquantes >5% |
| **Multi-facteurs, équilibré** | ezANOVA / aov | - 2+ facteurs<br>- Design complet<br>- Sphéricité OK |
| **Multi-facteurs, déséquilibré** | **lmer** ✅ | - 2+ facteurs<br>- Données manquantes<br>- Design incomplet<br>- Doublons |
| **Multi-facteurs, violations sphéricité** | lmer | - Epsilon < 0.75<br>- Correction GG/HF insuffisante |

### Votre cas → **lmer** car :
- ✅ Multi-facteurs (F×G)
- ✅ Déséquilibré (pas toutes les combinaisons)
- ✅ Doublons détectés
- ✅ Design incomplet (H pas croisé avec F×G)

---

## ✅ Solution proposée

### Étape 1 : Analyser la structure AVANT de décider

Ajouter dans `.multi_factor_analysis.R` autour de la ligne 2119 :

```r
} else {
  # ANALYSER LA STRUCTURE DES DONNÉES

  # 1. Vérifier équilibrage
  cell_counts <- table(interaction(factor_vars), id_var)
  is_balanced <- (length(unique(as.vector(cell_counts))) == 1)
  missing_pct <- sum(cell_counts == 0) / length(cell_counts)
  has_duplicates <- any(cell_counts > 1)

  # 2. Compter facteurs
  n_factors <- length(factor_vars)

  # 3. DÉCISION INTELLIGENTE
  if (n_factors == 1 && is_balanced && !has_duplicates) {
    # CAS 1: Friedman approprié
    method_choice <- "friedman"

  } else if (n_factors >= 2 || !is_balanced || has_duplicates || missing_pct > 0.05) {
    # CAS 2: Modèle mixte nécessaire
    method_choice <- "lmer"

  } else {
    # CAS 3: ezANOVA possible
    method_choice <- "ezanova"
  }
```

### Étape 2 : Exécuter LE bon test

```r
  # EXÉCUTER LE TEST CHOISI
  if (method_choice == "lmer") {
    k <- .vbse(
      paste0("Mixed model selected (design characteristics: ",
             n_factors, " factors, ",
             ifelse(is_balanced, "balanced", "unbalanced"),
             ifelse(has_duplicates, ", with duplicates", ""), ")."),
      paste0("Modèle mixte sélectionné (caractéristiques du design : ",
             n_factors, " facteurs, ",
             ifelse(is_balanced, "équilibré", "déséquilibré"),
             ifelse(has_duplicates, ", avec doublons", ""), ")."),
      verbose = verbose, k = k, cpt = "on"
    )

    # Construire formule lmer
    formula_lmer <- as.formula(paste0(
      deparse(formula[[2]]), " ~ ",
      paste(factor_vars, collapse = " * "),
      " + (1|", id_var, ")"
    ))

    # Exécuter lmer
    require(lme4)
    require(lmerTest)
    lmer_model <- lmer(formula_lmer, data = data)
    lmer_anova <- anova(lmer_model)

    # Afficher résultats
    k <- .vbse(
      paste0("Mixed model results (Type III):"),
      paste0("Résultats du modèle mixte (Type III) :"),
      verbose = verbose, k = k, cpt = "off"
    )

    for (i in 1:nrow(lmer_anova)) {
      effect_name <- rownames(lmer_anova)[i]
      f_val <- lmer_anova[i, "F value"]
      p_val <- lmer_anova[i, "Pr(>F)"]

      k <- .vbse(
        paste0("\t", effect_name, ": F = ", round(f_val, 2),
               ", p = ", .format_pval(p_val)),
        paste0("\t", effect_name, " : F = ", round(f_val, 2),
               ", p = ", .format_pval(p_val)),
        verbose = verbose, k = k, cpt = "off"
      )
    }

    # Stocker résultat
    robust_results$method <- "Mixed_Model_lmer"
    robust_results$test_result <- lmer_anova
    robust_results$model <- lmer_model

  } else if (method_choice == "friedman") {
    # [Code Friedman existant]

  } else {
    # [Code ezANOVA]
  }
}
```

### Étape 3 : Message clair et actionnable

Au lieu de :
```
Approche robuste suggérée... [liste de tout]
Note : pas encore automatique
```

Afficher :
```
9) Modèle mixte sélectionné (design : 2 facteurs, déséquilibré, avec doublons).
10) Résultats du modèle mixte (Type III) :
    F : F = 3.45, p = 0.042
    G : F = 5.67, p = 0.008
    F:G : F = 2.12, p = 0.098
    ==> Effets principaux de F et G significatifs.
```

---

## 🎯 Bénéfices de la solution

1. **Décision automatique intelligente** basée sur la structure des données
2. **Test exécuté** au lieu de suggestions
3. **P-values affichées** pour chaque effet
4. **Interprétation claire** pour débutant en statistiques
5. **Format cohérent** avec `.one_factor_analysis()`

---

## 📝 Résumé

| Aspect | Actuellement ❌ | Proposé ✅ |
|--------|----------------|-----------|
| Décision | Aucune | Automatique basée sur structure |
| Test | Aucun | Exécuté (lmer, Friedman, ou ezANOVA) |
| P-value | Absente | Affichée pour chaque effet |
| Message | "Pot-pourri" confus | Choix justifié + résultats |
| Utilisabilité | Nécessite expert | Accessible débutant |

**Voulez-vous que j'implémente cette logique de décision intelligente dans `.multi_factor_analysis.R` ?**
