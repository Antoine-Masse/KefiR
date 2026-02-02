# Analyse du problème : `.multi_factor_analysis()` pour mesures répétées

## Date: 2025-11-05

## 🐛 Problème identifié

### Test effectué:
```r
m.test(A~F*G, id="H", data=data)
# Routage : ANOVA multi-facteurs (repeated)
```

### Bilan reçu - Problèmes constatés:

1. **Test de Friedman non exécuté**
   - Suggestion : "Test de Friedman (friedman.test) pour k >= 3 conditions"
   - **Mais aucune p-value affichée**
   - Aucun résultat concret

2. **Erreur ezANOVA**
   ```
   Impossible d'effectuer ezANOVA pour le test de sphéricité : objet 'DV' introuvable
   ```

3. **Messages confus**
   ```
   Passage vers une approche d'analyse robuste...
   Approche robuste suggérée pour mesures répétées :
   - Test de Friedman...
   - Mesures répétées robustes du package {WRS2}
   - Modèle mixte avec estimation robuste (lmer avec REML)
   Note : L'implémentation de ces tests n'est PAS ENCORE complètement automatique
   ```

4. **Format incohérent avec `.one_factor_analysis()`**
   - `.one_factor_analysis()` : affiche le test exécuté + p-value claire
   - `.multi_factor_analysis()` : liste de suggestions sans action

---

## 🔍 Cause racine

**Fichier**: `R/.multi_factor_analysis.R` lignes 1900-1964

### NOTE PERSO #11 (code source):
```r
#* NOTE PERSO #11 : Implémentation automatique des tests robustes
#*
#* ➜ Problème identifié :
#*   Pour l'instant, la fonction suggère des démarches mais ne les applique PAS.
#*   Son ambition est de faire TOUS ces tests automatiquement, de contrôler
#*   les assomptions à chaque fois, et de changer de stratégie selon les contrôles.
#*
#*   Exemples d'implémentations manquantes :
#*   1. Friedman test (k >= 3, un seul facteur, appariement)
#*   2. Modèles mixtes (lmer) avec critères automatiques
#*   3. ANOVA par permutation (lmPerm::aovp)
#*   4. Tests robustes WRS2 (t2way, etc.)
#*   5. Scheirer-Ray-Hare pour 2-way non paramétrique
```

**Conclusion**: C'est un **travail en cours (NON RÉSOLU)** selon le code lui-même.

---

## 📊 Comparaison avec `.one_factor_analysis()`

### `.one_factor_analysis()` (CLAIR) ✅

```
1) Analyse des p-values du modèle et de ses coefficients.
   Test de Friedman appliqué.
   p-value : 0.0023
   ==> Différences significatives détectées.

2) Post-hocs: Tests de Conover avec correction de Bonferroni.
   Comparaison A vs B : p = 0.041
   Comparaison A vs C : p = 0.003
```

- Test exécuté
- P-value affichée
- Interprétation claire
- Post-hocs effectués

### `.multi_factor_analysis()` (CONFUS) ❌

```
9) Passage vers une approche d'analyse robuste en raison de :
   - Violations d'hypothèses, OU
   - Déséquilibre sévère, OU
   - Variable dépendante discrète, OU
   - Plan de mesures répétées déséquilibré, OU
   - Violations d'hypothèses ANCOVA
   Approche robuste suggérée pour mesures répétées :
   - Test de Friedman (friedman.test) pour k >= 3 conditions
   - Mesures répétées robustes du package {WRS2}
   - Modèle mixte avec estimation robuste (lmer avec REML)
   Note : L'implémentation de ces tests n'est PAS ENCORE complètement automatique
```

- **Aucun test exécuté**
- **Aucune p-value**
- **Juste des suggestions**
- **Message "pas encore automatique"**

---

## ✅ Solutions proposées

### Solution 1: PRIORITAIRE - Implémenter Friedman pour multi-facteurs

Pour le cas spécifique `F*G` avec mesures répétées, on peut :

1. **Détecter** : 2 facteurs within, k >= 3 niveaux par facteur
2. **Exécuter** : Test de Friedman sur interaction `F:G`
3. **Afficher** :
   ```
   9) Test de Friedman pour mesures répétées (facteurs F × G).
      Statistique de Friedman : χ² = 12.34, df = 5
      p-value : 0.029
      ==> Différences significatives détectées entre conditions.
   ```

### Solution 2: Améliorer le message "approche robuste"

Au lieu de :
```
Approche robuste suggérée...
Note : L'implémentation de ces tests n'est PAS ENCORE complètement automatique
```

Faire :
```
10) Application du test de Friedman pour mesures répétées multi-facteurs.
    [Exécution réelle du test]
    p-value : X.XX
    Interprétation : [...]
```

### Solution 3: Corriger l'erreur ezANOVA

L'erreur "objet 'DV' introuvable" vient probablement d'un problème de variables dans l'environnement.

Chercher dans le code où ezANOVA est appelé et corriger la préparation des données.

---

## 🎯 Recommandations immédiates

### Court terme (Critical):
1. ✅ Implémenter Friedman automatique pour cas multi-facteurs avec mesures répétées
2. ✅ Afficher p-value et résultat clair (comme `.one_factor_analysis()`)
3. ✅ Corriger erreur ezANOVA

### Moyen terme:
4. Implémenter modèles mixtes (lmer) avec critères automatiques
5. Implémenter Scheirer-Ray-Hare pour 2-way non paramétrique

### Long terme:
6. ANOVA par permutation (lmPerm)
7. Tests robustes WRS2 complets

---

## 📝 Plan d'action proposé

### Étape 1: Identifier où ajouter le code Friedman

Dans `.multi_factor_analysis.R` autour des lignes 1868-1893, après :
```r
k <- .vbse(
  "Switching to robust analysis approach...",
  "Passage vers une approche d'analyse robuste...",
  ...
)
```

**Ajouter**:
```r
# Cas spécial : mesures répétées multi-facteurs
if (paired && n_factors >= 1) {
  k <- .vbse(
    "Applying Friedman test for repeated measures (multi-factor design)...",
    "Application du test de Friedman pour mesures répétées (design multi-facteurs)...",
    verbose = verbose, k = k, cpt = "on"
  )

  # Exécuter Friedman
  friedman_result <- friedman.test(x ~ interaction(factor_vars) | id_var, data = data)

  k <- .vbse(
    paste0("Friedman test statistic: χ² = ", round(friedman_result$statistic, 2),
           ", df = ", friedman_result$parameter,
           "\n\tp-value: ", .format_pval(friedman_result$p.value)),
    paste0("Statistique du test de Friedman : χ² = ", round(friedman_result$statistic, 2),
           ", df = ", friedman_result$parameter,
           "\n\tp-value : ", .format_pval(friedman_result$p.value)),
    verbose = verbose, k = k, cpt = "off"
  )

  # Interprétation
  if (friedman_result$p.value < alpha) {
    k <- .vbse(
      "Significant differences detected between repeated measures conditions.",
      "Différences significatives détectées entre les conditions en mesures répétées.",
      verbose = verbose, k = k, cpt = "off"
    )
  } else {
    k <- .vbse(
      "No significant differences between repeated measures conditions.",
      "Pas de différences significatives entre les conditions en mesures répétées.",
      verbose = verbose, k = k, cpt = "off"
    )
  }

  # Stocker résultat
  robust_results$method <- "Friedman"
  robust_results$test_result <- friedman_result
  robust_results$posthoc_applicable <- (friedman_result$p.value < alpha)
}
```

---

## 🔧 Besoin d'aide ?

Voulez-vous que j'implémente ces corrections dans `.multi_factor_analysis.R` ?

Les priorités seraient :
1. Friedman automatique pour multi-facteurs repeated
2. Affichage clair avec p-value
3. Correction erreur ezANOVA

**Status**: Problème identifié, solutions proposées. Prêt pour implémentation.
