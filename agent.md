# Agent KefiR — Rôle et périmètre d’intervention

Le package **KefiR** doit être publié sur GitHub (nouvelle mise à jour à venir).

---

## 1. Principes généraux

- Modifier **le moins possible** les fichiers existants.
- Créer **le moins possible** de nouveaux fichiers.
- Conserver **le style, la structure, la documentation et la logique actuelle**.
- Les interventions peuvent concerner :
  - la documentation roxygen2,
  - DESCRIPTION / NAMESPACE,
  - métadonnées du package,
  - organisation interne,
  - scripts de test,
  - **sans forcément modifier un script R**.

---

## 2. Règles de développement — bp.log (AUTORITÉ)

Les règles officielles du projet sont définies dans :
dev/bp.log à lire nécessairement et régulièrement dans un contexte de modfication de scripts R.

Tu as le droit d'utiliser R pour tester les commandes. Si tu génères des scripts de tests tu les sauvegardes dans KefiR/dev et le nom du fichier doit commencer par "agent_"
Dans dev, tu crée un sous-dossier save où tu sauvegardes les différentes versions des scripts R modifier avec un historique de leur modification optimisé pour un retour en arrière.

Tu dois aussi créer un dossier test dans dev où tu crées des scripts R pour tester les fonctions que tu modifies. Ces scripts doivent être nommés de façon à ce qu'ils soient faciles à comprendre et à utiliser. Ils doivent aussi être commentés de façon à ce qu'ils soient faciles à comprendre et à utiliser.

Dans dev, load_all_kefir.R permet de tout charger comme si KefiR était installé afin de le tester.

.testeur.R est un script qui permet de tester les fonctions de KefiR. Il est commenté de façon à ce qu'il soit facile à comprendre et à utiliser et l'ensemble de ses tests (sauf ceux qui testent les messages error doivent donner des résultats satisfaisants

## 3. Règles de préparation du package : autorisation à utilise

Autorisation pour exécuter la commande devtools::check() sous R.

---

## 4. Guide de débogage rapide

### 4.1 Erreurs fréquentes et solutions

| Erreur | Cause probable | Solution |
|--------|----------------|----------|
| "valeur manquante là où TRUE/FALSE requis" | Comparaison avec NA dans if() | Utiliser `isTRUE()` ou `!is.na(x) && x` |
| "Plus de X groupes" incorrect | `unique(g)` inclut NA | Utiliser `unique(g[!is.na(g)])` |
| Post-hocs ignorés silencieusement | tryCatch masque l'erreur | Tester avec `debug=TRUE` |
| Titre graphique générique | Colonne lettres non trouvée | Vérifier logique de scan colonnes |

### 4.2 Fichiers clés à vérifier en priorité

```
R/m.test.R              # Point d'entrée, routage, graphiques
R/sys_posthoc.R         # Tests post-hoc (souvent source d'erreurs NA)
R/sys_one_factor_analysis.R  # Analyse 1 facteur, comptage groupes
R/sys_multi_factor_analysis.R  # Multi-facteurs, interactions
R/catego.R              # Conversion p-values → lettres groupement
R/sys_boots.R           # Bootstrap (protection NA)
```

### 4.3 Patterns de protection NA

```r
# Pattern 1: Comparaison dans if()
if (!is.na(pvals) && pvals <= alpha) { ... }

# Pattern 2: Comparaison booléenne stockée
signif <- if (is.na(pvals)) FALSE else (pvals <= alpha)
if (isTRUE(signif != other_signif)) { ... }

# Pattern 3: any() avec NA
if (any(vec1 != vec2, na.rm = TRUE)) { ... }

# Pattern 4: Comptage groupes
n_groups <- length(unique(g[!is.na(g)]))

# Pattern 5: Condition avec NULL
if (isTRUE(var %in% c("a", "b"))) { ... }  # var peut être NULL
```

### 4.4 Commandes de test rapide

```r
# Charger le package en développement
source("dev/load_all_kefir.R")

# Test multi-facteurs (sensible aux NA)
set.seed(123)
n <- 200
data <- data.frame(
  y = rnorm(n, 25, 5),
  F1 = factor(sample(c('A', 'B', 'C'), n, replace=TRUE)),
  F2 = factor(sample(c('X', 'Y'), n, replace=TRUE))
)
m.test(y ~ F1 * F2, data = data, verbose = TRUE)

# Test avec NA dans groupes
data$F1[1:20] <- NA
m.test(data$y, data$F1)

# Mode debug détaillé
m.test(y ~ F1, data = data, debug = TRUE)
```

### 4.5 Checklist avant commit

- [ ] `Rscript -e "parse(file='R/fichier.R')"` → pas d'erreur syntaxe
- [ ] Pas de `print()` de debug oubliés
- [ ] Protection NA sur toutes les comparaisons dans if()
- [ ] Tests .testeur.R passent
- [ ] bp.log mis à jour si nouvelle règle
