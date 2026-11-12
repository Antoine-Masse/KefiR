# KefiR 0.0.2.3

## Corrections de `m.test()` (audit de novembre 2026)

Détail, mesures avant et après : `dev/mtest_audit/RAPPORT.md`.

### Résultats qui changent pour les étudiants

* **Plan split-plot** `y ~ G*T + Error(id/T)` : l'interaction G:T n'est plus retirée sans message. Le facteur constant pour chaque sujet (G) est reconnu inter-sujet, et `$effects` reproduit `aov()`.
* **Post-hocs d'un split-plot** : ils passent par `emmeans` sur `lmer(y ~ G*T + (1 | id))`. Les tests appariés sur les cellules croisées apparaient des sujets différents.
* **Variances hétérogènes, données non normales** : la confirmation du test de Levene se fait maintenant par Brown-Forsythe **sur les variances** (`car::leveneTest(center = median)`), et non plus par `onewaytests::bf.test()`, qui compare des moyennes. Une hétéroscédasticité nette n'est plus écartée à tort, et l'ANOVA classique n'est plus retenue dans ce cas.
* **Données non normales, plus de deux groupes, variances hétérogènes confirmées** (asymétrie modérée) : `m.test` applique le F de Welch au lieu de Kruskal-Wallis, qui suppose des variances égales. Sous H0 (n = 5/5/40, sd = 4/4/1), les faux positifs passent de 38,5 % à 10,5 %.
* **Mesures répétées (`Error(id/T)`, 3 temps ou plus)** : la p-value des termes intra-sujet est **toujours corrigée de la sphéricité**, Greenhouse-Geisser si son epsilon est ≤ 0,75, Huynh-Feldt sinon (Girden, 1992). Le test de Mauchly reste affiché, à titre d'information : peu puissant avec peu de sujets, il ne décide plus de la correction. La p-value affichée, celle de `$effects` et la p-value globale sont les p-values corrigées.
  Les p-values des étudiants changent : en général elles augmentent (correction plus prudente), sauf quand F < 1. La correction s'applique au F du tableau `aov()` affiché.
* **Plan en tranches aux groupes inégaux** (`y ~ grp*temps + Error(id/temps)`, 5 + 9 sujets) : l'analyse aboutit (ANOVA à mesures répétées) au lieu de s'arrêter sur « Les variables explicatives de la formule ne correspondent pas aux colonnes de g ».
* Le récit du F de Welch ne le dit plus « robuste » : avec de petits groupes asymétriques, son risque de faux positif reste au-dessus d'alpha.
* **`y ~ T + (1|id)`** : l'analyse et les post-hocs sont appariés, comme avec `Error(id/T)`.
* **Mesures répétées, lignes non triées par sujet** : les données sont triées par sujet, et les post-hocs appariés comparent bien chaque sujet à lui-même.
* **MANOVA non normale** : la voie robuste utilise le lambda de Wilks sur les rangs (`rrcov::Wilks.test(method = "rank")`). L'appel prend moins d'une seconde au lieu de plus de 6 minutes. `$method` vaut `"Robust MANOVA (Wilks, ranks)"`.

### Plantages et refus corrigés

* `verbose = FALSE` sur un modèle mixte : plus d'erreur « $ operator is invalid for atomic vectors ».
* `m.test(formula = y ~ F, data = d)` est accepté.
* `m.test(x, g, id = id)` avec plus de deux modalités ne plante plus.

### Affichage

* Hors locale française, les messages s'affichent en anglais. Jusqu'ici, ils restaient toujours en français.
* En modèle mixte, les lettres d'une interaction nomment la cellule entière (« B t3 ») au lieu du seul premier facteur.
