# Périmètre de cor.bay()

Branche `claude/feat-cor-bay-59ccb0`, coupée de `feat/m.bay` (c8e8976). `m.test()`, `m.bay()` et leurs fichiers ne sont pas modifiés : `cor.bay()` réutilise le moteur de récit et de bilan de m.bay (`sys_bay_story.R`, `sys_bay_layout.R`, `sys_bay_bilan.R`) sans y toucher. Les décisions et leurs mesures sont dans `DECISIONS.md` (C1 à C9), les chiffres dans `VALIDATION.md`.

## 1. Ce qui existe dans KefiR

| Fonction | Ce qu'elle fait pour un lien quanti-quanti | Pour cor.bay |
|---|---|---|
| `cor.e()` (`R/cor.e.R`) | Génère une variable corrélée à une autre (simulation pour TP) ; aucun test | sert aux exemples et aux simulations |
| `corrigraph()` (`R/corrigraph.R:159, 233, 512`) | `cor.test()` de Pearson seul, seuil `alpha`, sur toutes les paires ; ni Spearman, ni Kendall, ni contrôle des valeurs éloignées | `cor.bay()` n'y est pas branché : à décider plus tard (corrigraph ferait des dizaines de tests, question de multiplicité) |
| `m.test()` | Pas de voie quanti-quanti | sans objet |
| `stats::cor.test()` | Pearson, Spearman, Kendall ; `alternative`, `conf.level`, formule `~ x + y` | l'interface de référence |

## 2. Interface

| Argument de `cor.test` | Dans `cor.bay` |
|---|---|
| `x`, `y` | identiques ; en plus, noms ou indices de colonnes de `data` (comme `m.bay`) |
| `formula = ~ x + y`, `data` | identiques ; `y ~ x` accepté aussi |
| `method` | `"auto"` (défaut, C3), `"pearson"`, `"kendall"`, `"spearman"` |
| `alternative` | identique ; change les hypothèses (C5) |
| `conf.level` | niveau de l'intervalle de crédibilité (HDI) |
| `exact`, `continuity` | acceptés par `...`, sans effet, avec un message |
| — | `verbose`, `detail`, `prior_h1`, `kappa`, `sensitivity`, `draws`, `seed`, comme `m.bay` |

Retour : objet `corbay` (hypothèses, `prob_link`, facteurs de Bayes, estimation, coefficients observés, méthode et raison, contre-épreuve, sensibilité, contrôles, récit). `print(x)`, `print(x, story = TRUE)`, `print(x, detail = TRUE)` comme `m.bay`.

## 3. Cas couverts

| Cas | Voie | Décision |
|---|---|---|
| Deux variables continues, sans valeur à part | Pearson, `BayesFactor::correlationBF`, a priori bêta étiré (Ly et al. 2016) | C1, C2 |
| Lien monotone non linéaire, rangs voulus | Kendall bayésien (van Doorn et al. 2018), sur demande ou par la règle | C3, C4 |
| Valeurs éloignées, paire influente (distance de Cook) | Kendall, choisi par la règle `auto` ; Pearson calculé en contre-épreuve | C3, C6 |
| Queues lourdes | mesuré : voir C6 (ce que la règle attrape et ce qu'elle laisse) | C6 |
| Échelle courte (notes 1 à 5, ex-aequo nombreux) | Kendall (tau-b) | C3 |
| Spearman | sur demande : Pearson bayésien sur les rangs (Makowski et al. 2020) | C4 |
| Hypothèse unilatérale | deux hypothèses (le signe testé, aucun lien) | C5 |
| Taille de l'échantillon | refus sous 3 paires ; avis sous 10 | C7 |
| Estimation de rho ou tau | médiane a posteriori et HDI (`bayestestR::hdi`) sous l'hypothèse d'un lien | C8 |
| Valeurs manquantes | paires incomplètes retirées, dit | — |

## 4. Hors périmètre, refusé ou non traité

- **Variable non numérique** : refus, avec renvoi vers `m.bay()` (catégories).
- **Variable constante** : refus.
- **Corrélation partielle, plusieurs variables, matrice de corrélation** : non traité (pas de `corrigraph` bayésien ici).
- **Spearman bayésien exact** (van Doorn et al. 2020, augmentation de données) : aucun paquet R établi ne le fournit ; voir C4.
- **Loi t bivariée (modèle robuste paramétrique)** : aucun paquet établi sans Stan ; voir C6.
- **Données appariées dans le temps (autocorrélation)** : non contrôlé ; le bilan ne le dit pas, c'est une limite à connaître.
