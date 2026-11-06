# Feuille de route de m.bay

Détail de l'existant : `PERIMETRE.md`. Choix : `DECISIONS.md`. Publication uniquement sur décision d'Antoine.

## Étapes

| # | Étape | État | Dépend de |
|---|---|---|---|
| 0 | Périmètre, comparaison des moteurs, méthodologie, bibliographie | Fait le 5 nov. 2026 (voir réserves du JOURNAL) | |
| 1 | Socle interne : probabilités, inclusion, sensibilité, ligne d'effets | Fait le 5 nov. 2026, testé ; manque le ROPE à plusieurs ddl | |
| 2 | Décisions D1, D2, D6, D9 ; défauts bloquants de la relecture à résoudre (multiplicité recoupée, `integrate` résolu à un facteur) | **Attend Antoine** | |
| 3 | Entrées : normalisation (D9) et refus explicites | À faire | D9 |
| 4 | `m.bay()` 2 groupes et k groupes, validé contre JZS | À faire | D1 |
| 5 | Multi-facteurs, inclusion, `$effects` | À faire | D1, D2 |
| 6 | ANCOVA, pentes hétérogènes | À faire | D1 |
| 7 | Apparié, mesures répétées, split-plot | À faire | D1 |
| 8 | MANOVA 1 facteur puis plusieurs facteurs | À faire | D7 |
| 9 | ROPE, probabilité de direction, taille d'effet, intervalle de crédibilité | À faire | D3 |
| 10 | Post-hoc | À faire | D11 |
| 11 | Bilan lisible, relecture par un étudiant ou industriel type | À faire | 4 à 10 |
| 12 | Voies Stan : pentes aléatoires, robuste, discret | À faire | D1, D6 |
| 3 bis | Variances inégales à un facteur, sans Stan | Calcul vérifié le 6 nov., **loi a priori rejetée** (B1, B2) ; à refaire | D6 |
| 13 | Table de simulation publique, vignette, documentation roxygen | À faire | tout |

## Couverture des plans de m.test()

| Plan | m.test | m.bay prévu | Étalon de comparaison possible |
|---|---|---|---|
| 2 groupes | Student, Welch, Wilcoxon | JZS ; robuste en plus | Intégrale JZS (validée) |
| k groupes | ANOVA, Welch, Kruskal, t1way, med1way | g-prior ; robuste en plus | Intégrale ZS (validée) |
| A*B, n facteurs | aov, Type I/II/III, SRH, t2way, aovp | Inclusion appariée | Double intégrale à 2 g |
| ANCOVA, pentes parallèles | Type III, rlm, HC3 | lmBF, covariable présente sous H0 et H1 | Liang (y~x validé) |
| Pentes hétérogènes | Retrait si p>=alpha | Terme dans l'espace de modèles, pas de retrait | À construire |
| Apparié | t apparié, Wilcoxon | JZS sur différences | Intégrale JZS |
| Mesures répétées | aov+Error, Friedman | Sujet aléatoire | **m.test : voir réserve** |
| Split-plot | Interaction perdue : **anomalie** | Modèle complet | **Pas m.test** : aov+Error direct |
| Mixte à intercept | lmer | whichRandom | lmer pour l'estimation |
| Pentes aléatoires | **Effacées** | Stan | À construire |
| MANOVA 1 et n facteurs | Wilks, Pillai, MCD | NIW analytique | d=1 anovaBF ; Monte-Carlo |
| Discret | Rangs selon route | Logistique, Poisson, cotes proportionnelles, tables | À construire |
| Covariable seule, ANCOVA appariée, MANOVA répétée, MANCOVA | Refusés | D8 | |

## Critères de fin (une voie n'est « validée » que si tous sont exécutés)

1. Comparaison à un calcul indépendant à a priori identiques, écart documenté.
2. Calibration sous l'a priori (moyenne, déciles) et, le cas échéant, SBC, R-hat, ESS.
3. Opérateurs sous effets fixes : faux positifs, puissance, gaussien et non gaussien, à côté de `m.test()`.
4. Sensibilité aux a priori exécutée et affichée.
5. Tests dans `tests/testthat`, sources ASCII, messages bilingues, références en commentaires, BP respectées.
6. Chaque résultat annonce ses a priori ; aucune p-value présentée comme probabilité d'hypothèse.
7. Relecture adverse sans défaut bloquant ouvert, et contrôle de lisibilité du bilan.
