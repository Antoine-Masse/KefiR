# Validation de g.bay()

11 novembre 2026. Scripts : `dev/gbay/scripts/sim_validation.R` (sortie `sim_validation.out`, `.rds`) et `dev/gbay/scripts/sim_bilan.R` (`sim_bilan.out`, `.rds`). R 4.6.1, BayesFactor 0.9.12-4.8, DescTools 0.99.60.

**Voie rapide des simulations.** `contingencyTableBF` et `GTest` coûtent environ 22 ms par appel. Les simulations appellent donc les fonctions internes que `contingencyTableBF` appelle elle-même (`BayesFactor:::contingencyJointMultinomial`, `contingencyIndepMultinomial`, `contingencyPoisson`, 0,15 ms) et une statistique G vectorisée. Chaque script vérifie d'abord, sur 20 tables tirées au hasard, l'égalité à 1e-10 avec `contingencyTableBF`, `g.bay()` et `DescTools::GTest`.

**Attention (11 novembre au soir).** Les § 2 (partie indépendance), 3 et 4 (partie indépendance) mesurent l'ancien moteur, `contingencyTableBF` (G3, G10). Ils sont gardés parce qu'ils fondent le remplacement. Le moteur actuel, le test log-linéaire de G11, est mesuré au § 9. Les § 1 (en partie), 2 (ajustement), 5 et 6 valent toujours : l'ajustement et l'estimation n'ont pas changé.

Sauf mention, « oui » veut dire probabilité a posteriori d'un lien (ou d'un écart) supérieure au seuil, à chances de départ égales ; « bilan oui » veut dire ce que dit vraiment le bilan (§ 3).

## 1. Accord numérique (tests testthat, `tests/testthat/test-gbay.R`)

- (Ancien moteur.) `g.bay()` = `contingencyTableBF` à 1e-10 sous les cinq plans. Remplacé par : log BF = (G − ddl·log N)/2 avec G de `GTest`, à 1e-8, sur quatre tables ; = −BIC/2 des `glm` logistiques à une ligne par individu (méthode de m.bay, D25) ; table, `x, y`, `y ~ g` et table transposée donnent le même nombre.
- Ajustement : forme close = intégration numérique à deux modalités (trois cas, a = 1/2, 1, 2) à 1e-7 ; = `multibridge::mult_bf_equality` pour des proportions égales à 1e-7 (test ignoré si `multibridge` n'est pas installé ; il l'a été ici dans une bibliothèque temporaire : log BF 0,38202 des deux côtés sur (3, 7, 12)).
- (Ancien moteur.) `g.bay(y ~ g)` = `m.bay(y ~ g)` à 1e-10 ; plus vrai tant que m.bay garde `contingencyTableBF`.
- Formes d'entrée : table, `x` et `y`, `~ A + B`, noms de colonnes, `Freq ~ A + B` donnent le même BF ; ajustement : vecteur, facteur, `~ A`, matrice à une ligne, `p` à normaliser.
- 2 × 2 à N = 4 000 : log du rapport des cotes à moins de 0,03 de `glm` binomial ; V de Cramér à moins de 0,01 de l'observé, dans son intervalle.
- Cases à N = 5 000 : rapport a posteriori à moins de 0,03 de O/E.

## 2. Taux de « oui » sous H0 (2 000 jeux par ligne, 1 000 à N = 1 000)

### Ajustement

P > 0,75 au plus : 6,4 % (5 modalités égales, N = 10) ; 3,4 % à N = 20 ; ≤ 1,4 % dès N = 50. P > 0,95 : 0,7 % au plus. Le test G (p < 0,05) : 3,7 à 11,5 % (11,5 % à N = 10, deux modalités, où il est trop libéral). La probabilité de l'ajustement dépasse 0,75 dans 81 à 100 % des jeux dès N = 50 : sous H0, g.bay soutient l'ajustement, ce que le test G ne peut pas faire.

| k | p | N = 10 | 20 | 50 | 200 | 1 000 |
|---|---|---|---|---|---|---|
| 2 | 0,5/0,5 | 1,8 / 0,1 | 4,2 / 0,2 | 1,4 / 0,3 | 0,8 / 0 | 0,2 / 0,1 |
| 3 | 1/3 chacune | 5,7 / 0,7 | 2,9 / 0,4 | 0,8 / 0,1 | 0,3 / 0 | 0,2 / 0,1 |
| 3 | 0,6/0,3/0,1 | 2,7 / 0,3 | 1,6 / 0,2 | 0,6 / 0,2 | 0,4 / 0,1 | 0 / 0 |
| 5 | 0,2 chacune | 6,4 / 0,6 | 3,4 / 0,4 | 1,3 / 0,3 | 0,2 / 0 | 0 / 0 |
| 5 | 0,4/0,25/0,15/0,1/0,1 | 3,3 / 0,5 | 2,6 / 0,4 | 0,8 / 0,1 | 0,1 / 0 | 0 / 0 |

(% de P > 0,75 / % de P > 0,95.)

### Indépendance (plan joint analysé en joint ; plan à lignes fixées analysé en indépendant et, à tort, en joint)

| Table | N | joint/joint | lignes/indép. | lignes/joint | test G |
|---|---|---|---|---|---|
| 2 × 2 égales | 20 | 10,5 / 1,3 | 4,4 / 0,3 | 11,7 / 1,0 | 5,4-6,9 |
| 2 × 2 égales | 50 | 5,9 / 0,5 | 2,4 / 0,4 | 6,2 / 0,5 | 6,0 |
| 2 × 2 égales | 200 | 2,7 / 0,3 | 1,8 / 0,1 | 2,8 / 0,2 | 5,0-5,5 |
| 2 × 2 inégales | 50 | 3,4 / 0,4 | 3,1 / 0,3 | 4,1 / 0,4 | 6,0-6,6 |
| 2 × 3 | 20 | 12,7 / 1,4 | 8,3 / 0,9 | 14,3 / 1,7 | 8,5-9,2 |
| 2 × 3 | 50 | 7,7 / 1,3 | 4,6 / 0,8 | 9,2 / 1,3 | 6,1-6,4 |
| 2 × 3 | 200 | 2,2 / 0,3 | 1,7 / 0,1 | 2,6 / 0,2 | 5,0-5,5 |
| 3 × 4 | 10 | 35,5 / 2,2 | 19,8 / 0,3 | 57,2 / 4,1 | 6,7-9,1 |
| 3 × 4 | 20 | 30,2 / 4,8 | 14,4 / 1,7 | 38,9 / 6,8 | 9,2-10,2 |
| 3 × 4 | 50 | 13,6 / 3,3 | 5,2 / 0,8 | 15,8 / 3,4 | 7,0 |
| 3 × 4 | 200 | 1,4 / 0,3 | 0,6 / 0 | 1,6 / 0,4 | 5,0-5,2 |

(% de P > 0,75 / % de P > 0,95 ; table complète dans `sim_validation.out`.)

**Lecture.**
- Dès N = 200, P > 0,75 dans 0,6 à 2,8 % des cas, P > 0,95 dans 0,4 % au plus : plus prudent que le test G à 5 %.
- Dans les tables creuses, le facteur de Bayes brut penche nettement vers un lien (jusqu'à 57 % de P > 0,75). C'est le défaut que traite G10.
- Le plan à marge fixée analysé en joint (ce que fait l'aiguillage automatique pour une table ou `x, y`) est plus libéral que le bon plan ; la comparaison des plans à l'étape 5 en tient compte (§ 3).

## 3. Ce que dit le bilan sous H0, avec et sans la prudence des tables creuses (`sim_bilan.R`, 1 000 jeux par ligne)

« Bilan oui » : P ≥ 0,75, facteur de Bayes ≥ 3, et l'ordre ne change pas quand a (1, 2, 4), les chances de départ (1:2 à 2:1) et le plan (joint, Poisson, indépendant lignes ou colonnes) varient.

| Table | N | attendu moyen | P > 0,75 brut | bilan oui (dont très probable), avant G10 | part de tables creuses | bilan oui après G10 (joint ; lignes fixées) |
|---|---|---|---|---|---|---|
| 2 × 2 | 10 | 2,5 | 8,2 | 2,4 (0,3) | 100 | 0 ; 0 |
| 2 × 2 | 20 | 5 | 9,9 | 5,7 (1,1) | 97 | 0,2 ; 0,3 |
| 2 × 2 | 50 | 12,5 | 6,6 | 6,6 (0,7) | 0,2 | 6,6 ; 3,2 |
| 2 × 2 | 200 | 50 | 2,2 | 2,2 (0,4) | 0 | 2,2 ; 1,7 |
| 2 × 3 | 20 | 3,3 | 12,3 | 5,2 (1,6) | 100 | 0 ; 0 |
| 2 × 3 | 50 | 8,3 | 8,6 | 6,4 (1,2) | 45 | 4,1 ; 2,0 |
| 2 × 3 | 100 | 16,7 | 5,2 | 4,1 (0,9) | 0 | 4,1 ; 2,0 |
| 3 × 4 | 20 | 1,7 | 31,8 | 7,3 (5,2) | 100 | 0 ; 0 |
| 3 × 4 | 50 | 4,2 | 13,3 | 6,5 (3,1) | 100 | 0 ; 0 |
| 3 × 4 | 100 | 8,3 | 4,9 | 1,9 (1,0) | 80 | 0,4 ; 0 |
| 3 × 4 | 200 | 16,7 | 1,5 | 0,7 (0) | 3 | 0,7 ; 0,8 |
| 4 × 5 | 20 | 1 | 66,8 | 6,2 (6,1) | 100 | 0 ; 0 |
| 4 × 5 | 50 | 2,5 | 41,6 | 11,9 (11,6) | 100 | 0 ; 0 |
| 4 × 5 | 100 | 5 | 13,8 | 3,7 (3,7) | 100 | 0 ; 0,2 |
| 4 × 5 | 200 | 10 | 2,2 | 0,3 (0,3) | 0 | 0,3 ; 0,3 |

**Lecture.** La sensibilité (a, chances, plans) retire déjà la plupart des faux « oui » bruts. Elle ne suffit pas dans les tables creuses, où les faux « oui » restent « très probables » (11,6 % en 4 × 5, N = 50). Après G10, le bilan dit « un lien existe » dans 6,6 % des cas au plus (2 × 2 à N = 50, plan joint), comparable au test G à 5 %, et dans moins de 1 % des cas en 3 × 4 et 4 × 5.

## 4. Puissance

### Comparée au test G (`sim_validation.R`, 1 000 jeux par ligne, P > 0,75 / P > 0,95 / test G p < 0,05)

| Cas | w | N = 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| indépendance 2 × 3 | 0,1 | 14 / 3 / 10 | 12 / 3 / 13 | 17 / 5 / 25 | 23 / 11 / 41 |
| indépendance 2 × 3 | 0,2 | 33 / 10 / 24 | 43 / 19 / 45 | 62 / 38 / 73 | 91 / 75 / 96 |
| indépendance 2 × 3 | 0,3 | 58 / 29 / 51 | 77 / 53 / 78 | 96 / 88 / 98 | 100 / 100 / 100 |
| ajustement k = 4 | 0,1 | 2 / 0 / 7 | 1 / 0 / 11 | 2 / 1 / 18 | 3 / 1 / 37 |
| ajustement k = 4 | 0,2 | 8 / 2 / 20 | 14 / 6 / 40 | 25 / 13 / 69 | 60 / 40 / 96 |
| ajustement k = 4 | 0,3 | 23 / 9 / 41 | 40 / 24 / 73 | 75 / 58 / 97 | 99 / 97 / 100 |

**Lecture.**
- Indépendance : à P > 0,75, la puissance suit celle du test G (un peu au-dessus à N = 50, un peu en dessous à N ≥ 200) ; à P > 0,95, elle est plus faible, comme attendu d'un seuil plus exigeant.
- Ajustement : nettement moins puissant que le test G. Contre une valeur ponctuelle, l'alternative uniforme de Dirichlet étale sa masse sur toutes les proportions possibles et favorise H0 face aux petits écarts (effet de Lindley ; Kass et Raftery 1995, § 5). Un écart w = 0,2 (« moyen » chez Cohen) demande environ N = 400 pour être dit avec P > 0,75 dans plus d'un cas sur deux. C'est le prix de taux de faux « oui » quasi nuls (§ 2) ; la page d'aide et le récit le disent par le facteur de Bayes, qui reste près de 1 (« non tranché ») quand les données sont trop peu nombreuses.

### Ce que dit le bilan sous un lien, avant et après G10 (`sim_bilan.R`, 500 jeux par ligne, % de « bilan oui »)

| Cas | N = 20 | 50 | 100 | 200 |
|---|---|---|---|---|
| 2 × 3, lien fort | 17 → 0 | 29 → 17 | 44 → 44 | 72 → 72 |
| 2 × 2, lien moyen | 14 → 0 | 29 → 29 | 48 → 48 | 75 → 75 |
| 3 × 4, lien moyen | 38 → 0 | 81 → 0 | 97 → 30 | 100 → 98 |

Dans les tables creuses, le bilan dit « plutôt un lien, sans conclure », avec la probabilité ; c'est le coût de G10.

## 5. Cases signalées sous H0 (G6, 400 jeux par ligne, plan joint)

| Table | N | tables avec au moins une case signalée | cases signalées |
|---|---|---|---|
| 2 × 2 | 50 | 6,0 % | 6,0 % |
| 2 × 2 | 200 | 6,8 % | 6,8 % |
| 2 × 3 | 50 | 6,0 % | 3,0 % |
| 2 × 3 | 200 | 8,0 % | 3,1 % |
| 3 × 4 | 50 | 26,3 % | 3,0 % |
| 3 × 4 | 200 | 35,0 % | 4,3 % |

**Lecture.** Une case isolée est signalée à tort dans 3 à 7 % des cas, ce qu'on attend d'un intervalle à 95 % lu case par case ; en 2 × 2 les quatre cases vont ensemble. Dans une grande table, une case au moins l'est souvent (26 à 35 % en 3 × 4). C'est pourquoi les cases ne sont lues dans le bilan que si un lien est établi (G6) ; quand l'indépendance est en tête, l'étape 6 dit « à titre indicatif » et le bilan ne les lit pas.

## 6. Couverture de l'intervalle à 95 % du rapport des cotes (2 × 2, OR vrai 3,5, 400 jeux)

N = 40 : 95,8 % ; N = 100 : 94,0 % ; N = 400 : 93,0 %.

## 7. Lisibilité (tests testthat)

Sept sorties (ajustement, 2 × 2, 2 × 3, 3 × 3, `y ~ g`, petits effectifs, `detail = TRUE`), en français et en anglais, avec `width` à 80 et 120 : aucune ligne au-delà de 80 colonnes (tabulation comptée 8), chaque ligne porte une marque (« ==> », « --> », « Note : », titre numéroté ou suite repliée), aucune p-value. À `width = 60`, rien au-delà de 60. `print(r, story = TRUE)` rejoue exactement le récit.

## 8. R CMD check (11 novembre, sur 1ecbda5)

`R CMD check --no-manual --ignore-vignettes`, suite de tests complète du paquet comprise : **Status: 1 NOTE**. La NOTE (« Problems with news in 'NEWS.md': No news entries found ») existe déjà sur c8e8976 : l'en-tête `# KefiR (development version, branch feat/m.bay)` n'est pas lu comme une version par R (vérifié avec `tools:::.build_news_db_from_package_NEWS_md` sur les deux versions du fichier). Exemples OK, dépendances des tests OK (le test `multibridge` passe par `getExportedValue`, ignoré si le paquet manque), code ASCII OK. Information sur les xrefs : `ez` absent de mbay-lib, comme pour m.bay.

## 9. Recalibration après G11 (test log-linéaire, BIC avec n = N), 11 novembre au soir

Script : `scripts/sim_loglin.R` (sortie `sim_loglin.out`, `.rds`). 2 000 tables par scénario (1 000 pour la puissance), plan joint. En tête du script, contrôle à 1e-8 : G vectorisé = `GTest` et log BF = `.gbay_lbf_ct()` (glm). « Bilan oui » : P ≥ 0,75, BF ≥ 3, et l'ordre tient de 1:2 à 2:1.

### Sous H0 : G11 contre l'ancien moteur (`contingencyTableBF`, a = 1)

| Table | Marges | N | G11 BF > 3 | **G11 BF > 20** | G11 bilan oui (très probable) | ancien BF > 20 : joint / Poisson / lignes fixées | test G p < 0,05 |
|---|---|---|---|---|---|---|---|
| 2 × 2 | égales | 30 | 1,95 % | 0,20 % | 1,95 % (0,20) | 0,55 / 0,80 / 0,40 % | 5,0 % |
| 2 × 2 | égales | 60 | 1,15 | 0,15 | 1,15 (0,15) | 0,40 / 0,55 / 0,25 | 4,9 |
| 2 × 2 | égales | 200 | 0,80 | 0,10 | 0,80 (0,10) | 0,55 / 0,60 / 0,40 | 5,1 |
| 2 × 2 | inégales | 30 | 1,65 | 0,20 | 1,65 (0,20) | 0,30 / 0,40 / 0,25 | 6,7 |
| 2 × 2 | inégales | 60 | 1,60 | 0,20 | 1,60 (0,20) | 0,25 / 0,40 / 0,15 | 6,1 |
| 2 × 2 | inégales | 200 | 0,50 | 0,05 | 0,50 (0,10) | 0,25 / 0,30 / 0,20 | 4,4 |
| 3 × 4 | égales | 30 | 0,30 | 0 | 0,30 (0) | 8,0 / 30,4 / 2,9 | 11,3 |
| 3 × 4 | égales | 60 | 0 | 0 | 0 (0) | 4,3 / 17,4 / 1,3 | 7,1 |
| 3 × 4 | égales | 200 | 0,05 | 0 | 0,05 (0) | 0,3 / 2,5 / 0,05 | 5,9 |
| 3 × 4 | inégales | 30 | 0 | 0 | 0 (0) | 2,1 / 10,4 / 1,1 | 7,0 |
| 3 × 4 | inégales | 60 | 0,10 | 0,05 | 0,10 (0,05) | 1,2 / 6,1 / 0,7 | 6,9 |
| 3 × 4 | inégales | 200 | 0 | 0 | 0 (0) | 0,05 / 1,2 / 0,05 | 7,0 |
| 4 × 5 | égales | 30 | 0 | 0 | 0 (0) | 21,1 / 90,4 / 5,6 | 12,9 |
| 4 × 5 | égales | 60 | 0 | 0 | 0 (0) | 11,6 / 76,9 / 2,7 | 10,4 |
| 4 × 5 | égales | 200 | 0 | 0 | 0 (0) | 0,9 / 13,2 / 0,1 | 6,4 |
| 4 × 5 | inégales | 30 | 0 | 0 | 0 (0) | 7,5 / 61,5 / 2,6 | 9,0 |
| 4 × 5 | inégales | 60 | 0 | 0 | 0 (0) | 2,4 / 42,8 / 0,7 | 8,1 |
| 4 × 5 | inégales | 200 | 0 | 0 | 0 (0) | 0,25 / 6,2 / 0,15 | 7,7 |

**Lecture.**
- La cible de la session parente, moins de 1 % de BF > 20 sous H0, est tenue partout : 0,2 % au plus, contre 90 % pour l'ancien moteur sous le plan Poisson (4 × 5, N = 30).
- Le bilan dit « un lien existe » sous indépendance dans 2 % des tables au plus (2 × 2, N = 30), « très probablement » dans 0,2 % au plus.
- Nos chiffres de l'ancien moteur rejoignent ceux de l'audit : 4 × 5 à marges égales, N = 60 : 11,6 % (joint) et 76,9 % (Poisson), contre 11,0 % et 82,7 % mesurés par la session parente à N = 50.
- Plus de dépendance au plan : G, ddl et N sont les mêmes sous tous les plans.

### Puissance (1 000 tables par ligne ; écart de taille w, w de Cohen)

| Table | w | N | G11 P > 0,75 | G11 P > 0,95 | ancien (joint) P > 0,75 | test G p < 0,05 |
|---|---|---|---|---|---|---|
| 2 × 2 | 0,3 | 30 | 25,8 % | 9,8 % | 44,2 % | 40,0 % |
| 2 × 2 | 0,3 | 60 | 46,0 | 20,2 | 65,9 | 66,4 |
| 2 × 2 | 0,3 | 200 | 93,6 | 82,5 | 97,8 | 99,1 |
| 2 × 2 | 0,5 | 30 | 68,0 | 39,4 | 83,1 | 80,6 |
| 2 × 3 | 0,3 | 60 | 27,1 | 12,3 | 63,9 | 56,7 |
| 2 × 3 | 0,3 | 200 | 82,7 | 63,3 | 96,5 | 98,2 |
| 2 × 3 | 0,5 | 60 | 86,5 | 69,6 | 97,4 | 96,8 |
| 3 × 4 | 0,3 | 60 | 2,4 | 0,5 | 54,6 | 45,2 |
| 3 × 4 | 0,3 | 200 | 15,7 | 9,4 | 83,6 | 93,1 |
| 3 × 4 | 0,1 | 200 | 0 | 0 | 5,1 | 14,8 |

(Table complète dans `sim_loglin.out`.)

**Lecture.**
- En 2 × 2, un lien moyen (w = 0,3) est dit à P > 0,75 dans 94 % des tables à N = 200, contre 99 % au test G. À N = 30-60, la puissance est environ deux tiers de celle du test G.
- Dans les grandes tables, la pénalité ddl·log N rend le BIC très prudent : en 3 × 4, un lien moyen n'est presque jamais dit à N = 60 (2,4 %), et rarement à N = 200 (15,7 %, contre 93 % au test G). Le bilan dit alors « indécis » ou « plutôt pas de lien » ; le V de Cramér et les cases (étape 6) restent lisibles. C'est le prix d'un taux de faux « oui » quasi nul (Raftery 1995 : le BIC favorise les modèles simples plus que les p-values).
- Suite possible, à décider à part : un a priori moins diffus que l'information unitaire sur l'interaction rendrait de la puissance aux grandes tables ; il faudrait le calibrer de même.

## 10. Choix du moteur d'indépendance (G12), 12 novembre

Scripts :
- `scripts/sim_engines.R` : BIC, BAS robuste, BAS hyper-g/n et marges neutres (a = 1). Sorties `sim_engines_h0.out` sous H0 ; `sim_engines_h1.out` sous un lien, sans BAS, dont `bas.glm` a planté pendant ces scénarios après avoir échoué sous H0.
- `scripts/sim_neutral.R` : marges neutres à a = 1, 2, 3, avec le garde-fou ; sortie `sim_neutral.out`.

1 000 tables par scénario. Lien diagonal de V de Cramér donné, marges égales.

### Sous H0 : P(lien) > 0,95, valeur la plus haute sur les 24 à 26 scénarios (2 × 2 à 4 × 5, N = 30, 60, 200, marges égales et inégales, lignes fixées)

| Moteur | Maximum | Où |
|---|---|---|
| BIC (G11) | 0,4 % | 2 × 2 |
| BAS robuste, n = N | 98,4 % | 4 × 5, N = 30 |
| BAS hyper-g/n, n = N | 83,0 % | 4 × 5, N = 30 |
| marges neutres, a = 1 | 1,1 % | 3 × 4, N = 30 |
| marges neutres, a = 2 | 0,7 % | 4 × 5, N = 60 |
| **marges neutres, a = 3** | **0,6 %** | 2 × 2 et 4 × 5 |

### Sous un lien modéré : indépendance « forte » (P(indépendance) ≥ 0,95), maximum sur les scénarios

| Moteur | V = 0,2 | V = 0,3 |
|---|---|---|
| BIC | 88,8 % (3 × 4, N = 50) | 66,3 % (3 × 4, N = 50) |
| marges neutres, a = 1 | 35,7 % (21,2 % avec le garde-fou) | 8,2 % |
| marges neutres, a = 2 | 8,9 % (7,4 %) | 0,8 % |
| **marges neutres, a = 3** | **2,9 % (2,3 %)** | **0 %** |

### Puissance de a = 3 (P > 0,75) et du test G (p < 0,05)

| Table | V | N = 50 | 100 | 200 |
|---|---|---|---|---|
| 2 × 2 | 0,2 | 24 % / 33 % | 40 / 53 | 71 / 81 |
| 3 × 3 | 0,2 | 29 / 34 | 51 / 58 | 78 / 90 |
| 3 × 4 | 0,2 | 28 / 31 | 43 / 52 | 70 / 85 |
| 2 × 2 | 0,3 | 49 / 59 | 76 / 87 | 98 / 99 |
| 3 × 3 | 0,3 | 64 / 63 | 88 / 93 | 99 / 100 |
| 3 × 4 | 0,3 | 56 / 60 | 85 / 88 | 99 / 100 |

**Lecture.** Les marges neutres à a = 3 sont le seul moteur qui tienne les deux critères à la fois : moins de 1 % de faux lien fort sous H0, au plus 5 % d'indépendance forte sous un lien modéré. Leur puissance vaut 73 à 100 % de celle du test G. La table de l'audit (12 11 6 / 9 10 14 / 7 12 19) passe de « indépendantes : 99,6 % (soutien très fort) » (BIC) à « indécis » (indépendance 53,9 %, BF 1,17, sensibilité 31 % à 80 %).

Les § 9 (BIC) et 2 à 4 (ancien moteur, partie indépendance) ne décrivent plus le moteur actuel. L'ajustement (§ 2) et l'estimation (§ 5, 6) n'ont pas changé.
