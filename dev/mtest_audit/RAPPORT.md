# Audit de `m.test()` — novembre 2026

Branche `claude/fix-m-test-audit-c7fc21`, partie de master (3675fc7, version 0.0.2.2).
Rien n'est poussé. Ni master, ni `m.bay`, ni les branches `feat/m.bay`, cor.bay et g.bay ne sont touchés.

**Méthode.** Chaque défaut est reproduit sur des données construites pour l'occasion, puis comparé à une fonction R de référence (`aov`, `car::leveneTest`, `pairwise.t.test`, `rrcov::Wilks.test`…). Chaque correction a son commit et son test dans `tests/testthat/test-m.test-audit.R`. Ce test échoue sur le code d'avant et passe après. Les taux de faux positifs viennent de simulations sous H0 (script `sim_h0.R`, section 3).

**Mesure du temps.** Un appel `m.test()` sur 50 observations coûte environ 1 s avec `boot = FALSE` et 47 s avec la valeur par défaut `boot = TRUE`, même quand `return = FALSE` ne demande qu'une p-value (machine partagée, 8 cœurs). Les simulations utilisent donc `boot = FALSE`. La section 3 vérifie que `boot` ne change pas la p-value globale.

## 1. Défauts signalés par les audits de m.bay (a à g)

| | Défaut | Reproduit ? | Avant | Après | Commit |
|---|---|---|---|---|---|
| a | Split-plot `y ~ G*T + Error(id/T)` : interaction retirée sans message | Oui | `aov` : G:T p = 0,030 ; `$effects` : G et T seulement (lmer `y ~ G + T + (1\|id)`) | `$effects` = G, T, G:T, p-values identiques à `aov` (écart < 1e-6) | a9ddfcb |
| b | Plantage avec `verbose = FALSE` | Oui | « $ operator is invalid for atomic vectors » | pas d'erreur | 573aa2c |
| c | `bf.test` (Brown-Forsythe sur les **moyennes**) pour confirmer Levene | Oui | Levene p = 0,00017 puis `bf.test` p = 0,93 : variances « homogènes », ANOVA classique | Brown-Forsythe sur les variances (`leveneTest(center = median)`) | 50e4411 |
| d | Post-hocs non appariés avec `(1\|id)` | Oui | t1-t2 : p = 1 (non apparié, `pool.sd`) | p = 0,024, identique à `Error(id/T)` | d3e1c55 |
| e | `m.test(formula = y ~ F, data = d)` refusé | Oui | « doit être une formule valide » | accepté, même p-value que `m.test(y ~ F, d)` | 17a0916 |
| f | Message anglais jamais affiché | Oui, en locale `C` | toute la sortie en français | sortie en anglais hors locale française | 2575f1b |
| g | MANOVA non normale : plus de 6 minutes | Oui, cause établie | `Wilks.test(method = "mcd")` : 3000 tirages de 0,13 s | Wilks sur les rangs, moins d'une seconde | 0e2765f |

### a. Interaction perdue en split-plot

Données : 12 sujets, G inter-sujet (2 niveaux), T intra-sujet (3 temps), effet seulement en B×t3. `summary(aov(y ~ G*T + Error(id/T)))` : G p = 0,199 ; T p = 0,0096 ; G:T p = 0,030.

Quatre causes en chaîne :
1. Faute de `within`, tous les facteurs étaient pris pour intra-sujet, G compris. Chaque sujet paraissait alors « incomplet » (6 conditions attendues, 3 observées).
2. Le plan partait donc vers `lmer`, dont la partie fixe était reconstruite en `within + between` : `y ~ G + T + (1|id)`, sans G:T.
3. Le contrôle « équilibrage RM strict » comptait les sujets de l'autre niveau du between comme incomplets : `table()` d'un facteur `id` garde tous ses niveaux.
4. Le contrôle de croisement attendait within × between cellules par sujet. Un sujet n'appartient pourtant qu'à un niveau du between.

Correction : un facteur constant pour chaque sujet est déclaré inter-sujet, et le message le dit. La partie fixe de `lmer` reprend les termes de la formule. Les deux contrôles raisonnent sujet par sujet. Le plan équilibré part maintenant vers `aov()` avec `Error()`, et `$effects` reproduit `aov` exactement. Un plan déséquilibré part vers `lmer`, qui garde G:T.

### b. Plantage avec `verbose = FALSE`

`valreg()` renvoie une liste `(valid, k)` seulement si `k != 0`. Avec `verbose = FALSE`, `.vbse()` n'incrémente pas `k`, et `valreg()` renvoie un simple logique : `valreg_output$valid` échouait. `.mixed_model_analysis()` accepte maintenant les deux formes. L'interface de `valreg()` n'est pas modifiée.

### c. Brown-Forsythe sur les moyennes

`car::leveneTest()` centre par défaut sur la **médiane** : c'est déjà le test de Brown-Forsythe. La « confirmation » appelait ensuite `onewaytests::bf.test()`, qui est le test de Brown-Forsythe sur les **moyennes** (une ANOVA de Welch modifiée). Sous H0 (moyennes égales), il ne confirme presque jamais, si bien que l'hétéroscédasticité était écartée avec le message « Levene probablement faussé par asymétrie/outliers ».

Correction : Levene centré sur la moyenne (`center = mean`), puis confirmation par Brown-Forsythe sur les variances (`center = median`), comme l'annoncent les messages et la sortie `code = TRUE`. Mesure sous H0 : section 3.

### d. `(1|id)` non apparié

L'`id` extrait de `y ~ T + (1|id)` était posé **après** le bloc qui déduit `paired = TRUE` de la présence d'`id`. La sortie annonçait « traitement comme mesures répétées », mais l'analyse et les post-hocs étaient ceux de groupes indépendants. Sur les mêmes données : t1-t2 p = 1 (avant), contre 0,024 avec `Error(id/T)` et avec un `pairwise.t.test(paired = TRUE)` corrigé par Holm.

### e. `formula =` nommé

`formula <- call_expr$formula` remplaçait l'argument évalué par l'expression non évaluée (classe `call`). La ligne est retirée. Une formule passée en texte (`"y ~ F"`) est convertie.

### f. Messages anglais

`.vbse()`, `.exit()` et `.dbg()` remplaçaient `ang` par `fr` dès que `length(ang) == 1`, c'est-à-dire pour tout message ordinaire. Le paquet parlait donc français dans toutes les locales. Vérifié en `LC_CTYPE = "C"` : tout était en français avant, tout est en anglais après. Le français ne sert plus que si l'anglais manque (`NULL`, vide ou `NA`). En locale française, rien ne change pour les étudiants.

### g. MANOVA non normale

Cause établie : la voie robuste appelait `rrcov::Wilks.test(method = "mcd")`. Celle-ci estime la loi de la statistique par simulation, avec `nrep = 3000` tirages par défaut. Mesures (3 réponses, 60 observations) : 7,6 s pour 50 tirages, 25 s pour 200, soit environ 6 min 20 s pour 3000. La p-value dépendait en outre du nombre de tirages (0,37 à 50, 0,29 à 200).

Décision : le lambda de Wilks sur les rangs (Nath et Pavur, 1985), `Wilks.test(method = "rank")` du même paquet, sans tirage. Calibration sous H0 non normale (réponses exponentielle, exponentielle au carré, log-normale ; N = 2000 par plan) :

| Effectifs | Wilks rangs | Pillai classique |
|---|---|---|
| 20/20/20 | 5,5 % | 3,6 % |
| 8/8/8 | 5,3 % | 3,3 % |
| 5/10/30 | 4,9 % | 5,6 % |

`$method` devient `"Robust MANOVA (Wilks, ranks)"` (changement visible : NEWS).

## 2. Défauts trouvés pendant l'audit

| | Défaut | Avant | Après | Commit |
|---|---|---|---|---|
| h | Mesures répétées, lignes non triées par sujet : post-hocs appariés entre sujets différents | t1-t2 p = 1 (données mélangées) au lieu de 0,024 | identique aux données triées | 204a39f |
| h' | `m.test(x, g, id = id)` avec plus de 2 modalités | plantage dans `terms(NULL)` | formule `x ~ g + Error(id/g)` construite | 204a39f |
| i | Split-plot : post-hocs par tests appariés entre cellules de sujets différents (A.t1 contre B.t3) | toutes les comparaisons à p ≈ 1 malgré G:T p = 0,030 | `emmeans` sur `lmer(y ~ G*T + (1\|id))` : B t1-B t3 p = 0,018 (Tukey), identique au calcul direct | 9c49338 |
| j | Lettres d'une interaction en modèle mixte : seul le premier facteur est nommé | « A : EMM = 9.165 » trois fois | « A t1 : EMM = 9.165 » | edd953f |
| k | Non normal, k > 2, hétéroscédasticité confirmée : Kruskal-Wallis « valide même avec hétéroscédasticité si formes similaires » | Kruskal-Wallis choisi 114 fois sur 200, 21 % de rejets sous H0 (calcul direct) | F de Welch (asymétrie modérée, aucun groupe constant) | 88f493d, dc74ae7 |

### h. Appariement par position

`pairwise.t.test(paired = TRUE)` apparie la i-ème observation de chaque groupe. Un commentaire de `m.test.R` (« NE PAS TRIER ») gardait l'ordre des lignes. Sur un plan à mesures répétées dont les lignes sont mélangées, les post-hocs comparaient donc des sujets différents, et la p-value globale (`aov` avec `Error()`) restait juste. Le tri par sujet est stable : il ne change pas l'ordre à l'intérieur d'un sujet.

### i. Post-hocs du split-plot

Une fois le split-plot reconnu (a), l'analyse globale passe par `aov()` avec `Error()`. Les post-hocs appelaient alors `pairwise.t.test(paired = TRUE)` sur les six cellules G×T. Ce test apparie par position la i-ème observation de A.t1 et celle de B.t3, c'est-à-dire deux sujets différents. Toutes les comparaisons sortaient à p ≈ 1. Ces plans passent maintenant par `emmeans` sur `lmer(y ~ G + T + G:T + (1 | id))`, comme la voie « modèle mixte » qui existait déjà, et la sortie dit pourquoi. Pour ce plan équilibré, le modèle reproduit exactement les F d'`aov()` (G p = 0,199 ; T p = 0,0096 ; G:T p = 0,030).

### k. Kruskal-Wallis sous hétéroscédasticité

Une fois c corrigé, l'hétéroscédasticité est bien détectée, mais le « diagnostic pragmatique » envoyait alors vers Kruskal-Wallis, avec le message « valide même avec hétéroscédasticité si formes similaires ». C'est faux : Kruskal-Wallis teste l'égalité des distributions, et des variances inégales suffisent à le faire rejeter. Calibration directe (N = 1000, n = 5/5/40, sd = 4/4/1) :

| Données | ANOVA classique | Kruskal-Wallis | `t1way` | F de Welch |
|---|---|---|---|---|
| normales | 50,4 % | 21,0 % | 9,6 % | 5,1 % |
| exponentielles centrées | 50,4 % | 33,3 % | 21,5 % | 17,2 % |

Décision : quand l'hétéroscédasticité est confirmée (Levene puis Brown-Forsythe) sur des données non normales, avec k > 2, une asymétrie modérée (max |skew| < 2), sans outlier extrême ni groupe constant, `m.test` applique le F de Welch, suivi de post-hocs t sans variance commune (comme la branche normale hétéroscédastique). C'est le mieux calibré des quatre, et il est recommandé par Delacre et al. (2019). Les autres cas gardent l'arbre existant (`med1way`, `t1way`, Kruskal-Wallis). Le premier essai échouait avec un groupe constant (test existant « un groupe constant ne fait pas planter l'analyse ») : la bascule exige maintenant des variances toutes positives (dc74ae7).

## 3. Taux de faux positifs sous H0

Script `sim_h0.R` (copie dans ce dossier, résumé par `summ.R`) : 200 jeux par scénario, `verbose = FALSE`, `boot = FALSE`, même graine par jeu avant et après correction. L'écart-type Monte-Carlo est d'environ 1,5 point autour de 5 %. « Avant » = master (3675fc7) ; « après » = 88f493d (dc74ae7 ne change que le cas d'un groupe constant).

| Scénario | Plan | Avant | Après | Commentaire |
|---|---|---|---|---|
| norm3 | 3 × 10, normal | 3,0 % | 3,0 % | |
| skew3 | 3 × 15, exponentiel | 3,0 % | 3,0 % | |
| lnorm4 | 4 × 12, log-normal | 3,0 % | 3,0 % | |
| **hetero_c** | n = 5/5/40, sd = 4/4/1, normal | **38,5 %** | **10,5 %** | c, puis k |
| **hetero_c2** | idem, exponentiel centré | **43,0 %** | **22,5 %** | c, puis k ; reste élevé (§ 4.1) |
| norm2 | 2 × 12, normal | 3,5 % | 3,5 % | |
| skew2 | 2 × 15, exponentiel | 4,5 % | 4,5 % | |
| hetero2 | n = 6/30, sd = 3/1 | 8,5 % | 8,5 % | stratégie en deux temps (§ 4.1) |
| paired2 | 15 sujets × 2 | 5,5 % | 5,5 % | |
| rm3 | 15 sujets × 3, `Error(id/T)` | 4,5 % | 4,5 % | post-hocs : 2 % |
| rm3_lmer | idem, `(1\|id)` | 8,5 % | 6,9 % (terme T) | identique à `aov` sur 25 jeux sur 25 : l'écart est du bruit |
| split | 12 sujets, G × T | **plantage 200/200** (b) | G 4,5 %, T 4,5 %, G:T 6,0 % | |
| twoway | 2 × 3 × 8, normal | A 2,5 %, B 3,5 %, A:B 6,0 % | inchangé | |
| twoway_sk | idem, exponentiel | A 3,5 %, B 7,0 %, A:B 5,5 % | inchangé | |
| ancova | 3 × 15, covariable | 7,5 % | 7,5 % | à 1,7 écart-type de 5 % : à surveiller, non conclu |

Après k, hetero_c se répartit en Welch 159 fois (8,8 % de rejets), `med1way` 36 fois (5,6 %) et ANOVA classique 5 fois (100 %).

Post-hocs (plus petite p-value des comparaisons) : 2 à 3 % dans les scénarios calibrés. Dans 0 à 4,5 % des jeux, une comparaison est significative alors que le test global ne l'est pas. Ce n'est pas un défaut : `m.test` affiche toujours les post-hocs, et ce choix est documenté dans le code.

## 4. Décisions et points ouverts

### 4.1 Stratégie en deux temps sous hétéroscédasticité (non modifiée)

L'inflation qui reste (hetero_c 10,5 %, hetero_c2 22,5 %, hetero2 8,5 %) vient de la logique « test de variance, puis test qui suppose ou non l'égalité ». Avec de petits groupes, le test de variance rate l'écart quand l'échantillon sous-estime par hasard la variance du petit groupe, et c'est là que l'ANOVA classique ou le Student rejettent le plus. Dans hetero_c2, l'ANOVA classique est encore choisie 30 fois sur 200 et rejette 77 % de ces fois. À deux groupes (n = 6/30, sd = 3/1, N = 4000), Student rejette 30 %, Welch 5,6 %, Wilcoxon 16 %, et `var.test` rate l'écart 11 % des fois.

La correction de fond serait d'appliquer Welch par défaut (t de Welch, F de Welch), sans test préalable de variance (Zimmerman, 2004 ; Delacre et al., 2017 et 2019). Elle change ce que les étudiants voient et apprennent (Fisher-Snedecor puis Student) : **c'est à Antoine de décider**, et ce n'est pas fait ici.

### 4.2 Autres constats, non corrigés

- **Lenteur avec `return = FALSE`.** Les post-hocs (bootstrap compris) sont calculés même quand seule la p-value est demandée : 47 s contre 0,85 s avec `boot = FALSE`.
- **`$groups` et `$p.value` vides en modèle mixte** (voie `lmer` existante et split-plot) : `.posthoc_mixed_model()` rend une liste par effet, que `m.test()` ne remonte pas. Les comparaisons s'affichent bien.
- **`$global_pvalue` d'un plan à plusieurs facteurs** = p-value du premier terme (0,199 pour G dans le split-plot, alors que G:T p = 0,030). Ce choix est fixé par un test existant, mais il trompe un lecteur pressé. `$effects` donne tous les termes.
- **Doublon en mesures répétées** : « ASSOMPTION DE BASE : Indépendance… » puis « ASSOMPTION 1/3 : Indépendance… ».
- **Note emmeans** « adjust = "tukey" was changed to "sidak" » imprimée même avec `verbose = FALSE`.
- **Documentation MANOVA** : elle annonce une « Permutation MANOVA » et un « Kruskal-Wallis multivarié » qui n'existent pas dans le code.
- **roxygen** : `mbay-lib` contient roxygen2 8.1.1, qui réécrit tout le NAMESPACE. Le `.Rd` de `.manova_analysis` a donc été corrigé à la main.

## 5. Vérifications

- `tests/testthat/test-m.test-audit.R` : 12 tests, 37 assertions. Chacun a été lancé sur le commit précédant sa correction (échec), puis après (succès).
- Non-régression au dernier commit : `test-m.test.R` (20 tests), `test-manova.R` (22 tests, 17 s au lieu de plus de 6 min) et `test-valreg.R` (6 tests), 0 échec.
- `R CMD check --no-manual --ignore-vignettes` (R 4.6.1, bibliothèque `mbay-lib`) sur dc74ae7 : 0 erreur, 0 avertissement, 2 NOTE (le dossier `.claude` inclus dans le build ; NEWS.md sans titre de version). Les deux sont corrigées au commit suivant : `^\.claude$` ajouté à `.Rbuildignore`, version 0.0.2.3 en titre de NEWS.md et dans DESCRIPTION. Relancé sur 2a1e8b1 : **Status: OK** (0 erreur, 0 avertissement, 0 NOTE), tests `testthat.R` compris.

## Références

- Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of variances. *Journal of the American Statistical Association*, 69(346), 364–367.
- Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by default use Welch's t-test instead of Student's t-test. *International Review of Social Psychology*, 30(1), 92–101.
- Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology*, 32(1), 13.
- Nath, R., & Pavur, R. (1985). A new statistic in the one-way multivariate analysis of variance. *Computational Statistics & Data Analysis*, 2(4), 297–315.
- Zimmerman, D. W. (2004). A note on preliminary tests of equality of variances. *British Journal of Mathematical and Statistical Psychology*, 57(1), 173–181.
