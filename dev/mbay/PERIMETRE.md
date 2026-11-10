# Périmètre de m.bay() établi à partir du code de m.test()

Worktree : `KefiR-mbay`, branche `feat/m.bay`, commit 3675fc7 (KefiR 0.0.2.2). Aucun fichier du dépôt modifié.

**Conventions.** Les références `fichier:ligne` sont relatives à `R/`, et `m.test.R` y est abrégé en `m:`, `sys_multi_factor_analysis.R` en `multi:`, `sys_one_factor_analysis.R` en `one:`, `sys_ancova_analysis.R` en `anc:`, `sys_manova_analysis.R` en `man:`, `sys_mixed_model_analysis.R` en `mix:`.
« [exéc. n] » renvoie à un cas réellement exécuté (scripts et sorties dans le scratchpad : `run_cases.R/.out` = cas 1 à 34 ; `run_cases2.R/.out` et `run_cases2b.R` = cas A à V ; `run_cases3.R` = traces). « [lu] » = établi par lecture seule. « non vérifié » = ni lu en détail ni exécuté.

**Limite de couverture.** `sys_multi_factor_analysis.R` (4615 lignes), `sys_posthoc.R` et `valreg.R` ont été lus par extraits ciblés, pas intégralement. Le cas Q (MANOVA non normale) n'a pas abouti : arrêté après plus de 6 minutes sans résultat, cause non établie.

---

## 0. Ce qu'il faut retenir pour décider

1. **Il n'existe aucune fonction de normalisation des entrées réutilisable** : tout est en ligne dans le corps de `m.test()` (BLOCS 3 à 10, m:262-1031). `m.bay` ne peut pas l'appeler sans extraction, donc sans toucher à `m.test`, ou devra le réécrire.
2. **Plusieurs formes d'entrée documentées échouent** : `formula=` nommé, `g` data.frame, `id` vecteur sans `data`. Ne pas les prendre comme contrat acquis (§2.3).
3. **Pas de transformation automatique** : ni Box-Cox ni log. Seul `log(y)` écrit par l'utilisateur est géré (§5.1). La prémisse « transformation Box-Cox/log » ne correspond pas au code.
4. **`$effects` ne porte que `term` et `p.value`**, avec un type de SC et une population de termes qui changent selon la route, et il est `NULL` à un facteur (§4.3). C'est le point de correspondance naturel avec une réponse bayésienne par terme, mais il ne fournit pas de spécification homogène.
5. **`$method` est peu fiable hors un facteur** : il annonce « Parametric ANOVA/t-test » alors que le test conduit est Friedman, Scheirer-Ray-Hare, aovp ou t2way [exéc. 17, F, G, I]. La couverture des plans ne se lit pas dans `$method`.
6. **Quatre anomalies de `m.test` touchent le plan et pèsent sur l'étalon de comparaison** (§4.4) : interaction silencieusement retirée en split-plot, plantage en `verbose=FALSE`, post-hocs non appariés avec `(1|id)`, et `bf.test` utilisé comme test de variance.

---

## 1. Paramètres de m.test() (m:208-213)

| Paramètre | Défaut | Sens et usage réel | Pour m.bay |
|---|---|---|---|
| `x` | NULL | Réponse : vecteur, matrice/data.frame (MANOVA), noms ou indices si `data` (m:115-120, m:335-408). Une formule passée ici est transférée (m:262-267) | garder |
| `g` | NULL | Facteur, vecteur, data.frame ; nom/indice si `data` (m:421-479). Un data.frame en 2e position devient `data` si `x` est une formule (m:270-276) | garder ; corriger le cas data.frame (§2.3) |
| `data` | NULL | Data.frame (validation m:411-415) | garder |
| `formula` | NULL | Formule. **Passée par nom elle échoue** : m:280-282 la remplace par l'appel non évalué, qui n'hérite pas de `formula`, puis m:597-600 sort [exéc. U] | garder, en la corrigeant côté m.bay |
| `paired` | FALSE | Force l'apparié ; activé aussi par `id`, `wt`, `Error(id)` (m:539, m:550, m:754). Sans `id`, l'appariement repose sur l'ordre des lignes (one:127-128) [exéc. 14] | garder |
| `id` | NULL | Sujet : nom, indice ou vecteur (m:548-589) | garder ; devient l'effet aléatoire |
| `wt` | NULL | Facteurs intra. Ajoutés à la formule (m:962-975), mais **jamais transmis aux analyses** : les appels m:1318-1324 et m:1354-1365 ne passent pas `within`. Sans effet sur `g` [exéc. 16] | adapter : lui donner un vrai rôle |
| `within` | NULL | Alias de `wt` (m:485-492) | garder comme alias |
| `between` | NULL | Facteurs inter ; transmis à multi et MANOVA seulement (m:1319, m:1355) | adapter : structure du modèle |
| `alpha` | 0.05 | Seuil de tous les tests et lettres ; `iter` en dépend (m:241) | **adapter** : pas de seuil de p. À remplacer par des seuils de BF/ROPE à décider |
| `control` | NULL | Niveau témoin (Dunnett), seulement vers `.posthoc` (m:1718, m:1752, m:1771) | adapter : contrastes contre le témoin |
| `verbose` | TRUE | Raisonnement pas à pas. Mis à FALSE si `return=FALSE` sans valeur explicite (m:226-233) | garder ; `.vbse` réutilisable [lu] |
| `plot` | TRUE | Graphes à lettres (m:1494, m:1836) | adapter : lettres sans objet, postérieurs à montrer |
| `return` | TRUE | **TRUE : objet complet. FALSE : un seul nombre**, la p-value globale, sans post-hoc (m:1406-1453) [exéc.]. Contrairement au roxygen m:148-149, qui parle de « raw model output » | adapter : quel scalaire renvoyer ? |
| `boot` | TRUE | Bootstrap des post-hocs ; n'est pas transmis aux analyses (m:1342-1344) | sans objet |
| `boot_type` | NULL | mean / median / meanbp, vers `.posthoc` | sans objet |
| `iter` | 0 | 0 devient `1/alpha*5`, soit 100 (m:241) | adapter : itérations MCMC |
| `conf` | 0.95 | Niveau des IC bootstrap et emmeans (m:1626, m:1670) | adapter : niveau de l'intervalle crédible |
| `maxcat` | 50 | Refus au-delà de ce nombre de groupes (m:1021) ; limite aussi les graphes | garder |
| `silent` | TRUE | Utilisé seulement m:1046 ; `.one_factor_analysis` a son propre défaut (one:62) | sans objet ou à simplifier |
| `code` | FALSE | Recette R ; force `verbose=FALSE` (m:235) | adapter : code de reproduction, à décider |
| `debug` | FALSE | Diagnostics ; force `verbose=FALSE` et `code=FALSE` (m:236-239) | garder |

Paramètres qui n'existent pas dans `m.test` et qu'il faudra ajouter : a priori, graine, nombre de chaînes, moteur. Hors périmètre du code lu.

---

## 2. Formes d'entrée et normalisation

### 2.1 Où se fait la normalisation

| Étape | Lignes | Effet |
|---|---|---|
| Formule en position `x` | m:262-277 | transfert vers `formula` |
| `data$` sans `data` | m:285-312 | retrouve le data.frame dans l'environnement appelant, nettoie la formule |
| `factor(F)` | m:316-319 | `.strip_factor_calls` |
| Variable absente | m:322-328 | erreur qui la nomme |
| `x` | m:335-408 | matrice à plusieurs colonnes : `check_manova` |
| `g` | m:421-479 | facteurs ; plusieurs colonnes : `check_anova` |
| `wt`, `id` | m:485-589 | colonnes temporaires dans `data` |
| `Y ~ F \| bloc` | m:630-647 | devient `Y ~ F + Error(bloc)` |
| `(1\|id)`, `(var\|id)` | m:648-686 | **extrait `id` et efface le terme**, y compris une pente aléatoire ; ne pose pas `paired` |
| `Error(id...)` | m:694-769 | id extrait si niveaux/n ≥ 5 %, sinon traité comme facteur de strate (m:722) |
| Contrôle du groupement | m:780-839 | numérique avec > 50 % de valeurs uniques **ou > 50 valeurs uniques** : refusé comme « continu » (m:793) |
| Extraction `x`, `g` | m:866-945 | prédicteur unique converti en facteur même s'il est numérique (m:937) |
| Construction de l'`Error` | m:950-998 | `Error(id/wt)` ou `Error(id)` |
| NA | m:1040-1068 | suppression listwise à un facteur ; conservés en multi et MANOVA |
| Petits groupes | m:1077-1177 | §5 |

Conséquences de lecture, non exécutées :
- un `id` numérique à plus de 50 sujets est refusé par m:793 s'il n'est pas déjà facteur ;
- un `Error(id)` avec moins de 5 % de niveaux par rapport à n (peu de sujets, beaucoup de mesures) n'est pas traité comme un sujet (m:722).

### 2.2 Formes testées

| Forme | Résultat | Réf. |
|---|---|---|
| `m.test(y ~ F, data=d)` | OK | exéc. 1, 2 |
| `m.test(x, g)` vecteurs | OK | exéc. 7 |
| `m.test("y", "F", data=d)` | OK | exéc. 10 |
| `m.test(d$y ~ d$F)` sans data | OK | exéc. 9 |
| `y ~ factor(F)` | OK | exéc. 21 |
| `log(y) ~ F` | OK, test sur le log | exéc. 22 |
| `cbind(y1,y2) ~ F` ; `A*B` | OK | exéc. 5, 28 |
| matrice `x`, facteur `g` | OK | exéc. 11 |
| `y ~ F + Error(id/F)` | OK, apparié | exéc. 6a |
| `y ~ F \| id` | OK, apparié | exéc. 13 |
| `y ~ F + (1\|id)` | OK, mais post-hocs non appariés (§4.4) | exéc. 6b |
| `y ~ T + (T\|id)` | Passe, **pente aléatoire effacée** (m:670-673) | exéc. 19 |
| `paired=TRUE` sans id, x,g | OK | exéc. 14 |
| `y ~ F, id="id"`, 2 niveaux | OK (multi délègue à one, multi:1438-1447) | exéc. 15, K, L |
| `y ~ F, id="id"`, 3 niveaux | OK | exéc. J |
| `within=` alias | OK | exéc. 17 |

### 2.3 Formes qui échouent

| Forme | Résultat | Cause |
|---|---|---|
| `formula = y ~ F` nommé | **Erreur** « doit être une formule valide » [exéc. U] | m:280-282 |
| `x` vecteur, `g` data.frame sans `data` | **Erreur** « colonnes non définies » [exéc. 8, A] | `.filter_small_groups` cherche les colonnes dans un data.frame qui ne contient que `.x_temp` (m:1102 ; `sys_filter_small_groups.R:84`) |
| `"y", c("A","B"), data=d` | **Erreur** [exéc. A2] | `.formulator_safe(NULL)` en BLOC 12 (m:1192) |
| `x`, `g`, `id` vecteur, 3 niveaux, sans formule | **Erreur** [exéc. T] | `.detect_model_type` reçoit `formula=NULL` (multi:459) |
| `x` sans `g` | Erreur attendue [exéc. 12] | m:1013 |
| `y ~ c1` covariable seule | Erreur trompeuse sur les « effectifs insuffisants » [exéc. 26] | m:937 factorise, chaque niveau a n=1 |
| `y ~ c1 + c2` | Erreur « au moins 2 groupes » [exéc. 27] | `.n_levels_g` = 1 (m:28, m:1027) |
| `y ~ G, id, wt="T"` | Pas d'erreur, mais **`T` ignoré dans le résultat** [exéc. 16] | §3 |

**Note.** D'après m:1192 et m:423-424, `g` data.frame avec `data` et sans formule suit le même chemin que A2 [lu, non exécuté].

### 2.4 Fonctions internes : contrat et réutilisabilité

| Fonction | Entrée | Sortie | Réutilisable telle quelle ? |
|---|---|---|---|
| `.normalize_from_formula` : **deux définitions**, m:61-95 et `sys_normalize_from_formula.R:1-60`. Celle de `sys_` l'emporte (formals vérifiés [exéc.]) | formule, data, id, wt_names, `na_strategy` ("preserve"/"listwise") | classe `norm_form` : `x`, `g` (toujours data.frame), `mf`, `keep`, `had_na`, `response`, `predictors` | **Appelée nulle part** (grep). `(1\|id)` n'est pas géré : colonne `"1 \| id"` dans `mf` [exéc.]. `g_out` garde id et wt (sys_normalize_from_formula.R:46), contrairement à `g_df` calculé puis inutilisé. Non, pas sans correction |
| `.drop_error_term(f)` (`sys_drop_error_term.R:5`) | formule | formule sans `Error()` ni `(a\|b)` | **Oui** [exéc.]. Développe `G*T` en `G + T + G:T` ; environnement non conservé en voie Error (`reformulate`, l.47) : passer `data` |
| `.strip_data_dollar_safe(f, data)` | formule, data | retire `data$`, à condition que l'objet s'appelle exactement comme l'argument | Oui avec réserve : `deparse(substitute(data))` (l.3) rend le résultat dépendant du nom de la variable appelante |
| `.normalize_formula_dollar(f, data)` (`sys_normalize_formula_dollar.R:1`) | idem | idem | Même réserve (l.4) |
| `.strip_factor_calls(f, data)` (`sys_normalize_formula_dollar.R:32`) | formule, data | `list(formula, data)` | **Oui** [lu, exéc. 21 via m.test] |
| `.formulator_safe(formula, data, mode)` (`sys_formulator_safe.R:1`) | formule, mode simple/linear/lmer | formule | Déconseillée : `linear` aplatit tout, `y~G*T+Error(id/T)` donne `y ~ G + T + id` [exéc.] ; passe en `lmer` dès qu'une barre est présente |
| `.formulator` (`sys_formulator.R:87`) | formule, mode aov/friedman/aovp/linear/simple/lmer | formule | Mode `aov` : convertit `(1\|r)` en `Error(r)`. Sens inverse de ce qu'il faut à m.bay |
| `.detect_model_type(formula, data, debug)` (`sys_detect_model_type.R:69`) | formule, data | **chaîne "TRUE"/"FALSE"**, malgré la doc qui annonce un booléen (l.27 ; class = character [exéc.]) ; `.exit` si aucun facteur | Oui, en comparant à "TRUE". Ne convertit les 0/1 qu'en local (l.117) |
| `.analyze_ancova_structure(formula, data, alpha, ...)` (`sys_analyze_ancova_structure.R:58`) | formule, data | `model_type`, `has_interaction`, `interaction_terms`, `factor_vars`, `covariate_vars`, `hierarchical_test`, `decision`, `k` | La partie descriptive (l.73-130) serait utile ; la décision est fréquentiste (l.225). À ne pas appeler telle quelle |
| `.n_levels_g(g)` (m:18) | facteur ou data.frame | produit des niveaux, covariables exclues | **Oui** |
| `.filter_small_groups(data, factor_vars, min_n, ...)` (`sys_filter_small_groups.R`) | data, facteurs | `data`, `removed_groups`, `removed_n`, `original_n`, `filtered`, `k` | Oui, si les colonnes sont dans `data` |
| `.effects_from_table`, `.first_factor_pvalue`, `.multi_effects_table` (`sys_effects_table.R:6, 30, 48`) | table ANOVA, modèle | §4.3 | Sans objet pour un BF, mais la structure `term` + valeur est à copier |
| `.letters_from_pairs(p_lower, position, alpha)` (`sys_effects_table.R:119`) | matrice de p | lettres via `catego()` | Sans objet tel quel ; adaptable si l'on veut des lettres |
| `.vbse`, `.msg`, `.exit`, `.dbg`, `.format_pval` | texte | sorties | **Oui** |
| `discret.test(vector)` (`discret.test.R`) | vecteur | TRUE si valeurs uniques < floor(√n) | **Oui** [exéc.] |

---

## 3. Plans traités et routes

**Routage de tête (m:1308-1376).** Drapeaux `paired`, `check_manova`, `check_anova`. `paired` + `Error` ou (`id` et > 2 niveaux) force `check_anova` (m:1230-1232). Un `y ~ F, id=` passe donc par multi, qui délègue à one pour 2 niveaux sans réplicats (multi:1438-1447).

**Post-hocs (m:1458-1827).** MANOVA → `.posthoc_MANOVA` ; `robust_results$method == "Mixed_Model_lmer"` → `.posthoc_mixed_model` ; ANCOVA → `.posthoc_ANCOVA` ; sinon `.posthoc`.

| Plan | Route | Tests | Résultat fréquentiste | Réf. |
|---|---|---|---|---|
| 2 groupes, normal, variances égales | one:1380-1442 | `var.test`, `t.test` | Student ; `$method "Student t-test"` | one:1414, 1419 ; exéc. 1 |
| 2 groupes, normal, variances inégales | one:1443-1469 | `t.test(var.equal=FALSE)` | Welch | one:1449 ; m:2024 [lu] |
| 2 groupes, non normal | one:817-1000 | `ks.test` (fiabilité), `wilcox.test` | Wilcoxon-Mann-Whitney | one:922 ; exéc. C |
| 2 groupes, discret et non normal | one:819-881 | `mood.test`, `ansari.test`, `wilcox.test` | Wilcoxon | one:826-862 [lu]. Discret mais jugé normal : Student [exéc. D] |
| k ≥ 3, normal, variances égales | one:1475-1494 | `aov` ; Student, Scheffé, Tukey, SNK, Waller-Duncan | ANOVA ; « Parametric ANOVA/t-test » | one:1480 ; exéc. 2 |
| k ≥ 3, normal, variances inégales | one:1495-1531 | `oneway.test`, `fanova.hetero` ; Games-Howell, Dunnett T3 | « Welch ANOVA » | one:1499 ; exéc. 22 |
| k ≥ 3, non normal | Arbre one:1147-1196 | outliers extrêmes > 5 % ou \|skew\| ≥ 2 : `med1way` ; formes semblables : Kruskal ; formes différentes + Fligner significatif : `t1way` ; sinon `med1way`. Post-hocs : Wilcoxon, Dunn, Brunner-Munzel, medpb2 | `$chosen_test` = med1way / t1way / kruskal | exéc. 30 (Kruskal), E (med1way) ; t1way [lu] |
| Facteur unique multi, non normal | multi:3508 | `kruskal.test` | `term = "global"` | [lu] |
| A\*B, équilibré, conditions OK | multi:1476-1531, 4196-4262 | `aov` ; SC Type I | `$effects` A, B, A:B | exéc. 3 |
| A\*B, déséquilibre ≤ 2 | multi:1575-1652 | `lm` contr.sum, GTest (`control_independence`) ; Type II ou III | `$effects` | exéc. 34 |
| A\*B, déséquilibre > 2 | multi:1557-1572, 3605 | Scheirer-Ray-Hare | `$effects` | exéc. H |
| A\*B, déséquilibre > 10 | multi:1540-1555, 3322-3370 | `WRS2::t2way` / `t3way`, tr adaptatif 5/10/20 % (20 % non lu) | `$effects` A, B, **AB** (nom différent) | exéc. I |
| A\*B, non normal ou hétéroscédastique | multi:2254, 2307, 3605 | Scheirer-Ray-Hare | `$effects` | exéc. F |
| ≥ 3 facteurs, non normal | multi:3704-3734 | `lmPerm::aovp` | `$effects` | exéc. G |
| ANCOVA, pentes parallèles, 1 ou plusieurs covariables | multi:459-504 → anc:49 | `lm` contr.sum, `car::Anova` Type III ; emmeans/Tukey | moyennes ajustées | exéc. 4, 25 |
| `F*cov`, interaction p ≥ alpha | `sys_analyze_ancova_structure.R:250-253`, anc:244-246 | `anova(réduit, complet)` puis **retrait** | formule réduite | exéc. 24 ; test-m.test.R:29-38 |
| `F*cov`, interaction p < alpha | anc:241-243, 259-262 | Pentes hétérogènes, covariable centrée, Type III | `$effects` F, cov, F:cov | exéc. 23 ; test-m.test.R:40-48 |
| ANCOVA 1 facteur, violations | anc:1661-1737 | `MASS::rlm` MM, Wald robustes | sans F global | test-m.test.R:138-146 [lu] |
| ANCOVA ≥ 2 facteurs, violations | anc:1739-1921 | HC3 (`sandwich`), repli `rlm`, puis classique (anc:1939-1971) | `$effects` | test-m.test.R:148-155 [lu]. Cas N : la route exacte n'est pas établie |
| Covariable seule | **non pris en charge** | `.exit` ou erreur | aucun | exéc. 26, 27 |
| ANCOVA appariée | **non pris en charge** | `.exit` | aucun | anc:82 ; exéc. P |
| MANOVA 1 facteur | man:900-1090 | Assomptions (Mahalanobis, Mardia, HZ, Box M), `manova` ; ANOVA protégées Bonferroni, discriminante | Wilks ; Pillai si violation (man:1201) | exéc. 5 |
| MANOVA ≥ 2 facteurs | man:918-938, 1063-1089 | `car::Manova` Type II | `$effects` par terme | exéc. 28 |
| MANOVA, ≥ 2 violations | man:876-877, 1104-1122 | `rrcov::Wilks.test` MCD, graine 2026 ; repli Pillai (man:1151-1166) | test sur groupes **croisés** | test-manova [lu] ; **cas Q non abouti** |
| MANOVA répétée, MANCOVA | **non pris en charge** | `.exit` | aucun | man:459, 510 ; exéc. 29, R |
| Apparié 2 niveaux | one:94-180 | normalité des différences, `t.test` ou `wilcox.test` appariés | `$method` exact | exéc. K, L, 14, 15 |
| RM 1 facteur, ≥ 3 niveaux | multi:606-, 1505-1522, 2018 | `aov+Error`, Mauchly (`rstatix::anova_test`), Student apparié | `$effects` par strate | exéc. 6a, J |
| RM non normal, 1 facteur | multi:2984-3023 ou 3400-3439 | `friedman.test`, Nemenyi | `term="global"` | exéc. 17 |
| Split-plot `G*T + Error(id/T)` | Auto-détection de tous les facteurs en intra (multi:666) → « problèmes » → mixte (multi:794-796) | `lmer` | **`y ~ G + T + (1\|id)`, interaction perdue** | exéc. B, B3, B4 ; mix:148-164 |
| Split-plot `G*T + (1\|id)` | `paired` reste FALSE | `aov+Error`, robuste, puis aucune méthode applicable | `$effects` G, T, G:T ; post-hocs non appariés | exéc. B2 |
| Split-plot `y ~ G, id, wt=T` | `data` reconstruit sans T (multi:397) | sans modèle pour T | Effet de T absent | exéc. 16 |
| Mixte avec réplicats | multi:731-800 → mix:77 | `lmerTest::lmer`, `valreg`, ANOVA Type III Satterthwaite ; emmeans (si effet significatif) | `$method "EMMs"` en mixte réel | exéc. B ; exéc. 20 ne passe pas par là |
| `y ~ F + (1\|id)` avec réplicats | `aov+Error` sur id, puis `.posthoc` non apparié | Wilcoxon non apparié | global 2.27e-08, `$p.value` 0.059 | exéc. 20 |
| Pentes aléatoires | **Effacées en entrée** (m:670-673) ; `random_slope` de `mix:85` : multi ne le passe qu'à FALSE (multi:3170, grep) | intercept seul | pente non estimée | exéc. 19 |
| Petits effectifs | Filtrage n < 3 (m:1077) ; tolérance Blanca exige n ≥ 15 (`sys_normality.R:199`) | selon n | exéc. 33, S |
| Hétéroscédasticité, 1 facteur | `.variance` : Bartlett si normal, sinon Levene puis `bf.test` | **Levene contredit par un test de moyennes** (§4.4) | exéc. 31 |

**Normalité (`sys_normality.R`).** Shapiro si n ≤ 50, Jarque-Bera si n ≤ 500, seuils d'asymétrie/kurtosis au-delà (l.164-170) ; Šidák par groupe ; retest `tolerance="extrem"` (Chaffin & Rhiel, Blanca, l.184-203).

---

## 4. Structure de retour

### 4.1 Forme selon la route

Toutes les sorties avec `return=TRUE` ont la classe `posthoc`, construite en m:2057-2073 (univarié) ou m:1568-1583 (MANOVA). Avec `return=FALSE` : un numeric seul (m:1452) [exéc. : 1.028506e-09, 1.835525e-05].

| Élément | y~F 2 niv. (1) | y~F 3 niv. (2) | y~A\*B (3) | y~F+cov (4) | cbind~F (5) | Error(id/F) (6a) | (1\|id) (6b) |
|---|---|---|---|---|---|---|---|
| `$groups` | df 2×3 | df 3×7 | df 4×7 | df 3×2 | df 2×3 **(VD)** | df 3×3 | df 3×7 |
| `$p.value` | num | matrice 2×2 | matrice 3×3 | matrice 3×3 **pleine** | named num (par VD) | matrice 2×2 | matrice 2×2 |
| `$global_pvalue` | 6.5e-04 | 1.0e-09 | 1.8e-05 (A) | 4.1e-09 (F) | 2.1e-12 | 3.3e-09 | 3.3e-09 |
| `$descriptive_stats` | group, mean, SD, SE, n, type | idem | idem (cellules) | group, adjusted_mean, SE, type | **absent** | idem | idem |
| `$method` | "Student t-test" | "Parametric ANOVA/t-test" | idem | "ANCOVA post-hoc (tukey adjustment)" | "MANOVA (Wilks' Lambda)" | "Parametric..." | "Parametric..." |
| `$effects` | absent | absent | term, p.value | term, p.value | `NULL` | term, p.value | term, p.value |
| `$bootstrap` | list(groups, p.value) | idem | idem | absent | absent | idem | idem |
| `$details` | | | | 5 éléments (pairwise_comparisons, emm_result, pairs_result, note, adjusted_means) | | | |
| Spécifique | | | | | test_statistics (Wilks, Pillai, HotellingLawley, Roy : statistic, F_value, df1, df2, p_value), discriminant_analysis (n_functions, prop_variance, structure_matrix), protected_anovas, effects_by_dv, pairwise_multivariate, groups_multivariate | | |

### 4.2 Détails de forme

- `$groups` : `categories` + une colonne de lettres par post-hoc. Le nom et le nombre de colonnes varient. Une **même lettre n'a pas le même sens entre colonnes** : Bootstrap et Student vont de a à c, Scheffé et Tukey en sens inverse [exéc. 2].
- `$groups` en multi-facteurs porte sur les **cellules croisées** (`a1.b1`...), pas sur les facteurs.
- `$groups_by_factor` n'apparaît qu'en ANCOVA à plusieurs facteurs (m:1685-1702 ; test-m.test.R:157-165).
- En MANOVA, `$groups` décrit les VD, pas les groupes. `y2` y est « ns » à p = 0.026 à cause de la correction de Bonferroni (`sys_posthoc_MANOVA.R:317`) [exéc. 5]. Les groupes sont dans `$groups_multivariate`.
- `$chosen_test` n'existe que si non NULL (m:2070), donc seulement pour la branche non normale à un facteur [exéc. 30, E].
- `$global_pvalue` en un facteur est le test global ; en multi, celui du **premier facteur de la formule** (`sys_effects_table.R:21-37`, BP-031, multi:4537-4541). Avec une interaction présente, ce scalaire porte sur le seul effet principal du premier facteur.

### 4.3 `$effects` et sys_effects_table.R

**Colonnes.** `term` (chr), `p.value` (num). Ni statistique, ni ddl, ni taille d'effet. Exception MANOVA multi-facteurs : `term, Wilks, p.Wilks, Pillai, p.Pillai` [exéc. 28].

**Absent ou NULL** : un facteur (one: ne renvoie pas `effects`, one:1551), MANOVA à un facteur [exéc. 5], `y ~ G, id, wt` [exéc. 16]. `.effects_from_table` retire `Residuals` et `(Intercept)` (l.14).

| Route | Fonction source | Type de SC | Termes | Réf. |
|---|---|---|---|---|
| ANOVA 2 facteurs équilibrée | `anova(model)` | **I** | A, B, A:B | `sys_select_ss_type.R:131` ; test-m.test.R:66-75 ; exéc. |
| 2 facteurs déséquilibré, interaction ns ou absente | `car::Anova` | **II** | A, B, A:B | `sys_select_ss_type.R:156, 176` ; exéc. (II) |
| 2 facteurs déséquilibré, interaction significative | `car::Anova` | **III** | | `sys_select_ss_type.R:165` [lu] |
| Autres ANOVA (dont ≥ 3 facteurs) | `car::Anova` | III par défaut | | `sys_select_ss_type.R` résultat initial ; exéc. (A+B+C = III) |
| ANCOVA | `car::Anova`, contr.sum | III en pratique ; I possible (`sys_select_ss_type.R:330`, non exécuté) | F, cov, cov2, F:cov | anc:1513 ; exéc. 4, 23, 25 |
| Mesures répétées | `summary(aovlist)` par strate | séquentiel (I) | par strate, empilés | `sys_effects_table.R:85-94` |
| Mixte lmer | `anova(type="III")` Satterthwaite | III | | `sys_effects_table.R:54-57` ; exéc. B |
| WRS2 t2way/t3way | `.p.value` | sans objet | **A, B, AB** | l.59-65 ; exéc. I |
| Kruskal, Friedman | `htest` | sans objet | **"global"** | l.67-69 ; exéc. 17 |
| aovp | `summary` | non vérifié | A, B, A:B... | l.71-77 ; exéc. G |
| Scheirer-Ray-Hare | table | SC sur rangs | A, B, A:B | l.79 ; exéc. F |
| rlm, HC3 | Wald ou `car::Anova(vcov.)` | III, repli II | | anc:1708-1721, 1779-1785 |
| MANOVA multi | `car::Manova` | **II** | | man:929 |

**Correspondance avec une réponse bayésienne « pour chaque terme ».**
- Un terme = une ligne : la structure `term` + valeur se transpose directement, en y ajoutant BF, sens du BF, probabilité postérieure de l'inclusion.
- Pas d'équivalent direct du **type de SC** en bayésien. La comparaison de modèles fixe ce que le test d'un terme signifie : retirer un terme en gardant ses termes d'ordre supérieur (analogue type III), ou l'exclure avec eux (analogue type II, marginalité). **Décision à prendre**, à documenter ; une valeur unique rendrait `m.bay` plus cohérent que `m.test`, qui change de type selon l'équilibre et une p-value.
- Pour un facteur, `m.bay` devrait produire `$effects` à une ligne, sans l'absence de `m.test`.
- Pour RM et split-plot, `m.test` empile des strates ; le bayésien aura un seul modèle mixte, et `Error` devra être traduit.
- Le terme de premier facteur ne doit pas être le seul scalaire résumé sans réflexion : `global_pvalue` porte sur un effet principal, alors que la question bayésienne naturelle est le modèle complet contre le nul, ou un terme désigné.

### 4.4 Échecs et surprises à l'exécution (rien n'a été corrigé)

**Touchent la validité des résultats**

1. **Split-plot : interaction retirée sans message.** `y ~ G*T + Error(id/T)` ajuste `y ~ G + T + (1|id)` : `$effects` ne contient pas G:T [exéc. B, B3]. Mécanisme : tous les facteurs sont déclarés intra (multi:666), le plan paraît déséquilibré (multi:794-796), puis `.mixed_model_analysis` construit la partie fixe avec `paste(collapse=" + ")` (mix:148-164). `between` n'est pas déduit.
2. **Plantage en `verbose=FALSE` sur le même plan.** « $ operator is invalid for atomic vectors » [exéc. 18, verbose=TRUE passe]. Pile jusqu'à `.mixed_model_analysis`. Cause lue : `valreg` ne renvoie une liste que si `k != 0` (valreg.R:84-85, 653-658), et `k` reste à 0 sans verbose ; mix:303 fait alors `valreg_output$valid` sur un logique. **Le résultat dépend de `verbose`** : à éviter dans m.bay.
3. **`(1|id)` : global sur le modèle à sujet, post-hocs non appariés.** `paired` reste FALSE (m:648-686 ne le pose pas). Mêmes données, t1-t2 : 0.139 avec `(1|id)`, 0.00208 avec `Error(id/F)` [exéc. 6a, 6b]. Tukey et Scheffé n'existent que si `paired == FALSE` (`sys_posthoc.R:626`). Cas extrême : global 2.27e-08, `$p.value` 0.059 [exéc. 20].
4. **`bf.test` utilisé comme test de variance.** `sys_variance.R:240` l'emploie pour confirmer Levene et commente « basé sur la médiane ». `onewaytests::bf.test` est une ANOVA de **moyennes** robuste à l'hétéroscédasticité (documentation lue en R). Vérification : données SD 0.5/2/5, Levene p = 2.7e-08, `bf.test` p = 0.80 [exéc.]. L'hétéroscédasticité forte est donc écartée quand les moyennes sont égales, ce qui donne une ANOVA classique [exéc. 31]. Hors périmètre m.bay, mais à signaler à Antoine.

**Autres**

5. `$method` trompeur : « Parametric ANOVA/t-test » avec colonnes Wilcoxon/Nemenyi et Friedman [exéc. 17], Scheirer [F], aovp [G], t2way [I], RM non normal.
6. Post-hocs incohérents avec le global : déséquilibre 3:1, global Scheirer-Ray-Hare, post-hocs Student/Tukey [exéc. H].
7. `.posthoc_mixed_model` ne fait rien si aucun effet n'est significatif (`sys_posthoc_mixed_model.R:95-102`), et renvoie `list(posthoc_results=NULL, k)` sans `groups` [lu].
8. **Cas Q non abouti** (MANOVA non normale avec `verbose=TRUE`, `boot=FALSE`) après plus de 6 minutes ; arrêté. Cause non établie (`rrcov`, `MVN` ou `.posthoc_MANOVA`, non vérifié). Le test `test-manova` équivalent (test-m.test.R:178-183) porte sur des données différentes.
9. Temps : 16 s pour y~F à 3 niveaux, 34,5 s au premier ANCOVA (cas 4), à cause du bootstrap et des chargements de packages (non décomposé) ; 0,8 s pour un ANCOVA ultérieur.

---

## 5. Comportements automatiques à trancher

**Marqueur.** **[p]** = la décision dépend d'une p-value ou d'un test. Le bayésien ne devrait pas reproduire ces points sans réflexion.

### 5.1 Comportements demandés explicitement

| Comportement | Ce que fait m.test | Réf. | Question pour m.bay |
|---|---|---|---|
| Box-Cox, log automatique | **N'existe pas.** Mention seulement dans un message à l'utilisateur (mix:314, 320). Pas de Box-Cox dans `valreg.R`, `bootreg.R`, `evolreg.R` (grep) | | Rien à reproduire. Modéliser (Student, log-normal) plutôt que transformer ? |
| Réponse `log(y)` de l'utilisateur (**BP-033**) | `x` déjà transformé ; la formule référence `y` seul pour éviter `log(log(y))` (multi:379-392 ; `dev/bp.log:882`) | test-m.test.R:77-87 | **Reprendre cette règle** : piège identique si m.bay reconstruit `data` |
| Retrait de l'interaction facteur×covariable **[p]** | `anova(réduit, complet)`, retrait si p ≥ alpha (`sys_analyze_ancova_structure.R:215-272` ; anc:244-246). Conséquence : la formule de l'utilisateur est modifiée | exéc. 24 | Comparer par BF ou LOO, et ne jamais retirer silencieusement ? Moyenner sur les modèles ? |
| Homogénéité des pentes testée même sans interaction écrite **[p]** | `lm` avec interaction ; p < alpha : `check_slopes_homogeneity <- FALSE` et voie robuste (anc:1337-1355, 1623) | | Terme testé par BF, sans changement de modèle |
| Filtrage des petits groupes | n < 3 : **supprimé avec message** (m:1077, `sys_filter_small_groups.R:91, 129`). À 2 groupes restants, le plan change [exéc. 33]. Porte sur les croisements en multi-facteurs | | **Point fort à discuter** : le bayésien sait traiter de petits groupes par rétrécissement, a priori et postérieur large. Garder un refus seulement si n = 0 ou 1 ? |
| Outliers | `rstatix::identify_outliers` ; > 5 % extrêmes : `med1way` (one:1037, 1148). Multi : règle le trim WRS2 (multi:3328-3339) ; **aucune suppression de données** [lu sur les extraits] | exéc. E | Modèle à queues lourdes (Student) par défaut ? Pas de suppression, diagnostic seulement |

### 5.2 Décisions à partir d'un test de p-value [p]

| Décision | Réf. | Effet |
|---|---|---|
| Normalité (Shapiro, JB, Šidák) | `sys_normality.R:164-170, 217-231` ; one:435, 461 ; multi:1941 | paramétrique ou non |
| Retest tolérant (seuils, sans p, conditionné au test précédent) | one:740-786 ; anc:1135-1163 ; multi:2189 | retour au paramétrique |
| Variance : `var.test`, Bartlett, Levene, BF | `sys_variance.R:63, 144, 209, 240` ; one:1414 | Student/Welch ; AOV/Welch ; robuste en multi (multi:2301-2307) |
| Diagnostic KS sur la forme | one:807-812, 1022, 1166 ; `sys_posthoc.R:871` | Kruskal ou t1way/med1way ; ajout de Brunner-Munzel |
| Fligner | one:1014, 1178 (seuil **0.05 en dur**, pas `alpha`) | t1way ou med1way |
| Mauchly | multi:2022 ; correction recommandée (multi:2135). Application effective de GG/HF : non vérifiée | sphéricité |
| Type de SC selon p d'interaction | `sys_select_ss_type.R:154` | Type II ou III |
| `control_independence` (GTest entre facteurs, Holm) : ajoute des interactions | multi:1635-1652 ; `sys_control_independence.R:52` | Modifie la formule ; choix final du modèle non lu en détail |
| Linéarité de la covariable (terme quadratique) | anc:1225-1238 ; man:630-633 | voie robuste |
| Indépendance covariable-facteur | anc:756, 770 | **avertissement seul**, pas de changement [exéc. O] |
| Normalité de la MANOVA, Box M (≥ 2 violations) | man:758-771, 829, 876 | robuste ; Wilks ou Pillai pour la p (man:1201) |
| Significativité pour les post-hocs mixtes | m:1619 ; mix:351, 433 ; `sys_posthoc_mixed_model.R:95` | **pas de post-hoc** si rien de significatif |
| ANOVA protégées | `sys_posthoc_MANOVA.R:317` | verdict significant/ns ; le test est choisi par VD (`dv_route`) |
| Lettres | `catego`, `.letters_from_pairs`, `cld` | groupes selon alpha |
| Repli lettres si `.posthoc` échoue | m:1789 | dérivées de la p globale |

**Point de méthode.** Sur ces décisions, le schéma « tester l'hypothèse, puis changer de méthode » est propre au fréquentiste. Reproduire la même chaîne avec des BF ajouterait les mêmes problèmes de sélection et invaliderait l'interprétation du postérieur. Alternative : un seul modèle, éventuellement robuste, avec vérification prédictive et comparaison de vraisemblances, le choix étant annoncé.

### 5.3 Décisions sans p-value à arbitrer

| Décision | Réf. | Question |
|---|---|---|
| Discret : uniques < floor(√n) | `discret.test.R` ; multi:420 (bascule en robuste) ; one:819 (2 groupes seulement) | Incohérence entre routes. Vraisemblance discrète ? |
| Moins de 5 valeurs distinctes en RM : robuste | multi:1344-1358 | idem |
| Ratio de déséquilibre > 2, > 10 | multi:1540, 1557 | Le bayésien gère le déséquilibre sans changer de méthode : à ne pas reproduire ? |
| Réplicats par sujet : mixte | m:1257-1283 ; multi:731-800 | Le mixte est la voie unique |
| `Error(id)` selon niveaux/n < 5 % | m:722 | Heuristique fragile ; préférer la déclaration explicite ? |
| Prédicteur numérique unique converti en facteur ; 0/1 converti | m:937 ; multi:370 ; `sys_detect_model_type.R:110` | `y ~ dose` est un facteur, `y ~ F + dose` une covariable. Reprendre ? |
| Centrage des covariables en interaction | anc:259-262 (BP-030) | Reprendre, avec le même rapport sur les valeurs d'origine |
| NA : listwise (un facteur), préservés (multi) | m:1040-1068 | Politique unique |
| `id` vu comme continu si > 50 valeurs | m:793 | Ne pas reproduire |
| Post-hocs toujours faits, même si global ns | m:1608, 1767 (BP-013), sauf mixte | Reprendre : cohérent avec le bayésien, sans correction de multiplicité à la fréquentiste |

---

## 6. Proposition de périmètre (pour décision, non validée)

Déduite des sections précédentes ; ce n'est pas une mesure.

- **Couvert par `m.test` et à reprendre** : 2 et k groupes, multi-facteurs avec interactions, ANCOVA (dont pentes hétérogènes et plusieurs covariables), apparié, RM, split-plot, mixte à intercept, MANOVA 1 et plusieurs facteurs.
- **Non couvert par `m.test`** (décision à prendre) : covariable seule, ANCOVA appariée, MANOVA répétée, MANCOVA, pentes aléatoires réelles. Le bayésien peut en couvrir certaines sans effort supplémentaire, notamment les pentes aléatoires.
- **Sans équivalent à reproduire** : tests robustes de la branche non normale (Kruskal, Brunner-Munzel, medpb2, t1way, bootstrap, permutation) ; le remplacement est une vraisemblance robuste ou de rangs, à choisir.
- **Étalon de comparaison** : `m.test` ne doit pas servir de référence pour le split-plot, ni pour `(1|id)` post-hocs, tant que les points §4.4 1 à 3 ne sont pas résolus.
- **À signaler à Antoine indépendamment de m.bay** : §4.4 points 1 à 4, et `formula=` nommé (§2.3).

## 7. Couverture de m.bay au 10 novembre 2026 (branche feat/m.bay)

Légende :
- **Couvert** : calculé, validé, raconté.
- **Refusé** : arrêt avec un message qui nomme le plan et renvoie à `m.test()`.
- **Non couvert** : ni calculé, ni refusé par un message propre à ce plan.

| Plan de `m.test()` | m.bay | Moteur, décision | Limites ou reste à faire |
|---|---|---|---|
| 2 et k groupes, réponse continue | Couvert | 4 lois comparées (normale ou t, variances égales ou non), D6 | — |
| Post-hoc à un facteur | Couvert | Regroupements, D11 ; Gibbs au-delà de 8 groupes | Calibration au-delà de 8 non simulée |
| Plusieurs facteurs, interactions | Couvert | `lmBF`, HOP, D2 | **Loi normale à variance commune seule** ; pas de post-hoc |
| ANCOVA, pentes parallèles ou différentes, plusieurs covariables | Couvert | `lmBF`, covariable centrée, D16 | Normale seule ; post-hoc à un facteur |
| Apparié, 2 niveaux | Couvert | `ttestBF(paired)`, D15 | Normale seule |
| Mesures répétées, 1 facteur intra | Couvert | `lmBF(whichRandom)`, D15 | Symétrie composée supposée, non testée |
| Plan en tranches (inter × intra) | Couvert | `lmBF(whichRandom)`, D15 | Pas de post-hoc |
| Modèle mixte avec répliques par sujet et condition | Refusé | — | À faire : `lmBF` avec sujet et sujet × condition |
| Pentes aléatoires | Refusé | — | Stan (rstanarm) seulement |
| ANCOVA appariée, MANCOVA, MANOVA répétée | Refusé | — | Non couverts par `m.test()` non plus |
| MANOVA, 1 ou plusieurs facteurs | Refusé | — | Formule Normale-Inverse-Wishart prévue (D1) |
| Réponse binaire ou catégorielle, 1 facteur | Couvert | `contingencyTableBF`, D17 voie 1 | Plusieurs facteurs : refusé |
| Comptages, 1 facteur | Couvert | Poisson et binomiale négative, BIC, D17 voie 2 | Plusieurs facteurs : refusé ; puissance modeste |
| Réponse ordonnée (notation bornée), 1 facteur | Couvert | `polr`, BIC, D17 voie 3 | Cotes proportionnelles non testées ; plusieurs facteurs : refusé |
| Tests de rangs (Wilcoxon, Kruskal) | Non couvert, raisons D17 | — | Leurs cas sont repris par les voies 3 et par les lois t |
| Covariable seule (régression) | Refusé | — | Hors `m.test()` |
| Lois t ou variances inégales à plusieurs facteurs | Non couvert (normale seule, annoncée) | — | **Priorité suivante** |
| Voie robuste ou Stan | Non couverte | — | rstanarm en Suggests (D1), non branché |

**Temps typique sur un jeu de TP** (N ≈ 40, valeurs par défaut, `scripts/resultats/agent_63.txt`) :

| Plan | Temps |
|---|---:|
| 2 groupes | 2,8 s |
| 4 groupes | 9,1 s |
| 2 × 2 facteurs | 1,1 s |
| Apparié | 0,5 s |
| Répété, 3 temps | 1,3 s |
| ANCOVA | 1,6 s |
| Oui/non | 0,1 s |
| Comptages | 0,2 s |
| Likert 1-5 | 0,3 s |

Tous restent sous 30 s. Le récit s'imprime étape par étape, donc la première ligne paraît dès l'étape 1. Toute voie future au-delà de 30 s devra afficher une ligne d'attente, comme le fait déjà le cas des plus de 50 modèles.

**Stan sur ce poste.** La compilation échoue : R 4.6 attend Rtools 4.5, et seul Rtools 4.4 est installé (`agent_62`). Un poste étudiant n'a en général pas Rtools. Les voies Stan (brms, rstanarm compilé) sont donc écartées pour l'usage en TP.

**Priorités, pour un TP de BUT** (complétées le 10 novembre avec la session parente).
0. **Post-hoc à plusieurs facteurs, en ANCOVA et en plan en tranches**, à priorité égale avec le point 1 : effets principaux quand l'interaction n'est pas probable, cases quand A:B l'est.
1. **Lois t ou variances inégales à plusieurs facteurs.** Les plans à deux facteurs dont les variances diffèrent selon les cases sont fréquents en TP. Aujourd'hui, ils sont répondus en loi normale à variance commune, avec un simple avertissement. La mise en défaut de D6 (variances inégales, 20 % de faux « oui » à P > 0,95 en voie normale seule à un facteur) laisse craindre le même défaut à plusieurs facteurs.
2. **Modèle mixte avec répliques**, plusieurs mesures par sujet et condition. C'est courant en TP de mesures répétées, et la réponse est aujourd'hui un refus.
3. **MANOVA.** Plus rare en BUT.
