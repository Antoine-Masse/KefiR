# Audit de `m.test()` — novembre 2026

Branche `claude/fix-m-test-audit-c7fc21`, partie de master (3675fc7, version 0.0.2.2).
Rien n'est poussé. Ni master, ni `m.bay`, ni les branches `feat/m.bay`, cor.bay et g.bay ne sont touchés.

**Méthode.** Chaque défaut est reproduit sur des données construites pour l'occasion, puis comparé à une fonction R de référence (`aov`, `car::leveneTest`, `pairwise.t.test`, `rrcov::Wilks.test`…). Chaque correction a son commit et son test dans `tests/testthat/test-m.test-audit.R`. Ce test échoue sur le code d'avant et passe après. Les taux de faux positifs viennent de simulations sous H0 (script `sim_h0.R`, section 3).

**Mesure du temps.** Un appel `m.test()` sur 50 observations coûte environ 1 s avec `boot = FALSE` et 47 s avec la valeur par défaut `boot = TRUE`, même quand `return = FALSE` ne demande qu'une p-value (machine partagée, 8 cœurs). Les simulations utilisent donc `boot = FALSE`. La section 3 vérifie que `boot` ne change pas la p-value globale.

## Synthèse

Dix-huit points, tous reproduits sur des données construites pour l'occasion, puis corrigés, chacun avec son commit et son test (a à g : audits de m.bay ; h à k : cet audit ; l à q : relectures de l'expert indépendant ; r : décision d'Antoine). Rien n'est fusionné : la fusion vers master revient à Antoine.

| | Défaut | Avant | Après | Ce que voit l'étudiant |
|---|---|---|---|---|
| a | Split-plot : interaction G:T retirée sans message | `$effects` = G, T ; G:T (p = 0,030) absente | G, T, G:T ; G identique à `aov`, G:T corrigée de la sphéricité (l) | l'interaction apparaît ; G est annoncé inter-sujet |
| b | Modèle mixte avec `verbose = FALSE` | plantage | résultat | — |
| c | Confirmation de Levene par un test de **moyennes** (`bf.test`) | n = 5/5/40, sd = 4/4/1 sous H0 : 38,5 % de faux positifs | 18 % (puis 10,5 % avec k) | Brown-Forsythe sur les variances ; l'ANOVA classique n'est plus retenue dans ce cas |
| d | `(1\|id)` analysé comme groupes indépendants | t1-t2 p = 1 | p = 0,024 (= `Error(id/T)`) | post-hocs appariés |
| e | `m.test(formula = y ~ F, data = d)` refusé | erreur | accepté | — |
| f | Messages toujours en français | hors locale française : tout en français | en anglais | rien en locale française |
| g | MANOVA non normale | plus de 6 min (MCD, 3000 tirages) | moins d'1 s ; 4,9 à 5,5 % sous H0 | `$method` = « Robust MANOVA (Wilks, ranks) » |
| h | Mesures répétées, lignes non triées par sujet | post-hocs entre sujets différents (p = 1 au lieu de 0,024) ; `m.test(x, g, id)` k > 2 plantait | identique aux données triées | post-hocs justes quel que soit l'ordre des lignes |
| i | Post-hocs du split-plot appariant des sujets de groupes différents | toutes les comparaisons à p ≈ 1 | `emmeans` sur `lmer` (B t1-B t3 p = 0,018) | nouvelle étape « comparaisons par emmeans » |
| j | Lettres d'une interaction en modèle mixte | « A » trois fois | « A t1 » | cellules nommées |
| k | Kruskal-Wallis sous hétéroscédasticité confirmée | 21 % (calcul direct) ; hetero_c 18 % | F de Welch : hetero_c 10,5 % | Welch au lieu de Kruskal-Wallis dans ce cas |
| l | Mesures répétées : p-value jamais corrigée de la sphéricité | voir § l | voir § l | p corrigée GG/HF toujours ; Mauchly en information |
| m | Welch annoncé « robuste » | récit trop affirmatif | « le choix le moins fragile », avec un avertissement pour les petits groupes asymétriques | message plus prudent |
| n | Résultat RM sans mention de la correction | « (p=0.22) » | « (p=0.22, corrigée de Greenhouse-Geisser) » | lien explicite avec l'étape de sphéricité |
| o | Plan en tranches déséquilibré (5 + 9 sujets) | plantage « Les variables explicatives… ne correspondent pas aux colonnes de g » | ANOVA à mesures répétées `aov()` + `Error()` ; grp identique à `aov` | le cas banal de TP aboutit ; la sortie signale les groupes inégaux |
| p | F corrigés de type III, tableau `aov` séquentiel | 5 + 9 sujets : « non corrigée » = type III (temps F = 3,17, p = 0,036) alors que le tableau affiche F = 5,14 (p = 0,0046) | ε appliqué aux F et ddl du tableau `aov` affiché | p non corrigée et p corrigée cohérentes avec le tableau |
| r | p globale d'un modèle à plusieurs termes = p du premier facteur (décision d'Antoine, 12 novembre) | plan en tranches 7 + 7 : 0,104 (grp) alors que temps p = 0,042 | min(1, 3 × 0,0421) = 0,126 (Bonferroni sur la plus petite p) | `$global_pvalue` et `return = FALSE` changent dès qu'il y a plusieurs termes ; le récit donne le calcul |
| q | Risque du Welch mal chiffré | « environ 10 à 17 % » | « environ 11 à 17 % », sources citées dans le code | chiffre exact |

**Format des nombres.** m.test écrit partout les p-values et les epsilons avec un point décimal (`.format_pval()`, `round()`), y compris en français. Les nouveaux messages gardent le point, par cohérence avec le reste de la sortie. m.bay utilise la virgule : à harmoniser plus tard si Antoine le souhaite.

**Points ouverts** (§ 4) : Welch par défaut sans test préalable de variance (décision d'Antoine) ; lenteur de `return = FALSE` ; mesure sous H0 après la correction de sphéricité et `R CMD check` final, en attente de mémoire libre.

## 1. Défauts signalés par les audits de m.bay (a à g)

| | Défaut | Reproduit ? | Avant | Après | Commit |
|---|---|---|---|---|---|
| a | Split-plot `y ~ G*T + Error(id/T)` : interaction retirée sans message | Oui | `aov` : G:T p = 0,030 ; `$effects` : G et T seulement (lmer `y ~ G + T + (1\|id)`) | `$effects` = G, T, G:T, p-values identiques à `aov` (écart < 1e-6) ; depuis l, T et G:T sont corrigées de la sphéricité | a9ddfcb |
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

Correction : un facteur constant pour chaque sujet est déclaré inter-sujet, et le message le dit. La partie fixe de `lmer` reprend les termes de la formule. Les deux contrôles raisonnent sujet par sujet. Le plan équilibré part maintenant vers `aov()` avec `Error()`, et `$effects` reproduit `aov` exactement. Un plan dont la structure intra-sujet est incomplète (sujet sans toutes les modalités, doublons) part vers `lmer`, qui garde G:T (test b). **Correction de cette phrase après la passe 2 de l'expert** : un plan complet mais aux groupes inter-sujets inégaux (5 + 9 sujets) ne partait pas vers `lmer`, il plantait (correctif o).

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

Décision : quand l'hétéroscédasticité est confirmée (Levene puis Brown-Forsythe) sur des données non normales, avec k > 2, une asymétrie modérée (max |skew| < 2), sans outlier extrême ni groupe constant, `m.test` applique le F de Welch, suivi de post-hocs t sans variance commune (comme la branche normale hétéroscédastique). C'est le mieux calibré des quatre (Delacre et al., 2019).

**Ce n'est pas un test robuste.** Avec de petits groupes asymétriques, le Welch reste fragile :
- ligne « exponentielles centrées » du tableau : 17,2 % de rejets sous H0 ;
- même plan, à l'intérieur de `m.test` (scénario hetero_c2, section 3) : 14,4 % ;
- relecture de l'expert indépendant, non refaite ici : avec 5 valeurs par groupe et une asymétrie de 1,4 à 1,75, Welch rejette 11 % des fois et Kruskal-Wallis 27 %.

Le récit ne dit donc plus « robuste à une non-normalité modérée ». Il dit « le choix le moins fragile ici » et avertit que, avec environ 5 valeurs par groupe asymétriques, le risque de faux positif reste au-dessus d'alpha (commit du correctif m, ci-dessous).

**Précision sur le message du commit 88f493d.** La parenthèse « 17 / 21.5 / 33 % with skewed data » se lit ainsi : plan n = 5/5/40, sd = 4/4/1, données exponentielles centrées, N = 1000 ; F de Welch 17,2 %, `t1way` 21,5 %, Kruskal-Wallis 33,3 %. Ce sont les trois tests de la ligne « exponentielles centrées » du tableau ci-dessus. Les autres cas gardent l'arbre existant (`med1way`, `t1way`, Kruskal-Wallis). Le premier essai échouait avec un groupe constant (test existant « un groupe constant ne fait pas planter l'analyse ») : la bascule exige maintenant des variances toutes positives (dc74ae7).

### l. Sphéricité jamais corrigée (signalé par l'expert indépendant)

L'ANOVA à mesures répétées (`aov()` + `Error()`) rendait toujours la p-value du F **non corrigé**. Quand Mauchly rejetait, la sortie décrivait les corrections de Greenhouse-Geisser et de Huynh-Feldt, puis écrivait « Les résultats ANOVA seront rapportés avec les p-values standard non corrigées ». Pire, cette étape ne pouvait pas se déclencher : les drapeaux étaient posés par `<<-` à l'intérieur d'un `tryCatch()`, ce qui écrit hors de la fonction. Avec peu de sujets, Mauchly est peu puissant, et la sortie affirmait souvent « sphéricité RESPECTÉE, aucune correction nécessaire ».

Correction (a00440c) :
- les termes intra-sujet reçoivent **toujours** une p-value corrigée : Greenhouse-Geisser si son epsilon est ≤ 0,75, Huynh-Feldt (plafonné à 1) sinon (Girden, 1992) ;
- les epsilons sont recalculés exactement par `car::Anova()` sur le plan de `rstatix` (qui les arrondit à 3 décimales : écart de 1e-4 sur p sinon) ;
- la p-value corrigée alimente le récit, `$effects` et la p-value globale ;
- Mauchly reste affiché comme information, avec la mention de sa faible puissance ;
- le test a (split-plot) compare maintenant G:T à la p-value corrigée ;
- depuis le correctif p (709080d), ε est appliqué aux F et ddl du **tableau `aov` affiché** (séquentiel), et non plus au F de type III de `car`. En plan équilibré, c'est le même F ; en plan déséquilibré, la p « non corrigée » annoncée est enfin celle du tableau.

Sens de la correction : la p corrigée est en général plus grande que la p non corrigée, mais pas toujours. Quand F < 1, réduire les degrés de liberté la diminue : dans l'exemple ci-dessous, 0,731 → 0,609.

Complément (5d4f7c2, correctif n) : la ligne de résultat le rappelle, « Aucun effet significatif de T (p=0.609, corrigée de Greenhouse-Geisser) », pour que l'étudiant fasse le lien avec l'étape de sphéricité.

L'expert indépendant a validé a00440c sur ses données (12 sujets × 6 temps, AR(1) ρ = 0,9) : ε GG = 0,361, identique à son calcul à la main ; p corrigée = 0,215, identique à `rstatix` p[GG] ; non corrigée 0,158.

Exemple du test l (12 sujets × 6 temps, AR(1) ρ = 0,9, graine 3) : ε GG = 0,481, Mauchly p = 0,019. F non corrigé p = 0,731 ; Greenhouse-Geisser p = 0,609. La référence est recalculée indépendamment (contrastes orthonormés, formule de l'epsilon).

Calibration sous H0 (12 sujets × 6 temps, 1000 jeux par scénario, `sim_sph.R`) :

| Mesure | Avant (master + a–k) | Après (a00440c) |
|---|---|---|
| AR(1) ρ = 0,9 : `m.test`, 1000 jeux | **8,1 %** (F non corrigé d'`aov` : 9,7 %) | **non mesuré ici**, voir ci-dessous |
| Symétrie composée (sphéricité vraie) : contrôle | non mesuré | non mesuré |
| Expert indépendant, AR(1) ρ = 0,9, 2000 jeux | F non corrigé : 8,7 % | Greenhouse-Geisser : 4,5 % |

**Mesure « après » non faite.** La simulation (`sim_sph.R`, 1000 jeux par scénario) a été interrompue deux fois : d'abord par la fin d'une session, puis par la saturation de la machine (0 Go libre sur 7,9 Go, 35 processus R). La règle commune aux sessions interdit ensuite toute simulation de plus d'une minute tant que la mémoire libre reste sous 1 Go. L'effet attendu est celui que l'expert a mesuré (8,7 % → 4,5 %), car la p-value de `m.test` est maintenant exactement la p-value de Greenhouse-Geisser ou de Huynh-Feldt (test l, et validation de l'expert sur ses données : p = 0,215 = `rstatix` p[GG]). **À relancer** sur un seul cœur quand la machine sera libre :

```
Rscript dev/mtest_audit/sim_sph.R <worktree> 1000 sph_after.rds 1 ar1 cs
```

(L'écart entre `m.test` à 8,1 % et le F d'`aov` à 9,7 % avant correction vient des jeux que l'arbre envoyait sur une autre voie, par exemple non paramétrique.)

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
- **`$global_pvalue` d'un plan à plusieurs facteurs : tranché par Antoine le 12 novembre 2026, voir § 4.3.** Règle antérieure, pour mémoire : C'est la p-value du **premier facteur principal de la formule** (`.first_factor_pvalue()`, `R/sys_effects_table.R` ; règle BP-031, fixée par le test « la p-value globale est celle du premier facteur »). Pour une mesure répétée, c'est la p-value corrigée de la sphéricité si ce facteur est intra-sujet. Exemples : split-plot de la section a, 0,199 (G) alors que G:T p = 0,030 ; plan en tranches équilibré de l'expert, 0,104 (grp) alors que temps p = 0,042. La règle ne décrit pas le plan entier ; l'expert propose le minimum corrigé ou le test du modèle entier. **Soumis à Antoine.** `$effects` donne tous les termes.
- **Doublon en mesures répétées** : « ASSOMPTION DE BASE : Indépendance… » puis « ASSOMPTION 1/3 : Indépendance… ».
- **Note emmeans** « adjust = "tukey" was changed to "sidak" » imprimée même avec `verbose = FALSE`.
- **Documentation MANOVA** : elle annonce une « Permutation MANOVA » et un « Kruskal-Wallis multivarié » qui n'existent pas dans le code.
- **roxygen** : `mbay-lib` contient roxygen2 8.1.1, qui réécrit tout le NAMESPACE. Le `.Rd` de `.manova_analysis` a donc été corrigé à la main.

### 4.3 Décision d'Antoine : p globale d'un modèle à plusieurs termes (correctif r, b824ef0)

**Décision d'Antoine Massé, 12 novembre 2026.** `$global_pvalue`, et la valeur rendue avec `return = FALSE`, servent à établir un lien entre variables, par exemple dans un réseau. Avec k termes testés (effets principaux et interactions), la p globale est `min(1, k × min(p))` : la correction de Bonferroni sur la plus petite p. Le taux de faux liens reste ainsi à alpha. Les p sont celles de `$effects`, corrigées de la sphéricité quand elle s'applique. Avec un seul terme, rien ne change : c'est l'appel de `corrigraph` (`m.test(x, g, return = FALSE)`, `R/corrigraph.R:456`).

Mise en œuvre, en un seul point (`R/m.test.R`, juste avant les post-hocs ; assistant `.global_pvalue_bonferroni()` dans `R/sys_effects_table.R`) : toutes les voies passent par là (ANOVA à plusieurs facteurs, mesures répétées, plan en tranches, modèle mixte, voies robustes à plusieurs termes, ANCOVA), avec `return = FALSE` comme avec `return = TRUE`. Le récit l'écrit : « p globale = 3 × 0.042 = 0.126 (plus petite p des 3 termes : temps, corrigée de Bonferroni) ».

Deux choix d'application, à valider par Antoine :
- **ANCOVA.** Un terme dont toutes les variables sont numériques (la covariable) est un ajustement, pas une comparaison de groupes : il n'est pas compté dans k. Sinon, une covariable très significative ferait passer tout modèle pour un « lien » avec le facteur. `y ~ x + F` garde donc la p de F, comme avant. Les interactions facteur × covariable sont comptées.
- **MANOVA.** À plusieurs termes, la règle porte sur la p par terme du test retenu (Wilks, ou Pillai si une assomption est violée). La MANOVA robuste (Wilks sur les rangs) teste les groupes croisés d'un seul bloc : un seul test, la règle est sans objet.

Mesures : plan en tranches équilibré de l'expert (7 + 7 sujets), p globale 0,104 (grp) → min(1, 3 × 0,0421) = 0,126. Le test r échoue sur le commit précédent (3 échecs sur 5) et passe après. Les tests qui figeaient l'ancienne règle (`test-m.test.R` : deux facteurs et réponse transformée ; `test-manova.R` : deux facteurs) suivent la nouvelle.

### 4.4 Durée d'un appel en plan en tranches (passe 2 de l'expert)

L'expert a vu plus de 10 minutes pour deux appels avec `boot = FALSE`. Mesuré séparément (`verbose = FALSE`, `return = TRUE`, machine encore chargée, 0,9 Go libre) : plan équilibré 7 + 7 sujets, 4,7 s puis 11,8 s ; plan 5 + 9 sujets après le correctif o, 2,9 s puis 2,3 s (avant : échec en 1,2 s). Profil de l'appel équilibré (9 s) : `.multi_factor_analysis()` 4,1 s, post-hocs `emmeans` du split-plot 1,8 s, sphéricité (`rstatix` + `car`) 0,3 s. Rien de pathologique : les 10 minutes viennent de la saturation de la machine (0 Go libre, échange disque) au moment de la mesure. L'estimation « environ 1 s » de l'introduction vaut pour un facteur à 50 observations ; un plan en tranches avec post-hocs `emmeans` coûte quelques secondes.

## 5. Vérifications

- `tests/testthat/test-m.test-audit.R` : 19 tests, 58 assertions, 0 échec au dernier commit (b824ef0). Chacun a été lancé sur le commit précédant sa correction (échec), puis après (succès).
- Non-régression au dernier commit : `test-m.test.R` (20 tests), `test-manova.R` (22 tests, 17 s au lieu de plus de 6 min) et `test-valreg.R` (6 tests) : 0 échec, 0 erreur (b824ef0).
- `R CMD check --no-manual --ignore-vignettes` (R 4.6.1, bibliothèque `mbay-lib`). Sur dc74ae7 : 2 NOTE (le dossier `.claude` inclus dans le build ; NEWS.md sans titre de version), corrigées au commit suivant. Sur 2a1e8b1 : Status: OK. Au dernier commit (l, m, n) : **non relancé**. Le check du 12 novembre a été tué par la saturation mémoire (« fatal error status code 4 » au chargement, sans message). La règle commune impose d'attendre que le verrou `VERROU_CHECK.txt` soit libre et que la mémoire libre dépasse 1 Go (0,08 Go au moment d'arrêter). À relancer.

## Références

- Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of variances. *Journal of the American Statistical Association*, 69(346), 364–367.
- Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by default use Welch's t-test instead of Student's t-test. *International Review of Social Psychology*, 30(1), 92–101.
- Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology*, 32(1), 13.
- Girden, E. R. (1992). *ANOVA: Repeated measures*. Sage University Paper Series on Quantitative Applications in the Social Sciences, 84. Newbury Park, CA: Sage.
- Nath, R., & Pavur, R. (1985). A new statistic in the one-way multivariate analysis of variance. *Computational Statistics & Data Analysis*, 2(4), 297–315.
- Zimmerman, D. W. (2004). A note on preliminary tests of equality of variances. *British Journal of Mathematical and Statistical Psychology*, 57(1), 173–181.
