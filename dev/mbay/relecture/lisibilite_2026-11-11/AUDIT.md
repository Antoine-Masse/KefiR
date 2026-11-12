# Audit de lisibilité de m.bay() — 11 novembre 2026

Lecteur visé : étudiant de BUT génie biologique (ANOVA, khi-deux, p-values ; jamais de bayésien).
Branche `feat/m.bay` (dernier commit c8e8976). Script : `run_cotes.R` ; sorties : `sorties_<plan>.txt`.

## Plans lancés

| Plan (fichier) | Appel | m.test | m.bay (durée) | Lignes m.bay / m.test |
|---|---|---|---|---|
| g2 | `y ~ g`, 2 × 10 | Student p = 0,13 | indécis, 54,9 % pour « aucune différence » (4,7 s) | 80 / 25 |
| g4 | `y ~ souche`, 4 × 8 | ANOVA p = 0,00094 | 97,6 % (13 s) | 92 / 35 |
| f2x3 | `y ~ milieu * temp`, effet de temp seul | temp p = 1e-07 | temp > 99,9 %, milieu non tranché (9,7 s) | 110 / 52 |
| f2x2int | `y ~ engrais * lumiere`, interaction | 3 termes significatifs | 3 termes « oui » (6,2 s) | 118 / 51 |
| ancova | `poids ~ regime + poids0` | regime p = 0,0006 | 97,8 % (1,3 s) | 84 / 116 |
| apparie | `y ~ t \| id`, 12 sujets | t apparié p = 0,00053 | 98,3 % (0,1 s) | 70 / 35 |
| repet3 | `y ~ t \| id`, 10 sujets × 3 temps | Friedman p = 0,061 | 86,9 % (2,4 s) | 81 / 60 |
| tranches | `y ~ groupe * temps + Error(id/temps)` | lmer additif, temps p = 0,0011 | 3 termes « oui » (43,2 s) | 159 / 108 |
| ouinon1 | `germe ~ lot` (oui/non, 3 lots) | **échec** (`quantile()` sur facteur) | non tranché, 79,5 % pour H0 (0,1 s) | 77 / — |
| comptages | `colonies ~ boite` (Poisson) | Welch p = 0,0065 | 95,4 % (8,8 s) | 88 / 34 |
| notes15 | `note ~ recette`, note = facteur **ordonné** 1–5 | **échec** (`quantile()` sur facteur ordonné) | 99,7 %, traité en **multinomial** (0,5 s) | 78 / — |
| notes15num | même jeu, note en **entiers** 1–5 (ajouté) | ANOVA p = 0,00013 | 99,4 %, modèle **ordinal** (2,5 s) | 66 / 52 |
| ouinon2x2 | `pousse ~ souche * milieu` | Scheirer-Ray-Hare puis **échec** (`x must be numeric`) | souche 86,8 % (0,2 s) | 67 / — |
| manova | `cbind(y1, y2) ~ g` | Wilks p = 3e-05 | > 99,9 % (3,5 s) | 69 / 57 |

Aucun appel n'a dépassé 60 s : `sensitivity` est resté par défaut partout.
Remarque : la première version du jeu `tranches` n'avait pas d'effet (erreur de simulation) ; elle a été corrigée et relancée, le fichier contient la version avec effet.

---

## 1. Cohérence du vocabulaire d'un plan à l'autre

### 1.1 L'état « on ne sait pas » a cinq noms
- Titre : « 7) Bilan - indécis (Aucune différence de y selon g : 54,9 %). » (g2) ; même mot dans ouinon1 et dans la première version de tranches.
- Verdict : « ==> Verdict : indécis. » (g2), mais « ==> Verdict : plutôt pas de différence, sans conclure. » (ouinon1, tranches v1).
- Ligne d'hypothèse : « 1. Aucune différence de y selon g : 54,9 % (non tranché) » (g2).
- Termes et paires : « milieu non tranché (53,8 %) » (f2x3), « Non tranchés : A/B, A/C, B/C. » (g4).
- Compte rendu : « Les données ne permettent pas de conclure » (g2, ouinon1).
Un étudiant ne sait pas si « indécis », « non tranché » et « sans conclure » sont la même chose. **Proposition : un seul mot, « non tranché », partout** (c'est déjà le plus fréquent et il est défini par la note des seuils).

### 1.2 L'écart-type de référence a six noms
- « soit 0,49 écart-type » (g2, g4, comptages, f2x2int — alors que f2x2int est un plan à 2 facteurs) ;
- « soit 2,37 écart-type résiduel (autres termes pris en compte) » (f2x3, étape 6) puis « Écart-type au sein des cases : 1,11 » (f2x3, bilan) : **deux noms pour la même grandeur dans la même sortie** ;
- « écart-type résiduel (covariables prises en compte) » (ancova) ;
- « écart-type résiduel au sein des sujets » (repet3, tranches) ;
- « soit -1,18 écart-type des différences (dz) » (apparie).
De plus, l'accord manque : « 2,37 écart-type », « 1,70 écart-type », « 3,69 écart-type résiduel » → « écarts-types ». Proposition : « écarts-types résiduels » partout dans l'étape 6, une seule ligne de définition au bilan.

### 1.3 Précision du calcul : quatre formulations
« Précision du calcul : environ 1,0 % sur le facteur de Bayes, négligeable. » (g2) ; « Précision du calcul : < 0,1 point sur les probabilités, négligeable. » (ancova) ; « Précision du calcul : < 0,1 point sur la probabilité. » (manova, sans « négligeable ») ; « Précision : calcul exact, sans tirage. » (apparie, ouinon1) ; « Précision : forme close, sans tirage. » (notes15num, ouinon2x2). « Forme close » est du jargon ; « calcul exact, sans tirage » dit la même chose.

### 1.4 « L'ordre ne change pas » : quatre formulations à l'étape 5
« --> L'ordre ne change pas ; seul le degré de soutien peut changer. » (g4, comptages) ; « --> La réponse ne change pas de côté. » (ancova, apparie, repet3, notes15) ; « --> L'ordre ne change pas. » (manova, notes15num) ; « --> Aucun terme ne change de côté. » (plans factoriels). Et deux formes de la ligne Sensibilité : « Sensibilité : « y diffère selon souche » va de 92,8 % à 98,9 % » (un facteur) contre « Sensibilité : la réponse en tête va de 93,1 % à 99,0 % » (ancova, apparie…). La première, qui nomme l'hypothèse, est plus claire.

### 1.5 Noms d'étapes
- Étape 3 : « Formes du bruit » (un facteur continu), « Forme du bruit » (autres plans continus), « Loi des données » (ouinon1, notes15, notes15num), « Loi des réponses » (ouinon2x2). « Loi des données » / « Loi des réponses » ne diffèrent que par le mot ; garder « Loi des données ».
- Étape 6 : « Quels groupes… » / « Quels niveaux… » / « Quelles conditions… » / « Quelles réponses diffèrent » / « De combien les deux niveaux diffèrent » / « De combien la part change » : justifié par le plan, acceptable.
- Référence de fonction placée à l'étape 4 dans un seul plan : « 4) Point de départ (a priori) [ttestBF(paired = TRUE) {BayesFactor}] » (apparie) ; ailleurs l'étape 4 n'a pas de crochet.

### 1.6 Notation « A/B » utilisée pour deux choses
« R1/R2 : différents à 41,2 % ; R2/R1 = 3,14 [0,80 ; 12,31]. » (notes15num) : la même barre désigne la paire puis un rapport de cotes, avec l'ordre inversé. Idem « souche : S2/S1 = 3,57 » (ouinon2x2).

### 1.7 Unités et effectifs
- « Moyennes observées, groupe : actif = 31,9 ; placebo = 30,7 (8 sujets x 3 mesures). » puis « temps : … (n = 16). » (tranches) : deux styles dans deux lignes voisines.
- ouinon1 : l'étape 6 parle de « écart de la part de oui, en points » mais le bilan dit « Écart le plus marqué : non, 26,7 % en L3 contre 43,3 % en L1. » (bascule sur « non »).
- notes15 (facteur ordonné) : « ==> R1 : 0 5, 1 4, 5 3, 9 2, 0 1 (sur 15). » est illisible quand les modalités sont des chiffres (lire « 0 fois la note 5, 1 fois la note 4… ») et l'ordre est décroissant. Au bilan : « 0 5 (< 0,1 %) » pour un effectif nul (→ « 0 % »), et « Écart le plus marqué : 2, 6,7 % en R3 contre 60,0 % en R1. » (« 2 » est la modalité, pas un nombre).

### 1.8 Échelle de soutien appliquée de deux façons
La note rappelle « facteur de Bayes ; 1 à 3 faible, 3 à 20 modéré, 20 à 150 fort ». Mais :
- tranches : « groupe:temps … : oui, soutien faible (66,0 % ; a priori 16,7 % ; facteur de Bayes 9,72) » — 9,72 est « modéré » selon la note ;
- f2x2int : le global « 97,0 % (soutien modéré) » (facteur de Bayes 16,1) est plus faible que chacun de ses termes « oui, soutien fort (96,2 % … 25,1) » ;
- comptages : 95,4 % est « soutien fort » (20,6), f2x2int : 97,0 % est « soutien modéré » (16,1).
L'étudiant lit le pourcentage, pas le facteur de Bayes : un 97 % « modéré » sous un 95 % « fort » le déroute.

### 1.9 Jargon et références internes laissés visibles
- « Note : queues lourdes avec des sujets : mesurées sans faux « oui » fort (loi normale gardée, D28). » (apparie, repet3, tranches) : **D28 est un numéro de décision interne**.
- « r va de 0,25 à 1 (rscaleCont avec lui) » (ancova) ; « r va de r/2 à 2r » (manova, alors que les autres plans écrivent « 0,25 à 1 ») ; « (HOP, termes indépendants, uniforme) » (ouinon2x2) ; « approximation de Laplace du modèle additif » (ouinon2x2) ; « epsilon de Greenhouse-Geisser 0,91 » (repet3, tranches) ; « sur dz » (apparie) ; « Dirichlet-multinomiale {KefiR} » en titre d'étape 6 (ouinon1, notes15).

---

## 2. Longueur

Sur les plans simples, m.bay est 2 à 3 fois plus long que m.test : g2 80 lignes contre 25, apparie 70/35, comptages 88/34, g4 92/35. Le plan en tranches atteint 159 lignes.

Coupables, du plus rentable au moins rentable :
1. **tranches, étape 6** : 15 paires de cases sur 3 lignes chacune (45 lignes), dont 9 « non tranchées », puis la liste complète redite au bilan (« Non tranchés : actif x T1/actif x T2, actif x T1/placebo x T1, … » sur 4 lignes). Les chaînes de regroupement (« actif x T1 = placebo x T2 = placebo x T1 = placebo x T3 = actif x T2 < actif x T3 ») débordent la largeur 80 et sont imprimées trois fois.
2. **Regroupement en double** (f2x3, f2x2int, tranches) : « ==> Regroupement en tête : 20C < 30C < 37C (77,6 %). » est suivi immédiatement de « ==> Si une différence existe, regroupements en tête : 20C < 30C < 37C (77,6 %) ; … ». La première ligne est redondante.
3. **Écart redit au bilan** : « ==> Écart estimé traite - temoin : 0,8 [-0,5 ; 2,2] (étape 6). » (g2, apparie, f2x3) reprend mot pour mot l'étape 6, quelques lignes plus haut.
4. **Étape 3 à un facteur** (g2, g4, comptages) : 4 lignes de probabilités des formes, une ligne « Variances », une ligne « La réponse porte surtout sur les moyennes », puis deux notes de 3 lignes. Pour l'étudiant, l'utile tient en deux lignes : la forme en tête et « variances égales / inégales ». La liste complète des quatre formes peut aller dans `detail = TRUE`.
5. **Sensibilité dite deux fois** : à l'étape 5 (« Sensibilité : … va de 92,8 % à 98,9 % ») et au bilan (« Selon le point de départ : 92,8 % à 98,9 % ; l'ordre ne change pas. »). Garder celle du bilan, réduire l'étape 5 à son titre quand l'ordre ne change pas.
6. **Notes répétées dans chaque sortie** : « Note : ces valeurs ne sont pas ajustées sur les données ; l'étape 5 les fait varier. » et « Note : l'écart estimé est ramené vers 0 par l'a priori, surtout quand les groupes sont petits. » (cette dernière dite pour des « groupes » même en apparié). Utiles une fois ; acceptables.
7. **Étape 1 des factoriels** : « 5 combinaisons de termes sont possibles (une interaction n'entre qu'avec ses facteurs : hérédité). » + « Une réponse par terme, obtenue en pesant ensemble ces 5 combinaisons. » : 4 lignes pour une information de méthode, à mettre dans `detail`.

---

## 3. Contradictions internes

1. **« les autres termes sont plus bas » alors qu'il n'y en a pas** — f2x2int : « ==> Verdict : une différence selon engrais, lumiere, engrais:lumiere ; les autres termes sont plus bas. » ; tranches : « une différence selon groupe, temps, groupe:temps ; les autres termes sont plus bas. » Dans f2x3 la même formule (« une différence selon temp ; les autres termes sont plus bas ») cache que milieu est « non tranché » et non « non ».
2. **« oui » à 66 % pour un terme, « non tranché » à 75 % près pour une paire** — tranches : titre « groupe:temps oui (66,0 %) », compte rendu « groupe:temps oui (66,0 %) », alors que la note du même bilan dit « différents à partir de 75 % de probabilité … sinon non tranchés ». Et ouinon1 : 79,5 % pour « aucune différence » est « non tranché ». Le seuil d'un « oui » change selon le niveau de lecture, sans que la sortie le dise.
3. **Toutes les paires « semblables » mais un regroupement avec une différence** — ouinon1 : « ==> Semblables : L1/L2, L1/L3, L2/L3 (modéré). » puis, deux lignes après, « regroupement le plus probable : L1 = L2 < L3 (46,2 %) », alors que l'étape 6 donne « L3 - L1 = 15,9 points [-6,9 ; 39,1] ». Le « si une différence existe » est là, mais l'étudiant lira « semblables » puis « < ».
4. **Titre d'étape 6 qui affirme un ordre pour un terme non tranché** — f2x3 : « 6) … - temp : 20C < 30C < 37C ; milieu : M2 > M1 (0,50). » alors que milieu est « non tranché (53,8 %) » et que l'étape 5 dit « milieu change de côté : 25,2 % à 62,2 % ».
5. **Global « modéré » sous des termes « fort »** — f2x2int : « 1. y diffère selon au moins un facteur : 97,0 % (soutien modéré) » ; « engrais … oui, soutien fort (96,2 % …) ». Voir 1.8.
6. **Hypothèse 2 qui oublie l'interaction** — f2x2int : « 2. Aucune différence de y selon engrais, lumiere : 3,0 % » (idem tranches, ouinon2x2). Correct (l'interaction n'existe pas sans ses facteurs) mais pas évident pour l'étudiant ; « selon aucun terme » serait plus sûr.
7. **Verdict nuancé, titre brut** — ouinon1 : titre « indécis », verdict « plutôt pas de différence, sans conclure », hypothèse « (non tranché) », compte rendu « Les données ne permettent pas de conclure ». Pas faux, mais quatre façons pour un même état (cf. 1.1).

Aucune contradiction trouvée entre la ligne « Verdict » et la phrase « pour un compte rendu » dans les plans à un facteur (g2, g4, comptages, apparie, repet3, ancova, manova) : ces bilans sont propres.

---

## 4. Comparaison avec m.test

| Plan | m.test | m.bay | Accord ? m.bay aide-t-il ? |
|---|---|---|---|
| g2 | p = 0,13, non significatif | non tranché, 45,1 % pour une différence | Accord. m.bay explique bien (« proche de 1 : les données départagent à peine les deux »). |
| g4, comptages, apparie, f2x2int, manova, notes15num | significatif | 95 à > 99,9 % | Accord. |
| f2x3 | milieu p = 0,16 | milieu non tranché 53,8 % | Accord sur le fond, mais le titre d'étape 6 « milieu : M2 > M1 » peut faire croire à un effet (contradiction 3.4). |
| ancova | moyennes ajustées 271,6 / 281,6 / 285,7, R1-R2 et R1-R3 significatifs | R1/R2 92,6 %, R1/R3 96,7 % | Accord, **mais** m.bay affiche « Moyennes brutes, regime : R1 = 274 ; R2 = 274 » à côté de « R2 - R1 = 8 ». La ligne « poids0 moyen par groupe » et « comparés à poids0 égal » expliquent, mais seul m.test montre les moyennes ajustées. |
| repet3 | **Friedman p = 0,061, non significatif** | **86,9 %, J14 différent** | **Désaccord non expliqué.** m.test conclut à une « non-normalité sévère » (Skewness = -3, Kurtosis = 13) ; m.bay dit « Aucune valeur éloignée ne met d'emblée la loi normale en doute ». Vérification hors paquet : résidus de `lm(y ~ t + id)` normaux (Shapiro p = 0,36), ANOVA à mesures répétées classique p = 0,0066. **C'est m.test qui se trompe** (diagnostic de normalité sur les résidus du modèle à strates) ; m.bay a raison mais ne peut pas le dire. À signaler côté m.test. |
| tranches | lit `groupe` comme intra-sujet (« 48 (attendu 96) »), passe à un lmer additif **sans interaction** ; temps seul significatif | groupe 79,6 %, temps 98,2 %, groupe:temps 66,0 % « oui » | Désaccord dû à m.test (plan mal lu). m.bay lit bien le plan, mais son « oui » à 66 % pour l'interaction est plus affirmatif que ses propres seuils (3.2). |
| ouinon1, notes15 (facteur ordonné) | **échec** : `quantile()` refuse un facteur | répond | m.test inutilisable sur une réponse facteur. |
| ouinon2x2 | Scheirer-Ray-Hare sur un 0/1, souche p = 0,0028, puis **échec** au post-hoc | souche 86,8 % | Accord sur souche ; m.test est mal outillé pour un oui/non à deux facteurs. |
| notes15 vs notes15num | — | facteur ordonné → **multinomial**, R1/R2 « différents à 87,1 % » ; entiers → **ordinal**, R1/R2 « non tranché, 41,2 % » | Mêmes données, deux codages, deux conclusions de paires. L'étudiant qui code proprement ses notes en `ordered()` obtient le modèle qui **ignore l'ordre**, sans avertissement à l'étape 1 (« Réponse catégorielle, 5 modalités : 1, 2, 3, 4, 5. »). Question de méthode (détection automatique), hors des retouches ci-dessous, mais à regarder. |
| comptages | Welch | lois continues | Accord. m.bay ne dit pas que des comptages pourraient être lus avec `response = "counts"`. |

---

## 5. Retouches proposées (texte et présentation seulement), par importance

**1. Supprimer la référence interne D28** — plans apparie, repet3, tranches (et tout plan à sujets).
- Actuel : « Note : queues lourdes avec des sujets : mesurées sans faux « oui » fort (loi normale gardée, D28). »
- Proposé : « Note : la loi normale est gardée ; des valeurs éloignées ne créent pas ici de faux « oui » fort. »

**2. Un « oui » de terme exige les mêmes 75 % qu'une paire, sinon « plutôt oui »** (libellé seul, sans changer les probabilités) — tranches (et tout factoriel).
- Actuel : titre « groupe:temps oui (66,0 %) » ; ligne « groupe:temps (…) : oui, soutien faible (66,0 % ; a priori 16,7 % ; facteur de Bayes 9,72). » ; compte rendu « groupe:temps oui (66,0 %) ».
- Proposé : « groupe:temps plutôt oui (66,0 %) » partout, et la note des seuils du bilan dit une seule règle : « oui à partir de 75 %, non jusqu'à 25 %, sinon plutôt oui / plutôt non ».

**3. Verdict factoriel : nommer les termes restants au lieu de « les autres termes sont plus bas »** — f2x2int, tranches, f2x3, ouinon2x2.
- Actuel (f2x2int) : « Verdict : une différence selon engrais, lumiere, engrais:lumiere ; les autres termes sont plus bas. »
- Proposé : « Verdict : une différence selon engrais, lumiere et engrais:lumiere. »
- Actuel (f2x3) : « Verdict : une différence selon temp ; les autres termes sont plus bas. »
- Proposé : « Verdict : une différence selon temp ; milieu et milieu:temp non tranchés. »

**4. Un seul mot pour l'indécision : « non tranché »** — g2, ouinon1 (et tout plan à verdict non tranché).
- Actuel : « 7) Bilan - indécis (Aucune différence de y selon g : 54,9 %). » / « ==> Verdict : indécis. » / « ==> Verdict : plutôt pas de différence, sans conclure. »
- Proposé : « 7) Bilan - non tranché (Aucune différence de y selon g : 54,9 %). » / « ==> Verdict : non tranché. » / « ==> Verdict : non tranché, penche vers aucune différence. »

**5. Une seule échelle de soutien, celle de la note** — tranches, f2x2int, comptages.
- Actuel : « oui, soutien faible (66,0 % ; … facteur de Bayes 9,72) » sous une note qui dit « 3 à 20 modéré » ; « 97,0 % (soutien modéré) » au-dessus de termes « soutien fort ».
- Proposé : le mot de soutien suit toujours l'échelle imprimée (9,72 → « modéré ») ; et pour le global factoriel, ajouter « (facteur de Bayes global, plus prudent que celui de chaque terme) » après « soutien modéré ».

**6. ouinon1 : ne pas dire « semblables » à côté d'un regroupement avec « < »** — ouinon1 (tout un-facteur où H0 mène).
- Actuel : « ==> Semblables : L1/L2, L1/L3, L2/L3 (modéré). » puis « ==> Si une différence existe, regroupement le plus probable : L1 = L2 < L3 (46,2 %). »
- Proposé : « ==> Aucune paire n'est dite différente (L1/L2, L1/L3, L2/L3 semblables, modéré). » puis « ==> S'il existait une différence, ce serait plutôt L3 au-dessus : L1 = L2 < L3 (46,2 %). »

**7. Réponse facteur à modalités chiffrées : écrire la modalité, dans l'ordre croissant, et prévenir que l'ordre n'est pas utilisé** — notes15 (et tout oui/non ou catégoriel à modalités numériques).
- Actuel : « ==> R1 : 0 5, 1 4, 5 3, 9 2, 0 1 (sur 15). » ; « 0 5 (< 0,1 %) » ; « Écart le plus marqué : 2, 6,7 % en R3 contre 60,0 % en R1. » ; étape 1 « Réponse catégorielle, 5 modalités : 1, 2, 3, 4, 5. »
- Proposé : « ==> R1 : « 1 » 0, « 2 » 9, « 3 » 5, « 4 » 1, « 5 » 0 (sur 15). » ; « « 5 » 0 (0 %) » ; « Écart le plus marqué : modalité « 2 », 6,7 % en R3 contre 60,0 % en R1. » ; ajouter en étape 1 : « --> L'ordre des modalités n'est pas utilisé ; pour des notes, les saisir en nombres entiers (modèle ordinal). »

**8. Écart-type : un seul nom, au pluriel, défini une fois** — f2x3, f2x2int, g2, g4, comptages, ancova, repet3, tranches, apparie.
- Actuel : « soit 2,37 écart-type résiduel (autres termes pris en compte). » ; « Écart-type au sein des cases : 1,11 » ; « soit 0,49 écart-type. » ; « soit -1,18 écart-type des différences (dz). »
- Proposé : étape 6 « soit 2,37 écarts-types résiduels. » ; bilan « Écart-type résiduel (au sein des cases, autres termes pris en compte) : 1,11 » ; apparié « soit -1,18 écart-type des différences avant/après. »

**9. Ne pas imprimer deux fois le même regroupement ; résumer les paires non tranchées quand elles sont nombreuses** — f2x3, f2x2int, tranches.
- Actuel : « ==> Regroupement en tête : 20C < 30C < 37C (77,6 %). » suivi de « ==> Si une différence existe, regroupements en tête : 20C < 30C < 37C (77,6 %) ; … » ; tranches : 15 paires sur 45 lignes.
- Proposé : supprimer la ligne « Regroupement en tête : … » ; au-delà de 6 paires, ne détailler à l'étape 6 que les paires différentes ou semblables et écrire « ==> 9 autres paires non tranchées (écarts avec detail = TRUE). »

**10. Notations et jargon** — notes15num, ouinon2x2, f2x3, ancova, manova, apparie, repet3, tranches.
- « R2/R1 = 3,14 [0,80 ; 12,31] » → « rapport des cotes R2 contre R1 = 3,14 [0,80 ; 12,31] » ; « souche : S2/S1 = 3,57 » → « souche : rapport des cotes S2 contre S1 = 3,57 ».
- f2x3, titre d'étape 6 : « milieu : M2 > M1 (0,50) » → « milieu (non tranché) : M2 - M1 = 0,50 ».
- « r va de 0,25 à 1 (rscaleCont avec lui) » → « r va de 0,25 à 1 (pour le facteur et la covariable) » ; manova « r va de r/2 à 2r » → « r va de 0,25 à 1 » ; « (HOP, termes indépendants, uniforme) » → « (trois façons de répartir les chances entre modèles) » ; « Précision : forme close, sans tirage. » → « Précision : calcul exact, sans tirage. » ; « epsilon de Greenhouse-Geisser 0,91 (proches d'être égales) » → « corrélations entre temps proches d'être égales (epsilon 0,91, 1 = parfait) ».

---

## Hors périmètre (méthode ou m.test), pour mémoire
- notes15 : un facteur `ordered()` part en multinomial, des entiers 1–5 en ordinal ; conclusions de paires différentes sur les mêmes données.
- m.test : échoue sur toute réponse facteur (oui/non, notes en facteur) ; lit mal `Error(id/temps)` avec un facteur entre sujets ; diagnostic de normalité faux en mesures répétées (repet3).
- ancova : m.bay ne montre pas les moyennes ajustées, que m.test montre.
- manova : pas de comparaison des groupes deux à deux, que m.test donne (A, B, C tous différents).

---

## Suites données (session m.bay, 11 novembre)

**Faites.**
1. Référence interne « D28 » retirée du récit : « la loi normale est gardée ; des valeurs éloignées ne créent pas ici de faux « oui » fort ».
2. Une seule règle de seuils pour les termes et les paires : « oui » à partir de 75 %, « non » jusqu'à 25 %, « plutôt oui » ou « plutôt non » entre les deux. Les probabilités ne changent pas.
3. Verdict à plusieurs termes : les autres termes sont nommés, avec leur réponse (« Verdict : une différence selon temp ; milieu non, milieu:temp non tranché. »).
4. Un seul mot pour l'indécision : « non tranché » (titre du bilan et verdicts ; « non tranché, penche vers … » au lieu de « plutôt …, sans conclure »).
5. Le qualificatif de soutien suit l'échelle du facteur de Bayes imprimée dans le bilan (3 à 20 modéré, 20 à 150 fort, au-delà très fort), quel que soit l'a priori. Seule la note de seuils l'annonçait jusqu'ici.
7. Fond : une réponse saisie en facteur ordonné (`ordered()`) prend la voie ordinale, comme des entiers. Les rangs sont nommés à l'étape 3. Sur notes15, les deux codages donnent désormais la même réponse.
8. « écarts-types » au pluriel à partir de 2.
9. « Regroupement en tête » n'est plus répété à l'étape 6 à plusieurs facteurs : le titre de l'étape le donne.
10. « Précision : calcul exact, sans tirage » ; « trois façons de répartir les chances entre les modèles » ; MANOVA : « r va de 0,25 à 1 » au lieu de « de r/2 à 2r ».

**Non reprises ici.**
- 6. « Semblables » à côté d'un regroupement avec « < » : le regroupement est donné « si une différence existe », ce qui est juste.
- Fin du 8 (un seul nom de l'écart-type dans le bilan), fin du 9 (résumé des paires au-delà de 6) et fin du 10 (libellés des rapports de cotes, epsilon) : à faire au prochain lot.

**Hors périmètre, transmis à la session parente.**
- m.test échoue sur les réponses en facteur ; il lit mal `Error(id/temps)` avec un facteur entre sujets ; son diagnostic de normalité est faux en mesures répétées.
- m.bay ne montre pas les moyennes ajustées en ANCOVA, ni les paires de groupes en MANOVA.
