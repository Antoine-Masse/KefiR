# Audit de lisibilité : sortie verbose de m.bay() comparée à m.test()

Données simulées communes (`set.seed(42)`, script `donnees.R`), six plans, `verbose = TRUE`, `LANGUAGE = "fr"`, m.test sur master et m.bay sur feat/m.bay, dans deux processus R séparés. Les sorties brutes sont dans ce dossier (`m.test_*.txt`, `m.bay_*.txt`). Aucun fichier du dépôt n'a été modifié.

## Mesures

| Plan | m.test lignes / >100 car. / mots | m.bay lignes / >100 car. / mots | Ligne la plus longue (test / bay) | Longueur moyenne (test / bay) |
|---|---|---|---|---|
| 2 groupes | 22 / 2 / 152 | 37 / 20 / 773 | 120 / 316 | 49 / 119 |
| 4 groupes | 62 / 2 / 486 | 51 / 30 / 1135 | 120 / 316 | 57 / 117 |
| 2 facteurs croisés | 71 / 2 / 515 | 37 / 19 / 699 | 115 / 252 | 54 / 103 |
| Apparié (id) | 35 / 2 / 256 | 39 / 16 / 687 | 147 / 255 | 53 / 101 |
| Répété 3 temps | 69 / 3 / 583 | 47 / 27 / 934 | 109 / 255 | 61 / 111 |
| ANCOVA | 126 / 3 / 847 | 46 / 27 / 958 | 108 / 290 | 48 / 115 |

Dans une console de 80 colonnes, m.bay occupe environ 69 à 102 lignes à l'écran (m.test : 26 à 147). Le « Verdict » de m.bay arrive vers la ligne d'écran 37 à 66, soit après la moitié de la sortie, sans rien qui le distingue visuellement.

Occurrences de jargon dans les six sorties m.bay : « a priori » 31, « facteur(s) de Bayes » 21, « écart-type » 18, « intervalle de crédibilité » 10, « a posteriori » 9, « cotes » 6, « Rouder » 4, « HOP / Beta(1,1) / hérédité » 4, « Zellner-Siow » 3, « dz » 3, « échantillonnage d'importance » 2, « rscaleCont » 2, « symétrie composée » 1.

## Verdict global

Sur le fond, m.bay est plus clair que m.test : chaque étape dit ce qu'elle a trouvé et ce qu'elle en fait, la conclusion nomme un verdict, chiffre les écarts et dit ce qu'elle ne prouve pas. Mais il ne suit pas encore la forme « étape par étape » de m.test. Ses lignes sont deux fois plus longues (100 à 120 caractères en moyenne, contre environ 50 ; jusqu'à 316). Il compte 2 à 5 fois plus de mots (773 contre 152 pour deux groupes). Il imprime deux fois les mêmes résultats (étapes 6 et 7). Il répète à chaque appel les mêmes paragraphes pédagogiques. Et la conclusion est fondue dans le texte. Un étudiant de BUT lit m.test d'un coup d'œil, au prix de nombreux tests qu'il ne comprend pas. Devant m.bay, il trouve des phrases justes mais trop longues, et doit remonter l'écran pour chercher le verdict. Il faut donc gagner en forme, pas en contenu. Les trois chantiers prioritaires sont la longueur des lignes, la redondance entre les étapes 6 et 7, et la visibilité du bilan.

## Conventions de m.test que m.bay ne suit pas

1. **Le résultat est dans le titre de l'étape.** m.test : `4) Test de Fisher-Snedecor [var.test()] - Variances homogènes (p = 0.47).` Chez m.bay, le titre décrit l'étape, et le résultat n'apparaît que plus bas : `5) Facteur de Bayes et précision numérique. [m.bay() {KefiR}]`.
2. **Des lignes courtes, une idée par ligne.** m.test : `==> Les résidus sont non-normaux (p-value : 0.012).` m.bay enchaîne plusieurs idées par des points-virgules : `==> 50,0 % de chances qu'une différence existe (pile ou face : ...) ; taille attendue r = 0,5 (échelle « medium » de Rouder et al. 2012) ; chacune des 4 formes du bruit part à 25,0 %.` (220 caractères).
3. **Peu de `==>`, et `-->` réservé aux décisions.** m.test : `--> Conservation de l'ANOVA classique (situation paramétrique).` m.bay met jusqu'à 7 `==>` dans une même étape (étape 6, quatre groupes). Il emploie aussi `-->` pour des explications, qui ne décident rien : `--> Toutes les paires viennent d'une seule loi sur les regroupements : ...`, `--> Pourquoi : une proportion ~ Beta(1,1) par ordre d'interaction (HOP, ...)`.
4. **Les explications portent l'étiquette « Note : ».** m.test : `Note : les p-value sont comparées à un alpha corrigé (Sidak) de 0.025.` m.bay glisse des lignes sans marque au milieu des constats : `Les valeurs éloignées sont cherchées une fois retiré le niveau propre de chaque sujet.` (insérée entre deux `==>`), ou `Écart estimé : médiane a posteriori, ...`.
5. **Sous-points a), b), c) et second niveau de retrait.** m.test : `a) Outliers détectés sur les résidus :` puis `\t\t==> 1 outliers (1.7%)...`. m.bay aligne ses 6 paires comme six `==>` au même niveau, et sa conclusion empile trois niveaux de retrait (tab, tab + 2, tab + 5 espaces) sans marque de sous-section.
6. **Une séparation nette entre le récit et le résultat.** m.test ferme le récit par une ligne isolée (`m.test() - version 02 bêta 2026 - ...`), puis donne un tableau compact. Chez m.bay, la conclusion est l'étape 7, indentée comme les autres, et la sortie se termine par `Formes du bruit, a priori en détail, mode de lecture : print(resultat, detail = TRUE).`
7. **Le compteur de progression.** m.test : `ASSOMPTION 2/3`, `ASSOMPTION 4/5`. m.bay n'en a pas, et son squelette change selon le plan : avec un facteur, l'a priori est à l'étape 3 et la forme du bruit à l'étape 4 ; dans les autres plans, c'est l'inverse. Le plan à deux facteurs n'a que 6 étapes.
8. **Les crochets uniquement sur la ligne-titre, avec la fonction.** m.test : `[shapiro.test()]`, `[HSD.test() of {agricolae}]`. m.bay en met aussi sur des lignes de constat (`... 0 au-delà de 3. [Tukey, quantile() {stats}]`, `... de 48,0 % à 90,9 %. [m.bay() {KefiR}]`). Il y place aussi un nom de méthode qui n'est pas une fonction (`[m.bay() {KefiR}, échantillonnage d'importance]`).

## Recommandations (classées, sans aucun changement de calcul)

### 1. Rendre le bilan visible en 10 secondes
Il faut sortir la conclusion du flux indenté, avec une ligne vide et un titre au bord gauche, et mettre le verdict en première ligne, à côté de la probabilité de l'hypothèse en tête.

Avant (4 groupes) :
```
7) Conclusion : les hypothèses classées, ce sur quoi chacune repose, son soutien.
	Verdict : une différence existe très probablement. Cela ne dit ni qu'elle est importante, ni pourquoi ; les paires sont plus bas.
	Probabilité de chaque hypothèse, compte tenu des données et du point de départ :
```
Après :
```

7) Bilan - y diffère selon g : 99,8 % (soutien très fort).
	==> Verdict : une différence existe très probablement.
	==> Regroupement le plus probable : a = c < b = d (58,6 %).
	Note : cela ne dit ni que la différence est importante, ni pourquoi.
```

### 2. Ne pas dépasser 100 caractères par ligne
Il faut découper une idée par ligne, et replier le texte trop long avec `strwrap(width = getOption("width") - 8)` et un retrait suspendu. Aujourd'hui, 16 à 30 lignes par sortie dépassent 100 caractères, contre 2 ou 3 pour m.test.

Avant (étape 5, deux groupes, 243 caractères) :
```
	==> Sensibilité : taille attendue r de 0,25 à 1, chances de départ d'une différence de 1 contre 2 à 2 contre 1, répartition des variances entre groupes. « Aucune différence de y selon g » va alors de 48,0 % à 90,9 %. [m.bay() {KefiR}]
```
Après :
```
	==> Sensibilité : « Aucune différence de y selon g » va de 48,0 % à 90,9 %
	    quand r va de 0,25 à 1 et les chances de départ de 1 contre 2 à 2 contre 1.
```

### 3. Ne pas imprimer deux fois les paires (étapes 6 et 7)
Avec 4 groupes, les 6 paires apparaissent deux fois (12 lignes de 130 à 170 caractères). L'étape 6 garde le détail. L'étape 7 ne garde que le regroupement et la liste courte des paires tranchées.

Avant (étape 7) :
```
	  a / b : différents, soutien très fort (différents à 99,4 %) ; écart estimé b - a = 2 [0,64 ; 3,3], delta = 1,1 [0,29 ; 1,8].
	  a / c : semblables, soutien modéré (différents à 19,5 %) ; écart estimé c - a = 0,03 [-1,2 ; 1,3], delta = 0,016 [-0,63 ; 0,68].
	  ... (4 autres lignes identiques à l'étape 6)
```
Après :
```
	==> Différents : a/b, a/d, b/c, c/d (soutien très fort).
	==> Semblables : a/c (modéré), b/d (faible). Écarts chiffrés : étape 6.
```
Le même doublon existe pour deux groupes et pour les données appariées : `Écart estimé B - A = 0,42 [-0,91 ; 1,8]` à l'étape 6, puis `B - A = 0,42 [-0,91 ; 1,8], delta = 0,18 [...]` à l'étape 7.

### 4. Ne pas répéter à chaque appel les paragraphes de mode d'emploi
Quatre blocs reviennent à l'identique dans les six sorties, pour environ 700 caractères : la définition de l'écart estimé et de l'intervalle de crédibilité (≈ 250 caractères), « Ce sont des probabilités d'hypothèses, pas des p-values ... (tirage au sort, ordre, opérateur) », « Ces valeurs ne sont pas ajustées sur les données ; l'étape 5 les fait varier ... » et « Aucune forme n'est choisie : la réponse de chacune est pondérée par sa probabilité. » Il suffit de garder une ligne « Note : » courte et de renvoyer à `detail = TRUE`, qui existe déjà.

Avant :
```
	Écart estimé : médiane a posteriori, légèrement ramenée vers zéro par l'a priori, donc un peu plus petite que l'écart observé de l'étape 7. Crochets : intervalle de crédibilité à 95 %, si une différence existe (la valeur s'y trouve avec 95 % de probabilité ; ce n'est pas un intervalle de confiance).
```
Après :
```
	Note : [ ; ] = intervalle de crédibilité à 95 % (95 % de chances d'y trouver l'écart).
```

### 5. Réserver chaque symbole à un seul rôle
`==>` sert à un constat chiffré, `-->` à une décision qui change la suite de l'analyse, et « Note : » à tout le reste. Les lignes sans marque au milieu d'une étape disparaissent.

Avant :
```
	--> Toutes les paires viennent d'une seule loi sur les regroupements : aucune correction paire par paire n'est nécessaire, et les réponses ne peuvent pas se contredire.
	Les valeurs éloignées sont cherchées une fois retiré le niveau propre de chaque sujet.
```
Après :
```
	Note : pas de correction paire par paire ; les paires ne peuvent pas se contredire.
	Note : valeurs éloignées cherchées après retrait du niveau propre de chaque sujet.
```

### 6. Mettre le résultat dans le titre de l'étape, comme m.test
Le titre doit aussi correspondre au contenu. « Facteurs de Bayes des 5 modèles » n'affiche aucun des 5 facteurs de Bayes.

Avant (deux facteurs) :
```
5) Facteurs de Bayes des 5 modèles. [lmBF() {BayesFactor}]
	==> Précision : recalculées avec une autre graine, les probabilités bougent d'au plus 0,28 point(s) de pourcentage, négligeable.
```
Après :
```
5) Solidité de la réponse [lmBF() {BayesFactor}] - B change de côté selon l'a priori.
	==> Précision du calcul : ± 0,3 point, négligeable.
```
De même : `4) Formes du bruit [m.bay() {KefiR}] - normale, variances égales en tête (51,9 %).`

### 7. Sortir le jargon de méthode du texte courant
« Zellner-Siow », « échantillonnage d'importance », « HOP, Taylor-Rodriguez et al. 2016 », « Beta(1,1) », « rscaleCont, sqrt(2)/4 », « symétrie composée », « égal à lmBF(y ~ temps + id) à 0,2 près en échelle log » et « cotes (rapport oui/non) » vont dans `detail = TRUE`. Le récit garde le nom de la fonction entre crochets sur la ligne-titre.

Avant (apparié) :
```
	--> Différence appariée : test t apparié bayésien de Rouder et al. (2009), r = 0,5 comme pour un effet sujet aléatoire (égal à lmBF(y ~ temps + id) à 0,2 près en échelle log).
```
Après :
```
4) Point de départ (a priori) [ttestBF(paired = TRUE) {BayesFactor}]
	==> 50 % de chances qu'une différence existe ; taille attendue r = 0,5.
```
Avant (ANCOVA) : `pente a priori de Cauchy d'échelle sqrt(2)/4 (rscaleCont, le « medium » de BayesFactor), partagée avec ses interactions avec un facteur`. Après : `==> La covariable x est présente dans tous les modèles.`

### 8. Formater les nombres avec une seule règle
Écarts observés :
- `b - a = 2 [0,64 ; 3,3]` à côté de `d - a = 2,5` et `c - a = 0,03` : 2 chiffres significatifs, donc des décimales qui changent d'une ligne à l'autre ;
- `delta = 0,016` contre `delta = 1` ;
- `50 %` (plans à sujets ou covariable) contre `50,0 %` (un facteur) ;
- moyennes `20,08` (apparié) contre `10,4` (ailleurs) ;
- même facteur de Bayes affiché `599000` à l'étape 5 et `10^6` à l'étape 7 (ANCOVA) ;
- précision `0,00045 point(s) de pourcentage` et `9,1e-08 point(s) de pourcentage`, avec une notation scientifique au milieu d'une phrase en français.

Proposition :
- les moyennes, les écarts et leurs bornes avec le même nombre de décimales (celui des moyennes, ici 1) ;
- les probabilités toujours à 1 décimale ;
- le facteur de Bayes écrit pareil partout (`599 000` ou `≈ 6 × 10^5`) ;
- une précision inférieure à 0,1 point écrite `< 0,1 point`.

Avant : `écart estimé c - a = 0,03 [-1,2 ; 1,3], soit 0,016 [-0,63 ; 0,68] écart-type`
Après : `écart estimé c - a = 0,0 [-1,2 ; 1,3], soit 0,02 [-0,63 ; 0,68] écart-type`

### 9. Accorder les pluriels au lieu d'écrire « (s) »
Le récit contient `1 valeur(s) au-delà de 3 écarts interquartiles`, `0,28 point(s)`, `Covariable(s) centrée(s) sur leur moyenne (x = 4,71)`, `Terme(s) qui changent de côté : B` et `covariable(s) x`. R fournit `ngettext()` pour accorder.

Avant : `--> 1 valeur(s) au-delà de 3 écarts interquartiles de leur case (étape 2) : ...`
Après : `--> 1 valeur au-delà de 3 écarts interquartiles (étape 2) : lire la réponse avec réserve.`

### 10. Supprimer les doubles négations et les gabarits qui ne vont pas au plan
- `a / c : semblables, soutien modéré (différents à 19,5 %)` : le lecteur lit « semblables » puis « différents ».
- Données appariées à 2 conditions : `Cela ne dit ni qu'elle est importante, ni quels groupes diffèrent, ni pourquoi.` Avec deux conditions, la question ne se pose pas.
- Deux groupes : `Moyennes observées (... ; ce bilan ne dit pas quels groupes diffèrent)`, alors que le bilan donne `B - A`.
- Un seul terme : `--> Aucun terme ne change de côté.`
- Deux facteurs : la question de l'interaction est formulée autrement que les autres (`Les écarts selon A varient-ils avec B (A:B) ?`), ce qui se défend, mais les autres questions commencent par `Différence de y selon ...`.

Avant : `a / c : semblables, soutien modéré (différents à 19,5 %)`
Après : `a / c : semblables, soutien modéré (probabilité qu'ils diffèrent : 19,5 %)`

### 11. Garder le même squelette d'étapes pour tous les plans
L'ordre attendu est 1 Plan, 2 Contrôles, 3 Forme du bruit, 4 Point de départ, 5 Solidité (facteur de Bayes, précision, sensibilité), 6 Paires, 7 Bilan. Aujourd'hui, avec un facteur, l'a priori est à l'étape 3 et le bruit à l'étape 4 ; dans les autres plans, c'est l'inverse. Le plan à deux facteurs passe de 5 à 6) Conclusion, sans annoncer que les paires ne sont pas calculées. On peut afficher l'étape sans rien calculer de plus.

Avant (deux facteurs) : `6) Conclusion : les hypothèses classées, ...`
Après :
```
6) Quels groupes diffèrent - non disponible à plusieurs facteurs.
7) Bilan - y diffère selon au moins un facteur : 99,0 % (soutien fort).
```

### 12. Présenter les hypothèses de la même façon dans tous les plans
Avec un facteur, la conclusion donne une liste numérotée et alignée (`1. Aucune différence de y selon g : 73,9 %`), qui se lit bien. Dans les autres plans, les deux hypothèses tiennent sur une seule ligne, séparées par un point-virgule et suivies d'une majuscule.

Avant (ANCOVA) :
```
	  y diffère selon g, à même x : > 99,9 % (soutien très fort) ; Aucune différence de y selon g : < 0,1 %.
```
Après :
```
	  1. y diffère selon g, à même x  : > 99,9 % (soutien très fort)
	  2. Aucune différence de y selon g :  < 0,1 %
```
Une remarque liée, à l'affichage seulement : en ANCOVA, m.test indique que la covariable x a un effet, alors que m.bay ne dit rien de x. L'étudiant se demandera « et x ? ». Une ligne `Note : x est présente dans tous les modèles ; son effet n'est pas la question posée.` suffit.

## Points où m.bay est déjà meilleur que m.test

- **Typographie française** : virgule décimale et espace avant `%` (`99,8 %`). m.test mélange le point décimal, `p-value`, `1 outliers`, `of {agricolae}` / `de {KefiR}`, et imprime des tableaux R en anglais (`Analysis of Variance Table`, `Signif. codes`) et des nombres bruts (`7.701097`, `3.242368e-06` en ANCOVA).
- **Un verdict explicite et des réponses chiffrées** : probabilité de chaque hypothèse, taille de l'écart avec son intervalle, regroupement le plus probable. m.test s'arrête à `p-value : 0.46` et à des lettres.
- **Pas de contradictions entre méthodes** : avec 4 groupes, m.test donne des lettres inversées d'une méthode à l'autre (`a`/`b` en bootstrap, `b`/`a` en Scheffé), puis un tableau bootstrap avec `ab`. m.bay donne une seule réponse cohérente.
- **Chaque étape se ferme sur une décision** (`--> Aucune valeur n'est retirée.`). Les étapes 4 et 5 de m.test pour deux groupes n'ont ni `==>` ni `-->`.
- **Une longueur stable** : 37 à 51 lignes logiques quel que soit le plan, contre 22 à 126 pour m.test, qui gonfle en ANCOVA. m.test répète aussi l'indépendance aux étapes 5 et 6 en mesures répétées, et numérote une étape `9)   ==> Test normalité strict RÉUSSI`.
- **L'honnêteté sur les limites** : sensibilité aux a priori, valeurs éloignées, causalité. m.bay dit quand la réponse peut s'inverser (`l'ordre peut s'inverser`, `B sera dit non tranché`), ce que m.test ne fait jamais.
- **Des hypothèses alignées en colonnes** dans le cas à un facteur (`y diffère selon g              : 26,1 %`), un modèle à étendre aux autres plans (recommandation 12).
