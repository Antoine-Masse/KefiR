# Méthodologie de m.bay

Pendant bayésien de `m.test`. Pour chaque terme : « Ce facteur a-t-il un effet ? Oui avec probabilité p ; non avec probabilité 1-p », a priori annoncés, sensibilité affichée.

Les clés entre crochets renvoient à BIBLIOGRAPHIE.md, qui indique si le contenu a été lu ou seulement les métadonnées vérifiées.

**Contrôles numériques faits pour ce document** (R 4.6.1, BayesFactor 0.9.12-4.8, 2026-11-05, un jeu de données simulé par cas). Ils portent uniquement sur l'algèbre et non sur les performances : voir §3 et §7.

> **État au 5 novembre 2026 : proposition de départ, contestée.** La relecture adverse (`RELECTURE_ADVERSE_2026-11-05.md`) y trouve 6 défauts bloquants, ouverts, listés dans `DECISIONS.md`. Les positions du §1 ne sont pas validées, en particulier le « π = 1/2 apparié » (§3), le ROPE à plusieurs ddl (§4), la branche gaussienne sans aiguillage (§7), la multiplicité entre termes et le réglage MANOVA (§6). Les parties §3 et §4 seront réécrites après arbitrage d'Antoine.

## 1. Positions en bref
1. **P(terme)** = probabilité issue du BF d'inclusion *apparié* (analogue du Type II), avec π = 1/2 **par terme**, annoncée comme telle.
2. **H0 ponctuelle** en première ligne ; **ROPE** toujours en seconde ligne, sous un autre nom (« effet négligeable »). Ne jamais les fusionner.
3. **MANOVA** : BF analytique Normal-Inverse-Wishart à a priori g par terme et Jeffreys sur Σ, intégré sur g. Réduit exactement à `anovaBF` à d = 1 (contrôlé).
4. **Pas de test préliminaire qui choisit la branche.** Le gaussien est toujours calculé ; le robuste s'y ajoute, et un désaccord s'affiche au lieu d'être arbitré en silence.
5. **Sensibilité** toujours imprimée : grille r × cotes a priori 3 × 3, plus la plage des p.
6. **Post-hoc** : partitions de moyennes pour k ≤ 5, sinon paires à cotes a priori ajustées.
7. **Promesse de calibration** limitée à ce qui est démontrable : calibration sous l'a priori. Le taux de faux positifs fréquentiste se mesure par simulation et se publie.

## 2. Notations et moteur commun

Modèle gaussien :

    y = μ + Σ_T X_T β_T + ε,   ε ~ N(0, σ²I),   p(μ, σ²) ∝ 1/σ²

Chaque terme T a son propre g :

    β_T | σ², g_T ~ N(0, σ² g_T I_{q_T}),   g_T ~ InvGamma(1/2, r_T²/2)

avec X_T en codage réduit orthonormal (q_T = a-1 pour a niveaux), comme [Rouder2012]. Marginalement, β_T/σ suit une Cauchy multivariée d'échelle r_T.

**Échelles par défaut** (celles de [BayesFactorManual], pour rester comparable à JASP et à la littérature) :

| Terme | r | Lecture |
|---|---|---|
| Facteur fixe | 1/2 | « medium » de `anovaBF` |
| Deux groupes, sur δ | √2/2 ≈ 0,707 | Même a priori que r = 1/2 en réduit (δ = √2 β*). Contrôlé : `ttestBF` et `anovaBF` donnent le même log BF (−0,6246) |
| Covariable continue | √2/4 ≈ 0,354 | « medium » de `lmBF` |
| Effet aléatoire | 1 | « nuisance » |

Piège d'affichage : « medium » n'a pas la même valeur selon la fonction. `m.bay` imprime toujours la valeur numérique et l'objet qu'elle met à l'échelle.

**Passage du BF à la probabilité.** Avec des cotes a priori O = π/(1-π) :

    P(T | D) = O·BF / (1 + O·BF)

Choix de calcul : voies analytiques ou `BayesFactor` en routine (rapide), Stan + bridge sampling [Carpenter2017, Burkner2017, Gronau2017, Gronau2020] seulement quand le modèle l'exige. Chaque p porte son erreur numérique. Essais `lmBF` : erreur relative de **23 %** sur un rapport à pentes aléatoires sans réplique (une observation par sujet × condition, cas dégénéré, §5.7) ; 0,9 % avec 3 répliques ; 0,9 % puis 0,4 % sur A + B avec 10⁴ puis 5·10⁵ itérations. Il faut donc viser une erreur cible, l'afficher, et refuser le cas sans réplique.

## 3. Définir P(terme a un effet | données)

### Trois définitions
Soit M+ et M− deux ensembles de modèles, avec et sans T.

| Définition | Ensembles | Sens | Défaut |
|---|---|---|---|
| **Toutes** | Tous les modèles respectant la marginalité | Cohérente, car marginale d'une seule loi sur les modèles. Mais l'effet de A y est mélangé avec des modèles où A intervient dans une interaction | Non, affichée sur demande |
| **Appariée** (« matched », Baws) | M+ : modèles avec T sans interaction d'ordre supérieur impliquant T. M− : les mêmes sans T [Mathot2017, bayestestRInc] | Interprétable comme « ajouter T apporte-t-il quelque chose, ses termes d'ordre inférieur étant présents ? » | **Oui** |
| **Modèle plein moins T** (`whichModels = "top"`) | Plein vs plein sans T | Analogue du Type III. Pour un effet principal avec interaction présente, viole la marginalité | Option avertie |

Formule de l'appariée :

    BF_T = [ Σ_{m∈M+} π(m) p(D|m) / Σ_{m∈M+} π(m) ] / [ Σ_{m∈M−} π(m) p(D|m) / Σ_{m∈M−} π(m) ]

### Pourquoi l'appariée
Argument : A + B contre B est exactement la comparaison du Type II pour A. Elle respecte la marginalité [Nelder1977, Peixoto1987, Peixoto1990]. [Langsrud2003] fait un argument parallèle côté fréquentiste, et [Rouder2016] plaide pour la comparaison de modèles nommés.

Limite franche : **les p appariés ne sont pas les marginales d'une même loi jointe**. P(A) et P(A:B) viennent de deux sous-espaces distincts et ne se combinent pas. C'est le prix à payer pour une réponse terme à terme lisible. Si l'utilisateur veut la cohérence, l'option `inclusion = "all"` avec a priori de Chipman la donne (ci-dessous).

### Rapport à `$effects` de m.test

`$effects` donne une p-value par terme mais **pas de façon homogène** : Type III pour les modèles mixtes (`sys_effects_table.R`), Type II pour `car::Manova` et en repli pour certaines ANCOVA, un seul test global pour Kruskal et Friedman.

- `$effects` : P(statistique au moins aussi extrême | H0, modèle). `m.bay` : P(T | D, a priori, espace de modèles). Pas de conversion possible entre les deux [WassersteinLazar2016, BergerSellke1987].
- `m.bay` conserve les **mêmes noms de termes** pour permettre une lecture côte à côte, et une colonne `analogue` (II, III) : un désaccord entre p-value et probabilité peut venir de l'hypothèse testée, pas seulement du paradigme.
- Cas où l'écart est attendu : un p ≈ 0,04 fréquentiste peut donner BF ≈ 1. La borne −e·p·ln p [Sellke2001] le montre : à p = 0,05, BF01 ≥ 0,41, donc P(H0|D) ≥ 0,29 avec π = 1/2. Le bilan doit l'expliquer.

### Influence de l'a priori sur les modèles

Énumération exacte des modèles respectant la marginalité, a priori uniforme sur les modèles (mêmes comptes et inclusions que `anovaBF` à 3 facteurs, contrôlé) :

| Facteurs | Modèles | P(effet principal) | P(interaction d'ordre 2) | P(ordre 3) |
|---|---|---:|---:|---:|
| 2 | 5 | 0,600 | 0,200 | – |
| 3 | 19 | 0,737 | 0,316 | 0,053 |

Donc « uniforme sur les modèles » **n'est pas** un a priori neutre : il favorise les effets principaux et pénalise les interactions sans le dire. Toute annonce « P(H0) = 1/2 » serait fausse avec ce réglage. C'est le point le plus trompeur des réglages courants.

| A priori | Ce qu'il dit | Décision |
|---|---|---|
| Uniforme sur les modèles | Inclusions de 0,6 ; 0,2 ; etc. | Rejeté comme défaut |
| **π = 1/2 par terme, BF apparié** | Chaque terme part de 1:1 | **Défaut**, jusqu'à 3 facteurs |
| Hérédité forte, inspirée de [Chipman1996] ; le réglage π(A:B | A, B) = 1/2 est propre à KefiR, non sourcé | Cohérent. Avec 2 facteurs : P(A) = P(B) = 1/2, P(A:B) = 1/8 | Sensibilité (`inclusion = "all"`) |
| Beta-binomial sur le nombre de termes [ScottBerger2010] | Ajustement automatique à la multiplicité | Défaut au-delà de 3 facteurs |

Justification du seuil : un π fixe n'ajuste pas la multiplicité [ScottBerger2010]. Avec ≤ 7 termes et un plan confirmatoire, la lisibilité l'emporte ; au-delà, l'exploration domine. Seuil de jugement, non démontré. Espace trop grand : exploration à la [Clyde2011, BAS]. Le modèle à probabilité médiane [BarbieriBerger2004] s'applique aux inclusions cohérentes, donc à `inclusion = "all"`, pas au défaut apparié.

## 4. H0 ponctuelle ou ROPE

**Elles répondent à deux questions** [Makowski2019] :

- *Existence* : l'effet est-il exactement nul ? Pertinent en plan expérimental, où une nullité vraie (ou quasi) est plausible : randomisation, manipulation peut-être inopérante.
- *Importance* : l'effet est-il plus grand qu'un seuil négligeable ? Pertinent en observationnel, où aucun effet n'est strictement nul, et à grand n, où de petits écarts font rejeter la nullité ponctuelle [MoreyRouder2011].

**Défaut** : le bilan répond à la question posée avec H0 ponctuelle, et ajoute toujours la ligne ROPE. Option `null = "rope"` pour inverser l'ordre. Sans H0 ponctuelle, un « non » ne serait pas accessible avec des a priori continus.

**Calcul.** Le ROPE moyenné sur H0 et H1 :

    P(ROPE | D) = P(H0 | D) + P(H1 | D) · P(|θ| < w | D, H1)

Le terme de droite vient des tirages a posteriori (`posterior()` ou Stan), sur la **loi entière** (`ci = 1` de [bayestestRrope]), pas sur l'intervalle à 95 % : ce dernier ne donne pas une probabilité.

**Mesure et largeur.** Pour tout terme, à l'échelle standardisée :

    f_T = sqrt( mean_j(α_j²) ) / σ,   α_j : effets de cellule centrés, σ résiduel

| Quantité | ROPE par défaut | Remarque |
|---|---|---|
| 2 groupes, δ | \|δ\| < 0,1 | [Kruschke2018] via [bayestestRrope] (±0,1·sd) |
| Terme à plusieurs ddl | f < 0,05 | Identique à \|δ\| < 0,1 pour 2 groupes égaux, puisque f = \|δ\|/2 |
| Pente de covariable | \|β·sd(x)/σ\| < 0,1 | Effet d'un écart-type de x |

Différences à annoncer : `bayestestR` standardise par sd(y), m.bay par σ résiduel. Avec une covariable forte, les deux divergent. 0,1 est une convention, pas un résultat. **La largeur doit être remplaçable par une valeur métier**, et une valeur métier vaut mieux.

Masse a priori à afficher : avec r = 0,707, la Cauchy met 2/π·atan(0,1/0,707) ≈ 0,089 dans le ROPE, donc P(ROPE) a priori ≈ 0,545 avec π = 1/2. Sans cela, le lecteur croit que le ROPE part de 0.

**Présentation** (chiffres fictifs) :

    A : a-t-il un effet ?          Oui p = 0,91 ; non 0,09   [H0 : effet exactement nul ; BF = 10 ; π = 1/2 ; r = 0,5]
        effet non négligeable ?    p = 0,62 ; négligeable 0,38 [ROPE : f < 0,05 ; masse a priori 0,55]

La confusion à éviter : p élevée à la première ligne et ROPE élevé à la seconde peuvent coexister. Le bilan le dit : « effet probablement réel, mais peut-être trop petit pour compter ».

## 5. Plans

Pour chaque plan, `P(T|D)` suit le §3 et la ligne ROPE le §4. Seul ce qui change est indiqué.

### 5.1 Deux groupes

- **Modèle** : y_i = μ ± δσ/2 + ε. Cauchy(0, 0,707) sur δ [Rouder2009].
- **Probabilités** : P(δ ≠ 0 | D) ; P(|δ| ≥ 0,1 | D) ; signe, P(δ > 0 | D, H1).
- **Limite** : à petit n, le résultat dépend de r. Le défaut ne remplace pas un a priori informé si l'utilisateur en a un [Gronau2020b].

### 5.2 Trois groupes ou plus

- **Modèle** : §2, un g pour le facteur. Forme exacte contrôlée : la formule bloc reproduit `anovaBF` déséquilibré (log BF 0,69632 contre 0,69632).
- **Probabilité** : P(au moins deux moyennes diffèrent | D). Terme global, ne dit pas lesquelles (§9).
- **Limite** : l'a priori lie tous les niveaux par un même g, et son contenu marginal sur une paire change avec le nombre de niveaux.

### 5.3 Plusieurs facteurs et interactions

- **Modèle** : un g par terme, interactions comprises [Rouder2012, Rouder2017]. `whichModels = "withmain"`.
- **Probabilités** : §3. Pour un effet principal : lecture « A ajoute, B présent, sans interaction A:B ». Si P(A:B) > 1/2, le bilan signale que la lecture de A seul est conditionnelle à un modèle peu crédible et renvoie aux effets simples.
- **Contrôle** : additif A + B déséquilibré, formule à deux g contre `lmBF` : log BF 1,084 contre 1,078 (erreur 0,4 % à 5·10⁵ itérations). Concordant à l'erreur Monte Carlo près.
- **Limite** : l'espace explose avec le nombre de facteurs (5, puis 19 modèles à 2 et 3 facteurs). Plan très déséquilibré ou cellules vides : modèles non identifiables, refuser de répondre sur les termes concernés.

### 5.4 ANCOVA

- **Modèle** : y = μ + X_A β + x γ + ε. Covariable centrée, r = 0,354 [RouderMorey2012].
- **Probabilité** : BF(A + x vs x). La covariable est présente dans les deux modèles, avec un a priori identique, comme nuisance. P(A | D, x), pas P(A).
- **Limites** : le BF sur A dépend de la covariable retenue. Une covariable choisie après coup invalide la lecture. Linéarité supposée. Covariable mesurée après traitement : aucun cadre statistique ne le corrige.

### 5.5 Pentes hétérogènes

- **Modèle** : ajout de A:x. **Probabilité** : BF(A + x + A:x vs A + x).
- **Si A:x est crédible** : « effet de A » n'existe plus comme quantité unique. Le bilan donne l'effet de A à des valeurs de x (moyenne, ±1 sd) par tirages a posteriori, non un P(A) unique.
- **Limite** : `lmBF` accepte A:x sans erreur (essai fait), mais **l'échelle a priori de ce terme reste à documenter** avant d'annoncer quoi que ce soit, le manuel consulté ne la précisant pas. Sinon brms avec a priori explicites.
- **Lien m.test** : `sys_ancova_analysis.R` teste l'interaction facteur × covariable de façon hiérarchique. m.bay la garde dans l'espace de modèles.

### 5.6 Apparié

- **Modèle** : sur les différences d_i : d = δσ_d + ε, Cauchy(0, 0,707) [Rouder2009].
- **Limite franche** : δ est standardisé par σ_d, pas par l'écart-type intra-groupe. Même r = cible différente selon la corrélation. Le ROPE doit préciser cette échelle, sinon un δ = 0,1 apparié et un δ = 0,1 indépendant ne sont pas comparables.

### 5.7 Mesures répétées et split-plot

- **Modèle** : sujet aléatoire présent partout (r = 1), facteurs fixes r = 1/2 [Rouder2017]. Split-plot : facteur inter, facteur intra, et leur interaction.
- **Défaut** : pentes aléatoires quand le plan a plusieurs observations par sujet × cellule. Avec **une seule** observation par sujet × condition, A:ID se confond avec le résidu : pas de pente estimable.
- **Probabilité** : §5.8 à plusieurs observations ; sinon BF(A + ID vs ID) apparié.
- **Limite** : pas de correction de sphéricité : l'a priori est composé-symétrique. Si la sphéricité est fausse, les BF héritent de cette hypothèse. `m.test` traite la question par ailleurs, le bilan doit le signaler. Voie alternative : MANOVA sur les contrastes intra (§6).

### 5.8 Mixte avec pentes aléatoires

- **Modèle** : y = Xβ + Z_id u + Z_{A:id} v + ε, composantes indépendantes dans BayesFactor (`whichRandom = c("id", "A:id")` : fonctionne à l'essai, mais le manuel ne décrit pas cet usage). Corrélations pentes-intercepts : seulement via brms, LKJ [Lewandowski2009].
- **Quelle comparaison ?** Le résultat dépend de la structure aléatoire sous H0 :
  - **Défaut : pentes présentes sous H0 et H1**, BF(A + id + A:id vs id + A:id). Cohérent avec [Barr2013], et [vanDoorn2021] rapporte que l'ignorance des pentes peut gonfler le BF.
  - **Sensibilité : moyenne sur structures** [Heck2021], qui soutient qu'aucune comparaison ne convient toujours.
- **Probabilité annoncée** : P(effet moyen de A ≠ 0 | D, pentes présentes). Ne pas présenter cela comme « A n'a d'effet sur personne » : H0 ici autorise des effets individuels de moyenne nulle.
- **Limites** : plusieurs effets aléatoires croisés ou peu de niveaux : erreur MC à surveiller (§2). Peu de sujets : Stan avec demi-t sur les écarts-types [Gelman2006]. Le débat n'est pas clos [vanDoorn2021, Heck2021].

MANOVA : §6. Non-normalité, hétéroscédasticité, petits effectifs : §7.

### 5.9 Données discrètes

| Donnée | Modèle | A priori | Calcul |
|---|---|---|---|
| Binaire, binomiale | Logistique | Cauchy(0, 2,5) sur prédicteurs mis à l'écart-type 0,5 [Gelman2008] (soit 1,25 sur prédicteurs réduits à 1) | brms + bridge sampling |
| Comptages | Poisson ou négatif binomial si surdispersion | Normale faible, échelle log, à calibrer par prédictif a priori [Schad2021] | Idem |
| Ordinale, peu de modalités | Cotes proportionnelles [McCullagh1980, BurknerVuorre2019] | Effets : même logique ; seuils Dirichlet | Idem |
| Tableau de contingence | Multinomial | Défaut de `contingencyTableBF` | `BayesFactor` |

**Probabilité** : celle du §3 avec un autre terme d'erreur. **Lien avec m.test** : quand le chemin discret de `m.test` traite les comptages comme continus ou en rangs, `m.bay` ne doit pas appliquer le gaussien sans le dire [LiddellKruschke2018].

**Limites** : la Cauchy sur un log-odds n'est pas le même objet que la JZS ; il n'y a pas d'équivalence entre r = 1/2 et 2,5. Séparation complète à petit n : l'a priori décide, il faut l'imprimer. ROPE sur log-odds : ±0,18 est un usage à vérifier, je ne l'ai pas sourcé et je ne le propose donc pas comme défaut. L'utilisateur fournit la valeur.

## 6. MANOVA

`BayesFactor` n'a rien pour cela (vérifié dans [BayesFactorManual]). Options réelles :

| Option | Pour | Contre | Décision |
|---|---|---|---|
| **A. NIW analytique** à g-prior, Jeffreys sur Σ | Exact, rapide, invariant, réduit au défaut univarié | Normalité multivariée, Σ commun, à écrire nous-mêmes | **Défaut** |
| B. BIC multivarié [Wagenmakers2007, Schwarz1978, KassWasserman1995] | Une ligne depuis le Λ de Wilks | A priori implicite peu lisible, biais vers H0 à n modeste | Contrôle de cohérence seulement |
| C. Stan `mvbind` + bridge [Burkner2017, Gronau2017] | Hétéroscédastique, t, mixte | Lent, MC. Sa prise en charge des a priori sur effets standardisés demande un code Stan dédié | Si A ne suffit pas |
| D. BF univariés combinés | Simple | Le produit suppose des DV indépendantes, et compte double une évidence corrélée. **Faux si r ≠ 0** | **Rejeté comme test global** ; affichés DV par DV à titre descriptif |
| E. BFpack [Mulder2021, Mulder2019, MulderGu2022] | Existant ; contraintes d'ordre | Fractionnaire, autre a priori, pas de moyenne de modèles ni de cohérence avec r | Comparaison externe |

### Formule retenue

Y (n × d), lignes N_d(XB, Σ). Y et X centrés, μ plat. Cas général (a priori par terme) :

    B_T | Σ, g_T ~ MN(0, g_T I_{q_T}, Σ),   g_T ~ InvGamma(1/2, r²/2),   p(Σ) ∝ |Σ|^{-(d+1)/2}

Avec V0 = diag(g_T I_{q_T}), A = V0⁻¹ + X'X, S = Y'Y − Y'X A⁻¹ X'Y :

    log m(Y | g) = C + (d/2)(−log|A| − log|V0|) − ((n−1)/2) log|S|

où **C est commun à tous les modèles**, car tous partagent Σ, μ et ν = n − 1. Le modèle nul : log m0 = C − ((n−1)/2) log|Y'Y|. Puis :

    BF_T = ∫ exp(log m(Y | g) − log m0) ∏ InvGamma(g_T) dg

Cas particulier d'un seul terme, en plan équilibré (n_c par cellule, X̃'X̃ = n_c I) :

    BF(g) = (1 + n_c g)^{-q d / 2} · ∏_{i=1}^{s} (1 − (n_c g / (1 + n_c g)) θ_i)^{-(n−1)/2}

où θ_i sont les valeurs propres de E0⁻¹H, avec E0 la SSCP totale corrigée et H la SSCP de l'hypothèse, c'est-à-dire **les mêmes objets que Pillai et Wilks**. Lien direct avec le `car::Manova` de `sys_manova_analysis.R`.

**Pourquoi ce choix** : pour toute combinaison a'Y, a'B | Σ, g a pour variance g·a'Σa : **chaque combinaison linéaire reçoit l'a priori univarié standardisé** avec le même r (partageant g). Cela donne une interprétation et une cohérence uniques. La structure conjuguée suit [Press2005, Rowe2002, Brown1998].

### Contrôles (un jeu de données par cas)

| Contrôle | Résultat |
|---|---|
| Formule NIW propre générale contre Monte Carlo brut (n = 6, d = 2, 200 000 tirages) | −19,42 contre −19,37 |
| Forme fermée contre limite NIW numérique | −1,076013 des deux côtés |
| d = 1 contre forme de [Liang2008] | −0,07531919 des deux côtés |
| d = 1 intégré sur g contre `anovaBF`, r = 0,25 ; 0,5 ; 1 | −0,05792 / −0,33385 / −0,94861 identiques à 10⁻⁶ |
| Invariance à Y → YA | Identique |
| BIC contre analytique, H0 vraie, n = 30, d = 3, r = 1/2 | log BF −8,46 contre −2,41 |

Dernière ligne : **une** simulation, mais le BIC y favorise H0 bien plus que l'analytique. Raison de plus pour ne pas en faire la réponse principale. À quantifier.

### Piège des a priori impropres
1. **μ plat et Jeffreys sur Σ sont permis, car communs à tous les modèles** : leur constante s'élimine. Conditions : même d, même ν = n − 1 partout, et même centrage [ConsonniVeronese2008].
2. **Un a priori impropre sur B_T est interdit** : constante arbitraire, BF arbitraire. Plus l'a priori propre est vague, plus H0 est favorisée [Lindley1957, Bartlett1957].
3. **Une Inverse-Wishart propre « vague » n'est pas anodine** : (ν0, S0) influence le BF, S0 impose une échelle et brise l'invariance à YA si S0 ∝ I. Si l'on préfère un a priori propre, S0 doit être calé empiriquement (avec le soupçon d'un double usage) ou informé, et on affiche l'effet de S0. D'où le choix de Jeffreys. Alternatives : [HuangWand2013], LKJ [Lewandowski2009], non conjuguées, donc voie Stan.
4. **Refus explicite** si n − 1 < d (E0 singulière). Si n − 1 − Σ q_T < d, E1 est singulière, des θ_i valent 1 et l'intégrale sur g peut diverger : refuser aussi, au lieu d'imprimer un nombre.

### Probabilités annoncées
- **Globale** : P(T a un effet sur le vecteur des DV | D).
- **Par DV (descriptive)** : BF univarié pour chaque DV, sans multiplicité ajustée, explicitement pas une décomposition de la globale.
- **ROPE multivarié** : Mahalanobis f²_T = tr(Σ⁻¹ B'B)/q ; seuil non standard. **Aucun défaut proposé**, à fournir ou à omettre. Une convention inventée serait pire.

### Limites

Σ commune par hypothèse, sans équivalent bayésien automatique d'un Box's M ; normalité multivariée ; le BF baisse mécaniquement avec d (pénalité q·d) quand les DV ajoutées sont sans rapport : observé −0,33 à d = 1 et −0,87 à d = 3, r = 0,5, un jeu de données. **r = 1/2 sur d grand n'a pas été calibré.** Facteurs multiples : intégration numérique sur g de dimension = nombre de termes, faisable jusqu'à 3 environ (quadrature), au-delà QMC ou Stan, à évaluer. Contrairement à d = 1, aucune référence multivariée de type `anovaBF` ne permet de valider le défaut : seuls les contrôles ci-dessus le soutiennent. **Cette construction est ma proposition, pas une méthode publiée telle quelle** : sa structure est classique, le réglage (g par terme à r partagé avec l'univarié) n'a pas de source vérifiée, il doit donc être présenté comme propre à KefiR et validé par simulation.

## 7. Branche non normale

### Faut-il un test préliminaire qui choisit la branche ?

**Non**, malgré la cohérence apparente avec `m.test`.

- Choisir par test préliminaire puis annoncer P(T | D, branche choisie) **ignore l'incertitude de sélection**. Ce que [Hoeting1999] reproche à toute sélection. Le p annoncé serait conditionnel à un modèle tenu pour certain, alors qu'il a été choisi sur les mêmes données.
- À petit n, un test de normalité n'a pas de puissance ; à grand n, il détecte l'insignifiant, alors que la conséquence réelle sur le BF peut être petite.
- [Oelrich2020] : les probabilités de modèle sont surconfiantes quand la variabilité est sous-estimée, ce qui est le cas des queues épaisses sous gaussien. Ce risque est réel.

**Mais la moyenne de vraisemblances n'est pas gratuite non plus** : lente, dépend des a priori de chaque branche, et mélange des paramètres d'échelles différentes (δ gaussien n'est pas δ Student). Le compromis proposé :

| Niveau | Déclenchement | Ce que l'on fait |
|---|---|---|
| 1 | Toujours | Gaussien, avec diagnostics (m.bay réutilise ceux de m.test, **à titre d'alerte, pas d'aiguillage**) |
| 2 | Alerte, ou `robust = TRUE` | Ajoute la ligne robuste : les deux p s'affichent, plus « concordant » ou « discordant » |
| 3 | `robust = "bma"` | P(T | D) = Σ_ℓ P(ℓ | D) P(T | D, ℓ), avec bridge sampling. Poids de vraisemblance affichés |

En cas de discordance, le bilan **ne désigne pas de gagnant** : il dit que la réponse dépend de la forme du bruit et que p n'est pas robuste. Ce sera moins confortable que `m.test`, c'est la réalité.

### Les modèles

**Student-t** [Lange1989, Kruschke2013, Burkner2017]. y ~ t_ν(Xβ, σ), même a priori standardisé sur β/σ, ν > 1 faiblement informé (voir [Kruschke2013] pour la forme, à recopier d'après la source avant implémentation). Pour rester à l'échelle JZS, un **code Stan dédié** est nécessaire : `brms` place par défaut ses a priori sur β et non sur β/σ. Calcul : bridge sampling, avec prudence sur l'erreur [Schad2023].

**Hétéroscédastique** : log σ_j par groupe, a priori hiérarchique faible. δ rapporté à σ_p = sqrt(mean σ_j²) : un choix, d'autres références de σ donneraient d'autres nombres.

**Rangs** [vanDoorn2020] : pour 2 groupes et apparié, disponible, à utiliser. **Pour k > 2 groupes, je n'ai trouvé aucune méthode publiée vérifiée d'équivalent bayésien de Kruskal-Wallis** (recherche Crossref sans résultat, qui n'exclut pas l'existence d'une méthode).

### Si m.test passerait en Kruskal

Proposition : **cotes proportionnelles** [McCullagh1980, BurknerVuorre2019], avec y traité comme ordonné.

- Pourquoi : n'utilise que l'ordre, donc invariant aux transformations monotones comme un test de rangs. Passe à plusieurs facteurs et mixte.
- **Limites** : ce n'est pas « le Kruskal bayésien », et le bilan ne doit pas employer ce terme. Un seuil par valeur distincte à données continues, donc coûteux, avec a priori sur les seuils à vérifier. Hypothèse PO elle-même testable. Pour Friedman : effet aléatoire sujet. Le lien entre Wilcoxon et test de score PO est de mémoire, à sourcer avant de l'écrire dans une documentation.

### Petits effectifs

**Dire la vérité** : à n < ~10 par groupe, le BF dépend beaucoup de r et de la forme de vraisemblance, et t, hétéroscédastique ou rangs n'ont aucune information pour se distinguer. Le bilan doit afficher la plage de sensibilité ; si elle traverse 1/2, il répond « les données ne tranchent pas », même si le point central est p = 0,8. Le seuil « ~10 » est un repère à valider par simulation.

## 8. A priori sur les hypothèses et sensibilité

**Défaut P(H0) = P(H1) = 1/2.** Argument : ce n'est pas neutre mais c'est lisible, et le BF se lit alors directement. La question « Oui avec probabilité p » **oblige** à le poser : il n'y a pas de p sans π [Kass1995, Morey2016b]. Présenter p sans π serait la faute principale de cet outil. L'utilisateur peut passer `prior_odds`.

Imprimé systématiquement (ici, pour BF = 10 au r central ; ce sont des calculs, pas des résultats) :

    Sensibilité de P(A a un effet | D)
                      cotes H1:H0
      r            1:4     1:1     4:1
      r/2 = 0,25   ...     ...     ...
      r   = 0,5    0,714   0,909   0,976
      2r  = 1,0    ...     ...     ...
    Plage : [min ; max]. Conclusion stable si tous les p sont du même côté de 1/2.

Les neuf cases coûtent peu : BF à 3 r, cotes par simple calcul. Pour Stan, trois exécutions : option `sensitivity = FALSE` si nécessaire, le bilan le mentionne alors.

Forme de la réponse si instable (exemple fictif) : « p = 0,91 (1:1, r = 0,5), mais 0,41 à 0,98 selon les a priori : conclusion dépendante des a priori ». Ne pas choisir un seuil de « preuve » : échelles de [Kass1995, Jeffreys1961] affichables en option, sans valeur normative. [Benjamin2018] relie 0,005 au BF, sans raison d'y adosser le défaut.

**Comportement à attendre** : doubler r baisse le BF quand l'effet est modeste et le monte pour de très gros effets [Lindley1957, Bartlett1957, Liang2008]. **La sensibilité n'est pas un remède au prior mal choisi**, elle le révèle. Note sur [TendeiroKiers2019] : le BF est sensible aux a priori et une « vérité » de H0 ou H1 n'est pas garantie : la grille ne le règle pas.

## 9. Post-hoc sans correction fréquentiste

**Pas de Bonferroni ni de Tukey : la multiplicité passe par l'a priori** [Westfall1997, ScottBerger2006, GopalanBerry1998].

### Option A (défaut, k ≤ 5) : partitions de moyennes

Une hypothèse par partition des k moyennes (Bell : 5, 15, 52 pour k = 3, 4, 5). Pour chacune, BF contre le nul par `lmBF` avec le facteur recodé. Puis :

    P(μ_i = μ_j | D) = Σ_{partitions réunissant i et j} P(partition | D)

- Pour : **cohérent et transitif**, pas de résultat « a = b, b = c, a ≠ c » impossible en probabilité. Inspiration : [GopalanBerry1998, vandenBerghDablander2025]. Cette dernière référence, non lue à ce stade, à consulter avant de finaliser l'a priori.
- A priori proposé : P(toutes égales) = 1/2, le reste uniforme sur les autres partitions. Compatible avec le π du terme global. Si π(T) = 1/2 : P(μ_i = μ_j) a priori = 1/2 + (1/2)·(B_{k-1} − 1)/(B_k − 1), soit **0,625 à k = 3** (2 partitions sur 5 réunissent i et j, dont le nul) et 0,5 + 0,5·4/14 ≈ 0,643 à k = 4. À afficher, car ce n'est pas 1/2.
- Limite : a priori sur les partitions lui-même discutable (uniforme contre processus de Dirichlet). Égalités exactes entre moyennes moins plausibles en observationnel : proposer alors le ROPE deux à deux, et la probabilité de signe.

### Option B (k > 5) : paires à cotes ajustées

Avec m = k(k-1)/2 paires, et l'idée que l'a priori sur la nullité globale fixe P(H0 famille) = π0 (ici 1/2). Sous indépendance a priori : π_ij = π0^{1/m}.

| k | m | P(H0_ij) a priori | Cotes H1:H0 par paire |
|---|---:|---:|---:|
| 3 | 3 | 0,794 | 0,26 |
| 4 | 6 | 0,891 | 0,12 |
| 5 | 10 | 0,933 | 0,07 |

**Cette forme simple est ma dérivation**, dans l'esprit de [Westfall1997], dont je n'ai pas lu les théorèmes (indépendance fausse en pratique pour des paires partageant des groupes). À présenter comme réglage de m.bay. Effet : pénalité croissante, proche d'un Bonferroni.

**Dans le bilan** : pas de lettres de groupes comme si c'était une décision : tableau des paires avec P(égales), P(différence hors ROPE), signe. Autre voie, pour l'estimation et non le test : rétrécissement hiérarchique [Gelman2012], qui répond à une autre question.

**Limite** : des BF par paires sans ajustement ne sont pas protégés du seul fait d'être bayésiens. Sans a priori tenant compte de m, aucune protection [ScottBerger2010].

## 10. Calibration

### Ce que l'on peut garantir

**Sous le modèle et l'a priori effectivement utilisés pour générer les données**, et calcul exact : parmi les analyses annonçant p, une proportion p ont H1 vraie. C'est le théorème de Bayes moyenné sur le prédictif a priori.

Vérifications exécutables :

1. **Moyenne a posteriori = a priori** (DAP) : simuler M, θ, y sous l'a priori, calculer P(M|y). La moyenne doit égaler π [Schad2023]. Schad et al. le précisent : la condition concerne la moyenne sur simulations, **pas un jeu de données unique**.
2. **Courbe de calibration** : p annoncée contre fréquence de H1.
3. **SBC** pour les paramètres [Talts2018, Cook2006, Modrak2025], et R-hat, ESS [Vehtari2021].
4. **Erreur BF** : même analyse répétée, écart-type de p ; signal immédiat si l'erreur MC dépasse un seuil fixé (ex. 5 %).

Ce sont **des tests de calcul**, pas des tests de vérité du modèle. Ils doivent couvrir chaque voie du §5, MANOVA comprise.

### Ce que l'on ne peut pas garantir

- **Taux de faux positifs fréquentiste.** Pour θ = 0 fixe, P(p > 0,95) n'est pas 5 % : il dépend de n, de r, du plan. Attendu, faible à grand n, puisque le BF est consistant sous H0.
- **Calibration sous a priori faux** : si les vrais effets ne ressemblent pas à la Cauchy, p peut être mal calibrée. [Oelrich2020] est utile sur la dispersion du BF sous spécification erronée.
- **Sens à un seul jeu de données** : p = 0,91 ne signifie pas « 91 % des études de ce type » sans admettre que π et r décrivent le domaine.
- **Arrêt optionnel** : [Rouder2014] le défend pour le BF, le débat suppose un a priori exact. Ne pas promettre la même propriété sous sensibilité large.

Le débat sur la calibration fréquentiste du BF est ouvert [Hoijtink2016, Morey2016c] : m.bay ne recalibre pas le BF, il **publie** ses opérateurs. Position : la p annoncée reste bayésienne, et une table simulée la complète.

### Table à publier (vignette)

Plans × n × d vrai (0 ; 0,2 ; 0,5 ; 0,8) × distribution (gauss, t3, asymétrique, hétéroscédastique), méthode inspirée de [Schonbrodt2018] :

| Sortie | Pour d = 0 | Pour d > 0 |
|---|---|---|
| P(p > 0,95) | « Faux positifs » | Puissance analogue |
| P(p < 0,05) | Soutien à H0 | Fausses nullités |
| P(0,25 < p < 0,75) | Indécis | Indécis |

Comparer avec `m.test` aux mêmes scénarios.

**Phrase à écrire dans la documentation** : « m.bay ne contrôle pas un taux d'erreur à long terme. Sa probabilité est conditionnelle aux a priori affichés et au modèle. La table de simulation fournie indique ses opérateurs fréquentistes pour les plans testés. »

## 11. Risques et points à trancher
1. **« Oui avec probabilité p » sera lu comme une vérité.** Sans π et modèle, trompeur. Imposer le rappel des a priori sur la même ligne.
2. **Simulation** (§10) pas faite ici : aucune affirmation de performance dans ce document.
3. **MANOVA §6** : algèbre contrôlée, calibration pas. Réglage propre à KefiR.
4. **Lectures manquantes** avant citation : [vandenBergh2020], [Hinne2020], [Westfall1997], [vandenBerghDablander2025], et tout le niveau « métadonnées » de la bibliographie.
5. **Choix laissés à Antoine** : seuil de 3 facteurs pour Beta-binomial ; `robust` niveau 2 par défaut ou à la demande ; partitions à k ≤ 5 ; acceptation d'un p non cohérent entre termes (matched) contre une cohérence moins lisible.
