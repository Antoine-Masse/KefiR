# g.bay() : décisions

`g.bay()` est le pendant bayésien du test G (`DescTools::GTest`) : lien entre variables qualitatives à partir d'effectifs. Même famille que `m.bay()` (branche `feat/m.bay`) : a priori fixés d'avance, jamais ajustés sur les données, analyse de sensibilité, paquets établis, bilan sans p-value. Les décisions sont numérotées G1, G2… pour ne pas se confondre avec les D1 à D28 de `dev/mbay/DECISIONS.md`.

Toutes les références ci-dessous ont été vérifiées le 11 novembre 2026 (notice d'éditeur, PubMed/PMC, dépôt institutionnel ou manuel CRAN), sauf mention contraire.

## G1. Ce que fait `GTest`, et ce que g.bay en reprend : DÉCIDÉ

Lecture du code de `DescTools::GTest` 0.99.60 :
- `x` vecteur (ou matrice à une ligne ou une colonne, ramenée à un vecteur) : **ajustement** à `p` (défaut : proportions égales), `rescale.p` pour normaliser `p`. G = 2 Σ O log(O/E), ddl = k − 1.
- `x` matrice ou table, ou `x` et `y` deux facteurs (`table(x, y)` après retrait des manquants) : **indépendance**, ddl = (r − 1)(c − 1).
- `correct = "williams"` (Williams 1976) divise G par un facteur q > 1 pour rapprocher sa loi d'un χ² aux petits effectifs ; `"yates"` (2 × 2 ou 2 modalités) rapproche chaque effectif de son attendu de 0,5.

Ce que g.bay garde : les mêmes arguments `x`, `y`, `p`, `rescale.p` et les mêmes cas, plus une formule avec `data`. Ce qu'il ne garde pas : `correct`. Les deux corrections servent à rapprocher une statistique de sa loi asymptotique ; le facteur de Bayes de G2 et G3 est **exact** à tout effectif, sans approximation à corriger. `correct` est accepté, ignoré, et le récit le dit.

## G2. Ajustement à des proportions données : DÉCIDÉ

**Modèle.** Effectifs n = (n₁…n_k) multinomiaux de proportions π.
- H0 : π = p (point, les proportions théoriques données).
- H1 : π ~ Dirichlet(a, …, a), a = 1 (toutes les proportions également possibles d'avance).

**Facteur de Bayes, forme close** (le coefficient multinomial se simplifie) :
BF₁₀ = [Γ(ka) / Γ(ka + N)] Π Γ(a + n_i) / Γ(a) ÷ Π p_i^{n_i}.

**Références.**
- Good, I. J. (1967). A Bayesian significance test for multinomial distributions. *JRSS B*, 29(3), 399-431 (avec discussion). Le test « Dirichlet contre point » d'une multinomiale ; Good met en plus une loi sur a, que l'on ne reprend pas (a fixé d'avance, principe de la famille).
- Sarafoglou, A., Haaf, J. M., Ly, A., Gronau, Q. F., Wagenmakers, E.-J., & Marsman, M. (2023). Evaluating multinomial order restrictions with bridge sampling. *Psychological Methods*, 28(2), 322-338. https://doi.org/10.1037/met0000411. Même modèle, a = 1 par défaut ; paquet CRAN `multibridge`, fonction `mult_bf_equality`.

**Paquet ou forme close ?** Aucun paquet établi ne donne ce BF pour des proportions théoriques quelconques : `multibridge::mult_bf_equality` ne traite que des proportions égales, `BayesFactor::proportionBF` que deux modalités, avec une autre loi a priori (logistique sur le logit). La forme close tient en une ligne. **Contrôles** (mesurés le 11 novembre) : égale à `multibridge::mult_bf_equality` à 1e-7 pour des proportions égales (log BF 0,38202 des deux côtés sur (3, 7, 12)) ; égale à l'intégration numérique (`integrate`) à deux modalités pour p quelconque (log BF −0,35933 des deux côtés pour (13, 4) contre p = (0,6 ; 0,4), a = 1 ; 3,01365 pour (7, 2) contre (0,3 ; 0,7), a = 1/2). `multibridge` n'est pas une dépendance : le test qui s'en sert est ignoré s'il n'est pas installé.

**Sensibilité.** a ∈ {1/2, 1, 2} (Jeffreys, uniforme, proportions attendues moins extrêmes), croisé avec les chances de départ 1 contre 2, 1 contre 1, 2 contre 1, comme dans m.bay (D4). Ici, contrairement à `contingencyTableBF`, a < 1 est permis.

**Estimation.** Loi a posteriori Dirichlet(n + a), tirée par `rgamma` (définition de la Dirichlet ; aucun moteur). Par modalité, le rapport π_i / p_i (« x fois plus fréquente que prévu ») avec son intervalle de crédibilité à 95 % (`bayestestR::hdi`) et P(π_i > p_i). Taille de l'écart : w de Cohen a posteriori, w = √Σ (π_i − p_i)² / p_i (Cohen, J. 1988, *Statistical Power Analysis for the Behavioral Sciences*, 2e éd., Erlbaum, chap. 7).

## G3. Indépendance dans une table r × c : REMPLACÉE par G11 (11 novembre, soir)

*Conservée pour l'historique : le moteur `contingencyTableBF` décrit ici n'est plus utilisé pour le facteur de Bayes.*


**Moteur.** `BayesFactor::contingencyTableBF` (Morey et Rouder ; Gunel et Dickey 1974 ; Jamil et al. 2017), a = 1. C'est le moteur déjà employé par m.bay (D17, voie 1), dont l'a priori a été retrouvé et testé à 1e-6.

**Références.**
- Gunel, E., & Dickey, J. M. (1974). Bayes factors for independence in contingency tables. *Biometrika*, 61(3), 545-557.
- Jamil, T., Ly, A., Morey, R. D., Love, J., Marsman, M., & Wagenmakers, E.-J. (2017). Default "Gunel and Dickey" Bayes factors for contingency tables. *Behavior Research Methods*, 49, 638-652. https://doi.org/10.3758/s13428-016-0739-8 (PMC5405059).

**Plan d'échantillonnage : ce que dit la littérature** (Jamil et al. 2017, lu le 11 novembre) :
- *Poisson* : rien n'est fixé, ni les cases ni le total (observation pure : on compte ce qui arrive pendant une durée).
- *Multinomial joint* : le total N est fixé (« interroger 100 personnes »), les deux classements sont aléatoires.
- *Multinomial indépendant* : une marge est fixée (« 30 plants par traitement ») ; le plan le plus courant en expérimentation.
- *Hypergéométrique* : les deux marges fixées (partage à la médiane) ; rare en pratique ; 2 × 2 seulement dans BayesFactor.
- Les BF diffèrent selon le plan, parfois d'une catégorie de Jeffreys (3,04 contre 9,19 sur une même table dans leur simulation), « même asymptotiquement ». Les auteurs ne disent pas quoi faire quand le plan est inconnu ; ils recommandent de vérifier la robustesse.

**Décision : aiguillage par la forme de l'entrée, réglable, et tous les plans montrés.**
1. `sampling = "auto"` (défaut) :
   - table, matrice, `x` et `y`, ou formule `~ A + B` (ou `n ~ A + B` avec une colonne d'effectifs) : **multinomial joint**. C'est la situation symétrique de `GTest`, où aucun des deux facteurs n'est désigné comme explicatif, et le cas le plus fréquent en TP (un échantillon de N individus, classés selon deux critères).
   - formule `y ~ g` (une réponse et un facteur explicatif) : **multinomial indépendant, marge de g fixée** (les groupes constitués d'avance). Même modèle que m.bay voie D17, donc les deux fonctions donnent le même nombre sur les mêmes données (vérifié par un test).
   - `fixed = "rows"` ou `"cols"` donné : multinomial indépendant sur cette marge.
2. `sampling` peut être fixé à `"jointMulti"`, `"indepMulti"`, `"poisson"` ou `"hypergeom"` (2 × 2) par qui connaît le recueil.
3. L'étape 5 donne le facteur de Bayes **sous chacun des autres plans** (poisson, joint, indépendant lignes, indépendant colonnes). Si le plan a été choisi automatiquement et qu'un autre plan plausible change la réponse de côté, la réponse est dite non tranchée, avec la raison (« la réponse dépend de la façon dont les données ont été recueillies »). Le plan hypergéométrique, rare, n'entre pas dans cette comparaison automatique.

**Sensibilité.** a ∈ {1, 2, 4} (BayesFactor refuse a < 1 : sous H0, l'a priori commun est Dirichlet(R(a − 1) + 1)), croisé avec les chances de départ 1:2, 1:1, 2:1 ; ce sont les valeurs de m.bay (D17). Jamil et al. : « si l'on attend de petits effets, on peut augmenter a ».

**Tables à plus de deux dimensions** : refusées avec un message (GTest les refuse aussi).

**Lignes ou colonnes vides** (une modalité déclarée mais jamais observée, fréquent avec un facteur aux niveaux inutilisés) : retirées avant le calcul, avec un « Changement de stratégie » dans le récit. Raison mesurée le 11 novembre : sur une table 2 × 3 à laquelle on ajoute une colonne vide, le log BF passe de 1,89 à 0,56 (joint), de 2,42 à 1,33 (Poisson), de 1,33 à −0,13 (indépendant, lignes fixées) ; une modalité vide ajoute des paramètres à H1 sans aucune donnée, et la réponse dépendrait de niveaux déclarés et non observés. `GTest` ignore aussi ces cases (son G est le même avec ou sans). S'il reste moins de deux lignes ou colonnes, refus.

**Ajustement, modalités à effectif nul** : gardées (une modalité attendue et jamais vue est une information sur l'écart) ; une modalité de proportion théorique nulle est refusée si son effectif est positif, retirée sinon.

## G4. Tables 2 × 2 : rapport des cotes : DÉCIDÉ

Rapport des cotes OR = π₁₁π₂₂ / (π₁₂π₂₁), calculé sur les tirages a posteriori des probabilités des cases donnés par `contingencyTableBF(..., posterior = TRUE)` sous le plan retenu. Il ne dépend que des probabilités conditionnelles, donc a la même lecture sous les quatre plans.
- Médiane et intervalle de crédibilité à 95 % (HDI de `bayestestR::hdi` sur log OR, ramené à l'échelle des cotes : l'HDI de log OR est l'intervalle le plus court sur l'échelle où la loi est la plus symétrique ; ses bornes exponentiées gardent 95 %).
- P(OR > 1).
- Jamil et al. (2017) donnent exactement cette lecture (log OR médian 0,61, IC à 95 % [0,31 ; 0,92] sur leur exemple), et elle rejoint l'estimation classique à grand effectif. Contrôle : test contre `glm` binomial à N = 4000.
- L'HDI est calculé sous H1 (estimation « si un lien existe »), comme dans m.bay ; le bilan le dit.

## G5. Petits effectifs : DÉCIDÉ

- Aucun effectif minimum. Le BF est exact (G2, G3) ; les petits effectifs le rapprochent de 1, ce qui rend la réponse moins tranchée, jamais fausse au sens d'une approximation.
- Le récit compte les cases d'effectif attendu < 5 et dit que le khi-deux et le test G seraient approximatifs là (correction de Williams), pas g.bay.
- Le seuil de m.bay (D14) s'applique : un BF entre 1/3 et 3 (« pas plus qu'une simple mention », Kass et Raftery 1995) donne « non tranché par les données », quelle que soit la probabilité a posteriori.
- Kass, R. E., & Raftery, A. E. (1995). Bayes factors. *JASA*, 90(430), 773-795. https://doi.org/10.1080/01621459.1995.10476572.
- Validation : taux de faux « oui » sous H0 mesurés aussi à N = 10 et 20 (VALIDATION.md).
- **Corrigé par G10 après mesure** : l'affirmation « les petits effectifs rendent la réponse moins tranchée, jamais fausse » tient pour l'ajustement, pas pour l'indépendance dans une table creuse. Le récit ne la dit plus quand la table est creuse.

## G6. Lecture des cases : rapport à l'indépendance a posteriori : DÉCIDÉ

**Ce qu'on remplace.** Les résidus standardisés (O − E)/√E, ou ajustés (Agresti, A. 2013, *Categorical Data Analysis*, 3e éd., Wiley, § 3.3), lus contre ±2.

**Quantité.** Pour chaque case, θ_ij = π_ij / (π_i+ π_+j) : combien de fois la case est plus (θ > 1) ou moins (θ < 1) fréquente que si les deux variables étaient indépendantes. Sa version observée est O/E, la quantité dont le logarithme fait la statistique G (G = 2 Σ O log(O/E)) : c'est donc la lecture du test G case par case. Calculée sur les tirages a posteriori de `contingencyTableBF` (plan retenu), avec :
- la médiane et l'intervalle de crédibilité à 95 % (HDI, `bayestestR::hdi`, sur log θ) ;
- P(θ_ij > 1).

**Précédents.** Albert, J. H. (1997). Bayesian testing and estimation of association in a two-way contingency table. *JASA*, 92(438), 685-693. https://doi.org/10.1080/01621459.1997.10474020 : les termes d'interaction log-linéaires a posteriori, et la détection des cases qui s'écartent du schéma. Kruschke, J. K. (2015). *Doing Bayesian Data Analysis* (2e éd.), Academic Press, chap. 24 : déviations d'interaction cellule par cellule sur un modèle de Poisson log-linéaire, lues par leurs intervalles de crédibilité. log θ_ij est, à une constante près par ligne et par colonne, ce terme d'interaction ; θ est plus lisible pour un étudiant.

**Règle de lecture.** Une case est dite « plus fréquente qu'attendu » si P(θ > 1) ≥ 0,975, « moins fréquente » si P(θ > 1) ≤ 0,025 (l'intervalle à 95 % exclut 1 de ce côté). Le seuil 0,975 plutôt que 0,95 : les cases sont lues toutes ensemble et l'on veut peu de cases signalées sous indépendance ; le taux de tables avec au moins une case signalée sous H0 est mesuré (VALIDATION.md). Les cases ne sont lues que « si un lien existe » ; quand l'indépendance est en tête, le bilan dit qu'aucune case n'est à interpréter. Pour l'ajustement (G2), même lecture modalité par modalité avec π_i / p_i.

## G7. Taille de l'effet : V de Cramér a posteriori : DÉCIDÉ

V = √(φ² / (min(r, c) − 1)), φ² = Σ (π_ij − π_i+π_+j)² / (π_i+π_+j), calculé sur chaque tirage a posteriori des probabilités des cases (Cramér, H. 1946, *Mathematical Methods of Statistics*, Princeton University Press, § 21.9). Médiane et HDI à 95 %.
- V est positif par construction : sa loi a posteriori n'atteint jamais 0, même sous indépendance. C'est une **taille**, pas une preuve de lien ; la preuve est le facteur de Bayes. Le bilan le dit.
- En 2 × 2, V = |φ|. Pour l'ajustement, w de Cohen (G2).
- Aucune étiquette « faible, moyen, fort » : les seuils de Cohen (0,1 ; 0,3 ; 0,5) ne valent que pour un seul degré de liberté et ne sont pas des repères de discipline.

## G8. Hypothèses et bilan : DÉCIDÉ

- Indépendance : « Lien entre A et B » contre « A et B indépendants ». Formule `y ~ g` : « La répartition de y dépend de g » contre « ne dépend pas ».
- Ajustement : « Les proportions de A s'écartent des proportions théoriques » contre « Les proportions de A suivent les proportions théoriques ».
- Probabilités a priori 1/2 chacune (`prior_h1`), classement par probabilité a posteriori, qualificatif de m.bay (seuils 0,6 ; 0,75 ; 0,95 ; 0,993, équivalents à 3, 20, 150 sur le BF à cotes égales, Kass et Raftery 1995), facteur de Bayes, jamais de p-value.
- Mêmes fonctions de mise en forme que m.bay (`sys_bay_story.R`, `sys_bay_layout.R`, `sys_bay_bilan.R`), appelées sans les modifier ; le bilan de g.bay a son propre gabarit (`R/sys_gbay.R`), car celui de m.bay parle de groupes et de différences.
- Même squelette en sept étapes : plan, contrôles, plan d'échantillonnage (loi des données), point de départ, solidité, où se trouve l'écart, bilan.
- `return = FALSE` : un nombre, P(lien) ou P(écart), sans les tirages a posteriori (rapide, pour les simulations).

## G10. Tables creuses : soutien d'un lien plafonné à « faible » : RETIRÉE avec G11 (11 novembre, soir)

*Le plafond corrigeait un symptôme du moteur de G3 ; G11 en retire la cause, et la recalibration ne montre plus de faux « oui » à corriger (VALIDATION.md, § 9). Texte conservé pour l'historique.*


**Défaut mesuré** (VALIDATION.md, § 2 et 3). Sous indépendance, dans une table où les effectifs attendus sont petits, le facteur de Bayes de `contingencyTableBF` à a = 1 penche vers un lien :
- P > 0,75 dans 32 à 36 % des tables 3 × 4 à N = 10-20 et dans 67 à 73 % des tables 4 × 5 à N = 10-20 (plan joint) ;
- une fois la sensibilité appliquée (a = 1, 2, 4 ; chances de départ ; plans), le bilan dit encore « un lien existe » dans 11,9 % des tables 4 × 5 à N = 50, dont 11,6 % « très probablement » ;
- dans les tables bien remplies (attendu moyen ≥ 8), le bilan dit « un lien existe » dans 6,6 % des cas au plus (2 × 2, N = 50), « très probablement » dans 1,2 % au plus.

L'a priori uniforme (Dirichlet(1) sur r × c probabilités) met beaucoup de masse sur des tables très déséquilibrées ; une table creuse, déséquilibrée par le hasard, y ressemble. La sensibilité de ce test à la concentration de la Dirichlet est connue : Good, I. J., & Crook, J. F. (1987). The robustness and sensitivity of the mixed-Dirichlet Bayesian test for "independence" in contingency tables. *The Annals of Statistics*, 15(2), 670-693. https://doi.org/10.1214/aos/1176350368.

**Décision.** Une table est dite **creuse** selon la règle de Cochran, sous sa forme courante : plus de 20 % des cases d'effectif attendu inférieur à 5, ou une case sous 1 (Cochran, W. G. 1954. Some methods for strengthening the common χ² tests. *Biometrics*, 10(4), 417-451. https://doi.org/10.2307/3001616 ; la formulation d'origine est un peu plus souple, voir Kroonenberg et Verbeek 2018, *The American Statistician*, 72(2), 175-183). Dans une table creuse :
- le soutien d'un lien est plafonné à « faible », comme m.bay plafonne déjà le soutien quand les variances diffèrent par case (D19) ; la probabilité reste affichée ;
- le verdict devient « plutôt un lien, sans conclure » ; le récit (étape 2, « Prudence ajoutée ») et le bilan disent pourquoi ; la phrase de compte rendu ajoute « table creuse, soutien faible » ;
- l'indépendance en tête n'est pas touchée : le défaut mesuré va dans le sens d'un lien.

**Pourquoi cette règle et pas une autre.** Elle est connue des étudiants (c'est la condition du khi-deux) et se décide avant de regarder la réponse. Avec elle, le bilan dit « un lien existe » sous indépendance dans 6,6 % des tables au plus, et dans 0,8 % au plus des tables 3 × 4 et 4 × 5 (VALIDATION.md, § 3).

**Coût mesuré.** Dans une table creuse, un vrai lien n'est jamais dit « existe » : sous un lien moyen en 3 × 4, la part de « oui » passe de 81 % à 0 % à N = 50, de 97 % à 30 % à N = 100, de 100 % à 98 % à N = 200 ; en 2 × 3 (lien fort), de 29 % à 17 % à N = 50, sans changement à N = 100. Le bilan dit alors « plutôt un lien, sans conclure », avec la probabilité.

**Ajustement.** Pas de plafond : sous H0, P > 0,75 dans 6,4 % des cas au plus, même à N = 10 (5 modalités), et P > 0,95 dans 0,7 % au plus.

**Suite possible, non faite.** Un a priori mieux adapté aux tables creuses (mélange sur la concentration, Good et Crook 1987) rendrait de la puissance. Il faudrait le coder à la main (pas de paquet établi) et le calibrer ; c'est une décision à part.

## G11. Indépendance : test log-linéaire de l'interaction, BIC avec n = N : DÉCIDÉ ET FAIT (11 novembre, soir), remplace G3 et G10

**Origine.** Point bloquant de la session parente, après l'avis d'un expert, recalculé par elle et par moi. Sous indépendance vraie, `contingencyTableBF` (a = 1) annonce un lien fort (BF > 20) dans :
- 82,7 % des tables 4 × 5 à N = 50 (plan Poisson), 11,0 % (joint), 2,7 % (indépendant, lignes fixées), d'après la mesure de la session parente ;
- 16 % des tables 4 × 5 à N = 50 à marges égales, et 4,1 % des tables 3 × 4 à N = 60, que la règle de Cochran ne classe pas comme creuses, d'après l'expert. Le plafond de G10 ne couvrait donc pas le défaut.

Sur une même table, le BF allait de 0,28 à 6,74 selon le plan.

**Cause.** L'a priori de Gunel et Dickey ne traite pas les marges de la même façon sous H0 et sous H1. Sous H1, les marges découlent d'une Dirichlet(a) sur les r × c cases ; sous H0, elles ont chacune leur propre Dirichlet (BayesFactor : `ldirich(ar - (J - 1))`, `ldirich(ac - (I - 1))`, code lu le 11 novembre). Le facteur de Bayes compare donc aussi deux lois des marges : des marges équilibrées, plus probables sous la Dirichlet des cases que sous celles des marges, comptent comme une preuve de lien. Le défaut grandit avec le nombre de cases et dépend du plan, qui change la loi des marges.

**Décision.** Le lien est l'interaction d'un modèle log-linéaire de Poisson, et les marges sont libres sous les deux hypothèses :
- H0 : log μ_ij = λ + α_i + β_j ;
- H1 : on ajoute γ_ij, soit (r − 1)(c − 1) paramètres.

Le facteur de Bayes est −BIC/2 de chaque modèle, avec n = N, le nombre d'individus et non le nombre de cases (approximation d'information unitaire : Schwarz 1978 ; Kass, R. E., & Wasserman, L. 1995, A reference Bayesian test for nested hypotheses and its relationship to the Schwarz criterion, *JASA*, 90(431), 928-934 ; Raftery 1995). Pour les modèles log-linéaires, Raftery (1986) l'écrit BIC = G² − ddl·log N (Raftery, A. E. 1986, Choosing models for cross-classifications, *American Sociological Review*, 51(1), 145-146 ; forme reprise dans le cours STAT 504 de Penn State, consulté le 11 novembre). D'où :

  log BF₁₀ = (G − ddl·log N) / 2, avec G la statistique du test G (`GTest`, sans correction).

**Calcul.** `stats::glm(n ~ lignes + colonnes, family = poisson)` sur les effectifs des cases. Sa déviance est G, puisque le modèle saturé a une déviance nulle ; la pénalité prend log N. `stats::BIC()` n'est pas appelé tel quel : sur une table agrégée, il compterait les cases (r × c) au lieu des N individus. Contrôles (tests) : égal à `(GTest(x)$statistic − ddl·log N)/2` à 1e-8 ; en 2 × 2 et 2 × 3 avec une réponse, égal à −BIC/2 calculé par `stats::BIC` sur les `glm` logistiques ou log-linéaires à une ligne par individu (la méthode de m.bay, D25).

**Pourquoi le BIC plutôt que `BAS::bas.glm`.**
- C'est la méthode de m.bay pour les réponses discrètes (D17 voies 2 et 3, D25, D26) : la famille reste cohérente, et la session parente demande la même correction à m.bay.
- BAS sur une table agrégée prendrait le nombre de lignes, c'est-à-dire de cases, comme taille d'échantillon dans ses a priori g (Bayarri et al. 2012). Il faudrait développer la table en N lignes et passer par un modèle multinomial, que `bas.glm` ne propose pas pour une réponse à plus de deux modalités.
- Le BIC se lit sur G, ce qui relie directement g.bay au test G que les étudiants connaissent.

**Pourquoi le plan d'échantillonnage n'a plus d'importance.** Sous les plans Poisson, multinomial et multinomial à une marge fixée, le noyau de vraisemblance du paramètre d'interaction et ses estimations sont les mêmes (Birch 1963 ; Agresti 2013, modèles log-linéaires) : G et ddl ne changent pas, et N non plus. Les arguments `sampling` et `fixed` sont donc retirés, ce qui supprime le problème « BF de 0,28 à 6,74 selon le plan ».

**Ce que le BIC coûte.** C'est une approximation d'ordre O(1). Elle penche vers H0, d'autant plus que la table est grande (pénalité ddl·log N). La puissance est mesurée au § 9 de VALIDATION.md. La sensibilité ne porte que sur les chances de départ : l'information unitaire n'a pas d'échelle à faire varier, comme dans D25.

**`contingencyTableBF` garde un seul usage : l'estimation.** Son tirage a posteriori, plan joint, donne une loi de Dirichlet(n + 1) sur les cases ; il sert aux cases, au V de Cramér et au rapport des cotes (G4, G6, G7). Cet usage est calibré : l'intervalle à 95 % du rapport des cotes couvre la vraie valeur dans 93 à 96 % des cas (VALIDATION.md, § 6). Son facteur de Bayes n'est plus utilisé.

**G10 retiré.** Le plafond des tables creuses corrigeait un symptôme du moteur de G3. Avec G11, aucun scénario ne dépasse 0,4 % de BF > 20 sous H0 (§ 9), et le plafond n'a plus d'objet.

**Ajustement (G2) inchangé.** H0 y est un vrai point (les proportions données) ; il n'y a pas de marges à traiter, et la calibration de G2 tient (§ 2).

**Pour m.bay.** La même correction est demandée à `sys_bay_contingency.R`. La formulation ci-dessus vaut pour les deux fonctions : test de l'interaction, marges libres, −BIC/2 avec n = N.

## G9. Validation : PRÉVUE, résultats dans VALIDATION.md

- Accord : `g.bay()` = `contingencyTableBF` sous chaque plan, à 1e-10 ; forme close d'ajustement = `multibridge` (proportions égales) et `integrate` (deux modalités) ; `y ~ g` = `m.bay(y ~ g)`.
- Sous H0 : taux de « oui » à P > 0,75 et P > 0,95, ajustement et indépendance, petits et grands effectifs.
- Puissance comparée au test G (`GTest`, p < 0,05) sous des écarts de taille croissante.
- Cases signalées sous H0 (G6).
