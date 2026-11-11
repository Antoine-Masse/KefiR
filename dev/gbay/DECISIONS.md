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

## G3. Indépendance dans une table r × c : DÉCIDÉ

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

## G9. Validation : PRÉVUE, résultats dans VALIDATION.md

- Accord : `g.bay()` = `contingencyTableBF` sous chaque plan, à 1e-10 ; forme close d'ajustement = `multibridge` (proportions égales) et `integrate` (deux modalités) ; `y ~ g` = `m.bay(y ~ g)`.
- Sous H0 : taux de « oui » à P > 0,75 et P > 0,95, ajustement et indépendance, petits et grands effectifs.
- Puissance comparée au test G (`GTest`, p < 0,05) sous des écarts de taille croissante.
- Cases signalées sous H0 (G6).
