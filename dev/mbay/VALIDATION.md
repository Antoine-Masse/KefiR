# Validation de m.bay

**Le §7 valide `m.bay()` à un facteur. Les §1 à 6 documentent les moteurs candidats et les essais antérieurs (5 et 6 novembre) ; ils ne valident pas `m.bay()`.** R 4.6.1, BayesFactor 0.9.12-4.8, Windows 11. Scripts dans `scripts/`, sorties dans `scripts/resultats/`. Détails : `rapport_moteurs.md`, `agent_2x_calibration_rapport.md`, `RELECTURE_ADVERSE_2026-11-05.md`.

Légende : **[recalculé]** = refait séparément par le coordinateur ; **[agent]** = exécuté par un agent, sortie présente, non refait.

## 1. Comparaisons à un calcul indépendant

| Résultat | Écart | Script | Statut |
|---|---|---|---|
| `ttestBF` contre intégrale JZS programmée | 1e-11 (autre jeu, formule réécrite) | `agent_09` | [recalculé] |
| `ttestBF` contre JZS, trois jeux | 8e-13 à 9e-11 | `agent_02` | [agent] |
| `anovaBF` contre intégrale Zellner-Siow, mêmes a priori | 4 à 5e-8 | `agent_02` | [agent] |
| `y~A+B` à 2 g, double intégrale | 1,5 % (propErr 0,9 %) | `agent_02` | [agent] |
| `y~x`, forme de Liang | 2e-9 | `agent_02` | [agent] |
| NIW propre, p=1 contre integrate et t multivariée | 4e-15, 7e-15 | `agent_05` | [agent] |
| NIW propre, p=2 contre Monte-Carlo | -41,4982 contre -41,5016 ± 0,0049 | `agent_05` | [agent] |
| Socle `.bay_inclusion` contre `bayestestR`, withmain, 2 jeux | 4e-16 sur log BF | `agent_10` | [agent] |

**Non comparable** : rstanarm + bridgesampling BF=0,207 contre anovaBF BF=3,46 ; a priori différents. Plus de 2 g sans vérification indépendante. Comparaison Stan à a priori identiques : non faite.

Écart 0,4547 / 0,4415 sur A:B apparié : différence de définition, pas de Monte-Carlo : `bayestestR` garde dans M+ les modèles qui violent la marginalité quand on lui passe `whichModels="all"`. Mécanisme vérifié à 5,6e-17 [agent].

## 2. Calibration

### Sous l'a priori (un facteur, test de calcul)

8 scénarios, N=2500 chacun, π=1/2 et 0,2, k=2 et 3, n=10 et 30 : aucune déviation détectée, moyennes |z| ≤ 1,27 [agent, `agent_2x_02_sim.R`].
Recalculé séparément, k=3, n=10, N=2000, `agent_14` [recalculé] :

| P annoncée, classe | Moyenne annoncée | Fréquence de H1 | n |
|---|---|---|---|
| [0 ; 0,2] | 0,186 | 0,171 | 258 |
| (0,2 ; 0,4] | 0,272 | 0,263 | 852 |
| (0,4 ; 0,6] | 0,493 | 0,512 | 215 |
| (0,6 ; 0,8] | 0,694 | 0,721 | 147 |
| (0,8 ; 1] | 0,967 | 0,960 | 526 |

Trouvé au passage : **`integrate()` échoue pour F très grand** et BayesFactor renvoie NA (reproduit dans `agent_14`, moyenne NA). 7 jeux sur 10 000 chez l'agent. `m.bay` devra gérer ce cas.

### Effets fixes (a priori faux à la génération), N=800, δ=(max-min)/σ [agent]

- δ=0 : P>0,95 ≤ 0,4 %, aov p<0,05 de 3,9 à 6,2 %.
- Soutien à H0 faible : P<0,05 observé seulement à k=3, n=100 (31,5 %). P minimale 0,284 (k=2, n=10) à 0,035 (k=3, n=100).
- δ=0,5, k=2, n=100 : P>0,95 62,7 %, p<0,05 94,9 %. Moins décisif que aov, à seuils non équivalents.

**Comparaison de seuils** : P>0,5 équivaut à p<0,084 et P>0,95 à p<0,0019 pour un facteur, k=3, n=15 (`agent_03`) : comparer P>0,5 à p<0,05 n'a pas de sens.

### Écarts au modèle

- t3 et log-normale, n=30, δ=0 : P>0,95 à 0 % [0 ; 0,5] [agent]. Petits n non testés.
- Hétéroscédasticité, H0 : n=(5,5,40), sd=(4,4,1), **P>0,95 dans 20,0 %** [18,6 ; 21,5], N=3000 ; aov 47,2 %, Welch 5,7 %, Kruskal 20,4 % [agent, `agent_21`] ; 20,2 % [18,2 ; 22,3], N=1500, autre graine [recalculé, `agent_23`]. Les 16 % (`agent_12`) et 25,5 % (relecteur) sont deux tirages de part et d'autre : mêmes données, même calcul, différence nulle [agent, `agent_21`]. Autre plan, SD 3 sur n=10 de 10,30,50 : P>0,95 8,0 % [5,4 ; 11,6], aov 26,3 % [agent].

### Sensibilité à r (0,25 ; 0,5 ; 1)

P>1/2 change de côté dans 4,5 à 9,8 % des jeux à δ=0 et 17 à 23 % à δ=0,5 [agent, N=600 ou 300]. Étape incomplète : k=2 à n=10 et 100 non exécutés.

### A*B, inclusion appariée : **non calibrée sous un générateur unique** [agent, 1521 jeux stratifiés]

- Générateur à hérédité forte : A:B annonce 0,350 pour une fréquence de 0,125 (z=51).
- Uniforme sur 5 modèles : A:B 0,379 contre 0,200 ; A et B z ≈ -8 à -10.
- Calibrées seulement chacune dans son monde : monde additif pour A et B, monde {A+B, A*B} pour A:B, et chaque monde rend les autres termes faux.
- Recherche post hoc sur 8854 mélanges : 2 seuls à p>0,05, collés à « nul ou A*B » (exact p=0,046). Argument algébrique, pas démonstration complète.

### Croisement 2×2 pur : P(A) apparié 0,215 alors que P(A:B)=1 [recalculé, `agent_11`]

Voir DECISIONS D2.

## 3. Socle `R/sys_bay_core.R`

`tests/testthat/test-bay-core.R` : FAIL 0, ERROR 0, 237 assertions, relancé par le coordinateur [recalculé]. ASCII : vérifié [agent]. Test de `.vbse` limité à cause du bug signalé au JOURNAL.
Non-régression `test-m.test.R` : 20 tests, 52 assertions, 0 échec [recalculé]. `test-manova.R`, `test-valreg.R` non lancés.

## 4. Convergence (R-hat, ESS)

Sans objet pour les voies BayesFactor analytiques. Stan non retenu à ce stade : seulement mesuré, `stan_glm` R-hat 1,014, ESS 383 ; `stan_lmer` R-hat 1,003, ESS 290 ; `stan_mvmer` R-hat 1,035, ESS 52 (insuffisant) [agent].

## 5. Non vérifié

ROPE f et MANOVA en √d (relecteur seul) ; toute voie non gaussienne ; hypothèses de sphéricité ; BAS ; brms avec Rtools.

## 6. Séance du 6 novembre

### Multiplicité 2³ sous H0, recoupée sur `anovaBF` (au moins un terme sur 7 au-dessus du seuil)

| Loi, seuil | N=16 (150 jeux, `agent_20`) | N=40 (226 jeux, `agent_20`) | N=40 (300 jeux, `agent_22`) |
|---|---:|---:|---:|
| Apparié P>0,5 | 0,740 [0,662 ; 0,808] | 0,580 [0,512 ; 0,645] | 0,603 [0,546 ; 0,659] |
| Apparié P>0,95 | 0,013 | 0,009 | 0,017 [0,005 ; 0,038] |
| Uniforme P>0,5 | 0,747 | 0,429 | 0,487 [0,429 ; 0,545] |
| Hérédité de D2 P>0,5 | non fait | non fait | 0,253 [0,205 ; 0,307] |
| Hérédité de D2 P>0,95 | non fait | non fait | 0,000 [0 ; 0,012] |
| aov, au moins un p<0,05 | 0,247 | 0,292 | 0,347 [0,293 ; 0,403] |

Le calcul maison du relecteur est confirmé. `agent_20` code l'hérédité en π^k (P(A:B)=1/16) : **ses lignes « hérédité » ne portent pas sur D2 et ne sont pas reprises**. `agent_22`, lancé par le coordinateur, emploie la loi de D2 : P(A:B)=1/8, P(A:B:C)=1/128. Erreur annoncée par BayesFactor : médiane 1,8 %, maximum 89 à 97 % sur un modèle ; effet sur les taux non mesuré. Forme finale « ligne 1, ligne 2 » : non simulée.

### `integrate()` à F très grand : résolu pour le calcul à un facteur

`.bay_bf_oneway_zs` : écart à `integrate` reparamétré et à une multiprécision inférieur à 3e-14 à F ≤ 1e6, inférieur à 8e-9 jusqu'à 1e16, 2,5e-3 au plus à 1e24-1e28 (cause supposée, l'arrondi de SSE, non testée) [agent, `agent_24_v1v2`]. Contre `anovaBF` : écart sous 1e-5 sur 6 plans × 2 r [tests relancés] ; à F=1e6, n=(5,5,40), BayesFactor rend 193,11 (erreur annoncée 2,0), zs 240,14 [agent]. Multi-facteurs : le problème reste entier.

### Voie à variances par groupe `.bay_bf_oneway_hetero` : **à l'essai, non validée**

Calcul indépendant.
- Grille déterministe, 9 jeux k=2 et 3 : aucun écart au-delà de 2 error_pct, 5 sur 9 au-delà de 1 [agent, `agent_24_v1v2`].
- Grille écrite par le coordinateur, k=2, 3 jeux : -0,7968/-0,7972 ; 0,1864/0,1864 ; -0,2840/-0,2838 [recalculé, `agent_29`].
- error_pct couvre correctement : rapport écart-type/error_pct 0,80 à 1,22 sur 30 graines [agent].

Simulation H0, n=(5,5,40), mêmes jeux que `agent_21`, N=1000 [agent, `agent_25_v3abc`] :

| | P>0,5 | P>0,95 | P médiane |
|---|---:|---:|---:|
| sd 4,4,1, hetero | 0,163 [0,141 ; 0,187] | **0,007 [0,003 ; 0,014]** | 0,290 |
| sd 4,4,1, gaussien | 0,544 | 0,207 [0,182 ; 0,233] | 0,552 |
| sd 1,1,4, hetero | 0,075 | 0,000 [0 ; 0,004] | 0,181 |
| sd 1,1,4, gaussien | 0,000 | 0,000 | 0,207 |

Recoupé : P>0,95 0/300 [0 ; 1,2 %], P>0,5 0,147, graine différente [recalculé, `agent_29`]. Le défaut de 20 % à P>0,95 disparaît sur **ce plan seulement**. P>0,5 reste à 16 %.

Mesures de faible résolution, à refaire avant toute conclusion [agent] :
- Homoscédastique, H0, N=60 par case : hetero et zs voisins, corrélation des log BF 0,84 à 1,00.
- Puissance, **N=16-17 par case** : aucune perte mesurable en homoscédastique. Avec variances inégales et grand n sur le groupe de petit écart-type, hetero nettement sous zs et Welch (k=3, d=1, n=30 : 0,25 contre 0,94 et 0,94 à P>0,95 ou p<0,05). Ce point est à examiner.
- Calibration sous a priori tronqué, N≈250 par case : un seul |z|>2 sur 18 classes ; (0,8 ; 1] : 0,958 annoncé contre 0,902 observé, n=(10,10), plage [1/20 ; 20]. Effet de la troncature : N=7-8, indicatif seulement.
- Sensibilité à r (0,25 ; 1), N=72 : P change de côté dans 3 à 7 % des jeux.
- Cas dégénérés, 10 graines : aucune erreur ; k=20 3 à 5,5 s par appel ; sd=1e-12 5,3 s.

Sur le modèle (méthodologie, `methodo_20261106/`, un seul jeu, Python) : un troisième groupe sans rapport (sd × 100) fait passer le BF de 46,4 à 2,9 avec τ arithmétique ; géométrique 43 à 58, harmonique jusqu'à 72. Beta(1,5 ; 1,5) au lieu du plat sur le rapport change le BF d'un facteur 1,02 à 1,23 ; Beta(5,5) peut faire changer de côté de 1. **La définition de τ est donc ouverte.**

### Défauts du modèle, relecture adverse du 6 novembre

Le calcul est juste, **le modèle ne l'est pas** (`RELECTURE_ADVERSE_2026-11-06.md`). Reproduit [recalculé, `agent_30`] :

| Cas | log BF |
|---|---:|
| n=(10,10,10,10,3), dernier groupe à 50σ : hetero, zs, Welch p | -6,92 ; +85,6 ; 3,3e-15 |
| Même plan, écart 5, 50, 500, 5000 | -1,50 ; -6,92 ; -12,25 ; -16,88 |
| 2 groupes seuls (d=1,3 ; n=10) | +1,64 |
| Les mêmes plus un 3e groupe de sd S = 1, 10, 100, 1e4, 1e6 | +2,63 ; +1,12 ; -1,19 ; -5,81 ; -10,42 |
| `zs` sur `c(.3,.3,.3,.1+.2,.1+.2,.1+.2)`, données constantes à l'arrondi près | Inf (à corriger) |

Par le relecteur seul, non refait : écart de 3σ sur ce plan, BF<1/3 dans 67 % des 150 jeux, Welch 83 % de détection ; n=(2,2,40), P(BF>19)=6 % (9/150) ; ESS et error_pct trompeurs à r=0,005 ; grille en g tronquée à k≥200.

### Tests

`test-bay-oneway.R` : 24 tests, 80 assertions, 0 échec ; `test-bay-core.R` : 237, 0 échec ; `test-m.test.R` : 52, 0 échec [recalculé, `agent_26_lance_tests_bay`]. Sources ASCII contrôlées [recalculé].

## 7. `m.bay()` à un facteur, 2 et k groupes, lois comparées (7 novembre)

Remplace le §6 pour cette voie. Tout se reproduit par les commandes du JOURNAL, qui dit aussi sur quelle version de code a tourné chaque lancement. Sorties : `scripts/resultats/agent_40.txt`, `agent_41_A.txt`, `agent_41_A2.txt`, `agent_41_B.txt`, `agent_41_C.txt`, `agent_41_D.txt`, `agent_42.txt`, `agent_43.txt`, `agent_44.txt`, tranches dans `agent_41_cas/`.

### Calcul indépendant

| Comparaison | Écart | Script |
|---|---|---|
| Les 8 modèles à k=2, 3 jeux, contre une quadrature écrite à part (`dnorm`, `dt`, `dcauchy`) | ≤ 0,02 sur log m(D), tous dans ±3 erreurs annoncées, normale exacte à 1e-4 | `agent_40` |
| Seconde forme sous loi t à variance commune, k=3, 2 jeux × 3 groupes, quadrature 4D | ≤ 0,022, ±2,4 erreurs annoncées | `agent_42` |
| Seconde forme, normale : `anovaBF` sur « s contre le reste » | < 1e-4 (test) | `test-bay-laws.R` |
| Première forme, normale : `.bay_bf_oneway_zs`, `anovaBF` ; k=2 : `ttestBF` | < 1e-4 (tests) | `test-bay-laws.R` |
| Noyau et constantes à k=3 contre intégration brute | ≤ 2e-3 (relecteur, sur le code d'avant le mélange) | `relecture/agent_rev3_03` |

**Non comparé** : t à variance par groupe à k ≥ 3, et la seconde forme à variance par groupe (normale ou t) : trop de dimensions pour une quadrature. Ces parties reposent sur la calibration ci-dessous.

### Calibration sous l'a priori

Générateur : l'a priori lui-même, loi comprise, deux formes d'effet comprises, N=600 par cas (`agent_41 A`). Test de calcul, pas de vérité du modèle : il est bon par construction si le calcul est juste, c'est ce qu'il contrôle.

| Plan | Moyenne de P(H1) | Fréquence de H1 | Moyennes de P(lois), attendu 0,25 environ | Lois écartées |
|---|---:|---:|---|---:|
| k=2, n=10 | 0.485 | 0.480 | 0.248 0.253 0.254 0.245 | 2 |
| k=2, n=30 | 0.489 | 0.493 | 0.253 0.240 0.250 0.258 | 3 |
| k=3, n=10 | 0.514 | 0.510 | 0.237 0.246 0.257 0.264 | 9 |
| k=3, n=30 | 0.519 | 0.515 | 0.247 0.229 0.251 0.274 | 2 |
| k=5, n=10 | 0.513 | 0.507 | 0.260 0.253 0.247 0.245 | 8 |
| k=5, n=30 | 0.481 | 0.495 | 0.226 0.256 0.271 0.248 | 4 |

Tranches de probabilité annoncée : voir `agent_41_A.txt`. Une seule |z| > 2,1 sur 36 tranches à N=600, à P > 0,95 (0,998 annoncé, 0,990 observé, k=3, n=30) et une à (0,8 ; 0,95] (k=2, n=10, z=−2,05). La tranche suspecte du premier lancement (k=2, n=30, (0,8 ; 0,95], 0,869 annoncé, 0,971 observé, z=3,64 avec la variance observée) ne se reproduit pas : autre graine, N=2000, 0,884 annoncé, 0,823 observé, z=−1,73, sens inverse (`agent_41_A2.txt`).

**Lois écartées (`ESS` < 30).** Au plus 1,5 % des jeux (9 sur 600 à k=3, n=10). À k=2 et k=3, tous ces jeux ont P=1 : effet énorme, sans conséquence sur la réponse. À k=5, P descend à 0,14 parmi eux : **la réponse peut alors en dépendre**, ce que le bilan annonce (« la comparaison est incomplète »). Le poids des lois restantes peut en être faussé.

### Effets fixes, N=300 par cas (`agent_41 B`)

Colonnes : `m.bay` ; gaussien seul (variance commune, `anovaBF`) ; tests classiques. **P > 0,95 et p < 0,05 ne sont pas des seuils équivalents** : P > 0,95 correspond plutôt à p ≈ 0,002 à ces effectifs (§2).

| Scénario | m.bay P>0,5 | m.bay P>0,95 | Gaussien P>0,95 | aov p<0,005 | Welch p<0,05 | Welch p<0,005 |
|---|---:|---:|---:|---:|---:|---:|
| hetero (5,5,40) sd(4,4,1) H0 | 0.267 | 0.013 | 0.190 | 0.260 | 0.037 | 0.003 |
| hetero (5,5,40) sd(1,1,4) H0 | 0.057 | 0.003 | 0.000 | 0.000 | 0.043 | 0.013 |
| hetero (10,30,50) sd(3,1,1) H0 | 0.110 | 0.020 | 0.107 | 0.167 | 0.053 | 0.020 |
| hetero (2,2,40) sd(4,4,1) H0 | 0.537 | 0.030 | 0.433 | 0.507 | 0.090 | 0.017 |
| (2x6) sd 1 H0 | 0.223 | 0.000 | 0.003 | 0.003 | 0.113 | 0.020 |
| hetero (5,5,40) sd(4,4,1) effet 4 groupe 1 | 0.823 | 0.153 | 0.780 | 0.830 | 0.300 | 0.053 |
| hetero (5,5,40) sd(4,4,1) effet 8 groupe 1 | 1.000 | 0.580 | 0.993 | 1.000 | 0.860 | 0.337 |
| hetero (30,30,30) sd(1,1,3) d=.8 entre 1 et 2 | 0.833 | 0.337 | 0.063 | 0.107 | 0.807 | 0.467 |
| normal k=2 n=10 d=0.0 | 0.097 | 0.000 | 0.003 | 0.003 | 0.033 | 0.003 |
| normal k=2 n=10 d=0.5 | 0.280 | 0.023 | 0.027 | 0.050 | 0.190 | 0.047 |
| normal k=2 n=10 d=0.8 | 0.533 | 0.033 | 0.037 | 0.087 | 0.373 | 0.083 |
| normal k=3 n=10 d=0.0 | 0.120 | 0.010 | 0.010 | 0.017 | 0.070 | 0.023 |
| normal k=3 n=10 d=0.5 | 0.303 | 0.017 | 0.023 | 0.033 | 0.167 | 0.043 |
| normal k=3 n=10 d=0.8 | 0.540 | 0.087 | 0.093 | 0.137 | 0.363 | 0.117 |
| normal k=2 n=30 d=0.0 | 0.047 | 0.000 | 0.000 | 0.003 | 0.027 | 0.003 |
| normal k=2 n=30 d=0.5 | 0.567 | 0.093 | 0.093 | 0.190 | 0.513 | 0.187 |
| normal k=2 n=30 d=0.8 | 0.897 | 0.377 | 0.377 | 0.520 | 0.850 | 0.517 |
| normal k=3 n=30 d=0.0 | 0.063 | 0.000 | 0.000 | 0.007 | 0.057 | 0.003 |
| normal k=3 n=30 d=0.5 | 0.523 | 0.100 | 0.107 | 0.183 | 0.483 | 0.160 |
| normal k=3 n=30 d=0.8 | 0.900 | 0.513 | 0.530 | 0.660 | 0.883 | 0.640 |
| normal k=2 n=100 d=0.0 | 0.050 | 0.003 | 0.003 | 0.013 | 0.050 | 0.013 |
| normal k=2 n=100 d=0.5 | 0.970 | 0.640 | 0.640 | 0.797 | 0.970 | 0.797 |
| normal k=2 n=100 d=0.8 | 1.000 | 0.993 | 0.993 | 0.997 | 1.000 | 0.997 |
| normal k=3 n=100 d=0.0 | 0.033 | 0.000 | 0.000 | 0.003 | 0.047 | 0.003 |
| normal k=3 n=100 d=0.5 | 0.960 | 0.717 | 0.720 | 0.867 | 0.963 | 0.853 |
| normal k=3 n=100 d=0.8 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| t3 k=2 n=10 d=0.0 | 0.087 | 0.000 | 0.000 | 0.013 | 0.050 | 0.007 |
| t3 k=2 n=10 d=0.5 | 0.457 | 0.043 | 0.027 | 0.060 | 0.307 | 0.060 |
| t3 k=2 n=10 d=0.8 | 0.737 | 0.177 | 0.120 | 0.200 | 0.503 | 0.190 |
| t3 k=3 n=10 d=0.0 | 0.087 | 0.000 | 0.000 | 0.003 | 0.040 | 0.000 |
| t3 k=3 n=10 d=0.5 | 0.363 | 0.020 | 0.010 | 0.027 | 0.217 | 0.037 |
| t3 k=3 n=10 d=0.8 | 0.703 | 0.207 | 0.150 | 0.190 | 0.510 | 0.233 |
| t3 k=2 n=30 d=0.0 | 0.047 | 0.003 | 0.000 | 0.003 | 0.060 | 0.003 |
| t3 k=2 n=30 d=0.5 | 0.730 | 0.247 | 0.157 | 0.247 | 0.567 | 0.243 |
| t3 k=2 n=30 d=0.8 | 0.997 | 0.860 | 0.610 | 0.723 | 0.893 | 0.720 |
| t3 k=3 n=30 d=0.0 | 0.030 | 0.003 | 0.000 | 0.000 | 0.030 | 0.000 |
| t3 k=3 n=30 d=0.5 | 0.713 | 0.300 | 0.153 | 0.243 | 0.527 | 0.257 |
| t3 k=3 n=30 d=0.8 | 0.990 | 0.870 | 0.583 | 0.667 | 0.880 | 0.693 |
| t3 k=2 n=100 d=0.0 | 0.030 | 0.000 | 0.000 | 0.000 | 0.040 | 0.000 |
| t3 k=2 n=100 d=0.5 | 0.993 | 0.947 | 0.677 | 0.800 | 0.937 | 0.800 |
| t3 k=2 n=100 d=0.8 | 1.000 | 1.000 | 0.987 | 0.990 | 0.997 | 0.990 |
| t3 k=3 n=100 d=0.0 | 0.007 | 0.003 | 0.000 | 0.000 | 0.033 | 0.003 |
| t3 k=3 n=100 d=0.5 | 1.000 | 0.983 | 0.720 | 0.827 | 0.950 | 0.843 |
| t3 k=3 n=100 d=0.8 | 1.000 | 1.000 | 0.980 | 0.987 | 0.997 | 0.987 |
| lognormal k=2 n=10 d=0.0 | 0.093 | 0.003 | 0.000 | 0.003 | 0.037 | 0.003 |
| lognormal k=2 n=10 d=0.5 | 0.450 | 0.060 | 0.023 | 0.043 | 0.240 | 0.040 |
| lognormal k=2 n=10 d=0.8 | 0.777 | 0.267 | 0.157 | 0.240 | 0.520 | 0.213 |
| lognormal k=3 n=10 d=0.0 | 0.100 | 0.000 | 0.000 | 0.007 | 0.037 | 0.000 |
| lognormal k=3 n=10 d=0.5 | 0.497 | 0.073 | 0.037 | 0.077 | 0.287 | 0.060 |
| lognormal k=3 n=10 d=0.8 | 0.853 | 0.350 | 0.147 | 0.243 | 0.573 | 0.260 |
| lognormal k=2 n=30 d=0.0 | 0.050 | 0.000 | 0.000 | 0.000 | 0.050 | 0.000 |
| lognormal k=2 n=30 d=0.5 | 0.857 | 0.437 | 0.183 | 0.283 | 0.587 | 0.280 |
| lognormal k=2 n=30 d=0.8 | 0.993 | 0.910 | 0.567 | 0.667 | 0.870 | 0.667 |
| lognormal k=3 n=30 d=0.0 | 0.050 | 0.000 | 0.003 | 0.013 | 0.063 | 0.003 |
| lognormal k=3 n=30 d=0.5 | 0.893 | 0.470 | 0.137 | 0.210 | 0.617 | 0.253 |
| lognormal k=3 n=30 d=0.8 | 1.000 | 0.963 | 0.587 | 0.707 | 0.950 | 0.793 |
| lognormal k=2 n=100 d=0.0 | 0.030 | 0.007 | 0.003 | 0.007 | 0.047 | 0.007 |
| lognormal k=2 n=100 d=0.5 | 1.000 | 0.980 | 0.643 | 0.773 | 0.933 | 0.773 |
| lognormal k=2 n=100 d=0.8 | 1.000 | 1.000 | 0.950 | 0.973 | 0.993 | 0.970 |
| lognormal k=3 n=100 d=0.0 | 0.013 | 0.000 | 0.000 | 0.010 | 0.080 | 0.007 |
| lognormal k=3 n=100 d=0.5 | 1.000 | 0.997 | 0.743 | 0.857 | 0.990 | 0.910 |
| lognormal k=3 n=100 d=0.8 | 1.000 | 1.000 | 0.990 | 1.000 | 1.000 | 1.000 |
| k=5 n=10 H0 | 0.103 | 0.000 | 0.000 | 0.003 | 0.067 | 0.007 |
| k=5 n=30 H0 | 0.043 | 0.000 | 0.000 | 0.003 | 0.037 | 0.003 |
| k=5 n=30 d=.8 | 0.900 | 0.537 | 0.480 | 0.593 | 0.883 | 0.583 |
| k=8 n=20 H0 | 0.050 | 0.000 | 0.000 | 0.007 | 0.047 | 0.003 |

IC binomial à 95 % sur N=300 : 0,05 → [0,03 ; 0,08] ; 0,01 → [0,002 ; 0,03].

**Lecture.**
- **Variances inégales, le critère décisif.** (5,5,40), sd (4,4,1), H0 : 1,3 % de P>0,95 contre 19 % en gaussien seul, 26 % pour aov à p<0,005, 0,3 % pour Welch à p<0,005. La voie reste plus exposée à (2,2,40) : 3 % à P>0,95, 54 % à P>0,5, avec l'avis « n < 3 ».
- **Puissance contre Welch à variances inégales.** (5,5,40), effet 4 ou 8 : au-dessus de Welch à seuil comparable. (30,30,30), sd (1,1,3), d=0,8 : 34 % à P>0,95, 6 % en gaussien, 47 % pour Welch à p<0,005 : **ici en dessous de Welch**.
- **Normale à variances égales, d=0,5 ou 0,8.** Le coût de la moyenne sur les lois est négligeable devant le gaussien seul.
- **t3 et lognormale.** Nette hausse de puissance sur le gaussien : n=30, k=2, d=0,8 : t3 86 % contre 61 %, lognormale 91 % contre 57 %. Faux positifs ≤ 1 % à P>0,95 à H0, sur tout le tableau hors (2,2,40), (10,30,50) et (2×6).
- **Non testé, à connaître.** Sous une loi t, le paramètre comparé est un centre, qui n'est pas la moyenne d'une loi asymétrique. Des groupes de formes différentes à moyennes égales (centres différents) n'ont pas été simulés : la réponse pourrait y annoncer un effet que les moyennes n'ont pas.

### Groupes de formes différentes, moyennes et variances égales (`agent_43`, N=200 par plan)

Groupe 1 lognormal, les autres normaux.

| Plan | m.bay P>0,5 | m.bay P>0,95 | Gaussien P>0,95 | Welch p<0,05 | Kruskal p<0,05 |
|---|---:|---:|---:|---:|---:|
| n=30, k=2 | 0,180 | 0,020 | 0,000 | 0,055 | 0,075 |
| n=100, k=2 | 0,455 | 0,095 | 0,000 | 0,035 | 0,140 |
| n=400, k=2 | 0,955 | **0,670** | 0,000 | 0,040 | 0,440 |
| n=30, k=3 | 0,160 | 0,005 | 0,000 | 0,050 | 0,100 |
| n=100, k=3 | 0,355 | 0,085 | 0,000 | 0,035 | 0,150 |

**Faux « Oui » si l'on s'intéresse aux moyennes**, vrai si l'on s'intéresse à une différence de centre : Kruskal réagit aussi. Le calcul n'est pas changé. Le bilan avertit désormais quand les queues lourdes dominent, et plus fort si la loi normale dit le contraire. Limite d'une méthode, à connaître, pas réglée.

### Variances étagées, sans effet (`agent_44`, N=100 par plan, 1500 tirages)

```
Sous H0, plans a ecarts-types etages : P>.5 P>.95 mediane | Welch p<.05 p<.005 | P(lois) moyennes | lois ecartees | err mediane
k=8 n=6 sd 1e-2 a 1e2 : 0.240 0.020 0.342 | 0.080 0.010 | 0 0.518 0 0.816 | 60 | 4.2
k=8 n=6 sd 1e-1 a 1e1 : 0.220 0.020 0.303 | 0.080 0.010 | 0 0.373 0 0.66 | 10 | 3.6
k=3 n=6 sd 1e-2 a 1e2 : 0.300 0.000 0.380 | 0.050 0.000 | 0 0.611 0 0.695 | 58 | 3.5
k=3 n=20 sd 1e-2 a 1e2 : 0.140 0.000 0.287 | 0.090 0.000 | 0 0.773 0 0.277 | 28 | 2.4
k=8 n=20 sd 1e-2 a 1e2 : 0.130 0.000 0.229 | 0.090 0.000 | 0 0.842 0 0.233 | 44 | 4.6
```

Colonnes de P(lois) : moyennes sur les jeux où la loi reste ; le relevé des lois écartées compte toutes les lois. **À ces variances et à 1500 tirages, t_het est écartée dans 5 à 44 % des jeux selon le plan, t_eq dans 5 à 19 %** : la comparaison y est incomplète et l'avis le dit. Sous H0, P>0,95 reste à 2 % au plus, voisin de Welch à p<0,005. Le jeu du relecteur (P=0,998) est dans cette queue : un cas défavorable, à cette résolution.

### Hyper-a priori sur r

Loi uniforme sur r/2, r, 2r, fixée d'avance : P>0,5 et P>0,95 inchangés à ±0,033 au plus sur les 66 scénarios ; P<0,05 peut bouger plus, jusqu'à 0,12 (lognormale, n=100, k=3, H0 : 0,483 contre 0,360). Colonnes « mix » de `agent_41_B.txt`. **Décision, D5 : r reste fixé à 1/2, la sensibilité montre r/2 et 2r.** Sur π : voir METHODOLOGIE §7. L'utilisateur n'a rien à régler dans les deux cas.

### Sensibilité aux a priori (`agent_41 C`, N=100 par cas)

| n | k | d | Tête change | Seuil franchi | Stable | P change de côté selon r | selon u |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 2 | 0.0 | 0.930 | 0.070 | 0.000 | 0.060 | 0.000 |
| 10 | 2 | 0.5 | 0.870 | 0.130 | 0.000 | 0.180 | 0.030 |
| 10 | 2 | 0.8 | 0.650 | 0.290 | 0.060 | 0.120 | 0.010 |
| 10 | 3 | 0.0 | 0.750 | 0.250 | 0.000 | 0.100 | 0.010 |
| 10 | 3 | 0.5 | 0.810 | 0.170 | 0.020 | 0.110 | 0.010 |
| 10 | 3 | 0.8 | 0.660 | 0.320 | 0.020 | 0.150 | 0.010 |
| 30 | 2 | 0.0 | 0.460 | 0.540 | 0.000 | 0.040 | 0.000 |
| 30 | 2 | 0.5 | 0.570 | 0.380 | 0.050 | 0.180 | 0.000 |
| 30 | 2 | 0.8 | 0.330 | 0.500 | 0.170 | 0.090 | 0.000 |
| 30 | 3 | 0.0 | 0.360 | 0.640 | 0.000 | 0.090 | 0.010 |
| 30 | 3 | 0.5 | 0.620 | 0.360 | 0.020 | 0.270 | 0.010 |
| 30 | 3 | 0.8 | 0.220 | 0.530 | 0.250 | 0.110 | 0.010 |
| 100 | 2 | 0.0 | 0.210 | 0.790 | 0.000 | 0.070 | 0.000 |
| 100 | 2 | 0.5 | 0.170 | 0.570 | 0.260 | 0.030 | 0.010 |
| 100 | 2 | 0.8 | 0.010 | 0.110 | 0.880 | 0.010 | 0.000 |
| 100 | 3 | 0.0 | 0.120 | 0.880 | 0.000 | 0.040 | 0.000 |
| 100 | 3 | 0.5 | 0.200 | 0.450 | 0.350 | 0.120 | 0.000 |
| 100 | 3 | 0.8 | 0.000 | 0.010 | 0.990 | 0.000 | 0.000 |

Avec la grille d'origine, cotes 1:4 à 4:1, l'hypothèse en tête changeait dans 98 % des jeux sans effet jusqu'à n=100 (premier lancement, chiffres remplacés). D'où la grille 1:2 à 2:1 et la ligne « il faudrait croire d'avance à l'autre à X contre 1 ». Sans effet, le verdict dit encore « ne tranchent pas » dans 93 % à n=10, 46 % à n=30, 21 % à n=100 : c'est vrai, H0 ne s'établit que lentement. Le choix de cette grille est un jugement, pas un résultat.

### Défauts de relecture (`agent_41 D`)

```
B1 : petit groupe eloigne, plus il s'eloigne, plus H1 doit gagner
n=(10,10,10,10,3) d=    5 : log BF     3.29, P(H1) 0.9642, P(variances inegales) 0.308, err 1.1%, zs 10.2, Welch p 0.00016
n=(10,10,10,10,3) d=   50 : log BF    12.24, P(H1) 1.0000, P(variances inegales) 0.255, err 1.1%, zs 85.6, Welch p 3.3e-15
n=(10,10,10,10,3) d=  500 : log BF    21.46, P(H1) 1.0000, P(variances inegales) 0.256, err 1.0%, zs 170.8, Welch p 1.8e-26
n=(10,10,10,10,3) d= 5000 : log BF    30.67, P(H1) 1.0000, P(variances inegales) 0.253, err 1.1%, zs 256.0, Welch p 1e-37
n=(30,30,2) d=    5 : log BF     2.38, P(H1) 0.9154, P(variances inegales) 0.622, err 1.0%, zs 12.1, Welch p 2.1e-30
n=(30,30,2) d=   50 : log BF     8.69, P(H1) 0.9998, P(variances inegales) 0.560, err 1.0%, zs 122.5, Welch p 1.2e-66
n=(30,30,2) d=  500 : log BF    15.56, P(H1) 1.0000, P(variances inegales) 0.559, err 1.0%, zs 255.8, Welch p 5.6e-103
n=(30,30,2) d= 5000 : log BF    22.47, P(H1) 1.0000, P(variances inegales) 0.561, err 1.1%, zs 389.3, Welch p 2.6e-139
n=(20,20,20,20,3) d=    5 : log BF     4.44, P(H1) 0.9883, P(variances inegales) 0.253, err 1.0%, zs 15.4, Welch p 4.5e-21
n=(20,20,20,20,3) d=   50 : log BF    13.51, P(H1) 1.0000, P(variances inegales) 0.245, err 1.1%, zs 163.7, Welch p 3.1e-48
n=(20,20,20,20,3) d=  500 : log BF    22.72, P(H1) 1.0000, P(variances inegales) 0.244, err 1.1%, zs 340.7, Welch p 1.6e-75
n=(20,20,20,20,3) d= 5000 : log BF    31.93, P(H1) 1.0000, P(variances inegales) 0.245, err 1.1%, zs 518.0, Welch p 8.8e-103
n=(10,10,10,10,10,10,10,3) d=    5 : log BF     3.50, P(H1) 0.9708, P(variances inegales) 0.225, err 1.4%, zs 10.5, Welch p 1.9e-06
n=(10,10,10,10,10,10,10,3) d=   50 : log BF    12.99, P(H1) 1.0000, P(variances inegales) 0.171, err 1.5%, zs 138.8, Welch p 3.4e-26
n=(10,10,10,10,10,10,10,3) d=  500 : log BF    22.26, P(H1) 1.0000, P(variances inegales) 0.170, err 1.4%, zs 286.6, Welch p 2.5e-46
n=(10,10,10,10,10,10,10,3) d= 5000 : log BF    31.47, P(H1) 1.0000, P(variances inegales) 0.169, err 1.5%, zs 434.0, Welch p 2.4e-66
n=(5,5,5,5,5,5,5,2) d=    5 : log BF     1.33, P(H1) 0.7913, P(variances inegales) 0.563, err 1.5%, zs 4.5, Welch p 0.013
n=(5,5,5,5,5,5,5,2) d=   50 : log BF     7.78, P(H1) 0.9996, P(variances inegales) NA, err 1.5%, zs 58.0, Welch p 8.2e-11
n=(5,5,5,5,5,5,5,2) d=  500 : log BF    14.70, P(H1) 1.0000, P(variances inegales) NA, err 1.5%, zs 122.4, Welch p 3.8e-20
n=(5,5,5,5,5,5,5,2) d= 5000 : log BF    21.46, P(H1) 1.0000, P(variances inegales) 0.334, err 1.3%, zs 186.8, Welch p 1.7e-29
n=(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,3) d=    5 : log BF     6.58, P(H1) 0.9986, P(variances inegales) NA, err 2.0%, zs 5.5, Welch p 0.0043
n=(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,3) d=   50 : log BF    15.78, P(H1) 1.0000, P(variances inegales) NA, err 2.2%, zs 286.9, Welch p 5e-39
n=(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,3) d=  500 : log BF    25.05, P(H1) 1.0000, P(variances inegales) NA, err 2.1%, zs 682.2, Welch p 2.5e-94
n=(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,3) d= 5000 : log BF    34.26, P(H1) 1.0000, P(variances inegales) NA, err 2.1%, zs 1078.3, Welch p 7e-150
B2 : un groupe disperse sans rapport ne doit pas changer la reponse sur les deux autres
deux groupes seuls : log BF 3.18, P(H1) 0.960
S=1 : log BF 5.55, P(H1) 0.996, err 0.2%
S=10 : log BF 2.63, P(H1) 0.933, err 0.9%
S=100 : log BF 2.42, P(H1) 0.918, err 0.9%
S=1000 : log BF 2.40, P(H1) 0.917, err 0.9%
S=10000 : log BF 2.40, P(H1) 0.917, err 0.9%
S=1e+06 : log BF 2.34, P(H1) 0.912, err 0.9%
c(0,1e-160,2e-160,1,1+1e-10,1) : log BF 68.942, P 1.0000, err 0.6%
donnees constantes a l'arrondi : y est constant : les probabilités sont indéterminées.
N = k+1, (1,1+1e-15,5) : log BF 3.029, P 0.9539, err 1.4%
N = k+1, (1,1,5) : Aucune variation dans aucun groupe : les lois normale et t ne s'appliquent pas.
une valeur aberrante, n=20, k=2 : log BF -1.683, P 0.1567, err 1.2%
groupe de n=1 : log BF 0.942, P 0.7196, err 1.0%
(2,2) : log BF 1.274, P 0.7814, err 0.8%
n=1e4 par groupe, k=2 : log BF -3.734, P 0.0233, err 0.0%
```

Pentes conformes au calcul à la main : d^(n_s + 2u − 2) une fois la seconde forme en jeu. Avant correction, k=20, n_s=3, d=10 donnait P(H1) = 0,17 (RELECTURE, B1). B2 : plancher à log BF ≈ 2,4, contre 3,18 pour les deux groupes seuls.

### Tests

Sur le code final : `test-bay-core` 237, `test-bay-oneway` 27, `test-m.bay` 57, `test-bay-laws` 77, `test-m.test` 52, `test-valreg` 12, `test-manova` 59 assertions, 0 échec, 0 erreur. Versions de code par lancement : JOURNAL.
