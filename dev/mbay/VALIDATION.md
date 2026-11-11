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

| t et normale à variance par groupe, k=3, H0, première forme, seconde forme pour chaque groupe s, et leur mélange, 3 jeux (graines 4501 à 4503) | log m : écart moyen ≤ 0,021, extrêmes -0,075 et +0,041 sur 3 graines KefiR à 4 000 tirages ; log BF mélange contre H0 : ≤ 0,061 (6 % sur le BF), sans biais de signe | référence `agent_45`, comparaison `agent_48` |

**Référence `agent_45`** : écrite par un agent qui n'avait pas lu `sys_bay_laws.R`, sur la seule spécification des a priori. Monte-Carlo naïf dans l'a priori (ϑ, ν, g et, pour t, les échelles λ_ij du mélange), μ, τ et l'effet intégrés exactement ; 10^8 tirages en t, 10^7 en normale ; erreur standard ≤ 0,0008. Contrôlée par quadratures déterministes (normale à variance par groupe, 5 formes ; t à variance par groupe sous H0, par une factorisation propre) à 1,5 erreur standard près, et par un Monte-Carlo entièrement direct, grossier. Pour t à variance par groupe sous H1, la référence est ce seul Monte-Carlo, de conditionnelle vérifiée. La référence porte sur y brut ; KefiR standardise y par l'écart-type intra-groupe poolé s_p, d'où le décalage (N−1) log s_p, que l'accord à 0,01 près confirme.

**Non comparé** : k ≥ 4 à variance par groupe ; la loi t à variance par groupe sous H1 n'a pas de quadrature déterministe.

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

### Seuil de l'avis « queues lourdes » (`agent_46`, N=300 par plan, 1500 tirages)

L'avis s'affiche quand P(loi t | données) ≥ 0,95.

| Données, k=2 et 3 | n=10 | n=30 | n=100 |
|---|---:|---:|---:|
| Normales, sans effet : fausse alerte | 0,3 à 0,7 % | 1,0 à 1,3 % | 0 à 0,3 % |
| Normales, d=0,5 : fausse alerte | 0 % | 0 à 0,3 % | 0 % |
| t3 : détection | 12 à 25 % | 52 à 75 % | 97 à 100 % |
| Lognormale contre normale : détection | 7 à 9 % | 29 à 43 % | 72 à 81 % |

À 0,9, la fausse alerte reste ≤ 1,3 % et la détection sur t3 gagne 4 à 9 points ; à 0,99, la fausse alerte tombe à 0 et la détection perd 5 à 15 points. 0,95 est retenu comme le seuil du « très probablement » du bilan. Contradiction (m.bay > 0,75 alors que la meilleure loi normale donne < 0,5, ou l'inverse) : jamais sur données normales, 0 à 9 % sur t3, jusqu'à 22 % sur lognormale à n=100, où la comparaison des centres diffère de celle des moyennes ; le bilan l'annonce alors en première ligne et plafonne le soutien.

### Tests

Sur le code final : `test-bay-core` 237, `test-bay-oneway` 27, `test-m.bay` 57, `test-bay-laws` 77, `test-m.test` 52, `test-valreg` 12, `test-manova` 59 assertions, 0 échec, 0 erreur. Versions de code par lancement : JOURNAL.

## 8. `m.bay()` à plusieurs facteurs, loi normale à variance commune, via BayesFactor (9 novembre)

Vraisemblances : `BayesFactor::lmBF` 0.9.12-4.8 (D1). KefiR : espace héréditaire, loi HOP (D2), probabilité par terme, sensibilité, bilan (`R/sys_bay_factorial.R`).

### Calcul indépendant (`agent_51`, agent qui n'a lu ni `R/`, ni `tests/`, ni ces documents)

Trois plans 2×2 (5 par case ; 3, 7, 4, 9 ; 10 par case), modèles A, B, A+B, A+B+A:B contre le modèle constant. Référence : forme fermée à g fixés, puis trapèzes en log g (1D, 2D, 3D), convergence ≤ 7e-10, recoupée par Monte-Carlo naïf dans l'a priori des g.

| Modèles | Écart `lmBF` − référence sur le log BF |
|---|---|
| À un terme (6) | < 1e-6 |
| A+B (3) | 0,0008 à 0,0015 |
| A+B+A:B (3) | 0,0001 à 0,023 |

Moyenne sur 10 graines : 5 cas sur 6 à moins d'une erreur standard, un à −2,5 (jeu déséquilibré, A+B, −0,006 sur le log BF), compatible avec le hasard sur 6 comparaisons. Codage lu dans le paquet installé : base orthonormée du complément de 1 (`fixedFromRandomProjection`), produits pour les interactions, r « medium » = 0,5 par terme.

**L'erreur annoncée par BayesFactor sous-estime son bruit** sur les modèles à plusieurs termes : écart-type entre graines 0,014 à 0,026 sur le log BF, contre 0,003 à 0,010 annoncés à 1e5 itérations. Recontrôlé par moi : à 1e4, 1e5 itérations, écart-type 0,022 et 0,020 (il ne baisse presque pas), erreur annoncée 0,018 et 0,013. KefiR garde 10 000 itérations et annonce au moins 3 % d'erreur dès qu'un modèle a plusieurs termes. Ce bruit (2 % sur le BF) ne change aucune conclusion à 3 décimales de probabilité près, sauf au voisinage immédiat d'un seuil.

**Non vérifié** : facteurs à plus de 2 niveaux (base de dimension ≥ 2), plus de deux facteurs, r ≠ 0,5, code C++ de l'a priori multi-termes.

### Concordance avec la voie interne à un facteur

`lmBF(y ~ g)` et la première forme de `normal_eq` : écart < 1e-4 à k = 2, 3, 5 et r = 0,25, 0,5, 1 (`test-bay-factorial.R`).

### Effets fixes, plans déséquilibrés, à côté d'aov de type II (`agent_50 B`)

Cases : 2×2 (4, 8, 12, 16) ; 2×3 (3, 9, 5, 7, 10, 4) ; 2×2×2 (3 à 8). Effets ±f (−f, 0, f pour 3 niveaux), σ = 1. N = 300 par cas (2×2×2 : 150). m.bay : loi HOP, r = 0,5. « Faux » = terme sans effet réel ; aov : `car::Anova(type = 2)`, contrastes somme nulle.

| Cas | m.bay : un faux terme P>0,5 / P>0,95 | aov : un faux p<0,05 | m.bay : vrai P>0,5 | aov : vrai p<0,05 |
|---|---:|---:|---:|---:|
| 2×2 sous H0 | 0,083 / 0,003 | 0,127 | – | – |
| 2×2, A et B, f = 0,25 | 0,017 / 0 | 0,070 | 0,355 | 0,302 |
| 2×2, A et B, f = 0,5 | 0,073 / 0 | 0,070 | 0,907 | 0,808 |
| 2×2, A:B seule, f = 0,5 | aucun faux terme (A et B agissent par A:B ; « oui » à A ou B : 0,653 / 0,180, voir ci-dessous) | 0,337 | 0,427 | 0,763 |
| 2×3 sous H0 | 0,097 / 0,003 | 0,123 | – | – |
| 2×3, A et B, f = 0,5 | 0,017 / 0,003 | 0,047 | 0,610 | 0,640 |
| 2×3, A, B et A:B, f = 0,5 | 0 / 0 | 0 | 0,656 | 0,623 |
| 2×2×2 sous H0 | 0,120 / 0 | 0,313 | – | – |
| 2×2×2, A, B, C, f = 0,5 | 0,140 / 0 | 0,200 | 0,933 | 0,887 |

**Lecture.** Sous H0, moins de « oui » à P>0,5 que de p<0,05 sur au moins un terme (8 à 12 % contre 12 à 31 %), et presque jamais P>0,95. Effets additifs : puissance comparable à aov, fausse interaction 1,7 à 7 %. **Interaction seule** : A et B y sont « oui » dans 65 % des jeux. C'est voulu par la question posée (« y diffère selon A, seul ou en interaction » : A joue bien, par l'interaction), et l'hérédité l'impose ; ce n'est pas une erreur de type I au sens d'aov, qui teste l'effet principal moyen. Le prix est ailleurs : **l'interaction réelle est moins souvent détectée** (P(A:B)>0,5 dans 43 % contre 76 % de p<0,05), parce qu'elle part de 1/6 a priori ; les seuils ne sont pas équivalents. À annoncer dans la documentation.

**Non simulé** : facteurs à plus de 3 niveaux, 4 facteurs, plans très déséquilibrés (rapport > 5 entre cases).

### Plans 2³ équilibrés, interaction réelle ou d'ordre 3 (`agent_49`, fin du 9 novembre)

Résultats complets dans `scripts/resultats/agent_49.txt` ; tableau et lecture dans DECISIONS, D2.
- Une interaction double réelle (A:B) fait monter les fausses interactions doubles sous HOP : 20,5 % (n=5) et 27 % (n=10), contre 13 % pour aov et 11 à 13 % pour 1/2 par terme.
- Avec toutes les interactions doubles réelles, A:B:C fausse sort dans 3 à 6 % des jeux.
- A:B:C seule, sans effet inférieur, est hors de l'espace héréditaire : puissance 0,33 contre 0,995 pour aov. Les seuls faux termes y sont A:B, A:C et B:C, puisque A, B et C agissent par A:B:C.
- Au-delà du qualificatif « faible » (P>0,75), les fausses interactions de la HOP tombent à 5 à 6,5 % (`agent_49b`, relecture des mêmes BF, `resultats/agent_49b.txt`).

### Calibration sous l'a priori (`agent_50 A`, N = 1 200 jeux par plan)

Modèle tiré selon la loi HOP, g ~ InvGamma(1/2, r²/2) par terme, effets sur le codage de BayesFactor (`fullDesignMatrix`), σ = 1. Probabilité annoncée par terme contre fréquence du terme vrai.

| Plan | Ensemble : annoncé / vrai | z | Aucun effet : annoncé / vrai |
|---|---|---:|---|
| 2×2, 5 par case | 0,384 / 0,390 | 0,87 | 0,340 / 0,323 |
| 2×3 déséquilibré | 0,396 / 0,399 | 0,55 | 0,324 / 0,323 |
| 2×2 déséquilibré | 0,392 / 0,399 | 1,20 | 0,328 / 0,323 |

Par tranche (0 à 0,1 ; 0,1 à 0,3 ; … ; 0,9 à 1), l'écart dépasse 0,03 seulement dans la tranche 0,5 à 0,7 : annoncé 0,585 à 0,589, observé 0,623 à 0,632, dans les trois plans (environ 230 termes par plan, soit un peu plus d'une erreur standard chacun). Légère prudence possible à ce niveau, sans effet sur les seuils de 0,75 et 0,95. À 600 jeux, le 2×2 déséquilibré donnait z = 2,48 : non confirmé à 1 200. Un jeu tiré avec g = 2·10⁴ faisait échouer l'`integrate` de `lmBF` (NA) : les modèles à un seul effet principal passent désormais par le calcul interne, et le jeu a été recalculé.

## 9. Post-hoc à un facteur (D11), 9 novembre

### Calcul indépendant (`agent_53`, agent qui n'a lu ni `R/`, ni `tests/`, ni `dev/mbay`)

Jeu : 4 groupes de 12, `set.seed(4)`, moyennes vraies 10 ; 11,2 ; 10,1 ; 10. L'agent a programmé la spécification de D11 : 15 partitions, BF de `lmBF`, a priori nul à 1/2 et Ewens(1) sur le reste, puis 20 000 tirages a posteriori et la HDI de bayestestR. Scripts et sorties : `relecture/agent_53/`.

| Paire | P(différents), référence | m.bay | Écart, référence (HDI) | m.bay (4 000 tirages) |
|---|---:|---:|---|---|
| L1-L2 | 0,9773 | 0,977 | 0,977 [0,287 ; 1,696] | 0,99 [0,28 ; 1,7] |
| L1-L3 | 0,2978 | 0,298 | −0,204 [−0,877 ; 0,457] | −0,22 [−0,88 ; 0,43] |
| L1-L4 | 0,4010 | 0,401 | −0,534 [−1,200 ; 0,156] | −0,55 [−1,2 ; 0,12] |
| L2-L3 | 0,9903 | 0,990 | −1,181 [−1,888 ; −0,429] | −1,2 [−1,9 ; −0,44] |
| L2-L4 | 0,9947 | 0,995 | −1,518 [−2,276 ; −0,769] | −1,5 [−2,3 ; −0,81] |
| L3-L4 | 0,3046 | 0,305 | −0,330 [−0,997 ; 0,346] | −0,33 [−1 ; 0,3] |

Regroupements de tête identiques : L1 = L3 = L4 < L2 (0,532), L1 = L3 | L2 | L4 (0,163), L1 | L2 | L3 = L4 (0,149). Les probabilités concordent à l'affichage près. Les écarts et les HDI concordent à la précision du tirage (4 000 contre 20 000 tirages), au signe de la paire près (m.bay donne groupe2 − groupe1). La cohérence P(i≠k) ≤ P(i≠j) + P(j≠k) tient sur les 24 triplets. **Limite** : les BF des partitions viennent de BayesFactor des deux côtés ; leur exactitude repose sur VALIDATION §8 (`agent_51`, calcul JZS indépendant à un terme).

**Dans `test-bay-posthoc.R`** :
- recalcul par `lmBF` direct à k = 3 ;
- nombres de Bell, a priori de somme 1, rapport des poids d'Ewens ;
- δ contre le d de Cohen à n = 400 (écart < 0,03), largeur de l'intervalle contre l'erreur standard asymptotique, borne contre l'intervalle de Student.

### Simulation à côté de Tukey HSD (`agent_52`, N = 300 jeux par scénario)

Moyennes vraies en écarts-types, σ = 1, r = 0,5, `prior_h1` = 1/2, Ewens α = 1. « Faux » : au moins une paire vraiment égale dite différente. « Puissance » : moyenne sur les paires vraiment différentes. Résultats : `scripts/resultats/agent_52.txt`.

| Scénario | Faux P>0,5 | Faux P>0,75 | Faux P>0,95 | Tukey faux p<0,05 | Puissance P>0,5 | Puissance P>0,75 | Tukey p<0,05 | Couverture HDI δ | Regroupement vrai en tête |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| k=3, n=10, H0 | 0,107 | 0,013 | 0 | 0,063 | – | – | – | – | 0,953 |
| k=3, n=20, H0 | 0,073 | 0,010 | 0,003 | 0,050 | – | – | – | – | 0,967 |
| k=4, n=10, H0 | 0,073 | 0,020 | 0,003 | 0,060 | – | – | – | – | 0,967 |
| k=5, n=10, H0 | 0,063 | 0,027 | 0,007 | 0,043 | – | – | – | – | 0,973 |
| k=3, n=15, un groupe à +1 | 0,110 | 0,033 | 0,007 | 0,023 | 0,758 | 0,588 | 0,625 | 0,903 | 0,630 |
| k=4, n=15, un groupe à +1 | 0,190 | 0,050 | 0,010 | 0,033 | 0,787 | 0,586 | 0,556 | 0,921 | 0,630 |
| k=4, n=15, deux blocs (0, 0, 0,8, 0,8) | 0,250 | 0,093 | 0,017 | 0,040 | 0,583 | 0,377 | 0,331 | 0,916 | 0,330 |
| k=3, n=15, trois niveaux (0 ; 0,6 ; 1,2) | – | – | – | – | 0,647 | 0,402 | 0,438 | 0,950 | 0,170 |

**Lecture.**
- Sous H0, une paire est dite différente au-delà de « faible » (P > 0,75) dans 1 à 3 % des jeux, contre 4 à 6 % de p < 0,05 chez Tukey. Le regroupement vrai (aucune différence) est en tête dans 95 à 97 % des jeux.
- **Avec des différences réelles, les paires vraiment égales sortent plus souvent qu'avec Tukey.** Au-delà de « faible », c'est 3 à 9 % des jeux contre 2 à 4 % ; à P > 0,5, de 11 à 25 %. Dans le bilan, une paire entre 0,5 et 0,6 est « non tranchée », car son soutien est nul. La puissance au-delà de « faible » est comparable à Tukey (0,38 à 0,59 contre 0,33 à 0,63). Les seuils ne sont pas équivalents.
- **La HDI de δ couvre la vraie valeur dans 90 à 95 % des cas**, sous le nominal quand seuls certains groupes diffèrent. Elle vient du modèle « tous les groupes peuvent différer », qui rapproche les écarts de zéro. C'est un intervalle de crédibilité, pas un intervalle de confiance. Le récit dit qu'il est « légèrement ramené vers zéro par l'a priori ».

### Audit de la session parente (9 novembre)

Données propres : 4 groupes de 12, moyennes vraies 0 ; 0,1 ; 1 ; 1,1.
- Cohérence avec TukeyHSD : a/c et b/c à 98,3 % (p = 0,0006), d/c à 80,8 % (p = 0,12), a/b à 24,4 %.
- Sous H0 à k = 4 (N = 40), la paire maximale dépasse 0,75 dans 5 % des jeux et 0,95 dans 0 %.
- Durée : 13 s.
- Quatre points demandés et traités : écart observé et écart estimé nommés, k > 6, couverture des intervalles, récit à variances inégales.

### Couverture des intervalles, toutes les paires (`agent_54`, N = 300 jeux par scénario)

Les intervalles viennent du modèle « tous les groupes peuvent différer ». Le tableau donne la part des jeux où l'intervalle à 95 % contient la vraie valeur, pour l'écart brut et pour δ. Résultats : `scripts/resultats/agent_54.txt`.

| Scénario | Paires égales : écart / δ (nb) | Paires différentes : écart / δ (nb) |
|---|---|---|
| k=2, n=10, δ = 0,5 | – | 0,947 / 0,950 (1) |
| k=3, n=10, H0 | 0,993 / 0,992 (3) | – |
| k=3, n=15, un groupe à +1 | 0,987 / 0,993 (1) | 0,922 / 0,913 (2) |
| k=4, n=20, deux blocs (0, 0, 0,8, 0,8) | 0,970 / 0,968 (2) | 0,925 / 0,927 (4) |
| k=4, n=10, graduel (0 ; 0,3 ; 0,6 ; 0,9) | – | 0,942 / 0,955 (6) |

**Lecture.**
- Le rétrécissement vers zéro sur-couvre les paires vraiment égales : 0,97 à 0,99.
- Il sous-couvre un peu les paires vraiment différentes quand d'autres sont égales : 0,91 à 0,93.
- Quand toutes les paires diffèrent (k = 2, graduel), la couverture est nominale : 0,94 à 0,96.

### Au-delà de 8 groupes : tirage de Gibbs

- Accord avec l'énumération exacte, en forçant le tirage : 1,4 point au plus à k = 6 et 7. Avec 3 000 tours, les 203 et les 877 partitions sont toutes visitées.
- Dans `test-bay-posthoc.R`, le seuil est abaissé à 3 par l'option `KefiR.bay_exact_max`. À k = 4, l'écart au calcul exact reste sous 3 points sur `p_differ` et sous 4 points sur `p_min` et `p_max` (repondération en α).
- Aussi testés : la constante d'Ewens contre la somme directe à k = 5, et les nombres de Bell jusqu'à 9.
- Récit à 9 groupes : `relecture/lecteur_recit/recits_posthoc.txt`, cas C. L'erreur de tirage annoncée est d'environ 4 points à 2 000 tours.

### Récit à variances inégales

Cas A de `relecture/lecteur_recit/recits_posthoc.txt` : 3 laboratoires, écarts-types 1, 1 et 3.
- L'étape 4 annonce des variances inégales (99,5 %).
- L'étape 6 donne les paires, puis « les lire avec réserve ».
- La conclusion le répète sous les paires.

C'est testé dans `test-bay-posthoc.R`.

**Non simulé** :
- variances inégales et queues lourdes : le post-hoc suppose la loi normale à variance commune, et le récit prévient alors ;
- α = 1/2 et 2, calculés en sensibilité (`p_min`, `p_max`), sans simulation ;
- plans déséquilibrés ;
- k de 6 à 8 : le calcul exact est vérifié, pas sa calibration ;
- au-delà de 8 groupes : la calibration du tirage n'est pas simulée.

## 10. Données appariées et mesures répétées (D15), 9 novembre

### Calcul indépendant (`agent_55`, agent qui n'a lu ni `R/`, ni `tests/`, ni `dev/mbay`)

Scripts et sorties : `relecture/agent_55/`. Mêmes jeux que les récits D, E et F (`relecture/lecteur_recit/recits_repete.txt`).

| Calcul | Référence | m.bay |
|---|---|---|
| A. Apparié, 13 paires, log BF10 (intégrale de Rouder et al. 2009 programmée par l'agent ; `ttestBF` à 1,7e-10 près) | 0,5474 | 0,548 (BF 1,73) |
| B1. 4 sessions, BF « session + sujet » contre « sujet seul » (intégrale en 2 dimensions sur g_fixe et g_sujet, programmée par l'agent ; `anovaBF` : 23,515) | 23,57 (log 3,1602) | 23,5 |
| B2. Paires S1-S2 / S1-S3 / S1-S4 / S2-S3 / S2-S4 / S3-S4 | 0,439 / 0,661 / 0,916 / 0,490 / 0,846 / 0,684 | 0,440 / 0,661 / 0,916 / 0,491 / 0,847 / 0,686 |
| B2. Trois regroupements de tête | 0,221 / 0,154 / 0,154 | 0,221 / 0,154 / 0,154 |
| C. Plan en tranches, log BF contre « id seul » : groupe / séance / groupe + séance / groupe * séance (`anovaBF`, 50 000 itérations) | −0,340 / 0,270 / −0,018 / 1,245 | −0,336 / 0,254 / −0,065 / 1,266 |

**B2, version corrigée.** Dans `ref.R`, la première version de l'agent mesurait les regroupements non nuls contre l'intercept seul, alors que « aucune différence » restait à 0. Sa propre spécification demandait « contre sujet seul ». Ses paires valaient donc les nôtres divisées par 1 − P(aucune différence) = 0,924, par exemple 0,991 au lieu de 0,916.
- Le calcul corrigé, avec la seule division ajoutée, est dans `ref_B2_corrige.R` et `sortie_B2_corrige.txt`.
- Il concorde à 0,3 point près.
- Le fichier de l'agent est laissé tel quel.

**Limite.** B2 et C reposent sur BayesFactor des deux côtés : ce sont des estimations par tirage, avec des erreurs de ±0,4 à ±2,5 %. Seuls A et B1 ont un calcul programmé indépendamment de BayesFactor.

### Calibration (`agent_56`, N = 200 jeux par scénario)

Écart-type sujet 2, écart-type résiduel 1, r = 0,5, `prior_h1` = 1/2. Le classique est le t apparié pour k = 2, et `aov(y ~ temps + Error(id/temps))` pour k ≥ 3. Résultats : `scripts/resultats/agent_56.txt`.

| Scénario | « Différence » P>0,5 / P>0,75 / P>0,95 | Classique p<0,05 | Couverture de δ (dz si k = 2) | Paire égale dite différente (P>0,75) | Regroupement vrai en tête |
|---|---|---:|---:|---:|---:|
| k=2, n=12, H0 | 0,090 / 0,030 / 0 | 0,040 | 0,970 | – | – |
| k=2, n=20, H0 | 0,090 / 0,005 / 0 | 0,040 | 0,990 | – | – |
| k=2, n=15, écart 0,7 (dz ≈ 0,49) | 0,635 / 0,295 / 0,070 | 0,415 | 0,935 | – | – |
| k=3, n=12, H0 | 0,055 / 0,030 / 0 | 0,040 | 0,983 | 0,025 | 0,960 |
| k=3, n=15, une condition à +0,8 | 0,655 / 0,385 / 0,095 | 0,555 | 0,952 | 0,015 | 0,385 |
| k=4, n=12, graduel (0 ; 0,3 ; 0,6 ; 0,9) | 0,520 / 0,275 / 0,065 | 0,440 | 0,944 | – | 0 |

**Lecture.**
- **Sous H0**, une différence est annoncée au-delà de « faible » (P > 0,75) dans 0,5 à 3 % des jeux, et au-delà de 0,95 jamais, contre 4 % de p < 0,05.
- **Avec un effet**, la méthode est plus prudente que le test classique : 30 à 39 % au-delà de 0,75, contre 42 à 56 % de p < 0,05. Le facteur de Bayes demande plus de données pour un même effet.
- **Les intervalles** couvrent la vraie valeur dans 94 à 99 % des cas. La sur-couverture sous H0 vient du rétrécissement vers zéro.
- **Le regroupement vrai** sort en tête dans 96 % des jeux sous H0. Il le fait rarement quand l'effet est graduel : à n = 12, les 4 niveaux distincts ne se séparent pas.

**Non simulé** :
- plans en tranches, et plusieurs facteurs intra ;
- sujets incomplets en mesures répétées ;
- écart à la symétrie composée, que le modèle suppose et que le récit nomme ;
- queues lourdes.

### Tests

`test-bay-repeated.R`, 7 tests :
- l'apparié égale `ttestBF(paired)` à 1e-8 ;
- les quatre écritures de l'apparié concordent ;
- une paire incomplète est retirée et dite, et l'ordre des lignes est sans effet ;
- les mesures répétées concordent avec `lmBF(+ id) / lmBF(id)` ;
- l'effet sujet rend la réponse plus nette qu'avec des groupes indépendants ;
- le plan en tranches concorde avec `anovaBF(whichRandom)` à 0,1 près en log ;
- un plan avec répliques est refusé comme modèle mixte.

`test-bay-posthoc.R` : à l'échelle 1e200, les intervalles restent finis. C'est le cas qui produisait le message « ci is too small » de bayestestR.

## 11. Analyse de covariance (D16), 9 novembre

### Calcul indépendant (`agent_57`, agent qui n'a lu ni `R/`, ni `tests/`, ni `dev/mbay`)

Mêmes jeux que les récits G, H et I (`relecture/lecteur_recit/recits_ancova.txt`). L'agent a appliqué la spécification de D16 :
- covariable centrée ;
- `lmBF` contre « covariable seule » ;
- a priori 1/2, 1/4, 1/4 sur les trois modèles quand l'interaction est écrite ;
- `posterior()` à 20 000 tirages et HDI.

Il l'a comparée à `car::Anova` type II et à `emmeans`. Scripts et sortie : `relecture/agent_57/`.

| Calcul | Référence | m.bay |
|---|---|---|
| G. log BF engrais / P(engrais) | 2,2151 / 0,9016 | 2,2181 / 0,9019 |
| G. Écarts ajustés E2−E1, E3−E1, E3−E2 (médiane) | 0,434 / 2,060 / 1,616 | 0,424 / 2,045 / 1,606 |
| G. HDI de E3−E1 | [0,456 ; 3,588] | [0,434 ; 3,621] |
| H. P(méthode) / P(méthode:note) | 0,999428 / 0,999321 | 0,999417 / 0,999308 |
| H. Écart projet − classique aux quartiles 11,36 / 13,58 / 15,59 | −0,981 / 0,495 / 1,836 | −0,995 / 0,495 / 1,829 |
| H. Zone où la HDI exclut 0 | [6,754 ; 11,630] et [13,846 ; 17,836] | [6,754 ; 11,630] et [13,957 ; 17,836] |
| I. log BF régime | 2,2419 | 2,2418 |
| I. Écart ajusté R2 − R1 | 4,318 [0,533 ; 8,286] | 4,482 [0,620 ; 8,203] (4 000 tirages) |
| I. Part de chaque groupe dans l'étendue commune | 0,067 et 0,067 | 0,067, plus faible part (R1) |

**Lecture.**
- **H, zone.** La zone bayésienne diffère d'un pas de grille sur une borne (0,11). Les bornes exactes de Johnson-Neyman fréquentiste sont 11,47 et 13,82, ce qui confirme le sens et la place des deux zones.
- **H, effets.** `car::Anova` donne méthode p = 0,30 et méthode:note p = 1,6e−6. Le « oui » à méthode de m.bay se lit « seul ou en interaction » : il vient des pentes, ce que le récit dit en changeant de stratégie.
- **I, écart.** À 20 000 tirages, m.bay donne 4,36 à 4,42 selon la graine, identique à un appel direct de `posterior()` avec la même graine. L'écart avec la référence est donc du tirage, plus grand quand les groupes se recouvrent peu, car l'écart y est mal déterminé. `emmeans` donne 4,575 [1,07 ; 8,08].

### Calibration (`agent_58`, N = 200 jeux par scénario)

3 groupes, covariable U(5, 15), pente 0,8, résidu 1. Résultats : `scripts/resultats/agent_58.txt`.

| Scénario | P(g) > 0,5 / > 0,75 / > 0,95 | car type II p < 0,05 | P(g:x) > 0,5 / > 0,75 | p(g:x) < 0,05 | Couverture de l'écart ajusté |
|---|---|---:|---|---:|---:|
| Pentes parallèles, H0, n = 15 | 0,055 / 0,010 / 0 | 0,025 | – | – | 0,990 |
| Pentes parallèles, un groupe à +1, n = 15 | 0,835 / 0,585 / 0,300 | 0,755 | – | – | 0,932 |
| `y ~ g * x` écrit, pentes égales, un groupe à +0,5, n = 20 | 0,235 / 0,120 / 0,025 | 0,345 | 0 / 0 | 0,070 | 0,972 |
| `y ~ g * x` écrit, une pente +0,4, un groupe à +0,5, n = 20 | 0,705 / 0,565 / 0,350 | 0,315 | 0,630 / 0,460 | 0,935 | 0,955 |

**Lecture.**
- **Sous H0**, une différence est annoncée au-delà de « faible » dans 1 % des jeux.
- **Avec un écart**, la méthode est plus prudente que le test de type II : 59 % au-delà de 0,75, contre 76 % de p < 0,05.
- **Pentes égales**, quand l'interaction est écrite : jamais de « pentes différentes », contre 7 % de p < 0,05.
- **Pente réellement différente** : elle n'est dite probable que dans 63 % des jeux, contre 94 % pour le test F, et au-delà de 0,75 dans 46 %. L'interaction part à 25 % (HOP) et son a priori de pente est petit (`rscaleCont` = √2/4, partagé avec la covariable). C'est un manque de puissance à n = 20. Quand l'interaction n'est pas reconnue, la réponse reste celle des pentes parallèles, et le récit dit pourquoi.
- **P(g)** dépasse le type II dans ce dernier scénario, car il se lit « seul ou en interaction ».

**Non simulé** :
- plusieurs facteurs avec covariables ;
- plusieurs covariables ;
- groupes qui se recouvrent peu ;
- non-linéarité.

### Tests

`test-bay-ancova.R`, 6 tests :
- le BF du facteur égale `lmBF(y ~ g + xc) / lmBF(y ~ xc)` à 0,05 près ;
- le centrage est vérifié : résultat identique quand la covariable est décalée de 1 000 ;
- les écarts ajustés concordent avec `emmeans` à 0,1 près, et le sens avec `car::Anova` type II ;
- avec des pentes différentes, les écarts aux quartiles concordent avec `emmeans(at = )` à 0,15 près, et la zone est trouvée ;
- un faible recouvrement donne une prudence ajoutée et un soutien plafonné ;
- la covariable seule et l'ANCOVA appariée sont refusées.

## 12. Réponse binaire ou catégorielle à un facteur (D17, voie 1), 9 novembre

### Concordance avec BayesFactor

Dans `test-bay-contingency.R` :
- **Forme close et package.** La forme close de `.bay_cont_lbf` égale `contingencyTableBF` à 1e-6 sur 4 tableaux (2×2 à 4×3), pour a = 1, 2 et 3.
  - L'a priori de la fonction a été retrouvé : Dirichlet(R(a−1)+1) sous H0. Il explique pourquoi elle refuse a < 1.
- **Regroupements.** Les probabilités de paires à k = 3 sont recalculées à la main, chaque regroupement passant par `contingencyTableBF` sur le tableau fusionné.
- **Formats de la réponse.** 0/1, facteur, texte et TRUE/FALSE donnent le même BF.
- **Tailles d'effet, à n = 2 000 par groupe.**
  - Différence de proportions : à 0,005 de l'observé.
  - Log du rapport de cotes : à 0,02 de `glm` binomial.

### Calibration (`agent_59`, N = 500 jeux par scénario)

Comparaison avec le khi-deux, en p simulé quand un effectif attendu est inférieur à 5. Résultats : `scripts/resultats/agent_59.txt`.

| Scénario | « Différence » P>0,5 / P>0,75 / P>0,95 | Khi-deux p<0,05 | Couverture de la différence de proportions |
|---|---|---:|---:|
| k=2, n=20, H0 (p = 0,3) | 0,114 / 0,024 / 0 | 0,044 | 0,952 |
| k=2, n=50, H0 (p = 0,5) | 0,074 / 0,016 / 0,002 | 0,042 | 0,960 |
| k=3, n=30, H0 (p = 0,2) | 0,052 / 0,014 / 0,002 | 0,058 | 0,946 |
| k=2, n=10, H0 (p = 0,2) | 0,166 / 0,026 / 0 | 0,004 | 0,966 |
| k=2, n=40, 0,3 contre 0,5 | 0,494 / 0,300 / 0,106 | 0,378 | 0,966 |
| k=3, n=30, un groupe à 0,8 | 0,776 / 0,610 / 0,370 | 0,700 | 0,945 |

**Lecture.**
- Sous H0, une différence est annoncée au-delà de « faible » dans 1,4 à 2,6 % des jeux, et au-delà de 0,95 dans 0,2 % au plus.
- Avec un écart, la méthode est un peu plus prudente que le khi-deux : 30 à 61 % au-delà de 0,75, contre 38 à 70 %.
- Les intervalles couvrent la vraie valeur dans 95 à 97 % des cas.
- À n = 10, le khi-deux simulé ne rejette presque jamais. m.bay dépasse 0,5 dans 17 % des jeux, mais jamais 0,95, et le bilan dit alors « non tranché ».

## 13. Comptages à un facteur (D17, voie 2), 10 novembre

### Concordance

`test-bay-counts.R` :
- **Probabilités.** P(différence) et P(binomiale négative) égalent le calcul à la main sur les −BIC/2 des quatre modèles, à 1e-8.
- **Rapports des moyennes.** Ils égalent exp(coef) de `glm.nb` à 1e-6.
- **Contre-vérification BAS.** Sur 3 jeux de Poisson, le signe du log-BF est le même qu'avec `bas.glm` (Poisson, a priori robuste).
- **Surdispersion.** Elle est trouvée (binomiale négative > 90 %) et dite en changement de stratégie.
- **Refus.** Ils sont vérifiés, et les lignes tiennent dans 80 colonnes.

### Calibration (`agent_60`, N = 300 jeux par scénario, 3 groupes)

Résultats : `scripts/resultats/agent_60.txt`.

| Scénario | m.bay P>0,5 / P>0,75 / P>0,95 | Poisson seul P>0,5 / P>0,75 / P>0,95 | glm.nb, rapport de vraisemblance p<0,05 | Binomiale négative en tête |
|---|---|---|---:|---:|
| H0, Poisson, µ = 5, n = 20 | 0,010 / 0 / 0 | 0,013 / 0 / 0 | 0,067 | 0,027 |
| H0, surdispersé (θ = 1), µ = 5, n = 20 | 0,020 / 0 / 0 | 0,500 / 0,450 / 0,333 | 0,070 | 1,000 |
| H0, surdispersé (θ = 0,5), µ = 3, n = 15 | 0,030 / 0,010 / 0,003 | 0,557 / 0,460 / 0,330 | 0,060 | 1,000 |
| H1, surdispersé (θ = 2), µ = 4, 4, 8 | 0,577 / 0,447 / 0,243 | 0,943 / 0,933 / 0,903 | 0,750 | 1,000 |
| H1, Poisson, µ = 4, 4, 6 | 0,717 / 0,623 / 0,413 | 0,720 / 0,633 / 0,440 | 0,783 | 0,017 |

**Lecture.**
- **Sous H0 surdispersé, Poisson seul invente une différence.** Il dépasse P > 0,95 dans un tiers des jeux. La comparaison des lois ramène ce taux à 0 à 0,3 %, et la binomiale négative est reconnue dans 100 % des jeux.
- **Sous Poisson vrai,** les deux voies se valent.
- **À effet réel et surdispersé,** m.bay est prudent : 45 % au-delà de 0,75, contre 75 % de p < 0,05 pour le rapport de vraisemblance. C'est le prix de l'a priori d'information unitaire. Le chiffre de Poisson seul (93 %) n'est pas une puissance, puisqu'il invente aussi sous H0.

### Approximation BIC contre BAS (10 novembre, audit de la session parente)

Quatre jeux de Poisson à 3 groupes, log BF de « différence » contre « aucune » :

| Cas | BIC (m.bay) | BAS g-prior (g = n, information unitaire) | BAS robust | BAS bic.prior |
|---|---:|---:|---:|---:|
| µ = 5, 5, 5, n = 20 | −4,04 | −4,06 | −4,06 | −4,04 |
| µ = 4, 4, 6, n = 20 | 3,69 | 3,55 | 3,46 | 3,69 |
| µ = 3, 5, 7, n = 15 | 11,46 | 11,13 | 10,91 | 11,46 |
| µ = 2, 2, 4, n = 30 | 1,41 | 1,33 | 1,28 | 1,41 |

L'écart à la vraisemblance marginale sous a priori g (BAS, approximation de Laplace) atteint 0,33 en log au plus, et 0,55 avec l'a priori robuste. On reste loin du seuil fixé, un facteur de Bayes de 3 (log 3 = 1,10). L'approximation est gardée. Script : `scratchpad/bicbas.R`, non versionné ; à reproduire avec `bas.glm(betaprior = g.prior(n))`.

### Puissance comparée (audit de la session parente, 100 jeux par scénario)

| Scénario | m.bay P>0,75 | m.bay P>0,95 | Test fréquentiste p<0,05 |
|---|---:|---:|---:|
| Poisson, rapport 1,5, n = 10 | 23 % | 10 % | 39 % (quasipoisson) |
| Binomiale négative k = 2, n = 10 | 10 % | – | 23 % |
| Binomiale négative k = 2, n = 20 | 12 % | – | 32 % |

Sous H0 Poisson : P > 0,75 dans 2 % des jeux, P > 0,95 dans 0 %. **Lecture.** L'a priori d'information unitaire est connu pour être conservateur. Il protège des fausses différences, surtout sous surdispersion (§13, agent_60), au prix d'une puissance environ moitié moindre que les tests fréquentistes. Le récit et le bilan ne le disent pas encore. À décider : une note de puissance, ou un a priori moins large.

## 14. Réponse ordonnée et routage des entiers (D17, voie 3), 10 novembre

`test-bay-ordinal.R` :
- **Routage** sur les 4 cas demandés par la session parente :
  - comptages surdispersés : comptages ;
  - tailles entières autour de 170 : continu ;
  - notes sur 20 : continu, dit « score borné » ;
  - Likert 1 à 5 : ordinal, annoncé.
- **Options.** `response = "ordinal"` et `response = "continuous"` forcent la voie.
- **Facteur de Bayes.** Il égale −BIC/2 de `MASS::polr` avec et sans le facteur, à 1e-8.
- **Rapports des cotes.** Ils égalent exp(coef) de `polr` à 1e-6.
- **Mise en forme.** Les lignes tiennent dans 80 colonnes.

### Calibration (`agent_61`, N = 300 jeux par scénario, 3 groupes)

Le Likert est tiré par seuils d'une variable latente normale. Résultats : `scripts/resultats/agent_61.txt`.

| Scénario | m.bay P>0,5 / P>0,75 / P>0,95 | polr, rapport de vraisemblance p<0,05 | Kruskal p<0,05 | Voie ordinale choisie seule |
|---|---|---:|---:|---:|
| 5 niveaux, n = 10, H0 | 0,047 / 0,027 / 0,003 | 0,077 | 0,053 | 1,000 |
| 5 niveaux, n = 30, H0 | 0,003 / 0 / 0 | 0,060 | 0,043 | 1,000 |
| 7 niveaux, n = 10, H0 | 0,063 / 0,027 / 0,010 | 0,083 | 0,060 | 1,000 |
| 7 niveaux, n = 30, H0 | 0,003 / 0 / 0 | 0,047 | 0,037 | 1,000 |
| 5 niveaux, n = 10, décalage 0,8 | 0,337 / 0,190 / 0,073 | 0,407 | 0,327 | 1,000 |
| 5 niveaux, n = 30, décalage 0,5 | 0,263 / 0,163 / 0,060 | 0,447 | 0,437 | 1,000 |
| 7 niveaux, n = 30, décalage 0,5 | 0,253 / 0,153 / 0,060 | 0,443 | 0,420 | 1,000 |
| 5 niveaux, n = 30, cotes non proportionnelles (même centre, dispersion × 2,5) | 0,007 / 0,003 / 0 | 0,080 | 0,053 | 1,000 |

**Lecture.**
- **Sous H0,** une différence est annoncée au-delà de « faible » dans 0 à 2,7 % des jeux.
- **Avec un effet,** la méthode est prudente, comme la voie 2 : 15 à 19 % au-delà de 0,75, contre 33 à 45 % pour les tests fréquentistes. C'est l'a priori d'information unitaire.
- **Cotes non proportionnelles** (même centre, dispersion différente) : ni m.bay, ni le rapport de vraisemblance, ni Kruskal ne trouvent de différence. Le récit ne l'invente pas ; il ne voit pas non plus l'écart de dispersion, puisque l'hypothèse n'est pas testée.

### BIC contre Laplace sous a priori propre (`agent_62`)

La référence exacte prévue, brms (`cumulative`) avec bridgesampling, n'a pas pu tourner : la compilation Stan échoue sur ce poste (R 4.6 attend Rtools 4.5). Elle est remplacée par l'approximation de Laplace sous l'a priori propre d'information unitaire, β ~ N(0, nV), que le BIC approche ([KassWasserman1995]).

| Cas | BIC | Laplace, information unitaire | Écart |
|---|---:|---:|---:|
| H0, n = 20 | −3,62 | −3,65 | 0,02 |
| décalage 0,6, n = 20 | −3,30 | −3,33 | 0,03 |
| 0 ; 0,5 ; 1, n = 15 | 3,63 | 3,47 | 0,16 |
| décalage 0,4, n = 40 | −3,97 | −3,99 | 0,02 |

L'écart atteint 0,16 en log au plus, loin de log 3. **Limite** : ce n'est pas une vraisemblance marginale exacte par tirage, seulement une approximation plus fine, sous le même a priori.

## 15. Variances inégales à plusieurs facteurs : le défaut mesuré (`agent_64`), 10 novembre

On rejoue la mise en défaut de D6 dans un plan 2 × 2 sous H0 (aucun effet), avec des effectifs et des écarts-types inégaux par case : n = (5, 5, 5, 40), écarts-types (4, 4, 4, 1). N = 200 jeux. Résultats : `scripts/resultats/agent_64.txt`.

| Méthode | Au moins un terme dit « oui » |
|---|---:|
| m.bay, P > 0,95 | 12,5 % |
| m.bay, P > 0,75 | 33,5 % |
| aov, p < 0,05 | 58,0 % |
| Welch sur les 4 cases, p < 0,05 | 5,0 % |

**Lecture.** Le moteur à plusieurs facteurs, qui suppose une variance commune, se trompe ici comme le faisait la voie normale seule à un facteur avant D6 : 12,5 % de faux « oui » à P > 0,95. Le critère de réussite de la prochaine étape est de ramener ce taux au niveau d'un facteur avec variances par groupe, soit 0,7 % à P > 0,95 (D6).

### Protection provisoire (10 novembre)

Les variances par case sont comparées par le moteur à un facteur, avec les cases pour groupes (lois normales à variance commune ou par groupe). Quand les variances inégales sont probables (≥ 0,75), le soutien de toute réponse en faveur d'un effet est plafonné à « faible », et le verdict reste prudent. Les probabilités ne changent pas.

| Plan | Mesure | Avant | Après |
|---|---|---:|---:|
| 2 × 2, H0, variances inégales (`agent_64`) | « oui fort » affiché, P > 0,95 sans plafond | 12,5 % | 0 % |
| 2 × 2, variances égales, effet de A de 0,8 (`agent_65`) | prudence déclenchée à tort | – | 5 % |
| 2 × 2, variances égales, effet de A de 0,8 (`agent_65`) | P(A) > 0,95 affiché « fort » | 17,5 % | 17,5 % |

**Lecture.** C'est une protection d'affichage, pas la correction. La probabilité de 12,5 % au-delà de 0,95 sous H0 demeure dans `$effects`. La correction attendue moyenne les lois à plusieurs facteurs, sur le critère de D6 : 0,7 % à P > 0,95, avec au plus un tiers de puissance perdue à variances égales.


## 16. Variances par case moyennées à plusieurs facteurs (D19, `agent_66`), 10 novembre

Correction de D19, option 1 (`R/sys_bay_fact_laws.R`), rejouée sur les plans de §15 avec les mêmes graines, **sans aucun plafond**. Résultats : `scripts/resultats/agent_66_*.txt`.

| Plan | Mesure | Variance commune seule | Lois moyennées |
|---|---|---:|---:|
| 2 × 2, H0, n = (5, 5, 5, 40), écarts-types (4, 4, 4, 1) | un terme à P > 0,95 | 12,5 % (25/200) | **0,9 %** (7/800) |
| idem | un terme à P > 0,75 | 33,5 % | 10,0 % (80/800) |
| 2 × 2, variances égales, effet de A de 0,8, n = 10 par case | P(A) > 0,95 | 17,5 % | **17,5 %** |
| idem | P(A) > 0,75 | 43,5 % | 45,0 % |
| idem | B ou A:B à P > 0,95 | – | 0 % |

Détail des 800 jeux sous H0 : 3/200 (graines d'agent_64), 1/300 et 3/300 (graines suivantes).

**Lecture.**
- Sous H0, le taux de faux « oui » à P > 0,95 passe de 12,5 % à 0,9 % (intervalle binomial à 95 % : 0,4 % à 1,8 %). C'est le niveau visé de 0,7 % (D6, un facteur), sans l'atteindre exactement ; l'écart reste dans l'intervalle.
- À variances égales, aucune puissance n'est perdue. P(variance par case) médiane : 0,22.
- La loi à variance par case est toujours calculée : 100 % des jeux, aucun échec d'échantillonnage, aucun plafond résiduel.
- Temps par appel (sans sensibilité, N = 40 ou 55) : médiane 1,8 à 1,9 s, maximum 4,4 s.
- Contrôle du noyau : avec w = 1, on retrouve `lmBF` et le calcul exact à un facteur à 0,05 log près (`test-bay-fact-laws.R`).

## 17. Post-hoc à plusieurs facteurs (D20), 10 novembre

Vérifié sur des cas construits (`test-bay-fact-laws.R`, `test-bay-repeated.R`, `test-bay-story.R`) :
- 3 × 2 avec un effet de A seul : la famille détaillée est A ; a1/a2 sont différents (P > 0,9) et a2/a3 penchent vers semblables ; A:B est nommé comme non comparé.
- 2 × 2 avec une seule case décalée : A:B probable, les 4 cases sont comparées.
- Split-plot (groupe inter × temps intra) : temps puis groupe, sujets en nuisance ; pas d'écart standardisé pour le facteur inter-sujets.
- ANCOVA à deux facteurs : les niveaux de A sont comparés à B et x égaux.

Temps : 6 à 9 s pour un appel complet avec sensibilité (3 × 2, 2 × 2, split-plot 16 sujets × 3 temps, ANCOVA 3 × 2 à N = 60). **Non fait** : calibration par simulation du post-hoc à plusieurs facteurs (fausses différences de paires sous H0).

## 18. Mesures répétées avec répliques (D21, `agent_67`), 10 novembre

10 sujets × 3 conditions × 4 répliques. Écarts-types : sujet 1, sujet × condition 0,7, résiduel 1. N = 200 jeux par scénario. Résultats : `scripts/resultats/agent_67_*.txt`.

| Scénario | m.bay P > 0,75 | m.bay P > 0,95 | Sujet seul (répliques indépendantes) P > 0,95 | ANOVA RM sur les moyennes, p < 0,05 |
|---|---:|---:|---:|---:|
| H0 | 1,0 % | 0 % | 7,0 % | 6,0 % |
| Effet (0 ; 0,5 ; 0,6) | 17,0 % | 4,0 % | 33,0 % | 26,5 % |

**Lecture.**
- Le terme sujet × condition est nécessaire : sans lui, 7 % de faux « oui » forts sous H0.
- Avec lui, aucun faux « oui » fort, au prix d'une puissance prudente, comme sur les autres voies.
- Temps par appel (sans sensibilité, N = 120) : médiane 4,3 s, maximum 6,1 s ; environ 5 s avec la sensibilité.

## 19. MANOVA à un facteur : choix du moteur (D22, `agent_68`), 11 novembre

3 groupes de 15 ; corrélations 0,3 entre réponses, sauf scénarios « cor » (0,8 entre y1 et y2, 0 ailleurs). 100 jeux par scénario. Résultats : `scripts/resultats/agent_68_*.txt`.

| Scénario | Pillai p < 0,05 | A BFpack P > 0,75 / 0,95 | B BIC P > 0,75 / 0,95 | **C lmBF empilé** P > 0,75 / 0,95 |
|---|---:|---:|---:|---:|
| H0, d = 3 | 3 % | 0 / 0 % | 0 / 0 % | 2 / 0 % |
| H0, d = 10 | 5 % | 0 / 0 % | 0 / 0 % | 3 / 0 % |
| H0, corrélations inégales | 3 % | 0 / 0 % | 0 / 0 % | 1 / 0 % |
| Effet 0,8 sur une réponse, d = 3 | 41 % | 6 / 2 % | 6 / 1 % | 16 / 4 % |
| Effet 0,8 sur une réponse, d = 10 | 21 % | 0 / 0 % | 0 / 0 % | 6 / 2 % |
| Effet 0,5 sur les 3 réponses | 27 % | 3 / 1 % | 2 / 1 % | 21 / 7 % |
| Effet 0,8 sur y1, corrélée 0,8 à y2 | 84 % | 29 / 14 % | 23 / 13 % | 15 / 4 % |

**Lecture.** C n'invente pas de différence, et il est le plus puissant, sauf quand l'effet passe par un contraste entre réponses corrélées : Pillai y exploite la covariance, C ne le peut pas. Un appel complet avec la sensibilité prend environ 15 s (3 réponses, N = 45).

## 20. MANOVA à deux facteurs (D23, `agent_69`), 11 novembre

2 × 2, 3 réponses (corrélations 0,3), n = 12 par case, 60 jeux par scénario, sans sensibilité. Résultats : `scripts/resultats/agent_69_*.txt`.

| Scénario | Un terme P > 0,75 | Un terme P > 0,95 | P(A) > 0,75 | P(A) > 0,95 | B ou A:B P > 0,95 | Pillai sur A, p < 0,05 |
|---|---:|---:|---:|---:|---:|---:|
| H0 | 1,7 % | 0 % | 0 % | 0 % | 0 % | 5,0 % |
| Effet de A (0,8 sur y1) | 38,3 % | 18,3 % | 38,3 % | 16,7 % | 1,7 % | 76,7 % |
| Effet de A réparti (0,45 sur les 3 réponses) | 16,7 % | 5,0 % | 16,7 % | 3,3 % | 1,7 % | 31,7 % |

**Lecture.** Aucun faux « oui » fort sous H0 ; puissance prudente, comme à un facteur (D22). Le rapport m.bay / Pillai est le même avec l'effet sur une réponse (0,50) et réparti (0,53) : la perte ne tient pas à l'effet porté par une seule réponse. Temps médian 21 s sans sensibilité, sous la charge de deux simulations ; 50 à 70 s avec la sensibilité, ligne d'attente affichée. Le test des 80 colonnes et R CMD check (Status OK) portent sur 3135136.

## 21. Queues lourdes à plusieurs facteurs (D24, `agent_70`), 11 novembre

2 × 2, n = 12 par case, bruit t3 (écart-type 1). « Avant » : D19 seul, 150 jeux. « Après » : D24, les 100 premiers jeux des mêmes graines. Résultats : `scripts/resultats/agent_70_*_avant.txt` et `*_apres.txt`. Dans le fichier « après », la colonne « m.bay facteurs (normale) » désigne le moteur à plusieurs facteurs, Student compris.

| Scénario | Mesure | Avant | Après | Voie à un facteur sur A | aov / aov sur rangs, p < 0,05 |
|---|---|---:|---:|---:|---:|
| H0 | P(A) > 0,95 | 0 % | 0 % | 0 % | 8 % / 6 % |
| Effet de A, 0,8 | P(A) > 0,75 | 70 % | 86 % | 88 % | 84 % / 96 % |
| Effet de A, 0,8 | P(A) > 0,95 | 44 % | 64 % | 71 % | |

Temps : environ 9 s par appel avec la sensibilité (2 × 2, N = 48).
