# Validation de m.bay

**`m.bay()` n'existe pas encore : rien ici ne le valide.** Sections 1 à 5 : 5 novembre ; section 6 : 6 novembre. Ce fichier documente ce qui a été exécuté sur les moteurs candidats et le socle. R 4.6.1, BayesFactor 0.9.12-4.8, Windows 11, 5 novembre 2026. Scripts dans `scripts/`, sorties dans `scripts/resultats/`. Détails : `rapport_moteurs.md`, `agent_2x_calibration_rapport.md`, `RELECTURE_ADVERSE_2026-11-05.md`.

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
