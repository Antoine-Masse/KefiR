# Validation de m.bay

**`m.bay()` n'existe pas encore : rien ici ne le valide.** Ce fichier documente ce qui a été exécuté sur les moteurs candidats et le socle. R 4.6.1, BayesFactor 0.9.12-4.8, Windows 11, 5 novembre 2026. Scripts dans `scripts/`, sorties dans `scripts/resultats/`. Détails : `rapport_moteurs.md`, `agent_2x_calibration_rapport.md`, `RELECTURE_ADVERSE_2026-11-05.md`.

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
- Hétéroscédasticité, H0 : n=(5,5,40), sd=(4,4,1), **P>0,95 dans 16 %** [12 ; 21], N=300 ; aov 46 %, Welch 7 %, Kruskal 18 % [recalculé, `agent_12`]. Le relecteur a mesuré 25,5 % sur un autre tirage : écart non élucidé. Autre plan, SD 3 sur n=10 de 10,30,50 : P>0,95 8,0 % [5,4 ; 11,6], aov 26,3 % [agent].

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

Multiplicité 2³ (calcul maison du relecteur, non recoupé avec anovaBF) ; ROPE f et MANOVA en √d (relecteur seul) ; toute voie non gaussienne ; hypothèses de sphéricité ; BAS ; brms avec Rtools.
