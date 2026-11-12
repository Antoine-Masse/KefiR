# Validation de cor.bay

R 4.6.1, BayesFactor 0.9.12-4.8, bayestestR 0.20.0, Windows 11, 11 novembre 2026. Scripts dans `scripts/`, sorties dans `scripts/resultats/` (`.rds`, `.log`). `proto.R` recopie les calculs du paquet sans le récit (pour aller vite) ; `sim_auto.R` appelle le code du paquet lui-même. Tous les taux sont des pourcentages de jeux simulés.

**R CMD check** (`--no-manual --ignore-vignettes`, suite de tests complète de KefiR comprise) sur le code de 06ad49b : **Status OK** (seule information : xref vers `ez`, absent de mbay-lib, comme pour m.bay). Tests de cor.bay : `tests/testthat/test-cor.bay.R`, 87 assertions.

Lecture : « P > 0,95 » est le seuil d'un « soutien fort » et « P > 0,75 » celui d'un « soutien modéré » (même échelle que m.bay). Ce ne sont pas des seuils de p-value ; la comparaison avec `cor.test` à p < 0,05 est donnée pour situer, pas pour égaler.

## 1. Accord numérique (`check_agree.R`, 200 jeux, n de 5 à 200, rho de -0,8 à 0,8)

| Calcul | Contre | Écart sur log BF |
|---|---|---|
| BF bilatéral = moyenne des deux moitiés de `correlationBF(nullInterval = c(0, 1))` | `correlationBF` bilatéral direct | médiane 0,0005, max 0,016 (erreur d'intégration propre à BayesFactor) |
| BF de Kendall (`integrate`, intervalle coupé autour du sommet) | méthode des trapèzes, 400 001 points sur [-1, 1] | 7e-12 |
| Probabilités des trois hypothèses | formule fermée à partir des deux moitiés | 1e-10 (test `test-cor.bay.R`) |
| A priori sur tau, kappa = 1 | pi/4 cos(pi tau / 2), éq. 10 de van Doorn et al. 2018 | 1e-12 (test) |

## 2. Sous H0 : faux « oui » par méthode (`sim_h0.R`, 500 jeux par case)

Aucun lien entre x et y. « Oui » = P(un lien, d'un signe ou de l'autre) au-delà du seuil. Données : normales indépendantes ; t à 3 ddl indépendantes ; **t3 bivariée** (x et y partagent une même échelle aléatoire : non corrélées mais dépendantes, le cas qui trompe Pearson) ; log-normales indépendantes ; échelle 1 à 5 indépendante ; normales avec **une paire aberrante** à (6, 6).

| Données | n | Pearson P > 0,75 | P > 0,95 | Kendall P > 0,75 | P > 0,95 | Spearman P > 0,95 | `cor.test` Pearson p < 0,05 | Kendall p < 0,05 |
|---|---|---|---|---|---|---|---|---|
| normales | 10 | 2,6 | 0,0 | 2,4 | 0,0 | 0,2 | 5,2 | 6,4 |
| normales | 20 | 2,0 | 0,2 | 3,0 | 0,0 | 0,4 | 4,4 | 5,6 |
| normales | 50 | 3,4 | 0,2 | 2,6 | 0,0 | 0,2 | 5,4 | 6,2 |
| normales | 100 | 2,2 | 0,0 | 2,4 | 0,2 | 0,2 | 5,4 | 4,6 |
| t3 indépendantes | 10 | 3,4 | 0,2 | 1,2 | 0,0 | 0,0 | 6,4 | 6,0 |
| t3 indépendantes | 20 | 4,4 | 0,2 | 4,4 | 0,6 | 0,8 | 6,0 | 5,8 |
| t3 indépendantes | 50 | 2,6 | 0,2 | 3,0 | 0,2 | 0,2 | 4,6 | 5,6 |
| t3 indépendantes | 100 | 2,8 | 0,0 | 1,4 | 0,0 | 0,0 | 5,4 | 5,4 |
| **t3 bivariée** | 10 | 11,4 | 3,0 | 3,4 | 0,0 | 0,4 | 16,4 | 9,2 |
| **t3 bivariée** | 20 | 15,0 | 5,6 | 3,6 | 0,4 | 0,6 | 19,2 | 5,8 |
| **t3 bivariée** | 50 | 24,6 | 10,2 | 6,2 | 0,4 | 0,4 | 31,4 | 9,8 |
| **t3 bivariée** | 100 | 21,4 | 11,2 | 3,8 | 0,4 | 0,4 | 28,6 | 7,8 |
| log-normales | 10 | 2,8 | 0,6 | 2,0 | 0,0 | 0,2 | 4,8 | 5,2 |
| log-normales | 20 | 4,2 | 1,2 | 2,6 | 0,2 | 0,2 | 5,2 | 5,2 |
| log-normales | 50 | 4,2 | 1,6 | 3,0 | 0,0 | 0,0 | 5,8 | 6,0 |
| log-normales | 100 | 2,6 | 1,2 | 1,4 | 0,0 | 0,0 | 4,2 | 5,6 |
| échelle 1 à 5 * | 10 | 2,4 | 0,0 | 1,2 | 0,0 | 0,0 | 4,0 | 3,8 |
| échelle 1 à 5 * | 20 | 2,4 | 0,0 | 2,6 | 0,0 | 0,0 | 5,2 | 4,8 |
| échelle 1 à 5 * | 50 | 3,4 | 0,2 | 3,4 | 0,6 | 0,6 | 5,8 | 4,4 |
| échelle 1 à 5 * | 100 | 1,2 | 0,2 | 1,6 | 0,2 | 0,2 | 4,2 | 4,0 |
| **une paire aberrante** | 10 | 90,0 | 18,8 | 2,4 | 0,0 | 0,2 | 95,6 | 8,8 |
| **une paire aberrante** | 20 | 95,6 | 57,2 | 4,8 | 0,2 | 0,4 | 98,0 | 7,0 |
| **une paire aberrante** | 50 | 88,4 | 57,6 | 3,6 | 0,4 | 0,4 | 94,2 | 6,6 |
| **une paire aberrante** | 100 | 68,0 | 33,6 | 2,4 | 0,2 | 0,0 | 81,6 | 5,2 |

**Lecture.**
- Sur des données bien faites (normales, t indépendantes, log-normales), aucune méthode ne dépasse 1,6 % de faux « oui » à P > 0,95, et 4,4 % à P > 0,75 ; la référence de m.bay est 0,7 % à P > 0,95 (D6, D26).
- **t3 bivariée** : Pearson bayésien dit « oui » à P > 0,95 dans 10 à 11 % des jeux à n ≥ 50, comme `cor.test` (29 à 31 % à p < 0,05) : le modèle normal bivarié est faux et les deux le subissent. Kendall reste à 0,4 %. C'est la raison de la voie Kendall pour les queues lourdes ; ce que la règle `auto` en attrape est en §6.
- **Une paire aberrante** : Pearson dit « oui » à P > 0,95 dans 19 à 58 % des jeux ; Kendall dans 0,4 % au plus.
- **Échelle 1 à 5** (* : relancé après la correction des ex-aequo, `sim_h0.R 500 likert _ties`, autre tirage) : Kendall au plus 0,6 % à P > 0,95. **Avant la correction**, il montait à 8,0 % à P > 0,75 et 2,0 % à P > 0,95 (n = 100) : T* avait un écart-type de 1,2 sous H0 (DECISIONS, C4). Spearman sur rangs y est prudent mais perd de la puissance (§5).

## 3. Calibration sous l'a priori (`sim_power.R`, partie 2)

Jeux tirés de l'a priori lui-même : moitié sans lien, moitié avec rho tiré de la bêta étirée (kappa = 1/3), données normales bivariées ; 1200 jeux par n. Si les probabilités annoncées sont justes, la fréquence d'un vrai lien dans chaque classe égale la probabilité annoncée. C'est exact par construction pour Pearson si le calcul est juste ; pour Kendall et Spearman sur rangs, c'est la mesure de l'approximation.

| Méthode | n | classe (0,2 ; 0,4] | (0,4 ; 0,6] | (0,6 ; 0,8] | (0,8 ; 1] |
|---|---|---|---|---|---|
| Pearson | 10 | 0,388 / 0,386 | 0,472 / 0,490 | 0,684 / 0,721 | 0,897 / 0,863 |
| Pearson | 20 | 0,345 / 0,350 | 0,475 / 0,443 | 0,700 / 0,694 | 0,922 / 0,923 |
| Pearson | 50 | 0,282 / 0,287 | 0,480 / 0,465 | 0,697 / 0,718 | 0,963 / 0,964 |
| Kendall | 10 | 0,390 / 0,418 | 0,466 / 0,485 | 0,676 / 0,693 | 0,871 / 0,849 |
| Kendall | 20 | 0,348 / 0,350 | 0,475 / 0,451 | 0,697 / 0,709 | 0,924 / 0,938 |
| Kendall | 50 | 0,288 / 0,295 | 0,487 / 0,505 | 0,694 / 0,651 | 0,960 / 0,962 |
| Spearman (rangs) | 10 | 0,387 / 0,414 | 0,466 / 0,488 | 0,692 / 0,687 | 0,883 / 0,844 |
| Spearman (rangs) | 20 | 0,346 / 0,353 | 0,475 / 0,465 | 0,687 / 0,667 | 0,928 / 0,937 |
| Spearman (rangs) | 50 | 0,281 / 0,288 | 0,487 / 0,519 | 0,697 / 0,713 | 0,963 / 0,962 |

(annoncé / fréquence observée ; classes de 50 à 680 jeux, donc un écart de 0,03 à 0,05 est dans le bruit.) **Les trois méthodes sont calibrées dès n = 10**, y compris l'approximation normale de Kendall que l'article ne garantit qu'à partir de n = 10 sous H0.

### Échelle 1 à 5 (`sim_calib_likert.R`, 1500 jeux par n) : Kendall, avant et après la correction des ex-aequo

Couples normaux latents (rho tiré de l'a priori, ou 0) coupés en 5 classes ; « vrai lien » = rho latent non nul.

| n | classe | avant : annoncé / observé | après : annoncé / observé |
|---|---|---|---|
| 20 | (0,6 ; 0,8] | 0,699 / 0,646 | 0,683 / 0,719 |
| 20 | (0,8 ; 1] | 0,928 / 0,842 | 0,917 / 0,941 |
| 50 | (0,6 ; 0,8] | 0,699 / 0,575 | 0,694 / 0,692 |
| 50 | (0,8 ; 1] | 0,961 / 0,910 | 0,950 / 0,963 |
| 100 | (0,4 ; 0,6] | 0,491 / 0,385 | 0,500 / 0,441 |
| 100 | (0,6 ; 0,8] | 0,698 / 0,520 | 0,695 / 0,656 |
| 100 | (0,8 ; 1] | 0,976 / 0,918 | 0,968 / 0,967 |

Avant, Kendall était trop sûr d'un lien (0,70 annoncé pour 0,52 réel à n = 100) ; après, les écarts restent dans le bruit des classes (90 à 400 jeux), sauf peut-être (0,4 ; 0,6] à n = 100 (0,50 / 0,44, 177 jeux), à surveiller.

## 4. Largeur de l'a priori : kappa = 1/3 (choisi) contre kappa = 1 (JASP) (`sim_kappa.R`, 300 jeux par case, Pearson)

| rho | n | P > 0,95, k = 1/3 | P > 0,95, k = 1 | P > 0,75, k = 1/3 | P > 0,75, k = 1 | « aucun lien » à 75 %, k = 1/3 | idem, k = 1 |
|---|---|---|---|---|---|---|---|
| 0 | 20 | 0,0 | 0,0 | 4,0 | 3,7 | 0,0 | 44,3 |
| 0 | 50 | 1,0 | 1,0 | 3,7 | 2,7 | 24,0 | 74,3 |
| 0,3 | 20 | 2,0 | 3,3 | 18,7 | 15,7 | 0,0 | 21,3 |
| 0,3 | 50 | 21,7 | 18,3 | 47,7 | 40,0 | 2,7 | 14,7 |
| 0,5 | 20 | 18,7 | 20,0 | 53,3 | 49,0 | 0,0 | 4,7 |
| 0,5 | 50 | 80,7 | 78,0 | 93,3 | 92,3 | 0,0 | 0,0 |

Pour dire « un lien », les deux largeurs se valent (écarts de 1 à 8 points). La différence est du côté de l'absence : avec kappa = 1, « aucun lien » atteint 75 % dans 44 % des jeux sans lien dès n = 20, **mais aussi dans 21 % des jeux où rho = 0,3** ; avec kappa = 1/3, presque jamais à n = 20. kappa = 1 attend des liens forts, donc il « innocente » vite un lien modéré. Pour des TP où rho = 0,3 est un lien réel et courant, c'est une erreur plus grave qu'un « indécis » : **on garde kappa = 1/3** (C2). Le prix est connu : sous 50 paires, cor.bay conclura rarement à l'absence de lien ; il dira « indécis », ce qui est la vraie réponse à cette taille.

## 5. Puissance comparée à `cor.test` (`sim_power.R`, 300 jeux par case)

P(lien positif) > 0,75 et > 0,95 ; `cor.test` bilatéral à p < 0,05. P = Pearson, K = Kendall, S = Spearman sur rangs.

| Données | rho | n | P > 0,75 | P > 0,95 | K > 0,75 | K > 0,95 | S > 0,95 | Pearson p < 0,05 | Kendall p < 0,05 | Spearman p < 0,05 |
|---|---|---|---|---|---|---|---|---|---|---|
| normales | 0,2 | 20 | 6,0 | 1,0 | 4,3 | 0,0 | 0,3 | 10,3 | 8,3 | 8,3 |
| normales | 0,2 | 50 | 19,7 | 3,7 | 17,3 | 3,0 | 6,0 | 30,3 | 27,7 | 29,0 |
| normales | 0,2 | 100 | 34,0 | 12,0 | 34,3 | 13,0 | 13,0 | 52,0 | 47,3 | 47,3 |
| normales | 0,3 | 20 | 12,7 | 1,3 | 12,3 | 1,0 | 1,3 | 22,7 | 20,7 | 22,0 |
| normales | 0,3 | 50 | 45,7 | 18,3 | 37,7 | 12,3 | 14,7 | 59,3 | 56,3 | 57,7 |
| normales | 0,3 | 100 | 77,0 | 50,0 | 73,0 | 43,3 | 44,0 | 85,0 | 81,7 | 80,7 |
| normales | 0,5 | 20 | 58,7 | 16,3 | 44,7 | 9,3 | 14,0 | 68,0 | 56,7 | 58,3 |
| normales | 0,5 | 50 | 94,0 | 75,7 | 91,7 | 66,0 | 69,7 | 97,0 | 96,0 | 96,3 |
| normales | 0,5 | 100 | 99,7 | 98,0 | 99,0 | 97,7 | 97,0 | 100,0 | 99,3 | 99,7 |
| t3 bivariée | 0,2 | 20 | 18,0 | 3,3 | 8,3 | 0,3 | 0,3 | 27,0 | 14,0 | 11,7 |
| t3 bivariée | 0,2 | 50 | 31,7 | 13,7 | 15,7 | 5,7 | 5,0 | 42,0 | 24,0 | 20,7 |
| t3 bivariée | 0,2 | 100 | 37,3 | 27,3 | 33,3 | 15,3 | 14,0 | 52,3 | 45,0 | 41,7 |
| t3 bivariée | 0,3 | 20 | 27,3 | 11,7 | 17,0 | 3,7 | 4,0 | 36,3 | 23,7 | 24,0 |
| t3 bivariée | 0,3 | 50 | 48,7 | 31,0 | 44,7 | 17,3 | 17,3 | 56,3 | 55,0 | 53,3 |
| t3 bivariée | 0,3 | 100 | 62,3 | 51,3 | 63,0 | 41,0 | 38,3 | 70,7 | 76,0 | 71,7 |
| t3 bivariée | 0,5 | 20 | 58,0 | 32,0 | 43,7 | 14,0 | 18,3 | 64,7 | 55,7 | 55,7 |
| t3 bivariée | 0,5 | 50 | 79,7 | 68,3 | 86,3 | 62,7 | 60,0 | 83,3 | 91,3 | 88,3 |
| t3 bivariée | 0,5 | 100 | 93,0 | 86,3 | 99,0 | 97,7 | 96,7 | 95,3 | 99,3 | 99,0 |
| monotone (y = exp(2 z)) | 0,2 | 20 | 6,3 | 0,3 | 7,3 | 0,3 | 0,3 | 11,3 | 12,0 | 13,0 |
| monotone | 0,2 | 50 | 7,0 | 1,7 | 16,0 | 4,7 | 5,3 | 12,3 | 24,7 | 25,0 |
| monotone | 0,2 | 100 | 11,3 | 2,0 | 27,0 | 9,7 | 11,0 | 18,0 | 39,0 | 40,3 |
| monotone | 0,3 | 20 | 8,7 | 1,3 | 13,3 | 0,7 | 2,3 | 14,3 | 21,7 | 20,7 |
| monotone | 0,3 | 50 | 17,7 | 4,3 | 40,0 | 12,0 | 12,3 | 26,3 | 54,3 | 55,3 |
| monotone | 0,3 | 100 | 26,7 | 9,7 | 66,7 | 38,7 | 38,7 | 39,0 | 79,3 | 80,3 |
| monotone | 0,5 | 20 | 21,7 | 2,7 | 45,3 | 13,0 | 19,7 | 32,7 | 59,0 | 59,3 |
| monotone | 0,5 | 50 | 50,7 | 19,3 | 88,0 | 66,0 | 65,0 | 63,7 | 91,7 | 92,3 |
| monotone | 0,5 | 100 | 71,7 | 43,0 | 99,0 | 96,3 | 96,7 | 82,3 | 99,7 | 99,7 |
| échelle 1 à 5 * | 0,2 | 20 | 6,7 | 1,3 | 3,7 | 0,0 | 1,3 | 13,0 | 10,7 | 10,3 |
| échelle 1 à 5 * | 0,2 | 50 | 14,7 | 2,3 | 16,3 | 3,7 | 2,3 | 20,7 | 20,7 | 21,0 |
| échelle 1 à 5 * | 0,2 | 100 | 28,3 | 10,7 | 27,3 | 10,3 | 9,3 | 44,0 | 41,7 | 41,7 |
| échelle 1 à 5 * | 0,3 | 20 | 14,3 | 0,7 | 11,3 | 0,7 | 1,3 | 23,0 | 22,7 | 22,7 |
| échelle 1 à 5 * | 0,3 | 50 | 36,0 | 10,3 | 36,7 | 9,7 | 10,3 | 45,3 | 44,0 | 43,0 |
| échelle 1 à 5 * | 0,3 | 100 | 65,7 | 38,0 | 66,0 | 38,7 | 34,7 | 77,0 | 75,7 | 75,3 |
| échelle 1 à 5 * | 0,5 | 20 | 44,7 | 14,7 | 41,0 | 7,0 | 14,7 | 54,3 | 52,0 | 52,7 |
| échelle 1 à 5 * | 0,5 | 50 | 85,7 | 67,0 | 88,0 | 60,3 | 64,3 | 92,3 | 91,7 | 91,3 |
| échelle 1 à 5 * | 0,5 | 100 | 98,7 | 94,7 | 98,7 | 94,7 | 97,3 | 99,7 | 99,3 | 99,3 |

**Lecture.**
- Sur des données normales, « soutien modéré » (P > 0,75) se déclare un peu moins souvent que p < 0,05 (5 à 14 points de moins), « soutien fort » (P > 0,95) nettement moins : c'est un seuil plus exigeant, et l'a priori ramène les petits effets vers 0. Pearson bat Kendall de 0 à 14 points, comme attendu ([vanDoorn2018], fig. 4).
- Lien monotone courbe : Kendall double ou triple la puissance de Pearson (rho = 0,5, n = 50 : 66 % contre 19 % à P > 0,95). La règle `auto` y prend Kendall (valeurs lointaines de exp(2 z)) : voir §6.
- Échelle 1 à 5 (* : colonnes Kendall relancées après la correction des ex-aequo, `sim_power.R 300 likert _ties`) : Kendall corrigé est un peu moins puissant que Pearson bayésien (rho = 0,5, n = 50 : 60 % contre 65 à 67 % à P > 0,95). Sa puissance d'avant (81 %) venait de sa surconfiance. La règle garde Kendall sur les échelles courtes : il ne dépend que de l'ordre des notes, pas des écarts supposés entre elles, comme la voie ordinale de m.bay (D17) ; le prix est d'environ 5 points de puissance.
- t3 bivariée : Pearson « trouve » plus souvent un lien faible que Kendall, mais sa calibration sous H0 sur ces données est à lire en §2 avant de s'en réjouir.

## 6. La règle `auto` du paquet (`sim_auto.R`, 300 jeux par case, code du paquet)

`sim_auto.R` appelle les fonctions du paquet (contrôles, choix, facteurs de Bayes, contre-épreuve), version finale de la règle (valeurs lointaines, distance de Cook dès 20 paires, échelle courte). Six scénarios lancés séparément (`sim_auto.R 300 <scénario> _v2_<scénario>`), sorties `sim_auto_v2_*.rds`. Données comme en §2 et §5, plus « t3 marges » (copule normale à marges t3) ; « paire aberrante » : sous H0 à (6, 6), sous H1 à (6, -6), contre le sens du lien. Colonnes : % de jeux où Kendall est choisi ; « un lien » à P > 0,75 et 0,95 (signe quelconque) ; « lien positif » à 0,95 ; « aucun lien » à 75 % ; `cor.test` à p < 0,05.

| Données | rho | n | Kendall choisi | lien > 0,75 | lien > 0,95 | positif > 0,95 | aucun lien ≥ 0,75 | Pearson p < 0,05 | Kendall p < 0,05 |
|---|---|---|---|---|---|---|---|---|---|
| normales | 0 | 10 | 10,7 | 2,0 | 0,0 | 0,0 | 0,0 | 2,7 | 4,3 |
| normales | 0 | 20 | 8,7 | 3,0 | 0,0 | 0,0 | 0,0 | 4,0 | 4,7 |
| normales | 0 | 50 | 1,0 | 1,0 | 0,0 | 0,0 | 25,0 | 2,3 | 2,3 |
| normales | 0 | 100 | 0,7 | 3,3 | 0,3 | 0,3 | 67,3 | 4,3 | 4,0 |
| normales | 0,3 | 50 | 0,3 | 45,0 | 21,7 | 20,0 | 2,3 | 57,3 | 53,0 |
| normales | 0,5 | 50 | 2,3 | 93,0 | 77,3 | 77,3 | 0,3 | 95,0 | 93,7 |
| t3 bivariée | 0 | 10 | 31,7 | 4,3 | 0,0 | 0,0 | 0,0 | 12,3 | 7,7 |
| t3 bivariée | 0 | 20 | 52,3 | 8,3 | 1,7 | 0,7 | 0,0 | 20,7 | 9,3 |
| t3 bivariée | 0 | 50 | 76,0 | 6,0 | 0,3 | 0,0 | 12,0 | 28,3 | 9,0 |
| t3 bivariée | 0 | 100 | 92,3 | 5,3 | 1,3 | 0,0 | 56,3 | 37,0 | 8,0 |
| t3 bivariée | 0,3 | 50 | 77,3 | 43,3 | 19,3 | 18,7 | 2,7 | 54,7 | 51,3 |
| t3 bivariée | 0,5 | 50 | 72,0 | 87,0 | 66,0 | 65,7 | 0,7 | 83,7 | 91,3 |
| t3 marges | 0 | 20 | 59,0 | 2,7 | 0,0 | 0,0 | 0,0 | 4,0 | 4,7 |
| t3 marges | 0 | 100 | 95,3 | 2,7 | 0,0 | 0,0 | 61,7 | 7,0 | 4,0 |
| t3 marges | 0,3 | 50 | 77,0 | 43,7 | 16,3 | 16,3 | 1,7 | 50,7 | 53,0 |
| t3 marges | 0,5 | 50 | 73,0 | 90,3 | 68,7 | 68,0 | 0,0 | 92,0 | 93,7 |
| monotone | 0 | 50 | 99,0 | 1,0 | 0,0 | 0,0 | 11,0 | 4,3 | 2,3 |
| monotone | 0,3 | 50 | 99,3 | 41,3 | 16,7 | 16,7 | 1,7 | 31,3 | 53,0 |
| monotone | 0,5 | 50 | 99,0 | 90,3 | 67,7 | 67,3 | 0,0 | 58,7 | 93,7 |
| échelle 1 à 5 | 0 | 20 | 100 | 3,0 | 0,0 | 0,0 | 0,0 | 4,3 | 5,0 |
| échelle 1 à 5 | 0 | 100 | 100 | 2,7 | 0,3 | 0,3 | 43,3 | 4,0 | 3,7 |
| échelle 1 à 5 | 0,3 | 50 | 100 | 35,3 | 12,3 | 11,7 | 0,0 | 49,3 | 45,7 |
| échelle 1 à 5 | 0,5 | 50 | 100 | 86,0 | 59,3 | 58,3 | 0,0 | 89,3 | 89,7 |
| paire aberrante | 0 | 10 | 91,7 | 8,7 | 0,3 | 0,3 | 0,0 | 95,3 | 7,3 |
| paire aberrante | 0 | 20 | 100 | 5,0 | 0,0 | 0,0 | 0,0 | 98,0 | 8,0 |
| paire aberrante | 0 | 50 | 100 | 1,7 | 0,0 | 0,0 | 14,0 | 96,0 | 5,7 |
| paire aberrante | 0 | 100 | 100 | 2,0 | 0,7 | 0,7 | 60,7 | 87,0 | 3,7 |
| paire aberrante | 0,3 | 100 | 100 | 56,7 | 25,0 | 25,0 | 5,7 | 5,0 | 69,3 |
| paire aberrante | 0,5 | 50 | 100 | 76,7 | 42,3 | 41,3 | 0,0 | 7,7 | 85,7 |
| paire aberrante | 0,5 | 100 | 100 | 98,3 | 96,0 | 96,0 | 0,0 | 10,7 | 99,3 |

Table complète (66 cases) : `Rscript analyse_auto.R`.

**Lecture.**
- **Faux « oui » sous H0** : au plus **1,7 %** à P > 0,95 sur les six scénarios et les quatre tailles (t3 bivariée, n = 20), 0,7 % ailleurs ; au plus 8,7 % à P > 0,75 (paire aberrante, n = 10). `cor.test` de Pearson, sur les mêmes jeux, dit « significatif » dans 87 à 98 % des jeux à paire aberrante et 21 à 37 % des jeux t3 bivariés dès n = 20.
- **t3 bivariée** : là où Pearson seul disait « oui » à P > 0,95 dans 10 à 11 % des jeux (§2), la règle ramène à 0,3 à 1,7 % : elle prend Kendall dans 52 à 92 % de ces jeux, et ceux qui restent à Pearson sont ceux dont les queues ne se voient pas.
- **Données normales** : Kendall n'est choisi que dans 0,3 à 11 % des jeux (8,7 à 10,7 % sous 20 paires, par la règle des valeurs lointaines seule ; 1 % et moins ensuite). Avec l'ancien critère MCD, c'était 42 à 65 % (`sim_auto_mcd.rds`) : c'est la mesure qui a fait changer C3. La puissance est celle de Pearson bayésien (§5).
- **Lien monotone courbe** : Kendall choisi dans 99 % des jeux ; puissance à P > 0,95 de 67,7 % à rho = 0,5, n = 50, contre 19 % pour Pearson seul (§5).
- **Paire aberrante contre le sens du lien** : `cor.test` de Pearson ne voit plus le lien (5 à 11 %) ; cor.bay le retrouve (96 % à P > 0,95, rho = 0,5, n = 100).
- **Contre-épreuve** : aucun désaccord compté (« contradiction », C9) dans les 66 cases. Les désaccords nets entre les deux méthodes existent (jusqu'à 49 % des jeux à paire aberrante, rho = 0,3, n = 100), mais toujours quand Kendall a été choisi à cause de paires à part : ils sont dits « expliqués », la réponse de Pearson n'est pas retenue, et le bilan tranche. Avant ce réglage, ces mêmes jeux auraient rendu le bilan « indécis » (37 à 58 % des cas, `sim_auto_mcd.rds`).
- **Absence de lien** : « aucun lien » à 75 % n'apparaît qu'à partir de 50 paires (11 à 25 % des jeux sans lien), 100 paires (43 à 67 %) : conséquence voulue de kappa = 1/3 (§4).

## 7. « Aucun lien » soutenu alors que l'intervalle de rho exclut 0 ? (`sim_null_hdi.R`, 1000 jeux par case, 12 novembre)

Demande de la session parente, pour la cohérence avec m.bay et g.bay : le bilan peut-il dire « aucun lien, soutien fort » alors que l'intervalle de crédibilité de rho, dans le même bilan, exclut 0 ? Liens modérés réels, données normales, code du paquet (règle `auto`, facteurs de Bayes, niveau de soutien de `.corbay_lead_info`, intervalle du bilan sous l'a priori d'estimation ; la dernière colonne reprend l'intervalle sous l'a priori du test). Soutien modéré = niveau 3 (P ≥ 0,75 et BF ≥ 3), fort = niveau 4 (P ≥ 0,95).

| rho | n | « aucun lien » en tête | … soutien modéré | … soutien fort | IC excluant 0 (tous jeux) | « aucun lien » modéré ou fort **et** IC excluant 0 | « aucun lien » en tête et IC excluant 0 | idem, IC sous l'a priori du test | P(aucun lien) maximale |
|---|---|---|---|---|---|---|---|---|---|
| 0,2 | 30 | 67,1 | 0,0 | 0 | 21,6 | **0** | 0,0 | 0 | 0,714 |
| 0,2 | 50 | 58,7 | 1,7 | 0 | 30,8 | **0** | 0,0 | 0 | 0,759 |
| 0,2 | 100 | 42,5 | 5,1 | 0 | 51,6 | **0** | 0,2 | 0 | 0,813 |
| 0,3 | 30 | 46,3 | 0,0 | 0 | 39,6 | **0** | 0,0 | 0 | 0,714 |
| 0,3 | 50 | 30,1 | 0,2 | 0 | 61,5 | **0** | 0,0 | 0 | 0,759 |
| 0,3 | 100 | 10,2 | 0,1 | 0 | 87,1 | **0** | 0,0 | 0 | 0,806 |

**Lecture.**
- **Jamais de « soutien fort » de l'absence de lien** jusqu'à 100 paires : avec kappa = 1/3, P(aucun lien) ne dépasse pas 0,81. Le facteur de Bayes en faveur de l'absence est plafonné à ces tailles (§4).
- **0 jeu sur 6000** avec « aucun lien » soutenu (modéré ou fort) et un intervalle excluant 0 ; le seuil de 5 % fixé par la session parente n'est pas approché. **Aucune ligne n'est ajoutée au bilan.**
- « Aucun lien » en tête avec un intervalle excluant 0 n'arrive que dans 2 jeux sur 1000 (rho = 0,2, n = 100), toujours sans soutien : le bilan y dit « indécis » et l'estimation y est donnée « si un lien existe ».
- À surveiller : à rho = 0,2 et n = 100, « aucun lien » en tête avec un soutien modéré dans 5,1 % des jeux (« plutôt pas de lien, sans conclure », jamais « très probablement ») ; l'intervalle contient alors 0, les deux lectures s'accordent. C'est le prix, déjà mesuré en §4, d'un a priori qui attend des liens d'ordre moyen.
