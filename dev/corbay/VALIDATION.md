# Validation de cor.bay

R 4.6.1, BayesFactor 0.9.12-4.8, bayestestR 0.20.0, Windows 11, 11 novembre 2026. Scripts dans `scripts/`, sorties dans `scripts/resultats/` (`.rds`, `.log`). `proto.R` recopie les calculs du paquet sans le récit (pour aller vite) ; `sim_auto.R` appelle le code du paquet lui-même. Tous les taux sont des pourcentages de jeux simulés.

Lecture : « P > 0,95 » est le seuil d'un « soutien fort » et « P > 0,75 » celui d'un « soutien modéré » (même échelle que m.bay). Ce ne sont pas des seuils de p-value ; la comparaison avec `cor.test` à p < 0,05 est donnée pour situer, pas pour égaler.

## 1. Accord numérique (`check_agree.R`, 200 jeux, n de 5 à 200, rho de -0,8 à 0,8)

| Calcul | Contre | Écart sur log BF |
|---|---|---|
| BF bilatéral = moyenne des deux moitiés de `correlationBF(nullInterval = c(0, 1))` | `correlationBF` bilatéral direct | médiane 0,0005, max 0,016 (erreur d'intégration propre à BayesFactor) |
| BF de Kendall (`integrate`, intervalle coupé autour du sommet) | méthode des trapèzes, 400 001 points sur [-1, 1] | 7e-12 |
| Probabilités des trois hypothèses | formule fermée à partir des deux moitiés | 1e-10 (test `test-cor.bay.R`) |
| A priori sur tau, kappa = 1 | pi/4 cos(pi tau / 2), éq. 10 de van Doorn et al. 2018 | 1e-12 (test) |

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
| échelle 1 à 5 | 0,2 | 20 | 6,7 | 1,3 | 13,0 | 3,0 | 1,3 | 13,0 | 10,7 | 10,3 |
| échelle 1 à 5 | 0,2 | 50 | 14,7 | 2,3 | 25,0 | 7,7 | 2,3 | 20,7 | 20,7 | 21,0 |
| échelle 1 à 5 | 0,2 | 100 | 28,3 | 10,7 | 44,3 | 21,3 | 9,3 | 44,0 | 41,7 | 41,7 |
| échelle 1 à 5 | 0,3 | 20 | 14,3 | 0,7 | 24,3 | 6,0 | 1,3 | 23,0 | 22,7 | 22,7 |
| échelle 1 à 5 | 0,3 | 50 | 36,0 | 10,3 | 47,3 | 24,7 | 10,3 | 45,3 | 44,0 | 43,0 |
| échelle 1 à 5 | 0,3 | 100 | 65,7 | 38,0 | 78,0 | 57,0 | 34,7 | 77,0 | 75,7 | 75,3 |
| échelle 1 à 5 | 0,5 | 20 | 44,7 | 14,7 | 57,0 | 25,3 | 14,7 | 54,3 | 52,0 | 52,7 |
| échelle 1 à 5 | 0,5 | 50 | 85,7 | 67,0 | 93,0 | 80,7 | 64,3 | 92,3 | 91,7 | 91,3 |
| échelle 1 à 5 | 0,5 | 100 | 98,7 | 94,7 | 99,3 | 97,3 | 97,3 | 99,7 | 99,3 | 99,3 |

**Lecture.**
- Sur des données normales, « soutien modéré » (P > 0,75) se déclare un peu moins souvent que p < 0,05 (5 à 14 points de moins), « soutien fort » (P > 0,95) nettement moins : c'est un seuil plus exigeant, et l'a priori ramène les petits effets vers 0. Pearson bat Kendall de 0 à 14 points, comme attendu ([vanDoorn2018], fig. 4).
- Lien monotone courbe : Kendall double ou triple la puissance de Pearson (rho = 0,5, n = 50 : 66 % contre 19 % à P > 0,95). La règle `auto` y prend Kendall (valeurs lointaines de exp(2 z)) : voir §6.
- Échelle 1 à 5 : Kendall (tau-b) est le plus puissant des trois méthodes bayésiennes, et à P > 0,75 il dépasse `cor.test` de Pearson ; Spearman sur rangs perd le plus (ex-aequo). La règle prend Kendall.
- t3 bivariée : Pearson « trouve » plus souvent un lien faible que Kendall, mais sa calibration sous H0 sur ces données est à lire en §2 avant de s'en réjouir.

## 2. Sous H0 : faux « oui » (`sim_h0.R`, 500 jeux par case)

(en cours)

## 6. La règle `auto` du paquet (`sim_auto.R`, 300 jeux par case)

(en cours)
