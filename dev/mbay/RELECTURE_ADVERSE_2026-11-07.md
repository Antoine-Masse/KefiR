# Relecture adverse du 7 novembre 2026 : m.bay à un facteur

Objet : `R/m.bay.R`, `R/sys_bay_input.R`, `R/sys_bay_laws.R`, `R/sys_bay_bilan.R`, `R/sys_bay_oneway.R`, tests associés, avant correction (commit `a315c6d`). Relecteur : agent, qui a reçu le code et non mes conclusions. Scripts et sorties : `dev/mbay/scripts/relecture/agent_rev3_*`. Rapport transcrit et résumé par le coordinateur. **État** = ce que j'ai fait du point ; une correction n'est dite faite que si un test l'exerce.

## Verdict du relecteur

Constantes et mesures justes entre les modèles (k=3 contre intégration brute, écart ≤ 2e-3, `agent_rev3_03`). Quatre défauts bloquants.

## BLOQUANT

| # | Défaut | Exemple exécuté | État |
|---|---|---|---|
| B1 | Un groupe très éloigné ne donne plus de preuve dès k ≥ n_s + 2u : BF ∝ max(constante, d^(n_s+2u−k)) | k=20, n_s=3, d=10, Welch p=6e-8, zs log BF 49,7 : **« Aucun effet » en tête, 83,4 %** ; k=8, d=1000 : P=0,647 | **Corrigé** : H1 est un mélange à parts égales de « tous les groupes » et « un seul groupe s'écarte ». Pente d^(n_s+2u−2), quel que soit k [main]. Test `B1 relecture`. L'échelle reste un réglage de KefiR |
| B2 | Groupe sans variation, ex-aequo massifs : plantage ou résultat tiré au sort. L'intégrale **diverge** si (k−1)u ≤ (n_i−1)/2 | `m.bay(c(3,3,3,3,3,1,2,3,4,5), rep(1:2,each=5))` : erreur ; Likert avec groupe constant : P=1,0000 selon la graine, alors que le gaussien donne log BF −0,28 | **Corrigé** : ces lois sont écartées avec un avis (`.bay_laws_allowed`), un poids non fini ou un ESS < 30 rend la loi caduque. **Limite** : seuil simple (valeur répétée sur plus de la moitié d'un groupe ou de l'échantillon) ; la voie discrète reste à faire |
| B3 | `.bay_law_kernel_grid` soustrait deux nombres presque égaux : écart 0,145 à c₂=1e14, Inf au-delà | `agent_rev3_14` | **Corrigé** par le centrage du relecteur ; test `B3 relecture` |
| B4 | Mauvaise analyse sans avis : `wt=` sans `id` traité en groupes indépendants ; `I(dose>2)` analysé sur la variable brute ; `g` plus court recyclé | `agent_rev3_16` | **Corrigé** : appariement par l'ordre comme m.test, refus d'un terme qui n'est pas une variable, refus de la longueur ; test `B4 relecture`. **Limite** : `log(x)` en prédicteur est refusé, il faut créer la colonne |

Le relecteur demandait d'arbitrer B1 auprès d'Antoine. Consigne d'autonomie : tranché ici, voir `DECISIONS.md` D6.

## IMPORTANT

| # | Défaut | État |
|---|---|---|
| I1 | Verdict de sensibilité dicté par une loi imposée que les données écartent (P(loi) 1e-24 à 1e-40) : « Les données ne tranchent pas » sous « > 99,9 % » | **Corrigé** : lois imposées sur une ligne à part, hors de la plage ; test |
| I2 | Cotes annoncées fausses si `prior_h1 ≠ 1/2` | **Corrigé** ; test |
| I3 | Borne droite de la grille en g : 3,7 d'écart sans signal, rapport extrême | **Corrigé** ; test `I3 relecture` |
| I4 | k=20, variances très différentes : ESS de 2 à 30, erreur sous-estimée | **En partie** : ESS < 30 écarte la loi, plafond de 40 000 tirages. Non refait à k=20 |
| I5 | Le bilan ne dit rien de la taille de l'effet | **En partie** : moyennes et plus grand écart, descriptifs. Taille d'effet et intervalle de crédibilité : jalon post-hoc |
| I6 | « Les données ne tranchent pas » sous H0 dans 78 % des jeux à N=60 | **Gardé, décision** : avec la cote 1:4 dans la grille, tout BF sous 4 fait changer l'hypothèse en tête, et c'est vrai : à ce N, l'absence d'effet n'est pas établie. À chiffrer (`agent_41` C) |

## MINEUR

Corrigés : `draws`, `seed`, `rscale`, `u` validés ; `y*1e±200` ; variable hors de `data` ; `id` composé ; `g` liste ou matrice ; `g` infini ; message « binaire » sur des données constantes ; articles et décimales des messages ; « chacune des 1 formes », lois annoncées seulement si elles sont retenues ; « proportions de précision » ; pilote de la voie t écarté de l'estimation.

Corrigés aussi : références de `sys_bay_oneway.R` (Liang, Trefethen, DOI contrôlés sur Crossref).

Non corrigés : 0,1+0,2 et 0,3 fusionnés en silence dans `.bay_oneway_prep` ; temps de calcul aux défauts à k ≥ 8 (35 s sans la seconde forme, à remesurer) ; `y~.` (« dans les données : .. ») ; `return=FALSE` non rappelé dans l'exemple.

## Non fait par le relecteur

Référence `normal_het` à k=3 (trop lent) ; marginale t complète à k ≥ 3 ; invariances y → −y et étiquettes (script écrit, non abouti) ; balayage de petits plans sous t, interrompu : **on ne sait pas si un plantage touche des données sans ex-aequo** ; citations non revérifiées. Sa relecture est antérieure au mélange de B1, **qui n'a donc pas été relu**.

## Seconde relecture, le même jour : le mélange et les corrections

Par un second agent, sur le commit `c349c30`. Scripts : `relecture/agent_rev4_*`. Constantes justes (formes contre ZS à 1e-15, mélange recalculé à la main). Aucun bloquant établi, un candidat instruit.

| # | Défaut | État |
|---|---|---|
| I1 | Le seuil ESS 30 rendait l'ensemble de lois dépendant de la graine : log BF 14,3 ou 6,2 sur un même jeu sans effet, sd de 0,01 à 100, P=0,998 | **Corrigé pour la graine** (pilote itératif : 4 lois, valeurs stables à 4000 et 16 000 tirages). **Instruit** : `agent_44`, 100 jeux, P>0,95 dans 2 % : ce jeu est un cas défavorable et non un biais systématique, à faible résolution. À 1500 tirages, t_het reste souvent écartée à ces variances |
| I2 | À k=20, les lois hétérogènes sont souvent écartées, donc le test de B1 y est faible | **Non corrigé.** Le test tient à k=8 et 12. À k=20, l'hypothèse « ce groupe est plus dispersé » peut manquer à la comparaison ; l'avis le dit |
| I3 | Le seuil d'ex-aequo écartait normal_het alors que son intégrale est finie | **Corrigé** : règles par divergence, dérivées à la main, avec marge. Non vérifiées par un tiers |
| I4 | La coupe à 5 groupes pour t, au-delà de k=8, perdait un vrai signal (log BF −9,6) | **Corrigé** : pas de coupe jusqu'à k=20. Au-delà, coupe sur les médianes, avec le risque inchangé |
| I5 | `g` data.frame recyclé ; `paired=1` ignoré | **Corrigé** |
| I6 | « comparée et non choisie » pour une seule loi, Dirichlet annoncé sans objet | **Corrigé** |

Mineurs corrigés : erreurs R brutes des entrées (attrapées et reformulées), `y` infini, `between` absent, `seed` maximal, `.bay_law_normal` sans calcul inutile, `dcauchy` sans avertissement, accord. Non corrigé : `h1_form="zs"` peut écarter une loi sans le dire ; lois écartées à r/2, 2r ou u/3, 3u non annoncées.

**Pas relu par un tiers** : verdict à 1:2–2:1, avis sur les queues lourdes, pilote itératif, règles d'ex-aequo finales.
