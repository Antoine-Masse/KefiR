# Relecture adverse du 9 novembre 2026 : m.bay à plusieurs facteurs

Relecteur : agent qui a lu le code (`R/sys_bay_factorial.R`, l'aiguillage de `R/m.bay.R`, `tests/testthat/test-bay-factorial.R`) mais aucun document de `dev/mbay/`. Scripts et sorties : `relecture/agent_rev4/` (rev4_a à rev4_e). Tout ce qui suit a été reproduit par exécution chez lui ; les corrections sont de moi, testées dans `test-bay-factorial.R` (bloc « relecture du 9 novembre »).

**Rien de bloquant.** La loi HOP a été recalculée à la main à 2 et 3 facteurs (écart 7e-18). Sur les 18 modèles héréditaires d'un plan à 3 facteurs, les log BF concordent avec `anovaBF(whichModels = "all")` (identiques à un terme, ≤ 0,08 sinon), et les probabilités par terme à 0,002 près. Écritures de formule, cases vides, une observation par case, déséquilibre, facteur à 3 niveaux, NA, graine et état global : sans défaut.

| # | Défaut | Gravité | Suite |
|---|---|---|---|
| 1 | Noms entre backticks (`` `my A` ``) : erreur R brute | important | **corrigé** : backticks retirés des termes |
| 2 | Facteurs confondus ou emboîtés (`A + C`, C = A renommé ; `A + N`, N emboîté) : chiffres trompeurs, sans avis | important | **corrigé** : refus, rang de la matrice additive ; plans emboîtés annoncés comme non traités |
| 3 | Précision annoncée non fiable (4 à 83 % sur un plan, probabilités stables à 0,02 ; `draws` plus grand n'aide pas ; colonnes de sensibilité ignorées) | important | **corrigé** : l'erreur de BayesFactor n'est plus utilisée ; le calcul central est refait avec une seconde graine, et l'écart maximal des probabilités est annoncé (`error_pct`, en points) ; avis au-delà de 2 points, sans conseil sur `draws` |
| 4 | « Écart-type au sein des cases » = résidu du modèle additif (1,603 contre 1,130 en croisement pur) | important | **corrigé** : écart-type autour des moyennes de case ; résiduel, et dit, s'il n'y a qu'une observation par case |
| 5 | `draws` < 10 000 sans effet | mineur | **documenté** dans `?m.bay` : au moins 10 000 itérations à plusieurs facteurs |
| 6 | Résultat central dépendant de `sensitivity` | mineur | **corrigé** : une graine par valeur de r, r central avec la graine de l'appel |
| 7 | « Interactions moins probables » dit sans interaction | mineur | **corrigé** |
| 8 | Instabilité de la réponse globale non dite | mineur | **corrigé** : ligne « l'ordre peut s'inverser » sous les hypothèses |
| 9 | 2⁴ : 279 s sans indication | mineur | **atténué** : message au-delà de 50 modèles ; la durée reste |
| 10 | `x, g` avec un data.frame : sans interaction, sans le dire | mineur | **corrigé** : avis « aucune interaction n'est testée » |

**Non vérifié par le relecteur** : mémoire à 4 facteurs à nombreux niveaux (3×3×3×3).
