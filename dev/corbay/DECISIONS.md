# Décisions de cor.bay()

Pendant bayésien de `stats::cor.test()` pour deux variables quantitatives. Mêmes principes que `m.bay` (dev/mbay/DECISIONS.md) : a priori fixés d'avance par la fonction, jamais réglés sur les données, avec sensibilité ; paquets R établis plutôt qu'un moteur maison ; bilan décisionnel lisible étape par étape, sans p-value. Chiffres des simulations : `VALIDATION.md` ; scripts : `scripts/`.

## Bibliographie vérifiée

Vérifiée le 11 novembre 2026 (notice de l'éditeur ou du dépôt institutionnel ; van Doorn et al. 2018 lu en entier sur le PDF de l'UvA).

- **[Ly2016]** Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology. *Journal of Mathematical Psychology*, 72, 19-32. https://doi.org/10.1016/j.jmp.2015.06.004 — test de corrélation de Jeffreys, a priori bêta étiré ; c'est le modèle `Jeffreys-beta*` de `BayesFactor::correlationBF`.
- **[vanDoorn2018]** van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Bayesian inference for Kendall's rank correlation coefficient. *The American Statistician*, 72(4), 303-308. https://doi.org/10.1080/00031305.2016.1264998 — éq. 6 (T*), 8 (relation de Greiner), 9 (a priori sur tau), 11 (vraisemblance normale de moyenne 3/2 tau racine(n), écart-type 1), 13 (facteur de Bayes) ; §2.4 : approximation normale suffisante sous H0 dès n = 10 d'après Ferguson, Genest et Hallin (2000), conservatrice sous H1.
- **[vanDoorn2020]** van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (2020). Bayesian rank-based hypothesis testing for the rank sum test, the signed rank test, and Spearman's ρ. *Journal of Applied Statistics*, 47(16), 2984-3006. https://doi.org/10.1080/02664763.2019.1709053 — Spearman exact par augmentation de données (échantillonneur de Gibbs) ; code dans JASP, pas dans un paquet CRAN.
- **[Makowski2020]** Makowski, D., Ben-Shachar, M. S., Patil, I., & Lüdecke, D. (2020). Methods and algorithms for correlation analysis in R. *Journal of Open Source Software*, 5(51), 2306. https://doi.org/10.21105/joss.02306 — paquet `correlation` : Spearman bayésien = `correlationBF` sur les rangs ; mêmes noms d'a priori que BayesFactor (medium = 1/3).
- **[Makowski2019]** Makowski, D., Ben-Shachar, M. S., & Lüdecke, D. (2019). bayestestR. *JOSS*, 4(40), 1541. https://doi.org/10.21105/joss.01541 — `hdi()`.
- **[KassRaftery1995]** Kass, R. E., & Raftery, A. E. (1995). Bayes factors. *JASA*, 90(430), 773-795. https://doi.org/10.1080/01621459.1995.10476572 — seuils 3, 20, 150, repris de m.bay.
- **[Ferguson2000]** Ferguson, T. S., Genest, C., & Hallin, M. (2000). Kendall's tau for serial dependence. *Canadian Journal of Statistics*, 28, 587-604 — cité par [vanDoorn2018] pour n ≥ 10 ; non relu directement.
- **[CookWeisberg1982]** Cook, R. D., & Weisberg, S. (1982). *Residuals and Influence in Regression*. New York : Chapman and Hall — distance de Cook, seuil usuel de 1 pour un point qui à lui seul déplace l'ajustement.

Paquets : BayesFactor 0.9.12-4.8, bayestestR 0.20.0, MASS, R 4.6.1.

## C1. Pearson : `BayesFactor::correlationBF` : DÉCIDÉ

Le test de [Ly2016] est celui de `correlationBF` (modèle « Jeffreys-beta* ») : on l'appelle tel quel. Pour avoir les trois hypothèses (positif, négatif, aucun), un seul appel avec `nullInterval = c(0, 1)` rend les deux moitiés, BF(rho > 0 contre 0) et BF(rho < 0 contre 0). Le BF bilatéral est leur moyenne (chaque moitié porte la moitié de la masse a priori). Écart avec l'appel bilatéral direct : voir VALIDATION §1.

## C2. A priori : DÉCIDÉ

- **Sur rho** : bêta étiré, (rho + 1)/2 ~ Bêta(1/kappa, 1/kappa), **kappa = 1/3**, le « medium » de BayesFactor et du paquet `correlation`. Il met 20,7 % de chances a priori sur |rho| > 0,5. JASP prend kappa = 1 (uniforme, [Ly2016]) : 50 % sur |rho| > 0,5, ce qui pénalise les liens modérés, fréquents en TP. Raison du choix : cohérence avec m.bay, qui prend les « medium » de BayesFactor (D5). JASP reste dans la sensibilité indirectement : kappa = 2/3 est la borne haute de la grille ; écart mesuré à kappa = 1 en VALIDATION §4.
- **Sur tau (Kendall)** : le même a priori sur rho, reporté par la relation de Greiner tau = 2/pi arcsin(rho) ([vanDoorn2018], éq. 9). Les deux méthodes partent donc du même point.
- **Sur les hypothèses** : P(aucun lien) = 1/2 ; le lien, 1/2, partagé à parts égales entre ses deux signes (1/4 chacun). C'est exactement l'hypothèse bilatérale de [Ly2016] découpée par signe : la probabilité d'un lien est la même qu'avec deux hypothèses.
- **Sensibilité** : kappa/2, kappa, 2 kappa, croisés avec des cotes de départ 1:2, 1:1, 2:1 (grille de m.bay, D4). Si l'hypothèse en tête change quelque part dans la grille, le bilan dit « non tranché ».
- `prior_h1` et `kappa` sont des arguments, comme `rscale` dans m.bay, documentés « à fixer d'avance » ; rien n'est estimé sur les données.

## C3. Choix automatique de la méthode : DÉCIDÉ (règle descriptive, sans p-value)

Consigne : « le choix se fait automatiquement, comme m.test choisit son test ». m.bay a écarté l'aiguillage par test préliminaire (D6) au profit d'une moyenne de modèles. Pour deux variables, la moyenne de modèles (normale bivariée contre t bivariée) demanderait la vraisemblance marginale d'une t bivariée : aucun paquet établi ne la donne sans Stan, indisponible ici. **On garde donc un aiguillage, mais par une règle descriptive annoncée, pas par une p-value**, et l'autre méthode est toujours calculée en contre-épreuve (C9).

**Règle** : tau de Kendall si au moins une condition, Pearson sinon :
1. une valeur au-delà de 3 écarts interquartiles des quartiles, sur x ou sur y (Tukey, valeur « lointaine » : probabilité 2e-6 par valeur sous la loi normale) ;
2. une paire influente : distance de Cook au-delà de 1 ([CookWeisberg1982]) dans la régression de y sur x ou de x sur y ; cherchée à partir de 20 paires ;
3. une échelle courte : 7 valeurs distinctes ou moins sur une variable (notes, Likert), où Pearson suppose une continuité absente et où les ex-aequo dominent.

**Changement mesuré (11 novembre).** La première version cherchait les « points hors du nuage » par une distance robuste (MCD de `MASS::cov.rob`, seuil khi-deux 0,999). Mesure : sur des données **normales**, ce critère se déclenchait dans 47 à 65 % des jeux (n = 10 à 100), et 40 à 55 % même avec repondération ou `robustbase::covMcd` (`scripts/check_mcd.R`, `check_mcd2.R`) : la règle envoyait vers Kendall la moitié des jeux sans raison, au prix de la puissance. La distance de Cook > 1 se déclenche dans 4,3 % des jeux normaux à n = 20, 1,3 % à n = 30, 0,3 % à n = 50, 0 à n = 100, et attrape 95 à 100 % des points plantés à (3, -3) ou (6, 6). Sous 20 paires (29 % de déclenchements à n = 10), elle n'est pas cherchée : seules les valeurs lointaines (règle 1) comptent.

Pourquoi Kendall plutôt que Spearman comme voie robuste : Kendall a un test bayésien publié avec sa propre vraisemblance ([vanDoorn2018]) ; le Spearman sur rangs est une approximation (C4). Mesures de la règle : C6 et VALIDATION §2 à §3.

## C4. Rangs : Kendall de van Doorn et al. (2018), Spearman sur rangs à la demande : DÉCIDÉ

- **Kendall** : aucun paquet CRAN ne l'exporte (le code est dans JASP et sur OSF). La méthode tient en une intégrale à une dimension (éq. 13), programmée en 30 lignes dans `sys_bay_cor.R` avec `stats::integrate`, à partir du tau-b de `stats::cor`. C'est la seule part « maison », assumée : elle recopie une formule publiée, et elle est vérifiée par un second calcul indépendant (méthode des trapèzes, test `test-cor.bay.R`) et par la calibration sous l'a priori (VALIDATION §3). Pour n grand, l'intégrande est une bosse étroite : l'intervalle est coupé autour de son sommet (test à n = 4000).
- **Ex-aequo : CORRIGÉ après mesure (11 novembre).** [vanDoorn2018] ne traite pas les ex-aequo. La première version prenait le tau-b et la variance « sans ex-aequo » de l'éq. 6 : sur une échelle de 1 à 5, T* avait alors un écart-type de **1,19 à 1,25 sous H0** au lieu de 1 (`scripts/check_ties.R`), et le bilan était trop sûr de lui (annoncé 0,70, fréquence réelle 0,52 à n = 100 ; VALIDATION §3). Correction : S = tau-b racine((n0 - n1)(n0 - n2)) est standardisé par sa variance sous H0 corrigée des ex-aequo (formule de Kendall, la même que `cor.test(method = "kendall", exact = FALSE)` : T* y est égal à sa statistique z, vérifié par un test), et la pente 3/2 racine(n) de l'éq. 11 est multipliée par le même rapport, pour que le paramètre reste le tau-b. **Sans ex-aequo, T* et la pente sont exactement ceux de l'article** (vérifié par un test). Après correction : écart-type 0,995 à 1,014 sous H0, calibration retrouvée (VALIDATION §3).
- **Spearman** : le test exact [vanDoorn2020] demande un échantillonneur de Gibbs à augmentation de données, sans paquet établi : on ne le reprogramme pas. `method = "spearman"` donne Pearson bayésien sur les rangs, l'approche du paquet `correlation` [Makowski2020] ; mesuré en VALIDATION §2-3. La règle `auto` ne le choisit jamais.

## C5. Hypothèse unilatérale : DÉCIDÉ

`alternative = "greater"` : deux hypothèses, lien positif (a priori `prior_h1`, sur rho > 0) et aucun lien (rho = 0) ; c'est le BF+0 de [Ly2016] et de JASP. `"less"` de même. Le bilan ne parle pas du signe non testé. Si le lien observé va dans l'autre sens, une ligne le dit : « il sort des hypothèses posées » ; le facteur de Bayes favorise alors « aucun lien », ce qui est la bonne lecture d'un test unilatéral mais pourrait tromper sans cette ligne. L'estimation est faite sous l'hypothèse testée (`nullInterval` pour Pearson, grille tronquée pour Kendall).

## C6. Queues lourdes et valeurs éloignées : MESURÉ

Mesures par méthode (VALIDATION §2, 500 jeux par case) :
- **Une paire aberrante sous H0** (normales indépendantes, une paire à (6, 6)) : Pearson bayésien dit « lien » à P > 0,95 dans 19 à 58 % des jeux (`cor.test` : 82 à 98 % à p < 0,05) ; Kendall dans 0,4 % au plus. C'est le cas que la règle 1 (valeur lointaine) attrape.
- **t3 bivariée sous H0** (x et y non corrélés mais dépendants par une échelle commune) : Pearson à P > 0,95 dans 3 à 11 % des jeux, de plus en plus quand n croît ; Kendall 0,4 % au plus. Les marges t3 ont souvent des valeurs au-delà de 3 écarts interquartiles, donc la règle prend Kendall dans la plupart de ces jeux ; ce qu'il en reste est mesuré par `sim_auto.R` (VALIDATION §6).
- **t3 indépendantes, log-normales** : toutes les méthodes restent sous 1,6 % à P > 0,95 ; Pearson n'y est pas trompé (sous indépendance, sa loi d'échantillonnage est robuste).
- **Puissance** (§5) : Kendall perd 0 à 14 points sur des données normales, en gagne 20 à 50 sur un lien monotone courbe, et sur la t3 bivariée il passe devant Pearson dès rho = 0,5.

Décision : la voie robuste est Kendall, choisie par la règle C3. **Pas de modèle t bivarié** : sa vraisemblance marginale n'existe dans aucun paquet établi utilisable sans Stan. Ce qui échappe à la règle (queues lourdes sans valeur au-delà de 3 écarts interquartiles ni paire influente) garde Pearson ; le taux de faux « oui » qui en résulte est mesuré en §6 plutôt que supposé.

## C7. Taille de l'échantillon : DÉCIDÉ

- Moins de 3 paires : refus (comme `cor.test` et `correlationBF`).
- Moins de 10 paires : une décision le dit à l'étape 2 (« la réponse s'appuiera sur le point de départ plus que sur les données ») ; les paires influentes ne sont pas cherchées sous 20 paires (C3). Kendall reste permis : sous H0, l'approximation normale est donnée suffisante dès n = 10 ([vanDoorn2018], §2.4) ; mesuré à n = 10 en VALIDATION §2.
- Pas de borne haute : `cor(method = "kendall")` est en n² ; 2000 paires prennent 3 s.

## C8. Estimation : DÉCIDÉ, a priori d'estimation CHANGÉ après mesure (11 novembre)

Sous l'hypothèse d'un lien (du signe testé en unilatéral), comme l'écart de m.bay : médiane a posteriori et intervalle HDI à `conf.level` (`bayestestR::hdi`).
- **A priori d'estimation** : bêta étirée alpha = beta = 1/2 (kappa = 2), le choix de Jeffreys pour l'estimation, que [vanDoorn2018] (§2.2) recommande « quand l'objet est l'estimation plutôt que le test » ; sur tau, il est uniforme. **Le test garde kappa = 1/3** (C2). Raison mesurée : avec kappa = 1/3 aussi pour l'estimation, un lien fort était trop ramené vers 0. Exemple : notes 1 à 5, n = 40, tau-b observé 0,79, r = 0,89. Avec kappa = 1/3, tau estimé 0,59 [0,40 ; 0,77] : l'intervalle excluait la valeur observée. Avec kappa = 2, 0,75 [0,54 ; 0,98]. Sur un lien modéré, l'écart est petit : r = 0,53, n = 30, Pearson 0,45 [0,19 ; 0,72] avec kappa = 1/3, 0,51 [0,23 ; 0,75] avec kappa = 2. Prix : l'estimation et le test ne partent pas du même a priori ; c'est dit en note à l'étape 6 et dans le détail.
- **Revu le 11 novembre au soir, à la demande de la session parente** (règle de famille : l'estimation vient du même modèle que le test, sauf argument fort chiffré). Mesure `scripts/check_estim.R` : 60 jeux par case, rho vrai fixé, intervalle à 95 %. On compte la couverture du coefficient vrai (tau vrai = 2/pi arcsin(rho)) et la part des jeux où l'intervalle contient la valeur observée.

| Méthode | rho vrai | n | Médiane, a priori du test (k = 1/3) | Médiane, Jeffreys (k = 2) | Couverture, k = 1/3 | Couverture, k = 2 | Observé dans l'intervalle, k = 1/3 | Idem, k = 2 |
|---|---|---|---|---|---|---|---|---|
| Kendall | 0,9 (tau 0,71) | 8 | 0,32 | 0,57 | 0 % | 100 % | 17 % | 98 % |
| Kendall | 0,9 | 15 | 0,43 | 0,63 | 13 % | 100 % | 22 % | 100 % |
| Kendall | 0,9 | 30 | 0,54 | 0,67 | 68 % | 98 % | 72 % | 100 % |
| Pearson | 0,9 | 8 | 0,66 | 0,84 | 57 % | 98 % | 88 % | 100 % |
| Pearson | 0,9 | 15 | 0,79 | 0,87 | 75 % | 95 % | 100 % | 100 % |
| Pearson | 0,9 | 30 | 0,85 | 0,88 | 87 % | 97 % | 100 % | 100 % |
| Kendall | 0,6 (tau 0,41) | 8 | 0,21 | 0,36 | 92 % | 98 % | 72 % | 100 % |
| Pearson | 0,6 | 8 | 0,38 | 0,55 | 88 % | 97 % | 93 % | 100 % |
| les deux | 0,3 | 8 à 30 | | | 90 à 95 % | 90 à 95 % | 97 à 100 % | 100 % |

  Avec l'a priori du test, un lien fort est rendu à la moitié de sa valeur (tau 0,32 pour 0,71 à n = 8), et l'intervalle de Kendall manque le vrai tau dans 87 à 100 % des jeux sous 15 paires. Il exclut même le tau observé dans 78 à 83 % des jeux. Or les liens forts sont courants en TP (étalonnage, mesures physiques). La bêta(3, 3) les juge très improbables : P(|rho| > 0,9) = 0,23 % a priori. **Décision : on garde Jeffreys pour l'estimation du bilan**, faute de quoi l'étudiant lirait un rho estimé démenti par son propre nuage de points. L'estimation sous l'a priori du test est gardée dans l'objet (`$estimate_test`) et dite par `print(x, detail = TRUE)`. Pour les liens faibles et modérés, les deux s'accordent. **Pour revenir à la règle de famille, une seule ligne** : `.corbay_kappa_est()` dans `sys_bay_cor.R`.
- Pearson et Spearman : 10 000 tirages de `correlationBF(posterior = TRUE)` (Metropolis-Hastings à proposition indépendante de BayesFactor), graine propre à l'appel (`.bay_with_seed`), état aléatoire global intact.
- Kendall : densité a posteriori de tau sur une grille de 4001 points autour du sommet (éq. 11 × éq. 9), puis 10 000 tirages pour `hdi`. L'écart-type de T* est pris égal à 1 ; il est plus petit quand |tau| est grand (éq. 14), donc l'intervalle de Kendall est un peu trop large pour un lien fort : prudent, comme l'article le dit de son test.
- Quand « aucun lien » est en tête, l'estimation est dite « si un lien existe ».
- Pas de moyenne entre hypothèses : P(aucun lien) est donnée à part, comme dans m.bay.

## C9. Contre-épreuve et désaccord : DÉCIDÉ

Avec `method = "auto"`, l'autre méthode (Kendall si Pearson est choisi, et inversement) est calculée et dite à l'étape 5. **Désaccord** : les deux méthodes mettent des hypothèses différentes en tête, chacune à 75 % au moins. Le bilan dit alors « indécis, la réponse dépend de la méthode », comme m.bay dit « la réponse dépend de la forme du bruit ». Un écart sans netteté des deux côtés n'est qu'une ligne : Kendall est moins puissant que Pearson sur des données normales ([vanDoorn2018], fig. 4), un désaccord faible est attendu.

**Réglage après mesure (11 novembre).** Quand Kendall a été choisi *à cause de* valeurs lointaines ou de paires influentes, la réponse de Pearson est précisément celle que ces paires peuvent porter : son désaccord est la raison du choix, pas une contradiction. Il est dit (« la réponse de Pearson peut être portée par les paires à part ; elle n'est pas retenue ») et ne rend pas le bilan indécis. Sans ce réglage, 37 à 58 % des jeux à paire aberrante opposée au lien finissaient « indécis » alors que Kendall avait la bonne réponse. Après réglage : aucune contradiction comptée sur les 66 cases de VALIDATION §6 ; la contradiction reste possible quand Pearson est choisi, ou quand Kendall l'est pour une échelle courte seulement.

## C10. Vocabulaire commun avec m.bay : VÉRIFIÉ ET ALIGNÉ (11 novembre, demande de la session parente)

m.bay et cor.bay ont été lancés sur des données comparables (deux groupes de 15 / 30 paires ; un cas indécis, un cas net), en français et en anglais. **Déjà identiques** : l'échelle de soutien (« soutien faible / modéré / fort / très fort », seuils 0,6 / 0,75 / 0,95 / 0,993 de `.bay_support`) et la note du facteur de Bayes (« 1 à 3 faible, 3 à 20 modéré, 20 à 150 fort, au-delà très fort ») ; « non tranché » et « non tranché par les données » ; la ligne « Selon le point de départ : … ; l'ordre ne change pas / peut s'inverser » ; le verdict réduit à sa première phrase ; la forme de la phrase de compte rendu (« <hypothèse> (probabilité a posteriori x %), facteur de Bayes y (<fonction> {KefiR}, <méthode à citer>) »). **Alignés** :
- compte rendu indécis : cor.bay donnait le BF bilatéral brut (« facteur de Bayes 0,663 ») ; il donne maintenant, comme m.bay, le BF de l'hypothèse en tête contre la deuxième, avec « en faveur de l'absence de lien » quand « aucun lien » est en tête ;
- note finale : « probabilités d'hypothèses, pas des p-values ; ni la force, ni la cause d'un lien ; la méthode avec print(resultat, detail = TRUE) », la mention « ni la force, ni la cause » seulement quand un lien est retenu, comme « ni l'importance, ni la cause d'une différence » de m.bay ;
- « Observé : » devient « Coefficients observés : », sur le modèle de « Moyennes observées : ».

**Écart gardé** : l'estimation est dite « Si un lien existe, rho estimé : … » quand « aucun lien » est en tête ; m.bay dit « Écart estimé » sans cette réserve. La réserve est exacte (l'estimation est faite sous l'hypothèse d'un lien) ; elle pourrait être reprise par m.bay. Le test `test-cor.bay.R` (« vocabulaire de la famille ») fige ces formes.
