# Journal de m.bay

## Séance 19 : 11 novembre 2026, audit de b447387, sphéricité, réponse ordonnée à plusieurs facteurs

**Fait.**
- **Repli trop tôt** (audit de b447387), deux causes corrigées à la source dans `.bay_wrap` :
  - la protection « a = b » soudait « facteur : > 99,9 % » en un bloc insécable de 18 caractères ; une coupure est désormais possible après le deux-points ;
  - le retrait de continuation (15 colonnes) réduisait aussi la largeur de la première ligne (12 colonnes) ; chaque ligne a maintenant sa propre largeur.
  Le cas oui/non 2 × 2 est ajouté au test des 80 colonnes, avec une vérification qu'aucune ligne numérotée n'est repliée trop tôt.
- **Arrêt de Gibbs sur l'erreur Monte-Carlo : essayé, non retenu** (VALIDATION §22, mesures du 11 novembre). À 10 groupes, l'erreur vaut 0,037 à 250 tours, 0,019 à 1000 et 0,017 à 2000. S'arrêter à 1000 tours gagnerait environ 5 s sur 73 s. Viser 2 points demanderait environ 5000 tours, donc plus de temps. Les 73 s viennent surtout des 4 lois avec la sensibilité. La ligne d'attente (« environ 70 s ») est la réponse retenue.
- **Ce que couvre m.test** (vérifié dans le code) : ni McNemar ni Cochran, donc le oui/non apparié n'est pas prioritaire et reste refusé ; Mauchly et Greenhouse-Geisser en mesures répétées ; Scheirer-Ray-Hare pour les données non normales à plusieurs facteurs.
- **D27, sphéricité.** Mesurée (`agent_74`) : pas de faux « oui » fort sous AR(1) ; pas de correction, ε dit dans le récit.
- **D26.** Réponse ordonnée à plusieurs facteurs (`polr`, BIC, HOP) ; `agent_75`.

**Vérifié.** R CMD check sur d7a295a : Status OK, dans un dossier propre. Un check lancé en même temps dans le même dossier avait échoué (« path does not exist ») : collision entre les deux runs, pas défaut du code. Depuis, un check à la fois, chacun dans son dossier. Tests bay : 801.

**Suite proposée.** Variances inégales et queues lourdes avec sujets ou covariables ; réponses discrètes et ordonnées avec sujets.

## Séance 18 : 11 novembre 2026, temps au-delà de 8 groupes, post-hoc à plusieurs facteurs calibré, réponses discrètes à plusieurs facteurs (D25)

**Fait.**
- **Audit de 8c0523a.** Ligne d'attente avec durée estimée avant le calcul des lois (le temps vient surtout des 4 lois, pas de Gibbs) ; tours adaptés mesurés et non retenus ; PERIMETRE §7 relu en entier (9922725).
- **D20 calibré** (`agent_72`, VALIDATION §23).
- **D25.** Oui/non et comptages à plusieurs facteurs entre groupes (`R/sys_bay_discrete_fact.R`) ; `agent_73` (VALIDATION §24).

**Vérifié.** R CMD check sur 9922725 et sur b447387 : Status OK. Tests bay ; deux tests qui attendaient le refus des réponses discrètes à plusieurs facteurs sont mis à jour.

**Non fait.** Variances inégales et queues lourdes avec sujets ou covariables ; réponses discrètes avec sujets.

## Séance 17 : 11 novembre 2026, audit de aed2319 et a206ba5, calibration au-delà de 8 groupes

**Fait.**
- **Audit.** Écarts de l'étape 6 tirés du mélange de lois (normale et Student), avec un garde-fou testé ; P(variance par case) affichée elle-même ; soutien de l'hypothèse globale sous un « oui » ; écart-type à 3 chiffres significatifs (D24, addendum).
- **Choix de la suite** : calibrer le post-hoc au-delà de 8 groupes avant les queues lourdes avec sujets ou covariables. C'est une voie déjà livrée et jamais simulée ; une erreur y toucherait les utilisateurs dès aujourd'hui, alors que l'autre voie demande un moteur neuf. `agent_71` lancé : accord Gibbs-exact à 8 groupes, H0 à 10 et 12 groupes, deux paquets à 10 groupes.

**Vérifié.** Tests bay (779) ; R CMD check sur 974f709 : Status OK. `agent_71` (VALIDATION §22) : aucune fausse paire sous H0 à 10 et 12 groupes ; Gibbs à 2,6 points du calcul exact en médiane ; pas de changement de méthode.

**Suite proposée.** Queues lourdes avec sujets ou covariables (moteur neuf).

## Séance 16 : 11 novembre 2026, audit de 9c44239 et queues lourdes à plusieurs facteurs (D24)

**Fait.**
- **Audit de 9c44239** (aed2319) : bilan mené par les termes dès qu'un terme reçoit un « oui » (tous les plans à plusieurs termes) ; lecture par réponse en MANOVA à plusieurs facteurs ; écart de puissance avec Pillai mesuré (même rapport, effet sur une réponse ou réparti).
- **D24.** Loi de Student à plusieurs facteurs (`R/sys_bay_fact_t.R`), moyennée avec les deux lois normales de D19.

**Vérifié.** R CMD check sur aed2319 et sur a206ba5 : Status OK. `agent_70` avant et après (VALIDATION §21). `test-bay-fact-t.R` (mesure contre `lmBF` et contre la forme close). Tests bay.

**Non fait.** Student par case ; queues lourdes avec sujets, covariables, ou en MANOVA.

## Séance 15 : 11 novembre 2026, audit de 8ffb16a et MANOVA à un facteur (D22)

**Fait.**
- **Audit de 8ffb16a.** « 10 sujets x 3 répliques » au lieu de « n = 30 » ; écarts standardisés sur l'écart-type d'une mesure d'un sujet dans une condition (résiduel + sujet × condition), décidé dans D21 (d36a62a).
- **D22.** Trois moteurs comparés (BFpack, BIC, lmBF empilé). `R/sys_bay_manova.R` : récit en 7 étapes, lecture par réponse, bilan.

**Vérifié.** `agent_68` (VALIDATION §19) ; `test-bay-manova.R` contre `lmBF` sur les données empilées ; tests bay.

- R CMD check sur 3135136 (D22 et les deux corrections d'audit) : Status OK.
- **D23.** MANOVA à plusieurs facteurs : `R/sys_bay_manova_fact.R`, deux formes par terme, HOP ; `agent_69` (VALIDATION §20).

**Non fait.** MANOVA avec covariables ou sujets ; queues lourdes à plusieurs facteurs.

## Séance 14 : 10 novembre 2026, audit de 23c941e et mesures répétées avec répliques (D21)

**Fait.**
- **Audit de 23c941e.** Quatre retouches : mots du compte rendu repris du bilan, écarts par niveau quand une interaction est non tranchée, repli des listes numérotées, réponse dans le titre de l'étape 6 (5eb84c5).
- **D21.** Plan « mixed » accepté sans pentes aléatoires : sujet et sujet × condition aléatoires, dans le calcul comme dans le post-hoc ; récit « avec répliques ».

**Vérifié.** `agent_67` (VALIDATION §18) ; `test-bay-mixed.R` contre `lmBF` ; tous les tests bay.

- `R CMD check --no-manual --ignore-vignettes` sur 8ffb16a (suite complète comprise) : Status OK.

**Non fait.** Pentes aléatoires.

## Séance 13 : 10 novembre 2026, variances par case (D19) et post-hoc à plusieurs facteurs (D20)

**Fait.**
- **D19, correction.** `R/sys_bay_fact_laws.R` : loi normale à une variance par case, moyennée avec la variance commune dans chaque facteur de Bayes. La prudence provisoire ne sert plus qu'en cas d'échec de l'échantillonnage.
- **D20.** `.bay_posthoc_factorial` et `.bay_story_posthoc_fact` : une famille par terme probable du plus haut ordre ; le bilan donne les paires par famille.
- **Retouches.** « -0,00 » s'écrit « 0,00 » (`.bay_fx`) ; accord « à B, x égaux ».

**Vérifié.**
- `agent_66` : sous H0, 0,9 % de faux « oui » à P > 0,95 (7/800, contre 12,5 % avant) ; puissance inchangée à variances égales ; 1,8 s médian à N = 40 (VALIDATION §16).
- Tous les tests bay passent (dont `test-bay-fact-laws.R`, nouveau) ; deux tests qui attendaient l'absence de post-hoc sont mis à jour.

- `R CMD check --no-manual --ignore-vignettes` (suite complète comprise) : Status OK, une information de xref sur `ez`.
- Relecture débutant (`relecture/lecteur_recit/RELECTURE_D19_D20.md`) : cinq retouches faites (terme non tranché détaillé « à titre indicatif », facteur entre sujets, « seul ou en interaction » dans le compte rendu, ordre des cases, phrase de l'étape 3 selon le résultat, note « à X égal »).

**Non fait.** Lois t à plusieurs facteurs ; variances par case avec sujets ou covariables ; calibration par simulation du post-hoc D20.

## Séance 12 : 10 novembre 2026, réponse ordonnée et routage (D17, voies 3 et 4)

**Fait.**
- **Routage des entiers.** Une échelle bornée va en voie ordinale (au plus 10 valeurs) ou reste continue, et le récit le dit. Les comptages ne passent en voie 2 qu'en l'absence de borne apparente.
- **Voie 3.** Cotes proportionnelles (`MASS::polr`, BIC), regroupements, rapports des cotes cumulées.
- **Voie 4 (rangs).** Non faite, faute de paquet, avec ses motifs (DECISIONS).
- **Tests.** Les tests du moteur continu dont les données ressemblaient à une échelle forcent désormais `response = "continuous"`.

**Vérifié.** `test-bay-ordinal.R` (2 tests, 4 cas de routage) ; tous les tests bay passent.

**Non vérifié.** La calibration de la voie ordinale n'est pas simulée.

## Séance 11 : 10 novembre 2026, comptages (D17, voie 2)

**Fait.**
- `response = "counts"` : Poisson et binomiale négative comparées, approximation BIC, réponse moyennée sur les deux lois.
- La surdispersion est dite en changement de stratégie.
- Post-hoc par regroupements, avec les rapports des moyennes.
- Même mise en forme que les autres plans.
- Formulation de l'étape 6 corrigée : « ramené vers 0 par l'a priori, surtout quand les groupes sont petits ».

**Vérifié.** `test-bay-counts.R` (5 tests) ; calibration `agent_60` (VALIDATION §13). Sous H0 surdispersé, Poisson seul dépasse 0,95 dans 33 % des jeux, m.bay dans 0 à 0,3 %.

**Audit de de64470 traité.**
- Les comptages sont choisis seuls, et la bascule est dite.
- Le nom de la réponse vient de l'objet passé.
- Les avertissements de `glm.nb` sont capturés.
- Une seule règle pour les facteurs de Bayes, partout.
- « reste > 99,9 % quel que soit le point de départ ».
- Le BIC est vérifié contre BAS : 0,55 en log au plus, sous le seuil de log 3.
- La puissance est documentée (VALIDATION §13).

**Reste.** D17, voies 3 et 4 : l'ordinal et les rangs.

## Séance 10 : 10 novembre 2026, mise en forme commune (D18)

**Cadre.** Consigne d'Antoine relayée par la session parente : des sorties « agréables à lire étape par étape comme m.test() ». J'ai appliqué les 12 recommandations de l'auditeur de lisibilité, avec la limite de 80 colonnes demandée par Antoine.

**Fait.**
- Récit :
  - lignes repliées à 80 colonnes ;
  - résultat dans le titre ;
  - « ==> », « --> » et « Note : » ;
  - jargon dans `detail = TRUE`.
- Même squelette de 7 étapes pour les 7 plans.
- Bilan commun visible d'emblée ; plus de doublon des paires ; une phrase prête pour un compte rendu.
- Les remarques du lecteur D17 sont appliquées :
  - effectifs en clair ;
  - « semblables » évité quand le bilan est indécis ;
  - un seul chiffre pour « aucune différence » ;
  - modalité qui porte l'écart ;
  - facteur de Bayes en faveur de l'absence.
- Bug corrigé : sur un refus de plan, l'étape 1 n'était plus imprimée, car elle restait dans le tampon.

**Mesures** (`relecture/lisibilite_2026-11-09/avant` et `apres`, `mesures.txt`) :
- lignes au-delà de 80 colonnes : de 24 à 38 par plan, à 0 ;
- longueur moyenne : de 104 à 122 caractères, à 55 à 61 ;
- ligne la plus longue : 308 caractères, à 80 ;
- mots : 8 à 28 % de moins ;
- lignes affichées : 37-51 à 60-80, car le repli à 80 colonnes allonge la sortie de lignes plus courtes.

**Vérifié.**
- Les 10 fichiers de tests de m.bay passent.
- Nouveau test de largeur : 7 plans, en français et en anglais, aux largeurs 60, 80 et 120.

**Relecture « après »** par un lecteur débutant indépendant (`relecture/lisibilite_2026-11-09/RELECTURE_APRES.md`).
- **Temps pour trouver le bilan** : 5 à 20 s, contre 20 à 40 s avant.
- **Notes avant → après** :
  - lisibilité : 2 → 4 ;
  - longueur : 2 → 4 ;
  - repérage du bilan : 2 → 5 ;
  - jargon : 2 → 3 ;
  - cohérence entre les plans : 3 → 4.
- **Corrigé à la suite**, sans toucher au calcul :
  - cas indécis : le titre du bilan et la phrase de compte rendu ne donnent plus l'hypothèse en tête comme acquise (« indécis (…) », « Les données ne permettent pas de conclure… ») ;
  - réponse par terme dans la phrase de compte rendu ;
  - explication d'une interaction dont le facteur de Bayes est supérieur à 1 mais la probabilité inférieure à 50 % (a priori de départ) ;
  - échelle du facteur de Bayes (Kass et Raftery 1995) ;
  - un seul vocabulaire pour la stabilité (« l'ordre ne change pas / peut s'inverser ») ;
  - règle fixe pour les paires : différents à partir de 75 %, semblables jusqu'à 25 % ;
  - définitions courtes de « forme du bruit » et « queues lourdes » ;
  - « pas de correction à la Tukey » ;
  - écarts-types par groupe et pourcentages dans les effectifs du discret.
- **Laissé** : moyennes ajustées en ANCOVA, lettres de groupes pour un graphique, moyenne des différences en apparié. Ce sont des sorties nouvelles, à décider plus tard.

## Séance 9 : 9 novembre 2026, réponse binaire ou catégorielle (D17, voie 1), relecture débutant

**Fait.**
- Tableau de contingence à un facteur (`R/sys_bay_contingency.R`) :
  - `contingencyTableBF`, dont l'a priori a été retrouvé et testé ;
  - sensibilité sur a = 1, 2, 4 ;
  - post-hoc par regroupements, en forme close ;
  - tailles d'effet : différence de proportions et log du rapport de cotes, ou écart total des proportions.
- La réponse peut être un facteur, du texte ou TRUE/FALSE.
- Relecture par un lecteur débutant des récits D15 et D16 (`relecture/lecteur_recit/RELECTURE_D15_D16.md`). Trois erreurs de sens sont corrigées :
  - « seul l'a priori le fait pencher » quand l'a priori est à 50 % ;
  - le cas H, où le groupe en tête change selon la covariable : la conclusion le dit ;
  - le cas I, où la confusion entre groupe et covariable est nommée avec les moyennes.
- Les autres remarques (jargon, échelle des qualificatifs, phrase modèle pour le compte rendu) attendent le rapport de l'auditeur de mise en forme, comme demandé par la session parente.

**Vérifié.** `test-bay-contingency.R` (6 tests) ; calibration `agent_59` (VALIDATION §12).

## Séance 8 : 9 novembre 2026, analyse de covariance (D16)

**Cadre.** Plan suivant du PERIMETRE. La session parente a validé le schéma et ajouté trois points : changer de stratégie quand les pentes diffèrent, centrer la covariable et dire son a priori, signaler le recouvrement. Contrôle demandé : `car::Anova` type II et `emmeans`.

**Fait.**
- Covariables centrées, présentes dans chaque modèle ; chaque réponse se lit contre « covariables seules » (`lmBF`, `rscaleCont` = √2/4, qui varie avec r dans la sensibilité). L'interaction groupe × covariable entre dans l'espace héréditaire sans exiger la covariable comme terme.
- Récit :
  - recouvrement des covariables entre groupes, avec « Prudence ajoutée » et soutien plafonné à faible quand il est mauvais ;
  - centre annoncé ;
  - linéarité nommée comme non testée ;
  - moyennes brutes étiquetées comme telles, avec la covariable moyenne par groupe.
- Post-hoc :
  - pentes parallèles : regroupements et écarts ajustés ;
  - pentes différentes probables : changement de stratégie, écarts aux quartiles et zone où les groupes diffèrent (Johnson-Neyman bayésien).
- Bibliographie : Aiken et West 1991, Johnson et Neyman 1936 et Preacher et al. 2006 ont été contrôlées par la session parente ; le DOI de Preacher et al. 2006 a été confirmé sur Crossref. Le doublon Westfall a été retiré.

**Vérifié.**
- `agent_57` concorde avec m.bay à la précision des tirages, pour les BF, les écarts et la zone.
- `agent_58` calibre la méthode : faux positifs rares, mais peu de puissance sur les pentes à n = 20.

**Non vérifié.**
- Relecture des récits G, H et I par un lecteur type.
- Plusieurs facteurs ou plusieurs covariables en simulation.

## Séance 7 : 9 novembre 2026, données appariées et mesures répétées (D15)

**Cadre.** Demande de la session parente : après le post-hoc, traiter les plans du PERIMETRE dans l'ordre suivant : appariés et mesures répétées, puis ANCOVA, puis le discret.

**Fait.**
- Apparié : `ttestBF(paired = TRUE)`, exact. La taille de l'écart est donnée en dz, avec son intervalle de crédibilité.
- Mesures répétées et plan en tranches : `lmBF(whichRandom = id)`, chaque modèle étant lu contre les sujets seuls. Un facteur intra reçoit le post-hoc par regroupements, avec sujet.
- Récit :
  - le plan reconnu, avec intra et inter ;
  - les sujets complets ;
  - les valeurs éloignées, cherchées une fois retiré le niveau de chaque sujet ;
  - la symétrie composée, nommée ;
  - l'étape « De combien les deux niveaux diffèrent » ou « Quelles conditions diffèrent ».
- Les plans avec répliques (modèle mixte) restent refusés.
- Post-hoc : y est standardisé avant le calcul. À l'échelle 1e200, les tirages débordaient et donnaient des intervalles NA (message de bayestestR vu dans la suite).

**Vérifié.**
- Calcul indépendant `agent_55` : l'apparié et 4 sessions concordent avec des intégrales programmées par l'agent. Le post-hoc concorde à 0,3 point près, après correction d'une erreur de dénominateur dans la référence.
- Calibration `agent_56`. Tout est dans VALIDATION §10.

**Non vérifié.**
- Relecture des récits D, E et F par un lecteur type.
- Calibration du plan en tranches et des sujets incomplets.

## Séance 6 : 9 novembre 2026, post-hoc à un facteur (D11)

**Cadre.** La session parente a demandé le post-hoc à un facteur puis à plusieurs : paquets établis, contrôle indépendant, intégration au récit, taille d'effet avec intervalle de crédibilité. Elle a ensuite audité le post-hoc non commité, en quatre points.

**Fait.**
- `R/sys_bay_posthoc.R` (nouveau) :
  - partitions d'égalité, avec un a priori nul plus Ewens (Gopalan et Berry 1998) et les BF de `lmBF` ;
  - écarts et δ par `posterior()` [BayesFactor] et `hdi()` [bayestestR] ; bayestestR passe en Imports.
- Récit : nouvelle étape 6, « Quels groupes diffèrent, et de combien ». La conclusion donne une ligne par paire, avec son qualificatif.
- Audit :
  - « écart observé » (conclusion) et « écart estimé » (post-hoc) sont nommés ;
  - calcul exact jusqu'à 8 groupes ; au-delà, tirage de Gibbs, annoncé comme changement de stratégie avec son erreur de tirage ;
  - au-delà de 15 paires, seules celles qui penchent vers « différents » sont détaillées ;
  - « À lire avec réserve » sous variances inégales, à l'étape 6 et dans la conclusion.

**Vérifié.**
- Calcul indépendant `agent_53` : concordance à l'affichage près.
- Simulations `agent_52` (à côté de Tukey) et `agent_54` (couverture des intervalles).
- Gibbs contre calcul exact à k = 6 et 7.
- Tout est consigné dans VALIDATION §9.

**Non vérifié.**
- Relecture des nouvelles phrases du post-hoc par un lecteur type.
- Calibration sous variances inégales et sous plans déséquilibrés.

## Séance 5 bis : 9 novembre 2026, audit des jalons 2 et 3

**Cadre.** Audit de la session parente : jalons 2 et 3 acceptés, quatre points à régler avant le post-hoc.

**Fait.**
- D2 recompté (`agent_49b`) avec la lecture « seul ou en interaction » : dans le cas A:B:C seule, A, B et C ne sont pas des termes faux. Les fausses interactions de la HOP tombent à 5–6,5 % au-delà de « faible » (P > 0,75).
- Étape 4 :
  - les variances inégales sont annoncées comme un constat, avec le groupe le plus dispersé ;
  - quand les écarts-types sont proches, le constat dit que les groupes diffèrent plutôt par leur forme.
- Règle du qualificatif (D14) : un facteur de Bayes entre 1/3 et 3 donne « non tranché par les données » et un verdict « indécis ». Un terme que les données poussent à l'opposé de l'a priori est « non tranché ».
- Langue : KefiR choisit la langue par `LC_CTYPE`, et `LANGUAGE=en` seul n'y suffit pas. En locale anglaise réelle, le récit sort en anglais (testé).
- Cas rares vus et testés :
  - groupe constant : loi écartée, puis « À vérifier » ;
  - échantillonnage d'importance trop pauvre (draws = 100) : loi écartée.
- Seconde relecture par un lecteur neuf, sur huit cas : 12 points, dont 10 corrigés (`relecture/lecteur_recit/RELECTURE.md`).

**Vérifié.**
- Suites relancées sur le code final, 0 échec : story 8/41, m.bay 11/61, factorial 8/56, oneway 7/27, core 32/237, laws 20/79, m.test 20/52, valreg 6/12, manova 22/59.
- `R CMD check` : Status OK.
- `m.test.R`, `sys_vbse.R` et `sys_msg.R` sont intacts.

**Non vérifié.** La version `recits_v4.txt` n'a pas été relue par un troisième lecteur.

## Séance 5 : 9 novembre 2026, jalon 3 (bilan en récit)

**Cadre.** Demande d'Antoine relayée par la session « KefiR officiel » : un bilan qui raconte son cheminement, comme le verbose de `m.test()`, et finit sur des hypothèses argumentées avec un indice de confiance, sans p-value. Décision D14.

**Fait.**
- `R/sys_bay_story.R` (nouveau) : le récit numéroté est imprimé par `.vbse()` au fil du calcul. « ==> » marque un constat, « --> » une décision, et la fonction et son paquet sont entre crochets.
- Les avertissements sont dits à l'étape où ils naissent. Le récit est gardé dans `$story` (en, fr) et `print(x, story = TRUE)` le rejoue.
- Étapes à un facteur : plan, contrôles (effectifs, ex-aequo, groupes constants, valeurs éloignées selon Tukey), point de départ, formes du bruit, facteur de Bayes, précision et sensibilité, conclusion.
- Étapes à plusieurs facteurs : plan et espace héréditaire, contrôles par case, forme du bruit (normale seule, dit), a priori HOP avec son motif, facteurs de Bayes [lmBF() {BayesFactor}], précision et sensibilité par terme, conclusion.
- Conclusion : le bilan court, en retrait, sans titre ni point de départ ni avertissements déjà dits. Chaque hypothèse et chaque terme reçoit une ligne « Repose sur » ou « Contre elle ».
- Un terme est dit « Non tranché » quand le soutien est nul ou que la sensibilité le renverse. Quand le facteur de Bayes est sous 3, la ligne dit que la réponse repose surtout sur l'a priori.
- Sous queues lourdes, l'écart est donné sur les médianes.
- Moyennes par case quand une interaction est testée.
- Relecture par un lecteur type (étudiante de BUT et ingénieur qualité) sur six cas : 14 points, dont une erreur, l'écart des « centres » qui était en fait celui des moyennes. 11 sont corrigés, 3 laissés avec leur motif (`relecture/lecteur_recit/RELECTURE.md`).

**Vérifié (comment).**
- `test-bay-story.R` (nouveau, 5 tests, 29 assertions) :
  - `verbose = FALSE` n'imprime rien ;
  - six étapes numérotées dans l'ordre ;
  - pas de titre ni de point de départ répétés dans la conclusion ;
  - une ligne d'argument par hypothèse et par terme ;
  - aucune p-value chiffrée ;
  - `print(story = TRUE)` est identique à l'impression directe ;
  - un plan refusé est nommé à l'étape 1 avant le refus.
- Suites relancées sur le code final, 0 échec : story 5/29, m.bay 11/61, factorial 8/56, oneway 7/27, core 32/237, laws 20/79, m.test 20/52, valreg 6/12.
- `R CMD check --no-manual --ignore-vignettes` : Status OK.
- `R/m.test.R`, `R/sys_vbse.R` et `R/sys_msg.R` sont identiques à 3675fc7.

**Non vérifié.**
- `test-manova.R` n'a pas été relancé (26 min ; aucun fichier qu'il touche n'a changé depuis le jalon 2, où il passait).
- La version corrigée du récit n'a pas été relue par un second lecteur.
- Le récit n'a pas été vu en session anglaise.
- Le récit des cas rares (loi écartée pour ESS, cas divergent) est écrit mais n'a été vu sur aucun jeu de données.

**Défaut signalé, non corrigé (hors périmètre).** `.vbse()` remplace tout texte anglais d'une seule chaîne par le texte français (`if (length(ang)==1) ang <- fr`). En session anglaise, le verbose de `m.test()` sort donc en français. `m.bay()` contourne ce défaut en choisissant la langue avant l'appel.

**Reproduction.** Six récits : `Rscript dev/mbay/relecture/lecteur_recit/lecteur_story.R (écrit dans le scratchpad, chemin à adapter)`, dont le texte est recopié dans `dev/mbay/relecture/lecteur_recit/recits_v2.txt`. Tests : `testthat::test_file("tests/testthat/test-bay-story.R")`.

## Séance 4 : 7 au 9 novembre 2026, jalon 2 (plusieurs facteurs), sous-session auditée par « KefiR officiel »

**Cadre.** Audit du jalon 1 accepté, avec six points (VALIDATION §7). Le moteur multi-facteurs interne avait été arrêté par le filtre de sécurité : non repris. Décision d'Antoine le 9 novembre (« Oui, via BayesFactor ») : D1 tranchée, puis son principe pour la suite (appeler des paquets établis). Simulations interrompues dans la nuit du 7 au 8 (machine), reprises depuis leurs tranches.

**Fait.**
- `R/sys_bay_factorial.R` : plusieurs facteurs croisés, loi normale à variance commune, vraisemblances `BayesFactor::lmBF` ; espace héréditaire, loi HOP, probabilité par terme, sensibilité (r × trois lois sur les modèles), bilan court. BayesFactor en Imports, `importFrom(lmBF, extractBF)`. Un facteur : calcul interne inchangé.
- D2 tranchée : HOP(1,1) (`agent_47`, `agent_49`). D1 consignée avec le principe d'Antoine.
- Validation : référence Monte-Carlo indépendante des lois hétérogènes (`agent_45`, `agent_48`) ; seuil des queues lourdes (`agent_46`) ; référence indépendante à deux g pour `lmBF` (`agent_51`) ; calibration et plans déséquilibrés (`agent_50`) ; relecture adverse (`RELECTURE_ADVERSE_2026-11-09.md`, dix points traités).
- Deux défauts de BayesFactor contournés : `integrate` en échec à très grand effet (NA), erreur annoncée non fiable sur les modèles à plusieurs termes.

**Vérifié (comment).** Voir VALIDATION §7 et §8 ; tests relancés sur le code final (7ad4dc2) : 8 fichiers, 126 tests, 582 assertions, 0 échec (dont test-m.test 20/52, test-manova 22/59, test-valreg 6/12) ; `R CMD check --no-manual --ignore-vignettes` : Status OK (une note de xref sur `ez`, absent de mbay-lib). `m.test()` intact : `git diff 3675fc7 HEAD -- R/m.test.R` vide ; aucun fichier de `R/` hors `m.bay.R` et `sys_bay_*` modifié.

**Non vérifié.** Facteurs à plus de 3 niveaux ; 4 facteurs ; plans très déséquilibrés ; r ≠ 0,5 dans la référence indépendante ; cas `agent_49` sur A:B:C (en cours) ; lisibilité du bilan multi-facteurs par un lecteur type ; queues lourdes et variances inégales à plusieurs facteurs (hors périmètre de ce jalon, annoncé dans le bilan).

**Écarts de méthode.** Deux fois l'argument `f=` capté par `fun` de `parLapplyLB` (agent_47, agent_50) : renommé. Les outils d'écriture transforment `\` et `\uXXXX` : patches par numéro de ligne et `chr(92)`. roxygen 8 réécrit NAMESPACE et d'autres Rd : rétablis à la main, seul `man/m.bay.Rd` régénéré.

**Reproduire**, depuis le worktree :
```bash
Rscript dev/mbay/scripts/agent_49_scenario_additif.R 300 200 5
Rscript dev/mbay/scripts/agent_50_multi_validation.R A 1200 3
Rscript dev/mbay/scripts/agent_50_multi_validation.R B 300 2
Rscript dev/mbay/scripts/agent_51_reference_2g.R
Rscript dev/mbay/scripts/agent_48_compare_mc.R
```

## Séance 3 : 6 et 7 novembre 2026, sous-session auditée par « KefiR officiel »

**Cadre.** Consigne de la session « KefiR officiel », en remplacement de la routine d'une heure. Démarrée à la fin de l'ancienne routine (`f2b8d9d`, worktree propre, vérifié). Interruption dans la nuit du 6 au 7 (PC endormi, cause non confirmée) : simulations arrêtées et reprises ; elles sont désormais reprenables, 3 ou 4 cœurs, chargement léger.

**Objectif.** Jalon 1 : `m.bay()` à un facteur, 2 et k groupes, de bout en bout.

**Fait.**
- `R/m.bay.R` (exporté), `R/sys_bay_input.R`, `R/sys_bay_laws.R`, `R/sys_bay_bilan.R` ; `.bay_bf_oneway_hetero` retirée ; `man/m.bay.Rd` ; `tests/testthat/test-m.bay.R`, `test-bay-laws.R`.
- Huit modèles comparés et moyennés, D6 faisable à k > 2. Deux formes d'effet à k ≥ 3, à la suite de la relecture.
- Plans autres qu'un facteur : refus explicite, classement déjà fait.
- Deux relectures adverses (`RELECTURE_ADVERSE_2026-11-07.md`). Scripts `agent_40` à `agent_44`, `agent_load_bay.R`. METHODOLOGIE §7, DECISIONS D1, D4, D5, D6, D9, VALIDATION §7.

**Vérifié (comment).**
- Tests relancés par moi sur le code final (`f8fe497`, code inchangé depuis) : `test-bay-core` 237, `test-bay-oneway` 27, `test-m.bay` 57, `test-bay-laws` 77, `test-m.test` 52, `test-valreg` 12, `test-manova` 59 : **0 échec, 0 erreur**. `m.test()` n'est pas touché (`git diff 3675fc7 HEAD -- R/m.test.R` vide).
- Calcul indépendant : `agent_40` (8 modèles à k=2 contre quadrature, ≤ 0,02), `agent_42` (seconde forme, t, k=3, ≤ 0,022).
- Calibration : `agent_41 A`, 6 plans, N=600, moyenne annoncée à 0,014 au plus de la fréquence, dans l'erreur de simulation. Cellule à z=3,64 du premier lancement : non reproduite, `agent_41 A2`.
- Effets fixes : `agent_41 B`, 66 scénarios, N=300 ; à variances inégales, (5,5,40) sd (4,4,1) sous H0 : 1,3 % de P>0,95 contre 19 % en gaussien seul.
- B1, B2 : `agent_41 D`. Sensibilité : `agent_41 C`, N=100. Formes différentes : `agent_43`. Variances étagées : `agent_44`.

**Lancements, à ne pas confondre.** `agent_41 B` scénarios 1 à 8, C, D : code final, copie `f8fe497`. `agent_41 B` scénarios 9 à 66 et A : `c349c30`, avant le pilote itératif et les règles d'ex-aequo finales ; ces changements ne touchent pas ces données, et les scénarios 1 à 8 ont redonné les mêmes chiffres à ±0,01. `agent_41 A2` : `a315c6d`, k=2, modèle identique. `agent_40` : avant le centrage de `.bay_law_kernel_grid`, sans effet hors c extrêmes (écart ≤ 2e-13 selon le relecteur). `agent_42`, `agent_43` : `c349c30` puis copie de travail. **Non relancés sur le code final.**

**Non vérifié.**
- Pas de calcul indépendant pour t à variance par groupe à k ≥ 3, ni pour la seconde forme à variance par groupe : calibration seule.
- Relecture par un tiers : pas sur le verdict 1:2–2:1, l'avis sur les queues lourdes, le pilote itératif, les règles d'ex-aequo finales (dérivées à la main).
- k=20 : lois hétérogènes souvent écartées, une graine seulement.
- Lisibilité du bilan : relue par moi seul, pas par un lecteur étudiant ou industriel.
- `devtools::check()` non lancé. `m.bay()` sans post-hoc, taille d'effet ni intervalle de crédibilité : jalon suivant.
- Bibliographie : Maier 2025 par résumeur automatique, définition de α selon δ non lue.

**Résultat principal.** Comparer les lois marche à k > 2, sans Stan : calibré, et les faux « Oui » du plan hétéroscédastique tombent de 19 % à 1,3 %. La relecture a trouvé un défaut plus grave que prévu : un petit groupe éloigné parmi 20 mettait H0 en tête. D'où la seconde forme d'effet, qui fait perdre l'identité avec `anovaBF` à k ≥ 3. **Défaut qui reste, annoncé dans le bilan** : à moyennes égales et formes différentes, la loi t compare des centres (P>0,95 dans 67 % à n=400).

**Décisions prises seul, réversibles.** Dirichlet(1,5) ; échelle par groupe ; deux formes à parts égales ; r fixe, l'hyper-a priori ne changeant rien à ±0,03 ; cotes 1:2–2:1 dans le verdict (jugement) ; seuils du qualificatif 0,75/0,95/0,993 ; ESS < 30 écarte une loi ; D1 : calcul interne à un facteur, BayesFactor seulement en test ; D9 option B.

**Écarts de méthode.** L'arrêt forcé de mes propres processus R de simulation a été refusé par le mode automatique ; je n'ai pas contourné, ils ont fini seuls. Machine saturée en mémoire au premier lancement (chaque processus chargeait tout KefiR) : d'où `agent_load_bay.R`. Un prototype a eu un bug de graine (`set.seed` dans une boucle, tous les jeux identiques) : chiffres jetés, jamais inscrits. Écrire avec l'outil Write transforme `\uXXXX` en accents : sources réécrites en ASCII, testé.

**Défauts de `m.test()` trouvés**, non corrigés, déjà signalés séance 1 : voir plus bas ; rien de nouveau cette séance.

**Reproduire**, depuis le worktree. Voir VALIDATION §7 pour les sorties attendues.
```bash
Rscript -e '.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib",.libPaths())); devtools::test(filter="bay")'
Rscript dev/mbay/scripts/agent_40_lois_calcul_independant.R
Rscript dev/mbay/scripts/agent_42_forme_un_groupe_t.R
Rscript dev/mbay/scripts/agent_41_lois_simulation.R A 600 3
Rscript dev/mbay/scripts/agent_41_lois_simulation.R B 300 3
Rscript dev/mbay/scripts/agent_41_lois_simulation.R C 100 3
Rscript dev/mbay/scripts/agent_41_lois_simulation.R D 1 1
Rscript dev/mbay/scripts/agent_43_formes_differentes.R 200 2
Rscript dev/mbay/scripts/agent_44_variances_extremes.R 100 3
```
Les tranches déjà calculées sont dans `dev/mbay/scripts/resultats/agent_41_cas/` (non commité) : les effacer pour tout recalculer. Comptez plusieurs heures.

**Prochain jalon, sans attendre.** Post-hoc et taille d'effet à un facteur, ou plusieurs facteurs, selon l'audit. Essai de moteur interne à plusieurs g concluant dans le scratchpad, à reprendre.

## Séance 2 : vendredi 6 novembre 2026

**Durée.** Lancée à 09:12, bloquée jusqu'à 13:40 (attente d'autorisations), budget d'une heure repris à 13:40 sur message de ta session « KefiR officiel », en mode Auto. Fin vers 14:45.

**Objectif.** Sans attendre tes décisions : voie à variances par groupe sans Stan, recoupement de la multiplicité 2³, échec d'`integrate()`, réécriture des §3 et §4.

**Fait.**
- `R/sys_bay_oneway.R` (`.bay_bf_oneway_zs`, `.bay_bf_oneway_hetero`, non exportés) et `tests/testthat/test-bay-oneway.R`.
- METHODOLOGIE §3 et §4 réécrits, marqués « en attente de D2 » et D3 ; §7 complété ; 0,707 n'est plus attribué à Rouder 2009.
- BIBLIOGRAPHIE : 10 références confrontées à leur source, section sur les variances inégales.
- Scripts `agent_20` à `agent_30`, `relecture/agent_rev2_*`, `methodo_20261106/`.

**Vérifié (comment).**
- Tests relancés par moi : `test-bay-oneway` 80 assertions, `test-bay-core` 237, `test-m.test` 52, 0 échec. Après ce lancement, seul un commentaire a été ajouté à `sys_bay_oneway.R` ; fichier réanalysé par `parse()`, tests non relancés. `test-manova`, `test-valreg` non lancés (aucun fichier existant modifié).
- Calcul de hetero : grille écrite à part, 3 jeux à k=2 (`agent_29`) : écart de l'ordre de l'erreur annoncée.
- 20 % de P>0,95 en gaussien hétéroscédastique : recalculé, N=1500, autre graine (`agent_23`). L'écart 16 % / 25,5 % est du hasard.
- Hetero sur ce plan : 0/300 à P>0,95 (`agent_29`), 0,7 % sur N=1000 chez l'agent.
- Multiplicité sous la loi de D2 : `agent_22`, lancé par moi, N=40.
- Défauts B1, B2 et I6 reproduits (`agent_30`).

**Non vérifié.**
- **La voie hetero n'est pas validée, et sa loi a priori est rejetée en l'état.**
- Puissance à N=16-17 par case, calibration à N≈250, sensibilité à N=72, troncature à N=7-8 : agent, faible résolution.
- Multiplicité : `agent_22` à N=40 seulement ; forme finale des deux lignes non simulée.
- Plusieurs défauts du relecteur (I1 à I5, I7, I8), non refaits. Corrections non programmées. Le contrôle de lisibilité du bilan n'a pas eu lieu : pas de bilan.
- Bibliographie : résumés lus par résumeur automatique ; Bertolino 2000 et Cano 2017 à lire ; attribution de 0,707 à confirmer dans le texte de Rouder 2009 ; 61 références sur 71 toujours au niveau métadonnées.

**Résultat principal.** Le calcul à un facteur sans Stan est possible et juste, mais le modèle proposé échoue : un petit groupe éloigné de 50σ donne log BF -6,9 alors que Welch donne p=3e-15, et plus il s'éloigne, plus le BF soutient H0. Il corrige le faux « Oui » du plan (5,5,40) mais crée un faux « Non » : le problème est déplacé, pas résolu. D6 reste ouvert.

**Écarts de méthode.** Ma consigne de simulation à l'agent était ambiguë (« π^ordre ») : d'où la mauvaise loi dans `agent_20`, repérée par l'agent de méthodologie et refaite (`agent_22`). J'ai aussi fait une recherche Gmail sur l'objet « KefiR m.bay » en début de séance, avant la consigne de ne pas le faire : seul le point d'hier en est sorti, aucune réponse de ta part. Pas de nouvelle recherche depuis.

**Décisions prises.** Aucune décision structurante. Lignes « hérédité » d'`agent_20` à ne pas citer.

**Arbitrage reçu le soir même.** D2, loi a priori des hypothèses : facteur 1/2, interaction 1/8, tranché par Antoine (« Oui, retiens-les ») sans connaître la simulation. Elle ne contredit pas ce choix : 25 % contre 49 % de jeux avec un faux terme à P>0,5 sous H0. Coût mesuré (`agent_31`, partiel) : interaction réelle plus souvent manquée à petit effectif, P>0,5 dans 37,5 % contre 51 % à n=10 par case, f=0,5. Voir `DECISIONS.md`. Le libellé de la question par terme, D1 et D9 restent ouverts.

**Prochaine séance, sans attendre.** Étudier l'orientation d'Antoine sur D6 (`DECISIONS.md`) : comparer et moyenner les lois d'échantillonnage, avec une loi propre sur les rapports d'échelle ; lire Maier 2025 et Dablander 2024 ; dire franchement si k>2 n'est pas calibré ou trop lent. Modèle à un g par groupe ou τ par contraste, testé d'abord contre B1 et B2. Corriger I6 et I7 de zs. Refaire puissance et calibration à N suffisant. Lire Bertolino 2000 et Cano 2017.

## QUESTIONS POUR ANTOINE (à jour au 6 novembre)

Les cinq questions du 5 novembre, plus bas, restent ouvertes. Éléments nouveaux :
- **D2.** Plan 2³ sous H0, N=40, au moins un terme à P>0,5 : hérédité 25 %, uniforme 49 %, apparié 60 %. L'écart vient en partie de l'a priori lui-même. À 3 facteurs, la loi retenue dans D2 place A:B:C à 1/128.
- **D6.** La voie gaussienne reste indéfendable en hétéroscédasticité (20 %) et la voie de remplacement écrite aujourd'hui est rejetée. Rien à arbitrer tant que je n'ai pas mieux.

## Séance 1 : jeudi 5 novembre 2026

**Durée.** Lancée vers 16:17 ; environ 70 minutes perdues en attentes d'autorisation. Sur ton accord (« tu peux continuer aujourd'hui »), la séance a dépassé l'heure prévue. Fin vers 18:35.

**Objectif.** Première séance : comparer les moteurs, définir le périmètre, écrire feuille de route et décisions, poser les questions. Puis avancer sur ce qui ne dépend pas de ces décisions : le socle interne.

**Fait.**
- Worktree `KefiR-mbay` et branche `feat/m.bay` créés depuis master ; paquets dans `mbay-lib` uniquement.
- `PERIMETRE.md` (environ 60 cas de `m.test()` exécutés), `METHODOLOGIE.md`, `BIBLIOGRAPHIE.md` (92 références), `DECISIONS.md`, `FEUILLE_DE_ROUTE.md`, `VALIDATION.md`, relecture adverse.
- Socle interne `R/sys_bay_core.R` et ses tests. Aucun fichier existant modifié, pas d'export, pas de `m.bay()`.
- Consigne de la tâche quotidienne modifiée à ta demande : mail quotidien, et retrait d'un titre erroné.

**Vérifié (comment).**
- `test-bay-core.R` : 237 assertions, 0 échec, relancé par le coordinateur.
- `test-m.test.R` : 20 tests, 52 assertions, 0 échec, exécuté jusqu'au bout. `test-manova.R` et `test-valreg.R` non lancés (aucun fichier existant modifié).
- JZS contre `ttestBF` recalculé : écart 1e-11.
- Calibration sous l'a priori à un facteur recalculée (N=2000) : voir VALIDATION.
- Croisement 2×2 et hétéroscédasticité recalculés : défauts confirmés.
- Anomalies `m.test` ci-dessous reproduites sur des données à moi.
- 3 DOI recontrôlés sur Crossref (Rouder 2009, Rouder 2012, Liang 2008) ; 8 références tirées au hasard confrontées par le relecteur : 1 erreur corrigée (Gelman 2008), 1 réglage reformulé (Chipman 1996).

**Non vérifié.**
- `m.bay()` n'existe pas : **aucune voie n'est validée**.
- 71 références sur 92 : métadonnées seules, la ligne « apport » vient de la mémoire de l'agent, à relire avant citation. 2 intitulés candidats non confirmés.
- Multiplicité 2³, ROPE f et MANOVA en √d : chiffres du relecteur seul.
- Écart 16 % contre 25,5 % de faux « Oui » en hétéroscédasticité : non élucidé.
- Voies Stan, brms avec Rtools, BAS, sphéricité.
- Sensibilité à r, k=2, n=10 et 100 : non exécutés.
- MANOVA non normale dans `m.test` : arrêtée après plus de 6 minutes, cause non établie.

**Résultat principal.** La méthodologie initiale ne tient pas. Le « Oui p / Non 1-p » par terme avec BF apparié est faux dès qu'une interaction joue : sur un croisement pur, P(A)=0,215 alors que P(A:B)=1,000 et que A agit manifestement. La simulation confirme : pas de générateur unique raisonnable sous lequel ces probabilités sont calibrées ensemble. Le calcul BayesFactor lui-même est sain (écarts minimes, calibration à un facteur). C'est la **définition de la question** qui est à revoir : voir D2.

**Décisions prises.** Aucune décision structurante. Tout reste « PROPOSÉ » dans `DECISIONS.md`. BayesFactor n'est pas ajouté à DESCRIPTION.

**Anomalies de `m.test()` hors m.bay, en production (non corrigées, on ne touche pas master).**
1. **Split-plot `y~G*T+Error(id/T)` : interaction retirée sans message.** Sur mes données, `aov` donne G:T p=0,039 ; `$effects` n'a que G et T.
2. Le même appel **plante avec `verbose=FALSE`** : « $ operator is invalid for atomic vectors ».
3. **`sys_variance.R:240` : `bf.test` compare des moyennes, pas des variances**, alors que le code l'utilise pour confirmer Levene. SD 0,5/2/5 : Levene p=4,6e-15, `bf.test` p=0,67, hétéroscédasticité écartée. Sur le plan n=(5,5,40), sd=(4,4,1), sous H0 : `m.test` choisit surtout l'ANOVA paramétrique, 26,7 % de p<0,05 (N=30 seulement, indicatif ; `agent_13`).
4. `(1|id)` : post-hocs non appariés (t1-t2 0,139 contre 0,002 avec Error ; chiffre de l'agent, lu et non refait).
5. `m.test(formula=y~F, data=d)` échoue : « doit être une formule valide ».
6. `.vbse` et `.exit` : `length(ang)==1` remplace `ang` par `fr`. D'après le code, un message anglais d'une seule chaîne ne s'affiche jamais (lu, non essayé en locale anglaise).

## QUESTIONS POUR ANTOINE (état du 5 novembre)

1. **Moteur (D1).** BayesFactor en Imports (+5 paquets, 17 Mo, sans compilateur), rstanarm seulement en Suggests pour les cas lourds, brms écarté, MANOVA en code interne : d'accord ?
2. **Question posée par terme (D2), la plus importante.** « Ce facteur joue-t-il un rôle, seul ou en interaction ? » en première ligne, avec l'effet principal en second : d'accord ? Et quelle loi a priori : uniforme sur les modèles (facteur 0,6, interaction 0,2) ou hérédité forte (facteur 1/2, interaction 1/8) ?
3. **`m.test` et `m.bay` partagent-ils la normalisation des entrées (D9) ?** Cela demande de toucher `m.test()`, sous non-régression. L'alternative est de recopier.
4. **Anomalies de `m.test()` ci-dessus.** Veux-tu que je les corrige, dans une tâche séparée ? Les points 1 et 3 touchent des résultats que les étudiants peuvent obtenir aujourd'hui.
5. **Aiguillage (D6).** Accepter que `m.bay` ne suive pas l'arbre de `m.test`, et ne produise pas de « Kruskal bayésien » pour k>2 ?

**Prochaine séance, sans attendre ces réponses.** Écrire et simuler la voie interne à variances par groupe (D6). Recouper la multiplicité 2³ sur `anovaBF`. Résoudre l'échec de `integrate()` à F très grand. Réécrire §3 et §4 de METHODOLOGIE en attente de D2.
