# Journal de m.bay

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
