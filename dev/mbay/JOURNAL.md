# Journal de m.bay

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
