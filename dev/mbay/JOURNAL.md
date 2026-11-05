# Journal de m.bay

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

## QUESTIONS POUR ANTOINE

1. **Moteur (D1).** BayesFactor en Imports (+5 paquets, 17 Mo, sans compilateur), rstanarm seulement en Suggests pour les cas lourds, brms écarté, MANOVA en code interne : d'accord ?
2. **Question posée par terme (D2), la plus importante.** « Ce facteur joue-t-il un rôle, seul ou en interaction ? » en première ligne, avec l'effet principal en second : d'accord ? Et quelle loi a priori : uniforme sur les modèles (facteur 0,6, interaction 0,2) ou hérédité forte (facteur 1/2, interaction 1/8) ?
3. **`m.test` et `m.bay` partagent-ils la normalisation des entrées (D9) ?** Cela demande de toucher `m.test()`, sous non-régression. L'alternative est de recopier.
4. **Anomalies de `m.test()` ci-dessus.** Veux-tu que je les corrige, dans une tâche séparée ? Les points 1 et 3 touchent des résultats que les étudiants peuvent obtenir aujourd'hui.
5. **Aiguillage (D6).** Accepter que `m.bay` ne suive pas l'arbre de `m.test`, et ne produise pas de « Kruskal bayésien » pour k>2 ?

**Prochaine séance, sans attendre ces réponses.** Écrire et simuler la voie interne à variances par groupe (D6). Recouper la multiplicité 2³ sur `anovaBF`. Résoudre l'échec de `integrate()` à F très grand. Réécrire §3 et §4 de METHODOLOGIE en attente de D2.
