# Relecture du récit de m.bay par un lecteur type (9 novembre 2026)

Le relecteur est un agent qui n'a lu que `recits.txt` : six sorties de `m.bay()` en `verbose = TRUE`, sans code ni documentation. Il a joué deux lecteurs : une étudiante de BUT génie biologique (2e année) et un ingénieur qualité habitué à Minitab. Pour chaque cas, il devait dire :
- ce qu'il a compris de la conclusion ;
- si la fonction a changé de stratégie, et pourquoi ;
- ce que vaut la confiance finale ;
- les phrases qui le bloquent.

`recits_v2.txt` montre les mêmes six cas après corrections. Les cas 2 et 3 se comprenaient déjà sans difficulté ; la structure en six étapes et l'annonce des a priori avant les données ont été jugées claires.

| # | Défaut relevé | Gravité | Suite |
|---|---|---|---|
| 1 | Cas 4 : « écart entre centres (M2 - M1 = 0,0133) ». Le chiffre est l'écart des moyennes, alors que les médianes diffèrent de 0,31. | **erreur** | **corrigé** : sous queues lourdes, l'écart est donné sur les médianes, avec l'écart des moyennes à côté |
| 2 | Cas 4 : « la loi normale donne une différence à 7,4 % » a été lu comme une p-value. 99,4 % à côté de « indécis » paraît contradictoire. | grave | **corrigé** : la phrase dit « selon la loi normale (qui compare les moyennes), la probabilité d'une différence est de 7,4 % ; … (qui comparent les centres) 99,4 % », puis « on ne retient pas 99,4 % tel quel ». « Une forme domine » remplace « aucune forme n'est choisie » au-delà de 95 %. |
| 3 | « oui/non avec aucun soutien » : le oui ou le non est lu seul. | grave | **corrigé** : « Non tranché (probabilité du oui …) » dès que le soutien est nul ou que la sensibilité renverse le terme |
| 4 | Cas 5 : « non avec soutien modéré » alors que les données ne divisent les cotes que par 1,21 | important | **corrigé** : pour un facteur de Bayes sous 3, la phrase devient « les données déplacent peu ses cotes : le non repose surtout sur l'a priori » |
| 5 | « Changement de stratégie » employé pour ce qui n'est qu'une prudence de formulation | important | **corrigé** : « Changement de stratégie » est réservé aux changements de calcul (loi écartée) ; « Prudence ajoutée » sert pour le plafonnement, la sensibilité et le non tranché |
| 6 | Jargon : « point de départ », « formes du bruit », « cotes », « Dirichlet », « centres », « seuil 95 % » (rapproché de 5 %) | important | **corrigé** : titre « Point de départ (a priori) » ; formes du bruit et centre définis à leur première mention ; « cotes (rapport oui/non) » ; Dirichlet retiré du récit (gardé en `detail`) ; « la méthode demande 95 % de probabilité » |
| 7 | Cas 1 : « Repose sur le plus grand écart » est lu comme un test sur cette paire ; « favorisent l'absence » est trop affirmatif pour un facteur de Bayes de 0,9 | important | **corrigé** : « l'ensemble des écarts entre groupes, le plus grand… » ; facteur de Bayes sous 3 : « à peine plus probables », « les données départagent à peine » |
| 8 | « Erreur numérique 1,7 % » est lu comme un risque d'erreur | moyen | **corrigé** : « imprécision du calcul par simulation… (ce n'est pas un risque d'erreur) » |
| 9 | Les chances a priori des termes (50 + 50 + 16,7 + 33,3) ne font pas 100 | moyen | **corrigé** : « pris seul (ces chances ne s'additionnent pas) » |
| 10 | Cas 6 : la réponse globale repose sur un terme non tranché ; les moyennes de case manquent alors qu'une interaction est testée | moyen | **corrigé** : « repose sur l'ensemble des termes ; aucun n'est tranché seul » ; ligne « Par case » quand une interaction est testée (jusqu'à 16 cases) |
| 11 | « note diffère-t-il » (accord) ; « dépend du protocole » lu comme le facteur « protocole » ; note de développeur sur lmBF | mineur | **corrigé** : « Différence de note selon groupe ? », « dépend de la façon dont l'étude a été menée », note retirée |
| 12 | « 84,4 % » au lieu de 72,0 + 12,5 = 84,5 % (arrondis) | mineur | laissé : le total est exact, c'est la somme des arrondis qui diffère |
| 13 | L'hypothèse classée première à 52,6 % est lue comme gagnante | mineur | laissé : le verdict « indécis » et « proche de pile ou face » le disent sur les lignes voisines |
| 14 | Données ordinales (1 à 4) : l'avertissement « discrète » n'a pas de suite dans la conclusion | hors périmètre | laissé : la voie discrète n'est pas encore validée (PERIMETRE) |

Non revérifié par un second lecteur : la version corrigée (`recits_v2.txt`).
