# Décisions structurantes de m.bay

État : **ARBITRÉ** = tranché par Antoine ; **DÉCIDÉ** = tranché ici par la session de développement, sur consigne d'autonomie d'Antoine du 6 novembre, réversible si des résultats le contredisent ; **PROPOSÉ** = en attente d'un jalon à venir, réversible sur cette branche ; **ACQUIS** = imposé par sa consigne ou sans alternative sérieuse. Chiffres : `VALIDATION.md`. Modèles et formules : `METHODOLOGIE.md`.

## D1. Moteur de calcul : ARBITRÉE par Antoine le 9 novembre 2026 pour plusieurs facteurs

**Plusieurs facteurs : BayesFactor en Imports.** Mots d'Antoine, transmis par « KefiR officiel » : « Oui, via BayesFactor ». Les vraisemblances marginales des modèles viennent de `BayesFactor::lmBF` ([BayesFactor], [Rouder2012]) ; KefiR n'écrit que l'assemblage (`R/sys_bay_factorial.R`) : espace des modèles héréditaires, loi a priori D2, probabilité par terme, sensibilité, bilan. Loi normale à variance commune seulement à ce jalon ; queues lourdes et variances inégales à plusieurs facteurs viendront plus tard, sans Stan si possible. **Un facteur** : le calcul interne déjà validé (huit modèles) est gardé ; à un facteur, `lmBF` et la première forme de `normal_eq` concordent à 1e-4 (k = 2, 3, 5 ; test `test-bay-factorial.R`).

**Principe pour la suite, ARBITRÉ par Antoine le 9 novembre** (ses mots, transmis par « KefiR officiel ») : « Oui, comme m.test() qui croise de nombreuses libraires R, m.bay() doit pouvoir faire de même. » Pour chaque étape suivante (ANCOVA, apparié, mesures répétées, mixte, MANOVA, discret, voies robustes), chercher d'abord un paquet R établi qui fait le calcul (BayesFactor, bayestestR, BAS ; rstanarm ou brms pour les voies lourdes) et l'appeler, plutôt que d'écrire un moteur. Le code interne reste pour ce qu'aucun paquet ne fait, comme la comparaison de lois à un facteur. Pour chaque paquet : Imports ou Suggests selon son poids et sa nécessité, motif consigné ici, version et référence citées, contrôle indépendant. Dépendance lourde (compilateur, Stan) : Suggests, avec un message clair si elle manque.

*Proposition d'origine, gardée pour mémoire :*

| Plan | Moteur proposé | Coût pour `install_github` |
|---|---|---|
| 2 groupes, k groupes, multi-facteurs, ANCOVA, apparié, mesures répétées, split-plot (gaussien) | **BayesFactor en Imports** | 5 paquets de plus que KefiR, environ 17 Mo, 35 s, sans compilateur |
| BF d'inclusion par terme | Code interne (`R/sys_bay_core.R`) | aucun |
| MANOVA | Code interne (formule Normal-Inverse-Wishart) | aucun |
| Pentes aléatoires, non gaussien, robuste | **rstanarm en Suggests**, appel conditionnel avec message d'installation | au moins 26 paquets, 207 Mo, 189 s, seulement pour qui le demande |

**Écartés.**
- **brms par défaut** : échoue sans compilateur (« 'make' not found » mesuré) ; un étudiant sous Windows n'a pas Rtools en général. Avec Rtools 4.4 hors PATH, la compilation a dépassé 6 minutes sans conclusion.
- **rstanarm en Imports** : 207 Mo au minimum imposés à tous les utilisateurs de KefiR, y compris ceux de `m.test()`.
- **BayesFactor en Suggests** : `install_github` ne l'installerait pas, `m.bay()` serait cassé d'origine pour son usage principal.
- **BAS** : installé, jamais comparé. Non écarté sur le fond, reste à évaluer pour l'exploration au-delà de 4 facteurs.
- **bayestestR en Imports** : son calcul d'inclusion est reproduit à la main ; pas de raison d'ajouter deux paquets.

**État au 7 novembre, PROPOSÉ.** La voie à un facteur n'emploie pas BayesFactor à l'exécution : calcul interne, sans compilateur et sans dépendance (`R/sys_bay_oneway.R`, `R/sys_bay_laws.R`). BayesFactor sert de contrôle dans les tests et reste hors de DESCRIPTION à ce stade. **Le jalon suivant dira si la proposition ci-dessus tient** : les variances par cellule et les lois t ne passent pas par `lmBF`, qui n'accepte pas de poids. Un essai de moteur interne à plusieurs g concorde avec `lmBF` sur 4 modèles à 2 facteurs et 3 facteurs (écarts de 0,01 à 0,09 sur log BF, de l'ordre des erreurs annoncées, scratchpad, un seul jeu) : à valider, pas validé.

**Réserves à connaître.**
- BayesFactor ne protège pas : `family=binomial()` et `var.equal=FALSE` sont acceptés et ignorés en silence, `cbind` et `(1+x|id)` donnent une erreur obscure, `x:id` aléatoire donne `BF=NA`. **`m.bay()` doit refuser lui-même ces cas.**
- L'erreur numérique annoncée monte jusqu'à 80 % sur certains modèles à 3 facteurs (n=100) : elle doit être lue, recalculée (`recompute`) ou signalée à l'utilisateur.
- 4 facteurs : 205 s à n=100, non terminé à n=500 au bout de 5 minutes. Temps à arbitrer.
- Sur les mêmes données, **rstanarm + bridgesampling donne BF=0,207 et anovaBF BF=3,46**, parce que les a priori diffèrent. Changer de moteur peut renverser la conclusion : on ne bascule jamais de moteur sans le dire, et le bilan annonce l'a priori de la voie employée. Le défaut de rstanarm n'est pas comparable au défaut JZS : une voie Stan devra fixer des a priori alignés sur le défaut, et cette équivalence reste à établir.

## D2. Définition de P(terme a un effet | données) : PROPOSÉ

BF d'inclusion **apparié** (analogue Type II) avec π = 1/2 par terme. Alternatives : tous modèles (cohérent, mais uniforme sur les modèles donne 0,6/0,2 à 2 facteurs et 0,737/0,316/0,053 à 3 facteurs, donc pas neutre) ; modèle plein moins T (analogue Type III).
Prix franc : les probabilités de termes ne sont pas les marginales d'une même loi jointe.

**Cette proposition est mise en défaut (bloquant).** Reproduit à la main (`scripts/agent_11_verif_croisement.R`) sur un croisement 2×2 pur : P(A) apparié = **0,215**, soit « Non 0,785 », alors que P(A:B) = 1,000 et que A agit manifestement, avec un signe opposé selon B. Le BF apparié A est de plus exactement le conditionnel de l'uniforme sur {B, A+B} (identité vérifiée) : « π = 1/2 » n'est donc pas un a priori distinct de celui qu'on prétendait rejeter.

**Nouvelle proposition.** La question d'Antoine est « ce facteur joue-t-il un rôle ? », pas « son effet principal moyen est-il non nul ? ». Répondre donc :
1. en première ligne, **P(A joue un rôle, seul ou en interaction)**, marginale d'une loi jointe unique et cohérente sur les modèles respectant la marginalité (type `all` de `.bay_inclusion`), avec la loi a priori annoncée ; dans le croisement ci-dessus elle vaut 1,000 ;
2. en seconde ligne, pour décrire, l'effet principal apparié, explicitement conditionnel et jamais présenté comme « Non, A n'a pas d'effet » ; l'effet de A selon les niveaux de B si l'interaction est crédible.

Loi jointe à arbitrer : uniforme sur les modèles (a priori 0,6 sur chaque facteur et 0,2 sur A:B à 2 facteurs) ou hérédité forte, chaque terme entrant à 1/2 si ses parents sont présents (A à 1/2, A:B à 1/8, A:B:C à 1/128 ; réglage KefiR, non publié tel quel). Sous cette dernière, P(A | pas de A:B) vaut 3/7, et non 1/2. Annoncée dans tous les cas. Non simulée à ce stade pour cette forme finale.

### D2, loi a priori sur les modèles : DÉCIDÉE le 9 novembre 2026, HOP(1,1)

**Décision : une proportion π_j ~ Beta(1,1) par ordre d'interaction, indépendantes, sous hérédité forte** : la HOP(1,1) de [TaylorRodriguez2016], « version hiérarchique » de la correction de multiplicité de [ScottBerger2010]. Pour un modèle M, p(M) = Π_j B(1 + s_j, 1 + e_j − s_j), où s_j est le nombre de termes d'ordre j dans M et e_j le nombre de termes d'ordre j dont tous les termes d'ordre inférieur sont dans M. **Loi construite, pas ajustée** : à π fixés, on tire chaque terme admissible ordre par ordre, ce qui donne une loi de somme 1 sur tout espace héréditaire ; l'intégrale sur les π garde cette propriété (contrôlé numériquement à 2 et 3 facteurs, et par la relecture du 9 novembre). A priori : A 1/2 ; A:B 1/6 ; A:B:C 1/32 ; aucun effet 1/3 à 2 facteurs, 1/4 à 3 facteurs ; **sachant A et B présents, A:B part à 1/2**. Programmée dans `R/sys_bay_factorial.R` (`.bay_fact_prior`), annoncée dans le bilan, avec en sensibilité 1/2 par terme [Chipman1996] et l'uniforme.

**Pourquoi pas la Beta commune (HUP), d'abord proposée.** Sachant A et B présents, elle met l'interaction à 3/4. Le scénario additif demandé par l'audit (`agent_49`, mêmes BF pour toutes les lois) la met en défaut :

| Cas (fausse interaction : P(A:B ou autre interaction) > 0,5) | 1/2 par terme | HUP, Beta commune | **HOP, Beta par ordre** | HOP(1,ch) | uniforme | aov p<0,05 |
|---|---:|---:|---:|---:|---:|---:|
| 2×2, n=10, A et B à f=0,5 | 0,043 | 0,413 | **0,057** | 0,043 | 0,057 | 0,043 |
| 2×2, n=25, A et B à f=0,5 | 0,073 | 0,480 | **0,077** | 0,073 | 0,077 | 0,047 |
| 2³, n=5, A, B, C à f=0,5 | 0,085 | 0,440 | **0,065** | 0,005 | 0,215 | 0,190 |
| 2×2 sous H0, n=10 : un terme > 0,5 | 0,163 | 0,120 | **0,120** | 0,040 | 0,227 | 0,160 |
| 2³ sous H0, n=10 : un terme > 0,5 | 0,180 | 0,045 | **0,050** | 0,020 | 0,265 | 0,315 |
| Puissance : 2×2, n=10, A:B réelle f=0,5, P(A:B)>0,5 | 0,397 | 0,503 | **0,443** | 0,307 | 0,497 | 0,863 |
| Puissance : 2³, n=5, A, B, C f=0,5, P(vrai)>0,5 | 0,858 | 0,888 | **0,938** | 0,808 | 0,982 | 0,862 |

La HOP garde la correction de multiplicité de la HUP sous H0 sans fabriquer d'interactions quand les effets sont additifs. HOP(1,ch) contrôle encore mieux les faux « oui », mais met un facteur seul à 1/3 (2 facteurs) ou 1/4 (3 facteurs) au lieu du 1/2 retenu à un facteur, et perd de la puissance. L'uniforme et « uniforme sur la taille » favorisent les effets principaux (0,6 à 0,74 a priori). **Coût de la HOP** : une interaction réelle à petit effectif est moins souvent détectée qu'avec la HUP (0,44 contre 0,50 ; aov 0,86, seuils non équivalents) ; même constat sur plans déséquilibrés (`agent_50 B`, VALIDATION §8).

**Restent à exécuter** (`agent_49` tourne) : 2³ à n = 10 additif, 2³ avec A:B réelle ou toutes les interactions doubles (fausse A:B:C), A:B:C seule. Ils seront ajoutés ici ; la décision ne sera revue que s'ils la contredisent.

**Libellé par terme**, tranché avec le libellé non causal de l'audit : « y diffère selon A, seul ou en interaction ».

*Historique : ce qui suit retrace les étapes antérieures (réglage 1/2, 1/8 du 6 novembre, puis HUP du 7 novembre). Gardé pour mémoire.*

**Candidate : π ~ Beta(1,1) intégrée sur les termes éligibles, sous hérédité forte.** Un terme n'est éligible que si tous ses parents sont dans le modèle ([Chipman1996]). Pour un modèle à s termes parmi e éligibles, p(M) = B(1+s, 1+e−s) ([ScottBerger2010], éq. 8-9, lues dans le PDF). La combinaison avec l'hérédité est propre à KefiR. A priori obtenus : à 2 facteurs, A 1/2, A:B 1/4, nul 1/3 ; à 3 facteurs, A 1/2, A:B 1/4, A:B:C 1/8, nul 1/4. Rien n'est réglé par l'utilisateur ni estimé sur les données.

**Comparaison** (`scripts/agent_47_loi_sur_modeles.R`, N=300 jeux par cas, vraisemblances `anovaBF` en loi normale à variance commune, mêmes BF pour les trois lois) :

| Cas | Critère | 1/2, 1/8 | Beta(1,1) | Uniforme | aov |
|---|---|---:|---:|---:|---:|
| 2³ sous H0, n=5 par case | au moins un terme P>0,5 / P>0,95 | 0,217 / 0,010 | **0,093 / 0,000** | 0,430 / 0,023 | p<0,05 : 0,270 |
| 2³ sous H0, n=10 | idem | 0,163 / 0,003 | **0,070 / 0,000** | 0,277 / 0,007 | 0,317 |
| 2×2 sous H0, n=10 | idem | 0,143 / 0,007 | **0,090 / 0,010** | 0,207 / 0,010 | 0,130 |
| 2×2 sous H0, n=25 | idem | 0,090 / 0,000 | **0,053 / 0,000** | 0,147 / 0,003 | 0,123 |
| Interaction réelle, n=10, f=0,25 | P(A:B)>0,5 / >0,95 | 0,057 / 0,003 | **0,107 / 0,007** | 0,090 / 0,007 | p<0,05 : 0,337 |
| Interaction réelle, n=10, f=0,50 | idem | 0,417 / 0,123 | **0,553 / 0,167** | 0,543 / 0,160 | 0,877 |
| Interaction réelle, n=25, f=0,25 | idem | 0,200 / 0,043 | **0,267 / 0,060** | 0,263 / 0,060 | 0,713 |
| Interaction réelle, n=25, f=0,50 | idem | 0,930 / 0,713 | **0,960 / 0,763** | 0,963 / 0,757 | 0,997 |

Sous H0, Beta(1,1) divise par deux environ les faux « oui » du réglage 1/2, 1/8. Face à une interaction réelle, elle la détecte aussi bien que l'uniforme et mieux que 1/2, 1/8 (+5 à +14 points). La raison : l'interaction part de 1/4 et non de 1/8, et la correction de multiplicité vient de l'intégration sur π, pas d'un a priori bas sur chaque terme. **Limites** : un seul plan à 3 facteurs, sous H0 ; aucun effet réel à 3 facteurs simulé ; loi normale à variance commune seulement ; l'écart à 2×2 sous H0 à P>0,95 (0,010 contre 0,007) est dans le bruit. À recontrôler sur le moteur multi-facteurs de m.bay quand il existera.

**Libellé par terme**, tranché avec le libellé non causal de l'audit : « y diffère selon A, seul ou en interaction ».

*Ce qui suit concerne l'ancien réglage 1/2, 1/8. Gardé pour mémoire.*

Définition à 3 facteurs retenue par cohérence (chaque terme entre à 1/2 si tous ses parents sont présents) : A:B à 1/8, A:B:C à 1/128. Ce prolongement est une lecture, pas ses mots : il n'a parlé que de facteur et d'interaction.

**La simulation ne contredit pas ce choix, mais il a un coût.**
- Plan 2×2×2 sous H0, N=40, 300 jeux (`agent_22`) : au moins un terme à P>0,5 dans 25 % [21 ; 31] des jeux, contre 49 % sous l'uniforme ; P>0,95 jamais, [0 ; 1,2 %]. Une partie de ce gain vient de l'a priori lui-même : ce n'est pas une preuve de justesse.
- Coût, interaction pure réelle, 2×2, N=200 par case (`agent_31`) :

| Cas | P(A:B)>0,5, hérédité | P(A:B)>0,5, uniforme | P(A:B)>0,95, hérédité | P(A:B)>0,95, uniforme | aov p<0,05 |
|---|---:|---:|---:|---:|---:|
| n=10 par case, f=0,25 | 0,045 | 0,115 | 0,000 | 0,010 | 0,340 |
| n=10 par case, f=0,50 | 0,375 | 0,510 | 0,125 | 0,185 | 0,855 |
| n=25 par case, f=0,25 | 0,180 | 0,225 | 0,025 | 0,025 | 0,715 |
| n=25 par case, f=0,50 | 0,930 | 0,955 | 0,720 | 0,790 | 1,000 |
| n=25 par case, f=0 | 0,000 | 0,005 | 0,000 | 0,000 | 0,045 |

Une interaction réelle à petit effectif est donc plus souvent manquée qu'avec l'uniforme : écart net à n=10, faible à n=25 et f=0,50. Le gros de l'écart avec aov existe sous les deux lois, à seuils non équivalents. Limites : lancement interrompu à 500 s, n=50 à f>0 non exécuté ; 3 facteurs à effet réel non simulés, alors que A:B:C part de 1/128 ; forme finale du bilan non simulée. Les probabilités des deux lois sont calculées sur les mêmes BF.

## D3. H0 ponctuelle ou ROPE : PROPOSÉ

H0 ponctuelle en première ligne, ROPE toujours affiché en seconde, avec sa masse a priori. Largeur par défaut |δ|<0,1 (f<0,05), convention remplaçable par une valeur métier.

**Mis en défaut (bloquant, relecteur, non reproduit par moi).** Pour un terme à 3 niveaux ou plus, f est positif par construction et son postérieur est biaisé vers le haut : vérité nulle, a=5, n=100, médiane de f à 0,105, soit deux fois le seuil ; la ligne ROPE recopie alors presque la ligne H0. **Proposition** : ROPE seulement sur des quantités à 1 ddl (paires, contrastes, pentes), pas de ROPE sur un terme à plusieurs ddl tant qu'une définition valide n'est pas trouvée. La masse a priori de 0,545 ne vaut que pour a=2.

## D4. Probabilité a priori des hypothèses et sensibilité : DÉCIDÉ, valeurs revues le 7 novembre

P(H0)=P(H1)=1/2 par défaut, annoncée dans le bilan. Grille r × cotes a priori (r/2, r, 2r ; 1:2, 1:1, 2:1) et Dirichlet u/3, 3u. Si la plage traverse 1/2, le bilan dit « les données ne tranchent pas » ; il donne aussi les cotes qu'il faudrait pour mettre les deux hypothèses à égalité, c'est-à-dire le facteur de Bayes.

**Changement.** La grille d'origine 1:4 à 4:1 a été abandonnée : l'hypothèse en tête y changeait dans 98 % des jeux sans effet, jusqu'à n=100 par groupe (`agent_41` C, premier lancement). Le verdict ne disait plus rien. **La nouvelle grille est un jugement, pas un résultat** : à 1:2, un BF de 2 suffit à trancher. Chiffres en `VALIDATION.md` §7.

## D5. Échelles a priori : DÉCIDÉ pour un facteur, r = 1/2 fixe

Celles de BayesFactor, pour rester comparable à JASP et à la littérature : facteur fixe r=1/2 ; deux groupes r=√2/2 sur δ (même a priori) ; covariable r=√2/4 ; effet aléatoire r=1.

**Hyper-a priori, testé le 7 novembre.** Loi uniforme sur r/2, r, 2r, fixée d'avance : P>0,5 et P>0,95 inchangés à ±0,033 sur 66 scénarios, P<0,05 jusqu'à 0,12 d'écart (`agent_41` B). Elle ne change donc presque rien : on garde le réglage le plus simple, avec r/2 et 2r dans la sensibilité. L'estimer sur les données (Bayes empirique) reviendrait à maximiser le BF : exclu, D6. Sur π : à un terme, un hyper-a priori se réduit à sa moyenne, donc sans objet. À plusieurs termes, il fait la correction de multiplicité, [ScottBerger2010] : à trancher au jalon suivant. **L'utilisateur n'a rien à régler.**

## D6. Pas d'aiguillage par test préliminaire : ARBITRÉE (orientation) ; voie à un facteur DÉCIDÉE

Contrairement à `m.test()`, pas de bascule automatique vers une méthode selon une p-value de normalité, de variance, etc. Le gaussien est toujours calculé ; les diagnostics deviennent des alertes ; la voie robuste s'y ajoute et un désaccord s'affiche. Conséquence à accepter : `m.bay()` ne reproduit **pas** l'arbre de décision de `m.test()`, et pas de « Kruskal bayésien » pour k>2 (aucune méthode publiée trouvée ; cotes proportionnelles proposées sous leur vrai nom).

**Mis en défaut (bloquant).** La voie gaussienne seule n'est pas défendable en cas d'hétéroscédasticité : H0 vraie, n=(5,5,40), sd=(4,4,1), **P>0,95 dans 16 % des jeux** (IC 95 % 12 à 21 %, N=300, `scripts/agent_12_verif_heterosc.R` ; 25,5 % chez le relecteur sur un autre tirage, écart non élucidé, tailles de simulation à comparer) ; aov 46 % de p<0,05, Welch 7 %. Et la voie de secours prévue ne se calcule pas sans Stan, et rstanarm n'a pas de variances hétérogènes prêtes à l'emploi. **Proposition** : voie interne sans compilateur pour les variances par groupe, et refuser de répondre « Oui » sans réserve quand le diagnostic d'hétéroscédasticité est net, tant que cette voie n'est pas validée. À ce jour, m.bay ne doit donc pas être publié sur ce cas.

**État au 6 novembre.** Taux corrigé : 20,0 % à P>0,95 (N=3000), le 16 % ci-dessus était un tirage bas. La voie `.bay_bf_oneway_hetero` existe, à un facteur seulement : 0,7 % à P>0,95 sur ce plan, à calcul vérifié. Ce n'est pas encore une solution : a priori propre à KefiR, τ à choisir, puissance mesurée sur N=16-17 seulement, et surtout **deux défauts bloquants reproduits** : un petit groupe éloigné pousse le BF vers H0 (n=(10,10,10,10,3), écart 50σ : log BF -6,9, Welch p=3e-15) et un groupe dispersé sans rapport fait basculer le BF, sans plancher (`RELECTURE_ADVERSE_2026-11-06.md`, B1, B2 ; `agent_30`). **Cette loi n'est donc pas proposée.** Pistes non programmées : un g par groupe, un τ propre à chaque contraste ; à défaut, refuser ou signaler min(n_i) ≤ k. Multi-facteurs non couverts.

**ARBITRÉE par Antoine, 6 novembre 2026 : orientation retenue.** Lois d'échantillonnage comparées selon les données, réponses moyennées sur leurs probabilités a posteriori ; a priori d'effet et d'hypothèse fixés d'avance, annoncés, avec sensibilité. Mots transmis par la session « KefiR officiel ».

**Faisabilité à k > 2 : établie pour un facteur (7 novembre).** Voie `R/sys_bay_laws.R`, huit modèles sans compilateur. Loi propre sur les rapports d'échelle : Dirichlet(1,5) sur les proportions de précision, d'après [Dablander2024] (forme) et [Maier2025] (valeur à k = 2). Détails dans METHODOLOGIE §7 et chiffres dans VALIDATION §7. **Réglage propre à KefiR à k ≥ 3** : aucune source consultée ne combine ces éléments pour tester des moyennes.

**Effet à deux formes, DÉCIDÉ le 7 novembre (k ≥ 3).** H1 mêle à parts égales « tous les groupes peuvent différer » (Zellner-Siow) et « un seul groupe s'écarte » (échelle r√2). Raison : avec la première seule, un petit groupe éloigné parmi 20 laissait H0 en tête à 83 %, Welch p=1e-111 (`RELECTURE_ADVERSE_2026-11-07.md`, B1). **Alternatives écartées** : un g par groupe, qui multiplie les dimensions et casse la lecture « toutes les différences » ; refuser ou seulement avertir, qui laisse une réponse fausse. **Prix** : à k ≥ 3, le BF global n'est plus celui d'`anovaBF`, même en gaussien ; le bilan donne la valeur « première forme seule ». Le poids 1/2 est un réglage, pas un résultat. Sur les 66 scénarios de `agent_41` B, P>0,5 et P>0,95 bougent de 0,07 au plus, le plus à k=5 quand un seul groupe s'écarte (P>0,95 : 0,537 contre 0,467).

**Limites à connaître.** Formes différentes à moyennes égales : la loi t compare des centres, donc un groupe asymétrique peut lire comme un effet (P>0,95 dans 67 % des jeux à n=400, `agent_43`) ; le bilan avertit, le calcul n'est pas changé. Valeurs répétées : les lois dont l'intégrale diverge sont écartées avec un avis (règles en tête de `.bay_laws_allowed`, dérivées à la main). Variances étagées sur 4 ordres de grandeur : voir VALIDATION §7.

Lecture initiale de cette orientation, pour mémoire :

**Orientation d'Antoine, 6 novembre** (transmise par sa session « KefiR officiel ») : « à toi de faire une fonction qui teste le plus rationnellement les lois les plus probables découlant des données ». Lecture de cette session, **à lui confirmer**, à étudier, rien n'est exécuté :
1. **Lois d'échantillonnage** (normale ou t, variance commune ou par groupe, plus tard discrètes) : `m.bay` les compare et moyenne ses réponses sur leurs probabilités a posteriori, sans aiguillage par tests préalables. Cela remplace l'hypothèse « le gaussien toujours, le robuste en plus, sans gagnant » et le niveau 3 de METHODOLOGIE §7 devient la voie étudiée par défaut.
2. **Condition technique** : la mesure plate sur les rapports d'échelle rend BF(hétéro contre homo) indéfini. Il faut une loi propre sur ces rapports. Piste : Maier et al. 2025, à 2 groupes (ρ ~ Beta(1,5 ; 1,5) d'après une extraction abîmée, à vérifier sur le PDF) ; le passage à k>2 est à établir (Dablander et al. 2024 emploie une Dirichlet sur les rapports, pour tester les variances). Calibration et puissance à mesurer. Le modèle devra aussi passer B1 et B2 de `RELECTURE_ADVERSE_2026-11-06.md`, que la moyenne de modèles ne corrige pas d'elle-même.
3. **Limite** : les a priori sur les paramètres d'effet (r) et sur les hypothèses ne sont jamais réglés sur les mêmes données. Maximiser le BF sur r ne peut que le gonfler. Ils restent fixés d'avance, annoncés, avec sensibilité. Seul le poids des lois d'échantillonnage est appris des données, par leur probabilité a posteriori sous un a priori lui aussi annoncé.
4. **Sortie de secours** : si cette voie n'est pas calibrée pour k>2 ou est trop lente, le journal le dit franchement.

## D7. MANOVA : PROPOSÉ

Forme analytique Normal-Inverse-Wishart, g par terme, Jeffreys sur Σ. **Réglage propre à KefiR, pas une méthode publiée telle quelle.** Algèbre contrôlée (réduction à `anovaBF` à d=1 ; invariance) ; calibration non faite. Le banc d'essai a mesuré un autre réglage (IW propre, S0 ∝ cov(Y), lambda fixe) : très conservateur, P>0,5 dans 43 à 72 % des jeux à effet connu contre 89,5 % de p<0,05 en manova. Le réglage retenu devra faire mieux par simulation avant d'être annoncé.

**Mis en défaut (bloquant, relecteur, non reproduit par moi).** L'algèbre tient, mais l'a priori impose sur la distance de Mahalanobis une échelle qui croît en √d (médiane 0,50 à 2,69) : effet réel f=0,41 sur une VD parmi 20, n=30, P médiane 0,016 alors que Pillai détecte 37 %. **Proposition** : ne pas livrer la MANOVA avant d'avoir un r qui dépend de d, calibré par simulation.

## D13. Multiplicité entre les termes : PROPOSÉ (ajouté après relecture)

**Bloquant, recoupé sur `anovaBF` le 6 novembre.** Plan 2³ sous H0, au moins un « Oui » (P>0,5) sur 7 termes dans 58 à 74 % des jeux selon N ; 25 % sous l'hérédité de D2 (N=40). π fixe par terme n'ajuste rien, et le seuil « Beta-binomial au-delà de 3 facteurs » tombe précisément sur le plan d'expériences type de l'industriel. **Proposition** : ne jamais utiliser P>0,5 comme verdict ; annoncer le nombre de termes testés et la probabilité que le modèle soit le nul ; Beta-binomial dès 3 facteurs, à valider.

## D8. Plans que m.test ne couvre pas : PROPOSÉ

Covariable seule, ANCOVA appariée, MANOVA répétée, MANCOVA, pentes aléatoires réelles : `m.test()` les refuse ou efface la pente. Le bayésien peut couvrir certaines de ces voies. Proposition : pentes aléatoires et covariable seule inclus ; les autres plus tard.

## D9. Reprise de la normalisation des entrées : PROPOSÉ

Toute la normalisation est en ligne dans `m.test()` (lignes 262-1031) : rien n'est réutilisable en bloc. Options :
- **A (proposée)** : extraire une fonction commune `.mtest_normalize_input()` et faire appeler `m.test()` et `m.bay()` dessus. Seule façon d'accepter durablement « exactement les mêmes entrées », mais touche `m.test()`, donc demande les tests de non-régression et ton accord.
- **B** : recopier dans `m.bay()` sans toucher `m.test()`. Sans risque immédiat, mais divergence garantie, et les bugs d'entrée s'y recopient.

**État au 7 novembre : option B prolongée, PROPOSÉ et réversible.** `R/sys_bay_input.R` réutilise les fonctions internes de `m.test()` sans les changer (`.strip_factor_calls`, `.msg`, `.exit`, `discret.test`). `m.test()` n'est pas touché. Les formes qui échouent dans `m.test()` (`formula=` nommé, `g` data.frame sans `data`) marchent ici. Divergence possible, à surveiller : les deux lectures de formule existent désormais.

## D10. Retour : ACQUIS sur le principe

Même structure que `m.test()` quand elle a un sens, classe distincte. `$effects` existe même à un facteur (une ligne), mêmes noms de termes, colonnes `term`, `bf10`, `log_bf`, `p_h1`, `p_h0`, `prior_h1`, `rscale`, `error_pct`, `analogue`. Aucune colonne ne se nomme `p.value`. Mode `return=FALSE` : valeur à décider (le « modèle complet contre nul » plutôt que le premier facteur, comme `$global_pvalue`).
Lettres de groupes : pas par défaut (tableau des paires à la place).

## D11. Post-hoc : PROPOSÉ

Partitions de moyennes jusqu'à k=5, puis paires à cotes a priori ajustées. Dérivation propre pour le second cas, à présenter comme réglage de KefiR.

## D12. Socle commun : ACQUIS

Fonctions internes de `R/sys_bay_core.R`, réutilisables par `cor.bay` et `g.bay`, sans développer ces deux fonctions avant ta demande.
