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

**Cas à 3 facteurs exécutés le 9 novembre** (`agent_49`, fin du calcul, N = 200 par cas, mêmes BF pour toutes les lois). On compte la proportion de jeux où au moins un terme faux dépasse P = 0,5, et la puissance P(vrai) > 0,5.

| Cas 2³ | 1/2 par terme | HUP | **HOP** | HOP(1,ch) | uniforme | aov p<0,05 |
|---|---:|---:|---:|---:|---:|---:|
| n=10, A, B, C à f=0,5 : un faux > 0,5 | 0,220 | 0,605 | **0,120** | 0,060 | 0,240 | 0,185 |
| n=5, A, B, C, A:B : un faux > 0,5 | 0,110 | 0,625 | **0,205** | 0,045 | 0,205 | 0,135 |
| n=10, A, B, C, A:B : un faux > 0,5 | 0,130 | 0,800 | **0,270** | 0,070 | 0,215 | 0,130 |
| n=5, toutes les interactions doubles : fausse A:B:C > 0,5 | 0,020 | 0,840 | **0,030** | 0,015 | 0,030 | 0,025 |
| n=10, toutes les interactions doubles : fausse A:B:C > 0,5 | 0,050 | 1,000 | **0,060** | 0,050 | 0,055 | 0,030 |
| n=10, A:B:C seule (sans effets principaux) : fausse A:B, A:C ou B:C > 0,5 | 0,220 | 0,495 | **0,340** | 0,100 | 0,505 | – |
| Même cas, puissance sur A:B:C | 0,215 | 0,490 | **0,325** | 0,100 | 0,490 | 0,995 |

**Au-delà de « faible »** (P > 0,75, `agent_49b`, mêmes BF), la part de jeux avec au moins une fausse interaction devient, sous HOP et contre 1/2 par terme :
- A, B, C réels : 5,0 % contre 6,5 % ;
- A:B réelle, n=5 : 5,5 % contre 3,0 % ;
- A:B réelle, n=10 : 6,5 % contre 4,0 % ;
- A:B:C seule : 23,5 % contre 17 %.

À P > 0,95, la part est de 0,5 % dans les trois premiers cas, et de 14,5 % pour A:B:C seule.

**Libellé corrigé le 9 novembre.** Une première version comptait comme faux, dans le cas A:B:C seule, A, B et C. Or la réponse par terme porte sur « seul ou en interaction » : A agit donc par A:B:C. Seuls A:B, A:C et B:C y sont faux. De même, dans le 2×2 « A:B seule » (`agent_50 B`), A et B ne sont pas des faux termes.

**Lecture.**
- Additif et toutes interactions doubles : la HOP reste au niveau de 1/2 par terme ou mieux. La HUP fabrique des interactions (60 à 100 %).
- **Une interaction réelle parmi plusieurs possibles** (A:B vraie, A:C et B:C fausses) : **la HOP fait moins bien** que 1/2 par terme et qu'aov à P > 0,5, avec 20 à 27 % de faux « oui » contre 11 à 13 %. L'écart se réduit à 5,5–6,5 % contre 3–4 % au-delà de « faible ». C'est la contrepartie de la proportion partagée par ordre : une interaction double nette relève π₂, donc l'a priori des autres interactions doubles. Ce coût est connu de la construction ([TaylorRodriguez2016] : échangeabilité à l'intérieur d'un ordre). Il n'existe pas sous H0 ni en additif, où la HOP fait mieux que 1/2 par terme (5 % contre 18 % sous H0 à n=10 ; 12 % contre 22 % en additif).
- **A:B:C seule, sans effet principal** : aucune loi héréditaire ne la représente, puisque A:B:C n'entre qu'avec A, B, C, A:B, A:C et B:C. Sa détection tombe à 0,33 (aov 0,995), et les interactions doubles, fausses, reçoivent un « oui » dans 34 % des jeux (23,5 % au-delà de « faible »). Ce cas sort de l'espace héréditaire par construction ([Chipman1996]). Il faut le dire dans la documentation plutôt que d'abandonner l'hérédité.

**Décision maintenue : HOP(1,1).** Le bilan de chiffres reste favorable, sous H0, en additif et en puissance additive. Les deux coûts sont documentés :
- sur-détection d'une seconde interaction quand une première est réelle ;
- aveuglement à une interaction d'ordre supérieur sans effets inférieurs.

1/2 par terme est déjà dans la sensibilité : un terme que la loi sur les modèles fait changer de côté est dit « non tranché ». La sensibilité ne couvre pas pour autant tout ce coût, car les fausses interactions à 20–27 % peuvent rester du même côté sous les trois lois. Non mesuré : la part de ces faux « oui » que la sensibilité requalifie en « non tranché ». **Point ouvert pour la session parente** : faut-il tester une HOP dont le paramètre de second ordre est moins partagé, ou une loi hybride ? Ni l'une ni l'autre n'est programmée.

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

## D11. Post-hoc : DÉCIDÉ le 9 novembre (un facteur ; plusieurs facteurs à suivre)

Proposition initiale, gardée pour mémoire : partitions de moyennes jusqu'à k=5, puis paires à cotes a priori ajustées.

**Décidé, selon le principe d'Antoine (paquets établis d'abord).**
- **Quels groupes diffèrent : partitions d'égalité.** Chaque façon de regrouper les k groupes en blocs de moyennes égales (« a = c < b ») est un modèle. Sa vraisemblance marginale vient de `lmBF(y ~ h)` [BayesFactor], où h est le facteur fusionné, ou du calcul interne identique à `lmBF` pour un terme seul (VALIDATION §8).
  - Le nombre de partitions suit les nombres de Bell : 5, 15, 52, 203, 877, 4 140 pour k = 3 à 8. Toutes sont calculées jusqu'à 8 groupes. Durées mesurées : 0,3 s à 6 groupes, 0,8 s à 7, 6,4 s à 8, 32 s à 9, d'où le seuil.
  - **Au-delà de 8 groupes (ajout du 9 novembre, après l'audit de la session parente).** Les partitions sont explorées par l'échantillonneur de Gibbs de [GopalanBerry1998], la méthode même de l'article dont vient l'a priori.
    - Chaque groupe est déplacé à son tour vers l'un des blocs des autres ou seul, avec le même a priori et les mêmes vraisemblances marginales (mises en cache).
    - 200 tours de chauffe, puis max(2 000, draws/2) tours.
    - Les probabilités de paires sont des fréquences de visite. Leur erreur de tirage est estimée par moyennes de lots (20 lots) et annoncée dans le récit.
    - La sensibilité à α repondère les mêmes tirages, avec la constante de normalisation d'Ewens exacte : Γ(α+k)/Γ(α) − α(k−1)!.
  - Écartées : les comparaisons par paires avec a priori de multiplicité ([Westfall1997], [ScottBerger2006]). Elles traitent chaque paire à part et peuvent donner a = b, b = c, a ≠ c, ce que la loi unique sur les partitions exclut. Garder la même loi de 3 à 12 groupes rend aussi les sorties comparables.
  - Accord avec l'énumération exacte, en forçant le tirage à 6 et 7 groupes : 1,4 point au plus sur les paires (3 000 tours). Durées : 20 s à 9 groupes, 40 s à 10, 93 s à 12.
  - Au-delà de 15 paires, le récit et le bilan ne détaillent que les paires qui penchent vers « différents » et résument les autres en une ligne. Toutes restent dans `$posthoc$pairs`.
  - La probabilité que deux groupes diffèrent est la somme des partitions qui les séparent. Une seule loi porte donc sur toutes les paires : pas de multiplicité ajoutée paire par paire, pas de résultats incohérents du type a = b, b = c, a ≠ c.
- **A priori sur les partitions.**
  - « Aucune différence » garde l'a priori de la réponse globale, 1 − `prior_h1`.
  - Les autres partitions se partagent `prior_h1` selon la loi d'Ewens (α = 1), celle du processus de Dirichlet : p ∝ Π (taille du bloc − 1)! ([GopalanBerry1998]).
  - Cette loi favorise les regroupements, ce qui corrige le nombre de paires. Elle ne dépend pas des données.
  - En sensibilité, α = 1/2 et α = 2.
- **De combien : taille d'effet et intervalle.**
  - Échantillons de la loi a posteriori du modèle « tous les groupes peuvent différer » : `posterior()` [BayesFactor].
  - Pour chaque paire : différence brute et différence standardisée δ = (m_j − m_i)/σ, médianes et intervalle de crédibilité à 95 % (HDI, `hdi()` [bayestestR]).
  - Les deux sont conditionnels à l'existence d'une différence (rétrécis vers zéro par l'a priori de Cauchy).
- **Deux chiffres pour un écart, nommés (ajout du 9 novembre).** La conclusion donne l'**écart observé** : différence des moyennes, ou des médianes, de l'échantillon. Le post-hoc donne l'**écart estimé** : médiane a posteriori, ramenée vers zéro par l'a priori. Le récit dit pourquoi le second est un peu plus petit.
- **Limite annoncée.** Le post-hoc suppose la loi normale à variance commune. Si l'étape 4 trouve des variances inégales probables ou des queues lourdes, le récit le dit à l'étape 6, et la conclusion le répète sous les paires (« À lire avec réserve »). La probabilité « aucune différence » du post-hoc peut alors différer de la réponse globale, qui moyenne les formes du bruit.
- **Paquets.** BayesFactor (déjà en Imports). bayestestR ≥ 0.13 passe en Imports : il est léger (insight, datawizard) et sert la HDI, référence [bayestestR].
- **Contrôles prévus.**
  - Probabilités de paires recalculées à la main sur un jeu, par un agent indépendant.
  - δ comparé au d de Cohen et à son intervalle fréquentiste à grand n.
  - Simulation : paires faussement « différentes » sous H0 et sous « un groupe s'écarte », puissance, à côté de Tukey HSD.

## D12. Socle commun : ACQUIS

Fonctions internes de `R/sys_bay_core.R`, réutilisables par `cor.bay` et `g.bay`, sans développer ces deux fonctions avant ta demande.

## D14. Bilan en récit : ACQUIS (demande d'Antoine, 9 novembre ; mise en œuvre décidée ici)

Demande, ses mots : « m.bay() doit faire comme m.test() proposer un bilan décisionnel adaptatif qui s'adapte en direct, présente ses éventuels changements de stratégies pour arriver enfin à la fin sur, au lieu de p-values cette fois, des hypothèses argumentés et associés à un indice de confiance. »

**Décidé.**
- `verbose = TRUE` raconte six étapes numérotées, imprimées au fil du calcul par `.vbse()` : plan, contrôles des données, a priori (ou forme du bruit puis a priori à plusieurs facteurs), formes du bruit, facteurs de Bayes et sensibilité, conclusion. « ==> » est un constat, « --> » la décision qui en découle, et la fonction et le paquet sont entre crochets.
- Chaque ligne ne dit que ce qui a été calculé dans l'appel : les avertissements sont dits à l'étape où ils naissent, pas répétés à la fin.
- Deux mentions distinctes, après la relecture par le lecteur type (`relecture/lecteur_recit/RELECTURE.md`, point 5) :
  - « Changement de stratégie : » quand le calcul change : loi écartée pour ex-aequo (intégrale divergente) ou pour échantillonnage trop pauvre ;
  - « Prudence ajoutée : » quand seule la conclusion se retient : contradiction entre lois, qui plafonne le soutien ; sensibilité qui peut renverser l'hypothèse en tête ; terme qui change de côté et devient « non tranché ».
- La conclusion est le bilan court, en retrait sous l'étape 6, sans titre, sans point de départ ni avertissements (dits plus haut). Chaque hypothèse porte une ligne « Repose sur » ou « Contre elle » :
  - un facteur : le plus grand écart observé en écarts-types, et le facteur de Bayes en clair (« 5 fois plus probable avec une différence que sans ») ;
  - plusieurs facteurs : pour la réponse globale, le terme le plus probable ; pour chaque terme, les données qui multiplient ou divisent ses cotes.
- `verbose = FALSE` n'imprime rien. L'objet garde le récit (`$story`, en et fr) et `print(x, story = TRUE)` le rejoue à l'identique.

**L'indice de confiance** est la probabilité a posteriori suivie du qualificatif de soutien (aucun, faible, modéré, fort, très fort ; seuils de Kass et Raftery). Le mot « confiance » n'est pas imprimé, pour qu'on ne le prenne pas pour un intervalle de confiance ou un niveau 1 − α ; le bilan dit « soutien ». Ce choix est réversible si Antoine tient au mot.

**Écart à `m.test()`.** `.vbse()` ne garde que son texte français quand on lui passe une chaîne (`if (length(ang)==1) ang <- fr`) : en session anglaise, `m.test()` raconte en français. `m.bay()` choisit la langue par `.msg()` avant d'appeler `.vbse()`. Le défaut de `.vbse()` est signalé, pas corrigé (il sert `m.test()`).

**Ajouts du 9 novembre, après l'audit de la session parente.**
- **Constats de l'étape 4.** La part des lois à variances inégales est annoncée comme un constat, l'équivalent du contrôle de Levene ou Bartlett dans `m.test()` :
  - « Variances inégales probables » à 75 % ou plus, avec le groupe le plus dispersé et son écart-type ;
  - « Variances égales probables » à 25 % ou moins ;
  - « les données ne tranchent pas » entre les deux.
  Les queues lourdes sont dites par la décision qui suit (« l'emportent », « les plus probables », « les formes normales… »).
- **Qualificatif quand les données apportent peu** (règle revue après la seconde relecture). Si le facteur de Bayes de la réponse en tête est entre 1/3 et 3 (« ne vaut guère plus qu'une mention », Kass et Raftery 1995) :
  - le qualificatif est « non tranché par les données » ;
  - le verdict est « indécis », dans un sens comme dans l'autre ;
  - la probabilité reste affichée.
  Pour un terme, s'y ajoute le cas où les données penchent dans le sens opposé à la réponse : il est « non tranché ; les données penchent vers le oui, l'a priori de l'autre côté ».
  Raisons :
  - un terme d'interaction part de 1/6, et un « non » à 93 % peut ne devoir presque rien aux données ;
  - la première version (« faible, venu de l'a priori ») était fausse quand l'a priori était de 1/2 et que les données avaient fait monter la probabilité ;
  - deux situations symétriques (facteur de Bayes 0,9 et 1,6) recevaient des verdicts différents.
  Testé dans `test-bay-story.R`.
- **Écarts-types proches mais dispersions inégales probables** (rapport < 1,3) : le constat dit que les groupes diffèrent plutôt par leur forme (asymétrie, valeurs éloignées) que par leur écart-type. La loi à variances inégales réagit aussi aux différences de forme.
- **Langue.** KefiR choisit la langue par `LC_CTYPE`, via `.msg()`. `LANGUAGE=en` seul laisse donc la sortie en français, pour `m.bay` comme pour `m.test`. En locale anglaise réelle (`Sys.setlocale("LC_CTYPE", "English_United States.utf8")`), le récit de `m.bay` est entièrement en anglais (testé). Le défaut propre à `.vbse()`, qui ne garde que le français dans `m.test`, reste signalé à part.

## D15. Données appariées et mesures répétées : DÉCIDÉ le 9 novembre

Ordre fixé par la session parente : appariés et mesures répétées, puis ANCOVA, puis le discret. Le principe d'Antoine s'applique : on appelle un paquet établi, ici BayesFactor, déjà en Imports.

**Plans couverts.** Ils sont reconnus par `.bay_parse_input`, comme dans `m.test()` : `y ~ temps | id`, `Error(id/temps)`, `id =`, `paired = TRUE` sans id (l'ordre des lignes apparie).
- **Apparié** : un facteur à 2 niveaux, une mesure par sujet et par niveau.
- **Mesures répétées** : un facteur intra à 3 niveaux ou plus, ou plusieurs facteurs, intra et inter (plan en tranches).
- **Refusés** avec un message clair :
  - plusieurs mesures d'un même sujet dans une même condition (modèle mixte) ;
  - pentes aléatoires ;
  - covariables avec sujets.

**Calcul.**
- **Apparié** : `ttestBF(x, y, paired = TRUE, rscale = r)` [BayesFactor ; Rouder et al. 2009]. C'est un calcul exact, par intégration, sans tirage, sur les différences au sein des sujets. Les sujets qui n'ont qu'une des deux mesures sont retirés, et le récit le dit.
  - Le même r = 0,5 sert pour l'apparié et pour l'effet sujet aléatoire. Sur 4 jeux, `ttestBF(paired, r = 0,5)` et `lmBF(y ~ temps + id, whichRandom = "id", rscaleFixed = 0,5) / lmBF(y ~ id)` diffèrent de 0,02 à 0,21 en log-BF ; avec r = 0,707 (défaut de `ttestBF`), l'écart est de 0,10 à 0,26. Les réponses à 2 temps et à 3 temps restent ainsi sur la même échelle.
- **Mesures répétées** : `lmBF(... + id, whichRandom = "id")` [BayesFactor ; Rouder et al. 2012] pour chaque modèle de l'espace héréditaire. Chaque facteur de Bayes se lit contre le modèle « sujets seuls », de sorte que l'effet sujet est présent dans tous les modèles, « aucune différence » compris.
  - La précision est mesurée en recalculant avec une autre graine, comme à plusieurs facteurs.
- **A priori.**
  - Un terme seul (apparié, ou un facteur intra) : `prior_h1`, comme à un facteur entre groupes. La sensibilité fait varier r (r/2 à 2r) et les chances de départ (1 contre 2 à 2 contre 1).
  - Plusieurs termes : HOP, comme D2.
- **Forme du bruit** : loi normale seule. Queues lourdes et variances inégales entre conditions ne sont pas validées avec des sujets, et le récit le dit. Les valeurs éloignées sont cherchées sur les écarts au niveau propre de chaque sujet.

**Post-hoc (D11 étendu).**
- **Apparié** : écart moyen (second niveau moins premier) et dz = écart / écart-type des différences, médianes et HDI à 95 %, tirés de `posterior(ttestBF)` [BayesFactor] et de `hdi()` [bayestestR].
- **Un facteur intra à 3 niveaux ou plus** : partitions d'égalité, comme D11, avec le BF de chaque regroupement h donné par `lmBF(y ~ h + id, whichRandom)` / `lmBF(y ~ id)`.
  - Le calcul est exact jusqu'à 5 niveaux (52 regroupements, 0,12 s chacun). Au-delà, on passe au tirage de Gibbs.
  - Les écarts viennent de `posterior(lmBF(y ~ temps + id))`. δ est rapporté à l'écart-type résiduel au sein des sujets, ce que le récit nomme.
- **Plusieurs facteurs avec sujets** : pas de post-hoc pour l'instant, comme à plusieurs facteurs sans sujets.

**Écart à `m.test()`.**
- `m.test()` contrôle la sphéricité (Mauchly) et corrige les degrés de liberté. Le modèle à effet sujet aléatoire suppose la symétrie composée : c'est la même hypothèse que l'ANOVA à mesures répétées non corrigée. Elle n'est pas testée ici ; le récit la nomme à l'étape 3 et renvoie à `m.test()`. Une voie qui la relâche (matrice de covariance libre, par exemple brms) reste à décider.
- Les post-hocs de `m.test()` avec `(1|id)` ne sont pas appariés (PERIMETRE §4.4). Ici, ils le sont toujours.

## D16. ANCOVA : DÉCIDÉ le 9 novembre (points de la session parente intégrés)

**Plan.** Un ou plusieurs facteurs entre groupes, une ou plusieurs covariables continues, sans sujets. Les plans suivants sont refusés : covariable seule, ANCOVA appariée, MANCOVA.

**Calcul**, par `lmBF` [BayesFactor ; Rouder et Morey 2012 pour les covariables].
- **Covariables centrées** avant tout calcul, comme `m.test()` depuis 5609fb7 ([AikenWest1991]). L'effet du groupe se lit alors à la moyenne de la covariable, et non à covariable = 0, souvent hors des données. Le récit annonce le centre.
- Les covariables sont présentes dans tous les modèles, « aucune différence » compris. Chaque terme se lit contre le modèle « covariables seules ». C'est l'ajustement de l'ANCOVA, analogue au type II de `car::Anova` pour le facteur.
- L'interaction groupe × covariable (pentes différentes) est un terme de l'espace héréditaire : elle n'entre qu'avec le facteur. Elle est testée seulement si la formule l'écrit (`y ~ g * x`), comme dans `m.test()`.
- **A priori.**
  - Facteur : r = 0,5.
  - Covariable : `rscaleCont` = √2/4, le défaut « medium » de BayesFactor, sur la pente standardisée. BayesFactor applique cette même échelle, avec un g commun, à la covariable et à ses interactions avec un facteur ; le récit le dit.
  - Sensibilité : r et `rscaleCont` sont divisés et multipliés par 2 ensemble, croisés avec les chances de départ (un terme) ou les trois lois sur les modèles (plusieurs termes).
- **Forme du bruit** : loi normale seule, comme à plusieurs facteurs. La linéarité de la covariable n'est pas testée, et le récit le dit.

**Recouvrement des covariables (prudence).** L'étendue de la covariable est comparée entre groupes, sur ce qu'elles ont en commun : l'intervalle entre le plus grand des minimums et le plus petit des maximums par groupe.
- Si cette partie commune contient moins de la moitié des observations de l'un des groupes, ou si elle est vide, l'ajustement extrapole. Le récit le dit en « Prudence ajoutée », et le qualificatif de la réponse sur le facteur est plafonné à « faible », comme sous contradiction entre lois.

**Post-hoc.**
- **Pentes parallèles**, quand le terme groupe × covariable n'est pas en tête (P < 0,5) ou absent. Les regroupements suivent D11 avec les covariables : BF de « h + covariables » contre « covariables seules ». Les écarts ajustés, à la moyenne de la covariable, sont tirés de `posterior(lmBF(y ~ g + x))`.
- **Pentes différentes probables** (P ≥ 0,5). C'est un changement de stratégie (« --> »), puisqu'un écart ajusté unique n'a plus de sens.
  - Le récit donne les écarts entre groupes au premier quartile, à la médiane et au troisième quartile de la covariable, à partir de `posterior(lmBF(y ~ g * x))`.
  - Il donne aussi la zone de la covariable, dans l'étendue observée, où l'intervalle à 95 % de l'écart exclut zéro. C'est l'équivalent bayésien de Johnson-Neyman ([JohnsonNeyman1936] ; [Preacher2006]).
  - Pas de regroupements dans ce cas : ils supposeraient un écart unique.

**Contrôle indépendant.** `car::Anova(lm, type = 2)` pour le sens des réponses, `emmeans` pour les écarts ajustés et leurs intervalles, et `emmeans(at = quartiles)` quand les pentes diffèrent.

## D17. Réponses discrètes : PROPOSÉ le 9 novembre (mesures faites, code à venir)

Demandé par la session parente dans l'ordre du PERIMETRE, avec trois points de cadrage : la surdispersion, les tableaux de contingence et l'ordinal.

**Mesures faites avant de décider.**
- **BAS 2.0.2**, sans compilateur. `bas.glm(force.heredity = TRUE)` regroupe les indicatrices d'un facteur en un seul terme (2 modèles, pas 4) et donne la vraisemblance marginale d'un GLM de Poisson ou binomial sous l'a priori robuste [Bayarri2012].
- **Binaire à un facteur.** Sur 5 jeux à 3 groupes, `contingencyTableBF` (Gunel et Dickey 1974, multinomiales indépendantes, lignes fixées) et BAS robuste vont dans le sens du χ². BAS est plus favorable à H0, de 1,5 à 2 en log-BF : par exemple +0,74 contre −0,98 à p = 0,027. Les a priori diffèrent ; aucune des deux voies n'est fausse.
- **Binomiale négative.**
  - brms reste impraticable sans compilateur : plus de 6 minutes de compilation sans fin, mesuré par agent_06.
  - `rstanarm::stan_glm.nb` est précompilé. Avec `bridgesampling::bridge_sampler`, il faut 8,8 s pour deux modèles à n = 60.
  - Les a priori par défaut de rstanarm sont larges, faits pour estimer : le log-BF obtenu (−4,97) penche fortement vers H0, ce qui est l'effet de Lindley. L'a priori des coefficients doit donc être fixé pour le test.
- **Ordinal.**
  - Aucun paquet CRAN ne propose les tests de rangs bayésiens de van Doorn, Ly, Marsman et Wagenmakers (2020, J. Appl. Stat. 47(16), 2984-3006). Il n'en existe que des scripts OSF et le module JASP, et seulement pour 2 groupes.
  - `rstanarm::stan_polr` refuse le modèle nul (`o ~ 1`). Sous H0, cependant, les probabilités de catégorie suivent l'a priori de Dirichlet de `stan_polr`, et la vraisemblance marginale est alors une Dirichlet-multinomiale exacte.

**Proposition.**
1. **Binaire ou catégoriel, un facteur** : `contingencyTableBF`. C'est le pendant direct du χ² et du test G de `m.test()`. Post-hoc par paires de groupes sur le même modèle.
2. **Binaire avec plusieurs facteurs ou une covariable** : `BAS::bas.glm(binomial, robust)`, espace héréditaire. Contrôle croisé avec la voie 1 sur le cas à un facteur.
3. **Comptages** : comparaison de lois, comme à un facteur continu, sans choix préalable.
   - Poisson et binomiale négative (`rstanarm`, Suggests), même a priori sur le log du rapport de taux, BF par `bridgesampling` pour les quatre modèles : {Poisson, BN} × {H0, H1}.
   - La part de la BN est le constat de surdispersion. Si elle l'emporte, le récit le dit en « Changement de stratégie ».
   - Sans rstanarm, voie BAS Poisson seule, avec le rapport déviance/ddl comme alerte, et « Prudence ajoutée » au-delà de 1,5.
   - À valider sous H0 avec des comptages surdispersés.
4. **Ordinal (peu de valeurs ordonnées, non comptages)** : logit cumulatif `stan_polr` (Suggests). Le modèle complet passe par bridge sampling, le nul par la Dirichlet-multinomiale exacte, sous le même a priori sur les seuils. À défaut de rstanarm, `contingencyTableBF` sans l'ordre, annoncé comme tel.
5. **A priori des coefficients** des voies 3 et 4 : à décider après lecture (Gelman et al. 2008 ; [Bayarri2012] pour la cohérence avec BAS), puis à vérifier par la calibration sous H0.

**Ordre de travail.** 1, puis 2, puis 3 ; 4 en dernier.

**Voie 1 FAITE le 9 novembre** (`R/sys_bay_contingency.R`).
- **Calcul**
  - Réponse acceptée : facteur, texte, TRUE/FALSE ou 0/1, avec un seul facteur.
  - BF : `contingencyTableBF(sampleType = "indepMulti", fixedMargin = "rows", priorConcentration = 1)`.
  - Son a priori, retrouvé et testé à 1e-6 : Dirichlet(a) par groupe sous H1, Dirichlet(R(a−1)+1) commun sous H0. Il n'est donc défini que pour a ≥ 1.
  - Sensibilité : a = 1, 2, 4, croisé avec les chances de départ 1 contre 2 à 2 contre 1.
- **Post-hoc**
  - Regroupements de D11. Le BF d'un regroupement est calculé en forme close (rapport de Dirichlet-multinomiales) ; à a = 1, il égale `contingencyTableBF` du tableau fusionné.
  - Taille d'effet, à partir des lois a posteriori de Dirichlet : à 2 modalités, différence de proportions et log du rapport de cotes ; au-delà, écart total des proportions (de 0 à 1).
- **Récit et refus**
  - Le récit dit que le calcul est exact, sans effectif minimum, contrairement au khi-deux.
  - Une réponse binaire avec plusieurs facteurs, covariables ou sujets est refusée en attendant la voie 2.

## D17, voie 2 : comptages, FAITE le 10 novembre

**Choix automatique (révisé le 10 novembre, règle d'Antoine : la fonction choisit seule).** Avec `response = "auto"`, valeur par défaut, et un facteur, la voie des comptages est prise si toutes ces conditions sont réunies :
1. des entiers positifs ou nuls ;
2. au moins 3 valeurs distinctes (2 valeurs donnent la contingence) ;
3. une moyenne de groupe inférieure à 10, car l'approximation normale de la loi de Poisson est médiocre en dessous ;
4. des zéros, ou plus de 7 valeurs distinctes : cela écarte les échelles de notation de 1 à 5 ou de 1 à 7, sans zéro.

**Annonce et options.**
- La bascule est annoncée à l'étape 2, en « --> Changement de stratégie ».
- `response = "counts"` force la voie des comptages, `response = "continuous"` la voie continue.
- C'est une règle raisonnée, pas un critère publié unique. Elle est à surveiller sur des données réelles, en particulier des notes avec zéro, par exemple de 0 à 10.

**Autres choix.**
- **Nom de la réponse** : dans la forme x, g, la réponse prend le nom de l'objet passé, par exemple « nb_poissons » plutôt que « x ».
- **Avertissements de `glm.nb`** : ils sont comptés et non imprimés. S'il y en a, le récit ajoute « Prudence ajoutée : ajustement instable ».
- **Facteurs de Bayes** : une seule règle, à trois chiffres significatifs avec les milliers séparés (« 12 300 ») et « > 10^9 » au-delà. La sensibilité s'écrit « reste > 99,9 % quel que soit le point de départ ». Les cas hors de ce cadre sont refusés : plusieurs facteurs, covariables, sujets, valeurs négatives ou non entières.

**Calcul** (`R/sys_bay_counts.R`).
- **Lois comparées, non choisies.** Deux lois, Poisson (`glm`) et binomiale négative (`MASS::glm.nb`, déjà en Imports), chacune avec et sans différence : quatre modèles.
- **Vraisemblances marginales.** Elles sont approchées par −BIC/2, l'approximation d'information unitaire ([Schwarz1978], [Raftery1995], [Wagenmakers2007]) : sans compilateur, sans tirage, sans échelle d'a priori à choisir.
- **Réponse.** Elle est moyennée sur les deux lois, comme à un facteur continu.
- **Surdispersion.** C'est la part de la binomiale négative, annoncée en « Changement de stratégie » quand elle passe en tête, avec le rapport de Pearson à côté.
- **Écartés.**
  - BAS : il est de type Poisson ou binomial seulement, sans binomiale négative, et un passage de Poisson à la binomiale négative aurait changé d'a priori. Il sert ici de contre-vérification du sens du BF sous Poisson (test).
  - rstanarm avec bridgesampling : précompilé, mais ses a priori par défaut écrasent le BF vers H0 (log BF −4,97 mesuré), ce qui imposait de choisir un a priori de test.
  - brms : compilateur indisponible.
- **Post-hoc.** Regroupements de D11, chacun ajusté sous les deux lois. Rapports des moyennes, avec un intervalle à 95 % par l'approximation de Laplace, sous la loi en tête.
- **Sensibilité.** Seulement les chances de départ, de 1 contre 2 à 2 contre 1 : l'approximation BIC n'a pas d'échelle.

**Coût annoncé.** L'a priori d'information unitaire est conservateur. À effet réel, la puissance est inférieure au test du rapport de vraisemblance (VALIDATION §13).

**Restent** (D17, voies 3 et 4) : l'ordinal, `stan_polr` ou la contingence sans l'ordre, et la voie par rangs. Pour les rangs, van Doorn et al. (2020) n'existe dans aucun paquet CRAN, seulement en scripts OSF et dans JASP.

## D17, voie 3 : réponse ordonnée, et routage des entiers, FAITS le 10 novembre

**Défaut corrigé** (audit de 141d5ae). Des notes sur 20 partaient en comptages, alors qu'une note bornée n'est pas un comptage : Poisson et binomiale négative n'ont pas de borne haute.

**Routage d'une réponse entière à un facteur**, avec `response = "auto"` (`.bay_int_route`) :
1. **Échelle bornée probable.** Le maximum observé est un sommet d'échelle usuel : 4, 5, 6, 7, 9, 10, 20 ou 100.
   - Avec au plus 10 valeurs distinctes, la réponse va en voie ordinale.
   - Au-delà (notes sur 20, pourcentages entiers), elle reste continue. Le récit le dit en « --> Score borné… gardé continu » et nomme l'option `response = "ordinal"`. Une échelle à beaucoup de valeurs se traite couramment comme continue, comme le fait aussi `m.test()`.
2. **Sinon, comptages** : les critères de la voie 2.
3. **Sinon, voie continue.**

**Annonce et options.**
- Chaque bascule est dite en « --> Changement de stratégie » avec sa raison.
- `response = "counts" | "ordinal" | "continuous"` force la voie.

**Une règle raisonnée, pas un critère publié unique.**
- Un comptage dont le maximum tombe sur 10 ou 20 serait pris pour une échelle. C'est une erreur peu grave : modèle ordinal ou continu, annoncée, et forçable.
- Des notes de 0 à 10 sans maximum atteint seraient prises pour des comptages.

**Calcul de la voie 3** (`R/sys_bay_ordinal.R`).
- **Modèle.** Cotes proportionnelles, `MASS::polr` (McCullagh 1980), vraisemblances marginales par −BIC/2 comme la voie 2.
- **Post-hoc.** Regroupements de D11, rapports des cotes cumulées avec intervalle de Laplace.
- **Limite.** L'hypothèse des cotes proportionnelles n'est pas testée (test de Brant), et le récit le dit.
- **Écarté.** `rstanarm::stan_polr` : précompilé, mais il refuse le modèle nul, et ses a priori ne sont pas faits pour le test.

**Voie 4, rangs : non faite, décision motivée.**
- La méthode de van Doorn et al. (2020) n'existe dans aucun paquet CRAN : seulement en scripts OSF, de licence non précisée, et dans le module JASP. Elle se limite à 2 groupes ou à des paires.
- Les cas qu'elle viserait sont déjà couverts :
  - données ordonnées : voie 3 ;
  - mesures continues à valeurs éloignées : la comparaison des lois à queues lourdes de la voie à un facteur (D6) ;
  - petits effectifs : calcul exact, sans approximation asymptotique.
- Réécrire l'échantillonneur de Gibbs des rangs serait du code maison sans paquet de référence pour le contrôler. C'est contraire au principe d'Antoine. La voie reste ouverte si un paquet paraît.

## D18. Mise en forme commune des sorties : DÉCIDÉE le 10 novembre

**Origine.** Consigne d'Antoine : les sorties doivent être « dans un format agréable à lire étape par étape comme m.test() ». L'auditeur de lisibilité a rendu 12 recommandations (`relecture/lisibilite_2026-11-09/RAPPORT.md`). Antoine a ajouté une correction : 80 caractères au plus par ligne, car les consoles des développeurs sont souvent réduites. Aucune de ces décisions ne touche au calcul.

**Moteur** (`R/sys_bay_story.R`, `R/sys_bay_layout.R`).
- Chaque ligne tient dans min(getOption("width"), 80) colonnes, retrait compris, une tabulation comptant pour 8. Le texte se replie avec `strwrap` et un retrait suspendu, titres compris. Les espaces insécables protègent les intervalles, « % », « : », « ; », les guillemets, « = » et « < ».
- Les lignes d'une étape sont gardées jusqu'à l'étape suivante, ce qui permet de mettre le résultat dans le titre, comme m.test : « 3) Formes du bruit [m.bay() {KefiR}] - normale, variances égales en tête (51,9 %). »
- Chaque marque a un seul rôle :
  - « ==> » pour un constat ;
  - « --> » pour une décision ;
  - « Note : » pour le reste.
- Aucune ligne n'est sans marque. Les crochets n'apparaissent que sur la ligne-titre.
- Le jargon de méthode passe en lignes « detail », gardées dans `$story` et montrées par `print(x, detail = TRUE)` : Zellner-Siow, HOP, Beta(1,1), rscaleCont, échantillonnage d'importance, symétrie composée, Gibbs, références.
- Si une erreur ou un refus interrompt l'analyse, l'étape en cours s'imprime quand même (`on.exit`).

**Squelette identique pour tous les plans.**
1. Plan
2. Contrôles
3. Forme du bruit, ou loi des données
4. Point de départ
5. Solidité de la réponse : facteur de Bayes, précision, sensibilité
6. Paires, ou « non disponible » / « non demandé »
7. Bilan

**Bilan.**
- Il est séparé par une ligne vide et porte un titre au bord gauche : « 7) Bilan - <hypothèse en tête> : xx % (soutien). »
- Viennent ensuite :
  - le verdict ;
  - les hypothèses numérotées et alignées ;
  - le facteur de Bayes, dit en clair (« N fois plus probables sous l'hypothèse 1 que sous la 2 ») ;
  - la stabilité ;
  - les réponses du plan (termes, pentes, liste brève des paires) ;
  - les observations ;
  - une phrase prête pour un compte rendu, avec la méthode à citer.
- Les paires sont détaillées à l'étape 6 seulement : pas de doublon.
- Une paire n'est dite « semblable » que si le bilan n'est pas indécis ; sinon, elle est « non tranchée » (relecture D17).
- Les regroupements sont donnés « si une différence existe » (probabilités conditionnelles), pour qu'un seul chiffre porte « aucune différence » (relecture D17).

**Retouches après l'audit de c50a311 par la session parente.**
- Au plus deux « Note : » par étape, placées en fin d'étape ; les autres passent dans `detail = TRUE`.
- Plus de doubles deux-points dans les notes.
- Hypothèses listées sans flèche, en retrait sous « ==> Probabilité de chaque hypothèse : ».
- Bilan en trois blocs séparés par une ligne vide :
  - hypothèses ;
  - termes, paires et regroupement ;
  - observations et compte rendu.

**Nombres, une seule règle.**
- Probabilités : 1 décimale.
- Facteur de Bayes : `.bay_bf_txt` partout.
- Précision : « < 0,1 » en dessous de 0,1.
- Écarts : les décimales des moyennes ; δ et dz : 2 décimales ; proportions : en points.
- Pluriels accordés.

**Discret (relecture D17).**
- Les effectifs s'écrivent « 18 oui, 12 non (sur 30) », la modalité d'intérêt en premier.
- Le bilan dit quelle modalité porte l'écart le plus marqué.
- Le facteur de Bayes inférieur à 1 est dit aussi en faveur de l'absence (« soit 2,16 en faveur de l'absence de différence »).

**Contrôle.** `test-bay-story.R` vérifie qu'aucune ligne ne dépasse 80 colonnes, sur 7 plans (les 6 de l'audit et le discret), en français et en anglais, avec `width` à 80 et à 120. Il vérifie aussi 60 colonnes avec `options(width = 60)`. Il contrôle enfin la marque de chaque ligne, le bilan après une ligne vide et le rejeu à l'identique. Les sorties avant et après sont dans `relecture/lisibilite_2026-11-09/`.

**Hors périmètre.** `m.test()` n'est pas touché. Ses défauts de forme relevés par l'auditeur (point décimal, tableaux R en anglais, lettres contradictoires) relèvent d'une correction à part.

## D19. Variances inégales à plusieurs facteurs : CORRECTION DÉCIDÉE ET FAITE (10 novembre, option 1)

**Constat** (VALIDATION §15, `agent_64`). Sous H0, en 2 × 2 avec des variances inégales, le moteur à variance commune dit « oui » à P > 0,95 dans 12,5 % des jeux.

**Provisoire, fait.**
- Le moteur à un facteur, sans deuxième implémentation, compare les variances par case.
- Si les variances inégales sont probables, l'étape 3 dit « Prudence ajoutée ». Le soutien des réponses en faveur d'un effet est alors plafonné à faible et le verdict reste prudent.
- Résultats : « oui fort » affiché sous H0 passé de 12,5 % à 0 % ; 5 % de déclenchement à tort à variances égales ; aucune puissance perdue.

**Correction, à faire.** Moyenner {variance commune, variance par case} sur chaque modèle de l'espace héréditaire.
- **Ce qui ne marche pas directement.** La reprise « presque telle quelle » de `sys_bay_laws` ne suffit pas. À un facteur, chaque moyenne de groupe a sa propre échelle (a_i ~ N(0, g σ_i²), correction du défaut B2). À plusieurs facteurs, un terme comme A s'étend sur des cases d'échelles différentes.
- **Options.**
  1. Matrice de plan pondérée par case, avec un a priori de Zellner-Siow sur chaque terme rapporté à une échelle commune des cases (moyenne harmonique des σ_c²), et un échantillonnage d'importance sur les proportions de précision (Dirichlet, comme à un facteur). Simple, mais il rouvre la question B2 : une case dispersée sans rapport change l'échelle des termes.
  2. Modèles qui sont des partitions des cases (nul, A seul, B seul, plein), calculés exactement par le moteur à un facteur ; le modèle additif A + B, qui n'en est pas une, est approché.
  3. Stan : écarté (D18, compilation).
- **Décision à prendre** après un essai de l'option 1 sur `agent_64` et `agent_65`, avec le critère de D6. Le post-hoc à plusieurs facteurs (PERIMETRE §7, priorité 0) passe avant si l'essai dépasse une séance.

**Correction faite (10 novembre au soir) : option 1.** `R/sys_bay_fact_laws.R`, sans toucher `sys_bay_laws.R` dont il réutilise les briques : statistiques par case, constante de mesure, proposition défensive en Dirichlet, ajustement en t multivariée, estimateur d'importance.
- **Modèle.** y ~ N(μ + x_c'β, σ_c²), σ_c² = 1 / (w_c τ̄), w_c = K ϑ_c, ϑ ~ Dirichlet(1,5) sur les K cases ; β_t | g_t ~ N(0, g_t / τ̄ I), g_t ~ InvGamma(1/2, r²/2), un g par terme. Le codage x_c est celui de BayesFactor (contrastes orthonormés, produits pour une interaction). Avec w = 1, on retrouve exactement l'a priori de `lmBF` : contrôlé à 0,02 log près sur 2 × 2 et 3 × 2 (`test-bay-fact-laws.R`).
- **Échelle commune.** 1/τ̄ est la moyenne harmonique des variances des cases. Un effet additif s'étend sur des cases de dispersions différentes ; il lui faut une seule échelle. Le prix est celui de B2 à un facteur : une case très dispersée déplace un peu cette échelle. Il est accepté ici, car les critères ci-dessous sont tenus.
- **Calcul.** Pour chaque modèle héréditaire, μ, τ̄ et β s'intègrent en forme close. L'intégrale sur (ϑ, g) se fait par échantillonnage d'importance autour du mode : t multivariée à 4 ddl à 85 %, partie défensive à 15 %, au moins 4000 tirages, et un réajustement si l'ESS < tirages/5. Les mêmes tirages servent aux trois r de la sensibilité. La loi à variance commune garde `lmBF`. Les deux lois sont moyennées à poids a priori égaux dans le facteur de Bayes de chaque modèle contre le nul. L'espace des modèles et son a priori HOP ne changent pas.
- **Conditions.** Au moins 3 valeurs par case, pas de sujets ni de covariables, et `normal_het` admis par `.bay_laws_allowed` (cases constantes exclues). Sinon, la loi commune est gardée et le récit le dit. Si l'échantillonnage échoue (ESS < 30), la prudence provisoire reprend, avec un plafond quand le rapport des écarts-types est ≥ 2.
- **Résultats** (`agent_66`, VALIDATION §16). Sous H0, en 2 × 2 à variances inégales, P > 0,95 sur un terme passe de 12,5 % à 1,5 % sur 200 jeux, sans aucun plafond (total sur 800 jeux : voir VALIDATION §16). À variances égales, P(A) > 0,95 reste à 17,5 % : aucune puissance perdue. Temps médian 1,8 s, maximum 4,1 s, à N = 40.
- **Non fait.** Lois t à plusieurs facteurs. Variances par case avec sujets ou covariables (plafond provisoire non appliqué non plus : seul le récit le signale).


## D20. Post-hoc à plusieurs facteurs : DÉCIDÉ ET FAIT (10 novembre)

**Demande** (PERIMETRE §7, priorité 0). En plan factoriel, ANCOVA à plusieurs facteurs ou split-plot : comparer les effets principaux quand aucune interaction n'est probable, et les cases quand A:B l'est.

**Règle** (`.bay_posthoc_factorial`, `R/sys_bay_posthoc.R`).
- On parcourt les termes de facteurs du plus haut ordre au plus bas. Un terme est détaillé si sa probabilité vaut au moins 1/2 et s'il n'est pas contenu dans un terme déjà détaillé.
- Une interaction détaillée compare ses cases ; un effet principal détaillé compare ses niveaux.
- Les autres facteurs (effets principaux), les covariables et les sujets entrent dans chaque regroupement comparé comme termes de nuisance. Chaque écart se lit donc « à autres facteurs égaux ».
- Chaque famille passe par `.bay_posthoc_oneway` : regroupements des niveaux, a priori d'Ewens, mêmes règles qu'à un facteur. Calcul exact jusqu'à 5 groupes avec nuisance (option `KefiR.bay_exact_max_rm`), Gibbs au-delà.
- Les termes non probables ne sont pas détaillés ; le récit les nomme. Les termes avec une covariable (pentes) ne sont pas comparés ici.

**Écart standardisé.** Il est rapporté à l'écart-type résiduel, autres termes pris en compte. Avec des sujets, c'est le résiduel au sein des sujets. Pour un facteur inter-sujets d'un split-plot, ce résiduel ne mesure pas l'écart entre sujets : l'écart standardisé n'est pas donné, seul l'écart brut l'est.

**Choix et limites.**
- Le seuil de 1/2 suit la lecture du bilan (le terme est « plutôt oui »). Il ne s'appuie sur aucune source : c'est un choix KefiR.
- Les a priori des familles sont indépendants de l'a priori HOP sur les modèles (prior_h1 = 1/2 par famille, comme à un facteur).
- Pas d'effets simples (A à chaque niveau de B) : les cases donnent la même information en un seul regroupement.
- Avec des variances inégales probables (D19), les paires supposent une variance commune et le récit le dit.

**Audit de 23c941e par « KefiR officiel » (10 novembre), quatre retouches faites.**
1. La phrase du compte rendu reprend les mots du bilan pour chaque terme (`.bay_term_answers`, une seule règle) : un terme « non tranché » au bilan l'est aussi dans le compte rendu.
2. Une interaction non tranchée (ou dont la probabilité est entre 25 et 50 %) sur un facteur comparé comme effet principal donne une note, gardée parmi les deux notes de l'étape. Les écarts de ce facteur à chaque niveau de l'autre sont tirés du modèle des cases (`lmBF`, `posterior`, `hdi`) et sont dans `detail = TRUE` et `$posthoc$by`.
3. Le repli d'un élément numéroté (« 1. … ») se fait sous son texte, après le numéro.
4. Le titre de l'étape 6 donne la réponse : « A : a2 > a1 (2,72) », ou le regroupement en tête au-delà de deux niveaux.

## D21. Mesures répétées avec répliques (modèle mixte) : DÉCIDÉ ET FAIT (10 novembre)

**Demande** (PERIMETRE §7, priorité suivante). Plusieurs mesures d'un même sujet dans une même condition. C'est courant en TP, et la réponse était jusqu'ici un refus.

**Choix.**
- Les effets aléatoires sont le sujet, puis le sujet croisé avec chaque combinaison de facteurs intra-sujet présente dans la formule (sujet × temps ; avec deux facteurs intra, sujet × A, sujet × B et sujet × A × B). Ils figurent dans chaque modèle, « aucune différence » compris. Chaque facteur de Bayes se lit contre le modèle des seuls effets aléatoires.
- Le calcul passe par `lmBF` avec `whichRandom` (Rouder et al. 2012, effets aléatoires d'échelle 1). Les termes sujet × condition y entrent comme facteurs aléatoires à part entière, ce qui équivaut au terme `cond:id` de BayesFactor. Contrôle : `test-bay-mixed.R`, contre `lmBF(y ~ cond + id + cond:id)` à 0,2 log près.
- Le post-hoc (un facteur, ou D20 à plusieurs) porte les mêmes effets aléatoires en nuisance.
- **Refusés.** Les pentes aléatoires (il faudrait Stan, qui ne compile pas ici) et les covariables avec des sujets.

**Pourquoi le terme sujet × condition.** Sans lui, les répliques comptent comme autant de sujets indépendants (pseudo-réplication). Dans `agent_67`, le modèle à sujet seul dit « oui » à P > 0,95 dans 7 % des jeux sous H0, et m.bay dans 0 % (VALIDATION §18).

**Limite.** La puissance est prudente : P > 0,95 dans 4 % des jeux pour un effet où l'ANOVA à mesures répétées sur les moyennes donne p < 0,05 dans 26,5 %. C'est le même écart que sur les autres voies entre P > 0,95 et p < 0,05.

**Audit de 8ffb16a par « KefiR officiel » (10 novembre), deux retouches faites.**
1. Les moyennes par niveau disent « 10 sujets x 3 répliques » (ou « x 3 mesures » pour un facteur entre sujets), au lieu de « n = 30 » qui faisait croire à 30 sujets.
2. **Unité des écarts standardisés avec répliques : décidé.** C'est l'écart-type d'une mesure d'un sujet dans une condition, soit le résiduel (entre répliques) plus le terme sujet × condition : σ²(1 + Σ g_sujet×condition), tiré du postérieur de `lmBF`. C'est la variabilité qu'un expérimentateur attend d'une mesure une fois le sujet connu. Le résiduel seul la sous-estimerait et grossirait les écarts standardisés. Le bilan donne cet écart-type sous ce nom, et l'étape 6 y renvoie.

## D22. MANOVA à un facteur : DÉCIDÉ ET FAIT (11 novembre), remplace la proposition D7

**Demande** (PERIMETRE §7). Plusieurs réponses comparées ensemble entre groupes. D7 avait arrêté un moteur maison Normal-Inverse-Wishart : l'échelle de son a priori croissait avec le nombre de réponses.

**Trois moteurs essayés** (`agent_68` ; 3 groupes de 15, 100 jeux par scénario, VALIDATION §19).
- **A.** `BFpack::BF` sur `aov` multivarié : facteur de Bayes fractionnaire (Mulder et al. 2021, JSS 100(18)). Paquet établi, installé dans mbay-lib pour l'essai seulement.
- **B.** BIC multivarié : loi normale à covariance libre, information unitaire (Kass et Raftery 1995).
- **C.** `lmBF` sur les réponses standardisées et empilées, y ~ g * réponse + individu aléatoire (Rouder et al. 2012). « Au moins une différence » réunit g seul (même écart partout) et g + g:réponse (écarts qui changent d'une réponse à l'autre), à poids égaux.

**Choix : C.**
- C est le seul à garder de la puissance : effet réparti sur 3 réponses, P > 0,75 dans 21 % des jeux contre 3 % (A) et 2 % (B), Pillai p < 0,05 dans 27 %. À 10 réponses, A et B ne détectent plus rien.
- Sous H0, C ne dépasse jamais 0,95 (0 % à d = 3 et d = 10) et reste à 1 à 3 % au-delà de 0,75.
- Il repose sur BayesFactor, déjà employé partout, sans nouvelle dépendance.

**Limite de C, mesurée.** L'effet individu suppose la même corrélation pour chaque paire de réponses (symétrie composée).
- Avec des corrélations très inégales (0,8 / 0 / 0), C reste calibré sous H0 (1 % à P > 0,75).
- Il perd en revanche de la puissance quand l'effet porte sur une réponse très corrélée à une autre : 15 % à P > 0,75, contre 84 % pour Pillai, qui exploite le contraste entre les deux réponses, et 29 % pour A.
- Le risque est donc de manquer une différence, pas d'en inventer. Quand les corrélations s'écartent de plus de 0,5, le récit et le bilan disent qu'un « non » est à lire avec réserve. Pas de plafond sur un « oui ».

**Lecture par réponse** (étape 6) : chaque réponse est lue seule, avec le moteur à un facteur à chances égales, pour voir où se trouve la différence. La réponse globale reste celle de l'étape 5.

**Refusés.** MANOVA à plusieurs facteurs, avec covariables ou avec sujets.

**Piste.** Moyenner C avec un moteur qui exploite la covariance (A ou B) quand les corrélations sont inégales, après calibration.

**Audit de D22 par « KefiR officiel » (11 novembre), deux corrections faites avant le commit.**
1. **Largeur.** Le repli d'un élément numéroté pouvait atteindre 83 colonnes (iris, noms de réponses longs) : la largeur disponible ne tenait pas compte du retrait de continuation. Corrigé à la source dans `.bay_wrap`, pour tous les plans : la largeur se calcule sur le plus large des deux retraits. iris est ajouté au test des 80 colonnes.
2. **Noms du résultat.** La MANOVA a désormais `laws`, `models` et `questions`, comme les autres plans. La probabilité d'un effet s'appelle `prob_effect` partout, MANOVA comprise ; aucun plan n'a d'élément `p_h1` à la racine (il est dans `$effects$p_h1`).

## D23. MANOVA à plusieurs facteurs entre groupes : DÉCIDÉ ET FAIT (11 novembre)

**Faisabilité.** Avec BayesFactor empilé, oui, proprement. Un plan 2 × 2 à 3 réponses compte 17 modèles et environ 20 s par appel sans sensibilité. Avec la sensibilité (3 valeurs de r), il faut 50 à 70 s : une ligne d'attente s'affiche alors (« Calcul de 17 modèles pour 3 valeurs de r, environ 70 s… »).

**Modèle** (`R/sys_bay_manova_fact.R`). Même empilement qu'à un facteur (D22) : réponses standardisées, individu aléatoire, `lmBF`.
- Chaque terme des facteurs entre sous l'une de deux formes : t (même écart sur chaque réponse) ou t + t:réponse (écarts qui changent avec la réponse).
- Les modèles sont les combinaisons héréditaires des termes (comme à une réponse), chaque terme présent sous l'une des deux formes.
- A priori : HOP sur les termes, puis 1/2 par forme. Chaque facteur de Bayes se lit contre « réponse + individu ».
- **Choix KefiR.** Les formes de deux termes sont indépendantes : A:B:réponse n'exige pas A:réponse. Cela garde 17 modèles au lieu de plusieurs dizaines.

**Réponses par terme.** P(terme) est la somme des modèles qui le contiennent, sous l'une ou l'autre forme. Le bilan et le compte rendu emploient les mêmes mots (`.bay_term_answers`). L'étape 6 donne, pour chaque terme, P(écarts qui changent avec la réponse | le terme agit).

**Calibration** (`agent_69`, 2 × 2, 3 réponses, n = 12 par case, 60 jeux par scénario ; VALIDATION §20).
- Sous H0, aucun terme ne dépasse P > 0,95 ; un terme dépasse 0,75 dans 1,7 % des jeux.
- Avec un effet de A de 0,8 sur une réponse, P(A) > 0,75 dans 38 % des jeux et > 0,95 dans 17 %, contre 77 % pour Pillai à p < 0,05. B ou A:B ne dépassent 0,95 que dans 1 jeu sur 60.

**Refusés.** MANOVA avec covariables ou avec sujets.

**Audit de 9c44239 par « KefiR officiel » (11 novembre) : calcul validé, trois retouches.**
1. **Bilan mené par les termes, pour tous les plans à plusieurs termes** (`.bay_bilan_common`). Dès qu'un terme reçoit un « oui », trois choses changent :
   - le titre donne les réponses par terme (« 7) Bilan - A oui (96,3 %), B non (19,4 %), A:B non (5,8 %) ») ;
   - le verdict dit « une différence selon A ; les autres termes sont plus bas » ;
   - la phrase de compte rendu commence par les termes, puis donne « au moins une différence : … ».
   « Indécis » n'est alors jamais dit. Sans aucun « oui », le libellé global reste.
2. **Lecture par réponse faite par m.bay.** À l'étape 6, chaque réponse est lue seule par le moteur factoriel (mêmes modèles, même a priori HOP). Au-delà de 10 modèles, seuls les termes au-dessus de 50 % sont lus, et le récit le dit. La lecture est dans `$responses` et dans le bilan.
3. **Puissance : l'écart avec Pillai ne tient pas à l'effet porté par une seule réponse.** Mesure : `agent_69`, scénario Arep, même plan, effet de A de 0,45 réparti sur les 3 réponses.
   - Effet sur une réponse : m.bay (P(A) > 0,75) / Pillai (p < 0,05) = 0,38 / 0,77 = 0,50.
   - Effet réparti : 0,17 / 0,32 = 0,53.
   - Le rapport est le même : c'est la prudence générale de P > 0,75 face à p < 0,05, comme sur les autres voies (D17, D21). La condition de l'avertissement n'est donc pas élargie.
   - Reste l'avertissement de D22, mesuré : avec des corrélations très inégales, le modèle empilé peut manquer une différence portée par un contraste entre réponses corrélées.

## D24. Queues lourdes à plusieurs facteurs entre groupes : DÉCIDÉ ET FAIT (11 novembre)

**Défaut mesuré d'abord** (`agent_70`, 2 × 2, n = 12 par case, bruit de Student à 3 ddl). Avec un effet de A de 0,8, le moteur à plusieurs facteurs (normales, D19 compris) donne P(A) > 0,95 dans 44 % des jeux, contre 69 % pour la voie à un facteur avec lois t appliquée à A seul. Sous H0 : 0 %.

**Correction** (`R/sys_bay_fact_t.R`). Une troisième loi, Student à échelle commune, est calculée pour chaque modèle héréditaire, puis moyennée à poids égaux avec les deux lois normales de D19, dans le facteur de Bayes de chaque modèle.
- **Mesure commune.** Sur les effets, β_t/σ suit une loi de Cauchy multivariée d'échelle r (la Zellner-Siow de `lmBF`, g intégré). p(ν) = exp(-(ν - 2)), comme à un facteur. Les données sont standardisées comme dans `.bay_laws_stats`.
- **Calcul.** Échantillonnage d'importance autour du mode (t multivariée, 4 ddl), 4000 tirages au moins ; les mêmes tirages servent aux trois r.
- **Contrôles** (`test-bay-fact-t.R`).
  - Le même code, avec la vraisemblance normale, redonne `lmBF` à 0,1 log près.
  - Il redonne aussi la forme close du modèle nul de `sys_bay_laws`, une fois ajouté log 2 : `sys_bay_laws` pose 1/τ sur la précision, ce code 1/σ, et les deux mesures diffèrent d'un facteur 2.
- **Conditions.** Celles de D19 (entre groupes, au moins 3 valeurs par case), plus `t_eq` admis par `.bay_laws_allowed` (exclu quand les ex-aequo sont trop nombreux). Un ESS sous 30 fait retirer la loi t, et le récit le dit.

**Résultat** (`agent_70`, rejoué sur les 100 premiers jeux des mêmes graines ; VALIDATION §21).
- Avec l'effet de A : P(A) > 0,95 passe de 44 % à 64 %, et P(A) > 0,75 de 70 % à 86 %. La voie à un facteur fait 71 % ; l'ANOVA sur rangs 96 % à p < 0,05.
- Sous H0 : 0 % au-delà de 0,75 et de 0,95.
- Une case à valeurs extrêmes peut aussi être lue comme « variance par case » : les deux lois se partagent alors les queues lourdes, sans changer la réponse.

**Non fait.** Student avec une échelle par case ; queues lourdes avec sujets ou covariables, et en MANOVA.

**Audit de aed2319 et a206ba5 par « KefiR officiel » (11 novembre), quatre retouches faites.**
1. **Écarts de l'étape 6 avec les queues lourdes : décidé.** Ils sont tirés du même mélange de lois que la réponse.
   - Les tirages du postérieur normal (`lmBF`) sont complétés par ceux du postérieur de Student (rééchantillonnage d'importance, `.bay_ft_gaps`), en proportion de P(Student) à l'étape 3.
   - La loi à variance par case est représentée par la normale commune ; c'est une approximation, que le récit ne détaille pas.
   - Quand Student l'emporte (P ≥ 0,5), l'écart standardisé est rapporté à l'échelle de la loi de Student, c'est-à-dire à la dispersion sans les queues. L'écart-type, lui, est gonflé par les queues.
   - Sur le jeu de l'audit, l'intervalle passe d'un [-0,13 ; 1,71] qui couvrait 0 à [0,19 ; 1,17].
   - **Garde-fou.** Si un « oui » à plus de 95 % s'affiche encore avec un intervalle qui couvre 0, une ligne explique que la probabilité compare des modèles entiers alors que l'intervalle estime un écart ramené vers 0. Testé (`test-bay-fact-t.R`).
2. **Étape 3, « 95,2 % » en double.** Ce n'était pas la même variable : « égales probables (95,2 %) » affichait 1 − P(variance par case), soit P(commune) + P(Student). La ligne donne désormais P(variance par case) elle-même (« inégales à 4,8 % »).
3. **Hypothèse globale sous un « oui ».** Elle garde le soutien de son propre facteur de Bayes sur l'échelle du bilan (« soutien faible »), jamais « non tranché ».
4. **Écart-type du bilan** à trois chiffres significatifs (`.bay_sig3`).
