# Relecture adverse de la méthodologie de m.bay

Date : 2026-11-05. Objet : `dev/mbay/METHODOLOGIE.md`, scripts `agent_02`, `agent_03`, `agent_05` et leurs sorties, `BIBLIOGRAPHIE.md`.

Scripts et sorties de cette relecture, tous dans `C:\Users\masse\AppData\Local\Temp\claude\dd1494da25364fcc\7e07ce24-9e29-4154-a698-e9d4fc5f31b8\scratchpad\relecture\` : `r1_algebre.R`, `r2_termes.R`, `r3_rope.R` (puis `r3b.R`, `r3c.R`), `r4_manova_part.R` (puis `r4b.R`), `r5_divers.R`, `r6.R` ; sorties `r1.out`, `r2.out`, `r3.out`, `r3b.out`, `r3c.out`, `r4.out`, `r4b.out`, `r5.out`, `r6.out`, `niw_rerun.out`.

Aucun fichier du dépôt n'a été modifié.

## Verdict

L'algèbre MANOVA résiste. Ce qui ne résiste pas, c'est ce que le chiffre est censé dire. Avec les réglages par défaut :

- « Oui avec probabilité p ; non avec probabilité 1-p » n'est pas une probabilité que l'utilisateur peut lire ainsi pour un effet principal dès qu'il y a d'autres facteurs ;
- la ligne ROPE ne porte presque aucune information propre dès 3 groupes ;
- la branche gaussienne annonce P > 0,95 sur du bruit pur dans 25 % des cas d'un plan hétéroscédastique banal, et le recours prévu dépend de Stan, non ajustable sans compilateur ;
- le « π = 1/2 par terme » retenu et l'« uniforme sur les modèles » rejeté sont, en fait, la même loi.

| Gravité | Nombre |
|---|---|
| Bloquant | 6 |
| Important | 13 |
| Mineur | 7 |

## Méthode de contrôle et ses limites

- Quadrature et tirages maison validés : forme fermée contre Monte Carlo indépendant (0,80332 contre 0,80385, se 0,008) ; d = 1 contre `anovaBF` déséquilibré, identique à 10⁻⁵ ; 2×2 contre `lmBF` (167,0/167,0 ; 53,7/53,8 ; 21,97/21,89) ; posterior de f contre `posterior()` (0,0257/0,0220 ; médiane 0,2045/0,2080).
- **Non recoupé** : le calcul 2³ (D, B5) étend ces formules à 7 termes par Monte Carlo a priori, sans comparaison directe à `anovaBF`.
- Taille des simulations : 100 à 300 jeux selon les cas, 1 500 pour la calibration des termes. Les écarts de quelques points ne sont pas significatifs. Un seul plan par scénario : ce sont des contre-exemples, pas une table de performance.
- Bibliographie : résumés et notices lus via un résumeur automatique, jamais le texte intégral.

---

# A. Défauts vérifiés par le calcul ou par lecture directe

## Bloquants

### B1. « Oui p / Non 1-p » est faux pour un effet principal dès qu'une interaction joue
**Preuve** (`r2.out` §A). Plan 2×2, n = 10 par cellule, croisement pur, effet de A de signe opposé selon B (p = 1,8·10⁻⁵ et 2,2·10⁻⁶ pour les effets simples).

| Quantité | Valeur |
|---|---|
| P(A) apparié, le défaut | **0,236**, soit « Non avec probabilité 0,764 » |
| P(A:B) apparié | 1,00000 |
| P(A), uniforme sur 5 modèles | 1,00000 |
| P(A), Chipman | 1,00000 |

Sous la marginalité, A:B implique A, donc P(A:B) ≤ P(A) pour toute loi jointe. Le couple (0,236 ; 1,000) n'est donc la marginale d'**aucune** loi. §3 parle de « prix à payer pour une réponse lisible », et §5.3 prévoit un signalement, mais la ligne imprimée reste fausse telle qu'elle se lit : A agit manifestement.

Sur 1 500 jeux générés selon l'uniforme, parmi ceux dont le modèle vrai contient A:B (donc A), P(A) apparié est < 0,5 dans **44,6 %** des cas.

**Correction.** Écrire la condition dans la phrase : « S'il n'y a pas d'interaction impliquant A : oui avec p ». Ne jamais imprimer « Non » pour un terme apparié, et supprimer la ligne quand une interaction parente n'est pas négligeable. Mieux : prendre comme défaut une inclusion cohérente (Chipman, déjà listée) et garder l'apparié comme lecture Type II nommée.

### B2. « π = 1/2 par terme » et « uniforme sur les modèles » sont la même loi, et la calibration promise est fausse au défaut
**Preuve algébrique.** Pour 2 facteurs, avec la loi uniforme sur 5 modèles :

- P(A | D, pas de A:B) = [m(A) + m(A+B)] / [m(1) + m(A) + m(B) + m(A+B)], exactement le BF apparié avec cotes 1:1 (`r2.out` §B : 0,235841 des deux côtés) ;
- P(A:B | D, A et B présents) = m(plein) / [m(plein) + m(A+B)], idem.

Le défaut de §3 **est** donc le conditionnel de la loi que le tableau de §3 marque « Rejeté comme défaut », avec l'argument « Toute annonce P(H0) = 1/2 serait fausse avec ce réglage ». Le 1/2 annoncé est une cote *conditionnelle au sous-espace*, pas P(A). La loi où P(A) vaut vraiment 1/2 est Chipman, et elle pondère autrement (1/4, 1/4, 1/4, 1/8, 1/8) : pas le même nombre.

**Preuve numérique** (`r2.out` §C, 1 500 jeux, tranches de p annoncé).

| Monde générateur | Lecture | Annoncé → réalisé | Moyenne annoncée / vraie |
|---|---|---|---|
| Uniforme | Apparié, tous jeux | 0,255 → 0,356 ; 0,371 → 0,459 ; 0,591 → 0,664 ; 0,800 → 0,884 | 0,526 / 0,597 |
| Uniforme | Apparié, modèles sans interaction | 0,256 → 0,248 ; 0,592 → 0,586 | 0,499 / 0,497 |
| Uniforme | Témoin : P(A) cohérent | 0,275 → 0,242 ; 0,589 → 0,581 | 0,603 / 0,597 |
| Chipman | Apparié | 0,591 → 0,672 ; 0,792 → 0,756 | 0,490 / 0,508 |

§10 écrit « parmi les analyses annonçant p, une proportion p ont H1 vraie ». Le défaut ne le tient que **conditionnellement à l'absence d'interaction**, ce que l'utilisateur ne sait pas. Le contrôle DAP de §10 passera ou non selon la condition choisie, et il faut la fixer avant de l'écrire.

**Correction.** Retirer le rejet de l'uniforme ou le défaut apparié, l'un des deux tombe. Annoncer sur la ligne « cote 1:1 *sachant* l'absence d'interaction d'ordre supérieur ».

### B3. La ligne ROPE est presque une copie de la ligne 1 dès a ≥ 3, à cause du spike et du biais positif de f
**Preuve** (`r3.out`, 40 jeux par cellule, médianes).

| a | n/groupe | f vrai | P(H1\|D) | P(f<0,05\|D,H1) | P(ROPE\|D) | Médiane f \| H1 |
|---|---:|---:|---:|---:|---:|---:|
| 5 | 100 | 0 | 0,016 | 0,044 | 0,984 | **0,105** |
| 10 | 100 | 0 | 0,001 | 0,004 | 0,999 | 0,108 |
| 10 | 30 | 0 | 0,022 | 0,000 | 0,978 | 0,170 |
| 5 | 30 | **0,20** | 0,449 | 0,002 | **0,552** | 0,233 |
| 3 | 30 | **0,20** | 0,460 | 0,020 | **0,549** | 0,224 |
| 2 | 30 | 0,20 | 0,412 | 0,130 | 0,646 | 0,165 |

Trois mécanismes :

1. **Biais de f** : f ≥ 0, et pour une vérité nulle f | H1 vaut environ √(q/(a·n)) (a = 5 : 0,237 ; 0,159 ; 0,098 ; 0,050 à n = 10 ; 30 ; 100 ; 500). Il faut 500 par groupe pour que la médiane atteigne le seuil 0,05, sur des données **sans aucun effet**. Le terme P(|θ|<w | D, H1) de la formule de §4 est donc ≈ 0.
2. **Spike** : P(ROPE) ≈ P(H0), la ligne 2 redit la ligne 1 sous un autre nom.
3. **Confusion d'évidence** : avec f = 0,20, soit **4 fois le seuil**, à n = 30, le bilan dit « négligeable 0,55 ». Ce 0,55 vient de l'a priori 0,50, pas d'une mesure de l'effet.

**Masse a priori** (`r1.out` §5) : P(f<0,05 | H1) = 0,089 ; 0,015 ; 0,003 ; 0,0009 ; 0,0003 pour a = 2 ; 3 ; 4 ; 5 ; 6. Le « 0,55 » de §4 ne vaut que pour a = 2, or l'exemple de bilan l'imprime avec r = 0,5 pour un terme quelconque. À 3 groupes, c'est 0,507. Sur données non nulles à effectifs modérés, un P(ROPE) d'environ 0,5 se lit donc « peu de chose appris », pas « une chance sur deux que ce soit négligeable ».

**Correction.**
- Séparer spike et continu : imprimer P(H0|D) et, à part, P(|effet|<w | D, H1) avec sa valeur a priori.
- Pour a > 2, passer aux contrastes par paire ou à un f sans biais (estimateur de type ω², ou f calculé sur le paramètre de non-centralité).
- Mieux, tester l'intervalle comme hypothèse, comme le propose [MoreyRouder2011], déjà en bibliographie.
- Dans le bilan, n'écrire « négligeable » que si P(ROPE|D) dépasse nettement sa valeur a priori.

### B4. Le niveau 1 gaussien annonce n'importe quoi en hétéroscédasticité, et le recours prévu n'est pas livrable
**Preuve** (`r3c.out` §4, 200 jeux).

| Scénario | P>0,5 | P>0,95 | P<0,25 | P médiane | p aov<0,05 | KW<0,05 |
|---|---:|---:|---:|---:|---:|---:|
| H0, n = (5, 5, 40), sd = (4, 4, 1) | **0,575** | **0,255** | 0,13 | 0,611 | 0,52 | 0,235 |
| H0, n = (5, 5, 40), sd = (1, 1, 4) | 0,000 | 0,000 | 0,905 | 0,207 | 0,00 | 0,00 |
| Log-normal, effet, données brutes | 0,27 | 0,015 | 0,265 | 0,351 | 0,175 | 0,49 |
| Même effet, données log | 0,51 | 0,075 | 0,17 | 0,513 | 0,41 | 0,41 |

Une seule valeur aberrante de +30 sur 45, moyennes −1, 0, 1 : P passe de **0,994 à 0,168** (« Non 0,83 »), alors que Kruskal-Wallis donne p = 0,003.

§1.4 et §7 font du gaussien la ligne toujours présente, et le niveau 2 s'ajoute « sur alerte ». Mais le niveau 2 exige « un code Stan dédié » (§7), et `resultats/rapport_moteurs.md` établit que brms échoue sans compilateur (l.37, l.203) et que l'essai avec Rtools 4.4 n'a produit aucun résultat en plus de 6 min (l.211). Les étudiants n'auront donc que la ligne fausse. Le point 5 de §11 laisse encore ouvert si le niveau 2 s'active par défaut. Aggravant : `ttestBF(var.equal=FALSE)` est accepté en silence et donne un BF à variance commune (l.114 du rapport).

**Correction.**
- Sur alerte, ne pas donner P gaussien comme réponse principale : l'imprimer « non interprétable » si le niveau 2 n'est pas disponible.
- Écrire une voie robuste **sans Stan** : par exemple BF sur scores normaux (normal scores) à valider par simulation, ou échantillonneur de Gibbs à mélange d'échelles en R pour t et σ_j.

§7 dit « le bilan ne désigne pas de gagnant », mais à ce stade il n'y a qu'un candidat.

### B5. π fixe par terme en plan factoriel : « Oui » sur du bruit dans plus d'une étude sur deux
**Preuve** (`r2.out` §D, 300 jeux, 2³ sous H0 globale, 7 termes, apparié).

| N | P>0,5 par terme | Au moins 1 « Oui » sur 7 | P>0,75, au moins 1 | P>0,95, au moins 1 | p<0,05, au moins 1 |
|---:|---:|---:|---:|---:|---:|
| 16 | 0,174 | **0,740** | 0,227 | 0,020 | 0,257 |
| 40 | 0,114 | **0,570** | 0,130 | 0,017 | 0,283 |

En moyenne 0,80 à 1,22 terme « probablement réel » sur du bruit. Le seuil de 7 termes de §3 (« la lisibilité l'emporte ») tombe exactement sur le plan d'expériences à 2 niveaux de l'industriel. Les termes d'interaction ont un X'X plus petit (N/4 et N/8 contre N/2) et atteignent 1/2 plus facilement (mécanisme déduit de la formule, pas isolé par terme ici).

**Correction.** Afficher au niveau de l'étude P(au moins un terme) et l'espérance du nombre de termes actifs a posteriori. Adopter un a priori hiérarchique sur l'inclusion **dès 2 facteurs**, après définition d'une version compatible avec la marginalité (I12). Donner la table de faux positifs par plan dans le bilan, pas seulement dans une vignette.

### B6. MANOVA : l'a priori gonfle la distance de Mahalanobis en √d et transforme un vrai effet en « Non 98 % »
**Preuve algébrique.** tr(Σ⁻¹B'B) | g ~ g·χ²_{qd}. « Chaque combinaison linéaire reçoit l'a priori univarié » (§6) est exact, et c'est le problème : la combinaison *maximale* n'a pas cet a priori.

**Preuve numérique** (`r4b.out` §2, 300 jeux). Effet sur **une** DV, moyennes −0,5 ; 0 ; 0,5 (f = 0,41), les autres DV étant du bruit indépendant.

| d | n/gr | Méd. Mahalanobis a priori | P médiane | P>0,5 | P<0,25 (« Non ») | Pillai <0,05 | H0 : P<0,05 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 30 | 0,50 | 0,977 | 0,927 | 0,013 | 0,913 | 0,000 |
| 3 | 30 | 0,99 | 0,849 | 0,737 | 0,130 | 0,803 | 0,590 |
| 5 | 30 | 1,30 | 0,619 | 0,553 | 0,303 | 0,693 | 0,707 |
| 10 | 30 | 1,88 | 0,133 | 0,240 | 0,627 | 0,510 | 0,897 |
| 20 | 30 | 2,69 | **0,016** | 0,083 | **0,827** | 0,370 | 0,960 |
| 10 | 15 | 1,88 | 0,085 | 0,117 | 0,780 | 0,173 | 0,637 |

Il y a un vrai effet, et le bilan dit « Non » avec 0,98. Pillai perd aussi de la puissance, mais il ne conclut pas à l'absence.

La calibration sous l'a priori exact, elle, tient (d = 3 : moyenne 0,493 ; d = 5 : 0,509) : la formule fait bien ce qu'elle dit. C'est l'a priori qui n'a pas de sens pour d ≥ 5. §6 attribue la baisse à une « pénalité mécanique » et s'en tient là.

**Preuve erronée de la méthodologie.** §6 écrit « DV ajoutées sans rapport : observé −0,33 à d = 1 et −0,87 à d = 3 ». Dans `niw.R` (l.14), `mu <- rbind(c(.5,0,.3),c(0,.4,0))` : les DV 2 et 3 **portent des effets réels**. L'exemple ne montre donc pas ce qui est affirmé (chiffres reproduits, `niw_rerun.out`).

**Correction.** Mettre à l'échelle r_d = r/√d, ou un a priori direct sur la distance de Mahalanobis, en conservant la réduction à d = 1. Valider sur la grille ci-dessus. Tant que ce n'est pas fait, refuser « Non » en MANOVA et imprimer les BF par DV au même niveau que le global.

## Importants

### I1. Option B de §9 : l'indépendance entre paires est logiquement impossible, et le passage k = 5 → 6 est une rupture
- μ1 = μ2 et μ2 = μ3 impliquent μ1 = μ3. Avec des événements indépendants de probabilité π : P(1=2, 2=3) = π² et P(les trois) = π³, or ces deux événements sont le même. C'est contradictoire pour π < 1. §9 dit « indépendance fausse en pratique » : elle est impossible, pas seulement fausse. Le bon compte si l'on veut une forme simple serait k − 1 (arbre couvrant), soit π0^(1/(k-1)) : 0,707 à k = 3 contre 0,794 annoncé.
- `r1.out` §6 : P(égales) a priori passe de **0,637 à k = 5** (partitions) à **0,955 à k = 6** (Option B). Les cotes a priori sont divisées par **12** pour un groupe de plus.
- `r4b.out` §5 : même paire, écart 0,8 σ, n = 20 par groupe, BF paire = 1,91 : P(μ1≠μ2) = 0,559 ; 0,494 ; 0,298 à k = 3, 4, 5, puis **0,083 ; 0,060 ; 0,046** à k = 6, 7, 8. Le bilan dit donc « égales à 92-95 % » pour un écart de 0,8 σ.

**Correction.** Une seule logique sur tout k : partitions avec a priori à processus de Dirichlet ou beta-binomial, exploration stochastique au-delà de k = 5.

### I2. Option A de §9 : aucune protection une fois le nul global écarté, et deux probabilités pour la même question
- `r1.out` §6 : P(égales) a priori **baisse** avec k au-delà de 4 (0,643 ; 0,637 ; 0,626 ; 0,615 ; 0,606 ; 0,598 à k = 4 à 9), et tend vers 1/2 car B_{k-1}/B_k → 0. La phrase « la multiplicité passe par l'a priori » ne vaut que pour le nul global.
- `r4b.out` §5 : k = 5, 4 groupes égaux et 1 décalé de 1,5 σ, n = 20, 150 jeux : le maximum sur les 6 paires nulles de P(différentes) dépasse **0,5 dans 61 %** des jeux et 0,75 dans 17 %.
- P(T|D) de §5.2 n'utilise que la partition « tous distincts » ; §9 somme sur toutes les partitions non nulles. Sur 150 jeux à k = 4 : écart médian 0,05, maximal 0,089, désaccord autour de 1/2 dans 3,3 % des jeux. Le bilan imprimerait deux réponses à « Le facteur a-t-il un effet ? ». « Compatible avec le π du terme global » n'est vrai qu'a priori.

**Correction.** Quand l'Option A tourne, P(T|D) doit être 1 − P(partition nulle | D).

### I3. r = 1/2 « medium » n'a rien de moyen et son sens dépend de a
**Preuve** : f² = (r²q/a)·F(q,1), vérifié par simulation (`r1.out` §5).

| a | Médiane de f a priori | P(f > 0,4) | P(f > 1) |
|---:|---:|---:|---:|
| 2 | 0,354 | 0,461 | 0,216 |
| 3 | 0,500 | 0,585 | 0,277 |
| 5 | 0,604 | 0,674 | 0,322 |
| 10 | 0,675 | 0,734 | 0,353 |

Dès 3 groupes, plus de la moitié de la masse dépasse le « grand » effet de Cohen (f = 0,4).

**Lecture dans d'autres mondes** (`r5.out` §2, 2 groupes, n = 30, 20 000 jeux) :

| Monde vrai | Annoncé → réalisé |
|---|---|
| Défaut (π = 0,5, Cauchy 0,707) | 0,927 → 0,955 ; 0,615 → 0,617 |
| π = 0,5, δ ~ N(0 ; 0,2²), effets petits | 0,976 → **0,853** ; 0,819 → **0,706** ; 0,221 → **0,450** |
| π = 0,1, δ ~ N(0 ; 0,5²), criblage | 0,925 → **0,611** ; 0,819 → **0,458** |
| π = 0,9, δ ~ N(0 ; 3²), procédé à σ faible | 0,221 → **0,459** ; 0,620 → 0,862 |

La grille de §8 (1:4 à 4:1) ne couvre pas π = 0,1 ni π = 0,9. Elle ne touche pas non plus à la **famille** de l'a priori.

**Correction.** Imprimer la lecture de r en langage d'effet (« médiane a priori de f = 0,5 »). Permettre `effet_attendu =` en unités de y, convertie en r. Élargir les cotes pour un criblage.

### I4. §1.4 interdit l'aiguillage que §5.3 et §5.5 pratiquent
§7 rejette le test préliminaire au nom de [Hoeting1999]. Or §5.3 (« Si P(A:B) > 1/2 … renvoie aux effets simples ») et §5.5 (« Si A:x est crédible ») sont des aiguillages à seuil sur les mêmes données. Et le BF apparié de A *suppose* l'absence d'interaction, c'est-à-dire exactement l'incertitude de sélection ignorée.

**Preuve** (`r5.out` §1, 2×2, n = 12, pI = 0,03) : P(A:B) = 0,690 ; P(A) = 0,266 (apparié), 0,327 (Chipman), 0,447 (uniforme). À P(A:B) = 0,49, le bilan ne dirait rien.

**Correction.** Remplacer le seuil par une moyenne sur l'espace de modèles pour les effets simples.

### I5. `agent_03` : l'échelle « d = 0,5 » est mal nommée et la simulation ne contrôle rien de promis
- Ligne 4 et ligne 7 : moyennes −0,5 ; 0 ; 0,5. Les groupes extrêmes sont à **1,0 σ** ; f = √(1/6) = **0,408**, soit un grand effet de Cohen. L'intitulé « d = 0,5 » de la sortie est faux pour un lecteur. Le commentaire « extremes voisins » n'a pas de sens.
- §10 prévoit une table à « d vrai (0 ; 0,2 ; 0,5 ; 0,8) » sans définir d pour k > 2 : même flou à venir.
- 200 répétitions : P>0,5 sous H0 entre 3,9 % et 11,5 %, P>0,95 sous effet entre 17,8 % et 30 % (déjà imprimé). La ligne compare P>0,5 à p<0,05, qui ne sont pas au même seuil : **P = 0,5 équivaut à p = 0,084** pour ce plan (sortie l.12, retrouvé par quadrature). D'où le 72 % contre 66 %, et le 7 % contre 4 %. Ce n'est pas un gain de puissance.
- Un seul n, un seul plan, aucune simulation sous l'a priori : la propriété promise n'est pas testée.
- `agent_05` §G teste un **autre** a priori : λ fixe et S0 = cov(Y), Y compris l'effet, donc double usage des données. Sa sortie (λ = 1/3 : P>0,5 à 0,435 contre 0,895 pour Pillai) n'est mentionnée nulle part dans METHODOLOGIE, alors que le §6 retient Jeffreys. Le passage de l'un à l'autre n'a aucune trace de script.

**Correction.** Définir l'effet par f, au moins 2 000 jeux, plusieurs n, table à seuils équivalents, et DAP par voie.

### I6. Les contrôles cités en §6 et §5 ne sont pas dans le dépôt
Les chiffres −19,42, −1,076013, −0,07531919, les log BF 0,69632 et 1,084 ne figurent dans aucun script de `dev/mbay/scripts/`. Ils viennent de `niw.R`, `chk.R` et `chk2.R`, dans un dossier temporaire de session. Je les ai rejoués pour `niw.R` : ils se reproduisent. `agent_05` contrôle une formule à a priori propre, pas la formule de §6.

**Correction.** Copier ces scripts et leurs sorties dans `dev/mbay/scripts/` avant toute citation.

### I7. Erreur Monte Carlo : §2 cite les chiffres rassurants, pas les pires
§2 cite 0,9 %, 0,4 % et 23 % (ce dernier imputé à un cas dégénéré). La sortie `agent_02` l.209 et l.212 donne **`propErr max = 0,805`** pour `anovaBF` à n = 100, 3 facteurs, et **0,753** à 4 facteurs (205 s, malgré la limite de 60 s, l.212). À 3 facteurs, c'est le plan du défaut, et le cas n'a rien de dégénéré.

La sortie est aussi tronquée : pas de `FIN`, cas n = 500 à 4 facteurs absent. §8 écrit « les neuf cases coûtent peu » ; à 3 facteurs et n = 500, un `anovaBF` prend 21 s (l.215), soit plus d'une minute pour 3 valeurs de r. À 4 facteurs, plus de 10 min.

**Correction.** Publier le pire cas. Dimensionner l'itération à une erreur cible sur les BF appariés (une somme de BF hérite du pire terme). Pour les plans équilibrés, envisager la quadrature en forme fermée.

### I8. Le ROPE n'est pas cohérent entre termes
- **Covariable contre facteur** : un facteur à 2 niveaux codé ±0,5 a sd(x) = 0,5. ROPE covariable |β·sd(x)/σ| < 0,1 donne |δ| < **0,2** ; ROPE facteur donne |δ| < **0,1**. Les a priori, eux, sont cohérents (r = 0,354 ↔ δ ~ Cauchy 0,707). Même variable, seuil double selon son codage.
- **§4 contre §6** : à d = 1 et a = 2, §4 donne f = |δ|/2 et §6 (f² = tr/q) donne |δ|/√2. La définition de §6 divise par q, celle de §4 par a.
- **Nombre de groupes** : un lot déviant de Δ donne f = Δ·√((1/a)(1−1/a)). Le seuil tolère donc Δ < 0,1 à a = 2 ; 0,167 à a = 10 ; 0,32 à a = 40.
- **σ résiduel** : même effet brut 0,15 : δ = 0,027 sans covariable, **0,291** avec (`r5.out` §5). « Négligeable » change avec une covariable et n'a pas de lien avec un seuil métier.
- **« Toujours »** : §1.2 dit « ROPE toujours » ; §6 ne donne « aucun défaut », §5.9 non plus.

**Correction.** Seuil de covariable à 0,05 sur la même métrique f ou seuil 0,1 partout, f unique, unités brutes pour un industriel.

### I9. Sphéricité : pas de correction, et la « voie alternative » de §5.7 ne marche pas telle quelle
**Preuve numérique** (`r6.out`, 100 jeux, T = 6, n = 20, H0).

| Covariance | P>0,5 | P>0,95 | F non corrigé <0,05 |
|---|---:|---:|---:|
| Composée symétrique | 0,01 | 0,00 | 0,00 |
| AR(1) 0,80 | **0,13** | 0,02 | 0,12 |
| AR(1) 0,95 | **0,12** | 0,02 | 0,09 |

**Preuve algébrique.** §5.7 renvoie à « MANOVA sur les contrastes intra (§6) ». Tester un effet intra revient à tester **μ = 0 sur les contrastes**. Or §6 pose μ plat et centre Y : le paramètre à tester est celui qu'on a intégré. La formule ne s'applique pas à ce cas. Il faudrait une autre dérivation, sans centrage, à ν = n et a priori g sur μ, donc une autre constante C. Seule l'interaction inter × intra est couverte par §6.

**Correction.** Écrire et contrôler cette dérivation. D'ici là, signaler l'absence de contrôle de sphéricité sur la ligne.

### I10. Pentes aléatoires : le cas continu échoue et n'est pas documenté
`agent_02` l.92-95 : `x:id` en `whichRandom` donne **`NA ± NA%`** avec un simple avertissement ; §5.8 ne parle que de A:id, facteur, avec un usage « non décrit par le manuel ». Cas fréquent en mixte. §5.8 dit aussi « défaut : pentes présentes » en s'appuyant sur [Barr2013], travail fréquentiste, et prend [Heck2021] comme simple sensibilité alors que ce dernier plaide pour un test par moyenne de modèles.

**Correction.** Refus explicite et voie Stan annoncée. Justifier le choix du défaut contre [Heck2021] au lieu de le citer à l'appui.

### I11. « Non » a un plafond fixé par n, et « P » n'est qu'un recodage de p à n donné
**Preuve** (`r1.out` §7, 3 groupes, r = 0,5) : même à F = 0 exactement.

| n/gr | P(H0) maximum | P(H1) à p = 0,05 | p correspondant à P(H1) = 0,5 |
|---:|---:|---:|---:|
| 3 | 0,720 | 0,676 | 0,179 |
| 15 | 0,859 | 0,598 | 0,084 |
| 30 | 0,909 | 0,532 | 0,058 |
| 100 | 0,965 | 0,359 | 0,027 |
| 1000 | 0,996 | 0,071 | 0,0037 |

- La corrélation rang(P, −p) = 1 dans `agent_03` n'est pas une trouvaille : à plan équilibré à 1 facteur, BF est une fonction de F. Le « Oui p » n'ajoute à la p-value que la dépendance à n et aux a priori. L'argument commercial doit le dire.
- « Stable » selon §8 exige BF hors [1/4 ; 4] sur la grille de cotes. À 5/gr, le BF minimal est 0,31 : un « non » stable y est **impossible**, quelles que soient les données.
- « Non 0,86 » à n = 45 est le **maximum atteignable** et doit être présenté comme tel.

**Correction.** Imprimer le plafond sur la ligne ; à petit n, expliquer qu'« absence de preuve » est la seule lecture.

### I12. Beta-binomial au-delà de 3 facteurs : non défini sous marginalité, et rupture au seuil
[ScottBerger2010] porte sur des inclusions échangeables ; sous marginalité, les termes ne le sont pas. §3 ne dit pas comment l'adapter, ni que P(A) a priori ne sera plus 1/2 alors que le bilan doit annoncer π. Comme pour les paires, un facteur de plus change brutalement le chiffre sur les mêmes données. Calcul non réalisé : c'est l'absence de définition qui est constatée.

**Correction.** Définir une hiérarchie (type Chipman à hyperparamètre Beta), et ne pas mettre de seuil.

### I13. Covariable : §2 et §5.4 pris à la lettre donnent un BF dépendant de l'unité
§2 écrit β_T | σ², g ~ N(0, σ²g I), g ~ IG(1/2, r²/2) pour « chaque terme », et §5.4 dit « Covariable centrée » (pas réduite).

**Preuve** (`r1.out` §8) : avec cette écriture, BF = 0,639 ; **0,001** ; 1,000 pour x, x·1000 et x/1000, contre 0,580 et 0,580 pour `lmBF`. La forme juste est celle de `agent_02` l.42 (g ~ IG(1/2, r²n/2), invariante, type Liang). Piège immédiat à l'implémentation, y compris pour une MANOVA avec covariable.

**Correction.** Écrire la covariable réduite dans la formule.

## Mineurs

- **M1. §6 point 4 : « l'intégrale sur g peut diverger » est faux.** Avec k = q + d − (n−1) valeurs θ = 1, l'intégrande se comporte en g^(−3/2 − (q−k)(d−k)/2) : l'intégrale converge toujours. Si k = min(q, d), **BF = 1 exactement** quelles que soient les données (`r1.out` §4 : logBF = 0 à g = 10² … 10⁸). Sinon, BF fini mais dégénéré (0,54 et 0,82). Refuser reste juste, pour la raison inverse : le nombre ne dépend pas ou guère des données.
- **M2. §5.2 : « son contenu marginal sur une paire change avec le nombre de niveaux » est faux** pour ce codage. Var(α_i − α_j)/g = 2 pour tout a (`r1.out` §5), donc Cauchy(0 ; 0,707) à tout a. Ce qui change, c'est f (I3).
- **M3. §3 : « A + B contre B est exactement la comparaison du Type II ».** Le BF apparié moyenne aussi A contre 1 : `agent_02` donne 20,49 (Type II exact) contre 18,7. Parenté, pas identité.
- **M4. §5.6 contre §5.7 :** deux chiffres pour un facteur intra à 2 niveaux : `ttestBF` 0,252 contre `lmBF` 0,322 ; 2,005 contre 2,167 (`r5.out` §3), et une échelle différente (σ_d contre σ résiduel, facteur √2). Le plan à 2 conditions n'est pas tranché.
- **M5. Libellé « NIW »** pour un a priori impropre (Jeffreys), et formulation qui varie entre §1.3, §6 et le tableau des options. La condition « même centrage » doit préciser : X centré *dans chaque modèle*.
- **M6. Vocabulaire contradictoire** : §1.5 « toujours imprimée » contre `sensitivity = FALSE` (§8) ; §7 « les données ne tranchent pas » contre §8 « conclusion dépendante des a priori » pour le même cas.
- **M7. Bibliographie :** voir section C.

---

# B. Jugements (non démontrés par le calcul)

1. **Le symbole p (important).** Dans un paquet dont `m.test` imprime des p-values, écrire « p = 0,91 » veut dire l'inverse pour les mêmes étudiants. Écrire « Prob. = 0,91 », jamais p.
2. **« Non » veut dire « exactement nul contre une Cauchy d'échelle r » (important).** Un industriel lira « sans effet ». Pour f = 0,20 à n = 30, P(H0) ≈ 0,55 (`r3.out`). Écrire « Pas de preuve d'un effet de la taille attendue » plutôt que « Non ».
3. **H0 ponctuelle par défaut en industrie (important).** À 10 000 par groupe et f = 0,03, la ligne de tête dit « Oui » avec 0,99 ou plus (`r3.out`). Un procédé observationnel à grand n devrait avoir le ROPE en tête. §4 le reconnaît en parlant d'observationnel, mais le défaut ne dépend de rien : c'est une absence d'aiguillage, et ici un aiguillage sur le *plan* (pas sur les données) serait légitime.
4. **Standardiser par σ est mauvais pour un métrologue (important).** Une machine qui s'améliore (σ baisse) rend les mêmes écarts « non négligeables » : écart de 0,002 mm, σ = 0,001 mm, tolérance ±0,5 mm donne δ = 2, « effet non négligeable ». §4 dit qu'une valeur métier vaut mieux ; le défaut devrait la demander quand les unités sont connues, et le bilan dire que 0,1 est une convention *à chaque fois*.
5. **« Stable » est mal choisi (mineur).** Une plage 0,51-0,94 est « stable » selon §8, alors que 0,51 n'affirme rien. Imprimer la plage sans le mot.
6. **Le niveau 3 `robust = "bma"` (important).** L'argument « a priori impropres communs » ne s'étend pas à gaussien contre hétéroscédastique : σ et (σ_1 … σ_k) ne sont pas le même paramètre. Il faudra des a priori propres sur log σ_j, et P(ℓ|D) héritera de leur largeur (Lindley, §6 l'écrit pour B). Rien dans §7 ne le dit. Non calculé ici.
7. **Refus de la cohérence (important).** §11.5 présente comme un choix l'acceptation d'un p incohérent contre une cohérence « moins lisible ». Vu B1 et B2, la lisibilité invoquée est une lisibilité de l'erreur.
8. **L'a priori commun aux 7 termes** ignore que les interactions d'ordre élevé sont souvent plus petites et plus rares. Un r unique pour les ordres 1, 2, 3 se défend mal en plan d'expériences (hiérarchie des effets, parcimonie).
9. **« Voici la sensibilité » rassure à tort** : la grille varie r et π, pas la vraisemblance, pas la famille de l'a priori, pas apparié contre cohérent. Or B1 montre que ce dernier axe peut faire passer de 0,24 à 1,00. Ajouter au moins la ligne cohérente à la grille.

---

# C. Bibliographie : 8 références tirées au hasard

Tirage reproductible (`set.seed(20261105)`, `tirage.R`) parmi 25 références dont l'apport est précis.

| Référence | Apport annoncé | Source consultée | Verdict |
|---|---|---|---|
| [Gelman2008] | §5.9 : Cauchy(0 ; 2,5) « sur prédicteurs standardisés » | Résumé arXiv 0901.4011 | **Inexact.** Échelle 2,5 confirmée, mais prédicteurs non binaires à **écart-type 0,5**, non 1. Sur des prédicteurs réduits à 1, l'équivalent est 1,25. Facteur 2 sur l'a priori en logistique : **à corriger** |
| [Chipman1996] | « Hérédité forte : π(A:B \| A, B) = 1/2 » (§3) | Résumé arXiv bayes-an/9510001, notices | **Partiel.** Le résumé traite d'a priori entre prédicteurs liés. Il ne parle ni d'hérédité forte ni de 1/2, et le mot « heredity » vient plutôt de Chipman, Hamada & Wu (1997). Le 1/2 est un réglage KefiR, non sourcé : le dire |
| [Barr2013] | « Omettre les pentes gonfle les faux positifs » | Résumé Europe PMC | **Compatible, non établi par le résumé**, qui parle de moins bonne généralisation. À vérifier dans le texte. Travail fréquentiste : son usage pour BF est une extrapolation |
| [Kruschke2013] | « Student-t à variances inégales, a priori sur ν » | Résumé Europe PMC | **Non établi par le résumé**, conforme à ma connaissance du modèle, non vérifiée ce jour. Le travail plaide pour l'estimation, pas pour des BF : le citer pour une branche BF demande une réserve |
| [Heck2021] | « Pas de comparaison unique : moyenne sur versions » | Notice et résumé par recherche ; PDF refusé | **Conforme**, mais METHODOLOGIE en fait une sensibilité, alors que l'article défend cette voie (I10). Date : parution vol. 6, 2023 ; la bibliographie écrit 2021 sans la mention « en ligne » utilisée ailleurs |
| [Langsrud2003] | Argument pour Type II | Notices (résumé partiel) | **Conforme** |
| [Peixoto1990] | Invariance des modèles bien formés au recodage | Recherche (résumé de seconde main) | **Conforme**, sur des polynômes |
| [Mulder2019] | BF fractionnaire, Cauchy matriciel par défaut | Résumé arXiv | **Conforme** |

Sur 8 : 1 erreur de fond, 1 attribution non sourcée, 2 non vérifiables au résumé, 4 conformes. Cela confirme l'avertissement de BIBLIOGRAPHIE.md (71 références sur 92 en « métadonnées ») : ce niveau ne permet pas de citer un point précis.

---

# D. Attaques tentées qui ont échoué

Consigné pour éviter de les refaire.

- **Formule §6** : log m(Y|g), forme équilibrée, réduction à d = 1 : corrects (Monte Carlo indépendant sous le postérieur de M0, 0,80332 contre 0,80385 ; `anovaBF` à 10⁻⁵ près sur plan déséquilibré).
- **C commun** : tient pour des modèles à g différents : C ne dépend que de n, d, et la constante du Jeffreys s'élimine, sous conditions : même d, même n, X centré dans chaque modèle, Σ de même rôle. L'argument ne s'étend pas à §5.7 (I9) ni à `robust = "bma"` (B.6).
- **DAP MANOVA** : moyenne 0,493 (se 0,009) à d = 3 ; 0,509 à d = 5.
- **Chiffres** : 0,089 ; 0,545 ; 0,625 ; 0,643 ; π0^(1/m) ; 0,6/0,2 ; 0,737/0,316/0,053 ; Sellke 0,41/0,29 : exacts.
- **Lot déviant noyé** : je pensais que le g partagé masquerait un lot isolé. Non : P médiane 0,996 à a = 20, Δ = 2σ, n = 10, et ≥ 0,82 à Δ = 1,5σ jusqu'à a = 20. Le masquage n'apparaît qu'à a = 40, Δ = 1,5σ (P médiane 0,36 ; p = 0,035).
- **Σ non commune en MANOVA** : très mauvais (P>0,5 à 0,427), mais Pillai fait pareil (0,423). Pas un défaut propre à m.bay.
- **t₂ sous H0** : P>0,5 à 0,05. Pas d'inflation.
- **Calibration dans le monde par défaut** : correcte en 2 groupes.

---

# E. Ordre de reprise proposé

1. Trancher B1 et B2, car ils décident de la phrase imprimée : défaut cohérent ou apparié conditionnel explicite.
2. Refaire le ROPE (B3, I8), avant toute table.
3. Rendre disponible une voie robuste sans Stan, ou un refus (B4).
4. Remettre r à l'échelle de d en MANOVA (B6), et écrire la dérivation sur μ (I9).
5. Unifier les post-hoc (I1, I2).
6. Faire ensuite la table de §10, avec f défini et les seuils équivalents (I5, B5).
