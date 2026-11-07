# Méthodologie de m.bay

Pendant bayésien de `m.test`. Pour chaque terme : « Ce facteur a-t-il un effet ? Oui avec probabilité p ; non avec probabilité 1-p », a priori annoncés, sensibilité affichée.

Les clés entre crochets renvoient à BIBLIOGRAPHIE.md, qui indique si le contenu a été lu ou seulement les métadonnées vérifiées.

**Contrôles numériques faits pour ce document** (R 4.6.1, BayesFactor 0.9.12-4.8, 2026-11-05, un jeu de données simulé par cas). Ils portent uniquement sur l'algèbre et non sur les performances : voir §3 et §7.

> **État au 5 novembre 2026 : proposition de départ, contestée.** La relecture adverse (`RELECTURE_ADVERSE_2026-11-05.md`) y trouve 6 défauts bloquants, ouverts, listés dans `DECISIONS.md`. Les positions du §1 ne sont pas validées, en particulier le « π = 1/2 apparié » (§3), le ROPE à plusieurs ddl (§4), la branche gaussienne sans aiguillage (§7), la multiplicité entre termes et le réglage MANOVA (§6). Les parties §3 et §4 seront réécrites après arbitrage d'Antoine.

## 1. Positions en bref
1. **P(terme)** *(en attente de D2)* : ligne 1, marginale d'une loi jointe unique, « le facteur joue un rôle, seul ou en interaction » ; ligne 2, effet principal apparié, conditionnel et jamais lu « non ». L'ancienne position (π = 1/2 par terme, apparié) ne tient pas : voir §3 et le JOURNAL du 5 nov.
2. **H0 ponctuelle** en première ligne ; **ROPE** en seconde ligne, seulement sur les quantités à 1 ddl *(en attente de D3)*. Ne jamais les fusionner.
3. **MANOVA** : BF analytique Normal-Inverse-Wishart à a priori g par terme et Jeffreys sur Σ, intégré sur g. Réduit exactement à `anovaBF` à d = 1 (contrôlé).
4. **Pas de test préliminaire qui choisit la branche.** Le gaussien est toujours calculé ; le robuste s'y ajoute, et un désaccord s'affiche au lieu d'être arbitré en silence.
5. **Sensibilité** toujours imprimée : grille r × cotes a priori 3 × 3, plus la plage des p.
6. **Post-hoc** : partitions de moyennes pour k ≤ 5, sinon paires à cotes a priori ajustées.
7. **Promesse de calibration** limitée à ce qui est démontrable : calibration sous l'a priori. Le taux de faux positifs fréquentiste se mesure par simulation et se publie.

## 2. Notations et moteur commun

Modèle gaussien :

    y = μ + Σ_T X_T β_T + ε,   ε ~ N(0, σ²I),   p(μ, σ²) ∝ 1/σ²

Chaque terme T a son propre g :

    β_T | σ², g_T ~ N(0, σ² g_T I_{q_T}),   g_T ~ InvGamma(1/2, r_T²/2)

avec X_T en codage réduit orthonormal (q_T = a-1 pour a niveaux), comme [Rouder2012]. Marginalement, β_T/σ suit une Cauchy multivariée d'échelle r_T.

**Échelles par défaut** (celles de [BayesFactorManual], pour rester comparable à JASP et à la littérature) :

| Terme | r | Lecture |
|---|---|---|
| Facteur fixe | 1/2 | « medium » de `anovaBF` |
| Deux groupes, sur δ | √2/2 ≈ 0,707 | Même a priori que r = 1/2 en réduit (δ = √2 β*). Contrôlé : `ttestBF` et `anovaBF` donnent le même log BF (−0,6246) |
| Covariable continue | √2/4 ≈ 0,354 | « medium » de `lmBF` |
| Effet aléatoire | 1 | « nuisance » |

Piège d'affichage : « medium » n'a pas la même valeur selon la fonction. `m.bay` imprime toujours la valeur numérique et l'objet qu'elle met à l'échelle.

**Passage du BF à la probabilité.** Avec des cotes a priori O = π/(1-π) :

    P(T | D) = O·BF / (1 + O·BF)

Choix de calcul : voies analytiques ou `BayesFactor` en routine (rapide), Stan + bridge sampling [Carpenter2017, Burkner2017, Gronau2017, Gronau2020] seulement quand le modèle l'exige. Chaque p porte son erreur numérique. Essais `lmBF` : erreur relative de **23 %** sur un rapport à pentes aléatoires sans réplique (une observation par sujet × condition, cas dégénéré, §5.7) ; 0,9 % avec 3 répliques ; 0,9 % puis 0,4 % sur A + B avec 10⁴ puis 5·10⁵ itérations. Il faut donc viser une erreur cible, l'afficher, et refuser le cas sans réplique.

## 3. Définir P(terme a un effet | données) *(en attente de D2)*

> **En attente de la décision D2 d'Antoine.** Loi jointe non choisie ; forme finale non simulée.

### Deux lignes, deux questions

La question d'Antoine est « ce facteur joue-t-il un rôle ? », pas « son effet principal moyen est-il non nul ? ». Le bilan donne donc deux lignes par facteur, et le mot « Non » ne s'imprime sur aucun effet principal apparié.

**Ligne 1 : le facteur joue-t-il un rôle, seul ou en interaction ?** Marginale d'une loi jointe unique π(m) sur les modèles respectant la marginalité :

    Prob(A joue un rôle | D) = Σ_{m contient A} π(m) p(D|m) / Σ_m π(m) p(D|m)

Sous la marginalité, tout modèle contenant une interaction impliquant A contient A : l'événement est bien « seul ou en interaction ». Les probabilités se combinent, et Prob(A:B | D) ≤ Prob(A | D) pour toute loi [relecteur, B1]. Dans le croisement 2×2 pur, cette ligne vaut 1,000 [recalculé, `agent_11`].

**Ligne 2 : l'effet principal apparié, à lire sous condition.**

    BF_A|c = [ Σ_{m∈M+} π(m) p(D|m) / Σ_{m∈M+} π(m) ] / [ Σ_{m∈M−} π(m) p(D|m) / Σ_{m∈M−} π(m) ]

M+ : modèles avec A sans interaction impliquant A ; M− : les mêmes sans A. Parent du Type II, sans lui être identique : 18,7 contre 20,49 pour le Type II exact [relecteur, M3].

Elle s'imprime toujours avec sa condition et la probabilité de cette condition, **sans seuil qui la masque** (un seuil serait l'aiguillage rejeté au §7, voir I4 de la relecture) :

    Effet principal, s'il n'y a pas d'interaction impliquant A (Prob. de cette condition : 0,00) : BF = 0,27

Elle donne « données compatibles avec un effet principal nul sous cette condition » et jamais « A n'a pas d'effet ». Dans le croisement ci-dessus, l'ancien défaut imprimait 0,215 et son complément en « Non » [recalculé, `agent_11`] (0,236 sur un autre tirage, [relecteur]) : c'est ce que cette ligne doit rendre impossible. Le BF affiché ci-dessus est fictif.

Si l'interaction a un poids non négligeable, le bilan donne aussi l'effet de A aux niveaux de B, moyenné sur les modèles, en quantités à 1 ddl (§4).

### Les deux lois candidates

**Arbitré par Antoine le 6 novembre : hérédité forte par défaut, l'uniforme reste en sensibilité** (`DECISIONS.md`, D2). Le libellé de la question par terme reste ouvert.

Modèles à 2 facteurs [main ; 0,6 et 0,2 : relecteur, agent]. **Le choix fixe tous les nombres : à annoncer sur chaque ligne.**

| Modèle ou quantité | Uniforme | Hérédité forte, π(A) = π(B) = 1/2, π(A:B \| A, B) = 1/2 |
|---|---:|---:|
| Nul | 1/5 | 1/4 |
| A | 1/5 | 1/4 |
| B | 1/5 | 1/4 |
| A + B | 1/5 | 1/8 |
| A + B + A:B | 1/5 | 1/8 |
| **P(A) a priori** | **0,6** | **1/2** |
| **P(A:B) a priori** | **0,2** | **1/8** |
| P(A \| pas de A:B) | 1/2 | **3/7 ≈ 0,429** |
| P(A:B \| A, B) | 1/2 | 1/2 |

**Trois points à ne pas perdre.**
- Le BF apparié avec cotes 1:1 **est** le conditionnel de l'uniforme [relecteur, B2 ; recalculé, D2]. Sous hérédité forte, la cote conditionnelle sur la ligne 2 n'est **pas** 1:1 mais 3:4 [main], et les poids de M+ sont inégaux (2:1), donc le BF_A|c change aussi. **Annoncer « π = 1/2 » sur la ligne 2 n'est donc vrai que sous l'uniforme.**
- L'uniforme ne part pas de 1/2 et favorise les effets principaux : 0,737 à 3 facteurs, 0,316 et 0,053 pour les interactions [agent, relecteur].
- L'hérédité forte n'est **pas publiée telle quelle** : l'idée vient de la littérature sur les prédicteurs liés [Chipman1996], le réglage est propre à KefiR.

**À 3 facteurs, l'hérédité forte n'est pas encore définie de façon unique.** À trancher :

| Variante | P(A) | P(A:B) | P(A:B:C) | P(nul) | Statut |
|---|---:|---:|---:|---:|---|
| π = 1/2 à chaque ordre (celle de D2) | 1/2 | 1/8 | 1/128 | 1/8 | [main], non programmée |
| π^k à l'ordre k (`agent_20`) | 1/2 | 1/16 | 1/4096 | 1/8 | [exécuté par l'agent] |

**La seconde donne 1/16 à 2 facteurs, pas le 1/8 annoncé dans D2** : les chiffres d'`agent_20` étiquetés « hérédité » ne portent donc pas sur la loi annoncée. L'uniforme y a pour a priori P(nul) = 1/19 ≈ 0,053.

### Multiplicité entre termes (D13)

Prob > 0,5 **n'est jamais un verdict**, ni dans le texte du bilan, ni dans un test unitaire. Le bilan annonce :
- le nombre de termes testés ;
- Prob(modèle nul | D), avec sa valeur a priori ;
- l'espérance du nombre de termes actifs a posteriori.

Motif, recoupé sur `anovaBF` le 6 novembre (2³ sous H0, 7 termes, au moins un terme à P>0,5) : apparié 0,74 à N=16 (150 jeux), 0,58 et 0,60 à N=40 (226 et 300 jeux) ; uniforme 0,75 à N=16, 0,43 et 0,49 à N=40 ; hérédité de D2, 0,25 à N=40 (300 jeux). P>0,95 : au plus 1,7 % [`agent_20` par un agent, `agent_22` lancé par le coordinateur]. L'hérédité π^k d'`agent_20` n'est pas la loi de D2 : ne pas citer ses 0,34 et 0,24. L'uniforme n'est pas un remède : chaque effet principal part de 0,737. Beta-binomial à partir de 3 facteurs : proposition de D13, **non définie sous marginalité** [relecteur, I12], à valider.

### Rapport à `$effects` de m.test

`$effects` donne une p-value par terme, de façon non homogène : Type III pour les mixtes, Type II pour `car::Manova` et en repli pour certaines ANCOVA, test global pour Kruskal et Friedman.

- `$effects` : Prob(statistique au moins aussi extrême | H0, modèle). `m.bay` : Prob(T | D, a priori, espace de modèles). Pas de conversion [WassersteinLazar2016, BergerSellke1987].
- Mêmes noms de termes, colonne `analogue`, pour la lecture côte à côte. **La ligne 2 est la seule à comparer à une p-value de Type II** ; la ligne 1 répond à une autre question.
- Un p ≈ 0,04 peut donner un BF ≈ 1 : à p = 0,05, BF01 ≥ 0,41, donc P(H0|D) ≥ 0,29 avec π = 1/2 [Sellke2001, formule de mémoire ; chiffres : relecteur].

Le symbole imprimé est « Prob. », jamais « p » [relecteur, jugement B.1].

### Options

- `whichModels = "top"` (plein moins T, analogue Type III) : option avertie, viole la marginalité pour un effet principal.
- Modèle à probabilité médiane [BarbieriBerger2004] : applicable à la ligne 1, pas à la ligne 2.

## 4. H0 ponctuelle, ROPE et autres indicateurs *(en attente de D2 et D3)*

### Six quantités à ne pas confondre

| Quantité | Répond à | Conditionnée à | Dépend de π ? | Ne dit pas |
|---|---|---|---:|---|
| **BF10** | Combien les données favorisent H1 contre H0 | Deux modèles, a priori sur les paramètres | Non | Si H1 est probable |
| **P(H1 \| D)** | Quelle probabilité que l'effet soit non nul | Modèles et π | **Oui** | La taille de l'effet |
| **P(\|θ\| < w \| D, H1)** | Si l'effet existe, est-il petit ? | H1 et w | Non | Si l'effet existe |
| **P(ROPE \| D)** = P(H0\|D) + P(H1\|D)·P(\|θ\|<w \| D, H1) | Quelle probabilité que l'effet soit négligeable | H0, H1, π et w | **Oui** | Elle mélange atome et partie continue : à ne pas lire seule |
| **pd** = max[P(θ>0 \| D, H1), P(θ<0 \| D, H1)] | Quel signe, si l'effet existe | H1 | Non | Ne teste pas H0. Un pd de 0,99 coexiste avec un P(H0\|D) élevé. Lien de Makowski p = 2(1−pd) : ne pas présenter comme une équivalence en dehors de ses conditions |
| **Taille d'effet et intervalle de crédibilité** | Quelle valeur, à quelle incertitude | H1 sauf mention, a priori | Non | **L'intervalle à 95 % n'est pas un test.** Il exclut 0 sans implication sur P(H0\|D). Sous H1 seul, il ignore l'atome en 0 |

Toute grandeur « sachant H1 » est imprimée comme telle.

### Domaine du ROPE : 1 ddl seulement

**ROPE seulement sur les quantités à 1 ddl** : deux groupes, paire de moyennes, contraste nommé, pente, effet simple. **Aucun ROPE sur un terme à plusieurs ddl** tant qu'aucune définition valide n'existe.

Motif : f = sqrt(mean(α_j²))/σ est positif par construction et biaisé. Sous vérité nulle à a = 5, n = 100, médiane de f | H1 = 0,105, le double du seuil 0,05 ; f vrai = 0,20 à n = 30 donne « négligeable 0,55 », dont l'essentiel vient de l'a priori [relecteur, B3, non reproduit]. Sur un terme à plusieurs ddl, le bilan renvoie à la liste des paires ou contrastes. Pistes non retenues à ce jour : f débiaisé, H0 par intervalle [MoreyRouder2011]. MANOVA : pas de ROPE.

### Largeurs

| Quantité | ROPE par défaut | Masse a priori sachant H1 | Remarque |
|---|---|---:|---|
| 2 groupes, δ | \|δ\| < 0,1 | 0,089 [relecteur] | ±0,1 : la moitié d'un petit effet de Cohen [Kruschke2018, texte lu] |
| Paire dans un facteur à k niveaux | \|α_i − α_j\|/σ < 0,1 | 0,089 pour tout k, à r = 1/2 [main ; variance : relecteur, M2] | Tant qu'il y a un g commun |
| Pente | **\|β·sd(x)/σ\| < 0,05** | 0,089 à r = √2/4 [main : (2/π)·atan(0,05/0,354)] | **Corrigé.** [Kruschke2018] donne ±0,05 sur variables standardisées (changement de ±0,1·S_y sur ±2·S_x). L'ancienne valeur 0,1 n'en vient pas. Ainsi, un facteur à 2 niveaux codé ±0,5 retrouve \|δ\| < 0,1 [main], ce qui corrige le premier point de I8 |

Trois réserves :
- **Échelle.** [Kruschke2018] et `bayestestR` standardisent par S_y ; m.bay par σ résiduel. Avec une covariable forte, les chiffres divergent (δ = 0,027 contre 0,291 pour un même effet brut [relecteur, I8]).
- **Règle.** Kruschke décide avec l'intervalle à 95 % et le ROPE, sans atome. La probabilité sur la loi entière vient de [bayestestRrope, Makowski2019]. Le mélange avec un atome en 0 est notre construction.
- **Convention.** 0,1 n'est pas un résultat. Kruschke lui-même prévient que ce seuil dépend de l'application. La largeur est remplaçable par une valeur métier, qui vaut mieux.

### Quelle ligne, dans quel ordre

**Défaut proposé** : H0 ponctuelle en première ligne, puis, pour chaque quantité à 1 ddl, atome et partie continue **séparés** :

    A vs B : effet non nul ?       Prob. = 0,91 ; a priori 0,50   [H0 : δ = 0 ; BF = 10 ; r = 0,707]
             si non nul, petit ?   Prob.(|δ| < 0,1 | H1) = 0,04 ; a priori 0,089
             signe, si non nul     Prob.(δ > 0 | H1) = 0,99
             taille, si non nul    δ = 0,55 [0,12 ; 0,98]

(chiffres **fictifs**, sauf 0,089 et le fait que BF = 10 donne 0,909 avec π = 1/2)

Le mot « négligeable » ne s'imprime que si la probabilité dépasse **nettement** sa valeur a priori ; sinon « peu d'information sur la taille » [relecteur, B3 ; aucune définition de « nettement » n'est validée]. Option `null = "rope"` pour inverser l'ordre. Ce choix peut légitimement dépendre du **plan** déclaré (observationnel à grand n), et non des données [relecteur, jugement B.3].

### Plafond de la ligne « non »

Le plafond de P(H0|D) dépend de n, même à F = 0 : 0,720 à n = 3 par groupe, 0,859 à 15, 0,909 à 30 (3 groupes, r = 0,5) [relecteur, I11]. Le bilan imprime ce plafond ; à petit n, « absence de preuve » est la seule lecture.

### Affirmations chiffrées à exécuter avant adoption

| Affirmation | État |
|---|---|
| 1/2, 1/8, 3/7 à 2 facteurs ; 1/8, 1/128, 1/8 à 3 facteurs | [main] |
| Ligne 1 et ligne 2 sous hérédité forte (forme finale) | [non exécuté] |
| Calibration des deux lignes, y compris DAP | [non exécuté]. Le relecteur a le témoin P(A) cohérent sous l'uniforme ; rien sous la forme finale |
| 0,089 pour une paire à k niveaux et pour la pente | [main] |
| Pente à 0,05 : cohérence avec le facteur | [main] ; aucune simulation |
| Masse et calcul de P(\|θ\|<w \| D, H1) à plusieurs facteurs | [non exécuté]. Quel modèle ? La moyenne sur les modèles contenant le terme reste à définir et à contrôler |
| Biais de f, 0,105, 0,55 | [relecteur], non reproduits |

## 5. Plans

Pour chaque plan, `P(T|D)` suit le §3 et la ligne ROPE le §4. Seul ce qui change est indiqué.

### 5.1 Deux groupes

- **Modèle** : y_i = μ ± δσ/2 + ε. Cauchy(0, 0,707) sur δ (échelle « medium » de [BayesFactorManual] ; [Rouder2009] pose la famille JZS, son échelle reste à confirmer dans le texte).
- **Probabilités** : P(δ ≠ 0 | D) ; P(|δ| ≥ 0,1 | D) ; signe, P(δ > 0 | D, H1).
- **Limite** : à petit n, le résultat dépend de r. Le défaut ne remplace pas un a priori informé si l'utilisateur en a un [Gronau2020b].

### 5.2 Trois groupes ou plus

- **Modèle** : §2, un g pour le facteur. Forme exacte contrôlée : la formule bloc reproduit `anovaBF` déséquilibré (log BF 0,69632 contre 0,69632).
- **Probabilité** : P(au moins deux moyennes diffèrent | D). Terme global, ne dit pas lesquelles (§9).
- **Limite** : l'a priori lie tous les niveaux par un même g, et son contenu marginal sur une paire change avec le nombre de niveaux.

### 5.3 Plusieurs facteurs et interactions

- **Modèle** : un g par terme, interactions comprises [Rouder2012, Rouder2017]. `whichModels = "withmain"`.
- **Probabilités** : §3. Pour un effet principal : lecture « A ajoute, B présent, sans interaction A:B ». Si P(A:B) > 1/2, le bilan signale que la lecture de A seul est conditionnelle à un modèle peu crédible et renvoie aux effets simples.
- **Contrôle** : additif A + B déséquilibré, formule à deux g contre `lmBF` : log BF 1,084 contre 1,078 (erreur 0,4 % à 5·10⁵ itérations). Concordant à l'erreur Monte Carlo près.
- **Limite** : l'espace explose avec le nombre de facteurs (5, puis 19 modèles à 2 et 3 facteurs). Plan très déséquilibré ou cellules vides : modèles non identifiables, refuser de répondre sur les termes concernés.

### 5.4 ANCOVA

- **Modèle** : y = μ + X_A β + x γ + ε. Covariable centrée, r = 0,354 [RouderMorey2012].
- **Probabilité** : BF(A + x vs x). La covariable est présente dans les deux modèles, avec un a priori identique, comme nuisance. P(A | D, x), pas P(A).
- **Limites** : le BF sur A dépend de la covariable retenue. Une covariable choisie après coup invalide la lecture. Linéarité supposée. Covariable mesurée après traitement : aucun cadre statistique ne le corrige.

### 5.5 Pentes hétérogènes

- **Modèle** : ajout de A:x. **Probabilité** : BF(A + x + A:x vs A + x).
- **Si A:x est crédible** : « effet de A » n'existe plus comme quantité unique. Le bilan donne l'effet de A à des valeurs de x (moyenne, ±1 sd) par tirages a posteriori, non un P(A) unique.
- **Limite** : `lmBF` accepte A:x sans erreur (essai fait), mais **l'échelle a priori de ce terme reste à documenter** avant d'annoncer quoi que ce soit, le manuel consulté ne la précisant pas. Sinon brms avec a priori explicites.
- **Lien m.test** : `sys_ancova_analysis.R` teste l'interaction facteur × covariable de façon hiérarchique. m.bay la garde dans l'espace de modèles.

### 5.6 Apparié

- **Modèle** : sur les différences d_i : d = δσ_d + ε, Cauchy(0, 0,707) (même remarque qu'au §5.1).
- **Limite franche** : δ est standardisé par σ_d, pas par l'écart-type intra-groupe. Même r = cible différente selon la corrélation. Le ROPE doit préciser cette échelle, sinon un δ = 0,1 apparié et un δ = 0,1 indépendant ne sont pas comparables.

### 5.7 Mesures répétées et split-plot

- **Modèle** : sujet aléatoire présent partout (r = 1), facteurs fixes r = 1/2 [Rouder2017]. Split-plot : facteur inter, facteur intra, et leur interaction.
- **Défaut** : pentes aléatoires quand le plan a plusieurs observations par sujet × cellule. Avec **une seule** observation par sujet × condition, A:ID se confond avec le résidu : pas de pente estimable.
- **Probabilité** : §5.8 à plusieurs observations ; sinon BF(A + ID vs ID) apparié.
- **Limite** : pas de correction de sphéricité : l'a priori est composé-symétrique. Si la sphéricité est fausse, les BF héritent de cette hypothèse. `m.test` traite la question par ailleurs, le bilan doit le signaler. Voie alternative : MANOVA sur les contrastes intra (§6).

### 5.8 Mixte avec pentes aléatoires

- **Modèle** : y = Xβ + Z_id u + Z_{A:id} v + ε, composantes indépendantes dans BayesFactor (`whichRandom = c("id", "A:id")` : fonctionne à l'essai, mais le manuel ne décrit pas cet usage). Corrélations pentes-intercepts : seulement via brms, LKJ [Lewandowski2009].
- **Quelle comparaison ?** Le résultat dépend de la structure aléatoire sous H0 :
  - **Défaut : pentes présentes sous H0 et H1**, BF(A + id + A:id vs id + A:id). Cohérent avec [Barr2013], et [vanDoorn2021] rapporte que l'ignorance des pentes peut gonfler le BF.
  - **Sensibilité : moyenne sur structures** [Heck2021], qui soutient qu'aucune comparaison ne convient toujours.
- **Probabilité annoncée** : P(effet moyen de A ≠ 0 | D, pentes présentes). Ne pas présenter cela comme « A n'a d'effet sur personne » : H0 ici autorise des effets individuels de moyenne nulle.
- **Limites** : plusieurs effets aléatoires croisés ou peu de niveaux : erreur MC à surveiller (§2). Peu de sujets : Stan avec demi-t sur les écarts-types [Gelman2006]. Le débat n'est pas clos [vanDoorn2021, Heck2021].

MANOVA : §6. Non-normalité, hétéroscédasticité, petits effectifs : §7.

### 5.9 Données discrètes

| Donnée | Modèle | A priori | Calcul |
|---|---|---|---|
| Binaire, binomiale | Logistique | Cauchy(0, 2,5) sur prédicteurs mis à l'écart-type 0,5 [Gelman2008] (soit 1,25 sur prédicteurs réduits à 1) | brms + bridge sampling |
| Comptages | Poisson ou négatif binomial si surdispersion | Normale faible, échelle log, à calibrer par prédictif a priori [Schad2021] | Idem |
| Ordinale, peu de modalités | Cotes proportionnelles [McCullagh1980, BurknerVuorre2019] | Effets : même logique ; seuils Dirichlet | Idem |
| Tableau de contingence | Multinomial | Défaut de `contingencyTableBF` | `BayesFactor` |

**Probabilité** : celle du §3 avec un autre terme d'erreur. **Lien avec m.test** : quand le chemin discret de `m.test` traite les comptages comme continus ou en rangs, `m.bay` ne doit pas appliquer le gaussien sans le dire [LiddellKruschke2018].

**Limites** : la Cauchy sur un log-odds n'est pas le même objet que la JZS ; il n'y a pas d'équivalence entre r = 1/2 et 2,5. Séparation complète à petit n : l'a priori décide, il faut l'imprimer. ROPE sur log-odds : ±0,18 est un usage à vérifier, je ne l'ai pas sourcé et je ne le propose donc pas comme défaut. L'utilisateur fournit la valeur.

## 6. MANOVA

`BayesFactor` n'a rien pour cela (vérifié dans [BayesFactorManual]). Options réelles :

| Option | Pour | Contre | Décision |
|---|---|---|---|
| **A. NIW analytique** à g-prior, Jeffreys sur Σ | Exact, rapide, invariant, réduit au défaut univarié | Normalité multivariée, Σ commun, à écrire nous-mêmes | **Défaut** |
| B. BIC multivarié [Wagenmakers2007, Schwarz1978, KassWasserman1995] | Une ligne depuis le Λ de Wilks | A priori implicite peu lisible, biais vers H0 à n modeste | Contrôle de cohérence seulement |
| C. Stan `mvbind` + bridge [Burkner2017, Gronau2017] | Hétéroscédastique, t, mixte | Lent, MC. Sa prise en charge des a priori sur effets standardisés demande un code Stan dédié | Si A ne suffit pas |
| D. BF univariés combinés | Simple | Le produit suppose des DV indépendantes, et compte double une évidence corrélée. **Faux si r ≠ 0** | **Rejeté comme test global** ; affichés DV par DV à titre descriptif |
| E. BFpack [Mulder2021, Mulder2019, MulderGu2022] | Existant ; contraintes d'ordre | Fractionnaire, autre a priori, pas de moyenne de modèles ni de cohérence avec r | Comparaison externe |

### Formule retenue

Y (n × d), lignes N_d(XB, Σ). Y et X centrés, μ plat. Cas général (a priori par terme) :

    B_T | Σ, g_T ~ MN(0, g_T I_{q_T}, Σ),   g_T ~ InvGamma(1/2, r²/2),   p(Σ) ∝ |Σ|^{-(d+1)/2}

Avec V0 = diag(g_T I_{q_T}), A = V0⁻¹ + X'X, S = Y'Y − Y'X A⁻¹ X'Y :

    log m(Y | g) = C + (d/2)(−log|A| − log|V0|) − ((n−1)/2) log|S|

où **C est commun à tous les modèles**, car tous partagent Σ, μ et ν = n − 1. Le modèle nul : log m0 = C − ((n−1)/2) log|Y'Y|. Puis :

    BF_T = ∫ exp(log m(Y | g) − log m0) ∏ InvGamma(g_T) dg

Cas particulier d'un seul terme, en plan équilibré (n_c par cellule, X̃'X̃ = n_c I) :

    BF(g) = (1 + n_c g)^{-q d / 2} · ∏_{i=1}^{s} (1 − (n_c g / (1 + n_c g)) θ_i)^{-(n−1)/2}

où θ_i sont les valeurs propres de E0⁻¹H, avec E0 la SSCP totale corrigée et H la SSCP de l'hypothèse, c'est-à-dire **les mêmes objets que Pillai et Wilks**. Lien direct avec le `car::Manova` de `sys_manova_analysis.R`.

**Pourquoi ce choix** : pour toute combinaison a'Y, a'B | Σ, g a pour variance g·a'Σa : **chaque combinaison linéaire reçoit l'a priori univarié standardisé** avec le même r (partageant g). Cela donne une interprétation et une cohérence uniques. La structure conjuguée suit [Press2005, Rowe2002, Brown1998].

### Contrôles (un jeu de données par cas)

| Contrôle | Résultat |
|---|---|
| Formule NIW propre générale contre Monte Carlo brut (n = 6, d = 2, 200 000 tirages) | −19,42 contre −19,37 |
| Forme fermée contre limite NIW numérique | −1,076013 des deux côtés |
| d = 1 contre forme de [Liang2008] | −0,07531919 des deux côtés |
| d = 1 intégré sur g contre `anovaBF`, r = 0,25 ; 0,5 ; 1 | −0,05792 / −0,33385 / −0,94861 identiques à 10⁻⁶ |
| Invariance à Y → YA | Identique |
| BIC contre analytique, H0 vraie, n = 30, d = 3, r = 1/2 | log BF −8,46 contre −2,41 |

Dernière ligne : **une** simulation, mais le BIC y favorise H0 bien plus que l'analytique. Raison de plus pour ne pas en faire la réponse principale. À quantifier.

### Piège des a priori impropres
1. **μ plat et Jeffreys sur Σ sont permis, car communs à tous les modèles** : leur constante s'élimine. Conditions : même d, même ν = n − 1 partout, et même centrage [ConsonniVeronese2008].
2. **Un a priori impropre sur B_T est interdit** : constante arbitraire, BF arbitraire. Plus l'a priori propre est vague, plus H0 est favorisée [Lindley1957, Bartlett1957].
3. **Une Inverse-Wishart propre « vague » n'est pas anodine** : (ν0, S0) influence le BF, S0 impose une échelle et brise l'invariance à YA si S0 ∝ I. Si l'on préfère un a priori propre, S0 doit être calé empiriquement (avec le soupçon d'un double usage) ou informé, et on affiche l'effet de S0. D'où le choix de Jeffreys. Alternatives : [HuangWand2013], LKJ [Lewandowski2009], non conjuguées, donc voie Stan.
4. **Refus explicite** si n − 1 < d (E0 singulière). Si n − 1 − Σ q_T < d, E1 est singulière, des θ_i valent 1 et l'intégrale sur g peut diverger : refuser aussi, au lieu d'imprimer un nombre.

### Probabilités annoncées
- **Globale** : P(T a un effet sur le vecteur des DV | D).
- **Par DV (descriptive)** : BF univarié pour chaque DV, sans multiplicité ajustée, explicitement pas une décomposition de la globale.
- **ROPE multivarié** : Mahalanobis f²_T = tr(Σ⁻¹ B'B)/q ; seuil non standard. **Aucun défaut proposé**, à fournir ou à omettre. Une convention inventée serait pire.

### Limites

Σ commune par hypothèse, sans équivalent bayésien automatique d'un Box's M ; normalité multivariée ; le BF baisse mécaniquement avec d (pénalité q·d) quand les DV ajoutées sont sans rapport : observé −0,33 à d = 1 et −0,87 à d = 3, r = 0,5, un jeu de données. **r = 1/2 sur d grand n'a pas été calibré.** Facteurs multiples : intégration numérique sur g de dimension = nombre de termes, faisable jusqu'à 3 environ (quadrature), au-delà QMC ou Stan, à évaluer. Contrairement à d = 1, aucune référence multivariée de type `anovaBF` ne permet de valider le défaut : seuls les contrôles ci-dessus le soutiennent. **Cette construction est ma proposition, pas une méthode publiée telle quelle** : sa structure est classique, le réglage (g par terme à r partagé avec l'univarié) n'a pas de source vérifiée, il doit donc être présenté comme propre à KefiR et validé par simulation.

## 7. Branche non normale

### Faut-il un test préliminaire qui choisit la branche ?

**Non**, malgré la cohérence apparente avec `m.test`.

- Choisir par test préliminaire puis annoncer P(T | D, branche choisie) **ignore l'incertitude de sélection**. Ce que [Hoeting1999] reproche à toute sélection. Le p annoncé serait conditionnel à un modèle tenu pour certain, alors qu'il a été choisi sur les mêmes données.
- À petit n, un test de normalité n'a pas de puissance ; à grand n, il détecte l'insignifiant, alors que la conséquence réelle sur le BF peut être petite.
- [Oelrich2020] : les probabilités de modèle sont surconfiantes quand la variabilité est sous-estimée, ce qui est le cas des queues épaisses sous gaussien. Ce risque est réel.

**Mais la moyenne de vraisemblances n'est pas gratuite non plus** : lente, dépend des a priori de chaque branche, et mélange des paramètres d'échelles différentes (δ gaussien n'est pas δ Student). Le compromis proposé :

| Niveau | Déclenchement | Ce que l'on fait |
|---|---|---|
| 1 | Toujours | Gaussien, avec diagnostics (m.bay réutilise ceux de m.test, **à titre d'alerte, pas d'aiguillage**) |
| 2 | Alerte, ou `robust = TRUE` | Ajoute la ligne robuste : les deux p s'affichent, plus « concordant » ou « discordant » |
| 3 | `robust = "bma"` | P(T | D) = Σ_ℓ P(ℓ | D) P(T | D, ℓ), avec bridge sampling. Poids de vraisemblance affichés |

En cas de discordance, le bilan **ne désigne pas de gagnant** : il dit que la réponse dépend de la forme du bruit et que p n'est pas robuste. Ce sera moins confortable que `m.test`, c'est la réalité.

### Les modèles

**Student-t** [Lange1989, Kruschke2013, Burkner2017]. y ~ t_ν(Xβ, σ), même a priori standardisé sur β/σ, ν > 1 faiblement informé (voir [Kruschke2013] pour la forme, à recopier d'après la source avant implémentation). Pour rester à l'échelle JZS, un **code Stan dédié** est nécessaire : `brms` place par défaut ses a priori sur β et non sur β/σ. Calcul : bridge sampling, avec prudence sur l'erreur [Schad2023].

**Hétéroscédastique, sans Stan, à un facteur : premier essai, retiré.** `.bay_bf_oneway_hetero` (mesure plate sur les rapports d'échelle, τ² commun) a été retirée du code : défauts B1 et B2 de `RELECTURE_ADVERSE_2026-11-06.md`. Modèle et chiffres dans `PROJET_METHODO_2026-11-06.md` et `VALIDATION.md` §6.

### Lois d'échantillonnage comparées, voie livrée pour 2 et k groupes (`R/sys_bay_laws.R`, 6-7 novembre)

Remplace `.bay_bf_oneway_hetero`, retirée du code (défauts B1 et B2 ; ses chiffres restent au §6 de `VALIDATION.md` à titre d'histoire). Huit modèles, à poids a priori égaux [Maier2025, à 2 groupes] :

    {normale, t} × {variance commune, variance par groupe} × {pas d'effet, effet}

**Modèle.** σ_i² = c_i / τ̄, avec c_i = 1/(k ϑ_i).

| Élément | Choix | Source |
|---|---|---|
| (μ, τ̄), partagés par tous les modèles | p ∝ 1/τ̄, μ plat | Mesure invariante à droite, emboîtement : [Dablander2024, §2.2, lu sur le PDF arXiv]. C'est ce qui rend le BF entre structures de variance défini |
| Proportions de précision ϑ, variance par groupe | Dirichlet(u, …, u), u = 1,5 ; c_i = 1 sinon | Forme : [Dablander2024]. Valeur : Beta(1,5 ; 1,5) de [Maier2025] à k = 2 (90 % de l'a priori sur un rapport d'écarts-types entre 1/3 et 3, d'après le résumé de leur texte). [Dablander2024] propose u = 1/2 comme défaut et montre que u ≤ 1/2 est nécessaire à la cohérence d'information. Ici u n'est pas testé mais sert de nuisance : u = 1,5 pénalise plus l'explication « un groupe éloigné n'est qu'une grande variance » (voir ci-dessous), et la sensibilité couvre 0,5 et 4,5 |
| Effet, H1 | a_i \| g ~ N(0, g σ_i²), g ~ InvGamma(1/2, r²/2), r = 1/2 | Zellner-Siow [Rouder2012] quand les c_i valent 1. **À k ≥ 3 avec variances inégales, réglage propre à KefiR.** L'échelle par groupe vient de B2 |
| Queues lourdes | y ~ t_ν de variance σ_i², ν − 2 ~ Exp(1) | [Maier2025] (densité exp(−(ν−2)), lue dans le résumé de leur texte) |

À c_i = 1 et loi normale, c'est exactement `anovaBF` et `.bay_bf_oneway_zs`, testé à 1e-5 et 1e-4.

**Calcul, sans compilateur et sans paquet supplémentaire.**
- Normale : conditionnellement à (c, g), tout s'intègre (μ, a, τ̄) ; il reste g par trapèzes en log g et ϑ par échantillonnage d'importance adaptatif. La proposition mêle une part défensive (a priori, variance propre, variance propre plus distance à la médiane pondérée des moyennes, version aplatie) à une t ajustée sur un pilote.
- t : la loi de θ sachant σ s'intègre en g et μ (forme fermée en tête de `.bay_t_prior_theta`), et la densité est directe. Mode, puis t multivariée à la mode, puis t ajustée sur pilote si l'ESS est faible.
- Les deux annoncent ESS et erreur relative. Le premier essai, sans pilote, avait un ESS de 2 à 5 pour un petit groupe très éloigné, avec des estimations fausses de plusieurs unités de log ; c'est corrigé, voir `VALIDATION.md`.

**Probabilité annoncée.** P(H1 | D) = Σ_ℓ π_ℓ π_1 m(D | ℓ, H1) / Σ_{ℓ,H} π_ℓ π_H m(D | ℓ, H). Le BF annoncé est le rapport des deux moyennes, tous modèles confondus : il n'est pas un BF d'une seule loi, et chaque loi garde le sien dans `$laws$log_bf`.

**Ce que le choix règle, et ce qu'il ne règle pas.**
- B2 (groupe dispersé sans rapport) : plancher, la réponse reste du côté de H1 quand S croît. L'échelle moyenne des précisions et le a_i par groupe y contribuent ensemble.
- B1 (petit groupe éloigné de d) : sous H0 hétéroscédastique, m ∝ d^−(n_s + 2u) ; sous H1 à g commun, m ∝ d^−k. Donc BF ∝ d^(n_s + 2u − k). Le signe a changé par rapport à la voie rejetée pour n_s + 2u > k. **Il reste un point faible**, annoncé dans le bilan, si n_s + 2u < k : un très petit groupe parmi beaucoup peut se lire comme une grande variance. Le remède proprement dit (un a priori par groupe) perdrait l'identité avec `anovaBF`. Ce n'est pas absurde sur le fond : avec 3 valeurs, distinguer « décalé » de « très dispersé » repose sur un t à 2 ddl.
- Le mélange ne remplace pas la vérification : il ne connaît que ces quatre lois. Forte asymétrie : non couverte (voir `VALIDATION.md`).

**Sans aiguillage.** Aucun test préalable : les lois sont comparées par leurs probabilités a posteriori. `r` et les cotes a priori restent fixés d'avance : les régler sur les données ne pourrait que gonfler le résultat.

**Rangs** [vanDoorn2020] : pour 2 groupes et apparié, disponible, à utiliser. **Pour k > 2 groupes, je n'ai trouvé aucune méthode publiée vérifiée d'équivalent bayésien de Kruskal-Wallis** (recherche Crossref sans résultat, qui n'exclut pas l'existence d'une méthode).

### Si m.test passerait en Kruskal

Proposition : **cotes proportionnelles** [McCullagh1980, BurknerVuorre2019], avec y traité comme ordonné.

- Pourquoi : n'utilise que l'ordre, donc invariant aux transformations monotones comme un test de rangs. Passe à plusieurs facteurs et mixte.
- **Limites** : ce n'est pas « le Kruskal bayésien », et le bilan ne doit pas employer ce terme. Un seuil par valeur distincte à données continues, donc coûteux, avec a priori sur les seuils à vérifier. Hypothèse PO elle-même testable. Pour Friedman : effet aléatoire sujet. Le lien entre Wilcoxon et test de score PO est de mémoire, à sourcer avant de l'écrire dans une documentation.

### Petits effectifs

**Dire la vérité** : à n < ~10 par groupe, le BF dépend beaucoup de r et de la forme de vraisemblance, et t, hétéroscédastique ou rangs n'ont aucune information pour se distinguer. Le bilan doit afficher la plage de sensibilité ; si elle traverse 1/2, il répond « les données ne tranchent pas », même si le point central est p = 0,8. Le seuil « ~10 » est un repère à valider par simulation.

## 8. A priori sur les hypothèses et sensibilité

**Défaut P(H0) = P(H1) = 1/2.** Argument : ce n'est pas neutre mais c'est lisible, et le BF se lit alors directement. La question « Oui avec probabilité p » **oblige** à le poser : il n'y a pas de p sans π [Kass1995, Morey2016b]. Présenter p sans π serait la faute principale de cet outil. L'utilisateur peut passer `prior_odds`.

Imprimé systématiquement (ici, pour BF = 10 au r central ; ce sont des calculs, pas des résultats) :

    Sensibilité de P(A a un effet | D)
                      cotes H1:H0
      r            1:4     1:1     4:1
      r/2 = 0,25   ...     ...     ...
      r   = 0,5    0,714   0,909   0,976
      2r  = 1,0    ...     ...     ...
    Plage : [min ; max]. Conclusion stable si tous les p sont du même côté de 1/2.

Les neuf cases coûtent peu : BF à 3 r, cotes par simple calcul. Pour Stan, trois exécutions : option `sensitivity = FALSE` si nécessaire, le bilan le mentionne alors.

Forme de la réponse si instable (exemple fictif) : « p = 0,91 (1:1, r = 0,5), mais 0,41 à 0,98 selon les a priori : conclusion dépendante des a priori ». Ne pas choisir un seuil de « preuve » : échelles de [Kass1995, Jeffreys1961] affichables en option, sans valeur normative. [Benjamin2018] relie 0,005 au BF, sans raison d'y adosser le défaut.

**Comportement à attendre** : doubler r baisse le BF quand l'effet est modeste et le monte pour de très gros effets [Lindley1957, Bartlett1957, Liang2008]. **La sensibilité n'est pas un remède au prior mal choisi**, elle le révèle. Note sur [TendeiroKiers2019] : le BF est sensible aux a priori et une « vérité » de H0 ou H1 n'est pas garantie : la grille ne le règle pas.

## 9. Post-hoc sans correction fréquentiste

**Pas de Bonferroni ni de Tukey : la multiplicité passe par l'a priori** [Westfall1997, ScottBerger2006, GopalanBerry1998].

### Option A (défaut, k ≤ 5) : partitions de moyennes

Une hypothèse par partition des k moyennes (Bell : 5, 15, 52 pour k = 3, 4, 5). Pour chacune, BF contre le nul par `lmBF` avec le facteur recodé. Puis :

    P(μ_i = μ_j | D) = Σ_{partitions réunissant i et j} P(partition | D)

- Pour : **cohérent et transitif**, pas de résultat « a = b, b = c, a ≠ c » impossible en probabilité. Inspiration : [GopalanBerry1998, vandenBerghDablander2025]. Cette dernière référence, non lue à ce stade, à consulter avant de finaliser l'a priori.
- A priori proposé : P(toutes égales) = 1/2, le reste uniforme sur les autres partitions. Compatible avec le π du terme global. Si π(T) = 1/2 : P(μ_i = μ_j) a priori = 1/2 + (1/2)·(B_{k-1} − 1)/(B_k − 1), soit **0,625 à k = 3** (2 partitions sur 5 réunissent i et j, dont le nul) et 0,5 + 0,5·4/14 ≈ 0,643 à k = 4. À afficher, car ce n'est pas 1/2.
- Limite : a priori sur les partitions lui-même discutable (uniforme contre processus de Dirichlet). Égalités exactes entre moyennes moins plausibles en observationnel : proposer alors le ROPE deux à deux, et la probabilité de signe.

### Option B (k > 5) : paires à cotes ajustées

Avec m = k(k-1)/2 paires, et l'idée que l'a priori sur la nullité globale fixe P(H0 famille) = π0 (ici 1/2). Sous indépendance a priori : π_ij = π0^{1/m}.

| k | m | P(H0_ij) a priori | Cotes H1:H0 par paire |
|---|---:|---:|---:|
| 3 | 3 | 0,794 | 0,26 |
| 4 | 6 | 0,891 | 0,12 |
| 5 | 10 | 0,933 | 0,07 |

**Cette forme simple est ma dérivation**, dans l'esprit de [Westfall1997], dont je n'ai pas lu les théorèmes (indépendance fausse en pratique pour des paires partageant des groupes). À présenter comme réglage de m.bay. Effet : pénalité croissante, proche d'un Bonferroni.

**Dans le bilan** : pas de lettres de groupes comme si c'était une décision : tableau des paires avec P(égales), P(différence hors ROPE), signe. Autre voie, pour l'estimation et non le test : rétrécissement hiérarchique [Gelman2012], qui répond à une autre question.

**Limite** : des BF par paires sans ajustement ne sont pas protégés du seul fait d'être bayésiens. Sans a priori tenant compte de m, aucune protection [ScottBerger2010].

## 10. Calibration

### Ce que l'on peut garantir

**Sous le modèle et l'a priori effectivement utilisés pour générer les données**, et calcul exact : parmi les analyses annonçant p, une proportion p ont H1 vraie. C'est le théorème de Bayes moyenné sur le prédictif a priori.

Vérifications exécutables :

1. **Moyenne a posteriori = a priori** (DAP) : simuler M, θ, y sous l'a priori, calculer P(M|y). La moyenne doit égaler π [Schad2023]. Schad et al. le précisent : la condition concerne la moyenne sur simulations, **pas un jeu de données unique**.
2. **Courbe de calibration** : p annoncée contre fréquence de H1.
3. **SBC** pour les paramètres [Talts2018, Cook2006, Modrak2025], et R-hat, ESS [Vehtari2021].
4. **Erreur BF** : même analyse répétée, écart-type de p ; signal immédiat si l'erreur MC dépasse un seuil fixé (ex. 5 %).

Ce sont **des tests de calcul**, pas des tests de vérité du modèle. Ils doivent couvrir chaque voie du §5, MANOVA comprise.

### Ce que l'on ne peut pas garantir

- **Taux de faux positifs fréquentiste.** Pour θ = 0 fixe, P(p > 0,95) n'est pas 5 % : il dépend de n, de r, du plan. Attendu, faible à grand n, puisque le BF est consistant sous H0.
- **Calibration sous a priori faux** : si les vrais effets ne ressemblent pas à la Cauchy, p peut être mal calibrée. [Oelrich2020] est utile sur la dispersion du BF sous spécification erronée.
- **Sens à un seul jeu de données** : p = 0,91 ne signifie pas « 91 % des études de ce type » sans admettre que π et r décrivent le domaine.
- **Arrêt optionnel** : [Rouder2014] le défend pour le BF, le débat suppose un a priori exact. Ne pas promettre la même propriété sous sensibilité large.

Le débat sur la calibration fréquentiste du BF est ouvert [Hoijtink2016, Morey2016c] : m.bay ne recalibre pas le BF, il **publie** ses opérateurs. Position : la p annoncée reste bayésienne, et une table simulée la complète.

### Table à publier (vignette)

Plans × n × d vrai (0 ; 0,2 ; 0,5 ; 0,8) × distribution (gauss, t3, asymétrique, hétéroscédastique), méthode inspirée de [Schonbrodt2018] :

| Sortie | Pour d = 0 | Pour d > 0 |
|---|---|---|
| P(p > 0,95) | « Faux positifs » | Puissance analogue |
| P(p < 0,05) | Soutien à H0 | Fausses nullités |
| P(0,25 < p < 0,75) | Indécis | Indécis |

Comparer avec `m.test` aux mêmes scénarios.

**Phrase à écrire dans la documentation** : « m.bay ne contrôle pas un taux d'erreur à long terme. Sa probabilité est conditionnelle aux a priori affichés et au modèle. La table de simulation fournie indique ses opérateurs fréquentistes pour les plans testés. »

## 11. Risques et points à trancher
1. **« Oui avec probabilité p » sera lu comme une vérité.** Sans π et modèle, trompeur. Imposer le rappel des a priori sur la même ligne.
2. **Simulation** (§10) pas faite ici : aucune affirmation de performance dans ce document.
3. **MANOVA §6** : algèbre contrôlée, calibration pas. Réglage propre à KefiR.
4. **Lectures manquantes** avant citation : [vandenBergh2020], [Hinne2020], [Westfall1997], [vandenBerghDablander2025], et tout le niveau « métadonnées » de la bibliographie.
5. **Choix laissés à Antoine** : seuil de 3 facteurs pour Beta-binomial ; `robust` niveau 2 par défaut ou à la demande ; partitions à k ≤ 5 ; acceptation d'un p non cohérent entre termes (matched) contre une cohérence moins lisible.
