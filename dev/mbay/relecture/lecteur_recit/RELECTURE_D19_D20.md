# Relecture D19-D20 : plusieurs facteurs (cas 1 à 4), par un lecteur débutant

Lecteur joué : étudiant de BUT génie biologique (bac+1/2). J'ai vu l'ANOVA, le khi-deux, les p-values et le seuil de 5 % en cours, mais jamais de statistique bayésienne. J'ai lu uniquement `recits_d19_d20.txt`.

---

## Cas 1 : deux facteurs (engrais, sol), variances très différentes, aucun effet

### 1. Ce que j'ai compris

On a testé si le rendement dépend de l'engrais, du sol ou de leur combinaison, avec des groupes très déséquilibrés (5 à 40 par case). Les moyennes sont toutes autour de 20 et la méthode ne tranche pas : 55,2 % pour « aucune différence », 44,8 % pour « une différence », et selon le point de départ ça peut s'inverser. Pour l'interaction, en revanche, on peut dire « non » (5,6 %).

Ce dont je ne suis pas sûr : ce qui a été fait avec les variances inégales, et si ça change quelque chose au résultat.

### 2. Ce qui m'a arrêté

- **« Variances par case : inégales probables (> 99,9 %) »** puis **« Les deux sont gardées : chaque réponse pèse la variance commune et une variance par case selon leur probabilité »**. Si les variances sont inégales à plus de 99,9 %, la « variance commune » pèse donc presque 0 ? Alors pourquoi dire « les deux sont gardées » ? J'ai lu la même phrase mot pour mot dans les cas 2 et 3, où les variances sont au contraire égales probables. J'en conclus que la phrase est automatique et je ne sais pas ce qui a vraiment été fait ici. En cours, avec des variances inégales, on nous dit de « changer de test » (Welch) : est-ce que c'est l'équivalent ? Rien ne me le dit.
- **Le titre « normale, variance commune ou par case, moyennées »** : « moyennées » quoi ? Les variances ? Les réponses ? Je comprends « on fait la moyenne de deux calculs », mais je ne sais pas si c'est ça.
- **« Écart-type au sein des cases : 2,92 »** (bilan). À l'étape 3, on m'a dit que les écarts-types vont de 0,82 à 7,6 et qu'ils sont différents. Donner un seul chiffre, 2,92, me semble contredire l'étape 3. Je ne le recopierais pas, mais si je le recopiais ce serait faux.
- **« aucune différence du tout : 33,3 % »** (étape 4) alors que chaque facteur part à 50 %. Je ne vois pas d'où sort 33,3 %. La note « ces chances ne s'additionnent pas » m'aide un peu, mais je ne saurais pas l'expliquer.
- **« facteur de Bayes 0,446 »** pour engrais. L'échelle donnée (« 1 à 3 faible… ») ne parle que des valeurs au-dessus de 1. Je devine qu'en dessous de 1 ça penche vers « non », mais il faut que je calcule 1/0,446 = 2,2 moi-même.
- **« engrais non tranché (30,8 %) »** mais **« engrais:sol : non, soutien modéré (5,6 %) »**. Avec les seuils donnés plus loin (cas 2 : « semblables jusqu'à 25 % »), je comprends la règle, mais dans ce cas 1 elle n'est écrite nulle part.
- **« La réponse globale peut s'inverser : le bilan ne l'affirmera pas. »** Bien compris, et utile.

### 3. Ce que je recopierais

La phrase « Les données ne permettent pas de conclure (probabilité d'une différence 44,8 %, facteur de Bayes 2,47 en faveur de l'absence de différence ; …) ». Je n'ajouterais pas l'écart-type 2,92. Je ne dirais rien sur les variances, parce que je ne saurais pas quoi dire.

---

## Cas 2 : deux facteurs (variété, serre), effet de la variété seule

### 1. Ce que j'ai compris

La hauteur dépend de la variété (99,4 %), pas clairement de la serre (non tranché, 29,9 %), et il n'y a pas d'interaction (6,2 %). V1 est plus petite que V2 (environ 3,1 de moins) et que V3 (2,5 de moins), et V2 et V3 sont semblables. Le classement le plus probable est V1 < V3 = V2.

### 2. Ce qui m'a arrêté

- **« ses niveaux sont comparés, à serre égal »**. Je devine : on compare les variétés dans une même serre, pour ne pas mélanger l'effet serre. Mais je ne sais pas si c'est « dans la serre S1 », « dans la serre S2 » ou « en moyenne sur les deux ». En ANOVA on parlait d'effet principal ; est-ce la même chose ?
- **« et dans aucune interaction probable »** : je comprends après coup (si l'interaction avait été probable, on n'aurait pas comparé les variétés seules), mais seulement grâce au cas 3. Ici, la phrase ne dit pas pourquoi c'est important.
- **« Regroupement en tête : V1 < V3 = V2 (75,3 %) »** puis juste dessous **« Si une différence existe, regroupements en tête : V1 < V3 = V2 (75,3 %) ; … »**. Même chiffre deux fois, une fois sans « si une différence existe », une fois avec. Est-ce que 75,3 % est la probabilité tout court ou « sachant qu'il y a une différence » ? Je ne sais pas lequel écrire.
- **« soit 1,13 écart-type résiduel (autres termes pris en compte) »**. « Résiduel » et « autres termes pris en compte » : je ne sais pas ce que ça change. Je comprends seulement « c'est un écart d'environ un écart-type ».
- **« Bilan - hauteur diffère selon au moins un facteur : 99,4 % (soutien très fort) »** avec un facteur de Bayes de 89,8. Or la note dit « 20 à 150 fort ». Donc 89,8 devrait être « fort », pas « très fort ». J'ai relu trois fois. Dans le cas 3, 62,2 est bien dit « fort ».
- **« Différents : V1/V2 (très fort) ; V1/V3 (fort) »** : les mots sans chiffres. À l'étape 6 j'ai 99,7 % et 99,3 % ; presque pareil, et pourtant « très fort » contre « fort ». Je suppose que c'est le facteur de Bayes de chaque paire, mais il n'est jamais donné.
- **« l'écart estimé est ramené vers 0 par l'a priori »** : je comprends l'idée (l'écart estimé 3,1 est un peu plus petit que l'écart observé 34,3 - 30,8 = 3,5). Utile.

### 3. Ce que je recopierais

La phrase de compte rendu, à laquelle j'ajouterais à la main : « V1 plus petite que V2 (écart 3,1 [1,5 ; 4,6]) et que V3 (2,5 [0,9 ; 4,0]) ; V2 et V3 semblables ». La phrase proposée ne contient pas les paires, alors que c'est ce que mon enseignant voudra lire.

---

## Cas 3 : deux facteurs (traitement, sexe), interaction

### 1. Ce que j'ai compris

Le traitement actif augmente le score chez les hommes (environ +8 à +9 points par rapport au placebo), mais pas chez les femmes. Les trois autres cases (F actif, F placebo, H placebo) sont semblables. Comme l'effet du traitement dépend du sexe, la méthode compare les 4 cases entre elles au lieu de comparer « actif contre placebo » en général. Ça, c'est clair et j'aurais aimé le voir aussi bien dit en ANOVA.

### 2. Ce qui m'a arrêté

- **« traitement (seul ou en interaction) : oui, soutien fort (98,7 %) »** et **« sexe (seul ou en interaction) : oui, soutien fort (99,0 %) »**. Si je lis vite, j'écris « le sexe a un effet sur le score ». Or d'après l'étape 6, les femmes et les hommes sous placebo sont semblables : le sexe ne compte que sous traitement actif. « (seul ou en interaction) » est là, mais il disparaît dans la phrase de compte rendu : « traitement oui (98,7 %), sexe oui (99,0 %), traitement:sexe oui ». Recopiée telle quelle, cette phrase dit trois effets, alors qu'il y en a un seul (le traitement chez H).
- **Six paires, dont des paires bizarres** : « F x actif/H x placebo », « H x actif/F x placebo ». Comparer une femme traitée à un homme sous placebo, ça ne répond à aucune question que je me pose. Ce que je veux : l'effet du traitement chez F, l'effet chez H, et peut-être F contre H sous chaque traitement. J'ai dû trier moi-même.
- **L'ordre des noms change** : l'étape 6 annonce « cases traitement x sexe » mais écrit « F x actif » (sexe d'abord) ; le bilan écrit « actif x F » (traitement d'abord). J'ai cru un moment que c'étaient des cases différentes.
- **« H x actif - F x actif = 9,4 [5,3 ; 13,2] »** : la ligne mélange « / », « - » et « x » ; c'est dur à lire. Le sens de la soustraction change d'une ligne à l'autre (parfois positif, parfois négatif), je dois relire à chaque fois.
- **« soit 1,76 écart-type »** : ici « écart-type » tout court, dans le cas 2 « écart-type résiduel (autres termes pris en compte) ». Est-ce le même ?
- **« Regroupement en tête : F x actif = H x placebo = F x placebo < H x actif (71,5 %) »** : même question qu'au cas 2 sur le « si une différence existe ».
- Les moyennes observées (actif x H = 59,6 ; placebo x H = 50,4, soit 9,2) et l'écart estimé (8,1) : j'ai compris grâce à la note sur le « ramené vers 0 ».

### 3. Ce que je recopierais

« L'effet du traitement dépend du sexe (interaction, 97,8 %) : chez les hommes, le score sous traitement actif est plus élevé que sous placebo (écart 8,1 [4,1 ; 11,9]) ; chez les femmes, pas de différence (semblables). » Je ne recopierais pas « sexe oui (99,0 %) » seul : je ne saurais pas le défendre à l'oral.

---

## Cas 4 : plan en tranches (groupe entre sujets, temps intra-sujet)

### 1. Ce que j'ai compris

Le poids augmente avec le temps (J0 < J7 < J14, presque sûr). Pour le groupe (traité contre témoin), la méthode ne tranche pas : 73,6 %, et selon le point de départ ça va de 43,7 % à 85,0 %. L'interaction n'est pas tranchée non plus. Une valeur très éloignée a été gardée, donc il faut lire avec réserve.

Ce dont je ne suis pas sûr : est-ce que le groupe traité pèse plus que le témoin ? Les moyennes disent 64,8 contre 60,1 (4,7 de plus), l'étape 6 dit 1,5, et le bilan dit « non tranché ».

### 2. Ce qui m'a arrêté

- **Étape 5 : « Prudence ajoutée : groupe sera dit non tranché. »** puis **étape 6 : « groupe probable (73,6 %) et dans aucune interaction probable : ses niveaux sont comparés »**, puis **bilan : « groupe … : non tranché (73,6 %) »**. Le même 73,6 % est appelé « probable » à l'étape 6 et « non tranché » à l'étape 7. Je ne sais pas ce que je dois croire. Si c'est non tranché, pourquoi comparer ses niveaux ?
- **« Écart estimé traite - temoin = 1,5 [-0,6 ; 6,1] »** alors que les moyennes observées donnent 64,8 - 60,1 = 4,7. Pour les autres cas, l'écart « ramené vers 0 » restait proche de l'observé (3,5 → 3,1). Ici il est divisé par trois. Et l'intervalle n'est pas centré (de -0,6 à 6,1 autour de 1,5). Je ne comprends pas, et je ne sais pas lequel des deux chiffres mettre dans mon rapport.
- **Pas de « différents à X % » pour traite/temoin**, alors que pour le temps chaque paire en a un. Je suppose que c'est 73,6 %, mais ce n'est pas écrit.
- **« Précision numérique : deux calculs avec des graines différentes diffèrent jusqu'à 5,1 points de probabilité. »** « Graines » : aucune idée de ce que c'est ici (pas des graines de plantes, je suppose). Et contrairement aux autres cas, il n'est pas dit si c'est grave ou « négligeable ». 5 points sur 73,6 %, ça me paraît beaucoup.
- **« Forme du bruit - normale, seule. »** « Seule » : seule par rapport à quoi ? Je comprends en lisant la note plus bas (pas de variances inégales avec des sujets), mais le titre seul est obscur.
- **« m.test() contrôle la sphéricité »**. « Sphéricité » : jamais vu. Et on me renvoie vers une autre fonction : dois-je la lancer ? Rien ne le dit.
- **« chaque sujet garde son niveau propre »** : j'ai compris (chaque animal a son poids de base), c'est bien dit.
- **« à temps égal »** pour le groupe : je comprends « en comparant les groupes au même jour », mais est-ce un jour précis ou la moyenne des trois jours ?
- **L'ordre des jours** : « J0 = 61,8 ; J14 = 63,1 ; J7 = 62,5 » et « J14/J7 : J7 - J14 ». L'ordre alphabétique met J14 avant J7, ça casse la lecture d'un suivi dans le temps.
- **« J14/J7 : différents à 89,7 % »** dit « modéré » au bilan, et **« 2,41 écart-type résiduel au sein des sujets »** : l'écart réel (1,2) est petit, mais 2,41 écarts-types paraît énorme. Je comprends que c'est parce que l'écart-type au sein des sujets est petit (0,44), mais c'est piégeux.
- **La note sur groupe:temps** (« les données augmentent les chances … mais ce terme partait à 16,7 %, d'où 37,8 % seulement ») : très claire, c'est la meilleure explication du point de départ de tout le texte.

### 3. Ce que je recopierais

« Le poids augmente avec le temps (J0 < J7 < J14, probabilité 88,8 % pour ce classement ; J14 - J0 = 1,2 [0,8 ; 1,5]). L'effet du groupe n'est pas tranché (73,6 %, sensible à l'a priori). » Je ne recopierais pas l'écart traite - temoin, car je ne comprends pas pourquoi il est si loin des moyennes.

---

## Les 8 remarques les plus importantes

1. **Cas 4 : « probable » à l'étape 6, « non tranché » au bilan, pour le même 73,6 %.** L'étape 5 annonce la prudence, l'étape 6 l'oublie. Proposition, étape 6 : « --> groupe non tranché (73,6 %, change de côté selon l'a priori) : ses niveaux sont comparés à titre indicatif, à temps égal. » Ou ne pas le comparer du tout, comme dans le cas 1.

2. **Cas 3 : la phrase de compte rendu dit « sexe oui (99,0 %) », alors que le sexe n'agit que sous traitement actif.** « (seul ou en interaction) » disparaît dans la phrase à recopier. Proposition : « traitement:sexe oui (97,8 %) : le traitement n'agit que chez H (H actif - H placebo = 8,1 [4,1 ; 11,9]) ; traitement et sexe ne s'interprètent pas seuls » ; et au bilan : « sexe (seul ou en interaction) : oui, mais à travers l'interaction ».

3. **Cas 1 : « inégales probables (> 99,9 %) » suivi de « Les deux sont gardées ».** La phrase est identique quand les variances sont égales à 90 % ; je ne sais pas ce qui a été fait. Proposition : « --> Les réponses tiennent compte de variances différentes par case (la variance commune ne pèse presque plus, < 0,1 %). Comparable à un test de Welch. » Et au bilan, remplacer « Écart-type au sein des cases : 2,92 » par « Écarts-types par case : 0,82 à 7,6 (inégaux) ».

4. **Cas 4 : écart estimé 1,5 contre 4,7 observé, sans explication.** La note générale « ramené vers 0 » ne suffit pas quand l'écart est divisé par trois. Proposition : « Note : écart observé 4,7, estimé 1,5 : avec 8 sujets par groupe et des sujets très différents entre eux, l'a priori ramène fortement l'écart vers 0. » (si c'est bien la raison).

5. **« à serre égal », « à temps égal » : non expliqué.** Proposition, la première fois : « ses niveaux sont comparés à serre égal (dans une même serre, en moyenne sur S1 et S2 ; possible car l'écart entre variétés ne dépend pas de la serre) ».

6. **Cas 3 : six paires dont des diagonales, et l'ordre des noms qui change.** Proposition : mettre en tête les comparaisons qui ont un sens : « Effet de traitement chez F : … ; chez H : … ; effet de sexe sous placebo : … ; sous actif : … », puis les deux paires restantes en dernier. Toujours écrire les cases dans l'ordre annoncé (« actif x F » partout, puisque le titre dit « traitement x sexe »).

7. **Cas 2 : « soutien très fort » pour un facteur de Bayes de 89,8, alors que l'échelle dit 20 à 150 = fort.** Le cas 3 (62,2 → « fort ») est cohérent, le cas 2 non. Proposition : « hauteur diffère selon au moins un facteur : 99,4 % (soutien fort) », ou dire dans la note sur quoi porte le mot (probabilité ou facteur de Bayes).

8. **« Regroupement en tête » répété deux fois, une fois sans « si une différence existe ».** Je ne sais pas si 75,3 % est conditionnel. Proposition : une seule ligne : « ==> Si une différence existe, regroupements en tête : V1 < V3 = V2 (75,3 %) ; V1 < V3 < V2 (24,2 %) ; … » ; et, plus largement, donner le facteur de Bayes (ou la probabilité) derrière les mots « très fort / fort / modéré » des paires au bilan, et rappeler qu'un facteur de Bayes sous 1 se lit à l'envers (0,446 = 2,2 fois en faveur de « non »).

Remarques secondaires : « graines » et « sphéricité » (cas 4) sans explication ; « normale, seule » obscur ; jours rangés J0, J14, J7 au lieu de J0, J7, J14 ; d'où vient 33,3 % à l'étape 4 ; la phrase de compte rendu du cas 2 ne contient pas les paires, qui sont pourtant la réponse utile.

---

## Suites données (session m.bay, 10 novembre)

Récits régénérés après correction dans `recits_d19_d20.txt`.
- **Cas 4, contradiction étape 6 / étape 5.** Un terme qui change de côté selon l'a priori est désormais dit « non tranché (…, étape 5) mais en tête : ses niveaux sont comparés à titre indicatif ». Pour un facteur entre sujets, une ligne explique que l'écart est fortement ramené vers 0 (a priori à l'échelle de la variation intra-sujet) et donne l'écart observé des moyennes (4,7).
- **Cas 3, phrase de compte rendu.** Un effet principal contenu dans une interaction probable s'écrit « oui, seul ou en interaction, (…) ».
- **Cas 3, ordre des noms de cases.** Les cases suivent l'ordre de la formule (« actif x F »), à l'étape 6 comme au bilan.
- **Étape 3, variances inégales.** La phrase de décision dépend maintenant du résultat : « calculée surtout avec une variance par case : une case très dispersée ne fait plus paraître les écarts plus grands qu'ils ne sont », ou « surtout avec une variance commune, comme d'habitude ».
- **« à serre égal ».** Une note l'explique : « chaque écart est mesuré dans un même niveau de serre, puis mis en commun ».

**Non repris ici.**
- L'échelle « soutien très fort » pour un facteur de Bayes de 89,8 relève de la mise en forme commune (D18).
- Les 6 paires de cases, dont des comparaisons sans intérêt (F actif contre H placebo) : on les garde toutes, car le regroupement les juge ensemble. Un filtre des paires « à un facteur près » est une piste.
- L'ordre J0, J14, J7 vient des niveaux du facteur dans les données.
- « graines » et « sphéricité ».
