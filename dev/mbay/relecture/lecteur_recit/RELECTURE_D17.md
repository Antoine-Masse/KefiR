# Relecture D17 : réponses discrètes (cas J, K, L), par un lecteur débutant

Lecteur joué : étudiant de BUT génie biologique (bac+1/2). J'ai vu le khi-deux, les p-values et le seuil de 5 % en cours, mais jamais de statistique bayésienne. J'ai lu uniquement `recits_contingence.txt`.

---

## Cas J : germination (oui/non) selon 3 traitements

### 1. Ce que j'ai compris

Les graines traitées au chaud germent plus dans l'échantillon (80 % contre 60 %), mais ce n'est pas assez pour conclure : la méthode dit « indécis ». Elle penche plutôt pour « pas de différence » (68,4 %), sans que ce soit sûr, puisque selon le point de départ ça peut passer sous 50 %.

Ce dont je ne suis pas sûr : est-ce que j'ai le droit d'écrire « les traitements n'ont pas d'effet » ? Les lignes « semblables » me poussent à le faire, le verdict « indécis » me l'interdit.

### 2. Ce qui m'a arrêté

- **« temoin 12/18 ; froid 12/18 ; chaud 6/24 (non/oui) »** (étape 2). J'ai d'abord lu « 12 sur 18 », comme une fraction. Je ne comprenais pas d'où sortait 18 si on avait 30 graines. C'est seulement à la fin, avec les pourcentages, que j'ai compris que c'était 12 « non » et 18 « oui ». En plus l'ordre est « non » puis « oui », alors que ce qui m'intéresse, c'est le « oui ».
- **« Facteur de Bayes, différence contre absence : 0,463 »**. Je ne sais pas si 0,463 c'est beaucoup ou peu. J'ai compris seulement à l'étape 7 avec « 2,16 fois plus probables si les proportions sont les mêmes ». Rien ne dit que 2,16 = 1/0,463, ni qu'en dessous de 1 le facteur penche pour « pas de différence ».
- **« Il n'y a pas de forme du bruit à comparer. »** Je ne sais pas ce qu'est la « forme du bruit ». La phrase n'a l'air de servir à rien pour moi.
- **« Dirichlet uniforme, Gunel et Dickey 1974 »**. Nom inconnu. La partie en clair (« toutes les proportions sont également possibles d'avance ») suffit ; la référence me fait peur sans m'aider.
- **« calcul en forme close, sans tirage ; aucune erreur de simulation »**. Je n'ai jamais entendu parler de tirage ni de simulation dans un test. Je comprends seulement que « c'est précis ».
- **« concentration de l'a priori de 0,5 à 2, croisée avec des chances de départ d'une différence de 1 contre 2 à 2 contre 1 »**. C'est la phrase la plus dure du texte. « Concentration » ne veut rien dire pour moi. « 1 contre 2 » : je ne sais pas si c'est 33 % ou 50 %.
- **« La réponse en tête va alors de 42,2 % à 89,4 % »**. « Réponse en tête » : laquelle ? J'ai deviné que c'était « pas de différence ». Et 42,2 %, ça veut dire que l'autre réponse passe à 57,8 % ; il faut que je fasse le calcul moi-même pour voir que l'ordre s'inverse.
- **Deux chiffres pour la même idée** : « temoin = froid = chaud (58,9 %) » à l'étape 6, et « ne diffèrent pas : 68,4 % » à l'étape 7. Pour moi c'est la même hypothèse (tout le monde pareil). Lequel je recopie ? Je pense que l'un des deux est faux.
- **« temoin < froid < chaud (6,8 %) »**. Témoin et froid ont exactement les mêmes résultats (60 % de oui chacun). Pourquoi proposer que froid soit au-dessus de témoin ?
- **« log du rapport de cotes »**. Jamais vu. « Rapport de cotes », je ne sais pas ce que c'est, et le log encore moins. Ensuite, à l'étape 7, ça s'appelle **« delta »**, sans dire que c'est la même chose.
- **« écart estimé chaud - temoin = 0,19 »** alors que les proportions sont données en pourcentages (60 % et 80 %). J'ai mis du temps à voir que 0,19 = 19 points. Et pourquoi 0,19 et pas 0,20 ? (je suppose que c'est la « médiane a posteriori », mot que je ne connais pas).
- **« temoin / chaud : semblables, soutien faible (différents à 35,1 %) »**. « Semblables » + « soutien faible » : ils sont semblables ou pas ? Et « temoin / froid : semblables, soutien modéré » à 81 % : alors que le verdict global est « indécis », ici on me dit « semblables ». J'ai peur d'écrire « le froid n'a pas d'effet ».
- **« Toutes les paires viennent d'une seule loi sur les regroupements : aucune correction paire par paire n'est nécessaire »**. En cours on m'a dit qu'il fallait corriger quand on fait plusieurs comparaisons (Bonferroni). Ici on me dit que non, mais je ne comprends pas pourquoi.
- **« (la valeur s'y trouve avec 95 % de probabilité ; ce n'est pas un intervalle de confiance) »**. Je ne connais pas la différence. Pour moi un intervalle de confiance, c'est justement « la valeur est dedans à 95 % ».

### 3. Ce qui me manque pour le compte rendu

- Une **phrase toute faite** à recopier, du genre : « Les données ne permettent pas de conclure sur un effet du traitement sur la germination (probabilité d'une différence : 32 %). »
- Ce qu'il faut dire de l'écart du chaud (+20 points) : j'ai le droit de le mentionner comme une tendance, ou non ?
- Le **nom de la méthode** et comment la citer (« facteur de Bayes de contingence » ? le nom de la fonction ?).
- La comparaison avec ce que j'aurais trouvé au **khi-deux**, puisque c'est ce que le prof attend peut-être : on me dit « pendant du test du khi-deux » mais pas ce qu'il aurait donné.
- Ce que je devrais faire pour trancher (combien de graines en plus ?).

---

## Cas K : réussite 0/1 selon 2 groupes de TD

### 1. Ce que j'ai compris

Le TD2 réussit un peu mieux dans l'échantillon (55,6 % contre 38,9 %, écart de 16 points), mais on ne peut pas conclure : verdict « indécis », 61,2 % pour « pas de différence », et ça peut s'inverser.

### 2. Ce qui m'a arrêté

- **« TD1 11/7 ; TD2 8/10 (0/1) »** : même souci de fraction que pour J ; et ici « 0/1 » ressemble encore plus à une division.
- **« TD1 : 0 61,1 %, 1 38,9 % »**. Les chiffres se collent : « 0 61,1 % » se lit comme un seul nombre. Et je dois me souvenir que 0 = échec et 1 = réussite (rien ne le dit, je le devine).
- **« différence des proportions de 1 »** : « proportion de 1 » est bizarre à lire ; j'aurais aimé « proportion de réussite ».
- **« Les proportions de reussite »** : on parle « des proportions » au pluriel alors qu'il n'y a qu'une chose qui compte (le taux de réussite). Le mot « reussite » sans accent fait négligé (je suppose que c'est le nom de la colonne).
- **« Facteur de Bayes ... : 0,633 »** puis **« 1,58 fois plus probables »** : même remarque que pour J, le lien n'est pas dit.
- **« soit 0,66 [-0,55 ; 2] en log du rapport de cotes »** : même remarque que J.
- Les étapes 2 à 5 sont **mot pour mot les mêmes** que pour J. Je les ai sautées à la deuxième lecture ; ça va, mais du coup j'ai failli rater les chiffres qui changent.

### 3. Ce qui me manque pour le compte rendu

- La phrase de conclusion prête à l'emploi.
- Est-ce que l'intervalle [-0,14 ; 0,45] qui contient 0 veut dire « pas significatif » comme en cours ? On me dit qu'il n'y a pas de seuil, mais je ne sais pas comment utiliser l'intervalle.
- Une indication de ce qu'il faudrait comme effectif pour voir un écart de 16 points.

---

## Cas L : couleur des fleurs (3 couleurs) selon 2 parcelles

### 1. Ce que j'ai compris

Les deux parcelles n'ont pas la même répartition des couleurs (99,1 %, soutien fort, et ça ne dépend pas du point de départ). P1 a surtout des fleurs blanches, P2 surtout des rouges. L'écart est de 0,37, mais je ne sais pas vraiment ce que mesure ce 0,37.

### 2. Ce qui m'a arrêté

- **« Cela ne dit ni qu'elle est importante, ni quels groupes diffèrent, ni pourquoi. »** Il n'y a que deux parcelles : si une différence existe, c'est forcément entre P1 et P2. La phrase m'a fait douter de ce que j'avais compris.
- **« Écart estimé P2 - P1 = 0,37 [0,19 ; 0,54] »** puis **« écart total des proportions sur les modalités (moitié de la somme des écarts, de 0 à 1) »**. Le « P2 - P1 » me fait croire que P2 a « plus » de quelque chose, mais de quoi ? Il y a trois couleurs. Je n'ai pas compris ce qu'est cette « moitié de la somme des écarts », ni si 0,37 c'est beaucoup.
- **Il manque le « quelle couleur »** : le texte ne dit nulle part que la différence vient surtout du blanc (57,5 % contre 20 %) et du rouge (17,5 % contre 50 %). Je l'ai trouvé moi-même dans les proportions observées, tout en bas. En khi-deux on regarde les contributions ou les résidus ; ici rien.
- **« Facteur de Bayes ... : 106 »** : là j'ai compris grâce à « 106 fois plus probables si les proportions diffèrent ». Mais je ne sais pas si 106 est « fort » dans une échelle connue.
- **« N = 80 »** dans le titre de l'étape 1 sans dire 40 par parcelle (je l'apprends à la toute fin).
- **« P1 23/10/7 »** : là, trois nombres, je n'ai pas confondu avec une fraction, mais j'ai dû remonter à « (blanc/rose/rouge) ».

### 3. Ce qui me manque pour le compte rendu

- Une phrase du type : « La répartition des couleurs diffère entre les parcelles (probabilité 99 %) : P1 compte plus de fleurs blanches, P2 plus de rouges. »
- Quelles couleurs « portent » la différence, avec un écart et un intervalle par couleur (par exemple : rouge, P2 - P1 = +32 points [..; ..]).
- Une façon de dire si 0,37 est un petit ou un grand écart.

---

## Mes 8 remarques classées (de la plus gênante à la moins gênante)

1. **Les effectifs « 12/18 », « 11/7 », « 0/1 » se lisent comme des fractions.**
   Reformulation : « temoin : 18 oui, 12 non (sur 30) ; froid : 18 oui, 12 non ; chaud : 24 oui, 6 non. » Mettre la modalité d'intérêt (oui, réussite) en premier.

2. **« Semblables » dans les paires contredit le verdict « indécis » (cas J).** Je risque d'écrire « le froid n'a pas d'effet ».
   Reformulation : « temoin / froid : plutôt semblables, sans certitude (différents à 18,8 %) » et « temoin / chaud : non tranché (différents à 35,1 %) ». Réserver « semblables » aux cas où le soutien est au moins modéré, et ne jamais l'écrire avec « soutien faible ».

3. **Deux probabilités différentes pour « aucune différence » (68,4 % et 58,9 %), cas J.**
   Reformulation : supprimer le pourcentage du regroupement à l'étape 7, ou ajouter : « (58,9 % ici, car le calcul répartit les chances de départ entre 5 regroupements au lieu de 2 ; pour la conclusion, retenez 68,4 %) ».

4. **Il manque une phrase de conclusion prête pour un compte rendu.**
   Reformulation, à la fin de l'étape 7 : « Pour un compte rendu : Les données ne permettent pas de conclure à un effet du traitement sur la germination (probabilité d'une différence : 32 % ; facteur de Bayes 0,46 ; test bayésien d'indépendance, Gunel et Dickey 1974). »

5. **Cas L : on ne dit pas quelles couleurs font la différence, et « ni quels groupes diffèrent » n'a pas de sens à 2 groupes.**
   Reformulation : « Verdict : la répartition des couleurs diffère très probablement entre P1 et P2. Cela ne dit pas pourquoi. » puis « Couleurs qui changent le plus : rouge (17,5 % -> 50 %), blanc (57,5 % -> 20 %). »

6. **La phrase de sensibilité est illisible (« concentration de l'a priori », « 1 contre 2 », « réponse en tête »).**
   Reformulation : « En changeant le point de départ (chance d'une différence de 33 % à 67 %, a priori plus ou moins resserré), la probabilité de "pas de différence" va de 42 % à 89 % : elle peut passer sous 50 %, donc l'ordre peut s'inverser. »

7. **« log du rapport de cotes », puis « delta » : mots inconnus, et le nom change entre l'étape 6 et l'étape 7.**
   Reformulation : garder seulement l'écart en points de pourcentage dans la conclusion (« chaud - temoin = +19 points [-2 ; +40] »), et renvoyer le log du rapport de cotes à une ligne à part, avec le même nom aux deux endroits : « (pour un article : log du rapport de cotes = 0,91 [-0,19 ; 2]) ».

8. **Le facteur de Bayes seul (0,463) ne se comprend pas ; le lien avec « 2,16 fois » n'est pas dit.**
   Reformulation à l'étape 5 : « Facteur de Bayes, différence contre absence : 0,463 (moins de 1 : les données penchent vers l'absence de différence ; 1/0,463 = 2,16 fois plus probables sans différence). »

Remarques mineures en plus : « Il n'y a pas de forme du bruit à comparer » (à supprimer pour ce cas) ; « forme close, sans tirage » (dire « calcul exact ») ; le mélange 0,19 / 80 % (choisir des points de pourcentage partout) ; « temoin < froid < chaud » proposé alors que témoin et froid sont identiques ; pourquoi pas de correction de Bonferroni (une demi-phrase de plus) ; la différence avec un intervalle de confiance n'est pas expliquée ; « 0 61,1 % » (écrire « 0 : 61,1 % » ou « échec 61,1 % ») ; les symboles « ==> » et « --> » ne sont jamais expliqués.

---

## Ce qui m'a aidé

- **« Pendant du test du khi-deux »** : tout de suite je sais à quoi ça correspond dans mon cours.
- **« contrairement au khi-deux, qui demande des effectifs attendus de 5 »** : je reconnais la règle du cours, et j'apprends qu'ici ce n'est pas un problème.
- **« ils sont 2,16 fois plus probables si les proportions sont les mêmes que si elles diffèrent »** et **« 106 fois plus probables »** : c'est la phrase qui m'a fait comprendre le facteur de Bayes.
- **« Verdict : indécis. Les données ne permettent de conclure ni à une différence, ni à son absence. »** : clair et honnête ; ça m'empêche d'écrire « pas de différence » comme on fait souvent avec p > 0,05.
- **« (non tranché) »** et **« (soutien fort) »** à côté des pourcentages : je sais tout de suite comment lire le chiffre.
- **« Ce sont des probabilités d'hypothèses, pas des p-values : aucun seuil de 5 % ne s'y applique. »** : utile, même si je ne sais pas encore par quoi remplacer le seuil.
- **« La cause d'une différence dépend de la façon dont l'étude a été menée »** : bon rappel pour la partie discussion du compte rendu.
- **Les proportions observées à la fin**, en pourcentages avec n : c'est la partie que je comprends le mieux et que je mettrais dans mon tableau.
- **La notation « temoin = froid < chaud »** : très lisible, on voit le classement d'un coup d'œil.
- **« fixé avant de voir les données »** et **« l'étape 5 les fait varier pour voir si la réponse tient »** : ça répond à la question « est-ce que la méthode triche ? ».
