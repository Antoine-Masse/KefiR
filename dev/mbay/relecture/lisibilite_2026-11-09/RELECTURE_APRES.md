# Relecture « après » — lecteur débutant (BUT génie biologique)

Lecteur : étudiant de BUT GB, bac+1/2. J'ai vu en cours le t de Student, l'ANOVA, le khi-deux, les p-values et Tukey, mais jamais de statistique bayésienne.
Fichiers lus : `avant/m.bay_*_fr.txt` et `apres/m.bay_*_fr.txt` (7 plans), rien d'autre. Je ne juge pas le calcul, seulement ce que je comprends.

---

## 1. Plan par plan (sorties « après »)

### Plan 1 : deux groupes

- **Temps pour trouver la conclusion** : environ 10 s. La ligne vide avant « 7) Bilan » se voit tout de suite. Ensuite, il me faut quelques secondes de plus pour comprendre que « Aucune différence… 73,9 % » ne veut **pas** dire « pas de différence ».
- **Ce que j'en comprends** : on ne peut pas dire si A et B diffèrent. Les données penchent un peu vers « pas de différence », mais pas assez pour l'affirmer.
- **Phrases qui m'arrêtent encore** :
  - « 5) Solidité de la réponse [m.bay() {KefiR}] - facteur de Bayes 0,354, la réponse peut s'inverser. » : quelle réponse ? Je n'ai encore lu aucune réponse à ce stade.
  - « Note : facteur de Bayes : combien de fois les données sont plus probables avec une différence que sans. » : j'ai lu juste avant « 0,354, soit 2,82 en faveur de l'absence ». Il me faut deux lectures pour faire le lien entre 0,354 et 2,82.
  - « Note : variés : r de 0,25 à 1, chances de départ de 1 contre 2 à 2 contre 1, répartition des variances. » : style télégraphique, je ne sais pas ce qui a été « varié » ni pourquoi.
  - « ==> Taille attendue r = 0,5 » : pour moi, r est un coefficient de corrélation. Ici, je ne sais pas ce que c'est.
  - « écart un peu ramené vers 0 » : ramené par qui, et pourquoi ?
  - « Formes du bruit » et « queues lourdes » : la sortie « avant » les expliquait, l'« après » ne le fait plus.
- **Ce qui manque pour mon compte rendu de TP** : les écarts-types par groupe (pour écrire « moyenne ± écart-type »). Surtout, la phrase « pour un compte rendu » s'écrit « Aucune différence de y selon g : probabilité a posteriori 73,9 % ». Si je la copie, j'écris qu'il n'y a pas de différence, alors que le verdict est « indécis ».

### Plan 2 : quatre groupes

- **Temps** : environ 5 à 8 s. Le titre du bilan suffit : « y diffère selon g : 99,8 % (soutien très fort) ».
- **Ce que j'en comprends** : les groupes diffèrent nettement. a et c vont ensemble, b et d aussi, et b et d sont plus hauts que a et c.
- **Phrases qui m'arrêtent** :
  - « ==> Semblables : a/c (modéré) ; b/d (faible). » : « semblables, faible » veut-il dire « probablement semblables » ou « on ne sait pas » ? Dans le plan 7, une paire à 35 % est classée « Non tranchés ». Ici, b/d, à 26,8 %, est classée « Semblables ». Je ne vois pas où passe la limite.
  - « Note : une seule loi pour toutes les paires : pas de correction, les réponses ne peuvent pas se contredire. » : on m'a appris qu'il faut corriger (Tukey, Bonferroni). « Une seule loi » ne m'explique pas pourquoi on peut s'en passer ici.
  - « ==> a/c : différents à 19,5 % » : la tournure est compréhensible, mais il me faut un temps pour la lire comme « probabilité que a et c diffèrent ».
- **Ce qui manque pour le TP** : les écarts-types par groupe, et des lettres de groupes (a, a, b, b) à mettre sur un graphique. « a = c < b = d » s'en approche beaucoup.

### Plan 3 : deux facteurs

- **Temps** : environ 15 à 20 s. Le titre du bilan (« au moins un facteur : 99,0 % ») ne dit pas lequel. Il faut descendre jusqu'aux lignes A, B, A:B.
- **Ce que j'en comprends** : A a un effet (A2 au-dessus de A1). Pour B et pour l'interaction, on ne peut pas trancher.
- **Phrases qui m'arrêtent** :
  - « 5 combinaisons de termes sont possibles (une interaction n'entre qu'avec ses facteurs : hérédité). » : « hérédité » me fait penser à la génétique.
  - « 2) Contrôles des données [quantile() {stats}] - 1 valeur éloignée, gardée. » alors que le corps dit « 6 au-delà de 1,5 écart interquartile ». Alors, 1 ou 6 ?
  - « 4) Point de départ (a priori) - aucune différence du tout : 33,3 %. » : dans tous les autres plans, la valeur est 50 %. Et « Note : ces chances ne s'additionnent pas » me perd.
  - « 5) Solidité de la réponse [lmBF() {BayesFactor}] - B change de côté selon l'a priori. » et « Note : variés : r de 0,25 à 1, trois lois sur les modèles. » : « change de côté » et « trois lois sur les modèles », je ne comprends ni l'un ni l'autre.
  - « ==> A:B (les écarts selon A changent avec B) : non tranché (39,3 % ; a priori 16,7 % ; facteur de Bayes 3,24). » : le facteur de Bayes est supérieur à 1 (il est donc en faveur de l'interaction ?), mais la probabilité est sous 50 %. La sortie « avant » expliquait ce paradoxe, l'« après » ne le fait plus.
  - « 6) Quels groupes diffèrent - non disponible à plusieurs facteurs. » : la case A1 x B2 (8,4) a l'air plus basse que les autres, et je ne peux pas le vérifier. Que dois-je faire à la place ?
- **Ce qui manque pour le TP** : une phrase « pour un compte rendu » **par facteur**. La phrase actuelle ne parle que d'« au moins un facteur », ce qui ne sert à rien dans un rapport de TP. Il manque aussi une piste pour comparer les cases.

### Plan 4 : apparié

- **Temps** : environ 5 à 8 s.
- **Ce que j'en comprends** : « apres » est plus élevé que « avant », d'environ 1,5 (entre 0,6 et 2,2), et c'est très probablement réel.
- **Phrases qui m'arrêtent** :
  - « ==> Taille attendue r = 0,5 sur dz (écart au sein du sujet rapporté à son écart-type). » : « dz » ne m'était pas connu. La parenthèse aide, mais la phrase reste lourde.
  - « 5) … - facteur de Bayes 75,4, l'ordre tient. » : « l'ordre » de quoi ? Au plan 2, la même idée s'écrit « réponse stable », et dans le corps « La réponse ne change pas de côté ».
  - « Chaque sujet est comparé à lui-même » revient trois fois (étapes 1, 3 et 4).
- **Ce qui manque pour le TP** : la moyenne des différences et leur écart-type (ce que j'écrirais pour un t apparié), et les écarts-types par condition.

### Plan 5 : mesures répétées (3 temps)

- **Temps** : environ 5 à 8 s.
- **Ce que j'en comprends** : y augmente de T1 à T2 puis à T3. T1/T3 est sûr, T1/T2 et T2/T3 le sont moins (« modéré »).
- **Phrases qui m'arrêtent** :
  - « Note : conditions supposées également corrélées au sein d'un sujet (non testé ; m.test() contrôle la sphéricité). » : je ne connais ni la sphéricité ni m.test(). Est-ce un problème pour mes données, oui ou non ?
  - « soit 0,88 écart-type résiduel au sein des sujets » : trop long, et répété trois fois.
  - « Sensibilité : la réponse en tête va de 99,9 % à > 99,9 % » : se lit, mais « Quel que soit le point de départ : > 99,9 % » (plan 6) est plus clair.
- **Ce qui manque pour le TP** : les écarts-types par temps, et des lettres pour un graphique.

### Plan 6 : ANCOVA

- **Temps** : environ 8 à 10 s.
- **Ce que j'en comprends** : à x égal, t1 est plus bas que t2 et t3, et t2 et t3 sont peut-être pareils.
- **Phrases qui m'arrêtent** :
  - « ==> Covariable x : étendue commune 1,15 à 9,03, au moins 80,0 % de chaque groupe. » : je ne sais pas si c'est bien ou mal, ni ce que je dois en faire.
  - « Note : variés : r de 0,25 à 1 (rscaleCont avec lui), chances de départ de 1 contre 2 à 2 contre 1. » : « rscaleCont » est du code, pas du français.
  - « facteur de Bayes 10^6 » : un étudiant de première année ne lit pas forcément 10^6. « Plus d'un million de fois » serait plus parlant.
  - « Sensibilité : la réponse en tête va de > 99,9 % à > 99,9 %. » : ça fait sourire.
  - « covariable(s) x » : le pluriel entre parenthèses fait automatique.
  - « Note : x est dans tous les modèles : les groupes sont comparés à x égal ; son effet propre n'est pas la question posée. » : utile, mais j'aimerais quand même savoir si x a un effet.
- **Ce qui manque pour le TP** : les **moyennes ajustées** (à x = 4,71). Seules les moyennes brutes sont données, alors que tout le reste de la conclusion porte sur les valeurs ajustées. Je ne sais donc pas quels chiffres mettre dans mon tableau.

### Plan 7 : réponse oui/non

- **Temps** : environ 12 à 15 s. Le titre du bilan tient sur deux lignes, et « germe » apparaît d'un coup (je suppose que c'est le nom de la variable).
- **Ce que j'en comprends** : on ne peut pas dire si le traitement change la proportion de germination. « chaud » a l'air plus haut (80 % contre 60 %), mais ce n'est pas tranché.
- **Phrases qui m'arrêtent** :
  - « 7) Bilan - Les proportions de germe ne diffèrent pas selon traitement : 68,4 % (non tranché). » : même problème qu'au plan 1, la phrase affirme une absence de différence alors que le verdict est « indécis ».
  - « Note : pas de forme du bruit à comparer avec des effectifs. » : je ne comprends pas cette phrase.
  - « [Dirichlet-multinomiale {KefiR}, hdi() {bayestestR}] » : jargon pur.
  - « Note : variés : l'étalement attendu des proportions, et les chances de départ (1 contre 2 à 2 contre 1). »
- **Ce qui manque pour le TP** : les proportions en % pour **chaque** groupe. Seuls les effectifs sont donnés, plus une seule paire en %. La sortie « avant » donnait « oui 60,0 % » pour chaque groupe.

---

## 2. Comparaison avant / après (note de 1 à 5)

| Critère | Avant | Après | Commentaire |
|---|:-:|:-:|---|
| Lisibilité | 2 | 4 | Les lignes coupées à environ 80 caractères et les titres d'étape résumés (« - rien à retirer », « - facteur de Bayes 623, réponse stable ») changent tout. Avant, les lignes de 300 caractères étaient illisibles dans la console. |
| Longueur | 2 | 4 | Environ 25 % de caractères en moins (par exemple 4489 contre 3313 pour le plan 1), même s'il y a plus de lignes. Les « Note : » se multiplient un peu (six dans le bilan du plan 6). |
| Repérage du bilan | 2 | 5 | La ligne vide avant « 7) Bilan », la conclusion dans le titre, le « Verdict : » puis les hypothèses numérotées : je trouve la réponse en moins de 10 s. Avant, il fallait 20 à 40 s et lire « Conclusion : les hypothèses classées, ce sur quoi chacune repose, son soutien ». |
| Jargon | 2 | 3 | Zellner-Siow, HOP, Beta(1,1), graine, Aiken et West, échantillonnage d'importance et « 9,1e-08 point(s) » ont disparu. Restent r, dz, rscaleCont, sphéricité, m.test(), Dirichlet-multinomiale, hérédité, 10^6 et « trois lois sur les modèles ». Les explications de « forme du bruit » et de « queues lourdes » ont aussi disparu, ce qui est un recul. |
| Cohérence entre plans | 3 | 4 | Les 7 étapes sont les mêmes partout (avant, le plan 3 s'arrêtait à 6), avec le même squelette de bilan et la même phrase pour le compte rendu. Restent les formulations différentes pour « la réponse tient » et les catégories Semblables / Non tranchés, qui n'ont pas la même limite d'un plan à l'autre. |

---

## 3. Les 10 remarques les plus importantes (par ordre de priorité)

1. **Cas « indécis » : la phrase du compte rendu et le titre du bilan affirment l'absence de différence** (plans 1 et 7). Je les recopierais tels quels, et mon rapport dirait le contraire du verdict.
   - Titre : « 7) Bilan - indécis : aucune différence légèrement favorisée (73,9 %), non tranché. »
   - Compte rendu : « Données non concluantes : l'absence de différence de y selon g est légèrement favorisée (probabilité a posteriori 73,9 %, facteur de Bayes 2,82 en faveur de l'absence), mais la conclusion peut s'inverser selon l'a priori (m.bay {KefiR}…). »

2. **Plan 3 : il faut une réponse par facteur dans le compte rendu, et une explication pour A:B.**
   - Compte rendu : « Effet de A : probabilité a posteriori 98,8 %, facteur de Bayes 82,4 ; B et l'interaction A:B : non tranchés (60,1 % et 39,3 %). »
   - Sous A:B : « Les données vont dans le sens d'une interaction (3,24 fois plus probables avec elle), mais pas assez pour compenser sa faible chance de départ (16,7 %). »
   - Étape 6 : remplacer « non disponible » par ce qu'il faut faire, par exemple « comparer les cases : analyser B à l'intérieur de chaque niveau de A ».

3. **Donner une seule fois une échelle de lecture du facteur de Bayes et des mots « faible / modéré / fort / très fort »**, au moment où le facteur de Bayes apparaît : « Lecture : 1 à 3 faible, 3 à 10 modéré, 10 à 100 fort, plus de 100 très fort ; en dessous de 1, il penche vers l'absence de différence (0,354 = 1/2,82). » Sans cette échelle, je ne sais pas pourquoi 98,7 % est « fort » et 99,8 % « très fort ».

4. **Un seul mot pour dire que la conclusion tient.** On lit aujourd'hui « réponse stable », « l'ordre tient », « La réponse ne change pas de côté », « l'ordre ne change pas » et « B change de côté ». Je propose, partout, « conclusion stable (ne dépend pas de l'a priori) » et « conclusion fragile : peut s'inverser selon l'a priori ».

5. **Mêmes catégories pour les paires dans tous les plans, avec leur limite.** Aujourd'hui, b/d à 26,8 % est classée « Semblables (faible) » au plan 2, et temoin/chaud à 35,1 % « Non tranchés » au plan 7. Je propose trois cases fixes : « Différents (≥ 75 %) », « Semblables (≤ 25 %) », « Non tranchés (entre les deux) », et j'écrirais « Non tranchés : b/d (différents à 26,8 %) ».

6. **Écrire la note de sensibilité en français, et les « chances » en pourcentages.** Pour « Note : variés : r de 0,25 à 1, chances de départ de 1 contre 2 à 2 contre 1, répartition des variances. », je propose « Vérifié en changeant le point de départ : taille d'effet attendue (petite à grande) et chance de départ d'une différence (33 % à 67 %). » Pour r : « taille d'effet attendue : moyenne (r = 0,5, valeur habituelle) », pour qu'on ne le confonde pas avec une corrélation.

7. **Remettre une explication courte de « forme du bruit » et de « queues lourdes »** (elle existait dans la sortie « avant ») : « Formes du bruit (comment les valeurs se dispersent autour de la moyenne de leur groupe) » ; « queues lourdes (plus de valeurs extrêmes que la loi normale) ». Au plan 7, supprimer « Note : pas de forme du bruit à comparer avec des effectifs. » ou la remplacer par « Avec des oui/non, il n'y a pas de dispersion à modéliser : seules les proportions comptent. »

8. **Chasser le code et le jargon qui restent** : « (rscaleCont avec lui) » à supprimer ; « m.test() contrôle la sphéricité » à remplacer par « non vérifié ici » ; « Dirichlet-multinomiale {KefiR} » à remplacer par « comparaison exacte des proportions {KefiR} » ; « hérédité » à remplacer par « une interaction n'est testée qu'avec ses deux facteurs » ; « 10^6 » à remplacer par « plus d'un million » ; « de > 99,9 % à > 99,9 % » à remplacer par « > 99,9 % quel que soit le point de départ » ; « trois lois sur les modèles » à remplacer par « trois façons de répartir les chances de départ entre les combinaisons ». Corriger aussi le titre de l'étape 2 du plan 3 (« 1 valeur éloignée » alors que le corps en compte 6 au seuil de 1,5) : « 6 valeurs éloignées dont 1 très éloignée, gardées ».

9. **Expliquer en une phrase les deux notes de l'étape 6 qui surprennent un étudiant formé à Tukey.**
   - « pas de correction » : « Pas besoin de correction type Tukey ou Bonferroni : toutes les paires sont estimées ensemble, à partir des mêmes regroupements, donc elles restent cohérentes entre elles. »
   - « écart un peu ramené vers 0 » : « l'écart estimé est un peu plus petit que l'écart observé : la méthode est volontairement prudente. »

10. **Ajouter ce qu'il faut pour remplir un tableau de TP** : moyenne ± écart-type par groupe ou condition (tous les plans) ; moyennes **ajustées** à x = 4,71 pour l'ANCOVA ; moyenne et écart-type des différences pour l'apparié ; % de oui pour chaque groupe au plan 7 ; si possible des lettres de groupes (« a = c < b = d » donne a : A, c : A, b : B, d : B) à reporter sur un graphique.

### Ce qu'il ne faut surtout pas changer

- **La ligne vide avant « 7) Bilan »** et le **résumé dans le titre de chaque étape** (« - rien à retirer », « - facteur de Bayes 623, réponse stable »). C'est ce qui me permet de lire la sortie en diagonale.
- **« Verdict : … »** suivi des **hypothèses numérotées avec leur %** (« 1. y diffère selon g : 99,8 % »).
- **« Selon le point de départ : X % à Y % ; l'ordre ne change pas / peut s'inverser »** : c'est simple, et je comprends tout de suite si la conclusion est solide.
- **La notation des regroupements « a = c < b = d »** : c'est plus clair que tout ce que j'ai vu avec Tukey.
- **La phrase « pour un compte rendu »** (à corriger dans les cas indécis, mais le principe est excellent).
- **« probabilités d'hypothèses, pas des p-values »** et **« cela ne dit ni que la différence est importante, ni pourquoi »** : ces deux phrases m'empêchent de mal interpréter.
- **L'explication des crochets** (« 95 % de chances d'y trouver l'écart »), les **« points »** pour les écarts de proportions, le rappel **« pendant du test du khi-deux »** et la ligne **« Écart le plus marqué : oui, 60,0 % en temoin contre 80,0 % en chaud »**.
- **Le découpage en 7 étapes identiques dans tous les plans** et les **lignes courtes**.
