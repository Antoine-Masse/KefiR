# Relecture « lecteur débutant » : cas D à I (mesures répétées, ANCOVA)

Lecteur joué : étudiant de BUT génie biologique (bac+1/2). J'ai eu un cours de statistiques classiques (t de Student, ANOVA, p-value, Tukey), mais jamais de bayésien. J'ai lu seulement `recits_repete.txt` et `recits_ancova.txt`.

---

## Cas D : tension avant / après, 14 patients (une mesure manquante)

**1. Ce que j'ai compris.** On ne peut pas dire si le traitement change la tension. Elle baisse un peu (environ 3 à 4 points), mais l'intervalle va de -7,5 à +0,58, donc il contient 0. En cours, j'aurais écrit « non significatif ».

**2. Ce qui m'a arrêté.**
- « mesuré sur les mêmes 14 sujets ; N = 27 » puis « Sujets : 13 ». Je n'ai compris qu'à l'étape 2 que 14 × 2 = 28, moins une mesure, ça fait 27. Ensuite on passe à 13 sujets et N = 26. Trois effectifs différents en dix lignes.
- « 1 sujet(s) sans mesure aux deux niveaux retiré(s) (3) ». Le « (3) » veut dire quoi ? Le patient numéro 3, ou 3 valeurs ?
- « queues lourdes et variances inégales entre conditions ne sont pas encore validées avec des sujets ». Je ne comprends pas « validées avec des sujets ». Qui valide quoi ? Je comprends seulement qu'on suppose la loi normale.
- « elle n'est pas mise à l'épreuve des autres » : les autres quoi ?
- « r = 0,5 sur la différence standardisée au sein des sujets (dz) » et surtout « égal à lmBF(y ~ temps + id) à 0,2 près en échelle log ». C'est illisible pour moi. Je ne sais pas ce qu'est r, et encore moins l'échelle log.
- « Facteurs de Bayes des 2 modèles » : quels deux modèles ? Je suppose « différence » et « pas de différence », mais ce n'est écrit nulle part.
- « Terme(s) qui changent de côté : moment (40,7 % à 77,6 %) ». J'ai fini par comprendre que ça passe sous et au-dessus de 50 %, mais il m'a fallu relire.
- **« Les données déplacent peu ses cotes (multiplient par 1,73) ; seul l'a priori (50,0 %) le fait pencher vers le oui. »** C'est la phrase qui m'a le plus gêné. On m'a dit plus haut que 50 %, c'est « pile ou face : la méthode ne prend pas parti ». Comment un pile ou face peut-il faire pencher vers le oui ? Pour moi, ce sont les données (×1,73) qui font passer de 50 % à 63,4 %. J'y vois une contradiction.
- « cotes » : je ne connais pas ce mot en statistique. Il est expliqué au cas E (« rapport oui/non »), mais pas ici, alors que D est le premier cas.
- « Écart estimé apres - avant = -3,4 », alors que les moyennes donnent 134,1 - 138,4 = -4,3. La phrase « légèrement ramenée vers zéro par l'a priori » m'a sauvé, mais j'ai d'abord cru à une erreur.
- « soit -0,45 [...] écart-type des différences (dz) » puis plus bas « delta = -0,45 » : dz et delta, c'est la même chose ? Deux noms pour un seul chiffre.
- « ce n'est pas un intervalle de confiance » : d'accord, mais alors je l'appelle comment dans mon compte rendu, et qu'est-ce qui change pour moi ?

**3. Ce qui me manque pour le compte rendu.** Une phrase modèle du genre : « La tension baisse de 3,4 mmHg en moyenne (IC95 % crédibilité [-7,5 ; 0,6]), mais les données ne permettent pas de conclure (probabilité d'un effet 63 %). » Il manque aussi l'unité (mmHg). Je ne sais pas non plus si « 63,4 % » est « presque » une preuve ou pas du tout, ni quoi conseiller (plus de patients ?).

---

## Cas E : temps de réaction à 4 sessions, 12 sujets

**1. Ce que j'ai compris.** Le temps de réaction change selon la session (95,9 %, soutien fort) : il diminue de S1 à S4. La seule différence nette est entre S4 et S1, et aussi entre S4 et S2. Pour les autres paires, on ne sait pas.

**2. Ce qui m'a arrêté.**
- « 2 au-delà de 3 » écarts interquartiles, puis « Aucune valeur n'est retirée », puis « la lire avec réserve ». Concrètement, je fais quoi ? Je dois le signaler dans mon compte rendu ? Je regarde ces valeurs ? Lesquelles sont-ce ?
- « symétrie composée » et « (m.test() contrôle la sphéricité par le test de Mauchly) ». Faut-il que je lance m.test() en plus ? La phrase me laisse croire que l'analyse est incomplète.
- « effet aléatoire d'échelle 1 (Rouder et al. 2012) », « Zellner-Siow », « chaque réponse se lit contre le modèle avec les sujets seuls » : trois termes de jargon de suite. Je les ai sautés.
- « Regroupements les plus probables (15 comparés) : S4 < S3 = S2 = S1 (22,1 %) ». Le « plus probable » n'a que 22 %, c'est peu. Je ne sais pas si je peux l'écrire comme conclusion.
- **Contradiction apparente.** Le regroupement le plus probable dit S3 = S1, mais la ligne des paires dit « S1 / S3 : différents, soutien faible (différents à 66,1 %) ». Juste au-dessus, on me promet que « les réponses ne peuvent pas se contredire ». Pour moi, elles se contredisent.
- « S1 / S3 : différents » avec l'intervalle [-24 ; 0,98], qui contient 0. Avec mon cours, un intervalle qui contient 0 veut dire « pas de différence prouvée ». Même chose pour S3 / S4 [-22 ; 2,4]. À l'inverse, au cas G, E2 / E3 a un intervalle qui exclut 0 [0,061 ; 3,2] pour « seulement » 83,9 %. Je ne sais pas lequel des deux croire.
- Les seuils : 44 % donne « non tranché », 66 % donne « soutien faible », 84,7 % et 91,6 % donnent « modéré », 95,9 % donne « fort ». Où est le tableau des seuils ? Sans lui, les mots ont l'air arbitraires.
- « aucune correction paire par paire n'est nécessaire » : bien, c'est ce qui remplace Tukey ? Si oui, il faut le dire, car c'est le mot que je connais.

**3. Ce qui me manque pour le compte rendu.** Des lettres comme avec Tukey (S1 a, S2 ab, S3 ab, S4 b) ou un équivalent que je peux mettre sur un graphique. Il me manque aussi l'unité (ms ?) et une phrase sur l'effet d'apprentissage (la baisse est régulière).

---

## Cas F : score, groupe (entre sujets) × séance (intra), 2 × 3

**1. Ce que j'ai compris.** On ne peut rien conclure : ni l'effet du groupe, ni l'effet de la séance, ni l'interaction ne sont tranchés. Pourtant, les moyennes par case montrent que le groupe « entraîné » part plus bas en s1 (48,4 contre 54,0), puis rattrape. J'aurais parié sur une interaction.

**2. Ce qui m'a arrêté.**
- « une réponse par terme, obtenue en pesant ensemble les 5 combinaisons de termes » : je ne vois pas quelles sont les 5 combinaisons.
- **« Chance de chaque terme, pris seul (ces chances ne s'additionnent pas) : groupe 50,0 %, seance 50,0 %, groupe:seance 16,7 % ; aucune différence du tout 33,3 % »**. Je n'ai rien compris. Pourquoi 16,7 % ? Pourquoi 33,3 % et pas 50 % comme dans les autres cas ? « Ne s'additionnent pas » me perd encore plus.
- « une proportion ~ Beta(1,1) par ordre d'interaction (HOP, Taylor-Rodriguez et al. 2016) » : jargon pur. Seule la fin, « ce qui protège des fausses interactions », est utile.
- « deux calculs avec des graines différentes diffèrent jusqu'à 2,3 points de probabilité ». Dans les autres cas, on me dit « négligeable ». Ici il n'y a aucun jugement : est-ce grave ? Et c'est quoi, une graine ?
- « croisé avec trois lois sur les modèles (HOP, 1/2 par terme, uniforme) » : jargon.
- Encore « seul l'a priori (50,0 %) le fait pencher vers le oui » (pour groupe et séance) : même contradiction qu'au cas D.
- « l'a priori (16,7 %) penche de l'autre côté : non tranché » : là, je comprends, parce que 16,7 % est bien en dessous de 50 %. Mais du coup, si l'interaction partait à 50 %, elle serait « oui » ? On dirait que c'est la méthode qui décide, pas les données.
- « Différence de score selon groupe, seul ou en interaction ? » : « seul ou en interaction » ne me parle pas. En ANOVA, j'avais un effet principal et une interaction séparés.
- Il n'y a pas d'étape « quels niveaux diffèrent » (on passe de 5 à 6 = conclusion). La phrase « ce bilan ne dit pas quels niveaux diffèrent » le reconnaît, mais ne me dit pas comment le savoir.

**3. Ce qui me manque pour le compte rendu.** Une phrase qui commente le profil des moyennes par case (le groupe entraîné part plus bas). Un conseil : « avec 10 sujets par groupe, l'étude manque de données pour trancher ». Le score n'a pas d'unité, et il manque l'écart des moyennes par case.

---

## Cas G : rendement selon l'engrais, à surface égale

**1. Ce que j'ai compris.** À surface égale, l'engrais E3 donne un meilleur rendement que E1 et E2 (environ +2 et +1,6), et E1 et E2 sont semblables. Les moyennes brutes faisaient croire que E2 était bien plus bas : c'est parce que ses parcelles sont plus petites (7,93 contre 10,70).

**2. Ce qui m'a arrêté.**
- « étendue commune aux groupes de 5,5 à 13 ; elle contient au moins 53,3 % de chaque groupe (E1 le moins) » : je comprends le chiffre, mais pas s'il est bon ou mauvais. Au cas I, 6,7 % déclenche une alerte. Où est la limite ?
- « Covariable(s) centrée(s) sur leur moyenne [...] (Aiken et West 1991) » : compréhensible, mais je ne vois pas pourquoi on me le dit.
- « pente a priori de Cauchy d'échelle sqrt(2)/4 (rscaleCont, le « medium » de BayesFactor), partagée avec ses interactions avec un facteur » : illisible.
- « Verdict : une différence existe probablement » (90,2 %, modéré) contre « très probablement » au cas E (95,9 %). La nuance est bien, mais toujours sans tableau des seuils.
- « les paires sont plus bas » : en fait, elles sont aussi plus haut (étape 6), en double. J'ai lu deux fois la même chose.
- « Moyennes observées [...] (brutes, non ajustées) » : il manque les **moyennes ajustées**. Pour mon graphique, j'ai besoin de E1, E2, E3 à surface 9,52, pas seulement des écarts.

**3. Ce qui me manque pour le compte rendu.** Les moyennes ajustées par engrais et la pente du rendement selon la surface (combien de rendement en plus par unité de surface). Il manque aussi les unités (t/ha ? m² ?).

---

## Cas H : score selon la méthode, avec la note d'entrée, pentes différentes

**1. Ce que j'ai compris.** Ce que je crois comprendre : l'effet de la méthode dépend de la note d'entrée. Les élèves faibles (note autour de 11) réussissent mieux en « classique », les forts (note 16) mieux en « projet ». À la note moyenne, il n'y a pas de différence nette. Mais le verdict me dit « score diffère selon methode : > 99,9 % », et j'ai d'abord compris « projet est meilleur », ce qui est faux.

**2. Ce qui m'a arrêté.**
- **« score diffère selon methode, à même note : > 99,9 % (soutien très fort) »**, alors que les moyennes sont 12,8 et 13,1 et qu'à note = 14 l'intervalle contient 0. Le chiffre « > 99,9 % » est en fait porté par l'interaction (« seul ou en interaction »), mais la ligne du verdict ne le dit pas. Un étudiant pressé recopie cette ligne et se trompe.
- « équivalent bayésien de Johnson et Neyman » : jamais entendu.
- « ils diffèrent de 6,8 à 12 et de 14 à 18 », alors que l'étendue commune est « de 7,67 à 17,5 ». Donc 6,8 et 18 sont hors des données ? Et à note = 14, l'intervalle [-0,11 ; 1,2] contient 0, alors que la zone « commence à 14 ». Ça m'a paru contradictoire (c'est sans doute un arrondi).
- « projet - classique : note = 11 : -0,99 » : il m'a fallu du temps pour lire que le signe moins veut dire « classique meilleur ». Une phrase en clair aiderait.
- « Les écarts selon methode varient-ils avec note » : bonne formulation, je la comprends mieux que « interaction ».
- « les paires sont plus bas » : il n'y a que deux groupes.
- « Par case : classique = 12,8 ; projet = 13,1. » répète exactement la ligne du dessus.
- Les cotes « multipliées par 1710 » et « par 4330 » : c'est énorme, mais je n'ai aucun repère (est-ce comme p < 0,001 ?).

**3. Ce qui me manque pour le compte rendu.** La phrase en clair « projet est meilleur pour les notes d'entrée au-dessus d'environ 14, classique en dessous d'environ 12 ». Les deux pentes (score gagné par point de note dans chaque méthode). Un graphique des deux droites serait idéal.

---

## Cas I : poids selon le régime, avec l'âge, groupes qui se recouvrent peu

**1. Ce que j'ai compris.** Le régime R2 semble donner un poids plus élevé (+4,5) à âge égal. Mais les deux groupes n'ont presque pas le même âge (26,8 ans contre 41,9 ans), donc on compare des gens qui n'existent pas dans les données. Il faut s'en méfier.

**2. Ce qui m'a arrêté.**
- « étendue commune aux groupes de 32,8 à 33,8 ; elle contient au moins 6,7 % de chaque groupe » : clair, et l'alerte qui suit est très utile.
- **« Verdict : une différence existe probablement »** avec « soutien faible ». « Probablement » et « faible » dans la même conclusion, ça m'a paru contradictoire. Le cas G a la même probabilité (90,2 %) et dit aussi « probablement », mais avec « modéré ». Je ne sais pas ce que je dois retenir.
- « ni quels groupes diffèrent » : il n'y a que deux groupes.
- Les groupes sont comparés à « age = 34,3 », justement là où presque personne n'a de données. Le texte le dit (« extrapolation »), mais je ne sais pas quoi écrire à la place : faut-il refaire l'étude ?
- « sensitivity = FALSE » : pourquoi l'avoir coupé ici ? Est-ce que je dois le remettre ?

**3. Ce qui me manque pour le compte rendu.** Une phrase du type : « Âge et régime sont confondus dans ce plan : on ne peut pas séparer l'effet de l'un de l'effet de l'autre. » Plus une recommandation (choisir des groupes du même âge). L'unité (kg) manque aussi.

---

## Mes 10 remarques classées (de la plus gênante à la moins gênante)

1. **Cas H, verdict « > 99,9 % » sur methode, alors que les moyennes sont presque égales.** Je conclus à tort que « projet » est meilleur.
   *Proposition* : « Verdict : l'effet de la méthode dépend de la note d'entrée (pentes différentes, > 99,9 %). Il n'y a pas de méthode meilleure pour tous : classique en dessous d'environ 12, projet au-dessus d'environ 14, pas de différence nette entre les deux. »

2. **« seul l'a priori (50,0 %) le fait pencher vers le oui » (cas D et F).** Cela contredit « pile ou face : la méthode ne prend pas parti ».
   *Proposition* : « Les données font passer la probabilité de 50 % à 63 % (cotes × 1,73) : c'est trop peu pour trancher. »

3. **Le regroupement contredit la ligne des paires (cas E : S3 = S1 dans le regroupement, mais S1 / S3 « différents »), juste après « les réponses ne peuvent pas se contredire ».**
   *Proposition* : pour toute paire entre 30 % et 70 % environ, écrire « non tranché », et non « différents, soutien faible ». Ou ajouter : « le regroupement le plus probable n'a que 22 % : aucun regroupement ne s'impose ».

4. **Les mots de soutien (faible, modéré, fort) n'ont pas de seuils visibles.** Je ne peux pas justifier ces mots dans mon compte rendu.
   *Proposition* : une ligne fixe en fin de conclusion, par exemple « Échelle : non tranché < 75 % ≤ faible < 90 % ≤ modéré < 95 % ≤ fort < 99 % ≤ très fort » (avec les vrais seuils).

5. **Probabilité contre intervalle qui contient 0 (cas E, S1 / S3 « différents » avec [-24 ; 0,98]).** Avec mon réflexe de cours, c'est contradictoire.
   *Proposition* : « L'intervalle contient 0 : la taille de l'écart reste incertaine, même si une différence est un peu plus probable qu'absente. »

6. **Pas de phrase modèle pour le compte rendu ni d'unités.** Je ne sais pas comment rédiger ni comment nommer l'intervalle.
   *Proposition* : terminer par « À écrire : La tension baisse de 3,4 (intervalle de crédibilité à 95 % : -7,5 à 0,6) ; probabilité d'une différence 63 %, non tranché (facteur de Bayes 1,73). »

7. **Cas F, les a priori « 50 % / 50 % / 16,7 % / 33,3 % (ces chances ne s'additionnent pas) » et « Beta(1,1) [...] HOP ».** C'est incompréhensible.
   *Proposition* : « Au départ, une interaction est jugée moins probable (16,7 %) qu'un effet simple (50 %), pour éviter d'en voir partout. Les données doivent donc la soutenir plus fort. »

8. **Cas I, « probablement » avec « soutien faible ».**
   *Proposition* : « Verdict : une différence semble exister, mais elle n'est pas fiable : les deux groupes n'ont pas le même âge (26,8 contre 41,9 ans), l'âge et le régime sont confondus. »

9. **Cas D, les effectifs 14 / 27 / 13 / 26 et « retiré(s) (3) ».**
   *Proposition* : « 14 patients, dont 1 (patient n° 3) sans mesure « après » : retiré. Analyse sur 13 paires (26 mesures). »

10. **Jargon technique dans le récit principal** (« lmBF(y ~ temps + id) à 0,2 près en échelle log », « Cauchy d'échelle sqrt(2)/4 (rscaleCont) », « Zellner-Siow », « Rouder et al. », « Johnson et Neyman », « graine », « symétrie composée »).
    *Proposition* : renvoyer tout cela vers `print(resultat, detail = TRUE)` et garder une phrase simple, par exemple « Taille d'effet attendue : moyenne (r = 0,5), réglage standard. » Ajouter aussi les **moyennes ajustées** en ANCOVA (cas G : « Moyennes à surface égale (9,52) : E1 = …, E2 = …, E3 = … ») : c'est ce que je mettrais dans mon graphique.

Remarques de moindre importance : « les paires sont plus bas » alors qu'il n'y a que 2 groupes (H, I) ; les doublons « Par case » (H) et étape 6 / conclusion (E, G) ; dz et delta pour un même chiffre (D) ; la précision de 2,3 points sans jugement (F) ; « m.test() contrôle la sphéricité » qui laisse croire qu'il faut lancer un autre test (E, F) ; « validées avec des sujets » (D, E, F).

---

## Ce qui m'a aidé et qu'il ne faut pas toucher

- La structure en étapes numérotées, avec « Une seule question : y a-t-il une différence... ? » au début : je sais tout de suite ce qu'on cherche.
- « Point de départ (a priori) » et « pile ou face : la méthode ne prend pas parti » : c'est la meilleure explication de l'a priori que j'aie lue.
- Le format « Réponse : Oui / Non tranché, probabilité du oui X %, a priori 50 % » : c'est facile à recopier.
- « Cela ne dit ni qu'elle est importante, ni pourquoi » : un rappel utile.
- « Ce sont des probabilités d'hypothèses, pas des p-values : aucun seuil de 5 % ne s'y applique » : indispensable pour quelqu'un comme moi.
- « multiplient ses cotes (rapport oui/non) par 23,5 (facteur de Bayes) » : enfin une définition du facteur de Bayes que je comprends (à placer dès le cas D).
- « légèrement ramenée vers zéro par l'a priori » : cela m'a évité de croire à une erreur de calcul.
- « Aucune valeur n'est retirée » et le comptage des valeurs éloignées selon Tukey : je connais.
- Cas H : « Changement de stratégie : [...] un écart ajusté unique n'aurait pas de sens ». C'est clair et honnête.
- Cas I : l'alerte sur le recouvrement et le « soutien plafonné à faible », expliqué avec sa raison.
- Cas G : la covariable moyenne par groupe à côté des moyennes brutes. C'est ce qui m'a fait comprendre pourquoi E2 « remonte ».
- « Les écarts selon methode varient-ils avec note ? » : bien plus parlant que « interaction ».
- « La réponse ne change pas de côté quand r ou les chances de départ changent » : rassurant et compréhensible.
