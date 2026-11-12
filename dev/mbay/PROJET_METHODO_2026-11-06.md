# Projet méthodologique de m.bay, 6 novembre 2026

Document de travail pour le coordinateur. Rien ici n'est arbitré par Antoine. **Tâche B déjà reportée dans METHODOLOGIE §3 et §4. Tâche A : la consistance d'information « attendue » au §A5 est fausse (RELECTURE_ADVERSE_2026-11-06.md, B1).** Aucun code ni autre document n'a été modifié.

**Marqueurs de statut des chiffres.**
- **[exécuté]** : calculé pour ce document. Scripts dans le dossier temporaire de session `C:\Users\masse\AppData\Local\Temp\claude\dd1494da25364fcc\09ecfbf3-1f60-4a99-9729-1c85cc9d368d\scratchpad\` (`het_bf.py`, `sens.py`, `chk_r.R`, `chk_py.py`, sorties `het_bf.out`) ; à copier dans `dev/mbay/scripts/` si on les cite (voir I6 de la relecture).
- **[main]** : dérivé à la main, sans exécution.
- **[relecteur]**, **[recalculé]**, **[agent]** : repris de `RELECTURE_ADVERSE_2026-11-05.md` ou `VALIDATION.md`, non refait ici.
- **[non exécuté]** : attendu ou conjecturé.

**Limite de lecture.** Les pages ont été lues par un résumeur automatique, sauf les PDF de Kruschke (2018) et de Barbieri et al. (2016), dont j'ai extrait et lu des passages. Je n'ai eu le texte intégral ni de Moreno et al. (1999), ni de Bertolino et al. (2000), ni de Wetzels et al. (2009), ni de Girón et del Castillo (2021).

---

# Tâche A. Variances par groupe sans Stan

## A1. Conclusion franche

1. **Le modèle n'existe pas tel quel dans ce que j'ai pu vérifier.** Sa structure est classique, son réglage est propre à KefiR. Aucune source consultée ne combine k ≥ 3 groupes, un g à Cauchy multivariée, un τ partagé et un a priori plat sur les rapports d'échelle.
2. **Il est calculable sans Stan, et mieux qu'attendu.** Sachant les rapports d'échelle ρ et g, tout est conjugué : c'est un `anovaBF` pondéré. Le BF se réduit à une intégrale de dimension k (k-1 rapports + g). À ρ = 1, ma formule reproduit `anovaBF` : log BF 4,641607 contre 4,641607 (k = 3) et −0,5881446 contre −0,5881446 (k = 2) [exécuté]. La quadrature concorde avec un Monte-Carlo indépendant (§A5) [exécuté].
3. **La justification par Berger, Pericchi et Varshavsky (1998) ne couvre pas le partage complet.** Elle couvre (μ, τ) et pas les k-1 rapports. Ce qui tient pour ces derniers est un autre argument, plus faible : limite de lois propres identiques sous H0 et H1, et orthogonalité. La constante s'élimine, mais **le choix de la mesure sur les rapports reste une convention**, dont j'ai mesuré l'effet (§A3).
4. **Il lève le blocage B4 sur le plan testé, puissance non mesurée.** H0, n = (5, 5, 40), sd = (4, 4, 1) : P > 0,95 dans 0/300 jeux, contre 16 à 20 % pour la voie gaussienne [exécuté ; `agent_12`, `agent_21`].
5. **La moyenne arithmétique pour τ² est le point faible.** Un troisième groupe très dispersé, sans rapport avec l'écart observé, fait passer le BF de 46 à 2,9 (un jeu) ; la moyenne géométrique reste entre 43 et 58 (§A5) [exécuté].
6. **Ce prior impropre interdit le BF contre le modèle à variance commune.** Il n'autorise donc pas `robust = "bma"` entre structures de variance. Pour cela, il faut une loi propre sur les rapports, comme Maier et al.

## A2. Littérature contrôlée

| Source | Vérification | Ce qu'elle fait | Écart avec le modèle KefiR |
|---|---|---|---|
| Moreno, Bertolino, Racugno (1999). Default Bayesian analysis of the Behrens-Fisher problem. *JSPI* 81(2), 323-333. DOI 10.1016/S0378-3758(99)00070-1 | Crossref ; texte non lu (ScienceDirect 403). Contenu de seconde main : Barbieri et al. (2016) | BF intrinsèques et fractionnaires pour éviter la difficulté des a priori impropres, **2 groupes** | Lois intrinsèques, pas de g JZS |
| Moreno, Girón, Martínez (1998). Model selection with vague prior information. *Rev. R. Acad. Cienc.* 92(4), 289-298 | EuDML, résumé. Pas de DOI affiché | Méthodes intrinsèque et fractionnaire, illustrées sur Behrens-Fisher | Idem |
| Bertolino, Racugno, Moreno (2000). Bayesian model selection approach to analysis of variance under heteroscedasticity. *JRSS D* 49, 495-502. DOI 10.1111/1467-9884.00249 | Crossref ; **résumé partiel, page refusée** | Le plus proche par le titre. Un résumé de recherche dit que l'accent porte sur homoscédasticité contre hétéroscédasticité | **À lire d'abord.** Je ne peux pas dire s'il teste les moyennes pour k > 2 |
| Cano, Carazo, Salmerón (2017). Objective Bayesian model selection approach to the two way analysis of variance. *Comput. Stat.* 33, 235-248. DOI 10.1007/s00180-017-0727-1 | Crossref ; contenu par résumé de recherche, PDF refusé | ANOVA à 2 facteurs, homo- et hétéroscédastique, a priori intrinsèques | Autre famille d'a priori. À lire |
| Girón, del Castillo (2021). A Bayesian solution to the Behrens-Fisher problem. *RACSAM* 115(4). DOI 10.1007/s13398-021-01098-0 | Dépôt RIUMA (notice et résumé) | BF en forme fermée par modèle hiérarchique, consistant, **2 groupes** selon le résumé | Autre a priori, non lu. Bon point de comparaison numérique |
| Maier, Bartoš, Quintana, Dablander, van den Bergh, Marsman, Ly, Wagenmakers (2025). Model-averaged Bayesian t tests. *PBR* 32(3), 1007-1031. DOI 10.3758/s13423-024-02590-5 | Europe PMC + texte PMC12092555 (formules abîmées à l'extraction) | **Le plus proche en esprit.** μ, σ² impropres communs ; ρ ~ Beta(1,5 ; 1,5) propre ; δ ~ Cauchy(0 ; 1/2) ; moyenne sur égales et inégales, normal et t. Dans JASP et RoBTT. **2 groupes**, le passage à k > 2 est annoncé comme travail futur | Loi **propre** sur le rapport, autre échelle de δ, échelle 1/2 et non √2/2 |
| Dablander, van den Bergh, Wagenmakers, Ly (2024). Default Bayes factors for testing the (in)equality of several population variances. *Bayesian Analysis* 19(3), 699-723. DOI 10.1214/23-BA1369 | Crossref ; HTML arXiv 2003.06278v2 | Teste les **variances**, k ≥ 2. Haar droit sur les communs, Dirichlet propre sur les rapports. Argument : constante commune, donc limite de lois propres | Autre question, mais c'est la source la plus claire sur le raisonnement. Son résultat suppose ϑ propre |
| Wetzels, Raaijmakers, Jakab, Wagenmakers (2009). How to quantify support for and against the null hypothesis. *PBR* 16(4), 752-760. DOI 10.3758/PBR.16.4.752 | Europe PMC, résumé ; **texte non lu** | Le résumé annonce le cas de groupes sans variance commune, par Savage-Dickey | **Modèle exact inconnu de moi** : des sources secondaires se contredisent (variances propres contre σ commune ou groupée). Lire p. 757 environ |
| Fu, Hoijtink, Moerbeek (2020). Sample-size determination for the Bayesian t test and Welch's test using the approximate adjusted fractional Bayes factor. *Behav. Res. Methods*. DOI 10.3758/s13428-020-01408-1 | Résumé par recherche ; page PMC bloquée. **Tome et pages non contrôlés** | AAFBF, 2 groupes, paquet `bain` | Fractionnaire approché, autre logique |
| Barbieri, Marin, Florin (2016). A fully objective Bayesian approach for the Behrens-Fisher problem using historical studies. arXiv:1611.06873 | Notice arXiv + PDF lu en extrait | Écrit ne connaître qu'une proposition BF entière antérieure : Moreno et al. (1999) | A priori construit sur données historiques |
| Berger, Pericchi, Varshavsky (1998). Bayes factors and marginal distributions in invariant situations. *Sankhyā A* 60, 307-321 | Résumé sur la page de Duke ; biblio : citation dans Hendriksen et al. **Pas de DOI trouvé, texte non lu** | Résultat valable pour des modèles invariants sous **le même groupe** ; condition portant sur le maximal invariant commun | Voir §A3 |
| Hendriksen, de Heide, Grünwald (2021). Optional stopping with Bayes factors. *Bayesian Analysis* 16(3), 961-989. DOI 10.1214/20-BA1234 | Dépôt CWI, résumé | Mots-clés : invariance de groupe, Haar droit sur les nuisances | Appui sur (μ, τ) seulement |
| Kass, Wasserman (1995). *JASA* 90(431), 928-934. DOI 10.1080/01621459.1995.10476592 | Résumé par recherche | Nuisance rendue orthogonale sous H0, **même loi marginale de nuisance sous les deux hypothèses**, a priori à information unitaire | Appui direct, mais asymptotique |
| Kass, Vaidyanathan (1992). Approximate Bayes factors and orthogonal parameters. *JRSS B* 54, 129-144. DOI 10.1111/j.2517-6161.1992.tb01868.x | Crossref seulement. **Contenu non vérifié**, cité de mémoire | Faible sensibilité du BF à la loi d'une nuisance orthogonale | À lire avant de le citer |
| Bayarri, Berger, Forte, García-Donato (2012). *Ann. Statist.* 40(3), 1550-1577. DOI 10.1214/12-AOS1013 | Crossref + notice arXiv | Critères de choix d'a priori (la notice ne nomme pas les critères) | Contenu précis non vérifié, sauf citation par Dablander pour le « predictive matching » |

**NON VÉRIFIÉES.**
- Dayal & Dickey (1976), « Bayes factors for Behrens-Fisher problems », *Sankhyā B* 38, 315-328 : vue seulement dans la bibliographie de Maier et al. Absente de Crossref. Elle y est présentée, avec Dickey, comme le travail pionnier : à chercher.
- Gönen, Johnson, Lu, Westfall (2005), *The American Statistician* 59(3), 252-257 : existence confirmée par des notices, mais **pas de DOI contrôlé**. Le résumé relève du t **à variance commune**. Ce n'est donc pas une source pour les variances inégales.
- Jeffreys (1961) traiterait le cas de deux moyennes à erreurs différentes : souvenir personnel, non contrôlé.
- Wang & Liu (2015) : seulement vu cité. Pas de lien établi avec les variances inégales.

**Recherche négative.** Aucun résultat pour une ANOVA bayésienne de Welch à k > 2 avec a priori de Cauchy, par deux requêtes seulement. Cela n'exclut pas son existence. Bertolino et al. (2000) et Cano et al. (2017) restent à lire.

## A3. Le partage de p(σ_i²) impropres tient-il ?

**Décomposition** [main]. Soit u_i = log(σ_i/σ_k), i < k, et τ une échelle homogène de degré 1 en σ. Alors ∏ d log σ_i = d log τ · ∏ du_i (jacobien 1, quel que soit le τ choisi). L'a priori du modèle s'écrit donc :

    p(μ, τ, u) ∝ (1/τ) dμ dτ · ∏ du_i

Deux facteurs, deux statuts.

**(μ, τ) : l'argument de groupe tient.** Le groupe affine commun y → cy + b agit sur H0 et H1, et (1/τ) dμ dτ est son Haar droit. C'est le cadre de Berger et al. (1998), repris par Hendriksen et al. et par Dablander et al. Il faut pour cela que a/τ soit un paramètre invariant : c'est le cas pour tout τ homogène.

**u : l'argument de groupe ne tient pas.** Un changement d'échelle par groupe y_ij → c_i y_ij ne laisse pas H0 invariante : des moyennes égales ne le restent pas. Il n'y a donc pas de groupe plus grand. Citer Berger et al. (1998) pour le partage de tous les p(σ_i²) serait **faux**, et la documentation ne doit pas le faire.

**Ce qui tient à la place.**
- **Limite de lois propres.** Soit une loi propre π_L(u) identique sous H0 et H1, troncature à |u| < L par exemple. Alors BF_L → BF si les deux intégrales convergent, et la constante disparaît pour une raison simple : u a le même sens sous H0 et H1, et π1(μ, τ, u, a) = π0(μ, τ, u) π(a | τ, u). C'est le raisonnement que Dablander et al. appliquent à (μ, τ). Il s'étend à u **une fois admis que u est commun**. Le résultat exact de ces auteurs exige toutefois leur partie propre.
- **Orthogonalité.** En loi normale, moyennes et variances sont orthogonales au sens de Fisher [main], ce qui correspond à la condition de Kass & Wasserman : même marginale de nuisance sous les deux hypothèses. L'influence de la loi sur u devrait alors décroître avec n. Argument asymptotique seulement.
- **Convergence** [main] : les intégrales sont finies si **n_i ≥ 2 et SS_i > 0 pour tout i**. Avec n_i = 1, sous H0, le facteur 1/|μ − ȳ_i| n'est pas intégrable : **BF indéfini**. m.bay doit refuser. À n_i = 2, le BF est défini mais sa loi est très plate.

**Ce qui ne tient pas : l'unicité.** Une autre mesure h(u) du donne un autre BF, car BF = E[BF(u) | y, H0], avec un postérieur dépendant de h. En k = 2, avec ρ = τ1/(τ1+τ2) en précisions, **le plat sur u est exactement Beta(0, 0)** sur ρ [main], limite impropre de la famille de Maier et al. (Beta(1,5 ; 1,5)) et de Dablander et al. (Beta(1/2 ; 1/2)).

Effet mesuré sur le BF10, k = 2, τ arithmétique, un jeu par ligne [exécuté, `sens.py`] :

| Plan | Beta(0,0), le plat | Beta(½,½) | Beta(1,5 ; 1,5) | Beta(5,5) |
|---|---:|---:|---:|---:|
| n = 5, 5 ; sd 1, 1 | 1,650 | 1,823 | 2,030 | 2,197 |
| n = 5, 5 ; sd 1, 3 | 0,717 | 0,756 | 0,803 | 0,793 |
| n = 5, 40 ; sd 4, 1 | 0,705 | 0,751 | 0,845 | 1,178 |
| n = 15, 25 ; sd 1, 2 | 38,79 | 39,42 | 39,67 | 35,77 |
| n = 40, 40 ; sd 1, 3 | 0,351 | 0,353 | 0,357 | 0,371 |

Le facteur reste de 1,02 à 1,23 entre le plat et Beta(1,5 ; 1,5). Le cas déséquilibré passe de 0,705 à 1,178 avec une loi informative, ce qui change le côté de 1. **Une ligne de sensibilité sur cette loi s'impose donc à petit n.** Cinq jeux : ordre de grandeur, pas une table.

**Conséquence non négociable** : le plat sur u est absent du modèle à variance commune, donc **BF(hétéro contre homo) n'est pas défini**. La remarque B.6 de la relecture adverse reste entière. Pour passer de l'aiguillage à une moyenne, il faudra une loi propre sur u.

## A4. Définitions de τ

Toutes reviennent à σ si les variances sont égales. Toutes laissent le prior sur u inchangé (§A3), et seule change la loi de a sachant u.

| τ² | Formule | k = 2 : δ standardise par | Si σ_1 → ∞ | Effet attendu [main, sauf §A5] |
|---|---|---|---|---|
| Arithmétique (proposé) | mean(σ_i²) | √((σ1²+σ2²)/2) | τ → ∞ | Échangeable, sans dépendance à n. **Un seul groupe dispersé élargit l'a priori sur tous les contrastes**, donc pousse vers H0 (Lindley) |
| Géométrique | exp(mean log σ_i²) | √(σ1σ2) | τ croît en σ_1^(1/k) | Échangeable. Effet ralenti, invariant à l'échelle propre de chaque groupe ; aucune lecture usuelle en taille d'effet |
| Harmonique | 1/mean(1/σ_i²) | moyenne harmonique | τ borné | Protège contre un groupe dispersé. Mais un groupe très précis (σ_i → 0) écrase τ et tout écart devient énorme : **pousse vers H1** |
| Pondérée par n | Σ n_i σ_i² / N | Proche de Maier et al. (formule reconstruite, extraction abîmée, à vérifier sur le PDF) | τ → ∞ | Prior dépendant du plan, non échangeable. Comparable à RoBTT seulement si cette lecture se confirme |
| Pondérée n-1 | Σ(n_i−1)σ_i²/Σ(n_i−1) | S_p classique | τ → ∞ | Idem. La lecture usuelle du d de Cohen |
| Information unitaire | n̄ · mean(σ_i²/n_i) | Proche de l'écart-type de Welch multiplié par √n̄ | τ → ∞ si n_1 petit | Le BF suit plutôt le t de Welch. Prior dépendant du plan |
| Groupe de référence | σ_ref² | Δ de Glass | Selon le groupe | Rompt l'échangeabilité. N'a de sens qu'avec un témoin désigné |
| Échelle par groupe | a_i ~ N(0, g σ_i²) | Pas d'échelle unique | Sans objet | Incompatible avec P = I − J/k sans redéfinir l'identifiabilité |

Mon avis, au vu du §A5 : à k = 2, la moyenne arithmétique se défend (lecture simple : d_av). **À k ≥ 3, je ne la prendrais pas comme défaut sans autre preuve**, car elle laisse un groupe sans rapport décider du BF d'une paire. La géométrique est la meilleure candidate à ce stade, au prix d'une lecture moins usuelle et de P < 0,25 presque absent sous H0 dans le premier scénario (2 %, contre 32 %). Choix à trancher par une simulation avec effet, pas sur ces seuls contre-exemples.

## A5. Contrôles exécutés [exécuté]

**Forme calculable.** Avec ρ_i = σ_i/τ, D = diag(ρ_i²/n_i), V_g = D + gP, Q_g = Σ SS_i/ρ_i² + ȳ'(V_g⁻¹ − V_g⁻¹11'V_g⁻¹/(1'V_g⁻¹1))ȳ, N = Σ n_i :

    BF(u, g) = [ |V_g|^(-1/2) (1'V_g⁻¹1)^(-1/2) Q_g^(-(N-1)/2) ] / [ même chose à g = 0 ]
    m0(u)    ∝ ∏ ρ_i^(-(n_i-1)) · |D|^(-1/2) (1'D⁻¹1)^(-1/2) Q_0^(-(N-1)/2)       (invariant si ρ → cρ)
    BF       = ∫ m0(u) [∫ BF(u, g) p(g) dg] du / ∫ m0(u) du

V_g est diagonale plus rang 1, donc chaque point coûte O(k) par Sherman-Morrison (non programmé ici, j'ai utilisé `solve` sur des piles de matrices).

| Contrôle | Résultat |
|---|---|
| δ a priori, k = 2 : médiane de \|δ\| contre √2/2 (400 000 tirages) | 0,7091 contre 0,7071. L'annonce Cauchy(√2/2) est **juste** (aussi [main]) |
| ρ = 1 contre `anovaBF` (k = 3, k = 2) | log BF identiques à 10⁻⁶ |
| Quadrature contre E[BF(σ) \| y, H0], tirages exacts sous H0 | k = 2 : 0,30190 / 0,30195 (se 0,00004). k = 3, H0 : 0,32459 / 0,32543 (se 0,00043). k = 3, effet : 3,616 / 3,573 (se 0,076) |
| Grille divisée par 2 | Inchangé à 5 chiffres |

**Piège de calcul.** L'estimateur simple E[BF(σ) | y, H0] est exact et sans MCMC, mais son ESS chute dès qu'il y a un effet (1 974 sur 20 000 sur la dernière ligne, [exécuté]). Sa variance peut devenir **infinie** quand l'écart inter dépasse SS_i (n_i·d² > SS_i, d étant l'écart de ȳ_i à μ) : argument heuristique sur le comportement à σ_i → 0, sans démonstration complète [main]. Il ne doit donc pas servir de moteur par défaut. Quadrature pour k ≤ 3 ou 4 ; au-delà, QMC ou échantillonnage d'importance à valider [non exécuté].

**Simulation H0, à comparer au 16 % de la voie gaussienne.**

Grille de 31 points par axe, r = 1/2, π = 1/2, mêmes jeux pour toutes les colonnes. Scénario du blocage B4 : n = (5, 5, 40), sd = (4, 4, 1), N = 300.

| τ² | Prob. médiane | P > 0,5 | **P > 0,95** | P < 0,25 |
|---|---:|---:|---:|---:|
| Arithmétique | 0,300 | 14,7 % | **0 %** | 32,3 % |
| Géométrique | 0,358 | 18,3 % | **0 %** | 2,0 % |
| Harmonique | 0,408 | 22,0 % | **0 %** | 0 % |
| Pondérée par n | 0,374 | 18,7 % | **0 %** | 1,3 % |
| Information unitaire | 0,192 | 9,0 % | **0 %** | 64,3 % |
| *Voie gaussienne, autre tirage* | *0,539* | *52,3 %* | *16 % (20 % à grand N)* | *[recalculé, `agent_12`, `agent_21`]* |

Scénario inverse : n = (5, 5, 40), sd = (1, 1, 4), N = 150. P > 0,95 à 0 % partout ; P > 0,5 de 5,3 % (pondérée par n) à 14,0 % (harmonique) ; P < 0,25 à 78,7 % (arithmétique), 52,7 % (géométrique), 8,0 % (harmonique), 84,0 % (pondérée), 74,7 % (information unitaire).

Lecture : **le modèle lève les faux « Oui » francs sur ce plan, pour tous les τ** (0/300 : borne supérieure à 95 % d'environ 1 %). Il reste 9 à 22 % de P > 0,5, ce qui confirme que ce seuil n'est pas un verdict. **La puissance n'est pas mesurée** : un 0 % sous H0 peut aussi venir d'une voie trop prudente, et aucune comparaison à Welch sous effet n'a été faite.

**Le groupe dispersé, test direct de τ.** Groupes 1 et 2 fixes (n = 20, écart d'environ 1 σ), groupe 3 multiplié par S sans toucher à sa moyenne. **Un seul jeu.** BF10 :

| S | Arithmétique | Géométrique | Harmonique |
|---:|---:|---:|---:|
| 1 | 46,4 | 46,3 | 46,2 |
| 3 | 40,0 | 43,2 | 43,7 |
| 10 | 24,8 | 55,4 | 59,9 |
| 30 | 9,49 | 57,6 | 69,1 |
| 100 | **2,90** | 48,1 | 71,7 |

Pondérée par n et information unitaire coïncident ici avec l'arithmétique (n égaux). **Sous τ arithmétique, un groupe sans rapport avec la différence observée fait passer le BF de 46 à 2,9.** Géométrique et harmonique ne s'effondrent pas, mais ne sont pas stables pour autant.

**Propriétés non établies** [non exécuté] : consistance sous H0 et H1, consistance d'information (BF → ∞ quand les écarts → ∞ à SS fixes, attendue point par point en u, non démontrée après intégration), calibration sous l'a priori (impossible à simuler avec un plat sur u sans troncature), non-normalité, voie multi-facteurs (tous les BF ci-dessus sont à un facteur).

## A6. Ce qui reste une convention de KefiR

Tout cela doit être écrit comme tel dans la documentation :
1. La mesure plate sur les rapports u (Beta(0,0) à k = 2).
2. La définition de τ, et donc le sens de δ.
3. r = 1/2 repris du cas homoscédastique, sans preuve de calibration ici.
4. Un seul g partagé pour tous les contrastes.
5. Pas de réduction à `anovaBF` : à données homoscédastiques, ce BF **diffère** du BF à σ commune, puisque les rapports restent intégrés. Il y a coïncidence seulement à ρ = 1 fixé.

## A7. Désaccords avec le modèle tel que formulé

1. **Sur la justification** : oui pour (μ, τ) ; non pour u, où Berger et al. (1998) ne s'applique pas, et il faut le dire (votre formulation le pressent, je le confirme par l'action de groupe).
2. **Sur « p(σ_i²) partagés entre H0 et H1 »** : exact, mais à condition de **refuser n_i < 2**, ce que le modèle ne disait pas.
3. **Sur τ arithmétique** : désaccord à k ≥ 3 (BF de 46 à 2,9 sous l'effet d'un groupe sans rapport, §A5). Elle exige au moins une simulation avec effet face à la géométrique.
4. **Sur la comparabilité avec JASP** : RoBTT prend Cauchy(0 ; 1/2) et non √2/2 ; une autre échelle de δ, et une loi propre. Le modèle KefiR **ne reproduira pas** JASP en variances inégales, même à 2 groupes. L'argument de D5 (« comparable à JASP ») ne s'étend pas ici.
5. **Sur le k = 2** : le calcul est juste, mais c'est un d_av. Un utilisateur qui connaît Cohen attendra S_p. L'écart est nul en plan équilibré et grand en déséquilibré [main].

---

# Tâche B. Texte prêt à remplacer les sections 3 et 4 de METHODOLOGIE

> Ci-dessous : texte de remplacement. Chiffres marqués selon la légende. Les marqueurs sont à retirer une fois les contrôles faits.

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

Motif : 2³ sous H0, 7 termes, au moins un terme à P > 0,5 dans 57 à 74 % des jeux en lecture appariée [relecteur] ; en moyenne 1,27 (N = 16) et 0,97 (N = 40) terme à P > 0,5 [agent, `agent_20`, sur `anovaBF`, 30 jeux à N = 16, effectif à N = 40 non relu ici]. La lecture cohérente ne supprime pas le problème : sous l'uniforme, à N = 16, chaque effet principal dépasse 0,5 dans 50 à 57 % des jeux, parce qu'il part de 0,737 [agent, `agent_20`, **30 jeux seulement**, à refaire avant citation]. Beta-binomial à partir de 3 facteurs : proposition de D13, **non définie sous marginalité** [relecteur, I12], à valider.

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

---

# Tâche C. Dix références « métadonnées seules »

**Tirage.** Lignes de BIBLIOGRAPHIE.md contenant exactement `niveau : métadonnées. ` (71 lignes, soit le compte annoncé), `random.Random(20261106).sample(x, 10)` en Python [exécuté]. Graine = date du jour, fixée avant de voir le résultat.

| Référence | Source consultée | Verdict |
|---|---|---|
| [Rouder2009] | Résumé Europe PMC ; sources secondaires ; **texte non obtenu** | **Non établi pour « r sur δ ».** Le résumé ne parle ni de JZS ni de r. **Attribution probablement à corriger** : §5.1 écrit « Cauchy(0, 0,707) [Rouder2009] » ; les sources secondaires rattachent cet article à une Cauchy d'échelle 1, et 0,707 au « medium » des logiciels. À confirmer dans le texte |
| [Wetzels2012] | Résumé (recherche) | **Conforme.** Précision utile : le test applique la sélection de variables à des g-priors variés |
| [OHagan1995] | Résumé (recherche), Crossref | **Imprécis.** BF fractionnaire : exact. « Fondement de BFpack » : la vignette indique que le défaut de BFpack est **ajusté** (Mulder 2014), avec O'Hagan en option. Pages 99-118 selon Crossref ; 99-138 sur un autre PDF (discussion comprise, à vérifier) |
| [Hoeting1999] | Résumé (dépôt Columbia) | **Conforme** : l'article parle d'inférences trop confiantes en ignorant l'incertitude de sélection. Un correctif de 2000 existe pour les formules |
| [Kruschke2018] | **PDF lu par extraits** | **Erreur de fond à l'usage.** ±0,1 sur δ : confirmé. **Mais pour un coefficient de régression, l'article donne ±0,05 sur variables standardisées**, alors que §4 prend 0,1 pour la pente, en s'autorisant de la même convention. L'article décide aussi par intervalle et ROPE, sans atome, et avertit que le seuil dépend de l'application |
| [Kelter2020] | Résumé Europe PMC | **Conforme.** Titre tronqué (« a non-technical introduction to Bayesian inference with JASP ») |
| [Schwarz1978] | Résumé (recherche) | **Conforme** |
| [KassWasserman1995] | Résumé (recherche) | **Conforme.** Précision : le lien normal-Schwarz est à erreur O(n^(-1/2)), et une Cauchy exige un terme correctif. Utile à §6 : le BIC ne fournit pas directement le BF de type JZS |
| [Rowe2002] | Crossref, notices et tables partielles | **Partiel.** Sous-titre manquant : *Models for Source Separation and Signal Unmixing*. Sujet conjugué et Wishart inverse confirmé, section NIW non localisée. 2002 selon Crossref, 2003 sur des notices |
| [vanDoorn2020] | Résumé arXiv 1712.06941 | **Conforme**, k > 2 absent du résumé aussi |

**Bilan.** Aucune erreur de métadonnées. 1 erreur de fond à l'usage (pente, Kruschke). 1 attribution probablement fausse (r = 0,707, Rouder, à confirmer). 2 imprécisions (O'Hagan et BFpack ; Rowe). 6 conformes. Hormis Kruschke, le contrôle s'est fait sur résumés, pas sur le texte : il ne valide pas un point précis.

**À ajouter à BIBLIOGRAPHIE.md si le coordinateur retient la tâche A** : Moreno1999, Bertolino2000, Maier2025, Dablander2024, Giron2021, Hendriksen2021, Wetzels2009, Cano2017, Barbieri2016, avec les niveaux du tableau A2.
