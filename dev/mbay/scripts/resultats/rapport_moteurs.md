# Moteurs candidats pour m.bay() : mesures (5 novembre 2026)

Tout ce qui suit a réellement tourné, sauf mention contraire dans « Ce qui n'a pas pu tourner ».
Scripts : `C:\Users\masse\Desktop\KefiR\KefiR-mbay\dev\mbay\scripts\agent_*.R`.
Sorties brutes : `...\dev\mbay\scripts\resultats\agent_*.txt` (sessionInfo complets à la fin de agent_04 et agent_05, reproduit plus bas ; agent_02 arrêté avant son sessionInfo).
Bibliothèque : `C:\Users\masse\Desktop\KefiR\mbay-lib` uniquement (247,2 Mo à la fin).

## Environnement et versions

R 4.6.1 (2026-06-24 ucrt), x86_64-w64-mingw32, Windows 11 x64 (build 26300), locale French_France.utf8.

| Paquet | Version |
|---|---|
| BayesFactor | 0.9.12-4.8 |
| bayestestR | 0.20.0 (insight 1.7.0, datawizard : voir agent_01) |
| BAS | 2.0.2 |
| bridgesampling | 1.2-1 |
| rstanarm | 2.32.2 |
| rstan / StanHeaders | 2.32.7 / 2.39.1.1 |
| brms | 2.23.0 |

Compilateur : `Sys.which("g++")` vide, `make` absent du PATH, `pkgbuild::has_build_tools()` = FALSE.
Rtools 4.4 existe dans `C:\RBuildTools\4.4` (variable RTOOLS44_HOME définie) mais hors PATH, et ce n'est pas la version prévue pour R 4.6.

## 1. Installation

Binaires CRAN pour R 4.6 Windows : disponibles pour tous, versions binaire = source pour les 8 candidats.
Fermeture récursive des Imports de KefiR : 41 Imports, 216 paquets hors base.

| Paquet | Dép. récursives | Au-delà de KefiR | Installés ici (paquet inclus) | Durée | Taille ajoutée | Compilateur pour ajuster ? |
|---|---|---|---|---|---|---|
| BayesFactor | 22 | 5 (coda pbapply hypergeo elliptic contfrac) | 6 | 35 s | 16,6 Mo | non |
| bayestestR | 2 | 2 (insight datawizard) | 3 | 13 s | 6,5 Mo | non |
| BAS | 0 | 0 | 1 | 7 s | 2,1 Mo | non (pas testé numériquement) |
| bridgesampling | 19 | 2 (Brobdingnag coda) | 2 | 9 s | 2,9 Mo | non |
| rstanarm | 117 | 37 | 26 | 189 s | 206,8 Mo | non (vérifié, g++ absent) |
| brms | 75 | 30 | 7 après rstanarm | 49 s après rstanarm | 12,3 Mo après rstanarm | **oui** (échec vérifié sans compilateur) |

Réserves :
- Durées et tailles sont **incrémentales**, dans cet ordre, sur cette machine. L'ordre favorise brms, qui hérite de rstan, StanHeaders, BH, etc. installés par rstanarm.
- rstanarm : 37 dépendances au-delà de KefiR, mais 13 déjà présentes dans la bibliothèque principale (297 paquets). Un étudiant qui n'a que KefiR paiera donc plus que 26 paquets, 189 s et 206,8 Mo : shiny, DT, etc. s'ajoutent.
- BH pèse à lui seul 123,2 Mo (rstanarm 21,3 ; BayesFactor 13,2 ; StanHeaders 11,4 ; brms 8,6).
- Union des six au-delà de KefiR : 52 paquets.

## 2a. ttestBF contre intégration JZS (Rouder 2009, agent_02)

| n1,n2 | d simulé | t | ttestBF | integrate() | écart relatif |
|---|---|---|---|---|---|
| 20,20 | 0 | -0,3847 | 0,32755749 | 0,32755749 | 8,1e-13 |
| 20,25 | 0,6 | 3,1433 | 12,619268 | 12,619268 | 2,0e-12 |
| 100,100 | 0,3 | 1,7695 | 0,66193194 | 0,66193194 | 9,1e-11 |

## 2b. anovaBF contre calcul indépendant

Mêmes a priori (Rouder 2012) : contrastes orthonormés, g ~ IG(1/2, r²/2), r = 1/2, mu plat, sigma Jeffreys. Comparaison strictement comparable.

| Cas | BayesFactor | integrate() | écart relatif | propErr |
|---|---|---|---|---|
| 3 gr, n=15, eff=0 | 0,20496348 | 0,20496347 | 5,4e-08 | 8,1e-05 |
| 3 gr, n=15, eff=0,5 | 0,92388933 | 0,92388929 | 3,6e-08 | 5,6e-05 |
| 3 gr, n=10,15,20, eff=0,5 | 22,397595 | 22,397594 | 4,6e-08 | 1,0e-04 |
| y~A+B (2 g, double intégrale) | 0,323062 | 0,327914 | 1,48e-02 | 8,8e-03 |
| y~x continu (Liang, r=√2/4) | 497,077 | 497,077 | 1,8e-09 | |

Un seul g : accord quasi parfait. Deux g : BayesFactor passe par Monte-Carlo, l'écart de 1,5 % est de l'ordre de 1,7 fois le propErr annoncé.

rstanarm + bridgesampling sur y~F (n=60, stan_glm, a priori par défaut, 2 chaînes, 2000 it) : **BF10 = 0,207, contre anovaBF = 3,46 sur les mêmes données**. Les a priori ne sont pas identiques (a priori par défaut faiblement informatifs de rstanarm, très larges sur l'échelle des effets, contre ZS r=1/2), donc ce n'est pas une vérification de BayesFactor. Mais c'est un résultat utile : le changement de moteur renverse la conclusion, un facteur 17, et change de côté de 1. `error_measures` affichait 0 %, que je ne prends pas pour une garantie (méthode « normal »).

## 2c. y~A*B, inclusion, plans

Données : A×B 2×3, 10 répétitions (`seed 221`). `whichModels="all"` : 2,48 s, 7 modèles + nul.

| Modèle | BF vs nul | propErr |
|---|---|---|
| A | 18,134 | 2,7e-09 |
| B | 0,3146 | 9,9e-05 |
| A:B | 0,4048 | 3,3e-04 |
| A+B | 6,445 | 2,7e-02 |
| A+A:B | 7,918 | 9,6e-03 |
| B+A:B | 0,1299 | 1,3e-02 |
| A+B+A:B | 2,931 | 1,3e-02 |

| Terme | Tous modèles (main) | bayestestR match=FALSE | Appariés (main) | bayestestR match=TRUE |
|---|---|---|---|---|
| A | 19,16 | exp(2,9527) = 19,16 | 18,70 | exp(2,9284) = 18,70 |
| B | 0,3576 | exp(-1,0282) = 0,3576 | 0,3533 | exp(-1,0405) = 0,3533 |
| A:B | 0,4396 | exp(-0,8219) = 0,4396 | **0,4547** | **exp(-0,8175) = 0,4415** |

Calcul à la main = bayestestR partout, **sauf A:B apparié** : (A+B+A:B)/(A+B) = 2,9306/6,4449 = 0,4547 contre 0,4415. Écart non élucidé dans le temps imparti. La définition de bayestestR pour l'interaction reste à lire avant de la citer. Anova classique sur ces données : A p=0,0018 ; B p=0,27 ; A:B p=0,30.

Autres plans (tous sans compilateur) :
- ANCOVA `lmBF(y~F+x)/lmBF(y~x)` : BF=3,95 (propErr 1,3 %) ; p classique du facteur 0,014.
- Apparié `lmBF(y~tm+id, whichRandom="id")/lmBF(y~id,...)` : BF=1,354 (propErr 0,8 %) ; `ttestBF(paired=TRUE)` : 1,274 (propErr 0,02 %) ; p=0,056. **Modèles et a priori différents** : ne pas présenter comme le même calcul.
- Split-plot `anovaBF(y~G*T+id, whichRandom="id")` : 3,65 s ; 4 modèles ; propErr jusqu'à 2,8 %. Interaction G:T appariée : BF=0,573 (propErr 2,0 %).

Temps et erreur, `anovaBF` par défaut (`withmain`), facteurs à 2 niveaux, machine chargée par d'autres calculs (les temps à 3-4 facteurs sont donc des majorants) :

| n | 1 facteur | 2 facteurs (4 modèles) | 3 facteurs (18 modèles) | 4 facteurs (166 modèles) |
|---|---|---|---|---|
| 30 | 0,04 s | 1,70 s ; propErr max 1,1 % | 5,08 s ; max 4,2 % | échec à 60,2 s (limite) |
| 100 | 0,11 s | 1,67 s ; max 2,3 % | 12,97 s ; **max 80,5 %**, médiane 2,0 % | 205 s ; **max 75,3 %**, médiane 4,7 % |
| 500 | 0,04 s | 2,08 s ; max 3,7 % | 21,39 s ; max 4,3 % | non terminé, arrêté au budget (> 5 min) |

La limite de 60 s n'est pas fiable : à n=100, 4 facteurs, le calcul a continué 205 s (le code compilé n'est pas interrompu par `setTimeLimit`). Sur des modèles isolés, propErr dépasse donc parfois très largement 1 % : m.bay devra le lire, relancer (`recompute`) ou avertir.

## 2d. Ce que BayesFactor ne sait pas faire (essayé)

| Essai | Résultat |
|---|---|
| `anovaBF(cbind(y1,y2)~F)`, `lmBF`, `generalTestBF` | erreur « colonnes non définies sélectionnées » |
| `lmBF(y~x+(1+x\|id))` | même erreur, peu parlante |
| `lmBF(y~x+id+x:id, whichRandom=c("id","x:id"))` | **pas d'erreur, BF = NA** (avertissement : 10000 NA retirés) |
| `lmBF(y~G+T+id+T:id, whichRandom=c("id","T:id"))` (facteur:id) | fonctionne : BF=0,28 (±1,9 %) |
| Variances hétérogènes | aucun argument dans `lmBF` ni `anovaBF` ; `ttestBF(var.equal=FALSE)` accepté en silence, retourne un BF JZS à variance commune |
| `lmBF(bin~F, family=binomial())` | **accepté en silence, ignoré** : BF=0,2516328, strictement identique au gaussien |
| y comptage `anovaBF(cnt~F)` | traité comme gaussien, sans avertissement |
| Binaire en table `contingencyTableBF` | fonctionne (BF=0,60), seule voie non gaussienne |
| NA dans y | erreur explicite |
| Une cellule à 1 observation | fonctionne, BF=0,386 |

m.bay doit donc refuser **lui-même** ces plans : BayesFactor ne protège pas (family et var.equal ignorés).

## 2e. Stan

Sans compilateur (g++ absent) :

| Appel | Temps | R-hat max | n_eff min |
|---|---|---|---|
| `stan_glm(y~F)`, n=60, 2 chaînes, 1000 it | 1,8 s | 1,014 | 383 |
| `stan_lmer(y~F+(1\|id))`, n=60, 2 chaînes, 1000 it | 5,8 s | 1,003 | 290 |
| `stan_mvmer(list(y~F+(1\|id), y2~F+(1\|id)))`, 1 chaîne, 500 it | 3,0 s | 1,035 (console : 1,09) | 52 (avertissements ESS) |
| `stan_mvmer(list(y~F, y2~F))` sans effet aléatoire | erreur « No random effects terms » | | |
| `stan_aov(prior=R2(.5))` | erreur (mode de R2 indéfini avec < 3 prédicteurs) ; relance sans résultat (voir plus bas) | | |
| `brm(y~F)` | **échec après 26,2 s** : « 'make' not found », puis « connexion incorrecte » | | |

Bridge sampling sur rstanarm, essayé :
- `bridge_sampler(stanreg)` sans `diagnostic_file` : erreur explicite, il faut l'option à l'ajustement.
- `bridge_sampler(fit$stanfit)` : erreur (« variable does not exist... gamma »).
- `stan_glm(..., diagnostic_file=...)` puis `bridge_sampler` : **fonctionne**, 1,3 s ; log m = -97,7803 contre -96,2068 pour le nul.

Multivarié sans compilateur : `stan_mvmer` seulement, et il exige un effet aléatoire. Pour MANOVA classique, ce n'est pas le bon outil. Le code Stan `brms` de `mvbind` se génère (1,3 s) mais n'a pas pu être ajusté sans compilateur.

## 2f. MANOVA conjuguée exacte (agent_05, sans aucune dépendance)

Vérifications du code :
- p=1 : formule = -36,9892961391 ; integrate() sur sigma² : -36,9892961391 (écart 4e-15) ; t multivariée fermée : écart 7e-15.
- p=2 : formule = -41,49818 ; Monte-Carlo direct (2e5 Sigma ~ IW) = -41,50163 ± 0,00485.

Données : p=3, 3 groupes, n=20/gr. H0 : logBF10 = -7,21 (manova Pillai p=0,73). H1 connu : logBF10 = 11,87 (p=2,7e-06). 5,4 ms par BF.

Sensibilité (logBF10), S0 = cov(Y) × (nu0-p-1) :

| lambda | H0, nu0=5 | H0, nu0=20 | H1, nu0=5 | H1, nu0=8 | H1, nu0=20 |
|---|---|---|---|---|---|
| 0,1 | -13,91 | -13,95 | 6,47 | 6,03 | 4,76 |
| 1/3 | -10,35 | -10,39 | 9,71 | 9,28 | 8,04 |
| 1 | -7,22 | -7,26 | 11,98 | 11,59 | 10,44 |
| 10 | -1,96 | -1,99 | 10,25 | 10,08 | 9,57 |

lambda (précision a priori sur les effets) domine ; nu0 compte moins. Avec S0 = I, nu0 pèse nettement plus sur H1 (6,4 -> 11,6 à lambda=0,1).

Invariance aux unités (Y1 ×1e6, Y2 ×1, Y3 ×1e-3), logBF10 :

| S0 | H0 avant -> après | H1 avant -> après |
|---|---|---|
| identité | -7,22 -> -7,56 | **11,93 -> 4,47** |
| diag(var Y) | -7,21 -> -7,21 | 11,87 -> 11,87 |
| cov(Y) | -7,22 -> -7,22 | 11,98 -> 11,98 |

Transformation linéaire non diagonale : cov(Y) invariant (11,979836 -> 11,979836), diag(var Y) non (11,87 -> 7,76).
Conclusion : S0 fixe = rejet ; **S0 proportionnel à cov(Y)**. Réserve : cet a priori dépend des données (empirique) et utilise cov(Y) totale, effet inclus.

Étalonnage, 200 jeux, même plan, S0=cov(Y), nu0=p+2 :

| lambda | effet | P(H1\|D) médiane | >0,5 | >0,95 | manova p<0,05 |
|---|---|---|---|---|---|
| 1/3 | 0 | 0,000 | 0,000 | 0,000 | 0,050 |
| 1/3 | connu | 0,336 | 0,435 | 0,235 | 0,895 |
| 1 | 0 | 0,002 | 0,005 | 0,000 | 0,065 |
| 1 | connu | 0,925 | 0,720 | 0,440 | 0,895 |

Le calcul est exact, mais très conservateur face à manova avec ces lambda. Le choix de lambda reste un vrai choix à arbitrer (un mélange sur g par intégration 1D, à la Zellner-Siow, serait possible ; **non testé ici**).

## 3. Simulation (agent_03, 400 jeux, 30,6 s, 77 ms/jeu)

y~F, 3 groupes, n=15/gr. « d=0,5 » interprété ainsi : moyennes -0,5 ; 0 ; +0,5 (sd=1). Graine 301. anovaBF, a priori 1/2 sur H1.

| | P(H1\|D) 5/25/50/75/95 % | P>0,5 | P>0,95 | aov p<0,05 |
|---|---|---|---|---|
| H0 | 0,146 / 0,168 / 0,209 / 0,307 / 0,552 | 0,070 [0,039 ; 0,115] | 0,010 | 0,040 [0,017 ; 0,077] |
| Effet | 0,218 / 0,447 / 0,777 / 0,939 / 0,997 | 0,720 [0,652 ; 0,781] | 0,235 [0,178 ; 0,300] | 0,660 [0,590 ; 0,725] |

Corrélation de rang entre P et -p sous effet : 1,000 (plan équilibré, un facteur). Équivalence observée : P>0,5 correspond à p jusqu'à 0,084 ; P>0,95 à p jusqu'à 0,0019.

## Tableau comparatif et recommandation par plan

| Plan | Moteur recommandé | Coût pour install_github |
|---|---|---|
| 2 groupes, 1 facteur, A*B, ANCOVA, apparié, split-plot, Gaussien | BayesFactor en Imports | 5 paquets en plus, ~17 Mo, ~35 s, pas de compilateur |
| Inclusion | à la main dans m.bay | 0 ; bayestestR facultatif (2 paquets, 6,5 Mo), A:B apparié à élucider |
| MANOVA | code conjugué interne | 0 dépendance ; à trancher : lambda ; S0 ∝ cov(Y) |
| Pentes aléatoires continues, hétéroscédasticité, non gaussien | rstanarm en Suggests, appel conditionnel | au moins 26 paquets, 189 s, 207 Mo si l'étudiant accepte ; BF de bridge très dépendant des a priori |
| brms | ne pas retenir par défaut | exige Rtools chez l'étudiant ; échec obscur sans lui |

Pourquoi pas BayesFactor en Suggests : `install_github` ne l'installerait pas, et m.bay serait cassé d'origine pour son usage principal. Pourquoi pas rstanarm en Imports : 207 Mo minimum imposés à tous, pour tous les utilisateurs de m.test aussi.

## Ce qui n'a pas pu tourner

- BAS : installé (0 dépendance, 2,1 Mo), **pas comparé numériquement**, faute de temps.
- `stan_aov` : l'appel de agent_04 a échoué (R2 mode). La relance (agent_07_stan_aov.R) n'a rien produit avant son plafond de 30 s, machine chargée. Fichier `resultats/agent_07_stan_aov.txt` vide : pas de chiffre pour stan_aov.
- brms sans compilateur : échec constaté. brms avec Rtools 4.4 (agent_06, mise au PATH du seul processus, g++ 13.2.0 trouvé, has_build_tools=TRUE) : **pas de résultat**, compilation de `brm(y~F)` encore en cours au bout de plus de 6 min, sur machine chargée, processus arrêté au budget. Ni succès ni échec établi.
- Comparaison anovaBF / brms+bridgesampling à a priori identiques : pas faite.
- BF multi-g à plus de 2 facteurs : pas de vérification indépendante.
- `anovaBF` n=500, 4 facteurs, et `whichModels="all"` à 3 facteurs : non atteints (agent_02 arrêté au budget). Fichier `agent_02` donc sans ligne FIN ni sessionInfo.



## sessionInfo (fin de agent_04_stan.R)

```
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26300)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8   
[3] LC_MONETARY=French_France.utf8 LC_NUMERIC=C                  
[5] LC_TIME=French_France.utf8    

time zone: Europe/Paris
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] brms_2.23.0            BayesFactor_0.9.12-4.8 Matrix_1.7-5          
[4] coda_0.19-4.1          rstanarm_2.32.2        Rcpp_1.1.2            

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1     dplyr_1.2.1          farver_2.1.2        
 [4] loo_2.10.1           S7_0.2.2             fastmap_1.2.0       
 [7] TH.data_1.1-5        tensorA_0.36.2.1     shinystan_2.7.0     
[10] promises_1.5.0       shinyjs_2.1.1        digest_0.6.39       
[13] estimability_2.0.0   mime_0.13            lifecycle_1.0.5     
[16] StanHeaders_2.39.1.1 processx_3.9.0       survival_3.8-6      
[19] magrittr_2.0.5       posterior_1.7.0      compiler_4.6.1      
[22] rlang_1.3.0          tools_4.6.1          igraph_2.3.4        
[25] bridgesampling_1.2-1 htmlwidgets_1.6.4    pkgbuild_1.4.8      
[28] plyr_1.8.9           RColorBrewer_1.1-3   dygraphs_1.1.1.6    
[31] multcomp_1.4-32      abind_1.4-8          miniUI_0.1.2        
[34] grid_4.6.1           stats4_4.6.1         xts_0.14.3          
[37] xtable_1.8-8         inline_0.3.21        ggplot2_4.0.3       
[40] emmeans_2.0.4        scales_1.4.0         gtools_3.9.5        
[43] MASS_7.3-65          mvtnorm_1.4-2        cli_3.6.6           
[46] reformulas_0.4.4     generics_0.1.4       otel_0.2.0          
[49] RcppParallel_6.2.1   reshape2_1.4.5       pbapply_1.7-5       
[52] minqa_1.2.8          rstan_2.32.7         stringr_1.6.0       
[55] shinythemes_1.2.0    splines_4.6.1        bayesplot_1.16.0    
[58] parallel_4.6.1       matrixStats_1.5.0    base64enc_0.1-6     
[61] vctrs_0.7.3          sandwich_3.1-3       boot_1.3-32         
[64] callr_3.8.0          crosstalk_1.2.2      glue_1.8.1          
[67] nloptr_2.2.1         codetools_0.2-20     DT_0.34.0           
[70] distributional_0.9.0 stringi_1.8.9        gtable_0.3.6        
[73] later_1.4.8          QuickJSR_1.11.0      lme4_2.1-0          
[76] tibble_3.3.1         colourpicker_1.3.0   pillar_1.11.1       
[79] htmltools_0.5.9      Brobdingnag_1.3-1    R6_2.6.1            
[82] Rdpack_2.6.6         shiny_1.14.0         lattice_0.22-9      
[85] markdown_2.1         rbibutils_2.4.1      backports_1.5.1     
[88] threejs_0.3.4        httpuv_1.6.17        MatrixModels_0.5-4  
[91] rstantools_2.7.1     gridExtra_2.3.1      nlme_3.1-169        
[94] checkmate_2.3.4      zoo_1.9-1            pkgconfig_2.0.3     
```
