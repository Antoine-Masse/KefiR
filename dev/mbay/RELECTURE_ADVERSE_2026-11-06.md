# Relecture adverse du 6 novembre 2026 : modèle à variances par groupe et zs

Objet : `R/sys_bay_oneway.R`, `tests/testthat/test-bay-oneway.R`, `dev/mbay/PROJET_METHODO_2026-11-06.md`. Scripts et sorties : `dev/mbay/scripts/relecture/agent_rev2_*`. Rien n'a été modifié dans `R/` ni `tests/`. Relecture arrêtée sur demande à 14:30 : ce qui n'est pas fini est noté comme tel en fin de document.

## Verdict

Les constantes et le calcul sont justes. **Le modèle a un défaut de fond, ni numérique ni propre à τ arithmétique : un petit groupe très éloigné fait pencher le BF vers H0, d'autant plus qu'il est loin.** Je ne le recommande pas comme défaut à variances inégales sans correction de la loi a priori. Il supprime bien les faux « Oui » du plan B4 à n_i ≥ 3.

## BLOQUANT

### B1. Plus un petit groupe s'éloigne, plus le BF soutient H0
`R/sys_bay_oneway.R:312-332` (g unique partagé), `:443-446` (modèle) ; `PROJET_METHODO:152`.

Exemple brut (`agent_rev2_08_complements.R`) :

```r
set.seed(7); y <- c(rnorm(40), rnorm(3, 50, 1))
g <- rep(c("a","b","c","d","e"), c(10,10,10,10,3))
.bay_bf_oneway_hetero(y, g)$log_bf   # -6.92, soit environ 1000 contre 1 pour H0 ; ESS 1879, erreur 2,3 %, aucun avis
.bay_bf_oneway_zs(y, g)              # +85.6
oneway.test(y ~ factor(g))$p.value   # 3e-15
```

**Loi d'échelle** [main, confirmée par le tableau]. Groupe de taille n_s à distance d, SS fixés. Sous H0, μ reste avec les autres groupes et σ_s absorbe l'écart : m(y|H0) ∝ d^(-n_s). Sous H1, le g unique doit couvrir d, ce qui gonfle tous les contrastes : m(y|H1) ∝ d^(-k). Donc **BF ∝ d^(n_s − k)**.

| Plan | Pente prévue | d = 10 | d = 100 | d = 1000 | d = 10⁴ | zs, d = 10⁴ |
|---|---:|---:|---:|---:|---:|---:|
| (30, 30, 5) | +2 | 2,63 | 7,19 | 11,79 | 16,40 | 476 |
| (20×4, 5) | 0 | −3,44 | −3,65 | −3,65 | −3,65 | 604 |
| (30, 30, 2) | −1 | −2,41 | −4,46 | −6,68 | −8,81 | 428 |
| (20×4, 3) | −2 | −5,99 | −11,79 | −16,91 | −21,52 | 570 |

(log BF, 4000 tirages au départ, `agent_rev2_04.out`. Pour (30, 30, 2), l'ESS tombe à 582 et 70 à d = 1000 et 10⁴ : ces deux cases sont moins précises.)

- **Pas un artefact de tirage** : quadrature 3D, (30, 30, 2), d = 100 : −4,460 contre −4,436 pour l'IS.
- **Plans équilibrés, n = 5** (d = 10, 100, 1000) : k = 4 : 1,3 / 3,6 / 5,9. k = 6 : −3,1 / −5,7 / −8,0. k = 8 : −5,8 / −14,6 / −21,5.
- **Fréquence** : (10, 10, 10, 10, 3), sd 1, écart réel de 3σ sur le dernier groupe, 150 jeux : BF < 1/3 dans **67 %**, BF > 3 dans 5 %. Welch p < 0,05 dans 83 %, zs BF > 3 dans 82 %.
- **La moyenne géométrique ne corrige pas** : (20×4, 3), d = 10⁴ : −21,73 ; (30, 30, 2) : −5,10.
- **À n_s = 2, k = 2, le BF ne croît presque plus** : (2, 2) à d = 10⁶ : référence 2,11 (BF ≈ 8) ; (1000, 2) : 2,59 (BF ≈ 13, ESS 20, imprécis). La consistance d'information, « attendue » au §A5, est fausse dès n_s ≤ k.

**Correction, à arbitrer par Antoine car elle change le modèle.**
1. g par groupe (a_i ~ N(0, g_i τ²)) : coût d^(-2) pour le seul groupe éloigné, soit BF ∝ d^(n_s − 2) [main]. Non programmé, et perd la réduction à JZS.
2. À défaut, refuser ou signaler tout plan avec min(n_i) ≤ k, et tester la croissance de `log_bf` en d.

### B2. Un groupe dispersé sans rapport fait basculer le BF, sans plancher
`R/sys_bay_oneway.R:316` (`tau2 <- rowMeans(s2)`).

Paire 1-2 à t ≈ 3,5, n = 20, groupe 3 étalé d'un facteur S autour de sa moyenne (`agent_rev2_04.out`) :

| S | 1 | 10 | 100 | 1000 | 10⁴ | 10⁶ |
|---|---:|---:|---:|---:|---:|---:|
| BF, τ arithmétique | 7,22 | 5,38 | **0,60** | 0,060 | 0,0060 | 6e-5 |
| log BF, τ géométrique | 1,98 | 2,68 | 2,42 | 1,76 | 1,01 | −0,52 |

Groupes 1 et 2 seuls : log BF = 2,98 (BF ≈ 20). BF ∝ 1/S à grand S, ESS stable à 3600 : aucun signal. La géométrique dérive aussi, plus lentement.

Fréquence, (30, 30, 30), écart 0,8 entre deux groupes, 150 jeux : sd (1, 1, 1) : BF > 3 dans 73 %, BF < 1/3 dans 4 %. sd (1, 1, 20) : 31 % et 27 %.

**Correction** : même arbitrage que B1. D'ici là, ne pas présenter ce BF comme un test des paires à k ≥ 3.

## IMPORTANT

### I1. L'ESS et error_pct peuvent mentir (r petit, n_i = 2)
`R/sys_bay_oneway.R:294-300`, `:336-343`, `:383-385`, `:447-448`.

Plan `st_of(c(2,30), c(0,1e4), c(1e-4,29))` (`agent_rev2_lib.R`), contre quadrature (`agent_rev2_05.out`, `agent_rev2_08.out`) :

| r | Tirages | log BF | Référence | Erreur annoncée | ESS |
|---|---:|---:|---:|---:|---:|
| 0,05 | 20 000 | 2,010 | 2,170 | 0,54 % | 12 668, sans avis (écart de 30 erreurs annoncées) |
| 0,005 | 2 000 | **−0,287** | **0,574** | 0,95 % | 1 694, sans avis (mauvais côté de 1) |
| 0,005 | 20 000 | −0,256 | 0,574 | 3,28 % | 887 |
| 0,005 | 200 000 | +1,332 | 0,574 | 79,8 % | 2 |

(r = 0,005 : `max_grow = 0`.) La masse manquée serait à s² ≈ n d², dans la queue de q [main]. **À r = 1/2, aucun biais non signalé trouvé** : (2,2) et (3,3), d = 10³, concordent (|z| ≤ 1,4). À (2,2), d = 10⁶, le signal LOW sort bien, mais l'erreur est sous-estimée : 6,8 % annoncés pour 39 % réels (seed 2, S = 4000).

- `rel_se` et ESS sont une seule information : rel_se² ≈ (S/ESS − 1)/S.
- Le commentaire `:447-448` (« weights are bounded ») ne garantit rien d'utile : la borne dépend des données et peut être énorme.

**Correction** : proposition en mélange défensif, avec une composante IG(ν_i/2, (SS_i + n_i d_i²)/2) par groupe ; diagnostic de queue ; ancrage par quadrature à k ≤ 3. Non programmé.

### I2. n_i = 2 : le faux « Oui » est déplacé, pas supprimé
`R/sys_bay_oneway.R:424` ; `PROJET_METHODO:22, 69, 138`. Sous H0, 150 jeux (`agent_rev2_06.out`) :

| Plan | sd | BF > 3 | BF > 19 (P > 0,95) | E[BF] | zs BF > 19 |
|---|---|---:|---:|---:|---:|
| (5, 5, 40) | 4, 4, 1 | 5,3 % | 0,7 % (1/150) | 0,97 | 16,7 % |
| (3, 3, 40) | 4, 4, 1 | 7,3 % | 0 % | 1,25 | 39 % |
| **(2, 2, 40)** | 4, 4, 1 | **18 %** | **6,0 % (9/150)** | **6,5** | 41 % |
| **(2×6)** | 1 | 12,7 % | 1,3 % | **6,5** | 0,7 % |

E[BF] est tiré par quelques valeurs extrêmes : indicatif à 150 jeux. Le gain sur la voie gaussienne reste net à (5, 5, 40) et ne s'étend pas à n_i = 2. Sur (2×6), à variances égales, ce modèle fait pire que zs.

Quasi-égalités, groupe n = 2 à 1,5σ d'un groupe n = 50 : eps = 1 : 0,52 ; 10⁻³ : 6,81 ; 10⁻⁸ : 18,32 (log BF). À eps = 0, refus.

**Correction** : refuser n_i = 2 par défaut, ou avis explicite ; ajouter ces plans à la simulation du projet.

### I3. Grille en g tronquée à gauche à grand k, sans signal
`R/sys_bay_oneway.R:309, 319`. Le commentaire borne le facteur a priori, pas l'intégrande ; le pic est en g ≈ r²/(k−1) [main].

Mêmes tirages, bornes élargies, n_i = 1000 (`agent_rev2_04c.out`) : k ≤ 100 : au plus 0,001. k = 200, r = 1/2 : 0,048. **k = 200, r = 1 : 4,09. k = 500, r = 1 : 71,1.**

**Correction** : borne gauche dépendant de k (`log(r^2 / (120 * k))`, par exemple), et test à k = 200.

### I4. Intégrale de H0 : blocage à très grand n
`R/sys_bay_oneway.R:266-279`. n = (10⁶, 10⁶, 10⁶), moyennes égales, SS ∝ (1, 4, 9) : **plus de 40 s**, arrêt par limite (`agent_rev2_07.out`). À 8·10⁵ : 0,22 s. Cause **probable, non démontrée** : `tol / length(a)` passe sous le bruit d'arrondi de `log_f`.

**Correction** : plancher de tolérance relatif à |log f|, plafond du nombre d'intervalles, test à n = 10⁶.

### I5. H0 : pic plus fin que l'arrondi de m
`R/sys_bay_oneway.R:70, 240`. Contre référence en coordonnées locales (`agent_rev2_03.out`) : écart 1,5·10⁻⁹ à 10¹⁵ σ, **+2,61 à 10¹⁸ σ**, +0,22 pour (2,2) à 10¹⁶. H0 surestimée, donc BF sous-estimé.

Cas brut : `y = c(0, 1e-160, 2e-160, 1, 1+1e-10, 1)`, deux groupes : **log BF = −308**, erreur 2,2 %, aucun avis (`agent_rev2_02_entrees.R`). Le lien de ce cas avec ce mécanisme est probable, non vérifié.

**Correction** : intégrer en u = μ − m_i par fenêtre, comme `ref_h0` dans `agent_rev2_03_h0.R`.

### I6. zs : Inf sur des données constantes à l'arrondi près
`R/sys_bay_oneway.R:66-69, 118-128`. Le seuil de bruit s'applique à SS, pas à SSB.

```r
.bay_bf_oneway_zs(c(.3,.3,.3,.1+.2,.1+.2,.1+.2), rep(1:2, each=3))  # Inf, "variance residuelle nulle"
```

Même effet à y = 10⁹ + N(0, 2·10⁻⁷) sous H0. Voisin : niveaux exacts + 1 ulp : Inf ; + 2 ulp : log BF = 274, sans signal.

**Correction** : traiter SSB sous le bruit comme nul, donc refus « constant » ; signaler un SSE/SST au niveau de l'arrondi.

### I7. zs : Inf à N = k + 1
`R/sys_bay_oneway.R:113, 123-128`. y = (1, 1 + eps, 5), groupes (1, 1, 2) : eps = 10⁻¹⁵ : **2,85** ; eps = 0 : **Inf** (avec avis). La divergence n'est que logarithmique [main]. L'argument du commentaire `:106` contre N = k vaut ici aussi.

**Correction** : refuser N = k + 1 avec SSE = 0.

### I8. Puissance sur le plan B4, « non mesurée » au `PROJET_METHODO:138`
n = (5, 5, 40), sd (4, 4, 1), 150 jeux, effet sur le groupe 1 :

| Effet | BF > 3 | BF > 19 | Welch p < 0,05 | Welch p < 0,005 |
|---:|---:|---:|---:|---:|
| 4 | 18 % | 4 % | 23 % | 5 % |
| 8 | 51 % | 13 % | 82 % | 29 % |

Prudence nette à effet fort. « Lève le blocage » doit s'accompagner de ces chiffres.

### I9. Tests faibles (lecture seule : suite non relancée, pas de mutation)
`tests/testthat/test-bay-oneway.R`.

- Aucune référence indépendante pour H1 à k ≥ 3. `:169` couvre k = 2 à σ fixé ; `:214` compare l'estimateur à lui-même et ne voit ni I1 ni I3.
- `:154` (H0 séparée) et `:288` (F très grand) : presque seulement `is.finite`.
- `:281` : un seul plan k = 2 : aucune monotonie en d, aucun petit n_s, donc rien sur B1.
- `:270` : proximité avec zs à 0,5 près, n = 300 seulement.
- `:224` : un seul cas, borne gauche non testée.
- `:162` : r → 0 à F modéré seulement.
- Absents : k grand, n = 10⁶, g avec NaN, y quasi constant, N = k + 1, `n_draws = Inf`.

## MINEUR

Preuves dans `agent_rev2_02_entrees.R`, sauf M7 (`agent_rev2_01`).

| # | Où | Défaut | Preuve | Correction |
|---|---|---|---|---|
| M1 | `:41-42` | NaN et Inf dans g forment un groupe, sans avis | `g = c(1,1,1,2,2,2,NaN,NaN)` : 1,257 contre 1,239 sans ces deux points | Tester `is.na(g)` avant `as.character`, refuser le non fini |
| M2 | `:36-41` | g liste ou matrice accepté ; g numérique fusionné | `0.1+0.2` et `0.3` confondus | `is.atomic`, `is.null(dim(g))` |
| M3 | `:60-65` | Sous-dépassement : groupe sd = 1 déclaré de variance nulle | valeur aberrante 10²⁰⁰ | SS à l'échelle de chaque groupe |
| M4 | `:151` | SSE sous-normale : erreur brute | `c(0,1e-160,2e-160,1,1,1)` : « 'to' doit être un nombre fini » | Refus propre |
| M5 | `:75, :150, :319` | r extrême : erreur brute | r = 10⁻³⁰⁰, 10³⁰⁰ | Borner r |
| M6 | `:454-463` | n_draws = Inf ou 10¹⁰ : erreur brute ; 150,7 et seed = 1,9 tronqués | seed 1,9 identique à seed 1 | Vérifier entier et borne |
| M7 | `:296-298` | Résultat dépendant des étiquettes, dans l'erreur MC | −0,5261 / −0,5265 | Trier les groupes avant tirage |
| M8 | `:9-20` | Références « a verifier », DOI manquants (BP-008) | lecture | Compléter |

Temps, 20 000 tirages : k = 3 : 1,5 s ; k = 30 : 6 à 10 s ; k = 60 : 48 à 56 s.

Règles maison : aucun caractère non ASCII dans les deux fichiers, messages bilingues, pas de `print` ni de `cat`.

## Projet méthodologique : attaques de fond

1. **§A1.4 « lève le blocage B4 ».** Vrai sur ce plan (zs reproduit à 16,7 %), mais voir I2 pour n_i = 2 et I8 pour la puissance. Hypothèse non démontrée : une part du gain viendrait du même mécanisme que B1, σ_i absorbant l'écart.
2. **§A4 et §A5, « géométrique meilleure candidate ».** Vue jusqu'à S = 100 seulement ; dérive au-delà et ne résout pas B1. Le choix de τ ne suffit pas.
3. **§A5, consistance d'information.** Fausse après intégration dès n_s ≤ k (B1). L'argument d'orthogonalité du §A3 est asymptotique en n, alors que n_s est ici fixe et petit.
4. **§A3, « convergence si n_i ≥ 2 ».** Exact, mais pas un critère d'emploi.
5. **§4, ROPE sur écart/σ** [main]. Non transposé aux variances inégales : avec l'échelle τ, une paire peut devenir « négligeable » parce qu'un autre groupe est dispersé (B2 se reporte sur le ROPE).
6. **Reproductibilité.** Les chiffres [exécuté] renvoient à un dossier temporaire (`:6`).

Refaits sans objection : tables a priori à 2 et 3 facteurs (1/2, 1/8, 3/7, 3:4, 1/128, 1/4096, 1/19, 14/19), 0,089 pour δ et la pente, cohérence ±0,5, Sellke 0,41 et 0,29 [main] ; plafonds 0,720, 0,859, 0,909 [exécuté].

## Essayé de casser, sans succès

- **Constantes, marginale après μ plat, remplacement de (I − J/k) par I** : refaits à la main, justes. Contrôle indépendant par ratio sous postérieur H0 (Gibbs) : k = 2 : 0,1025 (IS), 0,1045 (Gibbs), 0,1028 (quadrature) ; k = 3 : −0,7446 et −0,7435. Médiane de |δ| : 0,7077 pour √2/2.
- **Rater un pic de H0** : 3 pics jusqu'à 10¹⁵ σ, sd 10⁻⁶/10⁶, n_i = 10⁶ avec écart 1 : écart ≤ 1,5·10⁻⁹. Entre deux groupes éloignés, pas de pic intermédiaire [main].
- **Consistance à k fixe** : H0, n = 10² à 10⁵, pente −log n comme attendu ; H1, d = 0,3 : log BF 902 à n = 10⁵ ; n = (10, N) sous H0 : stable à −1,18.
- **Pas de grille en g** : écart h contre h/4 ≤ 1,7·10⁻⁷ pour k ≤ 60, à F nul ou très grand.
- **Invariances** : origine et échelle ; y × 10⁻³⁰⁰ exact pour zs.
- **Seuil de bruit à y ≈ 10⁹** : refuse seulement sous sd ≈ 3·10⁻⁷, soit environ 2 ulp ; sd = 10⁻⁶ reste traité. Le défaut est SSB (I6).
- **RNG** : état préservé, pas de `.Random.seed` créé, reproductible.
- **Calibration sans n_i = 2** : E[BF] de 0,97 et 1,25 sur (5,5,40) et (3,3,40) ; P(BF > 19) ≤ 0,7 %.

## Non fini

- Plan I1 à r = 0,1 et r = 1/2, et à S = 800 000 : non mesurés.
- H0 : « N 10⁶ chevauchants », k = 30, 200, 500, boucle aléatoire de 40 plans et balayage de sd_A : calculs interrompus (`agent_rev2_03.out` s'arrête avant).
- Consistance à n = 10⁶ et 10⁷ : interrompue, sans doute par I4.
- Suite de tests non relancée, aucune mutation.
- Correction par g multiples et proposition en mélange : non programmées.
- Non refaits : comparaison de zs à `anovaBF`, sources du §A2.
- Simulations à 150 jeux seulement (9/150 : environ 3 à 11 %). `agent_rev2_04_versionA.out` vient d'une première version du script (t ≈ 9,5 : BF de 4,7·10⁷ à 502).
