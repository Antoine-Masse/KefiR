# KefiR 0.0.3.1

## Lighter installation

* Multivariate normality in the MANOVA (Mardia, Henze-Zirkler) is now computed inside KefiR,
  from the closed forms of the original papers. The p-values are those of the MVN package
  (checked on 9 datasets and in the tests). MVN is no longer installed with KefiR.
* `onewaytests` removed: it was imported but no longer called.
* 17 packages fewer to install (225 to 208), 9 of them with compiled code.

# KefiR 0.0.3.0

## New: `cor.bay()`, the Bayesian counterpart of `cor.test()`

`cor.bay()` takes the inputs of `cor.test()` (`x`, `y`, or `~ x + y` with `data`; `method`;
`alternative`; `conf.level`) and gives the posterior probability of a positive link, a
negative link and no link, with the Bayes factor, the strength of the link with a credible
interval, and a sentence ready for a report. No p-value. It tells the analysis in the steps
of `m.bay()` (plan, checks, method, starting point, solidity, strength, bilan), in French or
English, within 80 columns.

* Pearson: `BayesFactor::correlationBF()`, stretched beta prior (Ly et al. 2016).
* Kendall: the Bayesian tau of van Doorn et al. (2018), with the same prior carried over.
* Spearman on request: Pearson's Bayesian test on the ranks.
* `method = "auto"` takes Kendall when a value is far (beyond 3 interquartile ranges), a point
  has a Cook's distance above 1, or a variable has 7 distinct values or fewer; Pearson otherwise.
  The other method is computed as a cross-check.
* Priors set in advance, never on the data; sensitivity over the width of the prior and the
  starting odds.

Choices and measures: `dev/corbay/DECISIONS.md` and `dev/corbay/VALIDATION.md`.

## New: `m.bay()`, the Bayesian counterpart of `m.test()`

`m.bay()` takes the same inputs as `m.test()` and answers "how probable is each hypothesis
given these data?" instead of giving a p-value. It tells the analysis in the same seven steps
for every plan (plan, checks, shape of the noise, starting point, solidity, which groups
differ, bilan), in French or English, within 80 columns, and states its priors and how far
the answer moves when they change.

Plans covered:

* One factor, 2 or k groups: four laws of the noise compared and averaged (normal or Student,
  equal or unequal variances); which groups differ through the groupings of their means.
* Several crossed factors between groups: hereditary models with a HOP prior; normal with a
  common variance, normal with one variance per cell and Student averaged; post-hoc on the
  probable terms.
* Analysis of covariance, with parallel or different slopes.
* Paired data, repeated measures, split-plot designs and repeated measures with replicates.
* Several responses (MANOVA), one or several factors between groups.
* Yes/no, categories, counts and ordered ratings with one factor; yes/no, counts and ordered
  ratings with several factors between groups.

Refused with a message, until validated: random slopes, covariates with subjects, MANOVA with
covariates or subjects, discrete responses with subjects, paired yes/no.

Discrete responses, in detail:

* Yes/no and categories with one factor: Dirichlet-multinomial with neutral margins, shared
  with `g.bay()` (D30); it can differ from the several-factor engine by up to 0.15 on the same
  table, said in the help (D34).
* Counts and ordered ratings, with one factor or several, and their post-hocs: deviance Bayes
  factor with a Zellner-Siow g, r = sqrt(2)/4 (D31). On 0-10 ratings with a 0.6 SD gap, a wrong
  "no" falls from 27 % to 7 % of datasets, with the level under H0 unchanged.

Small samples:

* A group of one value is analysed with the common-variance laws only, says so at step 2, and
  its support is capped at weak; the raw BayesFactor error is no longer printed (D32).
* With very few values (3 + 3), step 3 says that the shape of the noise cannot be recognised
  instead of naming a law ahead (D33).
* "Strong" and "very strong" support now start at Bayes factors 20 and 150, as printed in the
  bilan, instead of 19 and 142 (D33).

Every choice and its measure by simulation are in `dev/mbay/DECISIONS.md` and
`dev/mbay/VALIDATION.md`. `m.test()` is unchanged.

## New: `g.bay()`, the Bayesian counterpart of the G test

`g.bay()` takes the inputs of `DescTools::GTest()` (`x`, `y`, `p`, `rescale.p`, or a formula
with `data`) and answers "how probable is a link given these counts?" instead of giving a
p-value. Same seven steps and same bilan as `m.bay()`, in French or English, within 80
columns.

* Goodness of fit to given proportions: Dirichlet against the proportions given, closed-form
  Bayes factor (Good 1967; the model of `multibridge`).
* Independence in an r x c table: Dirichlet-multinomial with neutral margins (closed form,
  the BDe score of Heckerman et al. 1995; a = 2, sensitivity a = 1, 2, 4): under independence
  the margins follow the law that a link gives them, so balanced margins say nothing of a link,
  and the answer is the same whatever the sampling plan. Chosen among four engines measured on
  1 000 tables per scenario: under independence, a strong link (P > 0.95) in 0.7 % of the
  tables at most (1.2 % with very unequal margins). The default `contingencyTableBF()` Bayes
  factors (up to 90 % of false strong links), the BIC of the log-linear model (independence
  strongly supported under a moderate link in up to 89 %) and `BAS::bas.glm` were set aside.
  The same function (`.bay_dm_neutral_bf()`) is meant for `m.bay()`.
* Coherence guard: independence is never "strongly supported" while some cells depart from it
  with a 95 % interval that excludes 1; the bilan says why.
* 2 x 2 tables: odds ratio with its 95 % credible interval.
* Which cells carry the link: ratio of each cell to its value without a link, with a 95 %
  credible interval (in place of the standardised residuals); Cramer's V a posteriori
  (Cohen's w for a goodness of fit).
* No minimum count and no Williams or Yates correction; small counts make the answer less
  decisive.

Decisions, references and measures by simulation: `dev/gbay/DECISIONS.md`,
`dev/gbay/PERIMETRE.md`, `dev/gbay/VALIDATION.md`.

## Corrections de `m.test()` (audit de novembre 2026)

Détail, mesures avant et après : `dev/mtest_audit/RAPPORT.md`.

### Résultats qui changent pour les étudiants

* **p-value globale d'un modèle à plusieurs termes** (`$global_pvalue`, et la valeur rendue avec `return = FALSE`) : c'était la p-value du premier facteur de la formule. C'est maintenant la plus petite p-value des k termes testés (effets principaux et interactions, corrigées de la sphéricité s'il y a lieu), corrigée de Bonferroni : `min(1, k × min(p))`. Le récit l'écrit, par exemple « p globale = 3 × 0.042 = 0.126 ». Avec un seul terme, rien ne change (`corrigraph`). La covariable d'une ANCOVA n'est pas comptée. Décision d'Antoine Massé, 12 novembre 2026.

* **Plan split-plot** `y ~ G*T + Error(id/T)` : l'interaction G:T n'est plus retirée sans message. Le facteur constant pour chaque sujet (G) est reconnu inter-sujet, et `$effects` reproduit `aov()`.
* **Post-hocs d'un split-plot** : ils passent par `emmeans` sur `lmer(y ~ G*T + (1 | id))`. Les tests appariés sur les cellules croisées apparaient des sujets différents.
* **Variances hétérogènes, données non normales** : la confirmation du test de Levene se fait maintenant par Brown-Forsythe **sur les variances** (`car::leveneTest(center = median)`), et non plus par `onewaytests::bf.test()`, qui compare des moyennes. Une hétéroscédasticité nette n'est plus écartée à tort, et l'ANOVA classique n'est plus retenue dans ce cas.
* **Données non normales, plus de deux groupes, variances hétérogènes confirmées** (asymétrie modérée) : `m.test` applique le F de Welch au lieu de Kruskal-Wallis, qui suppose des variances égales. Sous H0 (n = 5/5/40, sd = 4/4/1), les faux positifs passent de 38,5 % à 10,5 %.
* **Mesures répétées (`Error(id/T)`, 3 temps ou plus)** : la p-value des termes intra-sujet est **toujours corrigée de la sphéricité**, Greenhouse-Geisser si son epsilon est ≤ 0,75, Huynh-Feldt sinon (Girden, 1992). Le test de Mauchly reste affiché, à titre d'information : peu puissant avec peu de sujets, il ne décide plus de la correction. La p-value affichée, celle de `$effects` et la p-value globale sont les p-values corrigées.
  Les p-values des étudiants changent : en général elles augmentent (correction plus prudente), sauf quand F < 1. La correction s'applique au F du tableau `aov()` affiché.
* **Mesures répétées et plans en tranches (`Error(id/T)`) : normalité des résidus.** Elle est testée sur les résidus de la strate intra-sujet (une valeur par observation), et non plus sur un mélange des « résidus » de toutes les strates, d'échelles différentes. Avec 12 sujets × 6 temps et des données normales, plus de la moitié des jeux partaient à tort vers le test de Friedman (« non-normalité sévère »). Ces jeux passent maintenant par l'ANOVA à mesures répétées corrigée de la sphéricité : **la méthode et la p-value changent pour eux.** Des données réellement non normales vont toujours vers Friedman.
* **Voie Friedman** : `$method` vaut `"Friedman test"` (et non plus `"Parametric ANOVA/t-test"`). Les post-hocs de Nemenyi apparient chaque sujet à lui-même ; ils échouaient ou appariaient des sujets différents quand les lignes étaient rangées par sujet.
* **Plan en tranches aux groupes inégaux** (`y ~ grp*temps + Error(id/temps)`, 5 + 9 sujets) : l'analyse aboutit (ANOVA à mesures répétées) au lieu de s'arrêter sur « Les variables explicatives de la formule ne correspondent pas aux colonnes de g ».
* Le récit du F de Welch ne le dit plus « robuste » : avec de petits groupes asymétriques, son risque de faux positif reste au-dessus d'alpha.
* **`y ~ T + (1|id)`** : l'analyse et les post-hocs sont appariés, comme avec `Error(id/T)`.
* **Mesures répétées, lignes non triées par sujet** : les données sont triées par sujet, et les post-hocs appariés comparent bien chaque sujet à lui-même.
* **MANOVA non normale** : la voie robuste utilise le lambda de Wilks sur les rangs (`rrcov::Wilks.test(method = "rank")`). L'appel prend moins d'une seconde au lieu de plus de 6 minutes. `$method` vaut `"Robust MANOVA (Wilks, ranks)"`.

### Plantages et refus corrigés

* `verbose = FALSE` sur un modèle mixte : plus d'erreur « $ operator is invalid for atomic vectors ».
* `m.test(formula = y ~ F, data = d)` est accepté.
* `m.test(x, g, id = id)` avec plus de deux modalités ne plante plus.

### Affichage

* Hors locale française, les messages s'affichent en anglais. Jusqu'ici, ils restaient toujours en français.
* En modèle mixte, les lettres d'une interaction nomment la cellule entière (« B t3 ») au lieu du seul premier facteur.
