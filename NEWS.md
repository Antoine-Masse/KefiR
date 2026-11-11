# KefiR (development version, branch feat/m.bay)

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

Every choice and its measure by simulation are in `dev/mbay/DECISIONS.md` and
`dev/mbay/VALIDATION.md`. `m.test()` is unchanged.

## New: `g.bay()`, the Bayesian counterpart of the G test

`g.bay()` takes the inputs of `DescTools::GTest()` (`x`, `y`, `p`, `rescale.p`, or a formula
with `data`) and answers "how probable is a link given these counts?" instead of giving a
p-value. Same seven steps and same bilan as `m.bay()`, in French or English, within 80
columns.

* Goodness of fit to given proportions: Dirichlet against the proportions given, closed-form
  Bayes factor (Good 1967; the model of `multibridge`).
* Independence in an r x c table: `BayesFactor::contingencyTableBF()` (Gunel and Dickey 1974;
  Jamil et al. 2017). The sampling plan is chosen from the form of the input (joint
  multinomial for a table or two variables, groups fixed for `y ~ g`, as in `m.bay()`) and can
  be set with `sampling =`; the answer under every other plan is shown, and called not settled
  if a plan reverses it.
* 2 x 2 tables: odds ratio with its 95 % credible interval.
* Which cells carry the link: ratio of each cell to its value without a link, with a 95 %
  credible interval (in place of the standardised residuals); Cramer's V a posteriori
  (Cohen's w for a goodness of fit).
* No minimum count and no Williams or Yates correction: the computation is exact.

Decisions, references and measures by simulation: `dev/gbay/DECISIONS.md`,
`dev/gbay/PERIMETRE.md`, `dev/gbay/VALIDATION.md`.
