# KefiR 0.0.2.2 (development version)

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

Every choice and its measure by simulation are in `dev/mbay/DECISIONS.md` and
`dev/mbay/VALIDATION.md`. `m.test()` is unchanged.
