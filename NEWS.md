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
