# sys_bay_dm_neutral.R
# Bayes factor of a link in an r x c table of counts: Dirichlet-multinomial with neutral margins,
# shared by g.bay() and m.bay() (categorical response with one factor). Self-contained: it uses
# nothing else from KefiR, so that both functions call the same code (dev/gbay/DECISIONS.md, G12).
#
#   H1: the r x c cell probabilities Dirichlet(a, ..., a).
#   H0: rows and columns independent, row margin Dirichlet(c a), column margin Dirichlet(r a):
#       the laws that H1 gives to the margins, so that the margins carry no evidence of a link.
#   By the aggregation property of the Dirichlet, the row margin cancels:
#       log BF10 = sum_i log DM(row_i; a) - log DM(column totals; r a),
#   the independent-multinomial Bayes factor of Gunel and Dickey (1974) with the common H0
#   proportions Dirichlet(r a) (instead of BayesFactor's Dirichlet(r (a - 1) + 1)). The same
#   number with the rows fixed, the columns fixed (it is symmetric) or only N fixed.
#   This is the BDe score of Bayesian networks, X -> Y against the empty graph (Heckerman,
#   Geiger and Chickering 1995), with an equivalent sample size r c a; uniform case: BDeu
#   (Buntine 1991).
#
# Default a = 2 (G12): prior median of Cramer's V about 0.2 to 0.4 depending on the table, a
# "medium" effect, as r = 1/2 in m.bay and kappa = 1/3 in cor.bay; sensitivity a = 1, 2, 4.
#
# References.
#   Gunel, E., & Dickey, J. M. (1974). Bayes factors for independence in contingency tables.
#     Biometrika, 61(3), 545-557.
#   Heckerman, D., Geiger, D., & Chickering, D. M. (1995). Learning Bayesian networks: the
#     combination of knowledge and statistical data. Machine Learning, 20, 197-243.
#   Buntine, W. (1991). Theory refinement on Bayesian networks. Proceedings of UAI 7, 52-60.

# Default concentration and sensitivity grid.
.bay_dm_neutral_a <- function() 2
.bay_dm_neutral_grid <- function() c(1, 2, 4)

# Log Bayes factor of a link (H1 against independence) for a table of counts with at least two
# non-empty rows and columns. tab: matrix or table of counts; a: concentration per cell.
.bay_dm_neutral_bf <- function(tab, a = .bay_dm_neutral_a()) {
	tab <- unclass(as.matrix(tab))
	ldm <- function(n, al) lgamma(sum(al)) - lgamma(sum(al) + sum(n)) + sum(lgamma(al + n) - lgamma(al))
	sum(apply(tab, 1L, ldm, al = rep(a, ncol(tab)))) - ldm(colSums(tab), rep(nrow(tab) * a, ncol(tab)))
}
