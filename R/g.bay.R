#' Bayesian counterpart of the G test: link between qualitative variables, from counts
#'
#' @description
#' `g.bay()` takes the inputs of [DescTools::GTest()] and answers the opposite question: not
#' "how surprising are these counts if there is no link?", but "how probable is a link given
#' these counts?". It covers the same two cases:
#' - **goodness of fit** of one variable to given proportions `p` (equal by default);
#' - **independence** of two variables in an r x c contingency table.
#'
#' The way is chosen alone, from the form of the input, and told step by step as [m.test()]
#' and [m.bay()] do. The bilan ranks the two hypotheses (link or independence; departure or
#' fit) by posterior probability, with a plain qualifier and the Bayes factor, never a p-value,
#' then says where the link lies (which cells), how large it is, and gives a sentence ready for
#' a report.
#'
#' @section Method (dev/gbay/DECISIONS.md, G1 to G11):
#' - **Goodness of fit**: multinomial counts; under "no departure" the proportions are exactly
#'   `p`, under "departure" they follow a uniform Dirichlet. The Bayes factor has a closed form
#'   (Good 1967; the model of `multibridge`, Sarafoglou et al. 2023).
#' - **Independence**: test of the interaction of a log-linear Poisson model, the margins
#'   free under both hypotheses. The Bayes factor is -BIC/2 with n = N, the number of
#'   individuals (unit-information prior; Raftery 1986, 1995), the same method as [m.bay()] for
#'   discrete responses (D25). It is read on the G statistic itself:
#'   log BF = (G - df log N) / 2. The answer does not depend on how the data were collected
#'   (total fixed, one margin fixed, or nothing fixed). The default Bayes factors of
#'   [BayesFactor::contingencyTableBF()] are not used: under independence they announced a
#'   strong link (Bayes factor above 20) in up to 21 % (joint multinomial) and 90 % (Poisson)
#'   of 4 x 5 tables with N = 30 (dev/gbay/VALIDATION.md, section 9).
#' - **2 x 2 tables**: odds ratio with its 95 % credible interval.
#' - **Cells**: for each cell, how many times it is more (or less) frequent than without a
#'   link, with a 95 % credible interval: the Bayesian counterpart of the standardised
#'   residuals. A cell is read when its interval excludes 1.
#' - **Size**: Cramer's V a posteriori (Cohen's w for a goodness of fit). It is positive even
#'   without a link: a size, not a proof.
#' - **Small counts**: no minimum count and no correction. Small counts make the answer less
#'   decisive: the Bayes factor stays near 1 or leans towards no link.
#' - **Power**: lower than the G test at p < 0.05, in exchange for almost no false "yes"
#'   (under independence, a Bayes factor above 20 in 0.2 % of the tables at most). In large
#'   tables the BIC is very cautious: in a 3 x 4 table with N = 200, a medium link (w = 0.3)
#'   is found in 16 % of the cases, against 93 % for the G test. Measures in
#'   dev/gbay/VALIDATION.md.
#' - **Sensitivity**: the starting chances (1 against 2 to 2 against 1); for a goodness of fit,
#'   also the concentration of the Dirichlet (1/2, 1, 2).
#'
#' Not handled: tables with three dimensions or more; paired data (McNemar); which groups
#' differ two by two (use [m.bay()] with `y ~ g`).
#'
#' @param x Counts: a vector (goodness of fit), a matrix, a table or a data.frame of counts
#'   (independence); or a factor; or a formula (`~ A`, `~ A + B`, `y ~ g`, `n ~ A + B` with a
#'   column of counts `n`, one row per combination); or a column name of `data`.
#' @param y Second variable, as in [DescTools::GTest()] (a vector or a column name of `data`).
#' @param p Proportions expected for a goodness of fit (equal by default).
#' @param rescale.p If `TRUE`, `p` is rescaled to sum to 1, as in [DescTools::GTest()].
#' @param data Data.frame, for a formula or column names.
#' @param correct Accepted for compatibility with [DescTools::GTest()], ignored: the
#'   computation is exact and needs no Williams or Yates correction.
#' @param prior_h1 Prior probability of a link (or of a departure). Default 1/2.
#' @param sensitivity If `TRUE` (default), the sensitivity to the priors is computed.
#' @param draws Posterior draws for the estimation (cells, V, odds ratio). Default 4000.
#' @param seed Seed of this call. The global random state is left untouched.
#' @param verbose If `TRUE`, tells the analysis in seven steps as [m.bay()] does: "==>" for a
#'   finding, "-->" for a decision, "Note:" for anything else, every line within 80 columns
#'   (or `getOption("width")` if smaller); then the bilan. `print(result, story = TRUE)`
#'   replays the steps, `print(result, detail = TRUE)` adds the method.
#' @param return If `TRUE`, the full object (invisibly when the bilan was printed); if
#'   `FALSE`, one number, the posterior probability of a link (or of a departure), computed
#'   without the posterior draws. This is not a p-value.
#' @param detail If `TRUE`, the bilan also gives the method in detail.
#'
#' @return With `return = TRUE`, a list of class `gbay`:
#'   \describe{
#'     \item{hypotheses}{The two hypotheses ranked: `label_en`, `label_fr`, `posterior`,
#'       `prior`, `null`.}
#'     \item{prob_effect, bf10, log_bf}{Posterior probability of a link (or departure), Bayes
#'       factor and its log.}
#'     \item{test}{`"independence"` or `"gof"`.}
#'     \item{G, df}{The G statistic and its degrees of freedom (independence: the Bayes factor
#'       is read on them).}
#'     \item{table, expected, p}{Counts, counts expected without a link (or under `p`), `p`.}
#'     \item{cells}{Per cell (or category): observed, expected, `ratio` to its value without a
#'       link (or to `p`) with `ratio_low`, `ratio_high` (95 % HDI), `p_more` (probability that
#'       the ratio exceeds 1) and `reading` (`"more"`, `"less"` or `"-"`).}
#'     \item{cramer_v, w, odds_ratio}{Size of the effect: median and 95 % HDI; the odds ratio
#'       (2 x 2) with `p_gt1`.}
#'     \item{sensitivity}{Range of the leading answer over the priors, and verdict.}
#'     \item{story, notes, method, priors, call, n, label}{Description of the analysis; `story`
#'       holds the steps in English and French.}
#'   }
#'
#' @references
#' Raftery, A. E. (1986). Choosing models for cross-classifications. American Sociological
#' Review, 51(1), 145-146.
#'
#' Raftery, A. E. (1995). Bayesian model selection in social research. Sociological
#' Methodology, 25, 111-163.
#'
#' Good, I. J. (1967). A Bayesian significance test for multinomial distributions. Journal of
#' the Royal Statistical Society B, 29(3), 399-431.
#'
#' Jamil, T., Ly, A., Morey, R. D., Love, J., Marsman, M., & Wagenmakers, E.-J. (2017).
#' Default "Gunel and Dickey" Bayes factors for contingency tables. Behavior Research Methods,
#' 49, 638-652 (the default Bayes factors not used here, see Method).
#'
#' Sarafoglou, A., Haaf, J. M., Ly, A., Gronau, Q. F., Wagenmakers, E.-J., & Marsman, M.
#' (2023). Evaluating multinomial order restrictions with bridge sampling. Psychological
#' Methods, 28(2), 322-338. \doi{10.1037/met0000411}
#'
#' Albert, J. H. (1997). Bayesian testing and estimation of association in a two-way
#' contingency table. Journal of the American Statistical Association, 92(438), 685-693.
#'
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. Journal of the American Statistical
#' Association, 90(430), 773-795.
#'
#' @seealso [DescTools::GTest()], [m.bay()], [BayesFactor::contingencyTableBF()]
#'
#' @examples
#' # Goodness of fit: 3:1 Mendelian ratio (the counts of GTest's help page).
#' g.bay(c(705, 224), p = c(3, 1) / 4, verbose = FALSE)$prob_effect
#'
#' # Independence in a 2 x 3 table (a table of counts).
#' tab <- matrix(c(20, 15, 25, 30, 10, 12), nrow = 2,
#'               dimnames = list(sex = c("F", "M"), colour = c("blue", "green", "red")))
#' r <- g.bay(tab)
#' r$cells[, c("row", "col", "ratio", "ratio_low", "ratio_high")]
#'
#' # 2 x 2: odds ratio with its credible interval.
#' g.bay(matrix(c(18, 7, 9, 16), 2), verbose = FALSE)$odds_ratio
#'
#' # Two factors in a data.frame, or a response and a factor.
#' d <- data.frame(treat = rep(c("cold", "warm"), each = 30))
#' set.seed(1)
#' d$germ <- ifelse(runif(60) < ifelse(d$treat == "warm", 0.8, 0.5), "yes", "no")
#' g.bay(germ ~ treat, data = d, verbose = FALSE)$prob_effect
#' g.bay(d$treat, d$germ, verbose = FALSE)$prob_effect
#'
#' # A data.frame of counts, one row per combination.
#' k <- data.frame(eye = c("blue", "blue", "brown", "brown"), hair = c("fair", "dark", "fair", "dark"),
#'                 n = c(32, 11, 10, 41))
#' g.bay(n ~ eye + hair, data = k, verbose = FALSE)$prob_effect
#'
#' @importFrom BayesFactor contingencyTableBF extractBF
#' @importFrom bayestestR hdi
#' @export
g.bay <- function(x, y = NULL, p = NULL, rescale.p = FALSE, data = NULL, correct = NULL,
                  prior_h1 = 1/2, sensitivity = TRUE, draws = 4000L, seed = 1L,
                  verbose = TRUE, return = TRUE, detail = FALSE) {
	x_name <- tryCatch(deparse1(substitute(x)), error = function(e) "x")
	y_name <- tryCatch(deparse1(substitute(y)), error = function(e) "y")
	clean <- function(s) {
		s <- sub("^.*[$]", "", s)
		if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", s) && nchar(s) <= 30L) s else NULL
	}
	x_name <- clean(x_name); y_name <- clean(y_name)
	if (is.null(x_name)) x_name <- "x"
	if (is.null(y_name)) y_name <- "y"
	call_expr <- match.call()
	if (isFALSE(return) && !("verbose" %in% names(call_expr))) verbose <- FALSE
	if (!is.numeric(prior_h1) || length(prior_h1) != 1L || !isTRUE(prior_h1 > 0 && prior_h1 < 1)) {
		.gbay_err("prior_h1 must be one number strictly between 0 and 1.", "prior_h1 doit \u00eatre un seul nombre strictement entre 0 et 1.")
	}
	if (!is.numeric(draws) || length(draws) != 1L || !isTRUE(draws >= 100 && draws <= 1e6)) {
		.gbay_err("draws must be one number between 100 and 1000000.", "draws doit \u00eatre un seul nombre entre 100 et 1000000.")
	}
	if (!is.numeric(seed) || length(seed) != 1L || !isTRUE(abs(seed) < .Machine$integer.max - 10)) {
		.gbay_err("seed must be one integer.", "seed doit \u00eatre un seul nombre entier.")
	}
	notes <- list(en = character(0), fr = character(0))
	story <- .bay_story(verbose)
	on.exit(.bay_story_flush(story), add = TRUE)
	note <- function(en, fr, tell = TRUE, kind = "fact") {
		notes$en <<- c(notes$en, en)
		notes$fr <<- c(notes$fr, fr)
		if (tell) .bay_tell_note(story, en, fr, kind)
	}
	inp <- .gbay_input(x, y, p, data, rescale.p, x_name, y_name, parent.frame())
	.gbay_run(inp, prior_h1, sensitivity, draws, seed, detail, verbose, return, story, note,
		function() notes, call_expr, !is.null(correct))
}

#' @export
print.gbay <- function(x, detail = isTRUE(x$detail), story = FALSE, ...) {
	told <- isTRUE(story) && length(x$story$en) > 0L
	if (told) .bay_story_replay(x$story)
	.gbay_print_bilan(x, detail = detail, narrated = told)
	invisible(x)
}
