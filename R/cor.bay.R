#' Bayesian counterpart of cor.test: probability of a link between two quantitative variables
#'
#' @description
#' `cor.bay()` takes the same inputs as [stats::cor.test()] and answers the opposite
#' question: not "how surprising is this correlation if there is no link?", but "how
#' probable is each hypothesis given these data?": a positive link, a negative link, no
#' link. Its main result is a bilan: the hypotheses ranked from the most to the least
#' probable, with their posterior probability, a plain qualifier and the Bayes factor, then
#' the strength of the link with a credible interval, and a sentence ready for a report.
#' A link observed without an experiment does not show a cause.
#'
#' A probability of a hypothesis depends on priors. The function sets them itself, in
#' advance, never on the data; the bilan states them and shows how far the answer moves when
#' they change. If they can change the leading hypothesis, it says the data do not settle it.
#'
#' @section Methods:
#' The DECISIONS file of the development folder (dev/corbay) gives the reasons and the
#' measures of each choice (C1 to C9).
#' - **Pearson** ([BayesFactor::correlationBF()]): bivariate normal model, Jeffreys's test
#'   with a stretched beta prior on rho (Ly, Verhagen and Wagenmakers 2016), width
#'   `kappa = 1/3` ("medium" of BayesFactor).
#' - **Kendall** (van Doorn, Ly, Marsman and Wagenmakers 2018): the standardised tau is
#'   normal around 3/2 tau sqrt(N); the same prior, carried over to tau. Based on ranks: a
#'   few far values cannot carry the answer, and any monotone link counts, straight or curved.
#' - **Spearman**, on request only: Pearson's Bayesian test on the ranks (the approach of the
#'   correlation package, Makowski et al. 2020). No established package gives the exact
#'   Bayesian Spearman of van Doorn et al. (2020).
#'
#' With `method = "auto"` (default), Kendall's tau is taken when a variable has a value
#' beyond 3 interquartile ranges, when a pair has a Cook's distance above 1 (from 20 pairs), or
#' when a variable has 7 distinct values or fewer; Pearson otherwise. No p-value chooses the
#' method. The other method is computed too, as a cross-check: if both answer firmly and
#' differently, the bilan does not decide.
#'
#' @param x Numeric vector, a column name or index of `data`, or a formula `~ x + y`.
#' @param y Numeric vector of the same length, or a column name or index of `data`.
#' @param data Data.frame, for the formula or the column names.
#' @param method `"auto"` (default), `"pearson"`, `"kendall"` or `"spearman"`.
#' @param alternative `"two.sided"` (default): three hypotheses, positive link, negative
#'   link, no link. `"greater"`: a positive link against no link; `"less"`: a negative link
#'   against no link. As in [stats::cor.test()], to be chosen before seeing the data.
#' @param verbose If `TRUE` (default), tells the analysis step by step as [m.bay()] does:
#'   plan, checks, method, starting point, solidity, strength of the link, bilan; "==>" for a
#'   finding, "-->" for a decision, "Note:" for anything else; every line within 80 columns
#'   (or `getOption("width")` if smaller). `print(result, story = TRUE)` replays the steps;
#'   `print(result)` gives the bilan alone; `print(result, detail = TRUE)` adds the method.
#' @param conf.level Level of the credible interval (highest density, [bayestestR::hdi()]).
#' @param prior_h1 Prior probability of a link. Default 1/2.
#' @param kappa Width of the stretched beta prior on rho. Default 1/3; 1 is uniform (JASP).
#' @param sensitivity If `TRUE` (default), the answer is recomputed with `kappa` halved and
#'   doubled and with starting chances from 1 against 2 to 2 against 1.
#' @param draws Posterior draws for the credible interval. Default 10000.
#' @param seed Seed of this call. The global random state is left untouched.
#' @param detail If `TRUE`, the bilan also gives the method in detail.
#' @param ... Other arguments of [stats::cor.test()] (`exact`, `continuity`): no role here,
#'   ignored with a message.
#'
#' @return A list of class `corbay`:
#'   \describe{
#'     \item{hypotheses}{Ranked hypotheses: `key` (`pos`, `neg`, `none`), `label_en`,
#'       `label_fr`, `posterior`, `prior`, `log_bf` (against no link), `null`.}
#'     \item{prob_link}{Posterior probability of a link (of the sign tested). Not a p-value.}
#'     \item{bf}{`log_bf10`, `log_bf_plus`, `log_bf_minus` (log Bayes factors against no link)
#'       and `error_pct`.}
#'     \item{estimate}{Under the hypothesis of a link, with Jeffreys's prior for estimation
#'       (stretched beta 1/2, recommended by van Doorn et al. 2018), which hardly pulls the
#'       estimate towards 0: `param` (`rho` or `tau`), `median`, `lo`, `hi`, `ci`.}
#'     \item{estimate_test}{The same under the prior of the test (`kappa`), shown with
#'       `detail = TRUE`; it can miss a strong link (dev/corbay/DECISIONS.md, C8).}
#'     \item{observed}{`r`, `tau`, `rho_s`: Pearson, Kendall and Spearman coefficients.}
#'     \item{method, auto, why}{The method, whether the rule chose it, and why.}
#'     \item{other}{The cross-check by the other method (`method = "auto"` only).}
#'     \item{contradiction}{`TRUE` if both methods answer firmly and differently.}
#'     \item{sensitivity}{Probability of the leading hypothesis over the priors, and verdict.}
#'     \item{checks}{Far values, influential pairs, ties.}
#'     \item{n, n_removed, alternative, priors, label, notes, call}{Description.}
#'     \item{story}{The steps told, in English (`en`) and French (`fr`).}
#'   }
#'
#' @references
#' Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor
#' hypothesis tests: explanation, extension, and application in psychology. Journal of
#' Mathematical Psychology, 72, 19-32. \doi{10.1016/j.jmp.2015.06.004}
#'
#' van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Bayesian inference for
#' Kendall's rank correlation coefficient. The American Statistician, 72(4), 303-308.
#' \doi{10.1080/00031305.2016.1264998}
#'
#' Makowski, D., Ben-Shachar, M. S., Patil, I., & Luedecke, D. (2020). Methods and algorithms
#' for correlation analysis in R. Journal of Open Source Software, 5(51), 2306.
#' \doi{10.21105/joss.02306}
#'
#' @seealso [stats::cor.test()], [m.bay()], [BayesFactor::correlationBF()]
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(30)
#' y <- 0.5 * x + rnorm(30)
#' r <- cor.bay(x, y)
#' r$prob_link
#'
#' # Formula, as cor.test; one direction tested.
#' d <- data.frame(size = x, mass = y)
#' cor.bay(~ size + mass, data = d, alternative = "greater", verbose = FALSE)
#'
#' # A far value: the rule takes Kendall's tau.
#' y[30] <- 12
#' cor.bay(x, y, verbose = FALSE, sensitivity = FALSE)$method
#'
#' @importFrom BayesFactor correlationBF extractBF
#' @importFrom bayestestR hdi
#' @export
cor.bay <- function(x, y = NULL, data = NULL, method = c("auto", "pearson", "kendall", "spearman"),
                    alternative = c("two.sided", "greater", "less"), verbose = TRUE, conf.level = 0.95,
                    prior_h1 = 1/2, kappa = 1/3, sensitivity = TRUE, draws = 10000L, seed = 1L,
                    detail = FALSE, ...) {
	method <- match.arg(method)
	alternative <- match.arg(alternative)
	call_expr <- match.call()
	if (!is.numeric(prior_h1) || length(prior_h1) != 1L || !isTRUE(prior_h1 > 0 && prior_h1 < 1)) {
		.exit("prior_h1 must be one number strictly between 0 and 1.",
			"prior_h1 doit \u00eatre un seul nombre strictement entre 0 et 1.")
	}
	if (!is.numeric(kappa) || length(kappa) != 1L || !isTRUE(kappa >= 0.01 && kappa <= 10)) {
		.exit("kappa must be one number between 0.01 and 10.", "kappa doit \u00eatre un seul nombre entre 0,01 et 10.")
	}
	if (!is.numeric(conf.level) || length(conf.level) != 1L || !isTRUE(conf.level > 0 && conf.level < 1)) {
		.exit("conf.level must be one number strictly between 0 and 1.",
			"conf.level doit \u00eatre un seul nombre strictement entre 0 et 1.")
	}
	if (!is.numeric(draws) || length(draws) != 1L || !isTRUE(draws >= 500 && draws <= 1e6)) {
		.exit("draws must be one number between 500 and 1e6.", "draws doit \u00eatre un seul nombre entre 500 et 1e6.")
	}
	if (!is.numeric(seed) || length(seed) != 1L || !isTRUE(abs(seed) < .Machine$integer.max - 10)) {
		.exit("seed must be one integer.", "seed doit \u00eatre un seul nombre entier.")
	}
	notes <- list(en = character(0), fr = character(0))
	story <- .bay_story(verbose)
	on.exit(.bay_story_flush(story), add = TRUE)
	note <- function(en, fr) {
		notes$en <<- c(notes$en, en)
		notes$fr <<- c(notes$fr, fr)
		.bay_tell_note(story, en, fr)
	}
	dots <- list(...)
	if (length(dots)) {
		nm <- paste(if (is.null(names(dots))) "..." else names(dots), collapse = ", ")
		note(paste0("Argument(s) without effect here, ignored: ", nm, "."), paste0("Argument(s) sans effet ici, ignor\u00e9(s) : ", nm, "."))
	}

	inp <- .corbay_input(x, y, data, substitute(x), substitute(y), parent.frame())
	xv <- inp$x; yv <- inp$y; xn <- inp$xn; yn <- inp$yn
	ok <- is.finite(xv) & is.finite(yv)
	n_removed <- sum(!ok)
	xv <- xv[ok]; yv <- yv[ok]
	n <- length(xv)
	.corbay_story_plan(story, xn, yn, n, alternative)
	if (n < 3L) {
		.exit(paste0("Only ", n, " complete pair(s): at least 3 are needed."),
			paste0("Seulement ", n, " paire(s) compl\u00e8te(s) : il en faut au moins 3."))
	}
	if (length(unique(xv)) < 2L || length(unique(yv)) < 2L) {
		cst <- if (length(unique(xv)) < 2L) xn else yn
		.exit(paste0(cst, " is constant: no link can be measured."),
			paste0(cst, " est constant : aucun lien ne peut \u00eatre mesur\u00e9."))
	}
	ck <- .corbay_checks(xv, yv, seed)
	.corbay_story_checks(story, ck, xn, yn, n_removed)
	ch <- .corbay_choose(ck)
	auto <- identical(method, "auto")
	meth <- if (auto) ch$method else method
	.corbay_story_method(story, meth, auto, ch$why, xn, yn)
	.corbay_story_prior(story, kappa, prior_h1, alternative, meth)

	lbf <- .corbay_lbf(meth, xv, yv, kappa)
	if (anyNA(lbf[c("plus", "minus")])) {
		.exit("The Bayes factor could not be computed on these data.",
			"Le facteur de Bayes n'a pas pu \u00eatre calcul\u00e9 sur ces donn\u00e9es.")
	}
	hyp <- .corbay_hyp(lbf, alternative, prior_h1, xn, yn)
	lead <- hyp$key[1L]
	lb10 <- switch(alternative, two.sided = lbf[["two"]], greater = lbf[["plus"]], less = lbf[["minus"]])
	other <- NULL
	contradiction <- FALSE
	if (auto) {
		om <- if (meth == "kendall") "pearson" else "kendall"
		ol <- .corbay_lbf(om, xv, yv, kappa)
		if (!anyNA(ol[c("plus", "minus")])) {
			oh <- .corbay_hyp(ol, alternative, prior_h1, xn, yn)
			other <- list(method = om, hypotheses = oh, log_bf10 = switch(alternative, two.sided = ol[["two"]], greater = ol[["plus"]], less = ol[["minus"]]))
			# Kendall chosen because of pairs apart: Pearson is the method those pairs can carry,
			# its disagreement is the reason of the choice, not a contradiction (DECISIONS, C9).
			explained <- meth == "kendall" && any(c("far", "influence") %in% ch$why)
			other$explained <- explained
			contradiction <- !explained && !identical(oh$key[1L], lead) && hyp$posterior[1L] >= 0.75 && oh$posterior[1L] >= 0.75
		}
	}
	sens <- if (isTRUE(sensitivity)) .corbay_sensitivity(meth, xv, yv, kappa, alternative, prior_h1, lead, xn, yn, lbf) else NULL
	# Estimation under Jeffreys's estimation prior, stretched beta 1/2 (DECISIONS, C8).
	est <- .corbay_estimate(meth, xv, yv, .corbay_kappa_est(), alternative, as.integer(draws), seed, conf.level)
	# The same estimate under the prior of the test, kept and shown with detail = TRUE.
	est_test <- .corbay_estimate(meth, xv, yv, kappa, alternative, as.integer(draws), seed, conf.level)
	observed <- c(r = stats::cor(xv, yv), tau = suppressWarnings(stats::cor(xv, yv, method = "kendall")),
		rho_s = suppressWarnings(stats::cor(xv, yv, method = "spearman")))
	if (alternative != "two.sided" && isTRUE(sign(observed[["r"]]) == (if (alternative == "greater") -1 else 1))) {
		note(paste0("The link observed goes the other way (", if (alternative == "greater") "negative" else "positive", "); it lies outside the hypotheses set by alternative = \"", alternative, "\"."),
			paste0("Le lien observ\u00e9 va dans l'autre sens (", if (alternative == "greater") "n\u00e9gatif" else "positif", ") ; il sort des hypoth\u00e8ses pos\u00e9es par alternative = \u00ab ", alternative, " \u00bb."))
	}
	out <- list(hypotheses = hyp, prob_link = sum(hyp$posterior[!hyp$null]),
		bf = list(log_bf10 = unname(lb10), log_bf_plus = unname(lbf[["plus"]]), log_bf_minus = unname(lbf[["minus"]]), error_pct = unname(lbf[["error_pct"]])),
		estimate = est, estimate_test = est_test, observed = observed, method = meth, auto = auto, why = ch$why, other = other,
		contradiction = contradiction, sensitivity = sens, checks = ck, n = n, n_removed = n_removed,
		alternative = alternative, priors = list(prior_h1 = prior_h1, kappa = kappa),
		label = paste0(xn, ", ", yn), notes = notes, call = call_expr, detail = isTRUE(detail))
	.corbay_story_solidity(story, out)
	.corbay_story_estimate(story, out)
	out$story <- .bay_story_lines(story)
	class(out) <- "corbay"
	if (isTRUE(verbose)) {
		.corbay_print_bilan(out, detail = detail, narrated = TRUE)
		return(invisible(out))
	}
	out
}

#' @export
print.corbay <- function(x, detail = isTRUE(x$detail), story = FALSE, ...) {
	told <- isTRUE(story) && length(x$story$en) > 0L
	if (told) .bay_story_replay(x$story)
	.corbay_print_bilan(x, detail = detail, narrated = told)
	invisible(x)
}

# Inputs as cor.test: x and y vectors; a formula ~ x + y (or y ~ x) with data; column names
# or indices of data. Returns the two numeric vectors and their names.
.corbay_input <- function(x, y, data, sx, sy, env) {
	nm_of <- function(s, def) {
		d <- tryCatch(deparse1(s), error = function(e) def)
		# A plain name (or data$name) only; an expression such as c(1, 2, 3) gives "x" or "y".
		if (length(d) != 1L || nchar(d) > 40L || !grepl("^[A-Za-z.][A-Za-z0-9._]*([$][A-Za-z.][A-Za-z0-9._]*)?$", d)) def else sub("^.*[$]", "", d)
	}
	col <- function(v, what) {
		if (is.null(data)) return(NULL)
		if (is.character(v) && length(v) == 1L) {
			if (!v %in% names(data)) .exit(paste0("Column ", v, " not found in data."), paste0("Colonne ", v, " absente de data."))
			return(list(v = data[[v]], n = v))
		}
		if (is.numeric(v) && length(v) == 1L && nrow(data) != 1L) {
			if (v < 1 || v > ncol(data)) .exit(paste0("Column ", v, " not found in data."), paste0("Colonne ", v, " absente de data."))
			return(list(v = data[[v]], n = names(data)[v]))
		}
		NULL
	}
	if (inherits(x, "formula")) {
		vars <- all.vars(x)
		if (length(vars) != 2L) {
			.exit("The formula must name two variables: ~ x + y.", "La formule doit nommer deux variables : ~ x + y.")
		}
		rhs <- all.vars(x[[length(x)]])
		vars <- if (length(x) == 3L) c(rhs, all.vars(x[[2L]])) else rhs
		mf <- tryCatch(stats::model.frame(stats::as.formula(paste("~", paste(vars, collapse = " + "))), data = if (is.null(data)) env else data,
			na.action = stats::na.pass), error = function(e) .exit(paste0("Formula not understood: ", conditionMessage(e), "."),
			paste0("Formule non comprise : ", conditionMessage(e), ".")))
		xv <- mf[[1L]]; yv <- mf[[2L]]; xn <- vars[1L]; yn <- vars[2L]
	} else {
		cx <- col(x, "x"); cy <- col(y, "y")
		xv <- if (is.null(cx)) x else cx$v
		yv <- if (is.null(cy)) y else cy$v
		xn <- if (is.null(cx)) nm_of(sx, "x") else cx$n
		yn <- if (is.null(cy)) nm_of(sy, "y") else cy$n
		if (is.null(yv)) .exit("y is missing: give two variables, or a formula ~ x + y.", "y manque : donner deux variables, ou une formule ~ x + y.")
	}
	for (v in list(list(xv, xn), list(yv, yn))) {
		if (!is.numeric(v[[1L]]) || is.logical(v[[1L]])) {
			.exit(paste0(v[[2L]], " is not numeric: cor.bay() links two quantitative variables (categories: m.bay())."),
				paste0(v[[2L]], " n'est pas num\u00e9rique : cor.bay() relie deux variables quantitatives (cat\u00e9gories : m.bay())."))
		}
	}
	if (length(xv) != length(yv)) .exit("x and y must have the same length.", "x et y doivent avoir la m\u00eame longueur.")
	if (identical(xn, yn)) { xn <- paste0(xn, "1"); yn <- paste0(yn, "2") }
	list(x = as.numeric(xv), y = as.numeric(yv), xn = xn, yn = yn)
}
