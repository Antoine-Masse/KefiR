#' Bayesian counterpart of m.test: probability of each hypothesis given the data
#'
#' @description
#' `m.bay()` takes the same inputs as [m.test()] and answers the opposite question:
#' not "how surprising are these data if there is no effect?", but "how probable is
#' each hypothesis given these data?". Its main result is a bilan: the hypotheses,
#' ranked from the most to the least probable, with their posterior probability and
#' a plain qualifier, then the answer for each term ("does the yield differ according
#' to the temperature? yes with this probability, no with that one"). A difference between
#' groups observed without randomisation does not show a cause.
#'
#' A probability of a hypothesis depends on priors. The bilan always states them and
#' always shows how far the answer moves when they change. If the priors can change
#' the leading hypothesis, it says that the data do not settle the question.
#'
#' `m.bay()` does not choose a method by preliminary tests as [m.test()] does. It
#' compares the plausible shapes of the noise (normal or heavy tails, equal or
#' unequal variances) and averages its answers over their posterior probabilities.
#'
#' @section Plans covered:
#' Available now: one factor, 2 groups or k groups, with the four shapes of noise
#' compared; several crossed factors (for example `y ~ A * B`), with the normal law and a
#' common variance only. With several factors, the marginal likelihoods come from the
#' BayesFactor package (Rouder et al. 2012) and each term is answered as "does y differ
#' according to A, alone or in interaction?". An interaction enters only with its factors,
#' and the prior on the models gives one proportion ~ Beta(1,1) per order of interaction
#' (Taylor-Rodriguez et al. 2016), which corrects for the number of terms; `prior_h1` and
#' `laws` have no role there. Two costs of this prior (DECISIONS, D2): when one interaction
#' is real, the other interactions of the same order are more often called "yes"; and an
#' interaction without its lower-order terms (A:B:C alone) is outside the hereditary models
#' and is poorly seen.
#'
#' Paired data and repeated measures (DECISIONS, D15), written as in [m.test()]
#' (`y ~ time | id`, `Error(id/time)`, `(1 | id)`, `id =`, `paired = TRUE`): the subject is a
#' random effect present in every model, "no difference" included. Two levels: Bayesian paired
#' t test on the differences, [BayesFactor::ttestBF()] with `paired = TRUE` (Rouder et al. 2009),
#' and the size of the difference as dz with its credible interval. More levels, or factors
#' between subjects too: [BayesFactor::lmBF()] with `whichRandom` (Rouder et al. 2012), each
#' model against the subjects alone; one factor within subject gets the post-hoc on
#' groupings. Normal law only. Several measures of one subject in one condition (replicates,
#' DECISIONS, D21): the subject and the subject crossed with the factors within subject are
#' random, so that the replicates do not count as more subjects; random slopes are refused.
#'
#' Analysis of covariance (DECISIONS, D16), `y ~ g + x` or `y ~ g * x`: covariates centred on
#' their mean, present in every model, each answer read against the covariates alone
#' ([BayesFactor::lmBF()], slope prior `rscaleCont` = sqrt(2)/4, varied with `r`). Groups that
#' hardly overlap on a covariate cap the support at "weak". Parallel slopes: adjusted gaps
#' and groupings as with one factor. Slopes that differ: gaps at the quartiles of the
#' covariate and the zone where the groups differ (Bayesian Johnson-Neyman). A covariate
#' alone, and covariates with subjects, are refused. The other plans of [m.test()] (mixed,
#' MANOVA, discrete data) stop with an explicit message until they are validated.
#'
#' @param x Response: numeric vector, column name or index if `data` is given, or a
#'   formula.
#' @param g Grouping: factor or vector, column name or index if `data` is given, or
#'   a data.frame of factors.
#' @param data Data.frame.
#' @param formula Formula, for example `y ~ F`. Same forms as [m.test()].
#' @param paired,id,wt,within,between Structure of the plan, as in [m.test()].
#' @param alpha,boot,boot_type,iter,conf,silent Accepted for compatibility with
#'   [m.test()]; they have no role here and are ignored, with a message if given.
#' @param control Reference level. Not used yet.
#' @param verbose If `TRUE`, tells the analysis step by step as [m.test()] does, in the same
#'   seven steps for every plan (plan, checks, shape of the noise, starting point, solidity,
#'   pairs, bilan): the result of each step in its title, the function and package in
#'   brackets, "==>" for a finding, "-->" for a decision, "Note:" for anything else; every
#'   line within 80 columns (or `getOption("width")` if smaller). If `FALSE`, nothing.
#'   `print(result, story = TRUE)` replays the steps; `print(result)` gives the bilan alone;
#'   `print(result, detail = TRUE)` adds the method in detail (priors, references).
#' @param plot Not used yet.
#' @param return If `TRUE`, the full object (invisibly when the bilan was printed);
#'   if `FALSE`, one number: the posterior probability that at least one effect
#'   exists. This is not a p-value.
#' @param maxcat Maximum number of levels of a factor.
#' @param code,debug Accepted for compatibility with [m.test()]. Not used yet.
#' @param prior_h1 Prior probability of an effect. Default 1/2.
#' @param rscale Prior scale r of a fixed effect (0.5 is the "medium" of the
#'   BayesFactor package; on a difference between 2 groups it gives a Cauchy of
#'   scale 0.707).
#' @param u Parameter of the Dirichlet prior on the proportions of precision between
#'   groups, in the laws with unequal variances. Default 1.5.
#' @param laws Sampling laws compared, among `"normal_eq"`, `"normal_het"`, `"t_eq"`
#'   and `"t_het"`.
#' @param sensitivity If `TRUE` (default), computes the sensitivity to the priors.
#'   The bilan says so when it is switched off.
#' @param draws Importance draws per model. Default 4000. With several factors, BayesFactor
#'   uses at least 10000 iterations; its noise hardly shrinks beyond, so a larger `draws`
#'   brings little there.
#' @param seed Seed of this call. The global random state is left untouched.
#' @param posthoc If `TRUE` (default), one factor: which groups differ and by how much. Every
#'   grouping of the groups into blocks of equal means is compared through [BayesFactor::lmBF()]
#'   (all of them up to 8 groups; beyond, explored by the Gibbs sampler of Gopalan and Berry
#'   1998), the probability that two groups differ is summed over the groupings that
#'   separate them, and each pair gets its difference and standardised difference with a 95 %
#'   credible interval ([BayesFactor::posterior()], [bayestestR::hdi()]). Normal law with a
#'   common variance; the story says so when the noise looks otherwise.
#' @param response `"auto"` (default): the function chooses. A factor, text, TRUE/FALSE or
#'   two values is categorical; with one factor, whole numbers from 0 with a group mean under
#'   10 and zeros (or more than 7 values) are counts; anything else is continuous.
#'   `"counts"` forces the count route: Poisson and negative binomial compared
#'   (overdispersion), Bayes factors by the unit-information approximation (BIC), rate
#'   ratios. A bounded rating scale (maximum at 4, 5, 6, 7, 9, 10, 20 or 100, at most 10
#'   values) takes the ordinal route: proportional odds, [MASS::polr()], cumulative odds
#'   ratios; `"ordinal"` forces it. `"continuous"` forces the continuous route (DECISIONS, D17).
#' @param detail If `TRUE`, the bilan also shows the shapes of the noise, the priors that
#'   are varied and how to read the numbers. Also `print(result, detail = TRUE)`.
#'
#' @return With `return = TRUE`, a list of class `mbay`:
#'   \describe{
#'     \item{hypotheses}{Ranked hypotheses: `label_en`, `label_fr`, `posterior`,
#'       `prior`.}
#'     \item{effects}{One row per term: `term`, `bf10`, `log_bf`, `p_h1`, `p_h0`,
#'       `prior_h1`, `rscale`, `error_pct`, `analogue`. No column is a p-value.}
#'     \item{laws}{Posterior probability of each sampling law, with the log Bayes
#'       factor within that law.}
#'     \item{models}{The eight models with their marginal likelihoods.}
#'     \item{sensitivity}{Grid over r and odds, range over the priors, verdict.}
#'     \item{priors}{The priors used.}
#'     \item{descriptive_stats}{n, mean and SD per group.}
#'     \item{prob_effect}{Posterior probability that at least one effect exists.}
#'     \item{method, plan, label, n, notes, call}{Description of the analysis.}
#'     \item{posthoc}{One factor: `pairs` (probability that two groups differ, its range over the
#'       prior on groupings, differences and standardised differences with 95 % credible
#'       intervals) and `partitions` (groupings ranked by posterior probability).}
#'     \item{story}{The steps told with `verbose = TRUE`, in English (`en`) and French
#'       (`fr`), kept even when nothing was printed.}
#'   }
#'
#' @seealso [m.test()]
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(y = c(rnorm(15, 10, 2), rnorm(15, 12, 2), rnorm(15, 10.5, 2)),
#'                 variety = rep(c("a", "b", "c"), each = 15))
#' r <- m.bay(y ~ variety, data = d, draws = 1000, sensitivity = FALSE)
#' r$effects[, c("term", "p_h1", "p_h0", "log_bf")]
#'
#' @importFrom BayesFactor lmBF extractBF posterior ttestBF
#' @importFrom bayestestR hdi
#' @export
m.bay <- function(x = NULL, g = NULL, data = NULL, formula = NULL,
                  paired = FALSE, id = NULL, wt = NULL, within = NULL, between = NULL,
                  alpha = NULL, control = NULL, verbose = TRUE, plot = FALSE,
                  return = TRUE, boot = NULL, boot_type = NULL, iter = NULL, conf = NULL,
                  maxcat = 50, silent = NULL, code = FALSE, debug = FALSE,
                  prior_h1 = 1/2, rscale = 1/2, u = 1.5,
                  laws = c("normal_eq", "normal_het", "t_eq", "t_het"),
                  sensitivity = TRUE, draws = 4000L, seed = 1L, detail = FALSE, posthoc = TRUE,
                  response = c("auto", "counts", "ordinal", "continuous")) {
	response <- match.arg(response)
	# Name of the response in the x, g form (audit of de64470, point 1).
	x_name <- if (!missing(x) && !inherits(x, "formula") && !is.character(x)) tryCatch(deparse1(substitute(x)), error = function(e) NULL) else NULL
	call_expr <- match.call()
	if (isFALSE(return) && !("verbose" %in% names(call_expr))) verbose <- FALSE

	# Checks on the priors. Nothing is ever tuned on the data (DECISIONS, D6).
	if (!is.numeric(prior_h1) || length(prior_h1) != 1L || !isTRUE(prior_h1 > 0 && prior_h1 < 1)) {
		.exit("prior_h1 must be one number strictly between 0 and 1.",
			"prior_h1 doit \u00eatre un seul nombre strictement entre 0 et 1.")
	}
	rscale <- .bay_check_rscale(rscale)
	if (!is.numeric(u) || length(u) != 1L || !isTRUE(is.finite(u) && u > 0)) {
		.exit("u must be one finite, strictly positive number.",
			"u doit \u00eatre un seul nombre fini et strictement positif.")
	}
	if (!isTRUE(rscale >= 1e-3 && rscale <= 1e3) || !isTRUE(u >= 1e-2 && u <= 1e4)) {
		.exit("rscale must lie between 0.001 and 1000, and u between 0.01 and 10000.",
			"rscale doit \u00eatre entre 0,001 et 1000, et u entre 0,01 et 10000.")
	}
	if (!is.numeric(draws) || length(draws) != 1L || !isTRUE(draws >= 100 && draws <= .bay_law_max_draws())) {
		.exit(paste0("draws must be one number between 100 and ", .bay_law_max_draws(), "."),
			paste0("draws doit \u00eatre un seul nombre entre 100 et ", .bay_law_max_draws(), "."))
	}
	if (!is.numeric(seed) || length(seed) != 1L || !isTRUE(abs(seed) < .Machine$integer.max - 10)) {
		.exit("seed must be one integer.", "seed doit \u00eatre un seul nombre entier.")
	}
	laws <- unique(laws)
	if (!length(laws) || !all(laws %in% .bay_law_names())) {
		.exit("laws must be among normal_eq, normal_het, t_eq, t_het.",
			"laws doit \u00eatre parmi normal_eq, normal_het, t_eq, t_het.")
	}

	ignored <- c(alpha = !is.null(alpha), boot = !is.null(boot), boot_type = !is.null(boot_type),
		iter = !is.null(iter), conf = !is.null(conf), silent = !is.null(silent),
		code = isTRUE(code), debug = isTRUE(debug), control = !is.null(control),
		plot = isTRUE(plot))
	notes <- list(en = character(0), fr = character(0))
	# The story of the call (verbose = TRUE): every note is also said at the step where it
	# arises, unless tell = FALSE (the step then says it in its own words).
	story <- .bay_story(verbose)
	# A refusal or an error still shows the step being told (its lines wait for the next step).
	on.exit(.bay_story_flush(story), add = TRUE)
	note <- function(en, fr, tell = TRUE, kind = "fact") {
		notes$en <<- c(notes$en, en)
		notes$fr <<- c(notes$fr, fr)
		if (tell) .bay_tell_note(story, en, fr, kind)
	}
	if (any(ignored)) {
		nm <- paste(names(ignored)[ignored], collapse = ", ")
		note(paste0("Argument(s) without effect here, ignored: ", nm, "."),
			paste0("Argument(s) sans effet ici, ignor\u00e9(s) : ", nm, "."))
	}

	# Whatever the input, no raw R error: what was not understood is said as such.
	input <- tryCatch(.bay_parse_input(x = x, g = g, data = data, formula = formula, paired = paired,
		id = id, wt = wt, within = within, between = between, maxcat = maxcat,
		env = parent.frame()), error = function(e) {
		msg <- conditionMessage(e)
		if (is.null(conditionCall(e))) stop(msg, call. = FALSE)
		.exit(paste0("Input not understood: ", msg, "."),
			paste0("Entr\u00e9es non comprises : ", msg, "."))
	})
	if (input$n_removed > 0L) {
		note(paste0(input$n_removed, " observation(s) with missing values removed."),
			paste0(input$n_removed, " observation(s) avec valeurs manquantes retir\u00e9e(s)."))
	}

	if (identical(input$response_name, "x") && length(x_name) == 1L && grepl("^[A-Za-z.][A-Za-z0-9._$]*$", x_name) && nchar(x_name) <= 30L) {
		nm <- sub("^.*[$]", "", x_name)
		input$response_name <- nm
		input$label <- sub("^x ~", paste(nm, "~"), input$label)
	}
	input$counts <- identical(response, "counts")
	# Counts chosen alone, before the plan is told (DECISIONS, D17, way 2).
	route <- if (identical(response, "auto") && identical(input$plan, "oneway") && !isTRUE(input$categorical))
		.bay_int_route(input$response, input$factors[[1L]]) else "none"
	input$counts_auto <- identical(route, "counts")
	if (input$counts_auto) input$counts <- TRUE
	input$ordinal <- identical(response, "ordinal") || identical(route, "ordinal")
	input$ordinal_auto <- identical(route, "ordinal")
	input$bounded <- identical(route, "bounded")
	# Several responses (DECISIONS, D22): their own route, which tells its own plan.
	if (identical(input$plan, "manova")) {
		return(.bay_m_manova(input, call_expr, function() notes, note, prior_h1, rscale, sensitivity, draws, seed, detail, verbose,
			return, story, posthoc))
	}
	# Step 1: the plan recognised, or refused and why.
	.bay_story_plan(story, input)

	# Plans not validated yet: an explicit refusal, never a wrong answer.
	if (identical(input$plan, "ancova") && !length(input$factors)) {
		.exit("A covariate alone (no factor) is a regression, not handled by m.bay() yet.",
			"Une covariable seule (sans facteur) est une r\u00e9gression, pas encore trait\u00e9e par m.bay().")
	}
	if (isTRUE(paired) && !identical(input$plan, "paired")) {
		.exit("paired = TRUE, but the measures do not form pairs (several measures per subject and level, or a factor between subjects); remove paired = TRUE to treat them as repeated measures.",
			"paired = TRUE, mais les mesures ne forment pas des paires (plusieurs mesures par sujet et niveau, ou un facteur entre sujets) ; retirer paired = TRUE pour les traiter en mesures r\u00e9p\u00e9t\u00e9es.")
	}
	if (identical(input$plan, "mixed") && length(input$slopes)) {
		.exit("Random slopes are not handled by m.bay() (they need Stan, not available here).",
			"Les pentes al\u00e9atoires ne sont pas trait\u00e9es par m.bay() (il faudrait Stan, indisponible ici).")
	}
	if (input$plan %in% c("paired", "repeated", "mixed") && length(input$covariates)) {
		.exit("Covariates with subjects (paired ANCOVA) are not handled by m.bay() yet.",
			"Des covariables avec des sujets (ANCOVA appari\u00e9e) ne sont pas encore trait\u00e9es par m.bay().")
	}
	if (!input$plan %in% c("oneway", "factorial", "paired", "repeated", "mixed", "ancova")) {
		plan_en <- c(factorial = "several factors", ancova = "covariates", paired = "paired data",
			repeated = "repeated measures", mixed = "mixed models", manova = "MANOVA")[input$plan]
		plan_fr <- c(factorial = "plusieurs facteurs", ancova = "les covariables", paired = "les donn\u00e9es appari\u00e9es",
			repeated = "les mesures r\u00e9p\u00e9t\u00e9es", mixed = "les mod\u00e8les mixtes", manova = "la MANOVA")[input$plan]
		.exit(paste0("m.bay() does not handle ", plan_en, " yet: this plan is not validated (use m.test())."),
			paste0("m.bay() ne traite pas encore ", plan_fr, " : ce plan n'est pas valid\u00e9 (utiliser m.test())."))
	}
	y <- input$response
	.bay_tell(story, "Checks on the data.", "Contr\u00f4les des donn\u00e9es.", "step")
	# Binary or categorical response (DECISIONS, D17, way 1): a contingency table with one factor.
	cat_resp <- isTRUE(input$categorical) || (is.numeric(y) && length(unique(y)) == 2L &&
		!isTRUE(diff(range(y)) <= 8 * .Machine$double.eps * max(abs(y))))
	if (cat_resp) {
		if (!identical(input$plan, "oneway")) {
			.exit("Binary or categorical response with several factors, covariates or subjects: not handled by m.bay() yet (one factor only).",
				"R\u00e9ponse binaire ou cat\u00e9gorielle avec plusieurs facteurs, des covariables ou des sujets : pas encore trait\u00e9e par m.bay() (un facteur seulement).")
		}
		yc <- if (is.factor(y)) y else factor(y)
		if (nlevels(yc) < 2L) .exit("The response has a single category.", "La r\u00e9ponse n'a qu'une seule modalit\u00e9.")
		return(.bay_m_contingency(input, input$factors[[1L]], yc, call_expr, function() notes, note, prior_h1, sensitivity,
			draws, seed, detail, verbose, return, story, posthoc, !("laws" %in% names(call_expr)), laws))
	}
	if (any(is.infinite(y))) {
		.exit("y contains infinite values.", "y contient des valeurs infinies.")
	}
	if (isTRUE(input$ordinal)) {
		if (!identical(input$plan, "oneway")) {
			.exit("Ordered response with several factors, covariates or subjects: not handled by m.bay() yet (one factor only).",
				"R\u00e9ponse ordonn\u00e9e avec plusieurs facteurs, des covariables ou des sujets : pas encore trait\u00e9e par m.bay() (un facteur seulement).")
		}
		return(.bay_m_ordinal(input, input$factors[[1L]], y, call_expr, function() notes, note, prior_h1, sensitivity,
			detail, verbose, return, story, posthoc, auto = isTRUE(input$ordinal_auto)))
	}
	if (isTRUE(input$bounded)) {
		vals <- sort(unique(y))
		.bay_tell(story, paste0("Bounded score (", length(vals), " values, from ", min(vals), " to ", max(vals), "): neither counts nor few ratings; kept continuous (response = \"ordinal\" for the ordinal route)."),
			paste0("Score born\u00e9 (", length(vals), " valeurs, de ", min(vals), " \u00e0 ", max(vals), ") : ni comptage, ni notation \u00e0 peu de valeurs ; gard\u00e9 continu (response = \u00ab ordinal \u00bb pour la voie ordinale)."), "decide")
	}
	if (isTRUE(input$counts)) {
		if (!identical(input$plan, "oneway")) {
			.exit("Counts with several factors, covariates or subjects: not handled by m.bay() yet (one factor only).",
				"Comptages avec plusieurs facteurs, des covariables ou des sujets : pas encore trait\u00e9s par m.bay() (un facteur seulement).")
		}
		if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
			.exit("response = \"counts\" needs non-negative whole numbers.", "response = \"counts\" demande des nombres entiers positifs ou nuls.")
		}
		if (all(y == y[1L])) .exit("All the counts are equal: the probabilities are undetermined.", "Tous les comptages sont \u00e9gaux : les probabilit\u00e9s sont ind\u00e9termin\u00e9es.")
		return(.bay_m_counts(input, input$factors[[1L]], round(y), call_expr, function() notes, note, prior_h1, sensitivity,
			detail, verbose, return, story, posthoc, auto = isTRUE(input$counts_auto)))
	}

	y_values <- length(unique(y))
	if (isTRUE(diff(range(y)) <= 8 * .Machine$double.eps * max(abs(y)))) {
		.exit("y is constant: the probabilities are undetermined.",
			"y est constant : les probabilit\u00e9s sont ind\u00e9termin\u00e9es.")
	}
	if (y_values == 2L) {
		.exit("Binary response: its own law is not validated yet, the normal and t laws do not apply.",
			"R\u00e9ponse binaire : sa loi propre n'est pas encore valid\u00e9e, les lois normale et t ne s'appliquent pas.")
	}
	if (isTRUE(input$discrete)) {
		note("The response looks discrete: the normal and heavy-tailed laws are only an approximation here.",
			"La r\u00e9ponse semble discr\u00e8te : les lois normale et \u00e0 queues lourdes ne sont ici qu'une approximation.")
	}

	if (input$plan %in% c("factorial", "paired", "repeated", "mixed", "ancova")) {
		return(.bay_m_factorial(input, call_expr, function() notes, note, prior_h1, rscale, laws, sensitivity,
			draws, seed, detail, verbose, return, !("prior_h1" %in% names(call_expr)), !("laws" %in% names(call_expr)),
			story, posthoc))
	}

	term <- names(input$factors)[1L]
	shown <- if (identical(term, ".g")) "g" else term
	gr <- input$factors[[1L]]
	st <- .bay_laws_stats(y, gr)
	chk <- .bay_story_checks(story, y, gr)
	if (min(st$n) < 3L) {
		small <- st$levels[st$n < 3L]
		note(paste0("Very small group(s) (n < 3), kept: ", paste(small, collapse = ", "),
				". Their variance is poorly known."),
			paste0("Tr\u00e8s petit(s) groupe(s) (n < 3), conserv\u00e9(s) : ", paste(small, collapse = ", "),
				". Leur variance est mal connue."))
	}

	laws_kept <- .bay_laws_allowed(st, laws, u)
	if (length(laws_kept) == 0L) {
		.exit("Repeated values: none of the requested laws applies. Allow normal_eq in laws.",
			"Valeurs r\u00e9p\u00e9t\u00e9es : aucune des lois demand\u00e9es ne s'applique. Autoriser normal_eq dans laws.")
	}
	lab <- .bay_law_labels()
	if (length(laws_kept) < length(laws)) {
		left <- setdiff(laws, laws_kept)
		note(paste0("Repeated values: for these data the integral of some laws diverges, they were left out: ",
				paste(lab$en[left], collapse = "; "), ". A law for discrete data is not validated yet."),
			paste0("Valeurs r\u00e9p\u00e9t\u00e9es : pour ces donn\u00e9es l'int\u00e9grale de certaines lois diverge, elles ont \u00e9t\u00e9 \u00e9cart\u00e9es : ",
				paste(lab$fr[left], collapse = " ; "), ". Une loi pour donn\u00e9es discr\u00e8tes n'est pas encore valid\u00e9e."), tell = FALSE)
		nl <- length(left)
		.bay_tell(story, paste0("Change of strategy: ", if (nl == 1L) "one law is" else paste(nl, "laws are"), " left out, without a finite value with a group without spread: ",
				paste(lab$en[left], collapse = "; "), "."),
			paste0("Changement de strat\u00e9gie : ", if (nl == 1L) "une loi \u00e9cart\u00e9e" else paste(nl, "lois \u00e9cart\u00e9es"), ", sans valeur finie avec un groupe sans dispersion : ",
				paste(lab$fr[left], collapse = " ; "), "."), "decide")
	} else if (chk$extremes[["mild"]] > 0L) {
		.bay_tell(story, "No value is removed: the heavy-tailed shapes are there for such values.",
			"Aucune valeur retir\u00e9e : les formes \u00e0 queues lourdes sont l\u00e0 pour ce genre de valeurs.", "decide")
	} else {
		.bay_tell(story, "Nothing to remove.", "Rien \u00e0 retirer.", "decide")
	}
	laws <- laws_kept
	draws <- as.integer(draws)
	seed <- as.integer(seed)

	# Step 3: the shapes of the noise, compared and not chosen.
	.bay_tell(story, "Shapes of the noise", "Formes du bruit", "step", fun = "m.bay() {KefiR}")
	.bay_tell(story, "the shape of the noise is how the values scatter in their group; heavy tails mean far values more frequent than under the normal law.",
		"la forme du bruit est la fa\u00e7on dont les valeurs s'\u00e9cartent dans leur groupe ; les queues lourdes, des valeurs \u00e9loign\u00e9es plus fr\u00e9quentes que sous la loi normale.", "note")
	.bay_tell(story, "Shape of the noise: how the values scatter around the centre of their group (normal or heavy tails, equal variances or not). Computed by importance sampling, eight models.",
		"Forme du bruit : la fa\u00e7on dont les valeurs s'\u00e9cartent du centre de leur groupe (normale ou \u00e0 queues lourdes, variances \u00e9gales ou non). Calcul par \u00e9chantillonnage d'importance, huit mod\u00e8les.", "detail")
	rs <- if (isTRUE(sensitivity)) c(rscale / 2, rscale, rscale * 2) else rscale
	tab <- .bay_laws_table(st, rs, u, draws, seed, laws = laws)
	ans <- .bay_laws_average(tab, rscale, prior_h1 = prior_h1)
	if (length(ans$laws_dropped)) {
		note(paste0("Law(s) left out because the computation is not reliable here, so the comparison is incomplete and the answer may change: ",
				paste(lab$en[ans$laws_dropped], collapse = "; "), ". Try a larger draws."),
			paste0("Loi(s) \u00e9cart\u00e9e(s) car le calcul n'est pas fiable ici : la comparaison est incompl\u00e8te et la r\u00e9ponse peut changer : ",
				paste(lab$fr[ans$laws_dropped], collapse = " ; "), ". Essayer un draws plus grand."), tell = FALSE)
		.bay_tell(story, paste0("Too few effective draws for: ", paste(lab$en[ans$laws_dropped], collapse = "; "), "."),
			paste0("Trop peu de tirages efficaces pour : ", paste(lab$fr[ans$laws_dropped], collapse = " ; "), "."))
		.bay_tell(story, "Change of strategy: these laws are left out; the answer may change (try a larger draws).",
			"Changement de strat\u00e9gie : ces lois sont \u00e9cart\u00e9es ; la r\u00e9ponse peut changer (essayer un draws plus grand).", "decide")
	}
	pl <- sort(ans$p_law, decreasing = TRUE)
	.bay_step_result(story, paste0(lab$en[names(pl)[1L]], " ahead (", .bay_pct(pl[1L], "."), ")"),
		paste0(lab$fr[names(pl)[1L]], " en t\u00eate (", .bay_pct(pl[1L], ","), ")"))
	.bay_tell(story, paste0("Probabilities: ", paste0(lab$en[names(pl)], " ", vapply(pl, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
		paste0("Probabilit\u00e9s : ", paste0(lab$fr[names(pl)], " ", vapply(pl, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))

	# What the comparison of the shapes finds about the variances, said as a finding (the
	# analogue of the check of variances in m.test()); the tails are said by the decision below.
	sds <- vapply(split(y, gr), function(v) if (length(v) > 1L) stats::sd(v) else NA_real_, numeric(1))
	if (any(names(ans$p_law) %in% c("normal_het", "t_het")) && any(names(ans$p_law) %in% c("normal_eq", "t_eq"))) {
		p_het <- sum(ans$p_law[intersect(names(ans$p_law), c("normal_het", "t_het"))])
		hi_sd <- which.max(sds); lo_sd <- which.min(sds)
		if (identical(hi_sd, lo_sd)) lo_sd <- setdiff(seq_along(sds), hi_sd)[1L]
		sd_txt <- function(dc) paste0(names(sds)[hi_sd], " ", gsub(".", dc, format(signif(sds[hi_sd], 3), trim = TRUE), fixed = TRUE),
			if (dc == ",") " contre " else " against ", names(sds)[lo_sd], " ", gsub(".", dc, format(signif(sds[lo_sd], 3), trim = TRUE), fixed = TRUE))
		alike <- isTRUE(sds[hi_sd] / sds[lo_sd] < 1.3)
		if (p_het >= 0.75 && alike) {
			.bay_tell(story, paste0("Unequal spreads probable (", .bay_pct(p_het, "."), "), standard deviations close (", sd_txt("."), ")."),
				paste0("Dispersions in\u00e9gales probables (", .bay_pct(p_het, ","), "), \u00e9carts-types proches (", sd_txt(","), ")."))
			.bay_tell(story, "The groups differ rather by their shape (asymmetry, far values) than by their spread.",
				"Les groupes diff\u00e8rent plut\u00f4t par leur forme (asym\u00e9trie, valeurs \u00e9loign\u00e9es) que par leur dispersion.", "note")
			.bay_tell(story, "No common spread is imposed.", "Aucune dispersion commune n'est impos\u00e9e.", "decide")
		} else if (p_het >= 0.75) {
			.bay_tell(story, paste0("Unequal variances probable (", .bay_pct(p_het, "."), "): ", names(sds)[hi_sd],
					" the most dispersed (standard deviation ", sd_txt("."), ")."),
				paste0("Variances in\u00e9gales probables (", .bay_pct(p_het, ","), ") : ", names(sds)[hi_sd],
					" le plus dispers\u00e9 (\u00e9cart-type ", sd_txt(","), ")."))
			.bay_tell(story, "No common variance is imposed.", "Aucune variance commune n'est impos\u00e9e.", "decide")
		} else if (p_het <= 0.25) {
			.bay_tell(story, paste0("Equal variances probable (", .bay_pct(1 - p_het, "."), ")."),
				paste0("Variances \u00e9gales probables (", .bay_pct(1 - p_het, ","), ")."))
		} else {
			lean_en <- if (p_het < 0.5) "lean towards equal variances" else "lean towards unequal variances"
			lean_fr <- if (p_het < 0.5) "penchent vers des variances \u00e9gales" else "penchent vers des variances in\u00e9gales"
			.bay_tell(story, paste0("Variances: the data ", lean_en, " without settling it (unequal at ", .bay_pct(p_het, "."), ")."),
				paste0("Variances : les donn\u00e9es ", lean_fr, " sans trancher (in\u00e9gales \u00e0 ", .bay_pct(p_het, ","), ")."))
		}
	}


	# Under heavy tails the groups are compared through their centres, which are not the
	# means of a skewed group: equal means, different shapes can read as a difference
	# (agent_43: one lognormal group, n = 400, P > 0.95 in 67 % of the data sets). When the
	# normal law and the answer fall on opposite sides, the bilan says so on its first
	# lines and caps the qualifier. Thresholds and false alarms: agent_46, VALIDATION.
	p_t <- sum(ans$p_law[intersect(names(ans$p_law), c("t_eq", "t_het"))])
	normal_kept <- intersect(names(ans$p_law), c("normal_eq", "normal_het"))
	contradiction <- FALSE
	p_normal <- NA_real_
	if (length(normal_kept)) {
		best_normal <- normal_kept[which.max(ans$p_law[normal_kept])]
		p_normal <- .bay_post_prob(log_bf = ans$log_bf_law[best_normal], prior_odds = prior_h1 / (1 - prior_h1))
		contradiction <- isTRUE((ans$p_h1 > 0.75 && p_normal < 0.5) || (ans$p_h1 < 0.25 && p_normal > 0.5))
	}
	heavy <- isTRUE(p_t >= .bay_heavy_threshold())
	if (heavy && !contradiction) {
		note("The noise looks heavy-tailed: the comparison bears on the centres of the groups, not exactly on their means.",
			"Le bruit semble \u00e0 queues lourdes : la comparaison porte sur le centre des groupes, pas exactement sur leur moyenne.", tell = FALSE)
	}
	if (heavy) {
		.bay_tell(story, paste0("Heavy tails prevail (", .bay_pct(p_t, "."), "): the answer bears on the centres of the groups, not their means."),
			paste0("Les queues lourdes l'emportent (", .bay_pct(p_t, ","), ") : la r\u00e9ponse porte sur les centres des groupes, pas sur leurs moyennes."), "decide")
	} else if (p_t > 0.5) {
		.bay_tell(story, paste0("Heavy tails the more probable (", .bay_pct(p_t, "."), "), below ", .bay_pct(.bay_heavy_threshold(), "."), ": the answer mixes means and centres."),
			paste0("Queues lourdes les plus probables (", .bay_pct(p_t, ","), "), sous ", .bay_pct(.bay_heavy_threshold(), ","), " : la r\u00e9ponse m\u00eale moyennes et centres."), "decide")
		.bay_tell(story, "The centre is the typical value of a group, close to its median.", "Le centre est la valeur typique d'un groupe, proche de sa m\u00e9diane.", "note")
	} else if (length(normal_kept)) {
		.bay_tell(story, paste0("The answer bears mainly on the means (normal shapes: ", .bay_pct(1 - p_t, "."), ")."),
			paste0("La r\u00e9ponse porte surtout sur les moyennes (formes normales : ", .bay_pct(1 - p_t, ","), ")."), "decide")
	}
	if (length(pl) > 1L && pl[1L] > 0.95) {
		.bay_tell(story, "One shape dominates, the answer is almost its own.", "Une forme domine, la r\u00e9ponse est presque la sienne.", "note")
	} else if (length(pl) > 1L) {
		.bay_tell(story, "No shape is chosen; each answer is weighted by its probability.",
			"Aucune forme n'est choisie ; chaque r\u00e9ponse est pond\u00e9r\u00e9e par sa probabilit\u00e9.", "note")
	}

	# Step 4: the priors, fixed before the computation.
	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	r_txt <- .bay_num_txt(rscale)
	user_prior <- "prior_h1" %in% names(call_expr)
	user_r <- "rscale" %in% names(call_expr)
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference", if (user_prior) " (given by the user)" else " (heads or tails)"),
		paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence", if (user_prior) " (donn\u00e9 par l'utilisateur)" else " (pile ou face)"))
	.bay_tell(story, paste0("Expected size r = ", r_txt$en, if (user_r) " (given by the user)" else "", "; each shape of the noise starts at ", .bay_pct(1 / length(laws), "."), "."),
		paste0("Taille attendue r = ", r_txt$fr, if (user_r) " (donn\u00e9e par l'utilisateur)" else "", " ; chaque forme du bruit part \u00e0 ", .bay_pct(1 / length(laws), ","), "."))
	.bay_tell(story, paste0("Prior on the effect: Zellner-Siow (Cauchy) of scale r = ", r_txt$en, ", the \"medium\" of Rouder et al. 2012; the ", length(laws), " shapes equally likely in advance."),
		paste0("A priori sur l'effet : Zellner-Siow (Cauchy) d'\u00e9chelle r = ", r_txt$fr, ", le \u00ab medium \u00bb de Rouder et al. 2012 ; les ", length(laws), " formes \u00e9galement probables d'avance."), "detail")
	if (length(st$n) > 2L) {
		.bay_tell(story, "More than 2 groups: two forms of a difference start equal, all groups may differ, or one group alone stands apart.",
			"Plus de 2 groupes : deux formes de diff\u00e9rence partent \u00e0 \u00e9galit\u00e9, tous les groupes peuvent diff\u00e9rer, ou un seul s'\u00e9carte.", "detail")
	}
	.bay_tell(story, if (isTRUE(sensitivity)) "These values are not tuned on the data; step 5 varies them." else "These values are not tuned on the data; not varied (sensitivity = FALSE).",
		if (isTRUE(sensitivity)) "Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier." else "Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; non vari\u00e9es (sensitivity = FALSE).", "note")

	# Step 5: the Bayes factor, its precision, what may change the strategy.
	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = "m.bay() {KefiR}")
	lb_txt <- list(en = .bay_bf_txt(ans$log_bf, "."), fr = .bay_bf_txt(ans$log_bf, ","))
	inv_en <- if (ans$log_bf < 0) paste0(", i.e. ", .bay_bf_txt(-ans$log_bf, "."), " in favour of no difference") else ""
	inv_fr <- if (ans$log_bf < 0) paste0(", soit ", .bay_bf_txt(-ans$log_bf, ","), " en faveur de l'absence de diff\u00e9rence") else ""
	.bay_tell(story, paste0("Bayes factor: ", lb_txt$en, inv_en, if (abs(ans$log_bf) < log(3)) " (close to 1: the data hardly tell the two apart)." else if (ans$log_bf >= 0) " (the data favour a difference)." else " (the data favour no difference)."),
		paste0("Facteur de Bayes : ", lb_txt$fr, inv_fr, if (abs(ans$log_bf) < log(3)) " (proche de 1 : les donn\u00e9es d\u00e9partagent \u00e0 peine les deux)." else if (ans$log_bf >= 0) " (les donn\u00e9es favorisent une diff\u00e9rence)." else " (les donn\u00e9es favorisent l'absence de diff\u00e9rence)."))
	.bay_tell(story, "The Bayes factor says how many times more probable the data are with a difference than without.",
		"Le facteur de Bayes dit combien de fois les donn\u00e9es sont plus probables avec une diff\u00e9rence que sans.", "note")
	low <- ans$models$ess[!is.na(ans$models$ess)] < 200
	if (any(low) || isTRUE(ans$error_pct > 5)) {
		note(paste0("Numerical precision is low: error on the Bayes factor about ",
				format(signif(ans$error_pct, 2), trim = TRUE), " %."),
			paste0("Pr\u00e9cision num\u00e9rique faible : erreur sur le facteur de Bayes d'environ ",
				gsub(".", ",", format(signif(ans$error_pct, 2), trim = TRUE), fixed = TRUE), " %."))
	} else {
		ep <- .bay_prec_txt(ans$error_pct, "%")
		.bay_tell(story, paste0("Precision of the computation: ", ep$en, " on the Bayes factor, negligible."),
			paste0("Pr\u00e9cision du calcul : ", ep$fr, " sur le facteur de Bayes, n\u00e9gligeable."))
	}
	if (contradiction) {
		.bay_tell(story, paste0("Normal law (means): a difference at ", .bay_pct(p_normal, "."), "; over the shapes (centres): ", .bay_pct(ans$p_h1, "."), "."),
			paste0("Loi normale (moyennes) : une diff\u00e9rence \u00e0 ", .bay_pct(p_normal, ","), " ; sur les formes (centres) : ", .bay_pct(ans$p_h1, ","), "."))
		.bay_tell(story, "Centres and means do not tell the same story; the normal law is implausible here, but it answers about the means.",
			"Centres et moyennes ne racontent pas la m\u00eame chose ; la loi normale est peu plausible ici, mais elle r\u00e9pond sur les moyennes.", "note")
		.bay_tell(story, "Caution added: support capped at weak; the verdict says that the answer depends on what is compared.",
			"Prudence ajout\u00e9e : soutien plafonn\u00e9 \u00e0 faible ; le verdict dit que la r\u00e9ponse d\u00e9pend de ce que l'on compare.", "decide")
	}

	effects <- .bay_effects_row(term = shown, log_bf = ans$log_bf, prior_h1 = prior_h1,
		rscale = rscale, error_pct = ans$error_pct, analogue = "global")
	hyp <- data.frame(
		label_en = c(paste0(input$response_name, " differs according to ", shown),
			paste0("No difference in ", input$response_name, " according to ", shown)),
		label_fr = c(paste0(input$response_name, " diff\u00e8re selon ", shown),
			paste0("Aucune diff\u00e9rence de ", input$response_name, " selon ", shown)),
		posterior = c(ans$p_h1, ans$p_h0),
		prior = c(prior_h1, 1 - prior_h1),
		null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	lead_is_h1 <- ans$p_h1 >= ans$p_h0

	sens <- NULL
	if (isTRUE(sensitivity)) {
		lb <- vapply(rs, function(r) .bay_laws_average(tab, r)$log_bf, numeric(1))
		names(lb) <- format(rs, trim = TRUE)
		# Odds 1:2 to 2:1 around the prior: with 1:4 to 4:1, the leading hypothesis changed
		# in 98 % of the data sets without effect, up to n = 100 per group, and the verdict
		# no longer told anything (agent_41 C, first run; DECISIONS, D4).
		odds_grid <- c(1/2, 1, 2) * prior_h1 / (1 - prior_h1)
		grid <- .bay_sensitivity(lb, prior_odds = odds_grid)
		# Dirichlet: u/3 and 3u, at the central r and odds. Then the noise imposed.
		us <- c(u / 3, 3 * u)
		lb_u <- vapply(us, function(uu) {
			laws_u <- .bay_laws_allowed(st, laws, uu)
			if (length(laws_u) == 0L) return(NA_real_)
			tu <- .bay_laws_table(st, rscale, uu, draws, seed, laws = laws_u)
			.bay_laws_average(tu, rscale)$log_bf
		}, numeric(1))
		names(lb_u) <- format(us, trim = TRUE)
		lb_laws <- ans$log_bf_law
		odds0 <- prior_h1 / (1 - prior_h1)
		# A law imposed against the data is not a plausible prior: it is shown apart and
		# stays out of the range (RELECTURE 2026-11-07, I1).
		all_p <- c(as.vector(grid$grid), if (any(laws %in% c("normal_het", "t_het"))) {
			.bay_post_prob(log_bf = lb_u, prior_odds = odds0)
		})
		all_p <- all_p[!is.na(all_p)]
		if (!lead_is_h1) all_p <- 1 - all_p
		range_lead <- range(all_p)
		sens <- list(grid = grid$grid, by_u = .bay_post_prob(log_bf = lb_u, prior_odds = odds0),
			by_law = .bay_post_prob(log_bf = lb_laws, prior_odds = odds0),
			all_groups = .bay_post_prob(log_bf = .bay_laws_average(tab, rscale, h1_form = "zs")$log_bf,
				prior_odds = odds0),
			r = rs, u = sort(c(us, u)), odds = colnames(grid$grid), range_lead = range_lead,
			# Prior odds that would make the two hypotheses level, against the lead.
			odds_to_reverse = exp(abs(ans$log_bf)),
			verdict = .bay_sensitivity_verdict(max(ans$p_h1, ans$p_h0), range_lead))
	}
	.bay_story_sensitivity(story, sens, hyp, het = any(laws %in% c("normal_het", "t_het")), log_bf = ans$log_bf)

	# Step 6: which groups differ, and by how much (DECISIONS, D11).
	ph <- NULL
	if (isTRUE(posthoc)) {
		ph <- .bay_posthoc_oneway(y, gr, rscale, prior_h1, draws, seed)
		p_het_all <- sum(ans$p_law[intersect(names(ans$p_law), c("normal_het", "t_het"))])
		ph$caution <- isTRUE(p_het_all >= 0.75) || isTRUE(p_t > 0.5)
		.bay_story_posthoc(story, ph, caution = ph$caution)
	} else {
		.bay_story_no_posthoc(story, "off")
	}

	desc <- data.frame(group = st$levels, n = st$n,
		mean = vapply(split(st$y_raw, st$g_raw), mean, numeric(1), USE.NAMES = FALSE),
		median = vapply(split(st$y_raw, st$g_raw), stats::median, numeric(1), USE.NAMES = FALSE),
		SD = vapply(split(st$y_raw, st$g_raw), function(v) {
			if (length(v) > 1L) stats::sd(v) else NA_real_
		}, numeric(1), USE.NAMES = FALSE), stringsAsFactors = FALSE)
	out <- list(
		hypotheses = hyp, effects = effects,
		laws = data.frame(law = names(ans$p_law), posterior = as.numeric(ans$p_law),
			log_bf = as.numeric(ans$log_bf_law[names(ans$p_law)]), stringsAsFactors = FALSE),
		models = ans$models, sensitivity = sens,
		priors = list(prior_h1 = prior_h1, rscale = rscale, u = u, law_prior = "equal",
			nu = "nu - 2 ~ Exp(1)"),
		questions = list(en = paste0("Does ", input$response_name, " differ according to ", shown, "?"),
			fr = paste0(input$response_name, " diff\u00e8re-t-il selon ", shown, " ?")),
		contradiction = contradiction, heavy = heavy, p_tails = p_t, posthoc = ph, detail = isTRUE(detail),
		descriptive_stats = desc, prob_effect = ans$p_h1,
		largest_gap = (max(desc$mean) - min(desc$mean)) / (st$unit * max(abs(st$y_raw))),
		method = "Averaged over sampling laws", plan = input$plan,
		label = if (identical(term, ".g")) paste0(input$response_name, " ~ g") else input$label,
		n = st$N, notes = notes, call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(ans$p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

#' @export
print.mbay <- function(x, detail = isTRUE(x$detail), story = FALSE, ...) {
	told <- isTRUE(story) && length(x$story$en) > 0L
	if (told) .bay_story_replay(x$story)
	if (identical(x$plan, "contingency")) {
		.bay_print_bilan_cont(x, detail = detail, narrated = told)
	} else if (identical(x$plan, "manova")) {
		.bay_print_bilan_manova(x, detail = detail, narrated = told)
	} else if (identical(x$plan, "ordinal")) {
		.bay_print_bilan_ord(x, detail = detail, narrated = told)
	} else if (identical(x$plan, "counts")) {
		.bay_print_bilan_counts(x, detail = detail, narrated = told)
	} else if (x$plan %in% c("factorial", "paired", "repeated", "mixed", "ancova")) {
		.bay_print_bilan_fact(x, detail = detail, narrated = told)
	} else {
		.bay_print_bilan(x, detail = detail, narrated = told)
	}
	invisible(x)
}

# Several factors: normal law with a common variance, likelihoods from BayesFactor
# (sys_bay_factorial.R). prior_h1 and laws do not apply and are said so when given.
.bay_m_factorial <- function(input, call_expr, get_notes, note, prior_h1, rscale, laws, sensitivity,
                             draws, seed, detail, verbose, return, no_prior, no_laws, story = NULL, posthoc = TRUE) {
	rm <- input$plan %in% c("paired", "repeated", "mixed")
	paired <- identical(input$plan, "paired")
	mixed <- identical(input$plan, "mixed")
	if (mixed) {
		# Replicated measures (D21): several measures of one subject in one condition. Random
		# terms: the subject, and the subject crossed with each combination of factors within
		# subject that is a term of the formula, so that the replicates are not counted as
		# independent subjects (pseudo-replication). Rouder et al. (2012), section 7.
		wi <- intersect(input$within, names(input$factors))
		tk <- .bay_term_key(gsub("`", "", input$terms, fixed = TRUE))
		rnd <- list(id = input$id)
		for (k in tk) {
			f <- strsplit(k, ":", fixed = TRUE)[[1]]
			if (length(f) && all(f %in% wi)) {
				rnd[[paste0("id:", k)]] <- interaction(c(list(input$id), as.list(input$factors[f])), drop = TRUE)
			}
		}
		input$id_rand <- as.data.frame(rnd, check.names = FALSE)
	}
	n_dropped <- 0L
	dropped_ids <- character(0)
	if (paired) {
		# A pair needs both measures: subjects with one only are removed (D15).
		pr <- .bay_paired_xy(input$response, input$factors[[1L]], input$id)
		keep <- as.character(input$id) %in% pr$id
		n_dropped <- nlevels(input$id) - length(pr$id)
		dropped_ids <- setdiff(levels(input$id), pr$id)
		if (length(pr$id) < 2L) {
			.exit("Fewer than 2 complete pairs: the paired comparison is impossible.",
				"Moins de 2 paires compl\u00e8tes : la comparaison appari\u00e9e est impossible.")
		}
		if (n_dropped > 0L) {
			note(paste0(n_dropped, " subject(s) without a measure at both levels removed."),
				paste0(n_dropped, " sujet(s) sans mesure aux deux niveaux retir\u00e9(s)."), tell = FALSE)
		}
		input$response <- input$response[keep]
		input$factors <- input$factors[keep, , drop = FALSE]
		input$id <- droplevels(input$id[keep])
		if (any(duplicated(paste(input$id, input$factors[[1L]])))) {
			.exit("Several measures for one subject at one level: the pairs are ambiguous (mixed model, not handled yet).",
				"Plusieurs mesures pour un sujet \u00e0 un m\u00eame niveau : les paires sont ambigu\u00ebs (mod\u00e8le mixte, pas encore trait\u00e9).")
		}
	}
	y <- input$response
	factors <- input$factors
	anc <- identical(input$plan, "ancova")
	cov_names <- if (anc) names(input$covariates) else character(0)
	terms <- .bay_space_terms(gsub("`", "", input$terms, fixed = TRUE), cov_names)
	sp <- .bay_fact_space(terms, cov_names)
	n_mod <- length(sp$models)
	single <- length(terms) == 1L
	overlap <- NULL
	fa <- NULL

	if (anc) {
		overlap <- .bay_story_ancova(story, input, y)
		if (!no_laws && !identical(laws, "normal_eq")) {
			note("With covariates only the normal law is computed; laws is ignored.",
				"Avec des covariables, seule la loi normale est calcul\u00e9e ; laws est ignor\u00e9.")
		}
	} else if (rm) {
		.bay_story_rm(story, input, y, paired, n_dropped, dropped_ids)
		if (!no_laws && !identical(laws, "normal_eq")) {
			note("With subjects only the normal law is computed; laws is ignored.",
				"Avec des sujets, seule la loi normale est calcul\u00e9e ; laws est ignor\u00e9.")
		}
	} else {
	# Step 2 (title said by m.bay): sizes per cell, ties, values far from their cell.
	cell_all <- table(interaction(factors, drop = FALSE))
	chk <- .bay_story_checks(story, y, interaction(factors, drop = TRUE), "cell", "case")
	if (any(cell_all == 0L)) {
		.bay_tell(story, paste0("Empty cells: ", sum(cell_all == 0L), " of ", length(cell_all), "."),
			paste0("Cases vides : ", sum(cell_all == 0L), " sur ", length(cell_all), "."))
	}
	if (all(cell_all <= 1L)) {
		.bay_tell(story, "One observation per cell: the residual of the additive model stands for the spread.",
			"Une observation par case : le r\u00e9sidu du mod\u00e8le additif tient lieu de dispersion.")
	}
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	# Step 3: one shape of the noise only.
	# Variances per cell (D19): one variance per cell is a second law, averaged with the common
	# variance in every Bayes factor (sys_bay_fact_laws.R). Needs 3 values per cell at least.
	het_cells <- FALSE
	cl <- interaction(factors, drop = TRUE)
	use_het <- !single && nlevels(cl) >= 2L && min(table(cl)) >= 3L && (no_laws || "normal_het" %in% laws)
	if (use_het) {
		st_c <- .bay_laws_stats(y, cl)
		use_het <- "normal_het" %in% .bay_laws_allowed(st_c, c("normal_eq", "normal_het"), 1.5)
	}
	fa <- .bay_factorial(input, rscale, sensitivity, as.integer(draws), as.integer(seed), note, prior_h1 = prior_h1, het = use_het)
	.bay_tell(story, "Shape of the noise", "Forme du bruit", "step")
	sds <- tapply(y, cl, stats::sd)
	sd_rng <- function(dec) paste0(gsub(".", dec, format(signif(min(sds), 2), trim = TRUE), fixed = TRUE),
		if (dec == ".") " to " else " \u00e0 ", gsub(".", dec, format(signif(max(sds), 2), trim = TRUE), fixed = TRUE))
	if (isTRUE(fa$het_used)) {
		ph <- fa$p_het
		.bay_step_result(story, if (is.na(fa$p_t)) "normal, common variance or one per cell, averaged" else "normal (common variance or one per cell) and Student t, averaged",
			if (is.na(fa$p_t)) "normale, variance commune ou par case, moyenn\u00e9es" else "normale (variance commune ou par case) et Student, moyenn\u00e9es")
		.bay_tell(story, paste0("Variances per cell: ", if (ph >= 0.75) "unequal probable" else if (ph <= 0.25) "equal probable" else "not settled",
				" (", if (ph >= 0.75) .bay_pct(ph, ".") else if (ph <= 0.25) .bay_pct(1 - ph, ".") else paste("unequal at", .bay_pct(ph, ".")),
				"); standard deviations from ", sd_rng("."), "."),
			paste0("Variances par case : ", if (ph >= 0.75) "in\u00e9gales probables" else if (ph <= 0.25) "\u00e9gales probables" else "non tranch\u00e9es",
				" (", if (ph >= 0.75) .bay_pct(ph, ",") else if (ph <= 0.25) .bay_pct(1 - ph, ",") else paste("in\u00e9gales \u00e0", .bay_pct(ph, ",")),
				") ; \u00e9carts-types de ", sd_rng(","), "."))
		if (ph >= 0.75) {
			.bay_tell(story, "So each answer is computed mostly with one variance per cell: a very dispersed cell no longer makes gaps look larger than they are.",
				"Chaque r\u00e9ponse est donc calcul\u00e9e surtout avec une variance par case : une case tr\u00e8s dispers\u00e9e ne fait plus para\u00eetre les \u00e9carts plus grands qu'ils ne sont.", "decide")
		} else if (ph <= 0.25) {
			.bay_tell(story, "So each answer is computed mostly with a common variance, as usual; one variance per cell keeps its small weight.",
				"Chaque r\u00e9ponse est donc calcul\u00e9e surtout avec une variance commune, comme d'habitude ; une variance par case garde son petit poids.", "decide")
		} else {
			.bay_tell(story, "Both are kept: each answer weighs the common variance and one variance per cell by their probability.",
				"Les deux sont gard\u00e9es : chaque r\u00e9ponse p\u00e8se la variance commune et une variance par case selon leur probabilit\u00e9.", "decide")
		}
		pt <- fa$p_t
		if (is.na(pt)) {
			.bay_tell(story, "heavy tails could not be computed here (too many ties or sampling failed).",
				"les queues lourdes n'ont pas pu \u00eatre calcul\u00e9es ici (trop d'ex-aequo ou tirage en \u00e9chec).", "note")
		} else {
			# D24: a third law, Student t with a common scale, averaged with the two normal laws.
			.bay_tell(story, paste0("Tails: ", if (pt >= 0.75) "heavy probable" else if (pt <= 0.25) "normal probable" else "not settled",
					" (Student t law at ", .bay_pct(pt, "."), "); it is averaged with the normal laws too."),
				paste0("Queues : ", if (pt >= 0.75) "lourdes probables" else if (pt <= 0.25) "normales probables" else "non tranch\u00e9es",
					" (loi de Student \u00e0 ", .bay_pct(pt, ","), ") ; elle est moyenn\u00e9e elle aussi avec les lois normales."))
			.bay_tell(story, "Heavy tails: Student t with a common scale, p(nu) = exp(-(nu - 2)), same priors on the effects (Cauchy of scale r); importance sampling around the mode (DECISIONS, D24).",
				"Queues lourdes : Student \u00e0 \u00e9chelle commune, p(nu) = exp(-(nu - 2)), m\u00eames a priori sur les effets (Cauchy d'\u00e9chelle r) ; \u00e9chantillonnage d'importance autour du mode (DECISIONS, D24).", "detail")
		}
		.bay_tell(story, "One variance per cell: Dirichlet(1.5) on the proportions of precision of the cells (Dablander et al. 2024); effects on the scale of the harmonic mean of the cell variances; importance sampling (DECISIONS, D19).",
			"Une variance par case : Dirichlet(1,5) sur les parts de pr\u00e9cision des cases (Dablander et al. 2024) ; effets \u00e0 l'\u00e9chelle de la moyenne harmonique des variances des cases ; \u00e9chantillonnage d'importance (DECISIONS, D19).", "detail")
	} else {
		.bay_step_result(story, "normal with a common variance, only", "normale \u00e0 variance commune, seule")
		.bay_tell(story, if (isTRUE(fa$het_failed)) "Several factors: the law with one variance per cell could not be computed; the common variance is kept." else "Several factors: only the normal law with a common variance is computed (fewer than 3 values in a cell).",
			if (isTRUE(fa$het_failed)) "Plusieurs facteurs : la loi \u00e0 une variance par case n'a pas pu \u00eatre calcul\u00e9e ; la variance commune est gard\u00e9e." else "Plusieurs facteurs : seule la loi normale \u00e0 variance commune est calcul\u00e9e (moins de 3 valeurs dans une case).")
		if (isTRUE(fa$het_failed) && max(sds) / min(sds) >= 2) {
			het_cells <- TRUE
			.bay_tell(story, paste0("Standard deviations per cell from ", sd_rng("."), ": caution added, support capped at weak."),
				paste0("\u00c9carts-types par case de ", sd_rng(","), " : prudence ajout\u00e9e, soutien plafonn\u00e9 \u00e0 faible."), "decide")
		} else {
			.bay_tell(story, "Heavy tails and unequal variances are not validated with several factors yet.",
				"Queues lourdes et variances in\u00e9gales ne sont pas encore valid\u00e9es \u00e0 plusieurs facteurs.", "note")
		}
	}
	if (!no_laws && any(c("t_eq", "t_het") %in% laws)) {
		note("With several factors the Student t laws are not computed; only the normal ones.",
			"\u00c0 plusieurs facteurs, les lois de Student ne sont pas calcul\u00e9es ; seulement les normales.")
	}
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(story, paste0(far, if (far == 1L) " value lies" else " values lie", " beyond 3 interquartile ranges (step 2): read the answer with caution."),
			paste0(far, if (far == 1L) " valeur" else " valeurs", " au-del\u00e0 de 3 \u00e9carts interquartiles (\u00e9tape 2) : lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(story, "No far value puts the normal law in doubt at first sight.",
			"Aucune valeur \u00e9loign\u00e9e ne met d'embl\u00e9e la loi normale en doute.", "decide")
	}
	}

	# Step 4: the priors, fixed by the method (DECISIONS, D2).
	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step",
		fun = if (paired) "ttestBF(paired = TRUE) {BayesFactor}" else NULL)
	lp <- if (single) log(c(1 - prior_h1, prior_h1)) else .bay_fact_prior(sp, "hop")
	incl <- vapply(seq_along(terms), function(i) sum(exp(lp[vapply(sp$models, function(m) i %in% m, logical(1))])), numeric(1))
	r_txt <- .bay_num_txt(rscale)
	if (single) {
		.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference"), paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence"))
		.bay_tell(story, paste0("Expected size r = ", r_txt$en, if (paired) " on dz (gap within subject over its standard deviation)." else "."),
			paste0("Taille attendue r = ", r_txt$fr, if (paired) " sur dz (\u00e9cart au sein du sujet rapport\u00e9 \u00e0 son \u00e9cart-type)." else "."))
	} else {
		.bay_step_result(story, paste0("no difference at all: ", .bay_pct(exp(lp[1L]), ".")), paste0("aucune diff\u00e9rence du tout : ", .bay_pct(exp(lp[1L]), ",")))
		.bay_tell(story, paste0("Chance of each term in advance: ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = "."), collapse = ", "), "."),
			paste0("Chance de chaque terme d'avance : ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = ","), collapse = ", "), "."))
		.bay_tell(story, paste0("Expected size r = ", r_txt$en, " per term."), paste0("Taille attendue r = ", r_txt$fr, " par terme."))
		.bay_tell(story, "these chances do not add up; the more terms, the less likely each one, which guards against false interactions.",
			"ces chances ne s'additionnent pas ; plus il y a de termes, moins chacun est probable, ce qui prot\u00e8ge des fausses interactions.", "note")
		.bay_tell(story, "Prior on the models: one proportion ~ Beta(1,1) per order of interaction (HOP, Taylor-Rodriguez et al. 2016); strong heredity. Effects: Zellner-Siow (Cauchy) of scale r per term (Rouder et al. 2012).",
			"A priori sur les mod\u00e8les : une proportion ~ Beta(1,1) par ordre d'interaction (HOP, Taylor-Rodriguez et al. 2016) ; h\u00e9r\u00e9dit\u00e9 forte. Effets : Zellner-Siow (Cauchy) d'\u00e9chelle r par terme (Rouder et al. 2012).", "detail")
	}
	if (anc) {
		.bay_tell(story, paste0(paste(cov_names, collapse = ", "), if (length(cov_names) == 1L) " is" else " are", " in every model; each answer is read against the covariates alone."),
			paste0(paste(cov_names, collapse = ", "), if (length(cov_names) == 1L) " est" else " sont", " dans tous les mod\u00e8les ; chaque r\u00e9ponse se lit contre les covariables seules."))
		.bay_tell(story, "Covariate slope: Cauchy prior of scale sqrt(2)/4 (rscaleCont, the \"medium\" of BayesFactor), shared with its interactions with a factor, varied with r in the sensitivity.",
			"Pente de la covariable : a priori de Cauchy d'\u00e9chelle sqrt(2)/4 (rscaleCont, le \u00ab medium \u00bb de BayesFactor), partag\u00e9 avec ses interactions avec un facteur, vari\u00e9 avec r dans la sensibilit\u00e9.", "detail")
	}
	if (rm) {
		.bay_tell(story, if (paired) "Each subject is compared with itself: the gap within subject is tested." else "Each subject keeps its own level, in every model, \"no difference\" included.",
			if (paired) "Chaque sujet est compar\u00e9 \u00e0 lui-m\u00eame : c'est l'\u00e9cart au sein du sujet qui est test\u00e9." else "Chaque sujet garde son niveau propre, dans tous les mod\u00e8les, \u00ab aucune diff\u00e9rence \u00bb compris.")
		.bay_tell(story, if (paired) "Paired: Bayesian paired t test (Rouder et al. 2009), r = 0.5 on dz, as for a random subject effect (equal to lmBF(y ~ time + id) within 0.2 on the log scale)." else "Subjects: random effect of scale 1 (Rouder et al. 2012); each answer is read against the model with the subjects only.",
			if (paired) "Appari\u00e9 : test t appari\u00e9 bay\u00e9sien (Rouder et al. 2009), r = 0,5 sur dz, comme pour un effet sujet al\u00e9atoire (\u00e9gal \u00e0 lmBF(y ~ temps + id) \u00e0 0,2 pr\u00e8s en \u00e9chelle log)." else "Sujets : effet al\u00e9atoire d'\u00e9chelle 1 (Rouder et al. 2012) ; chaque r\u00e9ponse se lit contre le mod\u00e8le avec les sujets seuls.", "detail")
	}
	.bay_tell(story, "These values are not tuned on the data; step 5 varies them.",
		"Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")
	if (!single && !any(sp$order > 1L)) {
		note("No interaction is tested (additive formula or x, g form); write y ~ A * B to test them.",
			"Aucune interaction n'est test\u00e9e (formule additive ou forme x, g) ; \u00e9crire y ~ A * B pour les tester.")
	}
	if (!no_prior && !single) {
		note("prior_h1 has no role with several factors: the prior on the models is fixed by the method.",
			"prior_h1 est sans effet \u00e0 plusieurs facteurs : la loi sur les mod\u00e8les est fix\u00e9e par la m\u00e9thode.")
	}

	# Step 5: solidity of the answer: Bayes factors of the models, precision, sensitivity.
	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step",
		fun = if (paired) "ttestBF(paired = TRUE) {BayesFactor}" else "lmBF() {BayesFactor}")
	if (n_mod > 50L) {
		.bay_tell(story, paste0(n_mod, " models, this can take a few minutes."),
			paste0(n_mod, " mod\u00e8les, cela peut prendre quelques minutes."), "note")
	}
	if (is.null(fa)) fa <- .bay_factorial(input, rscale, sensitivity, as.integer(draws), as.integer(seed), note, prior_h1 = prior_h1)
	if (single) {
		lb1 <- fa$log_bf_models[2L, if (isTRUE(sensitivity)) 2L else 1L]
		.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb1, "."), if (lb1 >= 0) " (the data favour a difference)." else " (the data favour no difference)."),
			paste0("Facteur de Bayes : ", .bay_bf_txt(lb1, ","), if (lb1 >= 0) " (les donn\u00e9es favorisent une diff\u00e9rence)." else " (les donn\u00e9es favorisent l'absence de diff\u00e9rence)."))
	} else {
		.bay_tell(story, paste0(n_mod, " models compared (combinations of terms)."), paste0(n_mod, " mod\u00e8les compar\u00e9s (combinaisons de termes)."))
	}
	if (paired) {
		.bay_tell(story, "Precision: exact computation, without draws.", "Pr\u00e9cision : calcul exact, sans tirage.")
	} else if (isTRUE(fa$spread <= 0.02)) {
		sp_txt <- .bay_prec_txt(100 * fa$spread, if (.bay_lang_fr()) "point" else "point")
		.bay_tell(story, paste0("Precision of the computation: ", sp_txt$en, " on the probabilities, negligible."),
			paste0("Pr\u00e9cision du calcul : ", sp_txt$fr, " sur les probabilit\u00e9s, n\u00e9gligeable."))
		.bay_tell(story, "Precision: the central computation redone with another seed, largest change of a probability.",
			"Pr\u00e9cision : le calcul central refait avec une autre graine, plus grand \u00e9cart d'une probabilit\u00e9.", "detail")
	}
	sens <- fa$sensitivity
	if (is.null(sens)) {
		.bay_step_result(story, "sensitivity not computed", "sensibilit\u00e9 non calcul\u00e9e")
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
		.bay_tell(story, "The answers are not tested against other starting points.",
			"Les r\u00e9ponses ne sont pas \u00e9prouv\u00e9es contre d'autres points de d\u00e9part.", "decide")
	} else {
		r1 <- .bay_num_txt(min(sens$r)); r2 <- .bay_num_txt(max(sens$r))
		moved <- which(sens$verdict_terms == "changes")
		if (single) {
			rg <- sens$range_lead
			stb <- .bay_stab_txt(sens$verdict)
			.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb1, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb1, ","), ", ", stb$fr))
			.bay_tell(story, paste0("Sensitivity: the leading answer ", .bay_range_txt(rg)$en, "."),
				paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate ", .bay_range_txt(rg)$fr, "."))
		} else if (length(moved)) {
			.bay_step_result(story, paste0(paste(terms[moved], collapse = ", "), if (length(moved) == 1L) " changes side with the prior" else " change side with the prior"),
				paste0(paste(terms[moved], collapse = ", "), if (length(moved) == 1L) " change de c\u00f4t\u00e9 selon l'a priori" else " changent de c\u00f4t\u00e9 selon l'a priori"))
		} else {
			.bay_step_result(story, "no term changes side", "aucun terme ne change de c\u00f4t\u00e9")
		}
		.bay_tell(story, paste0("r goes from ", r1$en, " to ", r2$en, if (anc) " (rscaleCont with it)" else "", if (single) ", the starting chances from 1 against 2 to 2 against 1." else ", with three priors on the models."),
			paste0("r va de ", r1$fr, " \u00e0 ", r2$fr, if (anc) " (rscaleCont avec lui)" else "", if (single) ", les chances de d\u00e9part de 1 contre 2 \u00e0 2 contre 1." else ", avec trois lois sur les mod\u00e8les."), "note")
		if (!single) {
			.bay_tell(story, "Priors on the models in the sensitivity: HOP, 1/2 per term (Chipman 1996), uniform.",
				"Lois sur les mod\u00e8les dans la sensibilit\u00e9 : HOP, 1/2 par terme (Chipman 1996), uniforme.", "detail")
		}
		if (length(moved) && !single) {
			rr <- sens$range_terms[moved, , drop = FALSE]
			for (m in seq_along(moved)) {
				.bay_tell(story, paste0(terms[moved[m]], " changes side: ", .bay_pct(rr[m, 1L], "."), " to ", .bay_pct(rr[m, 2L], "."), "."),
					paste0(terms[moved[m]], " change de c\u00f4t\u00e9 : ", .bay_pct(rr[m, 1L], ","), " \u00e0 ", .bay_pct(rr[m, 2L], ","), "."))
			}
			.bay_tell(story, paste0("Caution added: ", paste(terms[moved], collapse = ", "), if (length(moved) > 1L) " will be called not settled." else " will be called not settled."),
				paste0("Prudence ajout\u00e9e : ", paste(terms[moved], collapse = ", "), if (length(moved) > 1L) " seront dits non tranch\u00e9s." else " sera dit non tranch\u00e9."), "decide")
		} else if (single && identical(sens$verdict, "changes")) {
			.bay_tell(story, "Caution added: the answer can reverse; it will be called not settled.",
				"Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e.", "decide")
		} else {
			.bay_tell(story, if (single) "The answer does not change side." else "No term changes side.",
				if (single) "La r\u00e9ponse ne change pas de c\u00f4t\u00e9." else "Aucun terme ne change de c\u00f4t\u00e9.", "decide")
		}
		if (identical(sens$verdict, "changes") && !single) {
			.bay_tell(story, "The global answer can reverse: the bilan will not assert it.",
				"La r\u00e9ponse globale peut s'inverser : le bilan ne l'affirmera pas.", "decide")
		}
	}

	# Step 6 for one term with subjects: which conditions differ, and by how much (D11, D15).
	ph <- NULL
	if (anc && length(factors) == 1L && isTRUE(posthoc)) {
		g1 <- factors[[1L]]
		fn <- names(factors)[1L]
		int <- terms[grepl(":", terms, fixed = TRUE)]
		p_int <- if (length(int)) fa$effects$p_h1[match(.bay_term_key(int), .bay_term_key(fa$effects$term))] else numeric(0)
		covc <- as.data.frame(lapply(input$covariates, function(v) v - mean(v)), check.names = FALSE)
		if (length(int) && isTRUE(max(p_int, na.rm = TRUE) >= 0.5)) {
			top <- which.max(p_int)
			xn <- setdiff(strsplit(int[top], ":", fixed = TRUE)[[1]], fn)[1L]
			ph <- .bay_posthoc_slopes(y, g1, covc[[xn]], mean(input$covariates[[xn]]), covc[setdiff(names(covc), xn)], rscale, draws, seed)
			ph$xname <- xn
			.bay_story_slopes(story, ph, xn, int[top], p_int[top])
		} else {
			ph <- .bay_posthoc_oneway(y, g1, rscale, prior_h1, draws, seed, cov = covc)
			ph$caution <- FALSE
			ph$caution_overlap <- isTRUE(any(fa$overlap$poor))
			.bay_story_posthoc(story, ph, caution = ph$caution)
			if (length(int)) {
				.bay_tell(story, paste0("Slopes that differ are not the most probable (", int[1L], " ", .bay_pct(p_int[1L], "."), "): one adjusted gap per pair."),
					paste0("Pentes diff\u00e9rentes peu probables (", int[1L], " ", .bay_pct(p_int[1L], ","), ") : un \u00e9cart ajust\u00e9 par paire."), "decide")
			}
		}
	}
	if (rm && single && isTRUE(posthoc)) {
		g1 <- factors[[1L]]
		ph <- if (paired) .bay_posthoc_paired(y, g1, input$id, rscale, draws, seed) else
			.bay_posthoc_oneway(y, g1, rscale, prior_h1, draws, seed, id = if (mixed) input$id_rand else input$id)
		ph$caution <- FALSE
		.bay_story_posthoc(story, ph, subjects = TRUE)
	}
	# Several factors (PERIMETRE section 7): the probable term of highest order, its cells or its
	# levels, the other factors, covariates and subjects as nuisance terms.
	if (is.null(ph) && isTRUE(posthoc) && !single && !paired) {
		covc <- if (anc) as.data.frame(lapply(input$covariates, function(v) v - mean(v)), check.names = FALSE) else NULL
		# Undecided interactions, with the words of the bilan: their gaps per level are added.
		ta <- .bay_term_answers(list(effects = fa$effects, sensitivity = fa$sensitivity, overlap = fa$overlap,
			het_cells = exists("het_cells", inherits = FALSE) && isTRUE(het_cells)))
		und <- .bay_term_key(fa$effects$term)[grepl(":", fa$effects$term, fixed = TRUE) & fa$effects$p_h1 < 0.5 &
			(ta$open | fa$effects$p_h1 > 0.25)]
		ph <- .bay_posthoc_factorial(y, factors, fa$effects, rscale, 0.5, draws, seed, id = if (mixed) input$id_rand else if (rm) input$id else NULL, cov = covc,
			undecided = und)
		ph$caution <- isTRUE(fa$het_used) && isTRUE(fa$p_het >= 0.75)
		ph$open_terms <- if (!is.null(fa$sensitivity)) fa$effects$term[fa$sensitivity$verdict_terms == "changes"] else character(0)
		.bay_story_posthoc_fact(story, ph, caution = ph$caution)
	}
	if (is.null(ph)) .bay_story_no_posthoc(story, if (isTRUE(posthoc)) "factors" else "off")

	resp <- input$response_name
	fac_txt <- paste(names(input$factors), collapse = ", ")
	p_any <- 1 - fa$p_null
	hyp <- data.frame(
		label_en = c(paste0(if (single || (anc && length(factors) == 1L)) paste0(resp, " differs according to ", fac_txt) else paste0(resp, " differs according to at least one factor"),
				if (anc) paste0(", at equal ", paste(cov_names, collapse = ", ")) else ""),
			paste0("No difference in ", resp, " according to ", fac_txt)),
		label_fr = c(paste0(if (single || (anc && length(factors) == 1L)) paste0(resp, " diff\u00e8re selon ", fac_txt) else paste0(resp, " diff\u00e8re selon au moins un facteur"),
				if (anc) paste0(", \u00e0 m\u00eame ", paste(cov_names, collapse = ", ")) else ""),
			paste0("Aucune diff\u00e9rence de ", resp, " selon ", fac_txt)),
		posterior = c(p_any, fa$p_null), prior = c(1 - fa$prior_null, fa$prior_null),
		null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(
		hypotheses = hyp, effects = fa$effects,
		laws = if (isTRUE(fa$het_used)) { pt_ <- if (is.na(fa$p_t)) 0 else fa$p_t; data.frame(law = c("normal_eq", "normal_het", if (!is.na(fa$p_t)) "t_eq"), posterior = c(1 - fa$p_het - pt_, fa$p_het, if (!is.na(fa$p_t)) fa$p_t), stringsAsFactors = FALSE) } else NULL,
		models = data.frame(model = vapply(fa$space$models, function(m) {
				if (length(m)) paste(fa$space$terms[m], collapse = " + ") else "(null)"
			}, character(1)), log_bf = fa$log_bf_models[, if (isTRUE(sensitivity)) 2L else 1L],
			log_prior = fa$log_prior, stringsAsFactors = FALSE),
		n_models = length(fa$space$models),
		sensitivity = fa$sensitivity,
		priors = list(rscale = rscale, model_prior = if (single) "prior_h1" else "HOP(1,1)", law = if (isTRUE(fa$het_used)) c("normal_eq", "normal_het", if (!is.na(fa$p_t)) "t_eq") else "normal_eq", prior_h1 = if (single) prior_h1 else NULL),
		posthoc = ph, n_subjects = fa$n_subjects,
		questions = .bay_fact_questions(fa$effects$term, resp, adj = cov_names),
		overlap = fa$overlap, cov_centres = fa$cov_centres,
		cov_by_group = if (anc) sapply(input$covariates, function(v) tapply(v, interaction(input$factors, drop = TRUE), mean)) else NULL,
		het_cells = if (exists("het_cells", inherits = FALSE)) het_cells else FALSE,
		het_used = isTRUE(fa$het_used), p_t = fa$p_t,
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail),
		descriptive_stats = fa$means, cell_means = local({
			cl <- interaction(input$factors, sep = " x ", drop = TRUE)
			data.frame(cell = levels(cl), n = as.integer(table(cl)), mean = as.numeric(tapply(input$response, cl, mean)), stringsAsFactors = FALSE)
		}), sd_within = fa$sd_within, sd_is_within = fa$sd_is_within, balanced = fa$balanced,
		prob_effect = p_any,
		method = if (anc) "Analysis of covariance, normal law, covariates centred (BayesFactor)" else if (paired) "Paired data, Bayesian paired t test (BayesFactor)" else if (mixed) "Replicated repeated measures, normal law, random subject and subject x condition effects (BayesFactor)" else if (rm) "Repeated measures, normal law, random subject effect (BayesFactor)" else "Several factors, normal law with a common variance (BayesFactor)", plan = input$plan,
		label = input$label, n = length(input$response), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_fact(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_any)
	if (isTRUE(verbose)) invisible(out) else out
}
