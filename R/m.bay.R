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
#' `laws` have no role there. The other plans of [m.test()] (covariates, paired, repeated
#' measures, mixed, MANOVA, discrete data) stop with an explicit message until they are
#' validated.
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
#' @param verbose If `TRUE`, prints the bilan.
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
#' @importFrom BayesFactor lmBF extractBF
#' @export
m.bay <- function(x = NULL, g = NULL, data = NULL, formula = NULL,
                  paired = FALSE, id = NULL, wt = NULL, within = NULL, between = NULL,
                  alpha = NULL, control = NULL, verbose = TRUE, plot = FALSE,
                  return = TRUE, boot = NULL, boot_type = NULL, iter = NULL, conf = NULL,
                  maxcat = 50, silent = NULL, code = FALSE, debug = FALSE,
                  prior_h1 = 1/2, rscale = 1/2, u = 1.5,
                  laws = c("normal_eq", "normal_het", "t_eq", "t_het"),
                  sensitivity = TRUE, draws = 4000L, seed = 1L, detail = FALSE) {
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
	note <- function(en, fr) {
		notes$en <<- c(notes$en, en)
		notes$fr <<- c(notes$fr, fr)
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

	# Plans not validated yet: an explicit refusal, never a wrong answer.
	if (!input$plan %in% c("oneway", "factorial")) {
		plan_en <- c(factorial = "several factors", ancova = "covariates", paired = "paired data",
			repeated = "repeated measures", mixed = "mixed models", manova = "MANOVA")[input$plan]
		plan_fr <- c(factorial = "plusieurs facteurs", ancova = "les covariables", paired = "les donn\u00e9es appari\u00e9es",
			repeated = "les mesures r\u00e9p\u00e9t\u00e9es", mixed = "les mod\u00e8les mixtes", manova = "la MANOVA")[input$plan]
		.exit(paste0("m.bay() does not handle ", plan_en, " yet: this plan is not validated (use m.test())."),
			paste0("m.bay() ne traite pas encore ", plan_fr, " : ce plan n'est pas valid\u00e9 (utiliser m.test())."))
	}
	y <- input$response
	if (any(is.infinite(y))) {
		.exit("y contains infinite values.", "y contient des valeurs infinies.")
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
		note("The response looks discrete: the normal and t laws are only an approximation here.",
			"La r\u00e9ponse semble discr\u00e8te : les lois normale et t ne sont ici qu'une approximation.")
	}

	if (identical(input$plan, "factorial")) {
		return(.bay_m_factorial(input, call_expr, function() notes, note, prior_h1, rscale, laws, sensitivity,
			draws, seed, detail, verbose, return, !("prior_h1" %in% names(call_expr)), !("laws" %in% names(call_expr))))
	}

	term <- names(input$factors)[1L]
	shown <- if (identical(term, ".g")) "g" else term
	gr <- input$factors[[1L]]
	st <- .bay_laws_stats(y, gr)
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
	if (length(laws_kept) < length(laws)) {
		left <- setdiff(laws, laws_kept)
		note(paste0("Repeated values: for these data the integral of some laws diverges, they were left out: ",
				paste(.bay_law_labels()$en[left], collapse = "; "), ". A law for discrete data is not validated yet."),
			paste0("Valeurs r\u00e9p\u00e9t\u00e9es : pour ces donn\u00e9es l'int\u00e9grale de certaines lois diverge, elles ont \u00e9t\u00e9 \u00e9cart\u00e9es : ",
				paste(.bay_law_labels()$fr[left], collapse = " ; "), ". Une loi pour donn\u00e9es discr\u00e8tes n'est pas encore valid\u00e9e."))
	}
	laws <- laws_kept
	draws <- as.integer(draws)
	seed <- as.integer(seed)
	rs <- if (isTRUE(sensitivity)) c(rscale / 2, rscale, rscale * 2) else rscale
	tab <- .bay_laws_table(st, rs, u, draws, seed, laws = laws)
	ans <- .bay_laws_average(tab, rscale, prior_h1 = prior_h1)
	if (length(ans$laws_dropped)) {
		note(paste0("Law(s) left out because the computation is not reliable here, so the comparison is incomplete and the answer may change: ",
				paste(.bay_law_labels()$en[ans$laws_dropped], collapse = "; "), ". Try a larger draws."),
			paste0("Loi(s) \u00e9cart\u00e9e(s) car le calcul n'est pas fiable ici : la comparaison est incompl\u00e8te et la r\u00e9ponse peut changer : ",
				paste(.bay_law_labels()$fr[ans$laws_dropped], collapse = " ; "), ". Essayer un draws plus grand."))
	}
	low <- ans$models$ess[!is.na(ans$models$ess)] < 200
	if (any(low) || isTRUE(ans$error_pct > 5)) {
		note(paste0("Numerical precision is low: error on the Bayes factor about ",
				format(signif(ans$error_pct, 2), trim = TRUE), " %."),
			paste0("Pr\u00e9cision num\u00e9rique faible : erreur sur le facteur de Bayes d'environ ",
				gsub(".", ",", format(signif(ans$error_pct, 2), trim = TRUE), fixed = TRUE), " %."))
	}

	# Under heavy tails the groups are compared through their centres, which are not the
	# means of a skewed group: equal means, different shapes can read as a difference
	# (agent_43: one lognormal group, n = 400, P > 0.95 in 67 % of the data sets). When the
	# normal law and the answer fall on opposite sides, the bilan says so on its first
	# lines and caps the qualifier. Thresholds and false alarms: agent_46, VALIDATION.
	p_t <- sum(ans$p_law[intersect(names(ans$p_law), c("t_eq", "t_het"))])
	normal_kept <- intersect(names(ans$p_law), c("normal_eq", "normal_het"))
	contradiction <- FALSE
	if (length(normal_kept)) {
		best_normal <- normal_kept[which.max(ans$p_law[normal_kept])]
		p_normal <- .bay_post_prob(log_bf = ans$log_bf_law[best_normal], prior_odds = prior_h1 / (1 - prior_h1))
		contradiction <- isTRUE((ans$p_h1 > 0.75 && p_normal < 0.5) || (ans$p_h1 < 0.25 && p_normal > 0.5))
	}
	heavy <- isTRUE(p_t >= .bay_heavy_threshold())
	if (heavy && !contradiction) {
		note("The noise looks heavy-tailed: the comparison bears on the centres of the groups, not exactly on their means.",
			"Le bruit semble \u00e0 queues lourdes : la comparaison porte sur le centre des groupes, pas exactement sur leur moyenne.")
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
		stringsAsFactors = FALSE)
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
		contradiction = contradiction, heavy = heavy, detail = isTRUE(detail),
		descriptive_stats = desc, prob_effect = ans$p_h1,
		largest_gap = (max(desc$mean) - min(desc$mean)) / (st$unit * max(abs(st$y_raw))),
		method = "Averaged over sampling laws", plan = input$plan,
		label = if (identical(term, ".g")) paste0(input$response_name, " ~ g") else input$label,
		n = st$N, notes = notes, call = call_expr)
	class(out) <- "mbay"
	if (isTRUE(verbose)) .bay_print_bilan(out, detail = detail)
	if (isFALSE(return)) return(ans$p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

#' @export
print.mbay <- function(x, detail = isTRUE(x$detail), ...) {
	if (identical(x$plan, "factorial")) .bay_print_bilan_fact(x, detail = detail) else .bay_print_bilan(x, detail = detail)
	invisible(x)
}

# Several factors: normal law with a common variance, likelihoods from BayesFactor
# (sys_bay_factorial.R). prior_h1 and laws do not apply and are said so when given.
.bay_m_factorial <- function(input, call_expr, get_notes, note, prior_h1, rscale, laws, sensitivity,
                             draws, seed, detail, verbose, return, no_prior, no_laws) {
	if (!no_prior) {
		note("prior_h1 has no role with several factors: the prior on the models is fixed by the method.",
			"prior_h1 est sans effet \u00e0 plusieurs facteurs : la loi sur les mod\u00e8les est fix\u00e9e par la m\u00e9thode.")
	}
	if (!no_laws && !identical(laws, "normal_eq")) {
		note("With several factors only the normal law with a common variance is computed; laws is ignored.",
			"\u00c0 plusieurs facteurs, seule la loi normale \u00e0 variance commune est calcul\u00e9e ; laws est ignor\u00e9.")
	}
	if (!any(grepl(":", input$terms, fixed = TRUE))) {
		note("No interaction is tested (additive formula or x, g form); write y ~ A * B to test them.",
			"Aucune interaction n'est test\u00e9e (formule additive ou forme x, g) ; \u00e9crire y ~ A * B pour les tester.")
	}
	n_mod <- length(.bay_fact_space(gsub("`", "", input$terms, fixed = TRUE))$models)
	if (isTRUE(verbose) && n_mod > 50L) {
		message(.msg(paste0("Comparing ", n_mod, " models: this can take a few minutes."),
			paste0("Comparaison de ", n_mod, " mod\u00e8les : cela peut prendre quelques minutes.")))
	}
	fa <- .bay_factorial(input, rscale, sensitivity, as.integer(draws), as.integer(seed), note)
	resp <- input$response_name
	fac_txt <- paste(names(input$factors), collapse = ", ")
	p_any <- 1 - fa$p_null
	hyp <- data.frame(
		label_en = c(paste0(resp, " differs according to at least one factor"),
			paste0("No difference in ", resp, " according to ", fac_txt)),
		label_fr = c(paste0(resp, " diff\u00e8re selon au moins un facteur"),
			paste0("Aucune diff\u00e9rence de ", resp, " selon ", fac_txt)),
		posterior = c(p_any, fa$p_null), prior = c(1 - fa$prior_null, fa$prior_null),
		null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(
		hypotheses = hyp, effects = fa$effects, laws = NULL,
		models = data.frame(model = vapply(fa$space$models, function(m) {
				if (length(m)) paste(fa$space$terms[m], collapse = " + ") else "(null)"
			}, character(1)), log_bf = fa$log_bf_models[, if (isTRUE(sensitivity)) 2L else 1L],
			log_prior = fa$log_prior, stringsAsFactors = FALSE),
		n_models = length(fa$space$models),
		sensitivity = fa$sensitivity,
		priors = list(rscale = rscale, model_prior = "HOP(1,1)", law = "normal_eq"),
		questions = .bay_fact_questions(fa$effects$term, resp),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail),
		descriptive_stats = fa$means, sd_within = fa$sd_within, sd_is_within = fa$sd_is_within, balanced = fa$balanced,
		prob_effect = p_any,
		method = "Several factors, normal law with a common variance (BayesFactor)", plan = input$plan,
		label = input$label, n = length(input$response), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	if (isTRUE(verbose)) .bay_print_bilan_fact(out, detail = detail)
	if (isFALSE(return)) return(p_any)
	if (isTRUE(verbose)) invisible(out) else out
}
