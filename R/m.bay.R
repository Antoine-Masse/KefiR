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
#' and is poorly seen. The other plans of [m.test()] (covariates, paired, repeated
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
#' @param verbose If `TRUE`, tells the analysis step by step as [m.test()] does: numbered
#'   steps printed as the computation goes, "==>" for what was found, "-->" for what the
#'   function decides from it (a law left out, a lead that the priors can reverse), the
#'   function and package in brackets; the last step is the bilan. If `FALSE`, nothing.
#'   `print(result, story = TRUE)` replays the steps; `print(result)` gives the bilan alone.
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
	# The story of the call (verbose = TRUE): every note is also said at the step where it
	# arises, unless tell = FALSE (the step then says it in its own words).
	story <- .bay_story(verbose)
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

	# Step 1: the plan recognised, or refused and why.
	.bay_story_plan(story, input)

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
	.bay_tell(story, "Checks on the data.", "Contr\u00f4les des donn\u00e9es.", "step")
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
			draws, seed, detail, verbose, return, !("prior_h1" %in% names(call_expr)), !("laws" %in% names(call_expr)),
			story))
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
		.bay_tell(story, paste0("Change of strategy: law(s) left out because, with these repeated values, their integral diverges: ",
				paste(lab$en[left], collapse = "; "), ". The comparison goes on with the others."),
			paste0("Changement de strat\u00e9gie : loi(s) \u00e9cart\u00e9e(s) car, avec ces valeurs r\u00e9p\u00e9t\u00e9es, leur int\u00e9grale diverge : ",
				paste(lab$fr[left], collapse = " ; "), ". La comparaison continue avec les autres."), "decide")
	} else if (chk$extremes[["mild"]] > 0L) {
		.bay_tell(story, paste0("No value is removed: the ", length(laws_kept), " shapes of the noise stay in the running, and the heavy-tailed ones are there for such values."),
			paste0("Aucune valeur n'est retir\u00e9e : les ", length(laws_kept), " formes du bruit restent en lice, et celles \u00e0 queues lourdes sont l\u00e0 pour ce genre de valeurs."), "decide")
	} else {
		.bay_tell(story, paste0("Nothing to remove: the ", length(laws_kept), " shapes of the noise stay in the running."),
			paste0("Rien \u00e0 retirer : les ", length(laws_kept), " formes du bruit restent en lice."), "decide")
	}
	laws <- laws_kept
	draws <- as.integer(draws)
	seed <- as.integer(seed)

	# Step 3: the priors, fixed before the computation.
	.bay_tell(story, "Starting point (priors), set before seeing the data.", "Point de d\u00e9part (a priori), fix\u00e9 avant de voir les donn\u00e9es.", "step")
	r_txt <- .bay_num_txt(rscale)
	user_prior <- "prior_h1" %in% names(call_expr)
	user_r <- "rscale" %in% names(call_expr)
	.bay_tell(story, paste0(.bay_pct(prior_h1, "."), " chance that a difference exists", if (user_prior) " (given by the user)" else " (even odds: the method takes no side)",
			"; expected size r = ", r_txt$en, if (user_r) " (given by the user)" else " (\"medium\" scale of Rouder et al. 2012)",
			"; each of the ", length(laws), " shapes of the noise starts at ", .bay_pct(1 / length(laws), "."), "."),
		paste0(.bay_pct(prior_h1, ","), " de chances qu'une diff\u00e9rence existe", if (user_prior) " (donn\u00e9 par l'utilisateur)" else " (pile ou face : la m\u00e9thode ne prend pas parti)",
			" ; taille attendue r = ", r_txt$fr, if (user_r) " (donn\u00e9e par l'utilisateur)" else " (\u00e9chelle \u00ab medium \u00bb de Rouder et al. 2012)",
			" ; chacune des ", length(laws), " formes du bruit part \u00e0 ", .bay_pct(1 / length(laws), ","), "."))
	if (length(st$n) > 2L) {
		.bay_tell(story, "With more than 2 groups, two forms of a difference start equal: all groups may differ, or one group alone stands apart.",
			"Avec plus de 2 groupes, deux formes de diff\u00e9rence partent \u00e0 \u00e9galit\u00e9 : tous les groupes peuvent diff\u00e9rer, ou un seul groupe s'\u00e9carte.")
	}
	.bay_tell(story, if (isTRUE(sensitivity)) "These values are not tuned on the data; step 5 varies them to see whether the answer holds." else
			"These values are not tuned on the data; their influence is not checked (sensitivity = FALSE).",
		if (isTRUE(sensitivity)) "Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier pour voir si la r\u00e9ponse tient." else
			"Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; leur influence n'est pas contr\u00f4l\u00e9e (sensitivity = FALSE).", "decide")

	# Step 4: the shapes of the noise, compared and not chosen.
	.bay_tell(story, "Comparison of the shapes of the noise (how values scatter around the centre of their group: normal or heavy tails, equal variances or not).",
		"Comparaison des formes du bruit (la fa\u00e7on dont les valeurs s'\u00e9cartent du centre de leur groupe : normale ou \u00e0 queues lourdes, variances \u00e9gales ou non).", "step",
		fun = c("m.bay() {KefiR}, importance sampling", "m.bay() {KefiR}, \u00e9chantillonnage d'importance"))
	rs <- if (isTRUE(sensitivity)) c(rscale / 2, rscale, rscale * 2) else rscale
	tab <- .bay_laws_table(st, rs, u, draws, seed, laws = laws)
	ans <- .bay_laws_average(tab, rscale, prior_h1 = prior_h1)
	if (length(ans$laws_dropped)) {
		note(paste0("Law(s) left out because the computation is not reliable here, so the comparison is incomplete and the answer may change: ",
				paste(lab$en[ans$laws_dropped], collapse = "; "), ". Try a larger draws."),
			paste0("Loi(s) \u00e9cart\u00e9e(s) car le calcul n'est pas fiable ici : la comparaison est incompl\u00e8te et la r\u00e9ponse peut changer : ",
				paste(lab$fr[ans$laws_dropped], collapse = " ; "), ". Essayer un draws plus grand."), tell = FALSE)
		.bay_tell(story, paste0("The importance sampling is too poor for: ", paste(lab$en[ans$laws_dropped], collapse = "; "),
				" (too few effective draws)."),
			paste0("L'\u00e9chantillonnage d'importance est trop pauvre pour : ", paste(lab$fr[ans$laws_dropped], collapse = " ; "),
				" (trop peu de tirages effectifs)."))
		.bay_tell(story, "Change of strategy: these laws are left out; the comparison is incomplete and the answer may change (a larger draws may restore them).",
			"Changement de strat\u00e9gie : ces lois sont \u00e9cart\u00e9es ; la comparaison est incompl\u00e8te et la r\u00e9ponse peut changer (un draws plus grand peut les r\u00e9tablir).", "decide")
	}
	pl <- sort(ans$p_law, decreasing = TRUE)
	.bay_tell(story, paste0("Posterior probabilities: ", paste0(lab$en[names(pl)], " ", vapply(pl, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
		paste0("Probabilit\u00e9s a posteriori : ", paste0(lab$fr[names(pl)], " ", vapply(pl, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))

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
		.bay_tell(story, paste0("Heavy tails prevail (", .bay_pct(p_t, "."), "): the answer will bear on the centres of the groups, not exactly on their means."),
			paste0("Les queues lourdes l'emportent (", .bay_pct(p_t, ","), ") : la r\u00e9ponse portera sur les centres des groupes, pas exactement sur leurs moyennes."), "decide")
	} else if (p_t > 0.5) {
		.bay_tell(story, paste0("Heavy tails are the more probable (", .bay_pct(p_t, "."), "), not enough to speak of centres only (the method asks for ",
				.bay_pct(.bay_heavy_threshold(), "."), " of probability): the answer mixes means and centres (centre: the typical value of a group, close to its median)."),
			paste0("Les queues lourdes sont les plus probables (", .bay_pct(p_t, ","), "), pas assez pour ne parler que des centres (la m\u00e9thode demande ",
				.bay_pct(.bay_heavy_threshold(), ","), " de probabilit\u00e9) : la r\u00e9ponse m\u00eale moyennes et centres (centre : la valeur typique d'un groupe, proche de sa m\u00e9diane)."), "decide")
	} else if (length(normal_kept)) {
		.bay_tell(story, paste0("The normal shapes are the more probable (", .bay_pct(1 - p_t, "."), "): the answer bears mainly on the means."),
			paste0("Les formes normales sont les plus probables (", .bay_pct(1 - p_t, ","), ") : la r\u00e9ponse porte surtout sur les moyennes."), "decide")
	}
	if (length(pl) > 1L && pl[1L] > 0.95) {
		.bay_tell(story, "One shape dominates: the answer is almost its own (the others are kept with their small weight).",
			"Une forme domine : la r\u00e9ponse est presque la sienne (les autres restent avec leur faible poids).", "decide")
	} else if (length(pl) > 1L) {
		.bay_tell(story, "No shape is chosen: the answer of each one is weighted by its probability.",
			"Aucune forme n'est choisie : la r\u00e9ponse de chacune est pond\u00e9r\u00e9e par sa probabilit\u00e9.", "decide")
	}

	# Step 5: the Bayes factor, its precision, what may change the strategy.
	.bay_tell(story, "Bayes factor and numerical precision.", "Facteur de Bayes et pr\u00e9cision num\u00e9rique.", "step",
		fun = "m.bay() {KefiR}")
	lb_txt <- list(en = .bay_bf_txt(ans$log_bf, "."), fr = .bay_bf_txt(ans$log_bf, ","))
	.bay_tell(story, paste0("Bayes factor, difference against none: ", lb_txt$en, if (abs(ans$log_bf) < log(3)) " (close to 1: the data hardly tell the two apart)." else if (ans$log_bf >= 0) " (the data favour a difference)." else " (the data favour no difference)."),
		paste0("Facteur de Bayes, diff\u00e9rence contre absence : ", lb_txt$fr, if (abs(ans$log_bf) < log(3)) " (proche de 1 : les donn\u00e9es d\u00e9partagent \u00e0 peine les deux)." else if (ans$log_bf >= 0) " (les donn\u00e9es favorisent une diff\u00e9rence)." else " (les donn\u00e9es favorisent l'absence de diff\u00e9rence)."))
	low <- ans$models$ess[!is.na(ans$models$ess)] < 200
	if (any(low) || isTRUE(ans$error_pct > 5)) {
		note(paste0("Numerical precision is low: error on the Bayes factor about ",
				format(signif(ans$error_pct, 2), trim = TRUE), " %."),
			paste0("Pr\u00e9cision num\u00e9rique faible : erreur sur le facteur de Bayes d'environ ",
				gsub(".", ",", format(signif(ans$error_pct, 2), trim = TRUE), fixed = TRUE), " %."))
	} else {
		ep <- .bay_num_txt(ans$error_pct, 2)
		.bay_tell(story, paste0("Imprecision of the computation by simulation: about ", ep$en, " % on the Bayes factor, without effect on the conclusion (this is not a risk of error)."),
			paste0("Impr\u00e9cision du calcul par simulation : environ ", ep$fr, " % sur le facteur de Bayes, sans effet sur la conclusion (ce n'est pas un risque d'erreur)."))
	}
	if (contradiction) {
		.bay_tell(story, paste0("According to the normal law (which compares means), the probability of a difference is ", .bay_pct(p_normal, "."),
				"; averaged over the shapes, dominated here by heavy tails (which compare centres), it is ", .bay_pct(ans$p_h1, "."),
				". The two readings fall on opposite sides: the centres and the means do not tell the same story."),
			paste0("Selon la loi normale (qui compare les moyennes), la probabilit\u00e9 d'une diff\u00e9rence est de ", .bay_pct(p_normal, ","),
				" ; en moyenne sur les formes, domin\u00e9es ici par les queues lourdes (qui comparent les centres), elle est de ", .bay_pct(ans$p_h1, ","),
				". Les deux lectures tombent de part et d'autre : centres et moyennes ne racontent pas la m\u00eame chose."))
		.bay_tell(story, paste0("Caution added: ", .bay_pct(ans$p_h1, "."), " is not taken at face value; the support is capped at weak and the verdict says that the answer depends on what is compared (means or centres)."),
			paste0("Prudence ajout\u00e9e : on ne retient pas ", .bay_pct(ans$p_h1, ","), " tel quel ; le soutien est plafonn\u00e9 \u00e0 faible et le verdict dit que la r\u00e9ponse d\u00e9pend de ce que l'on compare (moyennes ou centres)."), "decide")
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
	.bay_story_sensitivity(story, sens, hyp, het = any(laws %in% c("normal_het", "t_het")))

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
	.bay_tell(story, "Conclusion: the hypotheses ranked, what each one rests on, its support.",
		"Conclusion : les hypoth\u00e8ses class\u00e9es, ce sur quoi chacune repose, son soutien.", "step")
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(ans$p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

#' @export
print.mbay <- function(x, detail = isTRUE(x$detail), story = FALSE, ...) {
	told <- isTRUE(story) && length(x$story$en) > 0L
	if (told) .bay_story_replay(x$story)
	if (identical(x$plan, "factorial")) {
		.bay_print_bilan_fact(x, detail = detail, narrated = told)
	} else {
		.bay_print_bilan(x, detail = detail, narrated = told)
	}
	invisible(x)
}

# Several factors: normal law with a common variance, likelihoods from BayesFactor
# (sys_bay_factorial.R). prior_h1 and laws do not apply and are said so when given.
.bay_m_factorial <- function(input, call_expr, get_notes, note, prior_h1, rscale, laws, sensitivity,
                             draws, seed, detail, verbose, return, no_prior, no_laws, story = NULL) {
	y <- input$response
	factors <- input$factors
	terms <- gsub("`", "", input$terms, fixed = TRUE)
	sp <- .bay_fact_space(terms)
	n_mod <- length(sp$models)

	# Step 2 (title said by m.bay): sizes per cell, ties, values far from their cell.
	cell_all <- table(interaction(factors, drop = FALSE))
	chk <- .bay_story_checks(story, y, interaction(factors, drop = TRUE), "cell", "case")
	if (any(cell_all == 0L)) {
		.bay_tell(story, paste0("Empty cells: ", sum(cell_all == 0L), " of ", length(cell_all), "."),
			paste0("Cases vides : ", sum(cell_all == 0L), " sur ", length(cell_all), "."))
	}
	if (all(cell_all <= 1L)) {
		.bay_tell(story, "One observation per cell: no spread within a cell, the residual of the additive model will stand for it.",
			"Une observation par case : pas de dispersion au sein d'une case, le r\u00e9sidu du mod\u00e8le additif en tiendra lieu.")
	}
	.bay_tell(story, "No value is removed.", "Aucune valeur n'est retir\u00e9e.", "decide")

	# Step 3: one shape of the noise only.
	.bay_tell(story, "Shape of the noise.", "Forme du bruit.", "step")
	.bay_tell(story, "With several factors, only the normal law with a common variance is computed: heavy tails and unequal variances are not validated there yet, so they are not compared.",
		"\u00c0 plusieurs facteurs, seule la loi normale \u00e0 variance commune est calcul\u00e9e : queues lourdes et variances in\u00e9gales n'y sont pas encore valid\u00e9es, elles ne sont donc pas compar\u00e9es.")
	if (!no_laws && !identical(laws, "normal_eq")) {
		note("With several factors only the normal law with a common variance is computed; laws is ignored.",
			"\u00c0 plusieurs facteurs, seule la loi normale \u00e0 variance commune est calcul\u00e9e ; laws est ignor\u00e9.")
	}
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(story, paste0(far, " value(s) lie beyond 3 interquartile ranges of their cell (step 2): the answer assumes a normal noise, read it with caution."),
			paste0(far, " valeur(s) au-del\u00e0 de 3 \u00e9carts interquartiles de leur case (\u00e9tape 2) : la r\u00e9ponse suppose un bruit normal, la lire avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(story, "No value lies far enough from its cell to put the normal law in doubt at first sight; it is not tested against the others.",
			"Aucune valeur n'est assez \u00e9loign\u00e9e de sa case pour mettre d'embl\u00e9e la loi normale en doute ; elle n'est pas mise \u00e0 l'\u00e9preuve des autres.", "decide")
	}

	# Step 4: the priors, fixed by the method (DECISIONS, D2).
	.bay_tell(story, "Starting point (priors), set by the method before seeing the data.", "Point de d\u00e9part (a priori), fix\u00e9 par la m\u00e9thode avant de voir les donn\u00e9es.", "step")
	lp <- .bay_fact_prior(sp, "hop")
	incl <- vapply(seq_along(terms), function(i) sum(exp(lp[vapply(sp$models, function(m) i %in% m, logical(1))])), numeric(1))
	r_txt <- .bay_num_txt(rscale)
	.bay_tell(story, paste0("Size of an effect: r = ", r_txt$en, " per term (Zellner-Siow). Chance of each term, taken alone (these chances do not add up): ",
			paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = "."), collapse = ", "),
			"; no difference at all ", .bay_pct(exp(lp[1L]), "."), "."),
		paste0("Taille d'un effet : r = ", r_txt$fr, " par terme (Zellner-Siow). Chance de chaque terme, pris seul (ces chances ne s'additionnent pas) : ",
			paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = ","), collapse = ", "),
			" ; aucune diff\u00e9rence du tout ", .bay_pct(exp(lp[1L]), ","), "."))
	if (any(sp$order > 1L)) {
		.bay_tell(story, "Why: one proportion ~ Beta(1,1) per order of interaction (HOP, Taylor-Rodriguez et al. 2016): the more terms, the less likely each one in advance, which guards against false interactions; an interaction enters only with its factors.",
			"Pourquoi : une proportion ~ Beta(1,1) par ordre d'interaction (HOP, Taylor-Rodriguez et al. 2016) : plus il y a de termes, moins chacun est probable d'avance, ce qui prot\u00e8ge des fausses interactions ; une interaction n'entre qu'avec ses facteurs.", "decide")
	} else {
		.bay_tell(story, "Why: one proportion ~ Beta(1,1) for the factors (HOP, Taylor-Rodriguez et al. 2016): the more factors, the less likely each model in advance.",
			"Pourquoi : une proportion ~ Beta(1,1) pour les facteurs (HOP, Taylor-Rodriguez et al. 2016) : plus il y a de facteurs, moins chaque mod\u00e8le est probable d'avance.", "decide")
		note("No interaction is tested (additive formula or x, g form); write y ~ A * B to test them.",
			"Aucune interaction n'est test\u00e9e (formule additive ou forme x, g) ; \u00e9crire y ~ A * B pour les tester.")
	}
	if (!no_prior) {
		note("prior_h1 has no role with several factors: the prior on the models is fixed by the method.",
			"prior_h1 est sans effet \u00e0 plusieurs facteurs : la loi sur les mod\u00e8les est fix\u00e9e par la m\u00e9thode.")
	}

	# Step 5: Bayes factors of the models, precision, sensitivity.
	.bay_tell(story, paste0("Bayes factors of the ", n_mod, " models."), paste0("Facteurs de Bayes des ", n_mod, " mod\u00e8les."), "step",
		fun = "lmBF() {BayesFactor}")
	if (n_mod > 50L) {
		.bay_tell(story, paste0(n_mod, " models: this can take a few minutes."),
			paste0(n_mod, " mod\u00e8les : cela peut prendre quelques minutes."), "info")
	}
	fa <- .bay_factorial(input, rscale, sensitivity, as.integer(draws), as.integer(seed), note)
	if (isTRUE(fa$spread <= 0.02)) {
		sp_txt <- .bay_num_txt(100 * fa$spread, 2)
		.bay_tell(story, paste0("Precision: computed again with another seed, the probabilities move by at most ", sp_txt$en, " percentage point(s), negligible."),
			paste0("Pr\u00e9cision : recalcul\u00e9es avec une autre graine, les probabilit\u00e9s bougent d'au plus ", sp_txt$fr, " point(s) de pourcentage, n\u00e9gligeable."))
	}
	sens <- fa$sensitivity
	if (is.null(sens)) {
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
		.bay_tell(story, "The answers are not tested against other starting points.",
			"Les r\u00e9ponses ne sont pas \u00e9prouv\u00e9es contre d'autres points de d\u00e9part.", "decide")
	} else {
		r1 <- .bay_num_txt(min(sens$r)); r2 <- .bay_num_txt(max(sens$r))
		.bay_tell(story, paste0("Sensitivity: r from ", r1$en, " to ", r2$en, ", crossed with three priors on the models (HOP, 1/2 per term, uniform)."),
			paste0("Sensibilit\u00e9 : r de ", r1$fr, " \u00e0 ", r2$fr, ", crois\u00e9 avec trois lois sur les mod\u00e8les (HOP, 1/2 par terme, uniforme)."),
			fun = "m.bay() {KefiR}")
		moved <- which(sens$verdict_terms == "changes")
		if (length(moved)) {
			rr <- sens$range_terms[moved, , drop = FALSE]
			.bay_tell(story, paste0("Term(s) that change side: ", paste0(terms[moved], " (", vapply(rr[, 1L], .bay_pct, character(1), dec = "."), " to ",
					vapply(rr[, 2L], .bay_pct, character(1), dec = "."), ")", collapse = ", "), "."),
				paste0("Terme(s) qui changent de c\u00f4t\u00e9 : ", paste0(terms[moved], " (", vapply(rr[, 1L], .bay_pct, character(1), dec = ","), " \u00e0 ",
					vapply(rr[, 2L], .bay_pct, character(1), dec = ","), ")", collapse = ", "), "."))
			many <- length(moved) > 1L
			.bay_tell(story, paste0("Caution added: ", paste(terms[moved], collapse = ", "), if (many) " will be called not settled, whatever their central probability." else " will be called not settled, whatever its central probability."),
				paste0("Prudence ajout\u00e9e : ", paste(terms[moved], collapse = ", "), if (many) " seront dits non tranch\u00e9s, quelle que soit leur probabilit\u00e9 centrale." else " sera dit non tranch\u00e9, quelle que soit sa probabilit\u00e9 centrale."), "decide")
		} else {
			.bay_tell(story, "No term changes side.", "Aucun terme ne change de c\u00f4t\u00e9.", "decide")
		}
		if (identical(sens$verdict, "changes")) {
			.bay_tell(story, "The global answer (at least one factor or none) can reverse: the conclusion will not assert it.",
				"La r\u00e9ponse globale (au moins un facteur ou aucun) peut s'inverser : la conclusion ne l'affirmera pas.", "decide")
		}
	}

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
		descriptive_stats = fa$means, cell_means = local({
			cl <- interaction(input$factors, sep = " x ", drop = TRUE)
			data.frame(cell = levels(cl), n = as.integer(table(cl)), mean = as.numeric(tapply(input$response, cl, mean)), stringsAsFactors = FALSE)
		}), sd_within = fa$sd_within, sd_is_within = fa$sd_is_within, balanced = fa$balanced,
		prob_effect = p_any,
		method = "Several factors, normal law with a common variance (BayesFactor)", plan = input$plan,
		label = input$label, n = length(input$response), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	.bay_tell(story, "Conclusion: the hypotheses ranked, what each one rests on, its support.",
		"Conclusion : les hypoth\u00e8ses class\u00e9es, ce sur quoi chacune repose, son soutien.", "step")
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_fact(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_any)
	if (isTRUE(verbose)) invisible(out) else out
}
