# sys_bay_manova.R
# Several responses, one factor between groups (m.bay, DECISIONS, D22). The responses are
# standardised and stacked: y ~ g * response + subject, the subject random (lmBF, whichRandom;
# Rouder et al. 2012). "The groups differ on at least one response" is the union of the
# models with g (g alone: the same gap on every response; g + g:response: gaps that differ
# from one response to another), against the model with the responses and subjects only.
# The subject effect carries the correlation between the responses of one individual, taken
# equal for every pair (compound symmetry): that is the approximation of this route, measured
# in agent_68 with very unequal correlations. Chosen over the fractional Bayes factor of
# BFpack and over the multivariate BIC, both nearly powerless at 3 responses (agent_68).
# Nothing is exported.

# Long format of the standardised responses.
.bay_manova_long <- function(Y, g) {
	Z <- scale(Y)
	d <- ncol(Y)
	n <- nrow(Y)
	data.frame(y = as.vector(Z), g = rep(g, d), vd = factor(rep(colnames(Y), each = n), levels = colnames(Y)),
		id = factor(rep(seq_len(n), d)))
}

# Log Bayes factors of "g" and "g + g:response" against "response + subject", for each r.
.bay_manova_lbf <- function(Y, g, rs, iterations, seeds) {
	long <- .bay_manova_long(Y, g)
	one <- function(fo, r) BayesFactor::extractBF(BayesFactor::lmBF(fo, data = long, whichRandom = "id",
		rscaleFixed = r, iterations = iterations, progress = FALSE), logbf = TRUE)$bf
	vapply(seq_along(rs), function(i) {
		.bay_with_seed(seeds[i], {
			b0 <- one(y ~ vd + id, rs[i])
			c(g = one(y ~ g + vd + id, rs[i]) - b0, gx = one(y ~ g + vd + g:vd + id, rs[i]) - b0)
		})
	}, numeric(2))
}

# P(at least one difference): prior_h1 split evenly over the two forms of a difference.
.bay_manova_post <- function(lb, prior_h1) {
	l1 <- .bay_log_sum_exp(lb) - log(2)
	p <- 1 / (1 + exp(-(l1 + log(prior_h1 / (1 - prior_h1)))))
	w <- exp(lb - max(lb))
	list(log_bf = l1, p_h1 = p, p_form = w / sum(w))
}

.bay_m_manova <- function(input, call_expr, get_notes, note, prior_h1, rscale, sensitivity, draws, seed, detail, verbose,
                          return, story, posthoc) {
	if (length(input$covariates) || !is.null(input$id)) {
		.bay_tell(story, paste0("Plan recognised: several responses (MANOVA) with covariates or subjects; N = ", NROW(input$response), "."),
			paste0("Plan reconnu : plusieurs r\u00e9ponses (MANOVA) avec des covariables ou des sujets ; N = ", NROW(input$response), "."), "step")
		.exit("Several responses: factors between groups only for now (no covariate or subject).",
			"Plusieurs r\u00e9ponses : facteurs entre groupes seulement pour l'instant (ni covariable, ni sujet).")
	}
	# Several factors (D23).
	if (length(input$factors) > 1L) {
		return(.bay_m_manova_fact(input, call_expr, get_notes, note, rscale, sensitivity, draws, seed, detail, verbose, return, story))
	}
	Y <- as.matrix(input$response)
	if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
	colnames(Y) <- make.unique(colnames(Y))
	g <- droplevels(as.factor(input$factors[[1L]]))
	fac <- names(input$factors)[1L]
	if (identical(fac, ".g")) fac <- "g"
	k <- nlevels(g)
	d <- ncol(Y)
	N <- nrow(Y)
	n <- as.integer(table(g))
	if (any(apply(Y, 2L, stats::sd) == 0)) {
		.exit("A response is constant: it cannot be standardised.", "Une r\u00e9ponse est constante : elle ne peut pas \u00eatre standardis\u00e9e.")
	}
	if (N - k < d) {
		.exit("Fewer individuals than responses within groups: the covariance of the responses cannot be estimated.",
			"Moins d'individus que de r\u00e9ponses au sein des groupes : la covariance des r\u00e9ponses ne peut pas \u00eatre estim\u00e9e.")
	}
	resp_txt <- paste(colnames(Y), collapse = ", ")

	# Step 1: the plan.
	.bay_tell(story, paste0("Plan recognised: ", d, " responses (", resp_txt, ") compared together between the ", k, " groups of ", fac,
			" (", paste(levels(g), collapse = ", "), "); N = ", N, " individuals."),
		paste0("Plan reconnu : ", d, " r\u00e9ponses (", resp_txt, ") compar\u00e9es ensemble entre les ", k, " groupes de ", fac,
			" (", paste(levels(g), collapse = ", "), ") ; N = ", N, " individus."), "step")
	.bay_tell(story, paste0("One question: do the groups differ on at least one of the responses?"),
		paste0("Une question : les groupes diff\u00e8rent-ils sur au moins une des r\u00e9ponses ?"), "decide")

	# Step 2: checks.
	.bay_tell(story, "Checks on the data", "Contr\u00f4les des donn\u00e9es", "step", fun = "quantile() {stats}")
	.bay_tell(story, paste0("Sizes: ", min(n), if (max(n) > min(n)) paste0(" to ", max(n)), " per group."),
		paste0("Effectifs : ", min(n), if (max(n) > min(n)) paste0(" \u00e0 ", max(n)), " par groupe."))
	far <- vapply(seq_len(d), function(j) {
		r <- Y[, j] - stats::ave(Y[, j], g)
		q <- stats::quantile(r, c(0.25, 0.75))
		sum(r < q[1L] - 3 * diff(q) | r > q[2L] + 3 * diff(q))
	}, numeric(1))
	.bay_step_result(story, if (sum(far)) paste0(sum(far), " far values kept") else "nothing to remove",
		if (sum(far)) paste0(sum(far), " valeurs \u00e9loign\u00e9es gard\u00e9es") else "rien \u00e0 retirer")
	.bay_tell(story, paste0("Values beyond 3 interquartile ranges from their group: ", paste0(colnames(Y), " ", far, collapse = ", "), "."),
		paste0("Valeurs au-del\u00e0 de 3 \u00e9carts interquartiles de leur groupe : ", paste0(colnames(Y), " ", far, collapse = ", "), "."))
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	# Step 3: shape of the noise.
	.bay_tell(story, "Shape of the noise", "Forme du bruit", "step")
	.bay_step_result(story, "normal on several responses, correlations taken equal", "normale \u00e0 plusieurs r\u00e9ponses, corr\u00e9lations prises \u00e9gales")
	E <- Y - apply(Y, 2L, function(v) stats::ave(v, g))
	R <- stats::cor(E)
	rr <- R[upper.tri(R)]
	f2 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
	.bay_tell(story, paste0("Each response is standardised; correlations between responses within groups: ",
			if (length(rr) == 1L) f2(rr, ".") else paste0("from ", f2(min(rr), "."), " to ", f2(max(rr), ".")), "."),
		paste0("Chaque r\u00e9ponse est standardis\u00e9e ; corr\u00e9lations entre r\u00e9ponses au sein des groupes : ",
			if (length(rr) == 1L) f2(rr, ",") else paste0("de ", f2(min(rr), ","), " \u00e0 ", f2(max(rr), ",")), "."))
	uneq <- length(rr) > 1L && diff(range(rr)) > 0.5
	if (uneq) {
		# agent_68: with very unequal correlations the stacked model stays calibrated under H0 but
		# can miss a difference that a contrast between correlated responses would show.
		.bay_tell(story, "The model takes these correlations equal; here they are far from it: it may then miss a difference, without inventing one; a \"no\" is to be read with caution.",
			"Le mod\u00e8le prend ces corr\u00e9lations \u00e9gales ; ici elles en sont loin : il peut alors manquer une diff\u00e9rence, sans en inventer ; un \u00ab non \u00bb est \u00e0 lire avec r\u00e9serve.", "decide")
	} else {
		.bay_tell(story, "The model takes these correlations equal, which they roughly are.",
			"Le mod\u00e8le prend ces corr\u00e9lations \u00e9gales, ce qu'elles sont \u00e0 peu pr\u00e8s.", "decide")
	}
	.bay_tell(story, "equal spreads in the groups are assumed, as in a classical MANOVA.",
		"des dispersions \u00e9gales dans les groupes sont suppos\u00e9es, comme dans une MANOVA classique.", "note")
	.bay_tell(story, "Responses standardised and stacked: y ~ g * response + individual (random), lmBF {BayesFactor}; the individual effect carries the correlation, the same for every pair (compound symmetry; DECISIONS, D22).",
		"R\u00e9ponses standardis\u00e9es et empil\u00e9es : y ~ g * r\u00e9ponse + individu (al\u00e9atoire), lmBF {BayesFactor} ; l'effet individu porte la corr\u00e9lation, la m\u00eame pour chaque paire (sym\u00e9trie compos\u00e9e ; DECISIONS, D22).", "detail")

	# Step 4: priors.
	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference"), paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence"))
	r_txt <- .bay_num_txt(rscale)
	.bay_tell(story, paste0("Two forms of a difference, at even odds: the same gap on every response, or gaps that change from one response to another; expected size r = ", r_txt$en, "."),
		paste0("Deux formes d'une diff\u00e9rence, \u00e0 poids \u00e9gaux : le m\u00eame \u00e9cart sur chaque r\u00e9ponse, ou des \u00e9carts qui changent d'une r\u00e9ponse \u00e0 l'autre ; taille attendue r = ", r_txt$fr, "."))
	.bay_tell(story, "these values are not tuned on the data; step 5 varies them.",
		"ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")

	# Step 5: solidity.
	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = "lmBF() {BayesFactor}")
	rs <- if (isTRUE(sensitivity)) c(rscale, rscale / 2, rscale * 2) else rscale
	iterations <- max(10000L, as.integer(draws))
	lbm <- .bay_manova_lbf(Y, g, rs, iterations, seed + c(0L, 1L, 2L)[seq_along(rs)])
	po <- .bay_manova_post(lbm[, 1L], prior_h1)
	lb <- po$log_bf
	lb2 <- .bay_manova_lbf(Y, g, rscale, iterations, seed + 7919L)
	spread <- abs(.bay_manova_post(lb2[, 1L], prior_h1)$p_h1 - po$p_h1)
	inv_en <- if (lb < 0) paste0(", i.e. ", .bay_bf_txt(-lb, "."), " in favour of no difference") else ""
	inv_fr <- if (lb < 0) paste0(", soit ", .bay_bf_txt(-lb, ","), " en faveur de l'absence de diff\u00e9rence") else ""
	.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb, "."), inv_en, "."), paste0("Facteur de Bayes : ", .bay_bf_txt(lb, ","), inv_fr, "."))
	sp_txt <- .bay_prec_txt(100 * spread, "point")
	.bay_tell(story, paste0("Precision of the computation: ", sp_txt$en, " on the probability."),
		paste0("Pr\u00e9cision du calcul : ", sp_txt$fr, " sur la probabilit\u00e9."))
	sens <- NULL
	if (isTRUE(sensitivity)) {
		phs <- c(1 / 3, 1 / 2, 2 / 3)
		grid <- as.vector(outer(seq_along(rs), phs, function(i, p) {
			vapply(seq_along(i), function(j) .bay_manova_post(lbm[, i[j]], p[j])$p_h1, numeric(1))
		}))
		lead_h1 <- po$p_h1 >= 0.5
		rng <- range(if (lead_h1) grid else 1 - grid)
		verdict <- .bay_sensitivity_verdict(max(po$p_h1, 1 - po$p_h1), rng)
		sens <- list(r = sort(rs), prior_h1 = phs, p = grid, range_lead = rng, verdict = verdict, verdict_terms = verdict,
			range_terms = matrix(range(grid), 1L))
		stb <- .bay_stab_txt(verdict)
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", ", stb$fr))
		.bay_tell(story, paste0("Sensitivity: the leading answer ", .bay_range_txt(rng)$en, "."),
			paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate ", .bay_range_txt(rng)$fr, "."))
		.bay_tell(story, if (identical(verdict, "changes")) "Caution added: the answer can reverse; it will be called not settled." else "The order does not change.",
			if (identical(verdict, "changes")) "Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e." else "L'ordre ne change pas.", "decide")
		.bay_tell(story, "r goes from r/2 to 2r, the starting chances from 1 against 2 to 2 against 1.",
			"r va de r/2 \u00e0 2r, les chances de d\u00e9part de 1 contre 2 \u00e0 2 contre 1.", "note")
	} else {
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", sensitivity not computed"), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", sensibilit\u00e9 non calcul\u00e9e"))
	}

	# Step 6: which responses carry the difference (one factor engine on each, for reading).
	per <- NULL
	if (isTRUE(posthoc)) {
		.bay_tell(story, "Which responses differ", "Quelles r\u00e9ponses diff\u00e8rent", "step", fun = "lmBF() {BayesFactor}")
		per <- data.frame(response = colnames(Y), log_bf = vapply(seq_len(d), function(j) .bay_bf_oneway_zs(Y[, j], g, rscale), numeric(1)),
			stringsAsFactors = FALSE)
		per$p_h1 <- 1 / (1 + exp(-per$log_bf))
		top <- per$response[which.max(per$p_h1)]
		.bay_step_result(story, paste0(top, " ahead (", .bay_pct(max(per$p_h1), "."), ")"), paste0(top, " en t\u00eate (", .bay_pct(max(per$p_h1), ","), ")"))
		for (j in seq_len(d)) {
			m <- tapply(Y[, j], g, mean)
			.bay_tell(story, paste0(per$response[j], ": difference between groups at ", .bay_pct(per$p_h1[j], "."), "; means ",
					paste0(names(m), " = ", .bay_fx(m, .bay_dg(m), "."), collapse = ", "), "."),
				paste0(per$response[j], " : diff\u00e9rence entre groupes \u00e0 ", .bay_pct(per$p_h1[j], ","), " ; moyennes ",
					paste0(names(m), " = ", .bay_fx(m, .bay_dg(m), ","), collapse = ", "), "."))
		}
		.bay_tell(story, paste0("If a difference exists, gaps that change from one response to another: ", .bay_pct(po$p_form[2L], "."), "."),
			paste0("Si une diff\u00e9rence existe, \u00e9carts qui changent d'une r\u00e9ponse \u00e0 l'autre : ", .bay_pct(po$p_form[2L], ","), "."))
		.bay_tell(story, "each response read alone (one factor, even odds), to see where the difference lies; the global answer is that of step 5.",
			"chaque r\u00e9ponse lue seule (un facteur, chances \u00e9gales), pour voir o\u00f9 se trouve la diff\u00e9rence ; la r\u00e9ponse globale est celle de l'\u00e9tape 5.", "note")
	} else {
		.bay_story_no_posthoc(story, "off")
	}

	resp <- paste(colnames(Y), collapse = ", ")
	hyp <- data.frame(
		label_en = c(paste0("At least one of the responses (", resp, ") differs according to ", fac), paste0("No difference in ", resp, " according to ", fac)),
		label_fr = c(paste0("Au moins une des r\u00e9ponses (", resp, ") diff\u00e8re selon ", fac), paste0("Aucune diff\u00e9rence de ", resp, " selon ", fac)),
		posterior = c(po$p_h1, 1 - po$p_h1), prior = c(prior_h1, 1 - prior_h1), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(hypotheses = hyp,
		effects = data.frame(term = fac, bf10 = exp(lb), log_bf = lb, p_h1 = po$p_h1, p_h0 = 1 - po$p_h1, prior_h1 = prior_h1,
			rscale = rscale, error_pct = 100 * spread, analogue = "MANOVA (Pillai)", stringsAsFactors = FALSE),
		# Same elements as the other plans (audit of D22): laws, models, questions.
		laws = NULL,
		models = data.frame(model = c("(null): response + individual", "g + response + individual", "g + response + g:response + individual"),
			log_bf = c(0, lbm[, 1L]), log_prior = log(c(1 - prior_h1, prior_h1 / 2, prior_h1 / 2)), stringsAsFactors = FALSE),
		questions = list(en = paste0("Do the groups of ", fac, " differ on at least one of ", resp_txt, "?"),
			fr = paste0("Les groupes de ", fac, " diff\u00e8rent-ils sur au moins une de ", resp_txt, " ?")),
		forms = data.frame(form = c("same gap on every response", "gaps that change with the response"), log_bf = lbm[, 1L],
			posterior_given_h1 = po$p_form, stringsAsFactors = FALSE),
		responses = per, correlations = R, unequal_correlations = uneq,
		sensitivity = sens, posthoc = NULL,
		descriptive_stats = data.frame(group = rep(levels(g), d), response = rep(colnames(Y), each = k), n = rep(n, d),
			mean = as.vector(apply(Y, 2L, function(v) tapply(v, g, mean))), stringsAsFactors = FALSE),
		priors = list(prior_h1 = prior_h1, rscale = rscale, forms = "1/2, 1/2", effects = "Zellner-Siow on the standardised responses"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = po$p_h1,
		method = "Several responses: standardised, stacked, lmBF with a random individual effect (BayesFactor)", plan = "manova",
		label = input$label, n = N, notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_manova(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(po$p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_print_bilan_manova <- function(obj, detail = FALSE, narrated = FALSE) {
	body <- function(emit, li, narrated) {
		per <- obj$responses
		if (!is.null(per)) {
			emit("fact", paste0("Each response read alone: ", paste0(per$response, " ", vapply(per$p_h1, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
				paste0("Chaque r\u00e9ponse lue seule : ", paste0(per$response, " ", vapply(per$p_h1, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))
		}
		if (isTRUE(obj$unequal_correlations)) {
			emit("decide", "Correlations between responses far from equal: the model may miss a difference; a \"no\" is to be read with caution.",
				"Corr\u00e9lations entre r\u00e9ponses loin d'\u00eatre \u00e9gales : le mod\u00e8le peut manquer une diff\u00e9rence ; un \u00ab non \u00bb est \u00e0 lire avec r\u00e9serve.")
		}
		e <- obj$effects
		if (nrow(e) > 1L) {
			# Several factors (D23): one line per term, with the words of the report sentence.
			ta <- .bay_term_answers(obj)
			for (i in seq_len(nrow(e))) {
				emit("fact", paste0(e$term[i], ": ", ta$long_en[i], " (", .bay_pct(e$p_h1[i], "."), "; prior ", .bay_pct(e$prior_h1[i], "."), "; Bayes factor ", .bay_bf_txt(e$log_bf[i], "."), ")."),
					paste0(e$term[i], " : ", ta$long_fr[i], " (", .bay_pct(e$p_h1[i], ","), " ; a priori ", .bay_pct(e$prior_h1[i], ","), " ; facteur de Bayes ", .bay_bf_txt(e$log_bf[i], ","), ")."))
			}
		}
		emit("block", "", "")
		ds <- obj$descriptive_stats
		grp <- if (!is.null(ds$group)) ds$group else ds$cell
		for (v in unique(ds$response)) {
			s <- ds$response == v
			dg <- .bay_dg(ds$mean[s])
			emit("fact", paste0("Observed means, ", v, ": ", paste0(grp[s], " = ", .bay_fx(ds$mean[s], dg, "."), collapse = "; "), "."),
				paste0("Moyennes observ\u00e9es, ", v, " : ", paste0(grp[s], " = ", .bay_fx(ds$mean[s], dg, ","), collapse = " ; "), "."))
		}
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)
	if (is.null(obj$forms)) {
		say("Each term of the factors in two forms at even odds: t (same gap on every response) or t + t:response; HOP prior on the terms; responses standardised, stacked, individual random effect (lmBF; DECISIONS, D23).",
			"Chaque terme des facteurs sous deux formes \u00e0 poids \u00e9gaux : t (m\u00eame \u00e9cart sur chaque r\u00e9ponse) ou t + t:r\u00e9ponse ; a priori HOP sur les termes ; r\u00e9ponses standardis\u00e9es, empil\u00e9es, effet individu al\u00e9atoire (lmBF ; DECISIONS, D23).")
		return(invisible(obj))
	}
	say(paste0("Forms of a difference, given a difference: same gap on every response ", .bay_pct(obj$forms$posterior_given_h1[1L], "."),
			"; gaps that change with the response ", .bay_pct(obj$forms$posterior_given_h1[2L], "."), "."),
		paste0("Formes d'une diff\u00e9rence, s'il y en a une : m\u00eame \u00e9cart sur chaque r\u00e9ponse ", .bay_pct(obj$forms$posterior_given_h1[1L], ","),
			" ; \u00e9carts qui changent avec la r\u00e9ponse ", .bay_pct(obj$forms$posterior_given_h1[2L], ","), "."))
	say("Responses standardised (mean 0, SD 1), stacked; individual random effect (compound symmetry); Zellner-Siow prior of scale r on g and g:response (Rouder et al. 2012).",
		"R\u00e9ponses standardis\u00e9es (moyenne 0, \u00e9cart-type 1), empil\u00e9es ; effet individu al\u00e9atoire (sym\u00e9trie compos\u00e9e) ; a priori de Zellner-Siow d'\u00e9chelle r sur g et g:r\u00e9ponse (Rouder et al. 2012).")
	invisible(obj)
}
