# sys_bay_contingency.R
# One factor and a binary or categorical response (DECISIONS, D17, way 1): the Bayesian
# counterpart of the chi-squared and G tests of m.test().
#   Whether the proportions differ: contingencyTableBF() [BayesFactor], independent
#     multinomials with the group margin fixed (Gunel and Dickey 1974), uniform Dirichlet
#     prior of concentration a on the proportions of each group. Exact, no draws.
#   Which groups differ: the partitions of D11; the Bayes factor of a grouping is that of the
#     same model on the merged table. Under this prior it is a ratio of Dirichlet-multinomial
#     marginal likelihoods, prod over blocks DM(n_b; a) / DM(n; a), identical to
#     contingencyTableBF (test-bay-contingency.R), computed here in closed form.
#   How much: draws from the Dirichlet posterior of each group; two categories: difference of
#     proportions and log odds ratio; more: total variation distance (half the sum of the
#     absolute differences of proportions).
#
# References.
#   Gunel, E., & Dickey, J. (1974). Bayes factors for independence in contingency tables.
#     Biometrika, 61(3), 545-557.
#   Jamil, T., Ly, A., Morey, R. D., Love, J., Marsman, M., & Wagenmakers, E.-J. (2017).
#     Default "Gunel and Dickey" Bayes factors for contingency tables. Behavior Research
#     Methods, 49, 638-652.

# Log Dirichlet-multinomial marginal likelihood, without the multinomial coefficient, which
# cancels in every Bayes factor below.
.bay_ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))

# Log Bayes factor "each row its own proportions" against "the same proportions", rows fixed.
# Prior of contingencyTableBF (checked to 1e-6, test-bay-contingency.R): each row Dirichlet(a)
# under H1; the common proportions Dirichlet(R (a - 1) + 1) under H0, so a >= 1. At a = 1 both
# are uniform, and the Bayes factor of a grouping of the rows is that of the merged table.
.bay_cont_lbf <- function(tab, conc = 1) {
	C <- ncol(tab); R <- nrow(tab)
	sum(apply(tab, 1L, .bay_ldm, a = rep(conc, C))) - .bay_ldm(colSums(tab), rep(R * (conc - 1) + 1, C))
}

.bay_cont_bf_pkg <- function(tab, conc = 1) {
	BayesFactor::extractBF(BayesFactor::contingencyTableBF(tab, sampleType = "indepMulti", fixedMargin = "rows",
		priorConcentration = conc), logbf = TRUE)$bf
}

.bay_cont_posthoc <- function(tab, conc, prior_h1, draws, seed, alphas = c(1, 1/2, 2)) {
	k <- nrow(tab)
	lev <- rownames(tab)
	pairs <- utils::combn(k, 2L)
	res <- data.frame(group1 = lev[pairs[1L, ]], group2 = lev[pairs[2L, ]], stringsAsFactors = FALSE)
	binary <- ncol(tab) == 2L
	prop <- tab[, ncol(tab)] / rowSums(tab)
	parts_tab <- NULL
	computed <- k <= .bay_posthoc_max_groups()
	if (computed) {
		parts <- .bay_partitions(k)
		lb <- vapply(parts, function(p) if (max(p) == 1L) 0 else .bay_cont_lbf(rowsum(tab, p, reorder = TRUE), conc), numeric(1))
		pm <- do.call(rbind, parts)
		post_of <- function(alpha) {
			lw <- lb + .bay_partition_prior(parts, prior_h1, alpha)
			exp(lw - .bay_log_sum_exp(lw))
		}
		sep_m <- pm[, pairs[1L, ], drop = FALSE] != pm[, pairs[2L, ], drop = FALSE]
		sep <- function(post) as.numeric(crossprod(sep_m, post))
		post <- post_of(alphas[1L])
		res$p_differ <- sep(post)
		ps <- matrix(vapply(alphas, function(a) sep(post_of(a)), numeric(ncol(pairs))), nrow = ncol(pairs))
		res$p_min <- apply(ps, 1L, min)
		res$p_max <- apply(ps, 1L, max)
		parts_tab <- data.frame(partition = vapply(parts, .bay_partition_label, character(1), lev = lev, means = as.numeric(prop)),
			code = vapply(parts, paste, character(1), collapse = ""), posterior = post,
			prior = exp(.bay_partition_prior(parts, prior_h1, alphas[1L])), log_bf = lb, stringsAsFactors = FALSE)
		parts_tab <- parts_tab[order(-parts_tab$posterior), , drop = FALSE]
		rownames(parts_tab) <- NULL
	} else {
		res$p_differ <- res$p_min <- res$p_max <- NA_real_
	}
	S <- max(4000L, as.integer(draws))
	P <- .bay_with_seed(seed, lapply(seq_len(k), function(i) {
		gm <- matrix(stats::rgamma(S * ncol(tab), shape = conc + tab[i, ]), S, ncol(tab), byrow = TRUE)
		gm / rowSums(gm)
	}))
	hdi_of <- function(v) {
		h <- bayestestR::hdi(as.numeric(v), ci = 0.95)
		c(stats::median(v), h$CI_low, h$CI_high)
	}
	st <- t(vapply(seq_len(ncol(pairs)), function(j) {
		a <- P[[pairs[1L, j]]]; b <- P[[pairs[2L, j]]]
		if (binary) {
			pa <- a[, 2L]; pb <- b[, 2L]
			c(hdi_of(pb - pa), hdi_of(stats::qlogis(pb) - stats::qlogis(pa)))
		} else {
			c(hdi_of(0.5 * rowSums(abs(b - a))), NA_real_, NA_real_, NA_real_)
		}
	}, numeric(6)))
	res$diff <- st[, 1L]; res$diff_low <- st[, 2L]; res$diff_high <- st[, 3L]
	res$delta <- st[, 4L]; res$delta_low <- st[, 5L]; res$delta_high <- st[, 6L]
	list(pairs = res, partitions = parts_tab, k = k, alphas = alphas, draws = S, computed = computed,
		method = if (computed) "exact" else "none", n_partitions = .bay_bell(k), sweeps = NA_integer_, mc_error = 0,
		max_exact = .bay_posthoc_max_groups(), delta_unit = if (binary) "logor" else "tvd",
		category = colnames(tab)[ncol(tab)], caution = FALSE)
}

# Whole analysis. g: the factor; yc: the response as a factor of categories.
.bay_m_contingency <- function(input, g, yc, call_expr, get_notes, note, prior_h1, sensitivity, draws, seed,
                               detail, verbose, return, story, posthoc, no_laws, laws) {
	tab <- table(g, yc)
	names(dimnames(tab)) <- NULL
	k <- nrow(tab); m <- ncol(tab)
	resp <- input$response_name
	fac <- names(input$factors)[1L]
	if (identical(fac, ".g")) fac <- "g"
	.bay_tell(story, paste0("Counts per group and category (", k, " groups x ", m, " categories): ",
			paste0(rownames(tab), " ", apply(tab, 1L, paste, collapse = "/"), collapse = "; "), " (", paste(colnames(tab), collapse = "/"), ")."),
		paste0("Effectifs par groupe et modalit\u00e9 (", k, " groupes x ", m, " modalit\u00e9s) : ",
			paste0(rownames(tab), " ", apply(tab, 1L, paste, collapse = "/"), collapse = " ; "), " (", paste(colnames(tab), collapse = "/"), ")."))
	if (any(tab == 0L)) {
		.bay_tell(story, paste0(sum(tab == 0L), " empty cell(s): no problem for this computation, which needs no minimum count."),
			paste0(sum(tab == 0L), " case(s) vide(s) : sans probl\u00e8me pour ce calcul, qui ne demande aucun effectif minimum."))
	}
	.bay_tell(story, "No approximation is made: the computation is exact, small counts do not make it wrong, only less decisive (unlike the chi-squared, which needs expected counts of 5).",
		"Aucune approximation : le calcul est exact, les petits effectifs ne le faussent pas, ils le rendent seulement moins tranch\u00e9 (contrairement au khi-deux, qui demande des effectifs attendus de 5).", "decide")
	.bay_tell(story, "Law of the data.", "Loi des donn\u00e9es.", "step")
	.bay_tell(story, paste0("Each group has its own proportions of the ", m, " categories (", if (m == 2L) "binomial" else "multinomial", " law); \"no difference\": the same proportions in every group. There is no shape of the noise to compare."),
		paste0("Chaque groupe a ses propres proportions des ", m, " modalit\u00e9s (loi ", if (m == 2L) "binomiale" else "multinomiale", ") ; \u00ab aucune diff\u00e9rence \u00bb : les m\u00eames proportions dans tous les groupes. Il n'y a pas de forme du bruit \u00e0 comparer."))
	if (!no_laws) {
		note("For a categorical response the laws of the noise do not apply; laws is ignored.",
			"Pour une r\u00e9ponse cat\u00e9gorielle, les lois du bruit ne s'appliquent pas ; laws est ignor\u00e9.")
	}

	.bay_tell(story, "Starting point (priors), set before seeing the data.", "Point de d\u00e9part (a priori), fix\u00e9 avant de voir les donn\u00e9es.", "step")
	ph <- .bay_num_txt(100 * prior_h1)
	.bay_tell(story, paste0(ph$en, " % chance that a difference exists", if (prior_h1 == 0.5) " (heads or tails: the method takes no side)" else "",
			"; in each group, every set of proportions equally likely in advance (uniform Dirichlet, Gunel and Dickey 1974)."),
		paste0(ph$fr, " % de chances qu'une diff\u00e9rence existe", if (prior_h1 == 0.5) " (pile ou face : la m\u00e9thode ne prend pas parti)" else "",
			" ; dans chaque groupe, toutes les proportions sont \u00e9galement possibles d'avance (Dirichlet uniforme, Gunel et Dickey 1974)."))
	.bay_tell(story, "These values are not fitted on the data; step 5 varies them to see whether the answer holds.",
		"Ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier pour voir si la r\u00e9ponse tient.", "decide")

	.bay_tell(story, "Bayes factor and sensitivity to the priors.", "Facteur de Bayes et sensibilit\u00e9 aux a priori.", "step",
		fun = "contingencyTableBF() {BayesFactor}")
	lb <- .bay_cont_bf_pkg(tab, 1)
	p_h1 <- .bay_cont_prob(lb, prior_h1)
	b1 <- .bay_num_txt(exp(lb))
	.bay_tell(story, paste0("Bayes factor, difference against none: ", .bay_bf_txt(lb, "."), "."),
		paste0("Facteur de Bayes, diff\u00e9rence contre absence : ", .bay_bf_txt(lb, ","), "."))
	.bay_tell(story, "Precision: closed-form computation, no draws; no simulation error.",
		"Pr\u00e9cision : calcul en forme close, sans tirage ; aucune erreur de simulation.")
	sens <- NULL
	if (isTRUE(sensitivity)) {
		# Concentration of the prior: 1 (uniform), 2 and 4 (proportions expected less extreme);
		# contingencyTableBF refuses a < 1, its prior under H0 being Dirichlet(R (a - 1) + 1).
		concs <- c(1, 2, 4); phs <- c(1 / 3, 1 / 2, 2 / 3)
		lbs <- vapply(concs, function(a) .bay_cont_bf_pkg(tab, a), numeric(1))
		grid <- outer(lbs, phs, function(l, p) .bay_cont_prob(l, p))
		lead_h1 <- p_h1 >= 0.5
		rng <- range(if (lead_h1) grid else 1 - grid)
		verdict <- .bay_sensitivity_verdict(max(p_h1, 1 - p_h1), rng)
		sens <- list(conc = concs, prior_h1 = phs, p = grid, range_lead = rng, verdict = verdict,
			verdict_terms = verdict, range_terms = matrix(range(grid), 1L))
		.bay_tell(story, paste0("Sensitivity: concentration of the prior from 1 to 4, crossed with starting chances of a difference from 1 against 2 to 2 against 1. The leading answer then goes from ",
				.bay_pct(rng[1L], "."), " to ", .bay_pct(rng[2L], "."), "."),
			paste0("Sensibilit\u00e9 : concentration de l'a priori de 1 \u00e0 4, crois\u00e9e avec des chances de d\u00e9part d'une diff\u00e9rence de 1 contre 2 \u00e0 2 contre 1. La r\u00e9ponse en t\u00eate va alors de ",
				.bay_pct(rng[1L], ","), " \u00e0 ", .bay_pct(rng[2L], ","), "."))
		if (identical(verdict, "changes")) {
			.bay_tell(story, "Caution added: the answer can reverse with the starting point; it will be called not settled.",
				"Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser selon le point de d\u00e9part ; elle sera dite non tranch\u00e9e.", "decide")
		} else {
			.bay_tell(story, "The answer does not change side.", "La r\u00e9ponse ne change pas de c\u00f4t\u00e9.", "decide")
		}
	} else {
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	ph_obj <- NULL
	if (isTRUE(posthoc)) {
		ph_obj <- .bay_cont_posthoc(tab, 1, prior_h1, draws, seed)
		.bay_story_posthoc(story, ph_obj)
	}

	hyp <- data.frame(
		label_en = c(paste0("The proportions of ", resp, " differ according to ", fac), paste0("The proportions of ", resp, " do not differ according to ", fac)),
		label_fr = c(paste0("Les proportions de ", resp, " diff\u00e8rent selon ", fac), paste0("Les proportions de ", resp, " ne diff\u00e8rent pas selon ", fac)),
		posterior = c(p_h1, 1 - p_h1), prior = c(prior_h1, 1 - prior_h1), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(hypotheses = hyp,
		effects = data.frame(term = fac, bf10 = exp(lb), log_bf = lb, p_h1 = p_h1, p_h0 = 1 - p_h1, prior_h1 = prior_h1,
			analogue = "chi-squared / G test", stringsAsFactors = FALSE),
		table = tab, proportions = prop.table(tab, 1L), sensitivity = sens, posthoc = ph_obj,
		priors = list(concentration = 1, prior_h1 = prior_h1, law = "multinomial"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = p_h1,
		method = "Contingency table, independent multinomials (BayesFactor, Gunel and Dickey 1974)", plan = "contingency",
		label = input$label, n = sum(tab), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	.bay_tell(story, "Conclusion: the hypotheses ranked, what each one rests on, its support.",
		"Conclusion : les hypoth\u00e8ses class\u00e9es, ce sur quoi chacune repose, son soutien.", "step")
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_cont(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_cont_prob <- function(lb, prior_h1) 1 / (1 + exp(-lb) * (1 - prior_h1) / prior_h1)

.bay_print_bilan_cont <- function(obj, detail = FALSE, narrated = FALSE) {
	lang_fr <- identical(.msg("en", "fr"), "fr")
	pre <- if (isTRUE(narrated)) "\t" else ""
	say <- function(en, fr) cat(pre, .msg(en, fr), "\n", sep = "")
	if (!isTRUE(narrated)) {
		say(paste0("Bilan of m.bay: ", obj$label, " (N = ", obj$n, ", contingency table)"),
			paste0("Bilan de m.bay : ", obj$label, " (N = ", obj$n, ", tableau de contingence)"))
	}
	v <- .bay_verdict(obj)
	say(v$en, v$fr)
	h <- obj$hypotheses
	e <- obj$effects
	lb_lead <- if (h$null[1L]) -e$log_bf else e$log_bf
	sup <- .bay_support_bf(h$posterior[1L], lb_lead)
	open <- identical(obj$sensitivity$verdict, "changes")
	say("Probability of each hypothesis, given the data and the starting point:",
		"Probabilit\u00e9 de chaque hypoth\u00e8se, compte tenu des donn\u00e9es et du point de d\u00e9part :")
	for (i in 1:2) {
		q <- if (i == 1L) paste0(" (", if (open) "not settled" else sup$en, ")") else ""
		qf <- if (i == 1L) paste0(" (", if (open) "non tranch\u00e9" else sup$fr, ")") else ""
		say(paste0("  ", i, ". ", h$label_en[i], ": ", .bay_pct(h$posterior[i], "."), q),
			paste0("  ", i, ". ", h$label_fr[i], " : ", .bay_pct(h$posterior[i], ","), qf))
	}
	bf_en <- .bay_bf_txt(abs(e$log_bf), "."); bf_fr <- .bay_bf_txt(abs(e$log_bf), ",")
	if (e$log_bf >= 0) {
		say(paste0("     Rests on the counts: they are ", bf_en, " times more probable if the proportions differ than if they are the same."),
			paste0("     Repose sur les effectifs : ils sont ", bf_fr, " fois plus probables si les proportions diff\u00e8rent que si elles sont les m\u00eames."))
	} else {
		say(paste0("     Rests on the counts: they are ", bf_en, " times more probable if the proportions are the same than if they differ."),
			paste0("     Repose sur les effectifs : ils sont ", bf_fr, " fois plus probables si les proportions sont les m\u00eames que si elles diff\u00e8rent."))
	}
	if (!is.null(obj$sensitivity)) {
		rg <- obj$sensitivity$range_lead
		if (open) {
			say(paste0("  Depending on the starting point, the first answer goes from ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), ": the order can reverse."),
				paste0("  Selon le point de d\u00e9part, la premi\u00e8re r\u00e9ponse va de ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " : l'ordre peut s'inverser."))
		} else {
			say(paste0("  Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order does not change."),
				paste0("  Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre ne change pas."))
		}
	}
	if (!is.null(obj$posthoc)) .bay_print_posthoc(obj$posthoc, say)
	pr <- obj$proportions
	say("Observed proportions per group (description of the sample):", "Proportions observ\u00e9es par groupe (description de l'\u00e9chantillon) :")
	for (i in seq_len(nrow(pr))) {
		say(paste0("  ", rownames(pr)[i], ": ", paste0(colnames(pr), " ", vapply(pr[i, ], .bay_pct, character(1), dec = "."), collapse = ", "), " (n = ", sum(obj$table[i, ]), ")"),
			paste0("  ", rownames(pr)[i], " : ", paste0(colnames(pr), " ", vapply(pr[i, ], .bay_pct, character(1), dec = ","), collapse = ", "), " (n = ", sum(obj$table[i, ]), ")"))
	}
	say("These are probabilities of hypotheses, not p-values: no threshold of 5 % applies. Whether a difference has a cause depends on how the study was run (randomisation, order, operator).",
		"Ce sont des probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values : aucun seuil de 5 % ne s'y applique. La cause d'une diff\u00e9rence d\u00e9pend de la fa\u00e7on dont l'\u00e9tude a \u00e9t\u00e9 men\u00e9e (tirage au sort, ordre, op\u00e9rateur).")
	if (length(obj$notes$en) && !isTRUE(narrated)) {
		say("Warnings:", "Avertissements :")
		for (i in seq_along(obj$notes$en)) say(paste0("  ", obj$notes$en[i]), paste0("  ", obj$notes$fr[i]))
	}
	invisible(obj)
}
