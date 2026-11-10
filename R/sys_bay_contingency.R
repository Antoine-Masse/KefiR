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
	cats <- rev(colnames(tab))
	for (i in seq_len(k)) {
		cn <- tab[i, cats]
		.bay_tell(story, paste0(rownames(tab)[i], ": ", paste0(cn, " ", cats, collapse = ", "), " (of ", sum(cn), ")."),
			paste0(rownames(tab)[i], " : ", paste0(cn, " ", cats, collapse = ", "), " (sur ", sum(cn), ")."))
	}
	if (any(tab == 0L)) {
		ne <- sum(tab == 0L)
		.bay_tell(story, paste0(ne, if (ne == 1L) " empty cell: no problem, no minimum count is needed." else " empty cells: no problem, no minimum count is needed."),
			paste0(ne, if (ne == 1L) " case vide : sans probl\u00e8me, aucun effectif minimum n'est demand\u00e9." else " cases vides : sans probl\u00e8me, aucun effectif minimum n'est demand\u00e9."))
	}
	.bay_step_result(story, "exact computation, no minimum count", "calcul exact, sans effectif minimum")
	.bay_tell(story, "small counts do not make the answer wrong, only less decisive (unlike the chi-squared).",
		"les petits effectifs ne faussent pas la r\u00e9ponse, ils la rendent moins tranch\u00e9e (contrairement au khi-deux).", "note")
	.bay_tell(story, "Law of the data", "Loi des donn\u00e9es", "step")
	.bay_step_result(story, if (m == 2L) "binomial" else "multinomial", if (m == 2L) "binomiale" else "multinomiale")
	.bay_tell(story, paste0("Each group has its own proportions of the ", m, " categories."),
		paste0("Chaque groupe a ses propres proportions des ", m, " modalit\u00e9s."))
	.bay_tell(story, "No difference: the same proportions in every group.",
		"Aucune diff\u00e9rence : les m\u00eames proportions dans tous les groupes.")
	.bay_tell(story, "no shape of the noise to compare with counts.", "pas de forme du bruit \u00e0 comparer avec des effectifs.", "note")
	if (!no_laws) {
		note("For a categorical response the laws of the noise do not apply; laws is ignored.",
			"Pour une r\u00e9ponse cat\u00e9gorielle, les lois du bruit ne s'appliquent pas ; laws est ignor\u00e9.")
	}

	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference"), paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence"))
	.bay_tell(story, "In each group, every set of proportions equally likely in advance.",
		"Dans chaque groupe, toutes les proportions sont \u00e9galement possibles d'avance.")
	.bay_tell(story, "these values are not tuned on the data; step 5 varies them.",
		"ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")
	.bay_tell(story, "Prior: uniform Dirichlet on the proportions of each group (Gunel and Dickey 1974, independent multinomials, group sizes fixed).",
		"A priori : Dirichlet uniforme sur les proportions de chaque groupe (Gunel et Dickey 1974, multinomiales ind\u00e9pendantes, effectifs des groupes fix\u00e9s).", "detail")

	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = "contingencyTableBF() {BayesFactor}")
	lb <- .bay_cont_bf_pkg(tab, 1)
	p_h1 <- .bay_cont_prob(lb, prior_h1)
	inv_en <- if (lb < 0) paste0(", i.e. ", .bay_bf_txt(-lb, "."), " in favour of no difference") else ""
	inv_fr <- if (lb < 0) paste0(", soit ", .bay_bf_txt(-lb, ","), " en faveur de l'absence de diff\u00e9rence") else ""
	.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb, "."), inv_en, "."), paste0("Facteur de Bayes : ", .bay_bf_txt(lb, ","), inv_fr, "."))
	.bay_tell(story, "The Bayes factor says how many times more probable the counts are with a difference than without.",
		"Le facteur de Bayes dit combien de fois les effectifs sont plus probables avec une diff\u00e9rence que sans.", "note")
	.bay_tell(story, "Precision: exact computation, without draws.", "Pr\u00e9cision : calcul exact, sans tirage.")
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
		stb <- .bay_stab_txt(verdict)
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", ", stb$fr))
		.bay_tell(story, paste0("Sensitivity: the leading answer goes from ", .bay_pct(rng[1L], "."), " to ", .bay_pct(rng[2L], "."), "."),
			paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate va de ", .bay_pct(rng[1L], ","), " \u00e0 ", .bay_pct(rng[2L], ","), "."))
		.bay_tell(story, "what varies is the expected spread of the proportions and the starting chances (1 against 2 to 2 against 1).",
			"on fait varier l'\u00e9talement attendu des proportions et les chances de d\u00e9part (1 contre 2 \u00e0 2 contre 1).", "note")
		.bay_tell(story, "Sensitivity: concentration a of the Dirichlet 1, 2 and 4, crossed with the prior odds 1:2, 1:1, 2:1.",
			"Sensibilit\u00e9 : concentration a de la Dirichlet 1, 2 et 4, crois\u00e9e avec les cotes a priori 1:2, 1:1, 2:1.", "detail")
		if (identical(verdict, "changes")) {
			.bay_tell(story, "Caution added: the answer can reverse; it will be called not settled.",
				"Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e.", "decide")
		} else {
			.bay_tell(story, "The answer does not change side.", "La r\u00e9ponse ne change pas de c\u00f4t\u00e9.", "decide")
		}
	} else {
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", sensitivity not computed"), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", sensibilit\u00e9 non calcul\u00e9e"))
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	ph_obj <- NULL
	if (isTRUE(posthoc)) {
		ph_obj <- .bay_cont_posthoc(tab, 1, prior_h1, draws, seed)
		.bay_story_posthoc(story, ph_obj)
	} else {
		.bay_story_no_posthoc(story, "off")
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
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_cont(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_cont_prob <- function(lb, prior_h1) 1 / (1 + exp(-lb) * (1 - prior_h1) / prior_h1)

.bay_print_bilan_cont <- function(obj, detail = FALSE, narrated = FALSE) {
	body <- function(emit, li, narrated) {
		ph <- obj$posthoc
		if (!is.null(ph) && ph$k > 2L) .bay_pairs_short(ph, emit, li$h$posterior[1L] < 0.6 || li$open, narrated)
		if (!is.null(ph) && ph$k == 2L) {
			pr <- ph$pairs
			pt <- function(v, dc) gsub(".", dc, formatC(100 * v, format = "f", digits = 1), fixed = TRUE)
			if (identical(ph$delta_unit, "logor")) {
				emit("fact", paste0("Estimated gap of the share of ", ph$category, ", ", pr$group2, " - ", pr$group1, ": ", pt(pr$diff, "."), " points [", pt(pr$diff_low, "."), "; ", pt(pr$diff_high, "."), "]."),
					paste0("\u00c9cart estim\u00e9 de la part de ", ph$category, ", ", pr$group2, " - ", pr$group1, " : ", pt(pr$diff, ","), " points [", pt(pr$diff_low, ","), " ; ", pt(pr$diff_high, ","), "]."))
			} else {
				emit("fact", paste0("Estimated total gap of the proportions: ", pt(pr$diff, "."), " points [", pt(pr$diff_low, "."), "; ", pt(pr$diff_high, "."), "]."),
					paste0("\u00c9cart total estim\u00e9 des proportions : ", pt(pr$diff, ","), " points [", pt(pr$diff_low, ","), " ; ", pt(pr$diff_high, ","), "]."))
			}
		}
		.bay_cont_observed(obj, emit)
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)
	say("Prior: in each group, uniform Dirichlet(1) on the proportions under a difference; under none, common proportions Dirichlet(R (a - 1) + 1) (contingencyTableBF, Gunel and Dickey 1974; Jamil et al. 2017).",
		"A priori : dans chaque groupe, Dirichlet(1) uniforme sur les proportions s'il y a une diff\u00e9rence ; sinon proportions communes Dirichlet(R (a - 1) + 1) (contingencyTableBF, Gunel et Dickey 1974 ; Jamil et al. 2017).")
	say("Sensitivity: concentration a of the Dirichlet from 1 to 4, crossed with the starting chances.",
		"Sensibilit\u00e9 : concentration a de la Dirichlet de 1 \u00e0 4, crois\u00e9e avec les chances de d\u00e9part.")
	say("Gaps: difference of the proportions (log odds ratio in $posthoc$pairs$delta) or total gap of the proportions; groupings in closed form (Dirichlet-multinomial).",
		"\u00c9carts : diff\u00e9rence des proportions (log du rapport de cotes dans $posthoc$pairs$delta) ou \u00e9cart total des proportions ; regroupements en forme close (Dirichlet-multinomiale).")
	invisible(obj)
}

# Observed counts, the category of interest first ("18 yes, 12 no (of 30)"), and the
# category that carries the largest gap (relecture D17, points 2 and 6).
.bay_cont_observed <- function(obj, emit) {
	emit("block", "", "")
	tab <- obj$table
	cats <- rev(colnames(tab))
	for (i in seq_len(nrow(tab))) {
		cn <- tab[i, cats]
		pc <- vapply(cn / sum(cn), .bay_pct, character(1), dec = ".")
		pf <- vapply(cn / sum(cn), .bay_pct, character(1), dec = ",")
		emit("fact", paste0(rownames(tab)[i], ": ", paste0(cn, " ", cats, " (", pc, ")", collapse = ", "), ", of ", sum(cn), "."),
			paste0(rownames(tab)[i], " : ", paste0(cn, " ", cats, " (", pf, ")", collapse = ", "), ", sur ", sum(cn), "."))
	}
	pr <- obj$proportions
	if (nrow(pr) >= 2L && ncol(pr) >= 2L) {
		rg <- apply(pr, 2L, function(v) diff(range(v)))
		cc <- names(which.max(rg))
		lo <- rownames(pr)[which.min(pr[, cc])]; hi <- rownames(pr)[which.max(pr[, cc])]
		emit("fact", paste0("Largest gap: ", cc, ", ", .bay_pct(min(pr[, cc]), "."), " in ", lo, " against ", .bay_pct(max(pr[, cc]), "."), " in ", hi, "."),
			paste0("\u00c9cart le plus marqu\u00e9 : ", cc, ", ", .bay_pct(min(pr[, cc]), ","), " en ", lo, " contre ", .bay_pct(max(pr[, cc]), ","), " en ", hi, "."))
	}
	invisible(NULL)
}
