# sys_bay_ordinal.R
# Ordered response with one factor (DECISIONS, D17, way 3): ratings 1 to 5, 1 to 7, 0 to 10.
#   Proportional odds (cumulative logit) model, MASS::polr(); each grouping against no
#   difference by the deviance Bayes factor with a Zellner-Siow g, r = sqrt(2)/4 (D31, extended
#   to one factor after passe 3 of the expert, I1), as for the counts (way 2): no draws.
#   Which groups differ: the partitions of D11, each grouping fitted by polr.
#   How much: cumulative odds ratio of each pair (above 1: higher values in the second group),
#   95 % interval from the Laplace approximation (log ratio approximately normal).
#   The proportional odds assumption is not tested; the story says so.
#
# References.
#   McCullagh, P. (1980). Regression models for ordinal data. JRSS B, 42, 109-142.
#   Venables, W. N., & Ripley, B. D. (2002). Modern Applied Statistics with S (MASS::polr).
#   Raftery, A. E. (1995). Bayesian model selection in social research. Sociological
#     Methodology, 25, 111-163.

# Usual tops of a rating or score scale.
.bay_scale_tops <- function() c(4, 5, 6, 7, 9, 10, 20, 100)

# Route of a whole-number response, chosen alone (D17, ways 2 and 3): "ordinal" (bounded
# scale with at most 10 values), "bounded" (bounded scale with more values, kept continuous
# and said), "counts" (no apparent bound), or "continuous".
.bay_int_route <- function(y, g) {
	if (!is.numeric(y) || any(!is.finite(y))) return("continuous")
	if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) return("continuous")
	nu <- length(unique(y))
	if (nu < 3L) return("continuous")
	if (max(y) %in% .bay_scale_tops()) return(if (nu <= 10L) "ordinal" else "bounded")
	if (.bay_is_counts(y, g)) return("counts")
	"continuous"
}

# Proportional odds checked (D26, audit of d7a295a): the polr model against a multinomial
# model of the same cells (nnet::multinom), where every level of the scale has its own effects,
# so without proportional odds. Both by -BIC/2, at even odds; returns P(not proportional), NA
# when nnet is missing or a fit fails. A comparison of the same kind as the Brant test, in the
# terms of this route (unit information).
.bay_po_check <- function(yo, h) {
	if (!requireNamespace("nnet", quietly = TRUE)) return(NA_real_)
	h <- droplevels(as.factor(h))
	if (nlevels(h) < 2L || nlevels(droplevels(yo)) < 3L) return(NA_real_)
	d <- data.frame(y = droplevels(yo), h = h)
	f1 <- tryCatch(.bay_quiet(MASS::polr(y ~ h, data = d)), error = function(e) NULL)
	f2 <- tryCatch(.bay_quiet(nnet::multinom(y ~ h, data = d, trace = FALSE)), error = function(e) NULL)
	if (is.null(f1) || is.null(f2)) return(NA_real_)
	l1 <- -stats::BIC(f1) / 2
	l2 <- -stats::BIC(f2) / 2
	if (!is.finite(l1) || !is.finite(l2)) return(NA_real_)
	1 / (1 + exp(l1 - l2))
}

# The story lines of that check.
.bay_po_tell <- function(story, p_np) {
	if (is.na(p_np)) {
		.bay_tell(story, "the gap between groups is assumed the same at every level of the scale (proportional odds, not checked here).",
			"on suppose que l'\u00e9cart entre groupes est le m\u00eame \u00e0 tous les niveaux de l'\u00e9chelle (cotes proportionnelles, non v\u00e9rifi\u00e9 ici).", "note")
	} else if (p_np >= 0.75) {
		.bay_tell(story, paste0("Caution added: the gaps between groups seem to change from one level of the scale to another (proportional odds doubtful, ", .bay_pct(1 - p_np, "."), "); read the answer with caution."),
			paste0("Prudence ajout\u00e9e : les \u00e9carts entre groupes semblent changer d'un niveau de l'\u00e9chelle \u00e0 l'autre (cotes proportionnelles douteuses, ", .bay_pct(1 - p_np, ","), ") ; lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(story, paste0("Proportional odds ", if (p_np <= 0.25) "plausible" else "not settled", " (", .bay_pct(1 - p_np, "."), ", against a model where every level has its own effects)."),
			paste0("Cotes proportionnelles ", if (p_np <= 0.25) "plausibles" else "non tranch\u00e9es", " (", .bay_pct(1 - p_np, ","), ", contre un mod\u00e8le o\u00f9 chaque niveau a ses propres effets)."))
	}
	.bay_tell(story, "Proportional odds checked by -BIC/2 of polr against nnet::multinom on the same groups or cells, at even odds (DECISIONS, D26).",
		"Cotes proportionnelles v\u00e9rifi\u00e9es par -BIC/2 de polr contre nnet::multinom sur les m\u00eames groupes ou cases, \u00e0 poids \u00e9gaux (DECISIONS, D26).", "detail")
	invisible(NULL)
}

.bay_polr_lm <- function(yo, h) {
	d <- data.frame(y = yo, h = h)
	one <- nlevels(droplevels(as.factor(h))) < 2L
	fit <- tryCatch(.bay_quiet(MASS::polr(if (one) y ~ 1 else y ~ h, data = d, Hess = TRUE)), error = function(e) NULL)
	ll <- if (is.null(fit)) NULL else stats::logLik(fit)
	list(fit = fit, lm = if (is.null(fit)) NA_real_ else -stats::BIC(fit) / 2,
		ll = if (is.null(ll)) NA_real_ else as.numeric(ll), k = if (is.null(ll)) NA_real_ else attr(ll, "df"), N = length(yo))
}

# Log Bayes factor of a grouping (fit f1) against no difference (fit f0): deviance Bayes factor
# with a Zellner-Siow g of scale r (D31 extended to one factor, I1 of the expert, passe 3).
.bay_polr_lbf <- function(f0, f1, r = .bay_dev_r()) {
	if (!is.finite(f0$ll) || !is.finite(f1$ll)) return(NA_real_)
	.bay_dev_lbf(2 * (f1$ll - f0$ll), f1$k - f0$k, f1$N, r)
}

.bay_ord_posthoc <- function(yo, g, prior_h1, alphas = c(1, 1/2, 2)) {
	g <- droplevels(as.factor(g))
	lev <- levels(g); k <- length(lev)
	pairs <- utils::combn(k, 2L)
	res <- data.frame(group1 = lev[pairs[1L, ]], group2 = lev[pairs[2L, ]], stringsAsFactors = FALSE)
	sc <- as.numeric(tapply(as.integer(yo), g, mean))
	parts_tab <- NULL
	if (k <= .bay_posthoc_max_groups()) {
		parts <- .bay_partitions(k)
		f00 <- .bay_polr_lm(yo, factor(rep(1L, length(yo))))
		lb <- vapply(parts, function(p) if (max(p) == 1L) 0 else .bay_polr_lbf(f00, .bay_polr_lm(yo, factor(p[as.integer(g)]))), numeric(1))
		lb[!is.finite(lb)] <- -Inf
		pm <- do.call(rbind, parts)
		post_of <- function(alpha) {
			lw <- lb + .bay_partition_prior(parts, prior_h1, alpha)
			exp(lw - .bay_log_sum_exp(lw[is.finite(lw)]))
		}
		sep_m <- pm[, pairs[1L, ], drop = FALSE] != pm[, pairs[2L, ], drop = FALSE]
		post <- post_of(alphas[1L])
		res$p_differ <- as.numeric(crossprod(sep_m, post))
		ps <- matrix(vapply(alphas, function(a) as.numeric(crossprod(sep_m, post_of(a))), numeric(ncol(pairs))), nrow = ncol(pairs))
		res$p_min <- apply(ps, 1L, min); res$p_max <- apply(ps, 1L, max)
		parts_tab <- data.frame(partition = vapply(parts, .bay_partition_label, character(1), lev = lev, means = sc),
			code = vapply(parts, paste, character(1), collapse = ""), posterior = post,
			prior = exp(.bay_partition_prior(parts, prior_h1, alphas[1L])), log_bf = lb, stringsAsFactors = FALSE)
		parts_tab <- parts_tab[order(-parts_tab$posterior), , drop = FALSE]
		rownames(parts_tab) <- NULL
	} else {
		res$p_differ <- res$p_min <- res$p_max <- NA_real_
	}
	f1 <- .bay_polr_lm(yo, g)$fit
	b <- c(0, stats::coef(f1)); V <- matrix(0, k, k)
	V[-1L, -1L] <- stats::vcov(f1)[names(stats::coef(f1)), names(stats::coef(f1)), drop = FALSE]
	st <- t(vapply(seq_len(ncol(pairs)), function(j) {
		i1 <- pairs[1L, j]; i2 <- pairs[2L, j]
		est <- b[i2] - b[i1]; se <- sqrt(max(V[i1, i1] + V[i2, i2] - 2 * V[i1, i2], 0))
		c(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
	}, numeric(3)))
	res$diff <- st[, 1L]; res$diff_low <- st[, 2L]; res$diff_high <- st[, 3L]
	res$delta <- res$delta_low <- res$delta_high <- NA_real_
	list(pairs = res, partitions = parts_tab, k = k, alphas = alphas, computed = !is.null(parts_tab),
		method = if (!is.null(parts_tab)) "exact" else "none", n_partitions = .bay_bell(k), sweeps = NA_integer_, mc_error = 0,
		max_exact = .bay_posthoc_max_groups(), delta_unit = "odds", caution = FALSE)
}

.bay_m_ordinal <- function(input, g, y, call_expr, get_notes, note, prior_h1, sensitivity, detail, verbose, return, story, posthoc,
                           auto = FALSE) {
	.bay_cw$n <- 0L
	resp <- input$response_name
	fac <- names(input$factors)[1L]
	if (identical(fac, ".g")) fac <- "g"
	g <- droplevels(as.factor(g))
	vals <- sort(unique(y))
	yo <- factor(y, levels = vals, ordered = TRUE)
	if (!is.null(input$ordered_levels)) {
		lv <- input$ordered_levels
		.bay_tell(story, paste0("Ordered response (ordered factor): ordinal model, only the order counts; ranks ", paste0(seq_along(lv), " = ", lv, collapse = ", "), "."),
			paste0("R\u00e9ponse ordonn\u00e9e (facteur ordonn\u00e9) : mod\u00e8le ordinal, seul compte l'ordre ; rangs ", paste0(seq_along(lv), " = ", lv, collapse = ", "), "."), "decide")
	} else if (auto) {
		.bay_tell(story, paste0("Change of strategy: bounded rating scale (", length(vals), " values, from ", min(vals), " to ", max(vals), "): ordinal model instead of the continuous laws."),
			paste0("Changement de strat\u00e9gie : \u00e9chelle de notation born\u00e9e (", length(vals), " valeurs, de ", min(vals), " \u00e0 ", max(vals), ") : mod\u00e8le ordinal au lieu des lois continues."), "decide")
		.bay_tell(story, "Ordinal route: whole numbers from 0, at most 10 values, maximum at a usual scale top (4, 5, 6, 7, 9, 10, 20, 100); response = \"continuous\" or \"counts\" forces another route.",
			"Voie ordinale : entiers \u00e0 partir de 0, au plus 10 valeurs, maximum \u00e0 un sommet d'\u00e9chelle usuel (4, 5, 6, 7, 9, 10, 20, 100) ; response = \u00ab continuous \u00bb ou \u00ab counts \u00bb force une autre voie.", "detail")
	}
	n <- as.integer(table(g))
	.bay_tell(story, paste0("Sizes: ", min(n), if (max(n) > min(n)) paste0(" to ", max(n)), " per group (N = ", length(y), ")."),
		paste0("Effectifs : ", min(n), if (max(n) > min(n)) paste0(" \u00e0 ", max(n)), " par groupe (N = ", length(y), ")."))
	.bay_step_result(story, "values kept as ordered categories", "valeurs gard\u00e9es comme cat\u00e9gories ordonn\u00e9es")
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	.bay_tell(story, "Law of the data", "Loi des donn\u00e9es", "step", fun = "polr() {MASS}")
	.bay_step_result(story, "proportional odds (cumulative logit)", "cotes proportionnelles (logit cumulatif)")
	.bay_tell(story, "Only the order of the values counts, not their gaps.", "Seul compte l'ordre des valeurs, pas leurs \u00e9carts.")
	p_np <- .bay_po_check(yo, g)
	.bay_po_tell(story, p_np)
	.bay_tell(story, "Proportional odds model of McCullagh (1980), MASS::polr().",
		"Mod\u00e8le \u00e0 cotes proportionnelles de McCullagh (1980), MASS::polr().", "detail")
	f0 <- .bay_polr_lm(yo, factor(rep(1L, length(y))))
	f1 <- .bay_polr_lm(yo, g)
	if (is.null(f0$fit) || is.null(f1$fit)) {
		.exit("The ordinal model could not be fitted (too many values for the sizes?); try response = \"continuous\".",
			"Le mod\u00e8le ordinal n'a pas pu \u00eatre ajust\u00e9 (trop de valeurs pour les effectifs ?) ; essayer response = \u00ab continuous \u00bb.")
	}
	if (.bay_cw$n > 0L) {
		.bay_tell(story, paste0("Caution added: the ordinal fit was unstable (", .bay_cw$n, if (.bay_cw$n == 1L) " warning" else " warnings", " of polr); read with caution."),
			paste0("Prudence ajout\u00e9e : l'ajustement ordinal a \u00e9t\u00e9 instable (", .bay_cw$n, if (.bay_cw$n == 1L) " avertissement" else " avertissements", " de polr) ; lire avec r\u00e9serve."), "decide")
	}

	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference"), paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence"))
	.bay_tell(story, paste0("Expected size r = ", .bay_num_txt(.bay_dev_r())$en, ", the same yardstick as for a continuous response."),
		paste0("Taille attendue r = ", .bay_num_txt(.bay_dev_r())$fr, ", le m\u00eame \u00e9talon que pour une r\u00e9ponse continue."))
	.bay_tell(story, "these values are not tuned on the data; step 5 varies them.", "ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")

	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step")
	lb <- .bay_polr_lbf(f0, f1)
	p_h1 <- 1 / (1 + exp(-lb) * (1 - prior_h1) / prior_h1)
	inv_en <- if (lb < 0) paste0(", i.e. ", .bay_bf_txt(-lb, "."), " in favour of no difference") else ""
	inv_fr <- if (lb < 0) paste0(", soit ", .bay_bf_txt(-lb, ","), " en faveur de l'absence de diff\u00e9rence") else ""
	.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb, "."), inv_en, "."), paste0("Facteur de Bayes : ", .bay_bf_txt(lb, ","), inv_fr, "."))
	.bay_tell(story, "Precision: exact computation, without draws.", "Pr\u00e9cision : calcul exact, sans tirage.")
	sens <- NULL
	if (isTRUE(sensitivity)) {
		phs <- c(1 / 3, 1 / 2, 2 / 3)
		r0 <- .bay_dev_r()
		lbs <- c(.bay_polr_lbf(f0, f1, r0 / 2), lb, .bay_polr_lbf(f0, f1, 2 * r0))
		grid <- as.vector(outer(lbs, phs, function(l, p) 1 / (1 + exp(-l) * (1 - p) / p)))
		rng <- range(if (p_h1 >= 0.5) grid else 1 - grid)
		verdict <- .bay_sensitivity_verdict(max(p_h1, 1 - p_h1), rng)
		sens <- list(prior_h1 = phs, p = grid, range_lead = rng, verdict = verdict, verdict_terms = verdict, range_terms = matrix(range(grid), 1L))
		stb <- .bay_stab_txt(verdict)
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", ", stb$fr))
		.bay_tell(story, paste0("Sensitivity: the leading answer ", .bay_range_txt(rng)$en, "."),
			paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate ", .bay_range_txt(rng)$fr, "."))
		.bay_tell(story, if (identical(verdict, "changes")) "Caution added: the answer can reverse; it will be called not settled." else "The order does not change.",
			if (identical(verdict, "changes")) "Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e." else "L'ordre ne change pas.", "decide")
		.bay_tell(story, "r goes from r/2 to 2r, the starting chances from 1 against 2 to 2 against 1.", "r va de r/2 \u00e0 2r, les chances de d\u00e9part de 1 contre 2 \u00e0 2 contre 1.", "note")
	} else {
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", sensitivity not computed"), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", sensibilit\u00e9 non calcul\u00e9e"))
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).", "Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	ph <- NULL
	if (isTRUE(posthoc)) {
		ph <- .bay_ord_posthoc(yo, g, prior_h1)
		.bay_story_posthoc(story, ph)
	} else {
		.bay_story_no_posthoc(story, "off")
	}

	hyp <- data.frame(
		label_en = c(paste0(resp, " differs according to ", fac), paste0("No difference in ", resp, " according to ", fac)),
		label_fr = c(paste0(resp, " diff\u00e8re selon ", fac), paste0("Aucune diff\u00e9rence de ", resp, " selon ", fac)),
		posterior = c(p_h1, 1 - p_h1), prior = c(prior_h1, 1 - prior_h1), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	tab <- table(g, yo); names(dimnames(tab)) <- NULL
	out <- list(hypotheses = hyp,
		effects = data.frame(term = fac, bf10 = exp(lb), log_bf = lb, p_h1 = p_h1, p_h0 = 1 - p_h1, prior_h1 = prior_h1,
			analogue = "proportional odds (polr)", stringsAsFactors = FALSE),
		table = tab, sensitivity = sens, posthoc = ph,
		descriptive_stats = data.frame(group = levels(g), n = n, median = as.numeric(tapply(y, g, stats::median)),
			mean = as.numeric(tapply(y, g, mean)), stringsAsFactors = FALSE),
		priors = list(prior_h1 = prior_h1, effects = "deviance Bayes factor, Zellner-Siow g, r = sqrt(2)/4"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = p_h1,
		method = "Ordered response: proportional odds (MASS::polr), deviance Bayes factor with a Zellner-Siow g", plan = "ordinal",
		label = input$label, n = length(y), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_ord(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_print_bilan_ord <- function(obj, detail = FALSE, narrated = FALSE) {
	body <- function(emit, li, narrated) {
		ph <- obj$posthoc
		if (!is.null(ph) && ph$k > 2L) .bay_pairs_short(ph, emit, li$h$posterior[1L] < 0.6 || li$open, narrated)
		if (!is.null(ph) && ph$k == 2L) {
			pr <- ph$pairs
			f <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
			emit("fact", paste0("Cumulative odds ratio ", pr$group2, "/", pr$group1, ": ", f(pr$diff, "."), " [", f(pr$diff_low, "."), "; ", f(pr$diff_high, "."), "]."),
				paste0("Rapport des cotes cumul\u00e9es ", pr$group2, "/", pr$group1, " : ", f(pr$diff, ","), " [", f(pr$diff_low, ","), " ; ", f(pr$diff_high, ","), "]."))
		}
		emit("block", "", "")
		d <- obj$descriptive_stats
		f1 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 1), fixed = TRUE)
		emit("fact", paste0("Observed medians: ", paste0(d$group, " = ", f1(d$median, "."), collapse = "; "), " (means ", paste(f1(d$mean, "."), collapse = "; "), ")."),
			paste0("M\u00e9dianes observ\u00e9es : ", paste0(d$group, " = ", f1(d$median, ","), collapse = " ; "), " (moyennes ", paste(f1(d$mean, ","), collapse = " ; "), ")."))
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)
	say("Cumulative odds ratio: above 1, the second group tends to higher values; Laplace approximation of the posterior.",
		"Rapport des cotes cumul\u00e9es : au-del\u00e0 de 1, le second groupe tend vers des valeurs plus hautes ; approximation de Laplace de la loi a posteriori.")
	invisible(obj)
}
