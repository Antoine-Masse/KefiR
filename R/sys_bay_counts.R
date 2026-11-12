# sys_bay_counts.R
# Counts with one factor (DECISIONS, D17, way 2), asked by response = "counts".
#   Two laws compared, not chosen: Poisson, and negative binomial for counts more spread
#   than Poisson allows (overdispersion). Each law with and without a difference according
#   to the factor: four models, glm(family = poisson) [stats] and glm.nb() [MASS].
#   Marginal likelihoods by the unit-information approximation, -BIC / 2 (Schwarz 1978;
#   Raftery 1995; Wagenmakers 2007): no compiler, no draws. The answer is averaged over the
#   two laws by their posterior probabilities, as for the continuous one-factor route.
#   The share of the negative binomial is the finding of overdispersion; the Pearson ratio
#   (chi-squared / degrees of freedom of the Poisson model) is said beside it.
#   Which groups differ: the partitions of D11, each grouping fitted under both laws.
#   How much: rate ratio of each pair, with a 95 % interval from the Laplace approximation
#   of the posterior (log rate ratio approximately normal), under the leading law.
#
# References.
#   Schwarz, G. (1978). Estimating the dimension of a model. Annals of Statistics, 6, 461-464.
#   Raftery, A. E. (1995). Bayesian model selection in social research. Sociological
#     Methodology, 25, 111-163.
#   Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values.
#     Psychonomic Bulletin & Review, 14, 779-804.

# Warnings of glm.nb (theta.ml, iteration limit) are counted, not printed: the console stays
# clean and the story says when the fit was unstable (audit of de64470, point 3).
.bay_cw <- new.env(parent = emptyenv())
.bay_cw$n <- 0L
.bay_quiet <- function(expr) {
	withCallingHandlers(expr, warning = function(w) {
		.bay_cw$n <- .bay_cw$n + 1L
		invokeRestart("muffleWarning")
	})
}

# Counts, chosen alone (DECISIONS, D17, way 2): non-negative whole numbers, at least 3
# distinct values, a group mean under 10 (the normal approximation of the Poisson law is
# poor below about 10), and zeros or more than 7 distinct values (which leaves out the
# rating scales 1 to 5 or 1 to 7, without zero).
.bay_is_counts <- function(y, g) {
	if (!is.numeric(y) || any(!is.finite(y))) return(FALSE)
	if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) return(FALSE)
	nu <- length(unique(y))
	if (nu < 3L) return(FALSE)
	if (min(tapply(y, g, mean)) >= 10) return(FALSE)
	any(y == 0) || nu > 7L
}

.bay_counts_fit <- function(y, h) {
	d <- data.frame(y = y, h = h)
	one <- nlevels(droplevels(as.factor(h))) < 2L
	fo <- if (one) y ~ 1 else y ~ h
	p <- .bay_quiet(stats::glm(fo, data = d, family = stats::poisson()))
	nb <- tryCatch(.bay_quiet(MASS::glm.nb(fo, data = d)), error = function(e) NULL)
	list(pois = p, nb = nb,
		lm_pois = -stats::BIC(p) / 2,
		lm_nb = if (is.null(nb) || !is.finite(nb$theta) || nb$theta > 1e6) NA_real_ else -stats::BIC(nb) / 2)
}

# Joint posterior over {law} x {no difference, difference}.
.bay_counts_post <- function(f0, f1, prior_h1) {
	lm <- c(pois0 = f0$lm_pois, pois1 = f1$lm_pois, nb0 = f0$lm_nb, nb1 = f1$lm_nb)
	lp <- log(c(0.5 * (1 - prior_h1), 0.5 * prior_h1, 0.5 * (1 - prior_h1), 0.5 * prior_h1))
	ok <- is.finite(lm)
	lw <- lm + lp
	w <- rep(0, 4L); names(w) <- names(lm)
	w[ok] <- exp(lw[ok] - .bay_log_sum_exp(lw[ok]))
	p_h1 <- w[["pois1"]] + w[["nb1"]]
	p_nb <- w[["nb0"]] + w[["nb1"]]
	# Bayes factor of a difference, averaged over the laws: posterior odds over prior odds.
	lb <- log(p_h1) - log1p(-p_h1) - (log(prior_h1) - log1p(-prior_h1))
	list(w = w, p_h1 = p_h1, p_nb = p_nb, log_bf = lb, nb_ok = all(ok))
}

.bay_counts_posthoc <- function(y, g, prior_h1, law, alphas = c(1, 1/2, 2)) {
	g <- droplevels(as.factor(g))
	lev <- levels(g); k <- length(lev)
	pairs <- utils::combn(k, 2L)
	res <- data.frame(group1 = lev[pairs[1L, ]], group2 = lev[pairs[2L, ]], stringsAsFactors = FALSE)
	means <- as.numeric(tapply(y, g, mean))
	parts_tab <- NULL
	if (k <= .bay_posthoc_max_groups()) {
		parts <- .bay_partitions(k)
		f0 <- .bay_counts_fit(y, factor(rep(1L, length(y))))
		lb <- vapply(parts, function(p) {
			if (max(p) == 1L) return(0)
			.bay_counts_post(f0, .bay_counts_fit(y, factor(p[as.integer(g)])), 0.5)$log_bf
		}, numeric(1))
		pm <- do.call(rbind, parts)
		post_of <- function(alpha) {
			lw <- lb + .bay_partition_prior(parts, prior_h1, alpha)
			exp(lw - .bay_log_sum_exp(lw))
		}
		sep_m <- pm[, pairs[1L, ], drop = FALSE] != pm[, pairs[2L, ], drop = FALSE]
		post <- post_of(alphas[1L])
		res$p_differ <- as.numeric(crossprod(sep_m, post))
		ps <- matrix(vapply(alphas, function(a) as.numeric(crossprod(sep_m, post_of(a))), numeric(ncol(pairs))), nrow = ncol(pairs))
		res$p_min <- apply(ps, 1L, min); res$p_max <- apply(ps, 1L, max)
		parts_tab <- data.frame(partition = vapply(parts, .bay_partition_label, character(1), lev = lev, means = means),
			code = vapply(parts, paste, character(1), collapse = ""), posterior = post,
			prior = exp(.bay_partition_prior(parts, prior_h1, alphas[1L])), log_bf = lb, stringsAsFactors = FALSE)
		parts_tab <- parts_tab[order(-parts_tab$posterior), , drop = FALSE]
		rownames(parts_tab) <- NULL
	} else {
		res$p_differ <- res$p_min <- res$p_max <- NA_real_
	}
	# Rate ratios under the leading law, Laplace approximation of the posterior.
	d <- data.frame(y = y, g = g)
	fit <- if (identical(law, "nb")) .bay_quiet(MASS::glm.nb(y ~ g - 1, data = d)) else stats::glm(y ~ g - 1, data = d, family = stats::poisson())
	b <- stats::coef(fit); V <- stats::vcov(fit)
	st <- t(vapply(seq_len(ncol(pairs)), function(j) {
		i1 <- pairs[1L, j]; i2 <- pairs[2L, j]
		est <- b[i2] - b[i1]; se <- sqrt(V[i1, i1] + V[i2, i2] - 2 * V[i1, i2])
		c(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
	}, numeric(3)))
	res$diff <- st[, 1L]; res$diff_low <- st[, 2L]; res$diff_high <- st[, 3L]
	res$delta <- res$delta_low <- res$delta_high <- NA_real_
	list(pairs = res, partitions = parts_tab, k = k, alphas = alphas, computed = !is.null(parts_tab),
		method = if (!is.null(parts_tab)) "exact" else "none", n_partitions = .bay_bell(k), sweeps = NA_integer_, mc_error = 0,
		max_exact = .bay_posthoc_max_groups(), delta_unit = "ratio", law = law, caution = FALSE)
}

.bay_m_counts <- function(input, g, y, call_expr, get_notes, note, prior_h1, sensitivity, detail, verbose, return, story, posthoc,
                          auto = FALSE) {
	.bay_cw$n <- 0L
	if (auto) {
		.bay_tell(story, "Change of strategy: count response (whole numbers from 0, zeros or many values, group means under 10); Poisson and negative binomial laws instead of the continuous ones.",
			"Changement de strat\u00e9gie : r\u00e9ponse de comptage (entiers \u00e0 partir de 0, z\u00e9ros ou valeurs nombreuses, moyennes sous 10) ; lois de Poisson et binomiale n\u00e9gative au lieu des lois continues.", "decide")
		.bay_tell(story, "Count detection: whole numbers >= 0, at least 3 distinct values, a group mean under 10 (normal approximation of Poisson poor below 10), zeros or more than 7 values (rating scales left out); response = \"continuous\" forces the continuous route.",
			"D\u00e9tection des comptages : entiers >= 0, au moins 3 valeurs distinctes, une moyenne de groupe sous 10 (approximation normale de Poisson m\u00e9diocre sous 10), des z\u00e9ros ou plus de 7 valeurs (\u00e9chelles de notation exclues) ; response = \u00ab continuous \u00bb force la voie continue.", "detail")
	}
	resp <- input$response_name
	fac <- names(input$factors)[1L]
	if (identical(fac, ".g")) fac <- "g"
	g <- droplevels(as.factor(g))
	k <- nlevels(g)
	n <- as.integer(table(g))
	mu <- as.numeric(tapply(y, g, mean)); vv <- as.numeric(tapply(y, g, stats::var))
	.bay_tell(story, paste0("Sizes: ", min(n), if (max(n) > min(n)) paste0(" to ", max(n)), " per group (N = ", length(y), ")."),
		paste0("Effectifs : ", min(n), if (max(n) > min(n)) paste0(" \u00e0 ", max(n)), " par groupe (N = ", length(y), ")."))
	z <- mean(y == 0)
	.bay_tell(story, paste0("Zeros: ", .bay_pct(z, "."), " of the counts."), paste0("Z\u00e9ros : ", .bay_pct(z, ","), " des comptages."))
	.bay_step_result(story, "counts kept as they are", "comptages gard\u00e9s tels quels")
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	.bay_tell(story, "Law of the counts", "Loi des comptages", "step", fun = c("glm(), glm.nb() {stats, MASS}", "glm(), glm.nb() {stats, MASS}"))
	f0 <- .bay_counts_fit(y, factor(rep(1L, length(y))))
	f1 <- .bay_counts_fit(y, g)
	po <- .bay_counts_post(f0, f1, prior_h1)
	phi <- sum(stats::residuals(f1$pois, type = "pearson")^2) / f1$pois$df.residual
	n_warn <- .bay_cw$n
	lead_nb <- isTRUE(po$p_nb >= 0.5)
	.bay_step_result(story, paste0(if (lead_nb) "negative binomial" else "Poisson", " ahead (", .bay_pct(if (lead_nb) po$p_nb else 1 - po$p_nb, "."), ")"),
		paste0(if (lead_nb) "binomiale n\u00e9gative" else "Poisson", " en t\u00eate (", .bay_pct(if (lead_nb) po$p_nb else 1 - po$p_nb, ","), ")"))
	phi_t <- .bay_num_txt(phi, 2)
	.bay_tell(story, paste0("Spread against Poisson (Pearson ratio): ", phi_t$en, " (1 under Poisson)."),
		paste0("Dispersion rapport\u00e9e \u00e0 Poisson (rapport de Pearson) : ", phi_t$fr, " (1 sous Poisson)."))
	.bay_tell(story, paste0("Probabilities: Poisson ", .bay_pct(1 - po$p_nb, "."), "; negative binomial ", .bay_pct(po$p_nb, "."), "."),
		paste0("Probabilit\u00e9s : Poisson ", .bay_pct(1 - po$p_nb, ","), " ; binomiale n\u00e9gative ", .bay_pct(po$p_nb, ","), "."))
	if (lead_nb) {
		.bay_tell(story, "Change of strategy: overdispersion, Poisson would overstate the evidence; the negative binomial leads.",
			"Changement de strat\u00e9gie : surdispersion, Poisson surestimerait les preuves ; la binomiale n\u00e9gative prend la t\u00eate.", "decide")
	} else {
		.bay_tell(story, "No overdispersion found: Poisson leads.", "Pas de surdispersion trouv\u00e9e : Poisson en t\u00eate.", "decide")
	}
	if (n_warn > 0L) {
		.bay_tell(story, paste0("Caution added: the negative binomial fit was unstable (", n_warn, if (n_warn == 1L) " warning" else " warnings", " of glm.nb); read the answer with caution."),
			paste0("Prudence ajout\u00e9e : l'ajustement de la binomiale n\u00e9gative a \u00e9t\u00e9 instable (", n_warn, if (n_warn == 1L) " avertissement" else " avertissements", " de glm.nb) ; lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	}
	.bay_tell(story, "the answer is weighted over the two laws by their probability.",
		"la r\u00e9ponse est pond\u00e9r\u00e9e sur les deux lois par leur probabilit\u00e9.", "note")
	.bay_tell(story, "Laws: Poisson, and negative binomial glm.nb() [MASS]; marginal likelihoods by -BIC/2, the unit-information approximation (Schwarz 1978; Raftery 1995; Wagenmakers 2007).",
		"Lois : Poisson, et binomiale n\u00e9gative glm.nb() [MASS] ; vraisemblances marginales par -BIC/2, approximation d'information unitaire (Schwarz 1978 ; Raftery 1995 ; Wagenmakers 2007).", "detail")

	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of a difference"), paste0(.bay_pct(prior_h1, ","), " de chances d'une diff\u00e9rence"))
	.bay_tell(story, "Each law starts at 50.0 %.", "Chaque loi part \u00e0 50,0 %.")
	.bay_tell(story, "prior on the effects: worth one observation (unit information), no scale to choose.",
		"a priori sur les effets : le poids d'une seule observation (information unitaire), sans \u00e9chelle \u00e0 choisir.", "note")

	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step")
	lb <- po$log_bf
	inv_en <- if (lb < 0) paste0(", i.e. ", .bay_bf_txt(-lb, "."), " in favour of no difference") else ""
	inv_fr <- if (lb < 0) paste0(", soit ", .bay_bf_txt(-lb, ","), " en faveur de l'absence de diff\u00e9rence") else ""
	.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb, "."), inv_en, "."), paste0("Facteur de Bayes : ", .bay_bf_txt(lb, ","), inv_fr, "."))
	.bay_tell(story, "Precision: exact computation, without draws.", "Pr\u00e9cision : calcul exact, sans tirage.")
	sens <- NULL
	if (isTRUE(sensitivity)) {
		phs <- c(1 / 3, 1 / 2, 2 / 3)
		grid <- vapply(phs, function(p) .bay_counts_post(f0, f1, p)$p_h1, numeric(1))
		lead_h1 <- po$p_h1 >= 0.5
		rng <- range(if (lead_h1) grid else 1 - grid)
		verdict <- .bay_sensitivity_verdict(max(po$p_h1, 1 - po$p_h1), rng)
		sens <- list(prior_h1 = phs, p = grid, range_lead = rng, verdict = verdict, verdict_terms = verdict, range_terms = matrix(range(grid), 1L))
		stb <- .bay_stab_txt(verdict)
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", ", stb$fr))
		.bay_tell(story, paste0("Sensitivity: the leading answer ", .bay_range_txt(rng)$en, "."),
			paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate ", .bay_range_txt(rng)$fr, "."))
		.bay_tell(story, if (identical(verdict, "changes")) "Caution added: the answer can reverse; it will be called not settled." else "The order does not change.",
			if (identical(verdict, "changes")) "Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e." else "L'ordre ne change pas.", "decide")
		.bay_tell(story, "the starting chances go from 1 against 2 to 2 against 1.", "les chances de d\u00e9part vont de 1 contre 2 \u00e0 2 contre 1.", "note")
	} else {
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", sensitivity not computed"), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", sensibilit\u00e9 non calcul\u00e9e"))
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).", "Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	ph <- NULL
	if (isTRUE(posthoc)) {
		ph <- .bay_counts_posthoc(y, g, prior_h1, if (lead_nb) "nb" else "pois")
		.bay_story_posthoc(story, ph)
	} else {
		.bay_story_no_posthoc(story, "off")
	}

	hyp <- data.frame(
		label_en = c(paste0(resp, " differs according to ", fac), paste0("No difference in ", resp, " according to ", fac)),
		label_fr = c(paste0(resp, " diff\u00e8re selon ", fac), paste0("Aucune diff\u00e9rence de ", resp, " selon ", fac)),
		posterior = c(po$p_h1, 1 - po$p_h1), prior = c(prior_h1, 1 - prior_h1), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(hypotheses = hyp,
		effects = data.frame(term = fac, bf10 = exp(lb), log_bf = lb, p_h1 = po$p_h1, p_h0 = 1 - po$p_h1, prior_h1 = prior_h1,
			analogue = "Poisson / negative binomial GLM", stringsAsFactors = FALSE),
		laws = data.frame(law = c("poisson", "negbin"), posterior = c(1 - po$p_nb, po$p_nb), stringsAsFactors = FALSE),
		models = po$w, dispersion = phi, sensitivity = sens, posthoc = ph, nb_warnings = n_warn,
		descriptive_stats = data.frame(group = levels(g), n = n, mean = mu, var = vv, stringsAsFactors = FALSE),
		priors = list(prior_h1 = prior_h1, effects = "unit information (BIC)", laws = "Poisson 1/2, negative binomial 1/2"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = po$p_h1,
		method = "Counts: Poisson and negative binomial compared (BIC approximation)", plan = "counts",
		label = input$label, n = length(y), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_counts(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(po$p_h1)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_print_bilan_counts <- function(obj, detail = FALSE, narrated = FALSE) {
	body <- function(emit, li, narrated) {
		ph <- obj$posthoc
		if (!is.null(ph) && ph$k > 2L) .bay_pairs_short(ph, emit, li$h$posterior[1L] < 0.6 || li$open, narrated)
		if (!is.null(ph) && ph$k == 2L) {
			pr <- ph$pairs
			f <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
			emit("fact", paste0("Rate ratio ", pr$group2, "/", pr$group1, ": ", f(pr$diff, "."), " [", f(pr$diff_low, "."), "; ", f(pr$diff_high, "."), "]."),
				paste0("Rapport des moyennes ", pr$group2, "/", pr$group1, " : ", f(pr$diff, ","), " [", f(pr$diff_low, ","), " ; ", f(pr$diff_high, ","), "]."))
		}
		emit("block", "", "")
		d <- obj$descriptive_stats
		f2 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
		emit("fact", paste0("Observed means: ", paste0(d$group, " = ", f2(d$mean, "."), collapse = "; "), "."),
			paste0("Moyennes observ\u00e9es : ", paste0(d$group, " = ", f2(d$mean, ","), collapse = " ; "), "."))
		emit("fact", paste0("Variances: ", paste0(d$group, " = ", f2(d$var, "."), collapse = "; "), " (equal to the means under Poisson)."),
			paste0("Variances : ", paste0(d$group, " = ", f2(d$var, ","), collapse = " ; "), " (\u00e9gales aux moyennes sous Poisson)."))
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)
	say(paste0("Dispersion (Pearson ratio of the Poisson model): ", format(signif(obj$dispersion, 3), trim = TRUE), "; law probabilities: Poisson ",
			.bay_pct(obj$laws$posterior[1L], "."), ", negative binomial ", .bay_pct(obj$laws$posterior[2L], "."), "."),
		paste0("Dispersion (rapport de Pearson du mod\u00e8le de Poisson) : ", gsub(".", ",", format(signif(obj$dispersion, 3), trim = TRUE), fixed = TRUE),
			" ; probabilit\u00e9s des lois : Poisson ", .bay_pct(obj$laws$posterior[1L], ","), ", binomiale n\u00e9gative ", .bay_pct(obj$laws$posterior[2L], ","), "."))
	say("Rate ratios: Laplace approximation of the posterior under the leading law (log ratio approximately normal).",
		"Rapports des moyennes : approximation de Laplace de la loi a posteriori sous la loi en t\u00eate (log du rapport approximativement normal).")
	invisible(obj)
}
