# sys_bay_cor.R
# Engine, story and bilan of cor.bay(), the Bayesian counterpart of cor.test() for two
# quantitative variables. Reasons and measures: dev/corbay/DECISIONS.md (C1 to C9).
#
# References (comments only):
#   Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor
#     hypothesis tests. Journal of Mathematical Psychology, 72, 19-32.
#     https://doi.org/10.1016/j.jmp.2015.06.004 (stretched beta prior on rho; correlationBF)
#   van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (2018). Bayesian inference for
#     Kendall's rank correlation coefficient. The American Statistician, 72(4), 303-308.
#     https://doi.org/10.1080/00031305.2016.1264998 (eq. 6, 9, 11-13: T* normal with mean
#     3/2 tau sqrt(n) and sd 1; prior on tau induced by Greiner's relation)
#   Makowski, D., Ben-Shachar, M. S., Patil, I., & Luedecke, D. (2020). Methods and algorithms
#     for correlation analysis in R. JOSS, 5(51), 2306. https://doi.org/10.21105/joss.02306
#     (Bayesian Spearman as correlationBF on ranks)
#   Kass, R. E., & Raftery, A. E. (1995). Bayes factors. JASA, 90(430), 773-795.

# ---------------------------------------------------------------------------
# Bayes factors: log BF against "no link" for rho > 0, rho < 0 and both
# ---------------------------------------------------------------------------

# Two halves from correlationBF (nullInterval), the two-sided BF as their mean
# (equal prior mass on each half). Error of the integration kept, in %.
.corbay_lbf_pearson <- function(x, y, kappa) {
	b <- BayesFactor::correlationBF(y = y, x = x, rscale = kappa, nullInterval = c(0, 1))
	e <- BayesFactor::extractBF(b, logbf = TRUE)
	l <- e$bf
	m <- max(l)
	c(plus = l[1L], minus = l[2L], two = m + log((exp(l[1L] - m) + exp(l[2L] - m)) / 2),
		error_pct = if (all(is.na(e$error))) NA_real_ else 100 * max(e$error, na.rm = TRUE))
}

.corbay_lbf_spearman <- function(x, y, kappa) .corbay_lbf_pearson(rank(x), rank(y), kappa)

# Standardised Kendall statistic T* (van Doorn et al. 2018, eq. 6), from tau-b.
.corbay_tstar <- function(x, y) {
	n <- length(x)
	tau <- suppressWarnings(stats::cor(x, y, method = "kendall"))
	tau * (n * (n - 1) / 2) / sqrt(n * (n - 1) * (2 * n + 5) / 18)
}

# Prior on tau induced by a stretched beta of width kappa on rho (eq. 9).
.corbay_prior_tau <- function(tau, kappa) {
	a <- 1 / kappa
	rho <- sin(pi * tau / 2)
	stats::dbeta((rho + 1) / 2, a, a) / 2 * (pi / 2) * cos(pi * tau / 2)
}

# Integrand of the marginal likelihood under H1, scaled by the likelihood under H0 (eq. 13).
.corbay_kendall_f <- function(ts, n, kappa) {
	m <- 1.5 * sqrt(n)
	l0 <- stats::dnorm(ts, 0, 1, log = TRUE)
	function(t) exp(stats::dnorm(ts, m * t, 1, log = TRUE) - l0) * .corbay_prior_tau(t, kappa)
}

.corbay_lbf_kendall <- function(x, y, kappa) {
	n <- length(x)
	ts <- .corbay_tstar(x, y)
	if (!is.finite(ts)) return(c(plus = NA_real_, minus = NA_real_, two = NA_real_, error_pct = NA_real_))
	f <- .corbay_kendall_f(ts, n, kappa)
	m <- 1.5 * sqrt(n)
	# The integrand is a narrow bump near T*/m for large n: the range is cut around it.
	pk <- min(max(ts / m, -1), 1)
	w <- 12 / m
	brk <- sort(unique(c(-1, max(-1, pk - w), 0, min(1, pk + w), 1)))
	err <- 0
	seg <- function(a, b) {
		if (b <= a) return(0)
		r <- stats::integrate(f, a, b, rel.tol = 1e-8, subdivisions = 500L)
		err <<- err + r$abs.error
		r$value
	}
	ints <- mapply(seg, brk[-length(brk)], brk[-1L])
	pos <- sum(ints[brk[-1L] > 0])
	neg <- sum(ints[brk[-length(brk)] < 0])
	c(plus = log(2 * pos), minus = log(2 * neg), two = log(pos + neg),
		error_pct = 100 * err / max(pos + neg, 1e-300))
}

.corbay_lbf <- function(method, x, y, kappa) {
	switch(method,
		pearson = .corbay_lbf_pearson(x, y, kappa),
		spearman = .corbay_lbf_spearman(x, y, kappa),
		kendall = .corbay_lbf_kendall(x, y, kappa))
}

# ---------------------------------------------------------------------------
# Hypotheses
# ---------------------------------------------------------------------------

# Prior: "no link" 1 - prior_h1; a link, prior_h1, split equally between its two signs
# (two-sided) or given to the one sign tested (alternative = "greater" or "less").
.corbay_hyp <- function(lbf, alternative, prior_h1, xn, yn) {
	keys <- switch(alternative, two.sided = c("none", "pos", "neg"), greater = c("none", "pos"), less = c("none", "neg"))
	pri <- c(none = 1 - prior_h1, pos = prior_h1, neg = prior_h1)
	if (identical(alternative, "two.sided")) pri[c("pos", "neg")] <- prior_h1 / 2
	pri <- pri[keys]
	ll <- c(none = 0, pos = unname(lbf["plus"]), neg = unname(lbf["minus"]))[keys]
	w <- log(pri) + ll
	post <- exp(w - max(w))
	post <- post / sum(post)
	lab_en <- c(none = paste0("No link between ", xn, " and ", yn), pos = paste0("Positive link between ", xn, " and ", yn),
		neg = paste0("Negative link between ", xn, " and ", yn))
	lab_fr <- c(none = paste0("Aucun lien entre ", xn, " et ", yn), pos = paste0("Lien positif entre ", xn, " et ", yn),
		neg = paste0("Lien n\u00e9gatif entre ", xn, " et ", yn))
	h <- data.frame(key = keys, label_en = unname(lab_en[keys]), label_fr = unname(lab_fr[keys]),
		posterior = unname(post), prior = unname(pri), log_bf = unname(ll), null = keys == "none", stringsAsFactors = FALSE)
	h <- h[order(-h$posterior), , drop = FALSE]
	rownames(h) <- NULL
	h
}

# ---------------------------------------------------------------------------
# Checks on the data and choice of the method
# ---------------------------------------------------------------------------

# Far values per variable (Tukey's rule, 1.5 and 3 interquartile ranges), influential pairs
# (Cook's distance), ties and distinct values. Descriptive only.
.corbay_checks <- function(x, y, seed) {
	tk <- function(v) {
		q <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
		i <- q[2L] - q[1L]
		if (!isTRUE(i > 0)) return(c(mild = 0L, far = 0L))
		c(mild = sum(v < q[1L] - 1.5 * i | v > q[2L] + 1.5 * i), far = sum(v < q[1L] - 3 * i | v > q[2L] + 3 * i))
	}
	far <- function(v) {
		q <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
		i <- q[2L] - q[1L]
		if (!isTRUE(i > 0)) return(rep(FALSE, length(v)))
		v < q[1L] - 3 * i | v > q[2L] + 3 * i
	}
	n <- length(x)
	# Influential pairs: Cook's distance above 1 in the regression of y on x or of x on y,
	# sought from 20 pairs (DECISIONS, C3: a robust distance from MASS::cov.rob flagged 40 to
	# 65 % of normal data sets; Cook's distance 4 % at n = 20, under 1 % from n = 30).
	infl <- NULL
	if (n >= 20L) {
		infl <- tryCatch({
			c1 <- stats::cooks.distance(stats::lm(y ~ x))
			c2 <- stats::cooks.distance(stats::lm(x ~ y))
			(is.finite(c1) & c1 > 1) | (is.finite(c2) & c2 > 1)
		}, error = function(e) NULL)
	}
	biv <- if (is.null(infl)) NA_integer_ else sum(infl)
	apart <- far(x) | far(y) | (if (is.null(infl)) FALSE else infl)
	list(n = n, tk_x = tk(x), tk_y = tk(y), influential = biv, apart = sum(apart),
		ties_x = 1 - length(unique(x)) / n, ties_y = 1 - length(unique(y)) / n,
		k_x = length(unique(x)), k_y = length(unique(y)))
}

# The rule of method = "auto" (DECISIONS, C3): ranks (Kendall) as soon as a few values could
# carry Pearson's answer or the scale is a short one; Pearson otherwise.
.corbay_choose <- function(ck) {
	why <- character(0)
	if (ck$tk_x[["far"]] + ck$tk_y[["far"]] > 0L) why <- c(why, "far")
	if (isTRUE(ck$influential > 0L)) why <- c(why, "influence")
	if (min(ck$k_x, ck$k_y) <= 7L) why <- c(why, "scale")
	list(method = if (length(why)) "kendall" else "pearson", why = why)
}

# ---------------------------------------------------------------------------
# Estimate of the link, with its credible interval
# ---------------------------------------------------------------------------

# Under the hypothesis of a link (both signs, or the sign tested). Pearson and Spearman:
# draws of rho from BayesFactor; Kendall: the posterior of tau on a grid, then draws.
.corbay_estimate <- function(method, x, y, kappa, alternative, draws, seed, ci) {
	ni <- switch(alternative, two.sided = NULL, greater = c(0, 1), less = c(-1, 0))
	if (method %in% c("pearson", "spearman")) {
		xx <- if (method == "spearman") rank(x) else x
		yy <- if (method == "spearman") rank(y) else y
		v <- .bay_with_seed(seed, {
			p <- NULL
			utils::capture.output(p <- suppressMessages(BayesFactor::correlationBF(y = yy, x = xx, rscale = kappa, nullInterval = ni,
				posterior = TRUE, iterations = draws, progress = FALSE)))
			as.numeric(p[, "rho"])
		})
	} else {
		n <- length(x)
		ts <- .corbay_tstar(x, y)
		f <- .corbay_kendall_f(ts, n, kappa)
		m <- 1.5 * sqrt(n)
		pk <- min(max(ts / m, -1), 1)
		lo <- max(-1, pk - 12 / m); hi <- min(1, pk + 12 / m)
		if (identical(alternative, "greater")) lo <- max(lo, 0)
		if (identical(alternative, "less")) hi <- min(hi, 0)
		if (hi <= lo) { lo <- if (identical(alternative, "greater")) 0 else -1; hi <- if (identical(alternative, "less")) 0 else 1 }
		g <- seq(lo, hi, length.out = 4001L)
		d <- f(g)
		d[!is.finite(d)] <- 0
		st <- g[2L] - g[1L]
		v <- .bay_with_seed(seed, sample(g, draws, replace = TRUE, prob = d) + stats::runif(draws, -st / 2, st / 2))
		v <- pmin(pmax(v, -1), 1)
	}
	h <- bayestestR::hdi(v, ci = ci)
	list(param = if (method == "kendall") "tau" else "rho", median = stats::median(v),
		lo = as.numeric(h$CI_low), hi = as.numeric(h$CI_high), ci = ci)
}

# ---------------------------------------------------------------------------
# Sensitivity: kappa / 2, kappa, 2 kappa; starting odds 1:2, 1:1, 2:1 (as m.bay, D4)
# ---------------------------------------------------------------------------

.corbay_sensitivity <- function(method, x, y, kappa, alternative, prior_h1, lead, xn, yn, lbf_central) {
	ks <- c(kappa / 2, kappa, 2 * kappa)
	odds <- c(1 / 2, 1, 2)
	o0 <- prior_h1 / (1 - prior_h1)
	grid <- matrix(NA_real_, 3L, 3L, dimnames = list(kappa = format(signif(ks, 3)), odds = c("1:2", "1:1", "2:1")))
	lead_same <- TRUE
	for (i in seq_along(ks)) {
		l <- if (i == 2L) lbf_central else .corbay_lbf(method, x, y, ks[i])
		for (j in seq_along(odds)) {
			o <- o0 * odds[j]
			h <- .corbay_hyp(l, alternative, o / (1 + o), xn, yn)
			grid[i, j] <- h$posterior[h$key == lead]
			if (!identical(h$key[1L], lead)) lead_same <- FALSE
		}
	}
	rg <- range(grid)
	verdict <- if (!lead_same) "changes" else if (.bay_support(rg[1L])$level != .bay_support(rg[2L])$level) "crosses" else "stable"
	list(grid = grid, kappa = ks, range_lead = rg, verdict = verdict)
}

# ---------------------------------------------------------------------------
# Words
# ---------------------------------------------------------------------------

.corbay_method_words <- function(method) {
	switch(method,
		pearson = list(en = "Pearson's correlation", fr = "corr\u00e9lation de Pearson", par = "rho", what_en = "a straight-line link",
			what_fr = "un lien en ligne droite", cite_en = "Bayesian Pearson correlation, Ly et al. 2016",
			cite_fr = "corr\u00e9lation de Pearson bay\u00e9sienne, Ly et al. 2016", fun = "correlationBF() {BayesFactor}"),
		spearman = list(en = "Spearman's correlation (Pearson on the ranks)", fr = "corr\u00e9lation de Spearman (Pearson sur les rangs)", par = "rho",
			what_en = "a monotone link (always in the same direction, straight or curved)",
			what_fr = "un lien monotone (toujours dans le m\u00eame sens, droit ou courbe)",
			cite_en = "Bayesian Pearson correlation on the ranks, Makowski et al. 2020",
			cite_fr = "corr\u00e9lation de Pearson bay\u00e9sienne sur les rangs, Makowski et al. 2020", fun = "correlationBF() {BayesFactor}"),
		kendall = list(en = "Kendall's tau", fr = "tau de Kendall", par = "tau",
			what_en = "a monotone link (always in the same direction, straight or curved)",
			what_fr = "un lien monotone (toujours dans le m\u00eame sens, droit ou courbe)",
			cite_en = "Bayesian Kendall's tau, van Doorn et al. 2018",
			cite_fr = "tau de Kendall bay\u00e9sien, van Doorn et al. 2018", fun = "cor.bay() {KefiR}"))
}

.corbay_num <- function(x, dg = 2L, dc = ".") .bay_fx(x, dg, dc)

.corbay_est_txt <- function(est, dc) {
	paste0(.corbay_num(est$median, 2L, dc), " [", .corbay_num(est$lo, 2L, dc), if (dc == ",") " ; " else "; ", .corbay_num(est$hi, 2L, dc), "]")
}

# ---------------------------------------------------------------------------
# Story, steps 1 to 6
# ---------------------------------------------------------------------------

.corbay_story_plan <- function(st, xn, yn, n, alternative) {
	.bay_tell(st, paste0("Plan recognised: two quantitative variables, ", xn, " and ", yn, "; N = ", n, " pairs."),
		paste0("Plan reconnu : deux variables quantitatives, ", xn, " et ", yn, " ; N = ", n, " paires."), "step")
	q <- switch(alternative,
		two.sided = list(en = paste0("One question: are ", xn, " and ", yn, " linked, and in which direction?"),
			fr = paste0("Une question : ", xn, " et ", yn, " sont-ils li\u00e9s, et dans quel sens ?")),
		greater = list(en = paste0("One question: do ", xn, " and ", yn, " rise together? (alternative = \"greater\")"),
			fr = paste0("Une question : ", xn, " et ", yn, " augmentent-ils ensemble ? (alternative = \u00ab greater \u00bb)")),
		less = list(en = paste0("One question: does ", yn, " fall when ", xn, " rises? (alternative = \"less\")"),
			fr = paste0("Une question : ", yn, " diminue-t-il quand ", xn, " augmente ? (alternative = \u00ab less \u00bb)")))
	.bay_tell(st, q$en, q$fr, "decide")
	.bay_tell(st, "a link observed without an experiment does not show a cause.",
		"un lien observ\u00e9 sans exp\u00e9rience ne montre pas une cause.", "note")
}

.corbay_story_checks <- function(st, ck, xn, yn, n_removed) {
	.bay_tell(st, "Checks on the data.", "Contr\u00f4les des donn\u00e9es.", "step", fun = "quantile() {stats}")
	.bay_tell(st, paste0("Pairs: ", ck$n, " complete", if (n_removed > 0L) paste0(", ", n_removed, " with a missing value removed") else "", "."),
		paste0("Paires : ", ck$n, " compl\u00e8tes", if (n_removed > 0L) paste0(", ", n_removed, " avec une valeur manquante retir\u00e9e", if (n_removed > 1L) "s" else "") else "", "."))
	pc <- function(v) format(round(100 * v), trim = TRUE)
	if (ck$ties_x > 0 || ck$ties_y > 0) {
		.bay_tell(st, paste0("Ties: ", pc(ck$ties_x), " % of the values of ", xn, ", ", pc(ck$ties_y), " % of ", yn, " repeat an earlier one (",
				ck$k_x, " and ", ck$k_y, " distinct values)."),
			paste0("Ex-aequo : ", pc(ck$ties_x), " % des valeurs de ", xn, ", ", pc(ck$ties_y), " % de ", yn, " r\u00e9p\u00e8tent une valeur d\u00e9j\u00e0 vue (",
				ck$k_x, " et ", ck$k_y, " valeurs distinctes)."))
	} else {
		.bay_tell(st, "Ties: none.", "Ex-aequo : aucun.")
	}
	.bay_tell(st, paste0("Far values: ", xn, " ", ck$tk_x[["far"]], ", ", yn, " ", ck$tk_y[["far"]], " beyond 3 interquartile ranges (",
			ck$tk_x[["mild"]], " and ", ck$tk_y[["mild"]], " beyond 1.5)."),
		paste0("Valeurs \u00e9loign\u00e9es : ", xn, " ", ck$tk_x[["far"]], ", ", yn, " ", ck$tk_y[["far"]], " au-del\u00e0 de 3 \u00e9carts interquartiles (",
			ck$tk_x[["mild"]], " et ", ck$tk_y[["mild"]], " au-del\u00e0 de 1,5)."))
	if (!is.na(ck$influential)) {
		.bay_tell(st, paste0("Influential pairs: ", ck$influential, " (Cook's distance above 1)."),
			paste0("Paires influentes : ", ck$influential, " (distance de Cook au-del\u00e0 de 1)."), fun = "cooks.distance() {stats}")
	} else {
		.bay_tell(st, "Influential pairs: not sought under 20 pairs.",
			"Paires influentes : non cherch\u00e9es sous 20 paires.")
	}
	.bay_tell(st, "Far values: Tukey's rule on each variable. Influential pairs: Cook's distance above 1 in the regression of y on x or of x on y (Cook and Weisberg 1982).",
		"Valeurs \u00e9loign\u00e9es : r\u00e8gle de Tukey sur chaque variable. Paires influentes : distance de Cook au-del\u00e0 de 1 dans la r\u00e9gression de y sur x ou de x sur y (Cook et Weisberg 1982).", "detail")
	nf <- ck$apart
	if (ck$n < 10L) {
		.bay_tell(st, "Few pairs: the answer will lean on the starting point more than on the data.",
			"Peu de paires : la r\u00e9ponse s'appuiera sur le point de d\u00e9part plus que sur les donn\u00e9es.", "decide")
	}
	.bay_step_result(st, if (nf > 0L) paste0(nf, if (nf == 1L) " pair apart, kept" else " pairs apart, kept") else "nothing to remove",
		if (nf > 0L) paste0(nf, if (nf == 1L) " paire \u00e0 part, gard\u00e9e" else " paires \u00e0 part, gard\u00e9es") else "rien \u00e0 retirer")
}

.corbay_story_method <- function(st, method, auto, why, xn, yn) {
	mw <- .corbay_method_words(method)
	.bay_tell(st, "Choice of the method.", "Choix de la m\u00e9thode.", "step", fun = "cor.bay() {KefiR}")
	why_en <- c(far = "far values", influence = "influential pairs", scale = "a short scale (7 values or fewer)")
	why_fr <- c(far = "valeurs \u00e9loign\u00e9es", influence = "paires influentes", scale = "une \u00e9chelle courte (7 valeurs ou moins)")
	if (auto) {
		if (identical(method, "kendall")) {
			.bay_tell(st, paste0("Seen at step 2: ", paste(why_en[why], collapse = ", "), "; a few values could carry Pearson's answer."),
				paste0("Vu \u00e0 l'\u00e9tape 2 : ", paste(why_fr[why], collapse = ", "), " ; quelques valeurs pourraient porter la r\u00e9ponse de Pearson."))
			.bay_tell(st, "Change of strategy: Kendall's tau, on the ranks, in place of Pearson's correlation.",
				"Changement de strat\u00e9gie : tau de Kendall, sur les rangs, \u00e0 la place de la corr\u00e9lation de Pearson.", "decide")
		} else {
			.bay_tell(st, "Neither far values, nor influential pairs, nor a short scale.",
				"Ni valeurs \u00e9loign\u00e9es, ni paires influentes, ni \u00e9chelle courte.")
			.bay_tell(st, "Pearson's correlation, the most precise when no value stands apart.",
				"Corr\u00e9lation de Pearson, la plus pr\u00e9cise quand aucune valeur ne se d\u00e9tache.", "decide")
		}
		.bay_tell(st, "the other method is computed too, as a cross-check, at step 5.",
			"l'autre m\u00e9thode est calcul\u00e9e aussi, en contre-\u00e9preuve, \u00e0 l'\u00e9tape 5.", "note")
	} else {
		.bay_tell(st, paste0("Method asked for: ", mw$en, "."), paste0("M\u00e9thode demand\u00e9e : ", mw$fr, "."), "decide")
		if (method == "pearson" && length(why)) {
			.bay_tell(st, paste0("Caution: ", paste(why_en[why], collapse = ", "), "; method = \"auto\" would take Kendall's tau."),
				paste0("Prudence : ", paste(why_fr[why], collapse = ", "), " ; method = \u00ab auto \u00bb prendrait le tau de Kendall."))
		}
	}
	.bay_tell(st, paste0("the method measures ", mw$what_en, "."), paste0("la m\u00e9thode mesure ", mw$what_fr, "."), "note_first")
	.bay_tell(st, "Rule of method = \"auto\": Kendall's tau if a variable has a value beyond 3 interquartile ranges, if a pair has a Cook's distance above 1 (from 20 pairs), or if a variable has 7 distinct values or fewer; Pearson otherwise. No p-value chooses the method.",
		"R\u00e8gle de method = \u00ab auto \u00bb : tau de Kendall si une variable a une valeur au-del\u00e0 de 3 \u00e9carts interquartiles, si une paire a une distance de Cook au-del\u00e0 de 1 (d\u00e8s 20 paires), ou si une variable a 7 valeurs distinctes ou moins ; Pearson sinon. Aucune p-value ne choisit la m\u00e9thode.", "detail")
	.bay_step_result(st, mw$en, mw$fr)
}

.corbay_story_prior <- function(st, kappa, prior_h1, alternative, method) {
	dc <- function(x, d) .bay_pct(x, d)
	res_en <- paste0(dc(prior_h1, "."), " chances of a link", if (alternative == "two.sided") paste0(" (", dc(prior_h1 / 2, "."), " positive, ", dc(prior_h1 / 2, "."), " negative)") else "")
	res_fr <- paste0(dc(prior_h1, ","), " de chances d'un lien", if (alternative == "two.sided") paste0(" (", dc(prior_h1 / 2, ","), " positif, ", dc(prior_h1 / 2, ","), " n\u00e9gatif)") else "")
	.bay_tell(st, "Starting point (priors).", "Point de d\u00e9part (a priori).", "step")
	.bay_step_result(st, res_en, res_fr)
	a <- 1 / kappa
	big <- 2 * stats::pbeta(0.25, a, a)
	kt <- .bay_num_txt(kappa, 3)
	.bay_tell(st, paste0("Expected strength: width kappa = ", kt$en, " on the correlation; a strong link (|rho| > 0.5) starts at ", dc(big, "."), "."),
		paste0("Force attendue : largeur kappa = ", kt$fr, " sur la corr\u00e9lation ; un lien fort (|rho| > 0,5) part \u00e0 ", dc(big, ","), "."))
	if (method == "kendall") {
		.bay_tell(st, "The same prior, carried over to Kendall's tau by Greiner's relation tau = 2/pi arcsin(rho).",
			"Le m\u00eame a priori, report\u00e9 sur le tau de Kendall par la relation de Greiner tau = 2/pi arcsin(rho).")
	}
	.bay_tell(st, "these values are not adjusted on the data; step 5 varies them.",
		"ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")
	.bay_tell(st, paste0("The prior on rho: stretched beta, (rho + 1) / 2 ~ Beta(1/kappa, 1/kappa), kappa = ", kt$en, " (\"medium\" of BayesFactor; JASP takes kappa = 1, uniform)."),
		paste0("L'a priori sur rho : b\u00eata \u00e9tir\u00e9e, (rho + 1) / 2 ~ B\u00eata(1/kappa, 1/kappa), kappa = ", kt$fr, " (\u00ab medium \u00bb de BayesFactor ; JASP prend kappa = 1, uniforme)."), "detail")
}

.corbay_story_solidity <- function(st, obj) {
	h <- obj$hypotheses
	mw <- .corbay_method_words(obj$method)
	lb <- obj$bf$log_bf10
	ll <- .corbay_lead_info(obj)$lb
	.bay_tell(st, "Solidity of the answer.", "Solidit\u00e9 de la r\u00e9ponse.", "step", fun = mw$fun)
	stb <- .bay_stab_txt(obj$sensitivity$verdict)
	.bay_step_result(st, paste0("Bayes factor ", .bay_bf_txt(ll, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(ll, ","), ", ", stb$fr))
	sg_en <- c(pos = "a positive link", neg = "a negative link", none = "no link")
	sg_fr <- c(pos = "un lien positif", neg = "un lien n\u00e9gatif", none = "aucun lien")
	.bay_tell(st, paste0("Bayes factor: ", .bay_bf_txt(ll, "."), " for ", sg_en[h$key[1L]], " against ", sg_en[h$key[2L]],
			if (obj$alternative == "two.sided") paste0(" (", .bay_bf_txt(lb, "."), " for a link of either sign against none).") else "."),
		paste0("Facteur de Bayes : ", .bay_bf_txt(ll, ","), " pour ", sg_fr[h$key[1L]], " contre ", sg_fr[h$key[2L]],
			if (obj$alternative == "two.sided") paste0(" (", .bay_bf_txt(lb, ","), " pour un lien de l'un ou l'autre sens contre aucun).") else "."))
	o <- obj$other
	if (!is.null(o)) {
		ow <- .corbay_method_words(o$method)
		same <- identical(o$hypotheses$key[1L], h$key[1L])
		.bay_tell(st, paste0("Cross-check, ", ow$en, ": ", tolower(substr(o$hypotheses$label_en[1L], 1, 1)), substring(o$hypotheses$label_en[1L], 2), " ",
				.bay_pct(o$hypotheses$posterior[1L], "."), if (same) ", the same answer." else ", another answer."),
			paste0("Contre-\u00e9preuve, ", ow$fr, " : ", tolower(substr(o$hypotheses$label_fr[1L], 1, 1)), substring(o$hypotheses$label_fr[1L], 2), " ",
				.bay_pct(o$hypotheses$posterior[1L], ","), if (same) ", m\u00eame r\u00e9ponse." else ", autre r\u00e9ponse."), fun = ow$fun)
		if (!same && isTRUE(o$explained)) {
			.bay_tell(st, "Pearson's answer can be carried by the pairs apart seen at step 2; it is not kept.",
				"La r\u00e9ponse de Pearson peut \u00eatre port\u00e9e par les paires \u00e0 part vues \u00e0 l'\u00e9tape 2 ; elle n'est pas retenue.", "decide")
		}
		if (isTRUE(obj$contradiction)) {
			.bay_tell(st, "Caution added: the two methods answer firmly and differently; the bilan will not decide.",
				"Prudence ajout\u00e9e : les deux m\u00e9thodes r\u00e9pondent nettement et diff\u00e9remment ; le bilan ne tranchera pas.", "decide")
		}
	}
	sens <- obj$sensitivity
	if (is.null(sens)) {
		.bay_tell(st, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).", "decide")
	} else {
		.bay_tell(st, paste0("Sensitivity: \"", h$label_en[1L], "\" ", .bay_range_txt(sens$range_lead)$en, "."),
			paste0("Sensibilit\u00e9 : \u00ab ", h$label_fr[1L], " \u00bb ", .bay_range_txt(sens$range_lead)$fr, "."))
		if (identical(sens$verdict, "changes")) {
			.bay_tell(st, "Caution added: the leading hypothesis can reverse; the bilan will not assert it.",
				"Prudence ajout\u00e9e : l'hypoth\u00e8se en t\u00eate peut s'inverser ; le bilan ne l'affirmera pas.", "decide")
		} else if (identical(sens$verdict, "crosses")) {
			.bay_tell(st, "The order does not change; only the degree of support can.",
				"L'ordre ne change pas ; seul le degr\u00e9 de soutien peut changer.", "decide")
		} else {
			.bay_tell(st, "The order does not change.", "L'ordre ne change pas.", "decide")
		}
		k1 <- .bay_num_txt(min(sens$kappa)); k2 <- .bay_num_txt(max(sens$kappa))
		.bay_tell(st, paste0("kappa goes from ", k1$en, " to ", k2$en, ", the starting chances from 1 against 2 to 2 against 1."),
			paste0("kappa va de ", k1$fr, " \u00e0 ", k2$fr, ", les chances de d\u00e9part de 1 contre 2 \u00e0 2 contre 1."), "note")
	}
	.bay_tell(st, "the Bayes factor says how many times the data are more probable with a link than without.",
		"le facteur de Bayes dit combien de fois les donn\u00e9es sont plus probables avec un lien que sans.", "note_first")
}

.corbay_story_estimate <- function(st, obj) {
	est <- obj$estimate
	mw <- .corbay_method_words(obj$method)
	.bay_tell(st, "Strength of the link.", "Force du lien.", "step",
		fun = paste0(if (obj$method == "kendall") "cor.bay() {KefiR}" else "posterior() {BayesFactor}", ", hdi() {bayestestR}"))
	.bay_step_result(st, paste0(mw$par, " = ", .corbay_est_txt(est, ".")), paste0(mw$par, " = ", .corbay_est_txt(est, ",")))
	ob <- obj$observed
	obs <- if (obj$method == "kendall") ob[["tau"]] else if (obj$method == "spearman") ob[["rho_s"]] else ob[["r"]]
	h0 <- isTRUE(obj$hypotheses$null[1L])
	.bay_tell(st, paste0(if (h0) "If a link exists, estimated " else "Estimated ", mw$par, ": ", .corbay_est_txt(est, "."), "; observed ", .corbay_num(obs, 2L, "."), "."),
		paste0(if (h0) "Si un lien existe, " else "", mw$par, " estim\u00e9 : ", .corbay_est_txt(est, ","), " ; observ\u00e9 ", .corbay_num(obs, 2L, ","), "."))
	if (obj$method == "kendall") {
		.bay_tell(st, paste0("For comparison, Pearson's r observed: ", .corbay_num(ob[["r"]], 2L, "."), " (tau is usually about 2/3 of rho)."),
			paste0("Pour comparer, r de Pearson observ\u00e9 : ", .corbay_num(ob[["r"]], 2L, ","), " (tau vaut en g\u00e9n\u00e9ral les 2/3 de rho environ)."))
	}
	pc <- format(round(100 * est$ci), trim = TRUE)
	.bay_tell(st, paste0("[; ] = ", pc, " % credible interval (", pc, " % chances to find ", mw$par, " in it)."),
		paste0("[ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 ", pc, " % (", pc, " % de chances d'y trouver ", mw$par, ")."), "note")
	.bay_tell(st, "the estimate is drawn towards 0 by the prior, the more so with few pairs.",
		"l'estimation est ramen\u00e9e vers 0 par l'a priori, d'autant plus que les paires sont peu nombreuses.", "note")
}

# ---------------------------------------------------------------------------
# Bilan
# ---------------------------------------------------------------------------

.corbay_verdict <- function(obj, li) {
	h <- obj$hypotheses
	p <- h$posterior[1L]
	if (isTRUE(obj$contradiction)) {
		return(list(en = "Verdict: undecided, the answer depends on the method (Pearson or ranks).",
			fr = "Verdict : ind\u00e9cis, la r\u00e9ponse d\u00e9pend de la m\u00e9thode (Pearson ou rangs)."))
	}
	if (abs(li$lb) < log(3)) p <- min(p, 0.59)
	open <- identical(obj$sensitivity$verdict, "changes")
	if (p < 0.6) {
		return(list(en = "Verdict: undecided. The data allow no conclusion, neither a link nor its absence.",
			fr = "Verdict : ind\u00e9cis. Les donn\u00e9es ne permettent de conclure ni \u00e0 un lien, ni \u00e0 son absence."))
	}
	if (isTRUE(h$null[1L])) {
		if (p >= 0.95 && !open) return(list(en = "Verdict: very probably no link of the size expected. This is not a proof of independence.",
			fr = "Verdict : tr\u00e8s probablement pas de lien de l'ordre attendu. Ce n'est pas une preuve d'ind\u00e9pendance."))
		return(list(en = "Verdict: rather no link, without concluding. This does not prove that the variables are independent.",
			fr = "Verdict : plut\u00f4t pas de lien, sans conclure. Cela ne prouve pas que les variables sont ind\u00e9pendantes."))
	}
	s_en <- if (h$key[1L] == "pos") "positive" else "negative"
	s_fr <- if (h$key[1L] == "pos") "positif" else "n\u00e9gatif"
	if (p >= 0.95 && !open) return(list(en = paste0("Verdict: a ", s_en, " link very probably exists. This says neither that it is strong, nor that it is a cause."),
		fr = paste0("Verdict : un lien ", s_fr, " existe tr\u00e8s probablement. Cela ne dit ni qu'il est fort, ni qu'il est une cause.")))
	if (p >= 0.75 && !open) return(list(en = paste0("Verdict: a ", s_en, " link probably exists. This says neither that it is strong, nor that it is a cause."),
		fr = paste0("Verdict : un lien ", s_fr, " existe probablement. Cela ne dit ni qu'il est fort, ni qu'il est une cause.")))
	list(en = paste0("Verdict: rather a ", s_en, " link, without concluding."), fr = paste0("Verdict : plut\u00f4t un lien ", s_fr, ", sans conclure."))
}

.corbay_lead_info <- function(obj) {
	h <- obj$hypotheses
	lb <- if (nrow(h) > 1L) log(h$posterior[1L] / h$posterior[2L]) - log(h$prior[1L] / h$prior[2L]) else NA_real_
	sup <- .bay_support_bf(h$posterior[1L], lb)
	if (isTRUE(obj$contradiction)) sup <- .bay_support(min(h$posterior[1L], 0.6))
	open <- identical(obj$sensitivity$verdict, "changes")
	if (open) sup <- list(en = "not settled", fr = "non tranch\u00e9", level = 1L)
	list(h = h, lb = lb, sup = sup, open = open)
}

.corbay_print_bilan <- function(obj, detail = FALSE, narrated = FALSE) {
	emit <- .bay_block_emitter()
	li <- .corbay_lead_info(obj)
	h <- li$h
	mw <- .corbay_method_words(obj$method)
	undecided <- isTRUE(obj$contradiction) || li$open || li$sup$level <= 1L || h$posterior[1L] < 0.6
	lead_en <- if (undecided) paste0("undecided (", h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), ")") else
		paste0(h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), " (", li$sup$en, ")")
	lead_fr <- if (undecided) paste0("ind\u00e9cis (", h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), ")") else
		paste0(h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), " (", li$sup$fr, ")")
	if (isTRUE(narrated)) {
		cat("\n")
		k <- obj$story$k + 1L
		emit("title", paste0(k, ") Bilan - ", lead_en, "."), paste0(k, ") Bilan - ", lead_fr, "."))
	} else {
		emit("title", paste0("Bilan of cor.bay (", obj$label, ", N = ", obj$n, ") - ", lead_en, "."),
			paste0("Bilan de cor.bay (", obj$label, ", N = ", obj$n, ") - ", lead_fr, "."))
	}
	v <- .corbay_verdict(obj, li)
	emit("fact", sub("^([^.]*\\.).*$", "\\1", v$en), sub("^([^.]*\\.).*$", "\\1", v$fr))
	w_en <- max(nchar(h$label_en)); w_fr <- max(nchar(h$label_fr))
	fit <- max(w_en, w_fr) + 32L <= .bay_width() - 12L
	pad <- function(x, w) if (fit) paste0(x, strrep("\u00a0", w - nchar(x))) else x
	emit("fact", "Probability of each hypothesis:", "Probabilit\u00e9 de chaque hypoth\u00e8se :")
	for (i in seq_len(nrow(h))) {
		emit("sub", paste0(i, ". ", pad(h$label_en[i], w_en), " : ", .bay_pct(h$posterior[i], "."), if (i == 1L) paste0(" (", li$sup$en, ")") else ""),
			paste0(i, ". ", pad(h$label_fr[i], w_fr), " : ", .bay_pct(h$posterior[i], ","), if (i == 1L) paste0(" (", li$sup$fr, ")") else ""))
	}
	sens <- obj$sensitivity
	if (!is.null(sens)) {
		rg <- sens$range_lead
		if (li$open) {
			emit("decide", paste0("Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order can reverse, not settled."),
				paste0("Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre peut s'inverser, non tranch\u00e9."))
		} else if (identical(.bay_pct(rg[1L], "."), .bay_pct(rg[2L], ".")) || rg[1L] > 0.999) {
			emit("fact", paste0("The answer ", .bay_range_txt(rg)$en, "."), paste0("La r\u00e9ponse ", .bay_range_txt(rg)$fr, "."))
		} else {
			emit("fact", paste0("Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order does not change."),
				paste0("Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre ne change pas."))
		}
	} else {
		emit("note", "sensitivity to the priors not computed (sensitivity = FALSE).", "sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}
	bf_en <- .bay_bf_txt(abs(li$lb), "."); bf_fr <- .bay_bf_txt(abs(li$lb), ",")
	emit("note", paste0("the data are ", bf_en, " times more probable under hypothesis 1 than under 2 (Bayes factor; 1 to 3 weak, 3 to 20 moderate, 20 to 150 strong, beyond very strong)."),
		paste0("les donn\u00e9es sont ", bf_fr, " fois plus probables sous l'hypoth\u00e8se 1 que sous la 2 (facteur de Bayes ; 1 \u00e0 3 faible, 3 \u00e0 20 mod\u00e9r\u00e9, 20 \u00e0 150 fort, au-del\u00e0 tr\u00e8s fort)."))
	emit("block", "", "")
	est <- obj$estimate
	h0 <- isTRUE(h$null[1L])
	if (!is.null(est)) {
		emit("fact", paste0(if (h0) "If a link exists, estimated " else "Estimated ", mw$par, ": ", .corbay_est_txt(est, "."), if (narrated) " (step 6)." else "."),
			paste0(if (h0) "Si un lien existe, " else "", mw$par, " estim\u00e9 : ", .corbay_est_txt(est, ","), if (narrated) " (\u00e9tape 6)." else "."))
	}
	if (!narrated) emit("fact", paste0("Method: ", mw$en, if (obj$auto) " (chosen by the rule of method = \"auto\")" else "", "."),
		paste0("M\u00e9thode : ", mw$fr, if (obj$auto) " (choisie par la r\u00e8gle de method = \u00ab auto \u00bb)" else "", "."))
	emit("block", "", "")
	ob <- obj$observed
	emit("fact", paste0("Observed: r = ", .corbay_num(ob[["r"]], 2L, "."), " (Pearson), tau = ", .corbay_num(ob[["tau"]], 2L, "."), " (Kendall), rho = ",
			.corbay_num(ob[["rho_s"]], 2L, "."), " (Spearman), N = ", obj$n, " pairs."),
		paste0("Observ\u00e9 : r = ", .corbay_num(ob[["r"]], 2L, ","), " (Pearson), tau = ", .corbay_num(ob[["tau"]], 2L, ","), " (Kendall), rho = ",
			.corbay_num(ob[["rho_s"]], 2L, ","), " (Spearman), N = ", obj$n, " paires."))
	if (length(obj$notes$en) && !isTRUE(narrated)) for (i in seq_along(obj$notes$en)) emit("fact", obj$notes$en[i], obj$notes$fr[i])
	est_en <- if (!is.null(est)) paste0("; ", mw$par, " = ", .corbay_est_txt(est, "."), ", ", format(round(100 * est$ci)), " % credible interval") else ""
	est_fr <- if (!is.null(est)) paste0(" ; ", mw$par, " = ", .corbay_est_txt(est, ","), ", intervalle de cr\u00e9dibilit\u00e9 \u00e0 ", format(round(100 * est$ci)), " %") else ""
	if (h0) { est_en <- ""; est_fr <- "" }
	h1 <- which(!h$null)[1L]
	p_link <- sum(h$posterior[!h$null])
	if (undecided) {
		emit("note", paste0("for a report: \"The data do not settle it (probability of a link ", .bay_pct(p_link, "."), ", Bayes factor ", .bay_bf_txt(obj$bf$log_bf10, "."),
				"); cor.bay {KefiR}, ", mw$cite_en, ".\""),
			paste0("pour un compte rendu : \u00ab Les donn\u00e9es ne permettent pas de conclure (probabilit\u00e9 d'un lien ", .bay_pct(p_link, ","), ", facteur de Bayes ",
				.bay_bf_txt(obj$bf$log_bf10, ","), ") ; cor.bay {KefiR}, ", mw$cite_fr, ". \u00bb"))
	} else {
		emit("note", paste0("for a report: \"", h$label_en[1L], " (posterior probability ", .bay_pct(h$posterior[1L], "."), "), Bayes factor ", bf_en, est_en,
				" (cor.bay {KefiR}, ", mw$cite_en, ").\""),
			paste0("pour un compte rendu : \u00ab ", h$label_fr[1L], " (probabilit\u00e9 a posteriori ", .bay_pct(h$posterior[1L], ","), "), facteur de Bayes ", bf_fr, est_fr,
				" (cor.bay {KefiR}, ", mw$cite_fr, "). \u00bb"))
	}
	emit("note", paste0("probabilities of hypotheses, not p-values; a link is not a cause; the method with print(result, detail = TRUE)."),
		paste0("probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values ; un lien n'est pas une cause ; la m\u00e9thode avec print(resultat, detail = TRUE)."))
	emit("flush", "", "")
	if (isTRUE(detail)) .corbay_print_detail(obj)
	invisible(NULL)
}

.corbay_print_detail <- function(obj) {
	say <- .bay_detail_head(obj)
	mw <- .corbay_method_words(obj$method)
	kt <- .bay_num_txt(obj$priors$kappa, 3)
	say(paste0("Hypotheses: no link (rho = 0), prior ", .bay_pct(1 - obj$priors$prior_h1, "."), "; a link, prior ", .bay_pct(obj$priors$prior_h1, "."),
			if (obj$alternative == "two.sided") ", split equally between its two signs." else ", on the sign tested only."),
		paste0("Hypoth\u00e8ses : aucun lien (rho = 0), a priori ", .bay_pct(1 - obj$priors$prior_h1, ","), " ; un lien, a priori ", .bay_pct(obj$priors$prior_h1, ","),
			if (obj$alternative == "two.sided") ", partag\u00e9 \u00e9galement entre ses deux signes." else ", sur le seul signe test\u00e9."))
	if (obj$method == "kendall") {
		say("For Kendall: T* = standardised tau-b, normal with mean 3/2 tau sqrt(N) and standard deviation 1 (van Doorn et al. 2018, eq. 6 and 11); Bayes factor by numerical integration (eq. 13); normal approximation checked from N = 10 (Ferguson et al. 2000).",
			"Pour Kendall : T* = tau-b standardis\u00e9, normal de moyenne 3/2 tau racine(N) et d'\u00e9cart-type 1 (van Doorn et al. 2018, \u00e9q. 6 et 11) ; facteur de Bayes par int\u00e9gration num\u00e9rique (\u00e9q. 13) ; approximation normale v\u00e9rifi\u00e9e d\u00e8s N = 10 (Ferguson et al. 2000).")
	} else {
		say(paste0(if (obj$method == "spearman") "For Spearman: ranks, then " else "", "correlationBF() {BayesFactor}, bivariate normal model, Jeffreys's test with a stretched beta prior (Ly et al. 2016); kappa = ", kt$en, "."),
			paste0(if (obj$method == "spearman") "Pour Spearman : rangs, puis " else "", "correlationBF() {BayesFactor}, mod\u00e8le normal bivari\u00e9, test de Jeffreys avec a priori b\u00eata \u00e9tir\u00e9 (Ly et al. 2016) ; kappa = ", kt$fr, "."))
	}
	if (!is.null(obj$sensitivity)) {
		g <- obj$sensitivity$grid
		say(paste0("Sensitivity grid (probability of the leading hypothesis): kappa ", paste(rownames(g), collapse = ", "), " by odds 1:2, 1:1, 2:1: ",
				paste(vapply(g, .bay_pct, character(1), dec = "."), collapse = ", "), "."),
			paste0("Grille de sensibilit\u00e9 (probabilit\u00e9 de l'hypoth\u00e8se en t\u00eate) : kappa ", gsub(".", ",", paste(rownames(g), collapse = " ; "), fixed = TRUE), " par cotes 1:2, 1:1, 2:1 : ",
				paste(vapply(g, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))
	}
	pr <- .bay_prec_txt(obj$bf$error_pct)
	say(paste0("Numerical error on the Bayes factor: ", pr$en, "."),
		paste0("Erreur num\u00e9rique sur le facteur de Bayes : ", pr$fr, "."))
	invisible(NULL)
}
