# sys_bay_manova_fact.R
# Several responses and several factors between groups (m.bay, DECISIONS, D23). Same
# stacking as with one factor (sys_bay_manova.R): standardised responses, y ~ terms + response
# + individual (random), lmBF. Each term of the factors enters in one of two forms: the same
# gap on every response (t), or gaps that change with the response (t + t:response). The
# models are the hereditary combinations of the terms (sys_bay_factorial.R), each term in
# either form; prior HOP on the terms, as with one response, then 1/2 per form. Nothing is
# exported.

# Expanded models: hereditary model m, and the form of each of its terms (0 same gap, 1 gaps
# that change with the response).
.bay_mf_models <- function(sp) {
	out <- list()
	for (j in seq_along(sp$models)) {
		m <- sp$models[[j]]
		if (!length(m)) {
			out[[length(out) + 1L]] <- list(base = j, terms = integer(0), forms = integer(0))
			next
		}
		grid <- as.matrix(expand.grid(rep(list(0:1), length(m))))
		for (r in seq_len(nrow(grid))) out[[length(out) + 1L]] <- list(base = j, terms = m, forms = as.integer(grid[r, ]))
	}
	out
}

# Log BF of every expanded model against "response + individual", for each r.
.bay_mf_lbf <- function(Y, factors, sp, em, rs, iterations, seeds) {
	safe <- paste0("f", seq_along(factors))
	names(safe) <- names(factors)
	n <- nrow(Y)
	d <- ncol(Y)
	long <- data.frame(y = as.vector(scale(Y)), vd = factor(rep(colnames(Y), each = n), levels = colnames(Y)),
		id = factor(rep(seq_len(n), d)))
	for (v in names(factors)) long[[safe[[v]]]] <- rep(factors[[v]], d)
	term_safe <- vapply(strsplit(sp$keys, ":", fixed = TRUE), function(f) paste(safe[f], collapse = ":"), character(1))
	vapply(seq_along(rs), function(i) {
		.bay_with_seed(seeds[i], {
			one <- function(rhs) BayesFactor::extractBF(BayesFactor::lmBF(stats::as.formula(paste("y ~", rhs)), data = long,
				whichRandom = "id", rscaleFixed = rs[i], iterations = iterations, progress = FALSE), logbf = TRUE)$bf
			b0 <- one("vd + id")
			vapply(em, function(x) {
				if (!length(x$terms)) return(0)
				tt <- c(term_safe[x$terms], if (any(x$forms == 1L)) paste0(term_safe[x$terms][x$forms == 1L], ":vd"))
				one(paste(c(tt, "vd", "id"), collapse = " + ")) - b0
			}, numeric(1))
		})
	}, numeric(length(em)))
}

# Answers per term from the expanded models: P(term), its Bayes factor, P(gaps that change
# with the response | term).
.bay_mf_answers <- function(sp, em, lb, lp_base) {
	lp <- vapply(em, function(x) lp_base[x$base] + length(x$terms) * log(1 / 2), numeric(1))
	post <- lb + lp
	post <- exp(post - .bay_log_sum_exp(post))
	pri <- exp(lp)
	eff <- do.call(rbind, lapply(seq_along(sp$keys), function(t) {
		has <- vapply(em, function(x) t %in% x$terms, logical(1))
		p1 <- sum(post[has]); q1 <- sum(pri[has])
		chg <- vapply(em, function(x) any(x$terms == t & x$forms == 1L), logical(1))
		data.frame(term = sp$terms[t], p_h1 = p1, p_h0 = 1 - p1, prior_h1 = q1,
			log_bf = log(p1 / (1 - p1)) - log(q1 / (1 - q1)), p_changes = if (p1 > 0) sum(post[chg]) / p1 else NA_real_,
			stringsAsFactors = FALSE)
	}))
	eff$bf10 <- exp(eff$log_bf)
	list(effects = eff, p_null = post[1L], prior_null = pri[1L], posterior = post, prior = pri)
}

.bay_m_manova_fact <- function(input, call_expr, get_notes, note, rscale, sensitivity, draws, seed, detail, verbose,
                               return, story) {
	Y <- as.matrix(input$response)
	if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
	colnames(Y) <- make.unique(colnames(Y))
	factors <- input$factors
	for (v in names(factors)) factors[[v]] <- droplevels(as.factor(factors[[v]]))
	terms <- gsub("`", "", input$terms, fixed = TRUE)
	sp <- .bay_fact_space(terms)
	em <- .bay_mf_models(sp)
	d <- ncol(Y)
	N <- nrow(Y)
	resp_txt <- paste(colnames(Y), collapse = ", ")
	fac_txt <- paste(names(factors), collapse = ", ")
	cells <- interaction(factors, drop = FALSE)
	if (any(table(cells) == 0L) && any(lengths(strsplit(terms, ":", fixed = TRUE)) > 1L)) {
		.exit("Empty cell: the interaction cannot be estimated, remove it from the formula.",
			"Case vide : l'interaction ne peut pas \u00eatre estim\u00e9e, la retirer de la formule.")
	}
	if (any(apply(Y, 2L, stats::sd) == 0)) {
		.exit("A response is constant: it cannot be standardised.", "Une r\u00e9ponse est constante : elle ne peut pas \u00eatre standardis\u00e9e.")
	}

	.bay_tell(story, paste0("Plan recognised: ", d, " responses (", resp_txt, ") compared together according to ", fac_txt,
			" (", paste0(names(factors), " ", vapply(factors, nlevels, integer(1)), " levels", collapse = ", "), "); N = ", N, " individuals."),
		paste0("Plan reconnu : ", d, " r\u00e9ponses (", resp_txt, ") compar\u00e9es ensemble selon ", fac_txt,
			" (", paste0(names(factors), " ", vapply(factors, nlevels, integer(1)), " niveaux", collapse = ", "), ") ; N = ", N, " individus."), "step")
	.bay_tell(story, paste0("Terms: ", paste(terms, collapse = ", "), "; for each, do the groups differ on at least one response?"),
		paste0("Termes : ", paste(terms, collapse = ", "), " ; pour chacun, les groupes diff\u00e8rent-ils sur au moins une r\u00e9ponse ?"), "decide")

	.bay_tell(story, "Checks on the data", "Contr\u00f4les des donn\u00e9es", "step", fun = "quantile() {stats}")
	cn <- table(interaction(factors, drop = TRUE))
	.bay_tell(story, paste0("Sizes: ", min(cn), if (max(cn) > min(cn)) paste0(" to ", max(cn)), " per cell."),
		paste0("Effectifs : ", min(cn), if (max(cn) > min(cn)) paste0(" \u00e0 ", max(cn)), " par case."))
	cl <- interaction(factors, drop = TRUE, sep = " x ", lex.order = TRUE)
	far <- vapply(seq_len(d), function(j) {
		r <- Y[, j] - stats::ave(Y[, j], cl)
		q <- stats::quantile(r, c(0.25, 0.75))
		sum(r < q[1L] - 3 * diff(q) | r > q[2L] + 3 * diff(q))
	}, numeric(1))
	.bay_step_result(story, if (sum(far)) paste0(sum(far), " far values kept") else "nothing to remove",
		if (sum(far)) paste0(sum(far), " valeurs \u00e9loign\u00e9es gard\u00e9es") else "rien \u00e0 retirer")
	.bay_tell(story, paste0("Values beyond 3 interquartile ranges from their cell: ", paste0(colnames(Y), " ", far, collapse = ", "), "."),
		paste0("Valeurs au-del\u00e0 de 3 \u00e9carts interquartiles de leur case : ", paste0(colnames(Y), " ", far, collapse = ", "), "."))
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	.bay_tell(story, "Shape of the noise", "Forme du bruit", "step")
	.bay_step_result(story, "normal on several responses, correlations taken equal", "normale \u00e0 plusieurs r\u00e9ponses, corr\u00e9lations prises \u00e9gales")
	E <- Y - apply(Y, 2L, function(v) stats::ave(v, cl))
	R <- stats::cor(E)
	rr <- R[upper.tri(R)]
	f2 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
	.bay_tell(story, paste0("Each response is standardised; correlations between responses within cells: ",
			if (length(rr) == 1L) f2(rr, ".") else paste0("from ", f2(min(rr), "."), " to ", f2(max(rr), ".")), "."),
		paste0("Chaque r\u00e9ponse est standardis\u00e9e ; corr\u00e9lations entre r\u00e9ponses au sein des cases : ",
			if (length(rr) == 1L) f2(rr, ",") else paste0("de ", f2(min(rr), ","), " \u00e0 ", f2(max(rr), ",")), "."))
	uneq <- length(rr) > 1L && diff(range(rr)) > 0.5
	.bay_tell(story, if (uneq) "The model takes these correlations equal; here they are far from it: it may then miss a difference, without inventing one; a \"no\" is to be read with caution." else "The model takes these correlations equal, which they roughly are.",
		if (uneq) "Le mod\u00e8le prend ces corr\u00e9lations \u00e9gales ; ici elles en sont loin : il peut alors manquer une diff\u00e9rence, sans en inventer ; un \u00ab non \u00bb est \u00e0 lire avec r\u00e9serve." else "Le mod\u00e8le prend ces corr\u00e9lations \u00e9gales, ce qu'elles sont \u00e0 peu pr\u00e8s.", "decide")
	.bay_tell(story, "equal spreads in the cells are assumed, as in a classical MANOVA.",
		"des dispersions \u00e9gales dans les cases sont suppos\u00e9es, comme dans une MANOVA classique.", "note")

	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	lp_base <- .bay_fact_prior(sp, "hop")
	.bay_step_result(story, paste0("no difference at all: ", .bay_pct(exp(lp_base[1L]), ".")), paste0("aucune diff\u00e9rence du tout : ", .bay_pct(exp(lp_base[1L]), ",")))
	incl <- vapply(seq_along(terms), function(i) sum(exp(lp_base[vapply(sp$models, function(m) i %in% m, logical(1))])), numeric(1))
	.bay_tell(story, paste0("Chance of each term in advance: ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = "."), collapse = ", "),
			"; each term, if present, takes one of two forms at even odds (same gap on every response, or gaps that change with the response)."),
		paste0("Chance de chaque terme d'avance : ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = ","), collapse = ", "),
			" ; chaque terme pr\u00e9sent prend l'une de deux formes \u00e0 poids \u00e9gaux (m\u00eame \u00e9cart sur chaque r\u00e9ponse, ou \u00e9carts qui changent avec la r\u00e9ponse)."))
	r_txt <- .bay_num_txt(rscale)
	.bay_tell(story, paste0("Expected size r = ", r_txt$en, " per term."), paste0("Taille attendue r = ", r_txt$fr, " par terme."))
	.bay_tell(story, "Prior on the terms: HOP (Taylor-Rodriguez et al. 2016), strong heredity, as with one response; 1/2 per form.",
		"A priori sur les termes : HOP (Taylor-Rodriguez et al. 2016), h\u00e9r\u00e9dit\u00e9 forte, comme avec une r\u00e9ponse ; 1/2 par forme.", "detail")

	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = "lmBF() {BayesFactor}")
	rs <- if (isTRUE(sensitivity)) c(rscale, rscale / 2, rscale * 2) else rscale
	n_calls <- length(em) * length(rs)
	# A waiting line when the computation is long (more than about 30 s).
	if (isTRUE(verbose) && n_calls > 25L) {
		message(.msg(paste0("Computing ", length(em), " models", if (length(rs) > 1L) " for 3 values of r" else "", ", about ", ceiling(n_calls * 1.2 / 10) * 10, " s..."),
			paste0("Calcul de ", length(em), " mod\u00e8les", if (length(rs) > 1L) " pour 3 valeurs de r" else "", ", environ ", ceiling(n_calls * 1.2 / 10) * 10, " s...")))
	}
	iterations <- max(10000L, as.integer(draws))
	lb <- .bay_mf_lbf(Y, factors, sp, em, rs, iterations, seed + c(0L, 1L, 2L)[seq_along(rs)])
	lb <- matrix(lb, ncol = length(rs))
	an <- .bay_mf_answers(sp, em, lb[, 1L], lp_base)
	eff <- an$effects
	.bay_tell(story, paste0(length(em), " models compared (combinations of terms, each in two forms)."),
		paste0(length(em), " mod\u00e8les compar\u00e9s (combinaisons de termes, chacun sous deux formes)."))
	sens <- NULL
	if (isTRUE(sensitivity)) {
		priors <- list(hop = lp_base, hip = .bay_fact_prior(sp, "hip"), unif = .bay_fact_prior(sp, "unif"))
		grid <- expand.grid(r = seq_along(rs), prior = names(priors), stringsAsFactors = FALSE)
		pt <- vapply(seq_len(nrow(grid)), function(i) .bay_mf_answers(sp, em, lb[, grid$r[i]], priors[[grid$prior[i]]])$effects$p_h1, numeric(nrow(eff)))
		pt <- matrix(pt, nrow = nrow(eff))
		pa <- vapply(seq_len(nrow(grid)), function(i) 1 - .bay_mf_answers(sp, em, lb[, grid$r[i]], priors[[grid$prior[i]]])$p_null, numeric(1))
		rng_terms <- t(apply(pt, 1L, range))
		verdict_terms <- vapply(seq_len(nrow(eff)), function(i) {
			lead <- eff$p_h1[i] >= 0.5
			r <- if (lead) rng_terms[i, ] else 1 - rev(rng_terms[i, ])
			.bay_sensitivity_verdict(max(eff$p_h1[i], eff$p_h0[i]), r)
		}, character(1))
		lead_any <- (1 - an$p_null) >= an$p_null
		rng_any <- range(if (lead_any) pa else 1 - pa)
		sens <- list(r = sort(rs), priors = names(priors), p_terms = pt, p_any = pa, range_terms = rng_terms,
			verdict_terms = verdict_terms, range_lead = rng_any, verdict = .bay_sensitivity_verdict(max(an$p_null, 1 - an$p_null), rng_any))
		moved <- which(verdict_terms == "changes")
		.bay_step_result(story, if (length(moved)) paste0(paste(eff$term[moved], collapse = ", "), " change side with the prior") else "no term changes side",
			if (length(moved)) paste0(paste(eff$term[moved], collapse = ", "), " changent de c\u00f4t\u00e9 selon l'a priori") else "aucun terme ne change de c\u00f4t\u00e9")
		for (i in moved) {
			.bay_tell(story, paste0(eff$term[i], " changes side: ", .bay_pct(rng_terms[i, 1L], "."), " to ", .bay_pct(rng_terms[i, 2L], "."), "."),
				paste0(eff$term[i], " change de c\u00f4t\u00e9 : ", .bay_pct(rng_terms[i, 1L], ","), " \u00e0 ", .bay_pct(rng_terms[i, 2L], ","), "."))
		}
		.bay_tell(story, "r goes from r/2 to 2r, with three priors on the models.", "r va de r/2 \u00e0 2r, avec trois lois sur les mod\u00e8les.", "note")
	} else {
		.bay_step_result(story, "sensitivity not computed", "sensibilit\u00e9 non calcul\u00e9e")
	}

	.bay_tell(story, "Which responses differ", "Quelles r\u00e9ponses diff\u00e8rent", "step")
	.bay_step_result(story, "by term, the form of the difference", "par terme, la forme de la diff\u00e9rence")
	for (i in seq_len(nrow(eff))) {
		.bay_tell(story, paste0(eff$term[i], ": if it acts, gaps that change from one response to another at ", .bay_pct(eff$p_changes[i], "."), "."),
			paste0(eff$term[i], " : s'il agit, \u00e9carts qui changent d'une r\u00e9ponse \u00e0 l'autre \u00e0 ", .bay_pct(eff$p_changes[i], ","), "."))
	}
	.bay_tell(story, "to see which response carries a term, run m.bay() on that response alone.",
		"pour voir quelle r\u00e9ponse porte un terme, lancer m.bay() sur cette r\u00e9ponse seule.", "note")

	p_any <- 1 - an$p_null
	hyp <- data.frame(
		label_en = c(paste0("At least one of the responses (", resp_txt, ") differs according to at least one factor"), paste0("No difference in ", resp_txt, " according to ", fac_txt)),
		label_fr = c(paste0("Au moins une des r\u00e9ponses (", resp_txt, ") diff\u00e8re selon au moins un facteur"), paste0("Aucune diff\u00e9rence de ", resp_txt, " selon ", fac_txt)),
		posterior = c(p_any, an$p_null), prior = c(1 - an$prior_null, an$prior_null), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(hypotheses = hyp, effects = eff[, c("term", "bf10", "log_bf", "p_h1", "p_h0", "prior_h1", "p_changes")],
		laws = NULL,
		models = data.frame(model = vapply(em, function(x) {
				if (!length(x$terms)) return("(null)")
				paste(paste0(sp$terms[x$terms], ifelse(x$forms == 1L, " (x response)", "")), collapse = " + ")
			}, character(1)), log_bf = lb[, 1L], posterior = an$posterior, prior = an$prior, stringsAsFactors = FALSE),
		questions = list(en = paste0("Do the groups differ on at least one of ", resp_txt, ", term by term?"),
			fr = paste0("Les groupes diff\u00e8rent-ils sur au moins une de ", resp_txt, ", terme par terme ?")),
		correlations = R, unequal_correlations = uneq, sensitivity = sens, posthoc = NULL,
		descriptive_stats = data.frame(cell = rep(levels(cl), d), response = rep(colnames(Y), each = nlevels(cl)),
			mean = as.vector(apply(Y, 2L, function(v) tapply(v, cl, mean))), stringsAsFactors = FALSE),
		priors = list(rscale = rscale, model_prior = "HOP(1,1), 1/2 per form"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = p_any,
		method = "Several responses and factors: standardised, stacked, lmBF with a random individual effect (BayesFactor)",
		plan = "manova", n_factors = length(factors), label = input$label, n = N, notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_manova(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_any)
	if (isTRUE(verbose)) invisible(out) else out
}
