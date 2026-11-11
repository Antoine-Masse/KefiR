# sys_bay_discrete_fact.R
# Discrete responses with several factors between groups (m.bay, DECISIONS, D25): yes/no
# (logistic regression) and counts (Poisson and negative binomial, averaged). Same logic as
# with one factor (D17): one glm per hereditary model of sys_bay_factorial.R, marginal
# likelihood by -BIC/2 (unit information; Schwarz 1978, Raftery 1995), HOP prior on the
# models, an answer per term. contingencyTableBF is not used: it handles a two-way table only
# (one factor against the response), not the three-way table of two factors. Nothing is
# exported.

# -BIC/2 of every hereditary model, one column per law.
.bay_df_lml <- function(y, factors, sp, kind) {
	safe <- paste0("f", seq_along(factors))
	names(safe) <- names(factors)
	d <- data.frame(y = y)
	for (v in names(factors)) d[[safe[[v]]]] <- factors[[v]]
	term_safe <- vapply(strsplit(sp$keys, ":", fixed = TRUE), function(f) paste(safe[f], collapse = ":"), character(1))
	one <- function(m, law) {
		fo <- stats::as.formula(paste("y ~", if (length(m)) paste(term_safe[m], collapse = " + ") else "1"))
		fit <- tryCatch(switch(law,
			binomial = .bay_quiet(stats::glm(fo, data = d, family = stats::binomial())),
			poisson = .bay_quiet(stats::glm(fo, data = d, family = stats::poisson())),
			negbin = .bay_quiet(MASS::glm.nb(fo, data = d)),
			polr = .bay_quiet(MASS::polr(fo, data = d))), error = function(e) NULL)
		if (is.null(fit)) return(NA_real_)
		if (identical(law, "negbin") && (!is.finite(fit$theta) || fit$theta > 1e6)) return(NA_real_)
		-stats::BIC(fit) / 2
	}
	laws <- if (kind == "binary") "binomial" else if (kind == "ordinal") "polr" else c("poisson", "negbin")
	out <- vapply(laws, function(l) vapply(sp$models, one, numeric(1), law = l), numeric(length(sp$models)))
	matrix(out, ncol = length(laws), dimnames = list(NULL, laws))
}

.bay_m_discrete_fact <- function(input, y, kind, call_expr, get_notes, note, sensitivity, detail, verbose, return, story, posthoc) {
	.bay_cw$n <- 0L
	factors <- input$factors
	for (v in names(factors)) factors[[v]] <- droplevels(as.factor(factors[[v]]))
	terms <- gsub("`", "", input$terms, fixed = TRUE)
	sp <- .bay_fact_space(terms)
	resp <- input$response_name
	fac_txt <- paste(names(factors), collapse = ", ")
	cl <- interaction(factors, sep = " x ", drop = TRUE, lex.order = TRUE)
	cn <- table(cl)
	bin <- identical(kind, "binary")
	ord <- identical(kind, "ordinal")
	if (bin) {
		yc <- if (is.factor(y)) droplevels(y) else factor(y)
		yes <- levels(yc)[2L]
		y01 <- as.integer(yc == yes)
		stat <- tapply(y01, cl, mean)
	} else if (ord) {
		# Ordered response (D26): the values keep their order only; mean rank of the level per cell.
		y01 <- if (is.ordered(y)) droplevels(y) else factor(y, levels = sort(unique(y)), ordered = TRUE)
		stat <- tapply(as.integer(y01), cl, mean)
	} else {
		y01 <- round(y)
		stat <- tapply(y01, cl, mean)
	}

	# Step 2 (title told by m.bay).
	.bay_tell(story, paste0("Sizes: ", min(cn), if (max(cn) > min(cn)) paste0(" to ", max(cn)), " per cell."),
		paste0("Effectifs : ", min(cn), if (max(cn) > min(cn)) paste0(" \u00e0 ", max(cn)), " par case."))
	f2 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
	if (bin) {
		.bay_tell(story, paste0("Share of \"", yes, "\" per cell: ", paste0(names(stat), " ", vapply(stat, .bay_pct, character(1), dec = "."), collapse = ", "), "."),
			paste0("Part de \u00ab ", yes, " \u00bb par case : ", paste0(names(stat), " ", vapply(stat, .bay_pct, character(1), dec = ","), collapse = ", "), "."))
		if (any(stat %in% c(0, 1))) {
			.bay_tell(story, "A cell is all yes or all no: the logistic fit is at its limit there; read the answer with caution.",
				"Une case est tout oui ou tout non : l'ajustement logistique y est \u00e0 sa limite ; lire la r\u00e9ponse avec r\u00e9serve.", "decide")
		}
	} else if (ord) {
		.bay_tell(story, paste0(nlevels(y01), " ordered values (", paste(levels(y01), collapse = " < "), "); mean rank per cell: ", paste0(names(stat), " ", f2(stat, "."), collapse = ", "), "."),
			paste0(nlevels(y01), " valeurs ordonn\u00e9es (", paste(levels(y01), collapse = " < "), ") ; rang moyen par case : ", paste0(names(stat), " ", f2(stat, ","), collapse = ", "), "."))
	} else {
		.bay_tell(story, paste0("Mean count per cell: ", paste0(names(stat), " ", f2(stat, "."), collapse = ", "), "."),
			paste0("Comptage moyen par case : ", paste0(names(stat), " ", f2(stat, ","), collapse = ", "), "."))
	}
	.bay_step_result(story, "nothing to remove", "rien \u00e0 retirer")
	.bay_tell(story, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")

	# Step 3: law.
	.bay_tell(story, if (bin) "Law of the answers" else if (ord) "Law of the ordered values" else "Law of the counts",
		if (bin) "Loi des r\u00e9ponses" else if (ord) "Loi des valeurs ordonn\u00e9es" else "Loi des comptages", "step",
		fun = if (bin) "glm(family = binomial) {stats}" else if (ord) "polr() {MASS}" else c("glm(), glm.nb() {stats, MASS}", "glm(), glm.nb() {stats, MASS}"))
	lml <- .bay_df_lml(y01, factors, sp, kind)
	n_warn <- .bay_cw$n
	p_nb <- NA_real_
	if (ord) {
		lm_m <- lml[, 1L]
		.bay_step_result(story, "proportional odds (cumulative logits)", "cotes proportionnelles (logits cumul\u00e9s)")
		.bay_tell(story, "Each model is a proportional-odds regression: only the order of the values counts.",
			"Chaque mod\u00e8le est une r\u00e9gression \u00e0 cotes proportionnelles : seul compte l'ordre des valeurs.")
		p_np <- .bay_po_check(y01, cl)
		.bay_po_tell(story, p_np)
	} else if (bin) {
		lm_m <- lml[, 1L]
		.bay_step_result(story, "binomial (logistic regression)", "binomiale (r\u00e9gression logistique)")
		.bay_tell(story, paste0("Each model is a logistic regression of the share of \"", yes, "\"."),
			paste0("Chaque mod\u00e8le est une r\u00e9gression logistique de la part de \u00ab ", yes, " \u00bb."))
	} else {
		nb_ok <- !anyNA(lml[, "negbin"])
		lp_m <- .bay_fact_prior(sp, "hop")
		if (nb_ok) {
			lm_m <- .bay_lse2(lml[, "poisson"], lml[, "negbin"]) - log(2)
			a <- c(.bay_log_sum_exp(lp_m + lml[, "poisson"]), .bay_log_sum_exp(lp_m + lml[, "negbin"]))
			p_nb <- exp(a[2L] - .bay_log_sum_exp(a))
		} else {
			lm_m <- lml[, "poisson"]
			p_nb <- 0
		}
		lead_nb <- p_nb >= 0.5
		.bay_step_result(story, paste0(if (lead_nb) "negative binomial" else "Poisson", " ahead (", .bay_pct(max(p_nb, 1 - p_nb), "."), ")"),
			paste0(if (lead_nb) "binomiale n\u00e9gative" else "Poisson", " en t\u00eate (", .bay_pct(max(p_nb, 1 - p_nb), ","), ")"))
		.bay_tell(story, paste0("Probabilities: Poisson ", .bay_pct(1 - p_nb, "."), "; negative binomial ", .bay_pct(p_nb, "."), "; each answer is weighted over the two laws."),
			paste0("Probabilit\u00e9s : Poisson ", .bay_pct(1 - p_nb, ","), " ; binomiale n\u00e9gative ", .bay_pct(p_nb, ","), " ; chaque r\u00e9ponse est pond\u00e9r\u00e9e sur les deux lois."))
		if (!nb_ok) {
			.bay_tell(story, "The negative binomial could not be fitted for every model: Poisson alone.",
				"La binomiale n\u00e9gative n'a pas pu \u00eatre ajust\u00e9e pour tous les mod\u00e8les : Poisson seule.", "decide")
		}
	}
	if (n_warn > 0L) {
		.bay_tell(story, paste0("Caution added: unstable fits (", n_warn, if (n_warn == 1L) " warning" else " warnings", " of glm); read the answer with caution."),
			paste0("Prudence ajout\u00e9e : ajustements instables (", n_warn, if (n_warn == 1L) " avertissement" else " avertissements", " de glm) ; lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	}
	.bay_tell(story, "Marginal likelihoods by -BIC/2, the unit-information approximation (Schwarz 1978; Raftery 1995), as with one factor (DECISIONS, D17, D25).",
		"Vraisemblances marginales par -BIC/2, approximation d'information unitaire (Schwarz 1978 ; Raftery 1995), comme \u00e0 un facteur (DECISIONS, D17, D25).", "detail")

	# Step 4: priors.
	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	lp <- .bay_fact_prior(sp, "hop")
	.bay_step_result(story, paste0("no difference at all: ", .bay_pct(exp(lp[1L]), ".")), paste0("aucune diff\u00e9rence du tout : ", .bay_pct(exp(lp[1L]), ",")))
	incl <- vapply(seq_along(terms), function(i) sum(exp(lp[vapply(sp$models, function(m) i %in% m, logical(1))])), numeric(1))
	.bay_tell(story, paste0("Chance of each term in advance: ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = "."), collapse = ", "), "."),
		paste0("Chance de chaque terme d'avance : ", paste0(terms, " ", vapply(incl, .bay_pct, character(1), dec = ","), collapse = ", "), "."))
	.bay_tell(story, "prior on the effects: worth one observation (unit information), no scale to choose.",
		"a priori sur les effets : le poids d'une seule observation (information unitaire), sans \u00e9chelle \u00e0 choisir.", "note")

	# Step 5: solidity.
	.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step")
	lb <- lm_m - lm_m[1L]
	eff <- .bay_fact_terms(sp, lb, lp)
	eff <- eff[match(.bay_term_key(terms), .bay_term_key(eff$term)), , drop = FALSE]
	eff$term <- terms
	rownames(eff) <- NULL
	p_null <- .bay_fact_null(lb, lp)
	.bay_tell(story, paste0(length(sp$models), " models compared (combinations of terms)."), paste0(length(sp$models), " mod\u00e8les compar\u00e9s (combinaisons de termes)."))
	.bay_tell(story, "Precision: closed form, without draws.", "Pr\u00e9cision : forme close, sans tirage.")
	sens <- NULL
	if (isTRUE(sensitivity)) {
		priors <- list(hop = lp, hip = .bay_fact_prior(sp, "hip"), unif = .bay_fact_prior(sp, "unif"))
		pt <- vapply(priors, function(q) {
			e <- .bay_fact_terms(sp, lb, q)
			e$p_h1[match(.bay_term_key(terms), .bay_term_key(e$term))]
		}, numeric(length(terms)))
		pt <- matrix(pt, nrow = length(terms))
		pa <- vapply(priors, function(q) 1 - .bay_fact_null(lb, q), numeric(1))
		rng_terms <- t(apply(pt, 1L, range))
		verdict_terms <- vapply(seq_along(terms), function(i) {
			lead <- eff$p_h1[i] >= 0.5
			r <- if (lead) rng_terms[i, ] else 1 - rev(rng_terms[i, ])
			.bay_sensitivity_verdict(max(eff$p_h1[i], eff$p_h0[i]), r)
		}, character(1))
		lead_any <- (1 - p_null) >= p_null
		rng_any <- range(if (lead_any) pa else 1 - pa)
		sens <- list(priors = names(priors), p_terms = pt, p_any = pa, range_terms = rng_terms, verdict_terms = verdict_terms,
			range_lead = rng_any, verdict = .bay_sensitivity_verdict(max(p_null, 1 - p_null), rng_any))
		moved <- which(verdict_terms == "changes")
		.bay_step_result(story, if (length(moved)) paste0(paste(terms[moved], collapse = ", "), " change side with the prior") else "no term changes side",
			if (length(moved)) paste0(paste(terms[moved], collapse = ", "), " changent de c\u00f4t\u00e9 selon l'a priori") else "aucun terme ne change de c\u00f4t\u00e9")
		.bay_tell(story, "three priors on the models (HOP, independent terms, uniform); no scale r with unit information.",
			"trois lois sur les mod\u00e8les (HOP, termes ind\u00e9pendants, uniforme) ; pas d'\u00e9chelle r avec l'information unitaire.", "note")
	} else {
		.bay_step_result(story, "sensitivity not computed", "sensibilit\u00e9 non calcul\u00e9e")
	}

	# Step 6: size of each probable main effect with two levels, from its own model fit.
	ratios <- NULL
	.bay_tell(story, if (bin) "By how much the share changes" else if (ord) "By how much the values shift" else "By how much the counts change",
		if (bin) "De combien la part change" else if (ord) "De combien les valeurs se d\u00e9placent" else "De combien les comptages changent", "step",
		fun = if (ord) "polr() {MASS}" else "glm() {stats}")
	main2 <- which(lengths(strsplit(terms, ":", fixed = TRUE)) == 1L & eff$p_h1 >= 0.5)
	main2 <- main2[vapply(terms[main2], function(v) nlevels(factors[[v]]) == 2L, logical(1))]
	if (length(main2)) {
		dd <- data.frame(y = y01, factors, check.names = FALSE)
		fo <- stats::as.formula(paste("y ~", paste0("`", names(factors), "`", collapse = " + ")))
		fit <- .bay_quiet(if (ord) MASS::polr(fo, data = dd, Hess = TRUE) else if (bin) stats::glm(fo, data = dd, family = stats::binomial()) else if (isTRUE(p_nb >= 0.5)) MASS::glm.nb(fo, data = dd) else stats::glm(fo, data = dd, family = stats::poisson()))
		co <- summary(fit)$coefficients
		for (i in main2) {
			v <- terms[i]
			nm <- paste0("`", v, "`", levels(factors[[v]])[2L])
			if (!nm %in% rownames(co)) nm <- paste0(v, levels(factors[[v]])[2L])
			if (!nm %in% rownames(co)) next
			est <- co[nm, 1L]; se <- co[nm, 2L]
			rr <- exp(c(est, est - 1.96 * se, est + 1.96 * se))
			ratios <- rbind(ratios, data.frame(term = v, group1 = levels(factors[[v]])[1L], group2 = levels(factors[[v]])[2L],
				ratio = rr[1L], low = rr[2L], high = rr[3L], stringsAsFactors = FALSE))
			lab_en <- if (bin) "odds ratio" else if (ord) "cumulative odds ratio" else "rate ratio"
			lab_fr <- if (bin) "rapport des cotes" else if (ord) "rapport des cotes cumul\u00e9es" else "rapport des moyennes"
			.bay_tell(story, paste0(v, ": ", lab_en, " ", levels(factors[[v]])[2L], "/", levels(factors[[v]])[1L], " = ", f2(rr[1L], "."), " [", f2(rr[2L], "."), "; ", f2(rr[3L], "."), "], at equal other factors."),
				paste0(v, " : ", lab_fr, " ", levels(factors[[v]])[2L], "/", levels(factors[[v]])[1L], " = ", f2(rr[1L], ","), " [", f2(rr[2L], ","), " ; ", f2(rr[3L], ","), "], \u00e0 autres facteurs \u00e9gaux."))
		}
		.bay_step_result(story, paste0(ratios$term[1L], ": ", ratios$group2[1L], "/", ratios$group1[1L], " = ", f2(ratios$ratio[1L], ".")),
			paste0(ratios$term[1L], " : ", ratios$group2[1L], "/", ratios$group1[1L], " = ", f2(ratios$ratio[1L], ",")))
		.bay_tell(story, "1 = no difference; [ ; ] = 95 % interval (Laplace approximation of the additive model).",
			"1 = pas de diff\u00e9rence ; [ ; ] = intervalle \u00e0 95 % (approximation de Laplace du mod\u00e8le additif).", "note")
	} else {
		.bay_step_result(story, "no probable main effect with two levels", "aucun effet principal probable \u00e0 deux niveaux")
		.bay_tell(story, "The cells are read in the bilan.", "Les cases se lisent au bilan.")
	}

	p_any <- 1 - p_null
	what_en <- if (bin) paste0("the share of \"", yes, "\"") else if (ord) paste0("the values of ", resp) else paste0("the counts of ", resp)
	what_fr <- if (bin) paste0("la part de \u00ab ", yes, " \u00bb") else if (ord) paste0("les valeurs de ", resp) else paste0("les comptages de ", resp)
	verb_en <- if (bin) " differs" else " differ"
	verb_fr <- if (bin) " diff\u00e8re" else " diff\u00e8rent"
	no_fr <- if (bin) paste0("Aucune diff\u00e9rence de ", what_fr) else if (ord) paste0("Aucune diff\u00e9rence des valeurs de ", resp) else paste0("Aucune diff\u00e9rence des comptages de ", resp)
	hyp <- data.frame(
		label_en = c(paste0(toupper(substr(what_en, 1L, 1L)), substring(what_en, 2L), verb_en, " according to at least one factor"), paste0("No difference in ", what_en, " according to ", fac_txt)),
		label_fr = c(paste0(toupper(substr(what_fr, 1L, 1L)), substring(what_fr, 2L), verb_fr, " selon au moins un facteur"), paste0(no_fr, " selon ", fac_txt)),
		posterior = c(p_any, p_null), prior = c(1 - exp(lp[1L]), exp(lp[1L])), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	out <- list(hypotheses = hyp, effects = eff,
		laws = if (bin) data.frame(law = "binomial", posterior = 1) else if (ord) data.frame(law = "polr", posterior = 1) else data.frame(law = c("poisson", "negbin"), posterior = c(1 - p_nb, p_nb)),
		models = data.frame(model = vapply(sp$models, function(m) if (length(m)) paste(terms[m], collapse = " + ") else "(null)", character(1)),
			log_bf = lb, log_prior = lp, stringsAsFactors = FALSE),
		questions = list(en = paste0("Does ", what_en, " differ, term by term?"), fr = paste0(what_fr, " diff\u00e8re-t-elle, terme par terme ?")),
		sensitivity = sens, posthoc = NULL, ratios = ratios, glm_warnings = n_warn,
		cell_stats = data.frame(cell = names(stat), n = as.integer(cn[names(stat)]), value = as.numeric(stat), stringsAsFactors = FALSE),
		kind = kind, priors = list(model_prior = "HOP(1,1)", effects = "unit information (BIC)"),
		contradiction = FALSE, heavy = FALSE, detail = isTRUE(detail), prob_effect = p_any,
		method = if (bin) "Yes/no with several factors: logistic regression, BIC approximation" else if (ord) "Ordered response with several factors: proportional odds (MASS::polr), BIC approximation" else "Counts with several factors: Poisson and negative binomial, BIC approximation",
		plan = "discrete_fact", label = input$label, n = length(y01), notes = get_notes(), call = call_expr)
	class(out) <- "mbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .bay_print_bilan_discrete_fact(out, detail = detail, narrated = TRUE)
	if (isFALSE(return)) return(p_any)
	if (isTRUE(verbose)) invisible(out) else out
}

.bay_print_bilan_discrete_fact <- function(obj, detail = FALSE, narrated = FALSE) {
	body <- function(emit, li, narrated) {
		e <- obj$effects
		ta <- .bay_term_answers(obj)
		for (i in seq_len(nrow(e))) {
			emit("fact", paste0(e$term[i], ": ", ta$long_en[i], " (", .bay_pct(e$p_h1[i], "."), "; prior ", .bay_pct(e$prior_h1[i], "."), "; Bayes factor ", .bay_bf_txt(e$log_bf[i], "."), ")."),
				paste0(e$term[i], " : ", ta$long_fr[i], " (", .bay_pct(e$p_h1[i], ","), " ; a priori ", .bay_pct(e$prior_h1[i], ","), " ; facteur de Bayes ", .bay_bf_txt(e$log_bf[i], ","), ")."))
		}
		emit("block", "", "")
		cs <- obj$cell_stats
		bin <- identical(obj$kind, "binary")
		ord <- identical(obj$kind, "ordinal")
		f2 <- function(v, dc) gsub(".", dc, formatC(v, format = "f", digits = 2), fixed = TRUE)
		emit("fact", paste0(if (bin) "Observed shares per cell: " else if (ord) "Mean rank per cell: " else "Mean counts per cell: ",
				paste0(cs$cell, " = ", if (bin) vapply(cs$value, .bay_pct, character(1), dec = ".") else f2(cs$value, "."), collapse = "; "), "."),
			paste0(if (bin) "Parts observ\u00e9es par case : " else if (ord) "Rang moyen par case : " else "Comptages moyens par case : ",
				paste0(cs$cell, " = ", if (bin) vapply(cs$value, .bay_pct, character(1), dec = ",") else f2(cs$value, ","), collapse = " ; "), "."))
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)
	say("Each hereditary model fitted by glm; marginal likelihood -BIC/2 (unit information); HOP prior on the models; counts: Poisson and negative binomial averaged at even odds (DECISIONS, D25).",
		"Chaque mod\u00e8le h\u00e9r\u00e9ditaire ajust\u00e9 par glm ; vraisemblance marginale -BIC/2 (information unitaire) ; a priori HOP sur les mod\u00e8les ; comptages : Poisson et binomiale n\u00e9gative moyenn\u00e9es \u00e0 poids \u00e9gaux (DECISIONS, D25).")
	invisible(obj)
}
