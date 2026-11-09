# sys_bay_factorial.R
# m.bay with several factors: normal law with a common variance only (DECISIONS, D1 of
# 9 November 2026). The marginal likelihood of each model comes from the BayesFactor
# package; KefiR only assembles: the space of hereditary models, the prior on these
# models (DECISIONS, D2), the probability of each term, the sensitivity and the bilan.
#
# References (comments only).
#   Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012). Default Bayes
#     factors for ANOVA designs. J. Math. Psych., 56, 356-374.
#     https://doi.org/10.1016/j.jmp.2012.08.001
#   Morey, R. D., & Rouder, J. N. BayesFactor: Computation of Bayes Factors for Common
#     Designs. R package, version 0.9.12-4.8 checked. https://doi.org/10.32614/CRAN.package.BayesFactor
#   Chipman, H. (1996). Bayesian variable selection with related predictors. Canadian
#     J. Statistics, 24, 17-36. https://doi.org/10.2307/3315687 (strong heredity)
#   Scott, J. G., & Berger, J. O. (2010). Bayes and empirical-Bayes multiplicity
#     adjustment in the variable-selection problem. Ann. Statist., 38, 2587-2619.
#     https://doi.org/10.1214/10-AOS792
#   Taylor-Rodriguez, D., Womack, A., & Bliznyuk, N. (2016). Bayesian variable selection
#     on model spaces constrained by heredity conditions. J. Comput. Graph. Statist., 25,
#     515-535. https://doi.org/10.1080/10618600.2015.1056793 (HOP prior, section 3.1)


# Hereditary model space. An interaction enters only with all its lower-order terms
# (strong heredity). Terms whose lower-order terms are missing from the formula are
# refused, rather than tested against a comparison nobody asked for.
.bay_fact_space <- function(terms) {
	keys <- .bay_term_key(terms)
	fac <- strsplit(keys, ":", fixed = TRUE)
	ord <- lengths(fac)
	parents <- lapply(seq_along(keys), function(i) {
		which(vapply(fac, function(f) length(f) < ord[i] && all(f %in% fac[[i]]), logical(1)))
	})
	short <- lengths(parents) < 2^ord - 2
	if (any(short)) {
		.exit(paste0("An interaction needs all its lower-order terms in the formula (write A * B): ",
				paste(terms[short], collapse = ", "), "."),
			paste0("Une interaction demande tous ses termes d'ordre inf\u00e9rieur dans la formule (\u00e9crire A * B) : ",
				paste(terms[short], collapse = ", "), "."))
	}
	models <- list(integer(0))
	for (i in order(ord)) {
		with_parents <- Filter(function(m) all(parents[[i]] %in% m), models)
		models <- c(models, lapply(with_parents, function(m) c(m, i)))
		if (length(models) > .bay_fact_max_models()) {
			.exit(paste0("More than ", .bay_fact_max_models(), " models to compare: write fewer interactions (for example A + B + C + A:B)."),
				paste0("Plus de ", .bay_fact_max_models(), " mod\u00e8les \u00e0 comparer : \u00e9crire moins d'interactions (par exemple A + B + C + A:B)."))
		}
	}
	list(terms = terms, keys = keys, order = ord, parents = parents, models = models)
}

.bay_fact_max_models <- function() 200L

# Log prior of each model.
#   "hop" : one pi ~ Beta(1,1) per order of interaction, integrated; within an order,
#           only the terms whose lower-order terms are all in the model can enter. This is
#           the HOP(1,1) of Taylor-Rodriguez et al. (2016), the hierarchical form of the
#           multiplicity correction of Scott & Berger (2010). Default (DECISIONS, D2).
#   "hip" : pi = 1/2 for each term that can enter (Chipman 1996).
#   "unif": uniform on the hereditary models.
# Each one is a proper law (a term-by-term draw, order by order): it sums to 1 on any
# hereditary space.
.bay_fact_prior <- function(sp, type = c("hop", "hip", "unif")) {
	type <- match.arg(type)
	n_mod <- length(sp$models)
	if (type == "unif") return(rep(-log(n_mod), n_mod))
	vapply(sp$models, function(m) {
		can <- which(vapply(seq_along(sp$keys), function(i) all(sp$parents[[i]] %in% m), logical(1)))
		if (type == "hip") return(length(can) * log(1 / 2))
		sum(vapply(unique(sp$order), function(o) {
			s <- sum(sp$order[m] == o)
			e <- sum(sp$order[can] == o)
			lbeta(1 + s, 1 + e - s)
		}, numeric(1)))
	}, numeric(1))
}

# Log BF of every model against the null model (intercept only), for each r.
# lmBF integrates over g by importance sampling when there are several terms: the seed
# of the call makes it reproducible.
.bay_fact_logbf <- function(y, factors, sp, rs, iterations, seed) {
	safe <- paste0("f", seq_along(factors))
	names(safe) <- names(factors)
	dat <- data.frame(.y = y, factors, check.names = FALSE)
	names(dat) <- c(".y", safe)
	term_safe <- vapply(strsplit(sp$keys, ":", fixed = TRUE), function(f) {
		paste(safe[f], collapse = ":")
	}, character(1))
	n_mod <- length(sp$models)
	lb <- matrix(0, n_mod, length(rs))
	err <- matrix(0, n_mod, length(rs))
	.bay_with_seed(seed, {
		for (j in seq_len(n_mod)[-1L]) {
			fo <- stats::as.formula(paste(".y ~", paste(term_safe[sp$models[[j]]], collapse = " + ")))
			for (ir in seq_along(rs)) {
				b <- BayesFactor::lmBF(fo, data = dat, rscaleFixed = rs[ir], iterations = iterations,
					progress = FALSE)
				x <- BayesFactor::extractBF(b, logbf = TRUE)
				lb[j, ir] <- x$bf
				err[j, ir] <- x$error
			}
		}
	})
	list(log_bf = lb, error = err)
}

# Probability of each term under one prior: marginal of a single joint law, so the
# terms can be read together (METHODOLOGIE, 3).
.bay_fact_terms <- function(sp, log_bf, log_prior) {
	tab <- lapply(sp$models, function(m) sp$keys[m])
	suppressWarnings(.bay_inclusion(log_bf, tab, type = "all", prior_prob = exp(log_prior)))
}

# Probability of the null model: no difference according to any factor.
.bay_fact_null <- function(log_bf, log_prior) {
	l <- log_bf + log_prior
	exp(l[1L] - .bay_log_sum_exp(l))
}

# Plain questions, without a cause. A main effect is read "alone or in interaction",
# because under heredity its probability includes the models where it interacts.
.bay_fact_questions <- function(terms, resp) {
	f <- strsplit(terms, ":", fixed = TRUE)
	en <- fr <- character(length(terms))
	for (i in seq_along(terms)) {
		if (length(f[[i]]) == 1L) {
			en[i] <- paste0("Does ", resp, " differ according to ", terms[i], ", alone or in interaction?")
			fr[i] <- paste0(resp, " diff\u00e8re-t-il selon ", terms[i], ", seul ou en interaction ?")
		} else {
			last <- f[[i]][length(f[[i]])]
			first <- paste(f[[i]][-length(f[[i]])], collapse = ":")
			en[i] <- paste0("Do the differences according to ", first, " vary with ", last, " (", terms[i], ")?")
			fr[i] <- paste0("Les \u00e9carts selon ", first, " varient-ils avec ", last, " (", terms[i], ") ?")
		}
	}
	list(en = en, fr = fr)
}

# Whole analysis with several factors. input: from .bay_parse_input; note: the notes
# collector of m.bay.
.bay_factorial <- function(input, rscale, sensitivity, draws, seed, note) {
	y <- input$response
	factors <- input$factors
	terms <- input$terms
	sp <- .bay_fact_space(terms)

	# An interaction with an empty cell cannot be estimated.
	for (i in which(sp$order > 1L)) {
		f <- strsplit(sp$keys[i], ":", fixed = TRUE)[[1]]
		cells <- table(interaction(factors[f], drop = FALSE))
		if (any(cells == 0L)) {
			.exit(paste0("Empty cell for ", terms[i], ": this interaction cannot be estimated, remove it from the formula."),
				paste0("Case vide pour ", terms[i], " : cette interaction ne peut pas \u00eatre estim\u00e9e, la retirer de la formule."))
		}
	}
	n_par <- 1L + sum(vapply(sp$keys, function(k) {
		prod(vapply(strsplit(k, ":", fixed = TRUE)[[1]], function(v) nlevels(factors[[v]]) - 1L, numeric(1)))
	}, numeric(1)))
	if (length(y) <= n_par) {
		.exit("No residual degree of freedom for the full model: fewer terms are needed.",
			"Aucun degr\u00e9 de libert\u00e9 r\u00e9siduel pour le mod\u00e8le complet : il faut moins de termes.")
	}

	rs <- if (isTRUE(sensitivity)) c(rscale / 2, rscale, rscale * 2) else rscale
	ic <- if (isTRUE(sensitivity)) 2L else 1L
	iterations <- max(10000L, as.integer(draws))
	bf <- .bay_fact_logbf(y, factors, sp, rs, iterations, seed)
	lp <- .bay_fact_prior(sp, "hop")
	eff <- .bay_fact_terms(sp, bf$log_bf[, ic], lp)
	# Same order as the formula.
	eff <- eff[match(.bay_term_key(terms), .bay_term_key(eff$term)), , drop = FALSE]
	eff$term <- terms
	eff$rscale <- rscale
	err_pct <- 100 * max(bf$error[, ic], na.rm = TRUE)
	eff$error_pct <- err_pct
	rownames(eff) <- NULL
	p_null <- .bay_fact_null(bf$log_bf[, ic], lp)
	prior_null <- exp(lp[1L])

	sens <- NULL
	if (isTRUE(sensitivity)) {
		priors <- list(hop = lp, hip = .bay_fact_prior(sp, "hip"), unif = .bay_fact_prior(sp, "unif"))
		cells <- expand.grid(r = seq_along(rs), prior = names(priors), stringsAsFactors = FALSE)
		p_terms <- sapply(seq_len(nrow(cells)), function(cc) {
			e <- .bay_fact_terms(sp, bf$log_bf[, cells$r[cc]], priors[[cells$prior[cc]]])
			e$p_h1[match(.bay_term_key(terms), .bay_term_key(e$term))]
		})
		p_terms <- matrix(p_terms, nrow = length(terms))
		p_any <- vapply(seq_len(nrow(cells)), function(cc) {
			1 - .bay_fact_null(bf$log_bf[, cells$r[cc]], priors[[cells$prior[cc]]])
		}, numeric(1))
		lead_any <- (1 - p_null) >= p_null
		rng_any <- range(if (lead_any) p_any else 1 - p_any)
		rng_terms <- t(apply(p_terms, 1L, range))
		verdict_terms <- vapply(seq_along(terms), function(i) {
			lead <- eff$p_h1[i] >= 0.5
			r <- if (lead) rng_terms[i, ] else 1 - rev(rng_terms[i, ])
			.bay_sensitivity_verdict(max(eff$p_h1[i], eff$p_h0[i]), r)
		}, character(1))
		sens <- list(r = rs, priors = names(priors), cells = cells, p_terms = p_terms, p_any = p_any,
			range_terms = rng_terms, verdict_terms = verdict_terms, range_lead = rng_any,
			verdict = .bay_sensitivity_verdict(max(p_null, 1 - p_null), rng_any))
	}
	if (isTRUE(err_pct > 5)) {
		note(paste0("Numerical precision is low: error on a Bayes factor up to ",
				format(signif(err_pct, 2), trim = TRUE), " %. Try a larger draws."),
			paste0("Pr\u00e9cision num\u00e9rique faible : erreur sur un facteur de Bayes jusqu'\u00e0 ",
				gsub(".", ",", format(signif(err_pct, 2), trim = TRUE), fixed = TRUE), " %. Essayer un draws plus grand."))
	}

	cells_n <- table(interaction(factors, drop = FALSE))
	balanced <- length(unique(as.integer(cells_n))) == 1L
	means <- lapply(names(factors), function(v) {
		data.frame(factor = v, level = levels(factors[[v]]),
			n = as.integer(table(factors[[v]])),
			mean = as.numeric(tapply(y, factors[[v]], mean)), stringsAsFactors = FALSE)
	})
	list(effects = eff, p_null = p_null, prior_null = prior_null, sensitivity = sens,
		space = sp, log_bf_models = bf$log_bf, error_models = bf$error, log_prior = lp,
		means = do.call(rbind, means), balanced = balanced, min_cell = min(cells_n),
		sd_within = stats::sigma(stats::lm(y ~ ., data = data.frame(y = y, factors))))
}

# Bilan with several factors. Same order as with one factor: verdict on the global
# question, then each term, its sensitivity, what was compared, the starting point.
.bay_print_bilan_fact <- function(obj, detail = FALSE) {
	lang_fr <- identical(.msg("en", "fr"), "fr")
	dec <- if (lang_fr) "," else "."
	say <- function(en, fr) cat(.msg(en, fr), "\n", sep = "")
	num <- function(x) gsub(".", dec, format(signif(x, 3), trim = TRUE), fixed = TRUE)
	say(paste0("Bilan of m.bay: ", obj$label, " (N = ", obj$n, ", several factors)"),
		paste0("Bilan de m.bay : ", obj$label, " (N = ", obj$n, ", plusieurs facteurs)"))
	v <- .bay_verdict(obj)
	say(v$en, v$fr)
	h <- obj$hypotheses
	sup <- .bay_support(h$posterior[1L])
	say(paste0("  ", h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), " (", sup$en, "); ",
			h$label_en[2L], ": ", .bay_pct(h$posterior[2L], "."), "."),
		paste0("  ", h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), " (", sup$fr, ") ; ",
			h$label_fr[2L], " : ", .bay_pct(h$posterior[2L], ","), "."))

	say("Answer per term, probability given the data and the starting point (prior in brackets):",
		"R\u00e9ponse par terme, probabilit\u00e9 compte tenu des donn\u00e9es et du point de d\u00e9part (a priori entre crochets) :")
	e <- obj$effects
	sens <- obj$sensitivity
	for (i in seq_len(nrow(e))) {
		s <- .bay_support(max(e$p_h1[i], e$p_h0[i]))
		lead_en <- if (e$p_h1[i] >= 0.5) "yes" else "no"
		lead_fr <- if (e$p_h1[i] >= 0.5) "oui" else "non"
		tail_en <- tail_fr <- ""
		if (!is.null(sens)) {
			rr <- sens$range_terms[i, ]
			if (identical(sens$verdict_terms[i], "changes")) {
				tail_en <- paste0(" Depending on the starting point ", .bay_pct(rr[1L], "."), " to ", .bay_pct(rr[2L], "."), ": not settled.")
				tail_fr <- paste0(" Selon le point de d\u00e9part ", .bay_pct(rr[1L], ","), " \u00e0 ", .bay_pct(rr[2L], ","), " : non tranch\u00e9.")
			}
		}
		say(paste0("  ", obj$questions$en[i], " Yes ", .bay_pct(e$p_h1[i], "."), " [", .bay_pct(e$prior_h1[i], "."), "], ",
				lead_en, " with ", s$en, ".", tail_en),
			paste0("  ", obj$questions$fr[i], " Oui ", .bay_pct(e$p_h1[i], ","), " [", .bay_pct(e$prior_h1[i], ","), "], ",
				lead_fr, " avec ", s$fr, ".", tail_fr))
	}
	if (!is.null(sens)) {
		if (all(sens$verdict_terms != "changes") && !identical(sens$verdict, "changes")) {
			say("  No answer changes side when r or the prior on the models changes.",
				"  Aucune r\u00e9ponse ne change de c\u00f4t\u00e9 quand r ou la loi sur les mod\u00e8les change.")
		}
	} else {
		say("  Sensitivity to the priors not computed (sensitivity = FALSE).",
			"  Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	m <- obj$descriptive_stats
	dg <- max(0L, 2L - floor(log10(max(abs(m$mean), 1e-300))))
	fmt <- function(x, dc) gsub(".", dc, formatC(x, format = "f", digits = dg), fixed = TRUE)
	say("Observed means per factor (description, not a test of which levels differ):",
		"Moyennes observ\u00e9es par facteur (description ; ce bilan ne dit pas quels niveaux diff\u00e8rent) :")
	for (f in unique(m$factor)) {
		mf <- m[m$factor == f, , drop = FALSE]
		say(paste0("  ", f, ": ", paste0(mf$level, " = ", fmt(mf$mean, "."), " (n = ", mf$n, ")", collapse = "; ")),
			paste0("  ", f, " : ", paste0(mf$level, " = ", fmt(mf$mean, ","), " (n = ", mf$n, ")", collapse = " ; ")))
	}
	say(paste0("  Within-cell standard deviation: ", fmt(obj$sd_within, "."), if (!isTRUE(obj$balanced)) "; unbalanced plan, these means mix the other factors." else "."),
		paste0("  \u00c9cart-type au sein des cases : ", fmt(obj$sd_within, ","), if (!isTRUE(obj$balanced)) " ; plan d\u00e9s\u00e9quilibr\u00e9, ces moyennes m\u00e9langent les autres facteurs." else "."))

	say("Shape of the noise: normal with a common variance only. With several factors, heavy tails and unequal variances are not compared yet.",
		"Forme du bruit : normale \u00e0 variance commune seulement. \u00c0 plusieurs facteurs, queues lourdes et variances in\u00e9gales ne sont pas encore compar\u00e9es.")
	say(paste0("Starting point, set by the method before seeing the data: expected size r = ", format(signif(obj$priors$rscale, 3), trim = TRUE),
			"; an interaction enters only with its factors; one chance in two per factor, interactions less likely ([] above)."),
		paste0("Point de d\u00e9part, fix\u00e9 par la m\u00e9thode avant de voir les donn\u00e9es : taille attendue r = ", num(obj$priors$rscale),
			" ; une interaction n'entre qu'avec ses facteurs ; une chance sur deux par facteur, les interactions moins probables ([] ci-dessus)."))
	say("These are probabilities of hypotheses, not p-values: no threshold of 5 % applies. Whether a difference has a cause depends on the protocol (randomisation, order, operator).",
		"Ce sont des probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values : aucun seuil de 5 % ne s'y applique. La cause d'une diff\u00e9rence d\u00e9pend du protocole (tirage au sort, ordre, op\u00e9rateur).")
	if (length(obj$notes$en)) {
		say("Warnings:", "Avertissements :")
		for (i in seq_along(obj$notes$en)) say(paste0("  ", obj$notes$en[i]), paste0("  ", obj$notes$fr[i]))
	}
	if (!isTRUE(detail)) {
		say("Bayes factors, priors in detail, how to read: print(result, detail = TRUE).",
			"Facteurs de Bayes, a priori en d\u00e9tail, mode de lecture : print(resultat, detail = TRUE).")
		return(invisible(obj))
	}

	say("Bayes factor per term (how much the data move the odds of the term):",
		"Facteur de Bayes par terme (de combien les donn\u00e9es d\u00e9placent les cotes du terme) :")
	say(paste0("  ", paste0(e$term, " ", vapply(e$log_bf, .bay_bf_txt, character(1), dec = "."), collapse = "; "), "."),
		paste0("  ", paste0(e$term, " ", vapply(e$log_bf, .bay_bf_txt, character(1), dec = ","), collapse = " ; "), "."))
	say("Priors in detail:", "A priori en d\u00e9tail :")
	say(paste0("  Effects: Zellner-Siow (Cauchy) priors of scale r = ", format(signif(obj$priors$rscale, 3), trim = TRUE),
			", one per term (Rouder et al. 2012), computed by the BayesFactor package."),
		paste0("  Effets : a priori de Zellner-Siow (Cauchy) d'\u00e9chelle r = ", num(obj$priors$rscale),
			", un par terme (Rouder et al. 2012), calcul\u00e9s par le paquet BayesFactor."))
	say(paste0("  Models: ", obj$n_models, " models respecting heredity; one proportion ~ Beta(1,1) per order of interaction (Taylor-Rodriguez et al. 2016, HOP), which corrects for the number of terms. No difference at all: ",
			.bay_pct(obj$hypotheses$prior[obj$hypotheses$null], "."), " a priori."),
		paste0("  Mod\u00e8les : ", obj$n_models, " mod\u00e8les respectant l'h\u00e9r\u00e9dit\u00e9 ; une proportion ~ Beta(1,1) par ordre d'interaction (Taylor-Rodriguez et al. 2016, HOP), qui corrige du nombre de termes. Aucune diff\u00e9rence du tout : ",
			.bay_pct(obj$hypotheses$prior[obj$hypotheses$null], ","), " a priori."))
	if (!is.null(sens)) {
		say(paste0("  Varied in the sensitivity: r from ", format(min(sens$r), trim = TRUE), " to ", format(max(sens$r), trim = TRUE),
				"; prior on the models: HOP, 1/2 per term (Chipman 1996), uniform."),
			paste0("  Vari\u00e9s dans la sensibilit\u00e9 : r de ", num(min(sens$r)), " \u00e0 ", num(max(sens$r)),
				" ; loi sur les mod\u00e8les : HOP, 1/2 par terme (Chipman 1996), uniforme."))
		for (i in seq_len(nrow(e))) {
			rr <- sens$range_terms[i, ]
			say(paste0("    ", e$term[i], ": ", .bay_pct(rr[1L], "."), " to ", .bay_pct(rr[2L], ".")),
				paste0("    ", e$term[i], " : ", .bay_pct(rr[1L], ","), " \u00e0 ", .bay_pct(rr[2L], ",")))
		}
	}
	say("How to read it:", "Comment lire :")
	say("  The probability of a factor includes the models where it acts in interaction: it answers \"does the factor play a role\", not \"is its average effect zero\".",
		"  La probabilit\u00e9 d'un facteur inclut les mod\u00e8les o\u00f9 il agit en interaction : elle r\u00e9pond \u00e0 \u00ab le facteur joue-t-il un r\u00f4le \u00bb, pas \u00e0 \u00ab son effet moyen est-il nul \u00bb.")
	say("  A p-value is the probability of data at least this extreme if there were no difference: neither of the above.",
		"  Une p-value est la probabilit\u00e9 de donn\u00e9es au moins aussi extr\u00eames s'il n'y avait pas de diff\u00e9rence : ni l'un ni l'autre ci-dessus.")
	invisible(obj)
}
