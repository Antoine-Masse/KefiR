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
.bay_fact_space <- function(terms, always = character(0)) {
	# always: variables present in every model (covariates, D16); a lower-order term made of
	# them only is not required.
	keys <- .bay_term_key(terms)
	fac <- strsplit(keys, ":", fixed = TRUE)
	ord <- lengths(fac)
	parents <- lapply(seq_along(keys), function(i) {
		which(vapply(fac, function(f) length(f) < ord[i] && all(f %in% fac[[i]]), logical(1)))
	})
	n_always <- vapply(fac, function(f) sum(f %in% always), numeric(1))
	short <- lengths(parents) < (2^ord - 2) - (2^n_always - 1)
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
# lmBF integrates over g by importance sampling when there are several terms. Each r has
# its own seed (seeds[ir]), so the central answer does not depend on whether the
# sensitivity is computed (RELECTURE 2026-11-09, point 6).
.bay_fact_logbf <- function(y, factors, sp, rs, iterations, seed, seeds = seed + seq_along(rs) - 1L, id = NULL,
                           cov = NULL, rconts = .bay_rcont() * rs / rs[1L]) {
	# Covariates (ANCOVA, D16): centred by the caller, present in every model; each Bayes
	# factor is read against the model with the covariates only; their slope prior is
	# rscaleCont, varied with r in the sensitivity.
	# Paired or repeated measures (DECISIONS, D15): the subject is a random effect present in
	# every model, "no difference" included; each Bayes factor is read against the model with
	# the subjects only. Two levels of one factor and complete pairs: ttestBF on the
	# differences, exact (no sampling), r on the standardised difference.
	paired <- !is.null(id) && length(factors) == 1L && nlevels(factors[[1L]]) == 2L
	safe <- paste0("f", seq_along(factors))
	names(safe) <- names(factors)
	dat <- data.frame(.y = y, factors, check.names = FALSE)
	names(dat) <- c(".y", safe)
	if (!is.null(id)) dat$.id <- id
	has_cov <- !is.null(cov) && ncol(cov) > 0L
	cov_txt <- ""
	if (has_cov) {
		cs <- paste0("c", seq_along(cov))
		names(cs) <- names(cov)
		for (i in seq_along(cov)) dat[[cs[i]]] <- cov[[i]]
		safe <- c(safe, cs)
		cov_txt <- paste(cs, collapse = " + ")
	}
	term_safe <- vapply(strsplit(sp$keys, ":", fixed = TRUE), function(f) {
		paste(safe[f], collapse = ":")
	}, character(1))
	n_mod <- length(sp$models)
	lb <- matrix(0, n_mod, length(rs))
	err <- matrix(0, n_mod, length(rs))
	for (ir in seq_along(rs)) {
		.bay_with_seed(seeds[ir], {
			lb0 <- 0
			if (has_cov) {
				lb0 <- BayesFactor::extractBF(BayesFactor::lmBF(stats::as.formula(paste(".y ~", cov_txt)), data = dat,
					rscaleCont = rconts[ir], iterations = iterations, progress = FALSE), logbf = TRUE)$bf
			} else if (!is.null(id) && !paired) {
				lb0 <- BayesFactor::extractBF(BayesFactor::lmBF(.y ~ .id, data = dat, whichRandom = ".id",
					iterations = iterations, progress = FALSE), logbf = TRUE)$bf
			}
			for (j in seq_len(n_mod)[-1L]) {
				m <- sp$models[[j]]
				if (paired) {
					pr <- .bay_paired_xy(y, factors[[1L]], id)
					x <- BayesFactor::extractBF(BayesFactor::ttestBF(x = pr$b, y = pr$a, paired = TRUE, rscale = rs[ir]), logbf = TRUE)
					lb[j, ir] <- x$bf
					err[j, ir] <- x$error
					next
				}
				# One main effect alone: the internal one-factor computation, equal to lmBF to
				# 1e-4 (test-bay-factorial.R) and robust where the integrate() of lmBF fails, at
				# very large effects (agent_50: g = 2e4 gave NA).
				if (is.null(id) && !has_cov && length(m) == 1L && sp$order[m] == 1L) {
					lb[j, ir] <- .bay_bf_oneway_zs(y, factors[[sp$keys[m]]], rs[ir])
					next
				}
				fo <- stats::as.formula(paste(".y ~", paste(c(if (has_cov) cov_txt, term_safe[m]), collapse = " + "), if (!is.null(id)) "+ .id" else ""))
				x <- tryCatch(suppressMessages(BayesFactor::extractBF(BayesFactor::lmBF(fo, data = dat,
					whichRandom = if (!is.null(id)) ".id" else NULL,
					rscaleFixed = rs[ir], rscaleCont = rconts[ir], iterations = iterations, progress = FALSE), logbf = TRUE)),
					error = function(e) list(bf = NA_real_, error = NA_real_))
				lb[j, ir] <- x$bf - lb0
				err[j, ir] <- x$error
			}
		})
	}
	list(log_bf = lb, error = err)
}

# Prior scale of a covariate slope: the "medium" of BayesFactor (rscaleCont).
.bay_rcont <- function() sqrt(2) / 4

# Terms of the model space: without the covariates alone, which are in every model.
.bay_space_terms <- function(terms, cov_names) {
	if (!length(cov_names)) return(terms)
	terms[!vapply(strsplit(terms, ":", fixed = TRUE), function(v) all(v %in% cov_names), logical(1))]
}

# Overlap of each covariate between groups (D16): the common range is from the largest of
# the group minima to the smallest of the maxima; poor when it is empty or holds less than
# half of the values of one group: the adjustment then extrapolates.
.bay_cov_overlap <- function(cov, groups) {
	if (is.null(cov) || !ncol(cov)) return(NULL)
	groups <- droplevels(as.factor(groups))
	out <- lapply(names(cov), function(v) {
		x <- cov[[v]]
		lo <- max(tapply(x, groups, min)); hi <- min(tapply(x, groups, max))
		fr <- if (lo <= hi) tapply(x >= lo & x <= hi, groups, mean) else stats::setNames(rep(0, nlevels(groups)), levels(groups))
		data.frame(covariate = v, low = lo, high = hi, min_share = min(fr), worst = names(fr)[which.min(fr)],
			poor = lo > hi || min(fr) < 0.5, stringsAsFactors = FALSE)
	})
	do.call(rbind, out)
}

# Complete pairs, matched by subject: a = first level, b = second.
.bay_paired_xy <- function(y, f, id) {
	lev <- levels(f)
	ia <- which(f == lev[1L]); ib <- which(f == lev[2L])
	ids <- intersect(as.character(id[ia]), as.character(id[ib]))
	list(a = y[ia][match(ids, as.character(id[ia]))], b = y[ib][match(ids, as.character(id[ib]))], id = ids)
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
.bay_fact_questions <- function(terms, resp, adj = character(0)) {
	f <- strsplit(terms, ":", fixed = TRUE)
	en <- fr <- character(length(terms))
	for (i in seq_along(terms)) {
		if (length(terms) == 1L) {
			en[i] <- paste0("Difference in ", resp, " according to ", terms[i], if (length(adj)) paste0(", at equal ", paste(adj, collapse = ", ")) else "", "?")
			fr[i] <- paste0("Diff\u00e9rence de ", resp, " selon ", terms[i], if (length(adj)) paste0(", \u00e0 m\u00eame ", paste(adj, collapse = ", ")) else "", " ?")
		} else if (length(f[[i]]) == 1L) {
			en[i] <- paste0("Difference in ", resp, " according to ", terms[i], if (length(adj)) paste0(" at equal ", paste(adj, collapse = ", ")) else "", ", alone or in interaction?")
			fr[i] <- paste0("Diff\u00e9rence de ", resp, " selon ", terms[i], if (length(adj)) paste0(" \u00e0 m\u00eame ", paste(adj, collapse = ", ")) else "", ", seul ou en interaction ?")
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
.bay_factorial <- function(input, rscale, sensitivity, draws, seed, note, prior_h1 = 0.5) {
	y <- input$response
	id <- input$id
	paired <- !is.null(id) && length(input$factors) == 1L && nlevels(input$factors[[1L]]) == 2L
	factors <- input$factors
	# Non-syntactic names come back from the formula with backticks (RELECTURE 2026-11-09, 1).
	cov <- input$covariates
	cov_names <- if (!is.null(cov)) names(cov) else character(0)
	if (length(cov_names)) cov <- as.data.frame(lapply(cov, function(v) v - mean(v)), check.names = FALSE) else cov <- NULL
	terms <- .bay_space_terms(gsub("`", "", input$terms, fixed = TRUE), cov_names)
	sp <- .bay_fact_space(terms, cov_names)

	# An interaction with an empty cell cannot be estimated.
	for (i in which(sp$order > 1L)) {
		f <- setdiff(strsplit(sp$keys[i], ":", fixed = TRUE)[[1]], cov_names)
		if (length(f) < 2L) next
		cells <- table(interaction(factors[f], drop = FALSE))
		if (any(cells == 0L)) {
			.exit(paste0("Empty cell for ", terms[i], ": this interaction cannot be estimated, remove it from the formula."),
				paste0("Case vide pour ", terms[i], " : cette interaction ne peut pas \u00eatre estim\u00e9e, la retirer de la formule."))
		}
	}
	# Confounded or nested factors (one determines the other): the terms cannot be told apart,
	# and the answers would be misleading (RELECTURE 2026-11-09, 2).
	fx <- stats::model.matrix(stats::as.formula(paste("~", paste0("`", names(factors), "`", collapse = " + "))),
		data = factors)
	if (qr(fx)$rank < ncol(fx)) {
		.exit("Some factors are confounded or nested (one determines another): they cannot be compared as crossed factors. Nested plans are not handled yet.",
			"Des facteurs sont confondus ou embo\u00eet\u00e9s (l'un d\u00e9termine l'autre) : on ne peut pas les comparer comme des facteurs crois\u00e9s. Les plans embo\u00eet\u00e9s ne sont pas encore trait\u00e9s.")
	}
	n_par <- 1L + sum(vapply(sp$keys, function(k) {
		prod(vapply(strsplit(k, ":", fixed = TRUE)[[1]], function(v) if (v %in% cov_names) 1 else nlevels(factors[[v]]) - 1L, numeric(1)))
	}, numeric(1)))
	if (!is.null(id)) n_par <- n_par + nlevels(id) - 1L
	n_par <- n_par + length(cov_names)
	if (length(y) <= n_par) {
		.exit("No residual degree of freedom for the full model: fewer terms are needed.",
			"Aucun degr\u00e9 de libert\u00e9 r\u00e9siduel pour le mod\u00e8le complet : il faut moins de termes.")
	}

	# Central r first, with the seed of the call; r/2 and 2r with their own seeds.
	rs <- if (isTRUE(sensitivity)) c(rscale, rscale / 2, rscale * 2) else rscale
	iterations <- max(10000L, as.integer(draws))
	bf <- .bay_fact_logbf(y, factors, sp, rs, iterations, seed, seeds = seed + c(0L, 1L, 2L)[seq_along(rs)], id = id,
		cov = cov, rconts = .bay_rcont() * rs / rscale)
	# One term alone: its prior is prior_h1, as with one factor between groups.
	single <- length(sp$keys) == 1L
	lp <- if (single) log(c(1 - prior_h1, prior_h1)) else .bay_fact_prior(sp, "hop")
	eff <- .bay_fact_terms(sp, bf$log_bf[, 1L], lp)
	# Same order as the formula.
	eff <- eff[match(.bay_term_key(terms), .bay_term_key(eff$term)), , drop = FALSE]
	eff$term <- terms
	eff$rscale <- rscale
	rownames(eff) <- NULL
	p_null <- .bay_fact_null(bf$log_bf[, 1L], lp)
	prior_null <- exp(lp[1L])

	# Numerical precision, measured rather than taken from BayesFactor, whose reported error
	# is unreliable on models with several terms, too low or too high, and hardly shrinks
	# with more iterations (agent_51; RELECTURE 2026-11-09, 3): the central computation is
	# redone with another seed, and the largest change of a probability is reported.
	multi <- lengths(sp$models) >= 2L
	spread <- 0
	if (any(multi) || (!is.null(id) && !paired) || !is.null(cov)) {
		bf2 <- .bay_fact_logbf(y, factors, sp, rscale, iterations, seed, seeds = seed + 7919L, id = id,
			cov = cov, rconts = .bay_rcont())
		e2 <- .bay_fact_terms(sp, bf2$log_bf[, 1L], lp)
		p2 <- e2$p_h1[match(.bay_term_key(terms), .bay_term_key(e2$term))]
		spread <- max(abs(p2 - eff$p_h1), abs(.bay_fact_null(bf2$log_bf[, 1L], lp) - p_null), na.rm = TRUE)
	}
	eff$error_pct <- 100 * spread
	if (anyNA(bf$log_bf)) {
		note("BayesFactor could not compute some models (very large effects); the terms that depend on them are left at NA.",
			"BayesFactor n'a pas pu calculer certains mod\u00e8les (effets tr\u00e8s grands) ; les termes qui en d\u00e9pendent restent \u00e0 NA.")
	}
	if (isTRUE(spread > 0.02)) {
		note(paste0("Numerical precision: two computations with different seeds differ by up to ",
				format(round(100 * spread, 1), trim = TRUE), " points of probability."),
			paste0("Pr\u00e9cision num\u00e9rique : deux calculs avec des graines diff\u00e9rentes diff\u00e8rent jusqu'\u00e0 ",
				gsub(".", ",", format(round(100 * spread, 1), trim = TRUE), fixed = TRUE), " points de probabilit\u00e9."))
	}

	sens <- NULL
	if (isTRUE(sensitivity)) {
		priors <- if (single) {
			list(p = lp, low = log(c(2 / 3, 1 / 3)), high = log(c(1 / 3, 2 / 3)))
		} else {
			list(hop = lp, hip = .bay_fact_prior(sp, "hip"), unif = .bay_fact_prior(sp, "unif"))
		}
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
		sens <- list(r = rs[order(rs)], priors = names(priors), cells = cells, p_terms = p_terms, p_any = p_any,
			range_terms = rng_terms, verdict_terms = verdict_terms, range_lead = rng_any,
			verdict = .bay_sensitivity_verdict(max(p_null, 1 - p_null), rng_any))
	}

	cells_n <- table(interaction(factors, drop = FALSE))
	balanced <- length(unique(as.integer(cells_n))) == 1L
	means <- lapply(names(factors), function(v) {
		data.frame(factor = v, level = levels(factors[[v]]),
			n = as.integer(table(factors[[v]])),
			mean = as.numeric(tapply(y, factors[[v]], mean)), stringsAsFactors = FALSE)
	})
	# Within-cell standard deviation, from the cell means (RELECTURE 2026-11-09, 4). With one
	# observation per cell there is none: the residual of the largest model compared is given.
	cell <- interaction(factors, drop = TRUE)
	within <- length(y) > nlevels(cell)
	sd_within <- if (!is.null(cov)) {
		# Residual once the covariates and the cells are accounted for.
		stats::sigma(stats::lm(y ~ ., data = data.frame(y = y, cell = cell, cov, check.names = FALSE)))
	} else if (!is.null(id)) {
		# Residual within subject: what the subjects' own levels do not explain.
		stats::sigma(stats::lm(y ~ id + cell))
	} else if (within) {
		stats::sigma(stats::lm(y ~ cell))
	} else {
		stats::sigma(stats::lm(stats::as.formula(paste("y ~", paste0("`", names(factors), "`", collapse = " + "))),
			data = cbind(y = y, factors)))
	}
	list(effects = eff, p_null = p_null, prior_null = prior_null, sensitivity = sens,
		space = sp, log_bf_models = bf$log_bf[, order(rs), drop = FALSE], error_models = bf$error, log_prior = lp,
		means = do.call(rbind, means), balanced = balanced, min_cell = min(cells_n),
		sd_within = sd_within, sd_is_within = within, spread = spread, single = single, paired = paired,
		n_subjects = if (!is.null(id)) nlevels(id) else NA_integer_, terms = terms,
		overlap = .bay_cov_overlap(input$covariates, interaction(factors, drop = TRUE)),
		cov_centres = if (length(cov_names)) vapply(input$covariates, mean, numeric(1)) else NULL)
}

# Bilan with several factors. Same order as with one factor: verdict on the global
# question, then each term, its sensitivity, what was compared, the starting point.
.bay_print_bilan_fact <- function(obj, detail = FALSE, narrated = FALSE) {
	lang_fr <- identical(.msg("en", "fr"), "fr")
	dec <- if (lang_fr) "," else "."
	# As the last step of the story: indented, without what the steps already said.
	pre <- if (isTRUE(narrated)) "\t" else ""
	say <- function(en, fr) cat(pre, .msg(en, fr), "\n", sep = "")
	num <- function(x) gsub(".", dec, format(signif(x, 3), trim = TRUE), fixed = TRUE)
	rm <- obj$plan %in% c("paired", "repeated")
	poor <- isTRUE(any(obj$overlap$poor))
	if (!isTRUE(narrated)) {
		if (rm) {
			say(paste0("Bilan of m.bay: ", obj$label, " (N = ", obj$n, ", ", obj$n_subjects, " subjects, ", if (obj$plan == "paired") "paired data)" else "repeated measures)"),
				paste0("Bilan de m.bay : ", obj$label, " (N = ", obj$n, ", ", obj$n_subjects, " sujets, ", if (obj$plan == "paired") "donn\u00e9es appari\u00e9es)" else "mesures r\u00e9p\u00e9t\u00e9es)"))
		} else {
			say(paste0("Bilan of m.bay: ", obj$label, " (N = ", obj$n, ", several factors)"),
				paste0("Bilan de m.bay : ", obj$label, " (N = ", obj$n, ", plusieurs facteurs)"))
		}
	}
	v <- .bay_verdict(obj)
	say(v$en, v$fr)
	h <- obj$hypotheses
	# Global Bayes factor of the leading answer: posterior odds over prior odds.
	lb_glob <- log(h$posterior[1L] / h$posterior[2L]) - log(h$prior[1L] / h$prior[2L])
	sup <- .bay_support_bf(h$posterior[1L], lb_glob)
	if (poor && sup$level > 2L) sup <- .bay_support(0.7)
	say(paste0("  ", h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), " (", sup$en, "); ",
			h$label_en[2L], ": ", .bay_pct(h$posterior[2L], "."), "."),
		paste0("  ", h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), " (", sup$fr, ") ; ",
			h$label_fr[2L], " : ", .bay_pct(h$posterior[2L], ","), "."))
	e <- obj$effects
	single <- nrow(e) == 1L
	top <- which.max(e$p_h1)
	if (length(top) && !single) {
		top_open <- !is.null(obj$sensitivity) && identical(obj$sensitivity$verdict_terms[top], "changes")
		if (isTRUE(obj$prob_effect >= 0.5) && top_open) {
			say(paste0("     Rests on the terms together; none is settled alone (the most probable: ", e$term[top], ", ", .bay_pct(e$p_h1[top], "."), ")."),
				paste0("     Repose sur l'ensemble des termes ; aucun n'est tranch\u00e9 seul (le plus probable : ", e$term[top], ", ", .bay_pct(e$p_h1[top], ","), ")."))
		} else if (isTRUE(obj$prob_effect >= 0.5)) {
			say(paste0("     Rests on the most probable term: ", e$term[top], " (yes ", .bay_pct(e$p_h1[top], "."), ")."),
				paste0("     Repose sur le terme le plus probable : ", e$term[top], " (oui ", .bay_pct(e$p_h1[top], ","), ")."))
		} else {
			say(paste0("     Rests on this: no term reaches 50 %; the most probable, ", e$term[top], ", is at ", .bay_pct(e$p_h1[top], "."), "."),
				paste0("     Repose sur ceci : aucun terme n'atteint 50 % ; le plus probable, ", e$term[top], ", est \u00e0 ", .bay_pct(e$p_h1[top], ","), "."))
		}
	}
	if (!single && !is.null(obj$sensitivity) && identical(obj$sensitivity$verdict, "changes")) {
		rg <- obj$sensitivity$range_lead
		say(paste0("  Depending on the starting point, the probability of the first answer goes from ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), ": the order can reverse."),
			paste0("  Selon le point de d\u00e9part, la probabilit\u00e9 de la premi\u00e8re r\u00e9ponse va de ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " : l'ordre peut s'inverser."))
	}

	if (single) say("Answer (yes, no or not settled; probability of a yes given the data and the starting point, and its prior):",
		"R\u00e9ponse (oui, non ou non tranch\u00e9 ; probabilit\u00e9 du oui compte tenu des donn\u00e9es et du point de d\u00e9part, et son a priori) :") else say("Answer per term (yes, no or not settled; probability of a yes given the data and the starting point, and its prior):",
		"R\u00e9ponse par terme (oui, non ou non tranch\u00e9 ; probabilit\u00e9 du oui compte tenu des donn\u00e9es et du point de d\u00e9part, et son a priori) :")
	e <- obj$effects
	sens <- obj$sensitivity
	for (i in seq_len(nrow(e))) {
		yes <- e$p_h1[i] >= 0.5
		s <- .bay_support_bf(max(e$p_h1[i], e$p_h0[i]), e$log_bf[i] * if (yes) 1 else -1)
		if (poor && s$level > 2L) s <- .bay_support(0.7)
		against <- !is.na(e$log_bf[i]) && ((e$log_bf[i] >= 0) != yes)
		open <- s$level == 1L || against || (!is.null(sens) && identical(sens$verdict_terms[i], "changes"))
		rr_en <- rr_fr <- ""
		if (!is.null(sens) && identical(sens$verdict_terms[i], "changes")) {
			rr <- sens$range_terms[i, ]
			rr_en <- paste0("; depending on the starting point ", .bay_pct(rr[1L], "."), " to ", .bay_pct(rr[2L], "."))
			rr_fr <- paste0(" ; selon le point de d\u00e9part ", .bay_pct(rr[1L], ","), " \u00e0 ", .bay_pct(rr[2L], ","))
		}
		p_en <- paste0("probability of a yes ", .bay_pct(e$p_h1[i], "."), ", prior ", .bay_pct(e$prior_h1[i], "."))
		p_fr <- paste0("probabilit\u00e9 du oui ", .bay_pct(e$p_h1[i], ","), ", a priori ", .bay_pct(e$prior_h1[i], ","))
		if (open) {
			say(paste0("  ", obj$questions$en[i], " Not settled (", p_en, rr_en, ")."),
				paste0("  ", obj$questions$fr[i], " Non tranch\u00e9 (", p_fr, rr_fr, ")."))
		} else {
			say(paste0("  ", obj$questions$en[i], " ", if (yes) "Yes" else "No", ", ", s$en, " (", p_en, ")."),
				paste0("  ", obj$questions$fr[i], " ", if (yes) "Oui" else "Non", ", ", s$fr, " (", p_fr, ")."))
		}
		if (!is.na(e$log_bf[i])) {
			up <- e$log_bf[i] >= 0
			bf_en <- .bay_bf_txt(abs(e$log_bf[i]), ".")
			bf_fr <- .bay_bf_txt(abs(e$log_bf[i]), ",")
			v_en <- if (up) "multiply" else "divide"
			v_fr <- if (up) "multiplient" else "divisent"
			w_en <- if (yes) "yes" else "no"
			w_fr <- if (yes) "oui" else "non"
			if (open && up != yes) {
				say(paste0("     The data lean towards ", if (up) "yes" else "no", " (odds ", if (up) "multiplied" else "divided", " by ", bf_en, "); the prior (", .bay_pct(e$prior_h1[i], "."),
						") leans the other way: not settled."),
					paste0("     Les donn\u00e9es penchent vers le ", if (up) "oui" else "non", " (cotes ", if (up) "multipli\u00e9es" else "divis\u00e9es", " par ", bf_fr, ") ; l'a priori (", .bay_pct(e$prior_h1[i], ","),
						") penche de l'autre c\u00f4t\u00e9 : non tranch\u00e9."))
			} else if (open && abs(e$log_bf[i]) < log(3)) {
				say(paste0("     The data hardly move its odds (", v_en, " by ", bf_en, "); only the prior (", .bay_pct(e$prior_h1[i], "."), ") makes it lean towards ", w_en, "."),
					paste0("     Les donn\u00e9es d\u00e9placent peu ses cotes (", v_fr, " par ", bf_fr, ") ; seul l'a priori (", .bay_pct(e$prior_h1[i], ","), ") le fait pencher vers le ", w_fr, "."))
			} else if (open) {
				say(paste0("     The data ", v_en, " its odds (yes against no) by ", bf_en, ": too little to settle it."),
					paste0("     Les donn\u00e9es ", v_fr, " ses cotes (rapport oui/non) par ", bf_fr, " : trop peu pour trancher."))
			} else if (abs(e$log_bf[i]) < log(3)) {
				say(paste0("     The data hardly move its odds (", v_en, " by ", bf_en, "): the \"", w_en, "\" rests mainly on the prior of ",
						.bay_pct(e$prior_h1[i], "."), "."),
					paste0("     Les donn\u00e9es d\u00e9placent peu ses cotes (", v_fr, " par ", bf_fr, ") : le \u00ab ", w_fr, " \u00bb repose surtout sur l'a priori de ",
						.bay_pct(e$prior_h1[i], ","), "."))
			} else {
				say(paste0("     Rests on the data, which ", v_en, " its odds (yes against no) by ", bf_en, " (Bayes factor)."),
					paste0("     Repose sur les donn\u00e9es, qui ", v_fr, " ses cotes (rapport oui/non) par ", bf_fr, " (facteur de Bayes)."))
			}
		}
	}
	if (!is.null(sens)) {
		if (all(sens$verdict_terms != "changes") && !identical(sens$verdict, "changes")) {
			if (single) say("  The answer does not change side when r or the starting chances change.",
				"  La r\u00e9ponse ne change pas de c\u00f4t\u00e9 quand r ou les chances de d\u00e9part changent.") else say("  No answer changes side when r or the prior on the models changes.",
				"  Aucune r\u00e9ponse ne change de c\u00f4t\u00e9 quand r ou la loi sur les mod\u00e8les change.")
		}
	} else if (!isTRUE(narrated)) {
		say("  Sensitivity to the priors not computed (sensitivity = FALSE).",
			"  Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	if (poor) say("  Support capped at weak: the covariate ranges of the groups hardly overlap, the adjustment extrapolates.",
		"  Soutien plafonn\u00e9 \u00e0 faible : les \u00e9tendues de la covariable se recouvrent peu entre groupes, l'ajustement extrapole.")
	if (!is.null(obj$posthoc)) .bay_print_posthoc(obj$posthoc, say)
	m <- obj$descriptive_stats
	dg <- max(0L, 2L - floor(log10(max(abs(m$mean), 1e-300))))
	if (isTRUE(obj$sd_within > 0)) dg <- min(10L, max(dg, 1L - floor(log10(obj$sd_within))))
	fmt <- function(x, dc) gsub(".", dc, formatC(x, format = "f", digits = dg), fixed = TRUE)
	if (identical(obj$plan, "ancova")) {
		cg <- obj$cov_by_group
		cg_txt <- function(dc) paste(vapply(colnames(cg), function(v) paste0(v, " : ", paste0(rownames(cg), " = ", gsub(".", dc, format(signif(cg[, v], 3), trim = TRUE), fixed = TRUE), collapse = if (dc == ",") " ; " else "; ")), character(1)), collapse = " | ")
		say(paste0("Raw observed means per factor (not adjusted; mean covariate per group, ", cg_txt("."), "):"),
			paste0("Moyennes observ\u00e9es par facteur (brutes, non ajust\u00e9es ; covariable moyenne par groupe, ", cg_txt(","), ") :"))
	} else 	if (!is.null(obj$posthoc)) say("Observed means per factor (description of the sample):",
		"Moyennes observ\u00e9es par facteur (description de l'\u00e9chantillon) :") else say("Observed means per factor (description, not a test of which levels differ):",
		"Moyennes observ\u00e9es par facteur (description ; ce bilan ne dit pas quels niveaux diff\u00e8rent) :")
	for (f in unique(m$factor)) {
		mf <- m[m$factor == f, , drop = FALSE]
		say(paste0("  ", f, ": ", paste0(mf$level, " = ", fmt(mf$mean, "."), " (n = ", mf$n, ")", collapse = "; ")),
			paste0("  ", f, " : ", paste0(mf$level, " = ", fmt(mf$mean, ","), " (n = ", mf$n, ")", collapse = " ; ")))
	}
	cm <- obj$cell_means
	if (!is.null(cm) && any(grepl(":", e$term, fixed = TRUE)) && nrow(cm) <= 16L) {
		say(paste0("  Per cell: ", paste0(cm$cell, " = ", fmt(cm$mean, "."), collapse = "; "), "."),
			paste0("  Par case : ", paste0(cm$cell, " = ", fmt(cm$mean, ","), collapse = " ; "), "."))
	}
	anc <- identical(obj$plan, "ancova")
	sd_en <- if (anc) "Residual standard deviation (covariates accounted for): " else if (rm) "Residual standard deviation within subject (each subject's own level removed): " else if (isTRUE(obj$sd_is_within)) "Within-cell standard deviation: " else "Residual standard deviation (one observation per cell): "
	sd_fr <- if (anc) "\u00c9cart-type r\u00e9siduel (covariables prises en compte) : " else if (rm) "\u00c9cart-type r\u00e9siduel au sein des sujets (niveau propre de chaque sujet retir\u00e9) : " else if (isTRUE(obj$sd_is_within)) "\u00c9cart-type au sein des cases : " else "\u00c9cart-type r\u00e9siduel (une observation par case) : "
	say(paste0("  ", sd_en, fmt(obj$sd_within, "."), if (!isTRUE(obj$balanced)) "; unbalanced plan, these means mix the other factors." else "."),
		paste0("  ", sd_fr, fmt(obj$sd_within, ","), if (!isTRUE(obj$balanced)) " ; plan d\u00e9s\u00e9quilibr\u00e9, ces moyennes m\u00e9langent les autres facteurs." else "."))

	if (!isTRUE(narrated)) say("Shape of the noise: normal with a common variance only. With several factors, heavy tails and unequal variances are not compared yet.",
		"Forme du bruit : normale \u00e0 variance commune seulement. \u00c0 plusieurs facteurs, queues lourdes et variances in\u00e9gales ne sont pas encore compar\u00e9es.")
	has_int <- any(grepl(":", e$term, fixed = TRUE))
	if (!isTRUE(narrated)) say(paste0("Starting point, set by the method before seeing the data: expected size r = ", format(signif(obj$priors$rscale, 3), trim = TRUE),
			if (has_int) "; an interaction enters only with its factors; one chance in two per factor, interactions less likely ([] above)." else "; one chance in two per factor ([] above)."),
		paste0("Point de d\u00e9part, fix\u00e9 par la m\u00e9thode avant de voir les donn\u00e9es : taille attendue r = ", num(obj$priors$rscale),
			if (has_int) " ; une interaction n'entre qu'avec ses facteurs ; une chance sur deux par facteur, les interactions moins probables ([] ci-dessus)." else " ; une chance sur deux par facteur ([] ci-dessus)."))
	say("These are probabilities of hypotheses, not p-values: no threshold of 5 % applies. Whether a difference has a cause depends on how the study was run (randomisation, order, operator).",
		"Ce sont des probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values : aucun seuil de 5 % ne s'y applique. La cause d'une diff\u00e9rence d\u00e9pend de la fa\u00e7on dont l'\u00e9tude a \u00e9t\u00e9 men\u00e9e (tirage au sort, ordre, op\u00e9rateur).")
	if (length(obj$notes$en) && !isTRUE(narrated)) {
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
	if (single) say(paste0("  Chance of a difference before the data: ", .bay_pct(obj$hypotheses$prior[!obj$hypotheses$null], "."), " (prior_h1)."),
		paste0("  Chance d'une diff\u00e9rence avant les donn\u00e9es : ", .bay_pct(obj$hypotheses$prior[!obj$hypotheses$null], ","), " (prior_h1).")) else
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
