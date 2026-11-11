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
	paired <- !is.null(id) && !is.data.frame(id) && length(factors) == 1L && nlevels(factors[[1L]]) == 2L
	safe <- paste0("f", seq_along(factors))
	names(safe) <- names(factors)
	dat <- data.frame(.y = y, factors, check.names = FALSE)
	names(dat) <- c(".y", safe)
	# Random terms: the subject, and with replicated measures subject x condition (D21).
	rand <- NULL
	if (is.data.frame(id)) {
		rand <- paste0(".id", c("", seq_along(id)[-1L]))
		for (j in seq_along(id)) dat[[rand[j]]] <- id[[j]]
	} else if (!is.null(id)) {
		rand <- ".id"
		dat$.id <- id
	}
	rand_txt <- paste(rand, collapse = " + ")
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
				lb0 <- BayesFactor::extractBF(BayesFactor::lmBF(stats::as.formula(paste(".y ~", rand_txt)), data = dat, whichRandom = rand,
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
				fo <- stats::as.formula(paste(".y ~", paste(c(if (has_cov) cov_txt, term_safe[m]), collapse = " + "), if (!is.null(id)) paste("+", rand_txt) else ""))
				x <- tryCatch(suppressMessages(BayesFactor::extractBF(BayesFactor::lmBF(fo, data = dat,
					whichRandom = rand,
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
.bay_factorial <- function(input, rscale, sensitivity, draws, seed, note, prior_h1 = 0.5, het = FALSE) {
	y <- input$response
	id <- input$id
	paired <- !is.null(id) && is.null(input$id_rand) && length(input$factors) == 1L && nlevels(input$factors[[1L]]) == 2L
	idr <- if (!is.null(input$id_rand)) input$id_rand else id
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
	bf <- .bay_fact_logbf(y, factors, sp, rs, iterations, seed, seeds = seed + c(0L, 1L, 2L)[seq_along(rs)], id = idr,
		cov = cov, rconts = .bay_rcont() * rs / rscale)
	# One term alone: its prior is prior_h1, as with one factor between groups.
	single <- length(sp$keys) == 1L
	lp <- if (single) log(c(1 - prior_h1, prior_h1)) else .bay_fact_prior(sp, "hop")
	# Variances per cell (D19): the Bayes factors are averaged over the common variance and one
	# variance per cell (sys_bay_fact_laws.R). If the importance sampling fails, the common
	# variance is kept and the caller caps the support.
	fl <- NULL
	p_het <- NA_real_
	if (isTRUE(het)) {
		# At least 4000 draws, so that the central answer does not move with draws below that.
		fl <- .bay_fl_het_all(y, factors, sp, rs, min(max(4000L, as.integer(draws)), 20000L), seed + 104729L)
		# Heavy tails (D24): Student t with a common scale, same measure, same models.
		if (!fl$failed && "t_eq" %in% .bay_laws_allowed(fl$stats, "t_eq", 1.5)) {
			ft <- .bay_ft_all(y, factors, sp, rs, min(max(4000L, as.integer(draws)), 20000L), seed + 130363L)
			l_t <- matrix(t(vapply(ft, `[[`, numeric(length(rs)), "log_ml")), ncol = length(rs))
			if (!anyNA(l_t) && all(is.finite(l_t)) && all(vapply(ft, `[[`, numeric(1), "ess") >= .bay_law_min_ess())) fl$l_t <- l_t
		}
		if (!fl$failed) {
			cb <- .bay_fl_combine(bf$log_bf, fl$eq0, fl$l_het, fl$l_t)
			p_laws <- .bay_fl_p_laws(c(list(cb$l_eq[, 1L], fl$l_het[, 1L]), if (!is.null(fl$l_t)) list(fl$l_t[, 1L])), lp)
			p_het <- p_laws[2L]
			p_t <- if (!is.null(fl$l_t)) p_laws[3L] else NA_real_
			bf$log_bf <- cb$log_bf
		}
	}
	if (!exists("p_t", inherits = FALSE)) p_t <- NA_real_
	het_used <- !is.null(fl) && !fl$failed
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
		bf2 <- .bay_fact_logbf(y, factors, sp, rscale, iterations, seed, seeds = seed + 7919L, id = idr,
			cov = cov, rconts = .bay_rcont())
		if (het_used) bf2$log_bf <- .bay_fl_combine(bf2$log_bf, fl$eq0, fl$l_het[, 1L, drop = FALSE],
			if (!is.null(fl$l_t)) fl$l_t[, 1L, drop = FALSE])$log_bf
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
			# Subjects per level (D21): n alone would read as n subjects with replicates.
			n_subj = if (!is.null(id)) as.integer(tapply(id, factors[[v]], function(x) length(unique(x)))) else NA_integer_,
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
		cov_centres = if (length(cov_names)) vapply(input$covariates, mean, numeric(1)) else NULL,
		het_used = het_used, het_failed = !is.null(fl) && fl$failed, p_het = p_het, p_t = p_t,
		het_ess = if (!is.null(fl)) min(fl$ess) else NA_real_)
}

# Answer per term, one rule for the bilan and the report sentence (audit of 23c941e, point 1):
# yes or no with its support, or "not settled" when the support is weak, the Bayes factor goes
# against the leading side, or the sensitivity reverses it.
.bay_term_answers <- function(obj) {
	e <- obj$effects
	sens <- obj$sensitivity
	out <- lapply(seq_len(nrow(e)), function(i) {
		yes <- e$p_h1[i] >= 0.5
		s <- .bay_support_bf(max(e$p_h1[i], e$p_h0[i]), e$log_bf[i] * if (yes) 1 else -1)
		if ((isTRUE(any(obj$overlap$poor)) || (isTRUE(obj$het_cells) && yes)) && s$level > 2L) s <- .bay_support(0.7)
		against <- !is.na(e$log_bf[i]) && ((e$log_bf[i] >= 0) != yes)
		open <- s$level == 1L || against || (!is.null(sens) && identical(sens$verdict_terms[i], "changes"))
		data.frame(open = open, yes = yes,
			short_en = if (open) "not settled" else if (yes) "yes" else "no",
			short_fr = if (open) "non tranch\u00e9" else if (yes) "oui" else "non",
			long_en = if (open) "not settled" else paste0(if (yes) "yes, " else "no, ", s$en),
			long_fr = if (open) "non tranch\u00e9" else paste0(if (yes) "oui, " else "non, ", s$fr),
			stringsAsFactors = FALSE)
	})
	do.call(rbind, out)
}

# Bilan with several factors. Same order as with one factor: verdict on the global
# question, then each term, its sensitivity, what was compared, the starting point.
.bay_print_bilan_fact <- function(obj, detail = FALSE, narrated = FALSE) {
	dec <- if (.bay_lang_fr()) "," else "."
	num <- function(x) gsub(".", dec, format(signif(x, 3), trim = TRUE), fixed = TRUE)
	e <- obj$effects
	sens <- obj$sensitivity
	anc <- identical(obj$plan, "ancova")
	single <- nrow(e) == 1L
	body <- function(emit, li, narrated) {
		ph <- obj$posthoc
		if (identical(ph$method, "slopes")) {
			gp <- ph$gaps
			flips <- any(gp$diff_high < 0) && any(gp$diff_low > 0)
			emit("fact", paste0("The gap between groups depends on ", ph$xname, if (flips) ": it changes sign, no group is ahead everywhere." else ": its size changes along it."),
				paste0("L'\u00e9cart entre groupes d\u00e9pend de ", ph$xname, if (flips) " : il change de sens, aucun groupe n'est devant partout." else " : son ampleur change selon sa valeur."))
			emit("note", if (narrated) paste0("gaps along ", ph$xname, ": step 6.") else "gaps: $posthoc$gaps.",
				if (narrated) paste0("\u00e9carts selon ", ph$xname, " : \u00e9tape 6.") else "\u00e9carts : $posthoc$gaps.")
		}
		if (nrow(e) > 1L) {
			ta <- .bay_term_answers(obj)
			for (i in seq_len(nrow(e))) {
				ans_en <- ta$long_en[i]
				ans_fr <- ta$long_fr[i]
				q <- strsplit(e$term[i], ":", fixed = TRUE)[[1]]
				q_en <- if (length(q) == 1L) paste0(e$term[i], " (alone or in interaction)") else paste0(e$term[i], " (the gaps along ", paste(q[-length(q)], collapse = ":"), " change with ", q[length(q)], ")")
				q_fr <- if (length(q) == 1L) paste0(e$term[i], " (seul ou en interaction)") else paste0(e$term[i], " (les \u00e9carts selon ", paste(q[-length(q)], collapse = ":"), " changent avec ", q[length(q)], ")")
				emit("fact", paste0(q_en, ": ", ans_en, " (", .bay_pct(e$p_h1[i], "."), "; prior ", .bay_pct(e$prior_h1[i], "."), "; Bayes factor ", .bay_bf_txt(e$log_bf[i], "."), ")."),
					paste0(q_fr, " : ", ans_fr, " (", .bay_pct(e$p_h1[i], ","), " ; a priori ", .bay_pct(e$prior_h1[i], ","), " ; facteur de Bayes ", .bay_bf_txt(e$log_bf[i], ","), ")."))
			}
		}
		if (nrow(e) > 1L) {
			odd <- which(!is.na(e$log_bf) & e$log_bf > log(1.5) & e$p_h1 < 0.5)
			if (length(odd)) {
				i <- odd[1L]
				emit("note", paste0("the data raise the odds of ", e$term[i], " (Bayes factor ", .bay_bf_txt(e$log_bf[i], "."), "), but it started at ", .bay_pct(e$prior_h1[i], "."),
						", hence ", .bay_pct(e$p_h1[i], "."), " only."),
					paste0("les donn\u00e9es augmentent les chances de ", e$term[i], " (facteur de Bayes ", .bay_bf_txt(e$log_bf[i], ","), "), mais ce terme partait \u00e0 ", .bay_pct(e$prior_h1[i], ","),
						", d'o\u00f9 ", .bay_pct(e$p_h1[i], ","), " seulement."))
			}
		}
		if (identical(ph$method, "factorial")) {
			for (j in seq_along(ph$families)) {
				f <- ph$families[[j]]
				lab_en <- if (f$cells) paste0("cells ", gsub(":", " x ", f$term, fixed = TRUE)) else f$term
				lab_fr <- if (f$cells) paste0("cases ", gsub(":", " x ", f$term, fixed = TRUE)) else f$term
				if (f$k > 2L) {
					emit("fact", paste0("Pairs, ", lab_en, ":"), paste0("Paires, ", lab_fr, " :"))
					.bay_pairs_short(f, emit, FALSE, narrated, note = j == length(ph$families))
				} else {
					pr <- f$pairs
					dg <- if (!is.null(f$dg)) f$dg else .bay_dg(obj$descriptive_stats$mean)
					ci <- function(dc) paste0(.bay_fx(pr$diff, dg, dc), " [", .bay_fx(pr$diff_low, dg, dc), if (dc == ",") " ; " else "; ", .bay_fx(pr$diff_high, dg, dc), "]")
					emit("fact", paste0("Estimated gap ", pr$group2, " - ", pr$group1, " (", lab_en, "): ", ci("."), if (narrated) " (step 6)." else "."),
						paste0("\u00c9cart estim\u00e9 ", pr$group2, " - ", pr$group1, " (", lab_fr, ") : ", ci(","), if (narrated) " (\u00e9tape 6)." else "."))
				}
			}
		} else if (!is.null(ph) && !identical(ph$method, "slopes")) {
			if (ph$k > 2L) .bay_pairs_short(ph, emit, li$h$posterior[1L] < 0.6 || li$open, narrated)
			if (ph$k == 2L) {
				pr <- ph$pairs
				dg <- .bay_dg(obj$descriptive_stats$mean)
				ci <- function(dc) paste0(.bay_fx(pr$diff, dg, dc), " [", .bay_fx(pr$diff_low, dg, dc), if (dc == ",") " ; " else "; ", .bay_fx(pr$diff_high, dg, dc), "]")
				emit("fact", paste0("Estimated gap ", pr$group2, " - ", pr$group1, ": ", ci("."), if (narrated) " (step 6)." else "."),
					paste0("\u00c9cart estim\u00e9 ", pr$group2, " - ", pr$group1, " : ", ci(","), if (narrated) " (\u00e9tape 6)." else "."))
			}
		}
		emit("block", "", "")
		m <- obj$descriptive_stats
		dg <- .bay_dg(m$mean)
		for (f in unique(m$factor)) {
			mf <- m[m$factor == f, , drop = FALSE]
			same_n <- length(unique(mf$n)) == 1L
			# With several measures per subject at a level, say subjects x measures (D21).
			per <- mf$n / mf$n_subj
			multi <- !anyNA(per) && any(per > 1) && length(unique(mf$n_subj)) == 1L && all(per == round(per)) && same_n
			nt <- function(dc) {
				if (!multi) return(paste0("n = ", mf$n[1L]))
				rep_w <- identical(obj$plan, "mixed") && isTRUE(mf$n_subj[1L] == obj$n_subjects)
				if (dc == ",") paste0(mf$n_subj[1L], " sujets x ", per[1L], if (rep_w) " r\u00e9pliques" else " mesures")
				else paste0(mf$n_subj[1L], " subjects x ", per[1L], if (rep_w) " replicates" else " measures")
			}
			txt <- function(dc) paste0(paste0(mf$level, " = ", .bay_fx(mf$mean, dg, dc), if (!same_n) paste0(" (n = ", mf$n, ")") else "", collapse = if (dc == ",") " ; " else "; "),
				if (same_n) paste0(" (", nt(dc), ")") else "")
			emit("fact", paste0(if (anc) "Raw means, " else "Observed means, ", f, ": ", txt("."), "."),
				paste0(if (anc) "Moyennes brutes, " else "Moyennes observ\u00e9es, ", f, " : ", txt(","), "."))
		}
		cm <- obj$cell_means
		if (!is.null(cm) && any(grepl(":", e$term, fixed = TRUE)) && !anc && nrow(cm) <= 16L) {
			emit("fact", paste0("Per cell: ", paste0(cm$cell, " = ", .bay_fx(cm$mean, dg, "."), collapse = "; "), "."),
				paste0("Par case : ", paste0(cm$cell, " = ", .bay_fx(cm$mean, dg, ","), collapse = " ; "), "."))
		}
		if (anc) {
			cg <- obj$cov_by_group
			for (v in colnames(cg)) {
				emit("fact", paste0("Mean ", v, " per group: ", paste0(rownames(cg), " = ", format(signif(cg[, v], 3), trim = TRUE), collapse = "; "), "."),
					paste0(v, " moyen par groupe : ", paste0(rownames(cg), " = ", gsub(".", ",", format(signif(cg[, v], 3), trim = TRUE), fixed = TRUE), collapse = " ; "), "."))
			}
			cv <- paste(colnames(cg), collapse = ", ")
			emit("fact", paste0(cv, " is in every model, so the groups are compared at equal ", cv, "; its own effect is not the question asked."),
				paste0(cv, " est dans tous les mod\u00e8les, donc les groupes sont compar\u00e9s \u00e0 ", cv, " \u00e9gal ; son effet propre n'est pas la question pos\u00e9e."))
			if (isTRUE(any(obj$overlap$poor))) {
				emit("decide", "support capped at weak, since the groups hardly overlap on the covariate and the adjustment extrapolates.",
					"soutien plafonn\u00e9 \u00e0 faible, car les groupes se recouvrent peu sur la covariable et l'ajustement extrapole.")
			}
		}
		if (isTRUE(obj$het_cells)) {
			emit("decide", "Support capped at weak: unequal variances per cell, which the common-variance model does not take into account.",
				"Soutien plafonn\u00e9 \u00e0 faible : variances in\u00e9gales par case, que le mod\u00e8le \u00e0 variance commune ne prend pas en compte.")
		}
		rm_ <- obj$plan %in% c("paired", "repeated", "mixed")
		mx_ <- identical(obj$plan, "mixed")
		sd_en <- if (anc) "Residual standard deviation (covariates accounted for): " else if (mx_) "Standard deviation of one measure of a subject in a condition (replicates and subject x condition): " else if (rm_) "Residual standard deviation within subject: " else if (isTRUE(obj$sd_is_within)) "Within-cell standard deviation: " else "Residual standard deviation: "
		sd_fr <- if (anc) "\u00c9cart-type r\u00e9siduel (covariables prises en compte) : " else if (mx_) "\u00c9cart-type d'une mesure d'un sujet dans une condition (r\u00e9pliques et sujet x condition) : " else if (rm_) "\u00c9cart-type r\u00e9siduel au sein des sujets : " else if (isTRUE(obj$sd_is_within)) "\u00c9cart-type au sein des cases : " else "\u00c9cart-type r\u00e9siduel : "
		emit("fact", paste0(sd_en, .bay_fx(obj$sd_within, dg + 1L, "."), if (isTRUE(obj$het_used) && isTRUE(!is.na(obj$p_t))) "; normal (common or per cell) and Student t, averaged." else if (isTRUE(obj$het_used)) "; common variance or one per cell, averaged." else "; normal law with a common variance only."),
			paste0(sd_fr, .bay_fx(obj$sd_within, dg + 1L, ","), if (isTRUE(obj$het_used) && isTRUE(!is.na(obj$p_t))) " ; normale (commune ou par case) et Student, moyenn\u00e9es." else if (isTRUE(obj$het_used)) " ; variance commune ou par case, moyenn\u00e9es." else " ; loi normale \u00e0 variance commune seule."))
		if (!is.null(ph) && isTRUE(ph$caution_overlap)) {
			emit("decide", "the adjusted gaps compare the groups where some have few values of the covariate.",
				"les \u00e9carts ajust\u00e9s comparent les groupes l\u00e0 o\u00f9 certains ont peu de valeurs de la covariable.")
		}
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)

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
