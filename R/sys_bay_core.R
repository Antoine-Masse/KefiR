# sys_bay_core.R
# Internal, engine-independent core of the Bayesian family (m.bay, later
# cor.bay and g.bay). Nothing here is exported and nothing calls a sampler:
# every function takes Bayes factors already computed by some engine.
# Design notes: dev/mbay/METHODOLOGIE.md, sections 2, 3, 4, 8 and 10.
#
# References (comments only, never in messages: dev/bp.log, 7.4 and BP-008).
#   Kass, R. E., & Raftery, A. E. (1995). Bayes factors. JASA, 90(430), 773-795.
#     https://doi.org/10.1080/01621459.1995.10476572
#   Morey, R. D., Romeijn, J.-W., & Rouder, J. N. (2016). The philosophy of
#     Bayes factors and the quantification of statistical evidence. Journal of
#     Mathematical Psychology, 72, 6-18.
#     https://doi.org/10.1016/j.jmp.2015.11.001
#   Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009).
#     Bayesian t tests for accepting and rejecting the null hypothesis.
#     Psychonomic Bulletin & Review, 16(2), 225-237.
#     https://doi.org/10.3758/PBR.16.2.225
#   Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012).
#     Default Bayes factors for ANOVA designs. Journal of Mathematical
#     Psychology, 56(5), 356-374. https://doi.org/10.1016/j.jmp.2012.08.001
#   Nelder, J. A. (1977). A reformulation of linear models. JRSS A, 140(1), 48.
#     https://doi.org/10.2307/2344517
#   Peixoto, J. L. (1987). Hierarchical variable selection in polynomial
#     regression models. The American Statistician, 41, 311-313.
#     https://doi.org/10.1080/00031305.1987.10475506
#   Chipman, H. (1996). Bayesian variable selection with related predictors.
#     Canadian Journal of Statistics, 24(1), 17-36.
#     https://doi.org/10.2307/3315687
#   Scott, J. G., & Berger, J. O. (2010). Bayes and empirical-Bayes multiplicity
#     adjustment in the variable-selection problem. Annals of Statistics, 38(5),
#     2587-2619. https://doi.org/10.1214/10-AOS792
#   Makowski, D., Ben-Shachar, M. S., Chen, S. H. A., & Luedecke, D. (2019).
#     Indices of effect existence and significance in the Bayesian framework.
#     Frontiers in Psychology, 10, 2767.
#     https://doi.org/10.3389/fpsyg.2019.02767
#   Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian
#     estimation. AMPPS, 1(2), 270-280.
#     https://doi.org/10.1177/2515245918771304
#   Mathot, S. (2017). Bayes like a Baws: interpreting Bayesian repeated
#     measures in JASP. Blog post, no DOI.
#     https://www.cogsci.nl/blog/interpreting-bayesian-repeated-measures-in-jasp


# ---------------------------------------------------------------------------
# Small internal helpers
# ---------------------------------------------------------------------------

# Numeric coercion that lets an all-NA logical vector through (bf = NA) and
# turns NaN into NA, so that callers only ever have to deal with one "missing".
.bay_as_num <- function(x, what) {
	if (is.logical(x) && all(is.na(x))) x <- as.numeric(x)
	if (!is.numeric(x)) {
		.exit(paste0(what, " must be numeric."),
			paste0(what, " doit \u00eatre num\u00e9rique."))
	}
	x_names <- names(x)
	x <- as.numeric(x)
	x[is.nan(x)] <- NA_real_
	names(x) <- x_names
	x
}

# log(sum(exp(x))) without overflow. Empty sum: -Inf. Any NA: NA.
.bay_log_sum_exp <- function(x) {
	if (length(x) == 0L) return(-Inf)
	if (anyNA(x)) return(NA_real_)
	top <- max(x)
	if (!is.finite(top)) return(top)
	top + log(sum(exp(x - top)))
}

# "1:4", "1:1", "4:1": odds printed as H1:H0.
.bay_odds_label <- function(prior_odds) {
	vapply(prior_odds, function(o) {
		if (isTRUE(o >= 1)) {
			paste0(format(signif(o, 3), trim = TRUE), ":1")
		} else {
			paste0("1:", format(signif(1 / o, 3), trim = TRUE))
		}
	}, character(1), USE.NAMES = FALSE)
}

.bay_check_odds <- function(prior_odds) {
	prior_odds <- .bay_as_num(prior_odds, "prior_odds")
	if (length(prior_odds) == 0L || anyNA(prior_odds) ||
		any(prior_odds <= 0) || any(is.infinite(prior_odds))) {
		.exit("prior_odds must be finite and strictly positive.",
			"prior_odds doit \u00eatre fini et strictement positif.")
	}
	prior_odds
}


# ---------------------------------------------------------------------------
# P(H1 | D) from a Bayes factor
# ---------------------------------------------------------------------------

# P(H1|D) = O*BF / (1 + O*BF), with O = pi / (1 - pi).
# Computed as plogis(log(O) + log(BF)), which neither overflows nor loses the
# small tail: BF = exp(1000) would give Inf/Inf on the natural scale. Give
# either bf or log_bf. The small complement P(H0|D) keeps its precision when
# requested as .bay_post_prob(log_bf = -log_bf, prior_odds = 1 / O), which is
# what .bay_effects_row does rather than computing 1 - p_h1.
# Kass & Raftery (1995); Morey et al. (2016): no p without pi.
.bay_post_prob <- function(bf = NULL, prior_odds = 1, log_bf = NULL) {
	if (is.null(bf) == is.null(log_bf)) {
		.exit("Give exactly one of bf and log_bf.",
			"Fournir exactement un argument parmi bf et log_bf.")
	}
	if (!is.null(bf)) {
		bf <- .bay_as_num(bf, "bf")
		if (any(bf < 0, na.rm = TRUE)) {
			.exit("A Bayes factor cannot be negative.",
				"Un facteur de Bayes ne peut pas \u00eatre n\u00e9gatif.")
		}
		log_bf <- log(bf)
	} else {
		log_bf <- .bay_as_num(log_bf, "log_bf")
	}
	prior_odds <- .bay_check_odds(prior_odds)
	if (length(log_bf) == 0L) return(numeric(0))
	n <- max(length(log_bf), length(prior_odds))
	if ((length(prior_odds) != 1L && length(prior_odds) != n) ||
		(length(log_bf) != 1L && length(log_bf) != n)) {
		.exit("Bayes factors and prior_odds have incompatible lengths.",
			"Les facteurs de Bayes et prior_odds ont des longueurs incompatibles.")
	}
	res <- stats::plogis(log(prior_odds) + log_bf)
	if (length(log_bf) == n) names(res) <- names(log_bf)
	res
}


# ---------------------------------------------------------------------------
# ROPE
# ---------------------------------------------------------------------------

# Prior mass of |theta| < w for a centred Cauchy of scale r (under H1), and the
# total prior mass of the ROPE once the point null is counted:
#   P(|theta| < w | H1) = (2/pi) * atan(w/r)
#   P(ROPE) = (1 - prior_h1) + prior_h1 * P(|theta| < w | H1)
# r = 0.707, w = 0.1, prior_h1 = 1/2 gives 0.0894 and 0.5447 (METHODOLOGIE, 4).
# IMPORTANT: one-dimensional parameter only (delta, a slope). For a term with
# several df the ROPE bears on f and this closed form does not apply as such.
# Rouder et al. (2009); Kruschke (2018); Makowski et al. (2019).
.bay_prior_mass_rope <- function(r, w, prior_h1 = 1/2) {
	r <- .bay_as_num(r, "r")
	w <- .bay_as_num(w, "w")
	prior_h1 <- .bay_as_num(prior_h1, "prior_h1")
	if (anyNA(r) || any(r <= 0) || any(is.infinite(r))) {
		.exit("r must be finite and strictly positive.",
			"r doit \u00eatre fini et strictement positif.")
	}
	if (anyNA(w) || any(w < 0)) {
		.exit("w must be positive or zero.",
			"w doit \u00eatre positif ou nul.")
	}
	if (anyNA(prior_h1) || any(prior_h1 < 0) || any(prior_h1 > 1)) {
		.exit("prior_h1 must lie between 0 and 1.",
			"prior_h1 doit \u00eatre compris entre 0 et 1.")
	}
	mass_h1 <- 2 / pi * atan(w / r)
	data.frame(rscale = r, rope_width = w, prior_h1 = prior_h1,
		mass_h1 = mass_h1, mass_total = (1 - prior_h1) + prior_h1 * mass_h1)
}

# P(ROPE|D) = P(H0|D) + P(H1|D) * P(|theta| < w | D, H1).
# The conditional comes from the whole posterior (ci = 1), never from a 95 %
# interval. When P(H1|D) = 0 the conditional is irrelevant: the result is 1 even
# if it is missing.
.bay_rope_prob <- function(p_h1, p_in_rope_given_h1) {
	p_h1 <- .bay_as_num(p_h1, "p_h1")
	p_in <- .bay_as_num(p_in_rope_given_h1, "p_in_rope_given_h1")
	if (any(p_h1 < 0 | p_h1 > 1, na.rm = TRUE) ||
		any(p_in < 0 | p_in > 1, na.rm = TRUE)) {
		.exit("Probabilities must lie between 0 and 1.",
			"Les probabilit\u00e9s doivent \u00eatre comprises entre 0 et 1.")
	}
	n <- max(length(p_h1), length(p_in))
	if (length(p_h1) == 0L || length(p_in) == 0L) return(numeric(0))
	if ((length(p_h1) != 1L && length(p_h1) != n) ||
		(length(p_in) != 1L && length(p_in) != n)) {
		.exit("p_h1 and p_in_rope_given_h1 have incompatible lengths.",
			"p_h1 et p_in_rope_given_h1 ont des longueurs incompatibles.")
	}
	res_names <- if (length(p_h1) == n) names(p_h1) else names(p_in)
	p_h1 <- rep_len(p_h1, n)
	p_in <- rep_len(p_in, n)
	res <- (1 - p_h1) + p_h1 * p_in
	res[which(p_h1 == 0)] <- 1
	names(res) <- res_names
	res
}


# ---------------------------------------------------------------------------
# One row per term: shared output format
# ---------------------------------------------------------------------------

# Columns are reusable by cor.bay and g.bay:
#   term       term name, same spelling as m.test()$effects
#   bf10       Bayes factor H1 against H0 (Inf if exp(log_bf) overflows)
#   log_bf     its log, the quantity to trust for very strong evidence
#   p_h1       P(H1|D), P(term has an effect | data, priors)
#   p_h0       P(H0|D), computed separately to keep the small tail
#   prior_h1   pi used for this term
#   rscale     numeric prior scale r
#   error_pct  numerical error on the BF, in PERCENT. BayesFactor reports a
#              proportion: the caller multiplies by 100.
#   analogue   closest frequentist hypothesis ("II", "III"), or NA. It does not
#              mean that the two answers are convertible (METHODOLOGIE, 3).
# On purpose, no column is a p-value and none is called p.value: p_h1 and p_h0
# are posterior probabilities of hypotheses.
.bay_effects_row <- function(term, log_bf = NULL, bf10 = NULL, prior_h1 = 1/2,
                             rscale = NA_real_, error_pct = NA_real_,
                             analogue = NA_character_) {
	if (is.null(bf10) == is.null(log_bf)) {
		.exit("Give exactly one of bf10 and log_bf.",
			"Fournir exactement un argument parmi bf10 et log_bf.")
	}
	if (!is.null(bf10)) {
		bf10 <- .bay_as_num(bf10, "bf10")
		if (any(bf10 < 0, na.rm = TRUE)) {
			.exit("A Bayes factor cannot be negative.",
				"Un facteur de Bayes ne peut pas \u00eatre n\u00e9gatif.")
		}
		log_bf <- log(bf10)
	}
	log_bf <- unname(.bay_as_num(log_bf, "log_bf"))
	term <- as.character(term)
	n <- length(term)
	if (n == 0L) {
		.exit("At least one term is required.",
			"Au moins un terme est requis.")
	}
	if (length(log_bf) != n) {
		.exit("term and the Bayes factors must have the same length.",
			"term et les facteurs de Bayes doivent avoir la m\u00eame longueur.")
	}
	prior_h1 <- unname(.bay_as_num(prior_h1, "prior_h1"))
	rscale <- unname(.bay_as_num(rscale, "rscale"))
	error_pct <- unname(.bay_as_num(error_pct, "error_pct"))
	if (anyNA(prior_h1) || any(prior_h1 < 0) || any(prior_h1 > 1)) {
		.exit("prior_h1 must lie between 0 and 1.",
			"prior_h1 doit \u00eatre compris entre 0 et 1.")
	}
	if (any(rscale <= 0, na.rm = TRUE) || any(error_pct < 0, na.rm = TRUE)) {
		.exit("rscale must be strictly positive and error_pct positive or zero.",
			"rscale doit \u00eatre strictement positif et error_pct positif ou nul.")
	}
	for (v in list(prior_h1, rscale, error_pct, analogue)) {
		if (length(v) != 1L && length(v) != n) {
			.exit("Each column must have length 1 or the number of terms.",
				"Chaque colonne doit \u00eatre de longueur 1 ou le nombre de termes.")
		}
	}
	prior_h1 <- rep_len(prior_h1, n)

	# pi = 0 or 1 leaves nothing to update: the data cannot change it.
	p_h1 <- prior_h1
	p_h0 <- 1 - prior_h1
	open <- prior_h1 > 0 & prior_h1 < 1
	if (any(open)) {
		odds <- prior_h1[open] / (1 - prior_h1[open])
		p_h1[open] <- .bay_post_prob(log_bf = log_bf[open], prior_odds = odds)
		p_h0[open] <- .bay_post_prob(log_bf = -log_bf[open], prior_odds = 1 / odds)
	}

	data.frame(
		term = term, bf10 = exp(log_bf), log_bf = log_bf,
		p_h1 = p_h1, p_h0 = p_h0, prior_h1 = prior_h1,
		rscale = rep_len(rscale, n), error_pct = rep_len(error_pct, n),
		analogue = rep_len(as.character(analogue), n),
		stringsAsFactors = FALSE, row.names = NULL
	)
}


# ---------------------------------------------------------------------------
# Inclusion Bayes factors
# ---------------------------------------------------------------------------

# Canonical form of a term, so that "B:A" and "A:B" are one term.
.bay_term_key <- function(term) {
	vapply(strsplit(term, ":", fixed = TRUE), function(f) {
		paste(sort(trimws(f), method = "radix"), collapse = ":")
	}, character(1), USE.NAMES = FALSE)
}

# Model table to a list of character vectors. Accepts a list of term vectors
# or strings such as "A + B + A:B". Null model: character(0), NULL, "1" or "".
.bay_parse_terms <- function(terms) {
	if (is.character(terms)) {
		terms <- strsplit(terms, "+", fixed = TRUE)
	}
	if (!is.list(terms)) {
		.exit("terms must be a list of term vectors, or a character vector.",
			"terms doit \u00eatre une liste de vecteurs de termes, ou du texte.")
	}
	lapply(terms, function(m) {
		m <- trimws(as.character(m))
		m <- unique(m[!m %in% c("", "1")])
		if (any(grepl("*", m, fixed = TRUE))) {
			.exit("Write interactions out in terms (A + B + A:B, not A*B).",
				"\u00c9crire les interactions en toutes lettres (A + B + A:B, pas A*B).")
		}
		m
	})
}

# Marginality: a term needs all the non-empty proper subsets of its factors in
# the model, unless they belong to the nuisance part present everywhere.
# Nelder (1977); Peixoto (1987).
.bay_violates_marginality <- function(model, nuisance) {
	available <- c(model, nuisance)
	for (tt in model) {
		f <- strsplit(tt, ":", fixed = TRUE)[[1]]
		if (length(f) < 2L) next
		for (k in seq_len(length(f) - 1L)) {
			parents <- utils::combn(f, k, paste, collapse = ":")
			if (!all(parents %in% available)) return(TRUE)
		}
	}
	FALSE
}

# Inclusion Bayes factor of each term, from a table of models.
#
# log_bf   : log BF of each model against ONE common reference model. The
#            reference does not matter (adding a constant changes nothing).
# terms    : terms of each model, null model included (see .bay_parse_terms).
# type     : "matched" (default, METHODOLOGIE 3) or "all".
# nuisance : terms present in every model and not tested, such as a subject id,
#            a covariate or a random slope "A:id". They are set aside first, so
#            A + id + A:id against id + A:id is the comparison for A and a slope
#            A:id does not count as a higher-order term of A (METHODOLOGIE
#            5.4 and 5.8).
# prior_prob : optional law on the models. Default: pi = 1/2 per term in
#            "matched", uniform on the models in "all".
#
# Definitions, for a term T:
#   "matched": M+ = models with T and without any higher-order term containing
#              T; M- = those same models with T removed. In a full marginal
#              space the two sets are in one-to-one correspondence, and the
#              comparison is the Type II one. Mathot (2017).
#   "all"    : M+ = all models with T, M- = all models without T. These are the
#              marginals of a single law on the models, but with a uniform law
#              they do NOT start at 1/2 (0.6 and 0.2 with 2 factors). Chipman
#              (1996) for a coherent alternative, Scott & Berger (2010) on
#              multiplicity.
#   BF_T = [sum_M+ pi(m) p(D|m) / sum_M+ pi(m)] /
#          [sum_M- pi(m) p(D|m) / sum_M- pi(m)]
#
# IMPORTANT: matched probabilities come from different subspaces, so they are
# not the marginals of one joint law and must not be combined across terms.
#
# Refusals and signals:
#   - "matched" refuses a space with a model violating marginality. Otherwise
#     M+ for A:B would mix in A:B, A + A:B and B + A:B against A + B alone,
#     which is not the comparison announced. "all" computes on the space given,
#     warns, and sets attr(res, "marginal_ok") to FALSE.
#   - a term in every model and not declared as nuisance is refused.
#   - a missing partner model in "matched" is dropped with a warning and counted
#     in n_unpaired; the BF would no longer compare like with like.
#   - a used model with an NA log BF gives NA for the terms that use it, rather
#     than a value silently computed without it.
.bay_inclusion <- function(log_bf, terms, type = c("matched", "all"),
                           prior_prob = NULL, nuisance = character(0)) {
	type <- match.arg(type)
	log_bf <- unname(.bay_as_num(log_bf, "log_bf"))
	models <- .bay_parse_terms(terms)
	n <- length(models)
	if (n != length(log_bf)) {
		.exit("log_bf and terms must describe the same models.",
			"log_bf et terms doivent d\u00e9crire les m\u00eames mod\u00e8les.")
	}

	if (is.null(prior_prob)) {
		weights <- rep(1 / n, n)
	} else {
		weights <- unname(.bay_as_num(prior_prob, "prior_prob"))
		if (length(weights) != n || anyNA(weights) || any(weights < 0) ||
			any(is.infinite(weights)) || sum(weights) <= 0) {
			.exit("prior_prob must give one finite weight >= 0 per model, not all 0.",
				"prior_prob doit donner un poids fini >= 0 par mod\u00e8le, pas tous nuls.")
		}
		weights <- weights / sum(weights)
	}

	# First spelling met is kept for display: same names as m.test()$effects.
	all_terms <- unique(unlist(models, use.names = FALSE))
	labels <- all_terms[!duplicated(.bay_term_key(all_terms))]
	names(labels) <- .bay_term_key(labels)
	models <- lapply(models, function(m) unique(.bay_term_key(m)))

	nuisance <- if (length(nuisance)) .bay_term_key(trimws(nuisance)) else character(0)
	in_all <- vapply(nuisance, function(tt) {
		all(vapply(models, function(m) tt %in% m, logical(1)))
	}, logical(1))
	if (!all(in_all)) {
		.exit(paste0("Nuisance terms must be in every model: ",
				paste(nuisance[!in_all], collapse = ", "), "."),
			paste0("Les termes de nuisance doivent figurer dans tous les mod\u00e8les : ",
				paste(nuisance[!in_all], collapse = ", "), "."))
	}
	models <- lapply(models, function(m) sort(setdiff(m, nuisance), method = "radix"))

	model_keys <- vapply(models, paste, character(1), collapse = "+")
	if (anyDuplicated(model_keys)) {
		.exit("The same model appears twice in the table.",
			"Le m\u00eame mod\u00e8le figure deux fois dans le tableau.")
	}
	if (!"" %in% model_keys) {
		.exit("The null model (nuisance terms only) must be in the table.",
			"Le mod\u00e8le nul (termes de nuisance seuls) doit figurer dans le tableau.")
	}

	tested <- unique(unlist(models, use.names = FALSE))
	if (length(tested) == 0L) {
		.exit("No term to test.", "Aucun terme \u00e0 tester.")
	}
	everywhere <- vapply(tested, function(tt) {
		all(vapply(models, function(m) tt %in% m, logical(1)))
	}, logical(1))
	if (any(everywhere)) {
		.exit(paste0("Term in every model, declare it as nuisance: ",
				paste(tested[everywhere], collapse = ", "), "."),
			paste0("Terme pr\u00e9sent dans tous les mod\u00e8les, le d\u00e9clarer en nuisance : ",
				paste(tested[everywhere], collapse = ", "), "."))
	}

	violates <- vapply(models, .bay_violates_marginality, logical(1),
		nuisance = nuisance)
	if (any(violates)) {
		shown <- paste(model_keys[violates], collapse = ", ")
		if (type == "matched") {
			.exit(paste0("Matched inclusion refused, marginality is violated by: ",
					shown, "."),
				paste0("Inclusion appari\u00e9e refus\u00e9e, marginalit\u00e9 viol\u00e9e par : ",
					shown, "."))
		}
		warning(.msg(paste0("Marginality is violated by: ", shown, "."),
				paste0("Marginalit\u00e9 viol\u00e9e par : ", shown, ".")),
			call. = FALSE)
	}

	# Display order: main effects, then increasing order of interaction.
	term_order <- lengths(strsplit(tested, ":", fixed = TRUE))
	tested <- tested[order(term_order)]

	# Log of the prior-weighted mean of p(D|m) on a set of models. Zero-weight
	# models are out of the prior, so they cannot contribute even if BF = Inf.
	log_mean <- function(idx) {
		idx <- idx[weights[idx] > 0]
		if (length(idx) == 0L) return(NA_real_)
		.bay_log_sum_exp(log(weights[idx]) + log_bf[idx]) -
			.bay_log_sum_exp(log(weights[idx]))
	}

	n_terms <- length(tested)
	res_log_bf <- numeric(n_terms)
	res_prior <- numeric(n_terms)
	n_plus <- integer(n_terms)
	n_minus <- integer(n_terms)
	n_unpaired <- integer(n_terms)

	for (i in seq_len(n_terms)) {
		tt <- tested[i]
		has_t <- vapply(models, function(m) tt %in% m, logical(1))
		if (type == "all") {
			idx_plus <- which(has_t)
			idx_minus <- which(!has_t)
		} else {
			f <- strsplit(tt, ":", fixed = TRUE)[[1]]
			higher <- Filter(function(u) {
				fu <- strsplit(u, ":", fixed = TRUE)[[1]]
				length(fu) > length(f) && all(f %in% fu)
			}, tested)
			parents <- character(0)
			for (k in seq_len(length(f) - 1L)) {
				parents <- c(parents, utils::combn(f, k, paste, collapse = ":"))
			}
			no_higher <- vapply(models, function(m) !any(higher %in% m), logical(1))
			has_parents <- vapply(models, function(m) {
				all(parents %in% c(m, nuisance))
			}, logical(1))
			cand_plus <- which(has_t & no_higher)
			cand_minus <- which(!has_t & no_higher & has_parents)
			# Pair by the model obtained once T is removed.
			key_plus <- vapply(models[cand_plus], function(m) {
				paste(setdiff(m, tt), collapse = "+")
			}, character(1))
			key_minus <- model_keys[cand_minus]
			keep_plus <- key_plus %in% key_minus
			keep_minus <- key_minus %in% key_plus
			idx_plus <- cand_plus[keep_plus]
			idx_minus <- cand_minus[keep_minus]
			n_unpaired[i] <- sum(!keep_plus) + sum(!keep_minus)
		}
		n_plus[i] <- length(idx_plus)
		n_minus[i] <- length(idx_minus)
		mass_plus <- sum(weights[idx_plus])
		mass_minus <- sum(weights[idx_minus])

		value <- log_mean(idx_plus) - log_mean(idx_minus)
		# Inf - Inf or -Inf - (-Inf): the ratio is undetermined.
		res_log_bf[i] <- if (is.nan(value)) NA_real_ else value
		res_prior[i] <- if (type == "matched" && is.null(prior_prob)) {
			1 / 2
		} else if (mass_plus + mass_minus > 0) {
			mass_plus / (mass_plus + mass_minus)
		} else {
			NA_real_
		}
	}

	if (any(n_unpaired > 0L)) {
		shown <- paste(labels[tested[n_unpaired > 0L]], collapse = ", ")
		warning(.msg(paste0("Models without a partner were left out for: ", shown, "."),
				paste0("Mod\u00e8les sans partenaire \u00e9cart\u00e9s pour : ", shown, ".")),
			call. = FALSE)
	}

	# No model of the pair has a prior weight: nothing can be said on this term.
	no_prior <- is.na(res_prior)
	res <- .bay_effects_row(
		term = unname(labels[tested]), log_bf = res_log_bf,
		prior_h1 = ifelse(no_prior, 1 / 2, res_prior),
		analogue = if (type == "matched") "II" else NA_character_
	)
	res$prior_h1[no_prior] <- NA_real_
	res$p_h1[no_prior] <- NA_real_
	res$p_h0[no_prior] <- NA_real_
	res$n_plus <- n_plus
	res$n_minus <- n_minus
	res$n_unpaired <- n_unpaired
	attr(res, "type") <- type
	attr(res, "nuisance") <- nuisance
	attr(res, "marginal_ok") <- !any(violates)
	res
}


# ---------------------------------------------------------------------------
# Sensitivity grid
# ---------------------------------------------------------------------------

# Grid r x prior odds of P(H1|D), with its range and a stability flag.
# log_bf_by_r: log BF at each r, named by r if possible (METHODOLOGIE, 8).
# stable is TRUE only if every cell is known and strictly on the same side of
# 1/2. A cell at NA leaves stability unproved, hence FALSE, and is counted in
# n_missing; range is then that of the known cells.
.bay_sensitivity <- function(log_bf_by_r, prior_odds = c(1/4, 1, 4)) {
	log_bf_by_r <- .bay_as_num(log_bf_by_r, "log_bf_by_r")
	prior_odds <- .bay_check_odds(prior_odds)
	nr <- length(log_bf_by_r)
	if (nr == 0L) {
		.exit("At least one r is required.", "Au moins un r est requis.")
	}
	r_names <- names(log_bf_by_r)
	if (is.null(r_names)) r_names <- paste0("r", seq_len(nr))
	no <- length(prior_odds)
	grid <- matrix(
		.bay_post_prob(log_bf = rep(unname(log_bf_by_r), times = no),
			prior_odds = rep(prior_odds, each = nr)),
		nrow = nr, ncol = no,
		dimnames = list(r = r_names, odds = .bay_odds_label(prior_odds))
	)
	known <- grid[!is.na(grid)]
	n_missing <- sum(is.na(grid))
	all_h1 <- n_missing == 0L && all(known > 1/2)
	all_h0 <- n_missing == 0L && all(known < 1/2)
	list(
		grid = grid,
		range = if (length(known)) range(known) else c(NA_real_, NA_real_),
		stable = isTRUE(all_h1) || isTRUE(all_h0),
		side = if (isTRUE(all_h1)) "h1" else if (isTRUE(all_h0)) "h0" else "mixed",
		n_missing = n_missing,
		prior_odds = prior_odds
	)
}


# ---------------------------------------------------------------------------
# Prior announcement
# ---------------------------------------------------------------------------

# Short statement of pi and r, to go with the probabilities: a p without its
# priors is the main fault this tool can commit (METHODOLOGIE, 8 and 11).
# Returns list(en, fr), to pass as .vbse(lab$en, lab$fr, ...). Lines are broken
# with "\n\t" to stay at 80 characters or fewer once .vbse has added its prefix.
# target says what r applies to, because "medium" differs from one BayesFactor
# function to the next (METHODOLOGIE, 2).
.bay_prior_label <- function(prior_h1 = 1/2, rscale = NULL,
                             target = c("fixed", "delta", "covariate", "random"),
                             per_term = TRUE) {
	prior_h1 <- .bay_as_num(prior_h1, "prior_h1")
	if (length(prior_h1) != 1L || is.na(prior_h1) || prior_h1 <= 0 || prior_h1 >= 1) {
		.exit("prior_h1 must be one value strictly between 0 and 1.",
			"prior_h1 doit \u00eatre une valeur strictement entre 0 et 1.")
	}
	if (!is.null(rscale)) {
		rscale <- .bay_as_num(rscale, "rscale")
		if (anyNA(rscale) || any(rscale <= 0)) {
			.exit("rscale must be strictly positive.",
				"rscale doit \u00eatre strictement positif.")
		}
	}
	targets <- c("fixed", "delta", "covariate", "random")
	if (identical(target, targets)) target <- rep_len("fixed", length(rscale))
	target <- rep_len(target, length(rscale))
	if (!all(target %in% targets)) {
		.exit("target must be fixed, delta, covariate or random.",
			"target doit valoir fixed, delta, covariate ou random.")
	}
	num <- function(x, dec) format(signif(x, 3), trim = TRUE, decimal.mark = dec)
	pi_txt <- function(dec) {
		if (isTRUE(all.equal(prior_h1, 1/2))) "1/2" else num(prior_h1, dec)
	}
	odds_txt <- .bay_odds_label(prior_h1 / (1 - prior_h1))
	target_en <- c(fixed = "fixed factor", delta = "on delta",
		covariate = "covariate", random = "random effect")
	target_fr <- c(fixed = "facteur fixe", delta = "sur delta",
		covariate = "covariable", random = "effet al\u00e9atoire")

	build <- function(dec, head, odds_word, by_term, tgt, sep) {
		first <- paste0(head, pi_txt(dec), if (isTRUE(per_term)) by_term,
			" (", odds_word, " H1:H0 = ", odds_txt, ")")
		if (is.null(rscale)) return(first)
		r_txt <- paste0("r = ", vapply(rscale, num, character(1), dec = dec),
			" (", tgt[target], ")")
		# Greedy wrap: 64 + the two tabs of .vbse cpt = "sub" (8 columns each) = 80.
		lines <- character(0)
		current <- first
		for (piece in r_txt) {
			if (nchar(current) + nchar(sep) + nchar(piece) > 64L) {
				lines <- c(lines, current)
				current <- piece
			} else {
				current <- paste0(current, sep, piece)
			}
		}
		paste(c(lines, current), collapse = "\n\t")
	}
	list(
		en = build(".", "Priors: P(H1) = ", "odds", " per term", target_en, "; "),
		fr = build(",", "A priori : P(H1) = ", "cotes", " par terme", target_fr, " ; ")
	)
}
