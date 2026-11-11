# sys_bay_posthoc.R
# Which groups differ, and by how much, for one factor (DECISIONS, D11).
#   Which: every partition of the k groups into blocks of equal means is a model; its
#     marginal likelihood is that of lmBF(y ~ h) [BayesFactor], h the merged factor (computed
#     by .bay_bf_oneway_zs, equal to lmBF for one term). "No difference" keeps the prior of
#     the global answer; the other partitions share prior_h1 by the Ewens law (Dirichlet
#     process, alpha = 1; Gopalan and Berry 1998). The probability that two groups differ is
#     the sum over the partitions that separate them: one single law for all pairs.
#   How much: draws from the posterior of the model "all groups may differ", posterior()
#     [BayesFactor]; for each pair the raw and standardised differences, medians and 95 %
#     highest density intervals, hdi() [bayestestR].
#   Many groups: up to .bay_posthoc_max_groups() groups (8: 4140 partitions, about 6 s) every
#     partition is computed. Beyond, the partitions are explored by the Gibbs sampler of
#     Gopalan and Berry (1998, section 3), each group moved in turn among the blocks of the
#     others or alone, with the same prior and the same marginal likelihoods (cached); the
#     pair probabilities are visit frequencies, with a batch-means Monte Carlo error. The
#     sensitivity to alpha reweights the same draws (Ewens normalising constant exact).
# Normal law with a common variance only; the story says so when step 4 found otherwise.
#
# References.
#   Gopalan, R., & Berry, D. A. (1998). Bayesian multiple comparisons using Dirichlet
#     process priors. JASA, 93(443), 1130-1139.
#   Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012). Default Bayes
#     factors for ANOVA designs. J. Math. Psych., 56, 356-374.
#   Makowski, D., Ben-Shachar, M. S., & Luedecke, D. (2019). bayestestR. JOSS, 4(40), 1541.

# Option KefiR.bay_exact_max, for the tests and the validation only.
.bay_posthoc_max_groups <- function(rm = FALSE) {
	if (rm) as.integer(getOption("KefiR.bay_exact_max_rm", 5L)) else as.integer(getOption("KefiR.bay_exact_max", 8L))
}

# Log Bayes factor of one grouping h against "no difference". Without subjects: the internal
# one-factor computation, equal to lmBF. With subjects (D15): lmBF(y ~ h + id) over
# lmBF(y ~ id), subject random; lb0 is the second, computed once.
.bay_partition_bf <- function(y, h, rscale, id = NULL, lb0 = 0, iterations = 10000L, cov = NULL) {
	if (is.null(id) && is.null(cov)) return(.bay_bf_oneway_zs(y, h, rscale))
	d <- .bay_nuisance_data(y, id, cov)
	d$h <- h
	BayesFactor::extractBF(BayesFactor::lmBF(stats::as.formula(paste("y ~ h +", .bay_nuisance_rhs(d))), data = d,
		whichRandom = .bay_rand_names(d), rscaleFixed = rscale, rscaleCont = .bay_rcont(),
		iterations = iterations, progress = FALSE), logbf = TRUE)$bf - lb0
}
# Model with the nuisance terms only (subjects, D15; or covariates, D16), against the intercept.
.bay_subject_lb0 <- function(y, id, iterations = 10000L, cov = NULL) {
	if (is.null(id) && is.null(cov)) return(0)
	d <- .bay_nuisance_data(y, id, cov)
	BayesFactor::extractBF(BayesFactor::lmBF(stats::as.formula(paste("y ~", .bay_nuisance_rhs(d))), data = d,
		whichRandom = .bay_rand_names(d), rscaleCont = .bay_rcont(),
		iterations = iterations, progress = FALSE), logbf = TRUE)$bf
}
.bay_nuisance_data <- function(y, id, cov) {
	d <- data.frame(y = y)
	# id: the subject, or a data frame of random factors (subject, then subject x condition
	# for replicated measures, D21), named id, id2, id3...
	if (is.data.frame(id)) {
		d$id <- id[[1L]]
		for (j in seq_along(id)[-1L]) d[[paste0("id", j)]] <- id[[j]]
	} else if (!is.null(id)) {
		d$id <- id
	}
	if (!is.null(cov)) for (i in seq_along(cov)) d[[paste0("c", i)]] <- cov[[i]]
	d
}
.bay_nuisance_rhs <- function(d) paste(setdiff(names(d), c("y", "h", "g")), collapse = " + ")
# Random terms of a nuisance data frame (NULL when none).
.bay_rand_names <- function(d) {
	r <- grep("^id[0-9]*$", names(d), value = TRUE)
	if (length(r)) r else NULL
}
# The subject factor itself.
.bay_id1 <- function(id) if (is.data.frame(id)) id[[1L]] else id

# Number of partitions of k items (Bell number), for the story.
.bay_bell <- function(k) {
	row <- 1
	for (i in seq_len(k - 1L)) {
		nxt <- numeric(i + 1L)
		nxt[1L] <- row[i]
		for (j in seq_len(i)) nxt[j + 1L] <- nxt[j] + row[j]
		row <- nxt
	}
	row[length(row)]
}

# Log of the Ewens normalising constant over the partitions other than "no difference":
# sum of alpha^b prod((size - 1)!) over all partitions is Gamma(alpha + k) / Gamma(alpha);
# the single block weighs alpha (k - 1)!.
.bay_ewens_lz <- function(k, alpha) {
	a <- lgamma(alpha + k) - lgamma(alpha)
	b <- log(alpha) + lgamma(k)
	a + log1p(-exp(b - a))
}

# Gibbs sampler over the partitions (Gopalan and Berry 1998), for many groups.
# Returns the visited partitions (matrix, one row per sweep) and the cached log Bayes factors.
.bay_partition_gibbs <- function(y, g, rscale, prior_h1, sweeps, burn, seed, alpha = 1, id = NULL, cov = NULL) {
	lb0 <- .bay_subject_lb0(y, id, cov = cov)
	gi <- as.integer(g)
	k <- nlevels(g)
	cache <- new.env(hash = TRUE)
	canon <- function(p) match(p, unique(p))
	lz <- .bay_ewens_lz(k, alpha)
	target <- function(p) {
		key <- paste(p, collapse = ".")
		lb <- cache[[key]]
		if (is.null(lb)) {
			lb <- if (max(p) == 1L) 0 else .bay_partition_bf(y, factor(p[gi]), rscale, id, lb0, cov = cov)
			assign(key, lb, envir = cache)
		}
		lp <- if (max(p) == 1L) log(1 - prior_h1) else {
			sz <- tabulate(p)
			log(prior_h1) + length(sz) * log(alpha) + sum(lgamma(sz)) - lz
		}
		lb + lp
	}
	.bay_with_seed(seed, {
		p <- seq_len(k)
		out <- matrix(0L, sweeps, k)
		for (s in seq_len(burn + sweeps)) {
			for (i in seq_len(k)) {
				others <- unique(p[-i])
				cand <- lapply(c(others, max(p) + 1L), function(b) { q <- p; q[i] <- b; canon(q) })
				cand <- cand[!duplicated(vapply(cand, paste, character(1), collapse = "."))]
				lw <- vapply(cand, target, numeric(1))
				w <- exp(lw - max(lw))
				p <- cand[[sample.int(length(cand), 1L, prob = w)]]
			}
			if (s > burn) out[s - burn, ] <- p
		}
		out
	}) -> draws
	keys <- ls(cache)
	list(draws = draws, log_bf = stats::setNames(vapply(keys, function(x) cache[[x]], numeric(1)), keys))
}

# All partitions of k items, as restricted growth strings (first item in block 1).
.bay_partitions <- function(k) {
	out <- list(1L)
	if (k == 1L) return(out)
	for (i in 2:k) {
		nxt <- list()
		for (p in out) {
			for (b in seq_len(max(p) + 1L)) nxt[[length(nxt) + 1L]] <- c(p, b)
		}
		out <- nxt
	}
	out
}

# Log prior of each partition: "no difference" at 1 - prior_h1, the others sharing prior_h1
# by the Ewens law of parameter alpha: p ~ alpha^b * prod((size - 1)!).
.bay_partition_prior <- function(parts, prior_h1, alpha = 1) {
	lw <- vapply(parts, function(p) {
		sz <- tabulate(p)
		length(sz) * log(alpha) + sum(lgamma(sz))
	}, numeric(1))
	null <- vapply(parts, function(p) max(p) == 1L, logical(1))
	lp <- numeric(length(parts))
	lp[null] <- log(1 - prior_h1)
	w <- lw[!null]
	lp[!null] <- log(prior_h1) + w - .bay_log_sum_exp(w)
	lp
}

# Label of a partition, blocks ordered by their mean: "a = c < b".
.bay_partition_label <- function(p, lev, means) {
	blocks <- split(seq_along(p), p)
	bm <- vapply(blocks, function(ix) mean(means[ix]), numeric(1))
	blocks <- blocks[order(bm)]
	paste(vapply(blocks, function(ix) paste(lev[ix][order(means[ix])], collapse = " = "), character(1)), collapse = " < ")
}

.bay_posthoc_oneway <- function(y, g, rscale, prior_h1, draws, seed, alphas = c(1, 1/2, 2), id = NULL, cov = NULL,
                                keep_draws = FALSE) {
	rm <- !is.null(id)
	if (!is.null(cov) && !ncol(cov)) cov <- NULL
	# Computed on the standardised response (Bayes factors and delta do not depend on the
	# scale), the gaps put back in the units of y: no overflow at extreme scales (1e200).
	sc <- max(abs(y - stats::median(y)))
	if (!is.finite(sc) || sc <= 0) sc <- 1
	dg0 <- .bay_dg(tapply(y, g, mean))
	y <- (y - mean(y)) / sc
	g <- droplevels(as.factor(g))
	lev <- levels(g)
	k <- length(lev)
	means <- as.numeric(tapply(y, g, mean))
	pairs <- utils::combn(k, 2L)
	res <- data.frame(group1 = lev[pairs[1L, ]], group2 = lev[pairs[2L, ]], stringsAsFactors = FALSE)
	parts_tab <- NULL
	p_sens <- NULL
	method <- if (k <= .bay_posthoc_max_groups(rm || !is.null(cov))) "exact" else "gibbs"
	mc_error <- 0
	n_sweeps <- NA_integer_
	if (method == "exact") {
		parts <- .bay_partitions(k)
		lb <- .bay_with_seed(seed, {
			lb0 <- .bay_subject_lb0(y, id, cov = cov)
			vapply(parts, function(p) {
				if (max(p) == 1L) return(0)
				.bay_partition_bf(y, factor(p[as.integer(g)]), rscale, id, lb0, cov = cov)
			}, numeric(1))
		})
		pm <- do.call(rbind, parts)
		post_of <- function(alpha) {
			lp <- .bay_partition_prior(parts, prior_h1, alpha)
			lw <- lb + lp
			exp(lw - .bay_log_sum_exp(lw))
		}
		post <- post_of(alphas[1L])
		sep_m <- pm[, pairs[1L, ], drop = FALSE] != pm[, pairs[2L, ], drop = FALSE]
		sep <- function(post) as.numeric(crossprod(sep_m, post))
		res$p_differ <- sep(post)
		p_sens <- matrix(vapply(alphas, function(a) sep(post_of(a)), numeric(ncol(pairs))), nrow = ncol(pairs))
		parts_tab <- data.frame(partition = vapply(parts, .bay_partition_label, character(1), lev = lev, means = means),
			code = vapply(parts, paste, character(1), collapse = ""),
			posterior = post, prior = exp(.bay_partition_prior(parts, prior_h1, alphas[1L])), log_bf = lb,
			stringsAsFactors = FALSE)
	} else {
		n_sweeps <- max(2000L, as.integer(draws) %/% 2L)
		gb <- .bay_partition_gibbs(y, g, rscale, prior_h1, n_sweeps, 200L, seed, alphas[1L], id = id, cov = cov)
		D <- gb$draws
		nb <- apply(D, 1L, max)
		null <- nb == 1L
		# Reweighting to another alpha: only the Ewens part of the prior changes.
		w_of <- function(a) {
			lw <- ifelse(null, 0, nb * (log(a) - log(alphas[1L])) - .bay_ewens_lz(k, a) + .bay_ewens_lz(k, alphas[1L]))
			w <- exp(lw - max(lw))
			w / sum(w)
		}
		sep_m <- D[, pairs[1L, ], drop = FALSE] != D[, pairs[2L, ], drop = FALSE]
		res$p_differ <- colMeans(sep_m)
		p_sens <- matrix(vapply(alphas, function(a) as.numeric(crossprod(sep_m, w_of(a))), numeric(ncol(pairs))), nrow = ncol(pairs))
		# Monte Carlo error: batch means over 20 batches, largest over the pairs.
		bt <- ceiling(seq_len(nrow(D)) * 20 / nrow(D))
		bm <- apply(sep_m, 2L, function(v) tapply(v, bt, mean))
		mc_error <- max(apply(bm, 2L, stats::sd) / sqrt(20))
		codes <- apply(D, 1L, paste, collapse = ".")
		fr <- table(codes)
		ucode <- names(fr)
		parts <- lapply(strsplit(ucode, ".", fixed = TRUE), as.integer)
		lp <- vapply(parts, function(p) {
			if (max(p) == 1L) return(log(1 - prior_h1))
			sz <- tabulate(p)
			log(prior_h1) + length(sz) * log(alphas[1L]) + sum(lgamma(sz)) - .bay_ewens_lz(k, alphas[1L])
		}, numeric(1))
		parts_tab <- data.frame(partition = vapply(parts, .bay_partition_label, character(1), lev = lev, means = means),
			code = gsub(".", "", ucode, fixed = TRUE),
			posterior = as.numeric(fr) / nrow(D), prior = exp(lp), log_bf = unname(gb$log_bf[ucode]),
			stringsAsFactors = FALSE)
	}
	res$p_min <- apply(p_sens, 1L, min)
	res$p_max <- apply(p_sens, 1L, max)
	parts_tab <- parts_tab[order(-parts_tab$posterior), , drop = FALSE]
	rownames(parts_tab) <- NULL

	# How much: posterior of the model where all groups may differ.
	dd <- .bay_nuisance_data(y, id, cov)
	dd$g <- g
	fo <- if (rm || !is.null(cov)) stats::as.formula(paste("y ~ g +", .bay_nuisance_rhs(dd))) else y ~ g
	po <- .bay_with_seed(seed, BayesFactor::posterior(BayesFactor::lmBF(fo, data = dd, rscaleFixed = rscale,
			rscaleCont = .bay_rcont(), whichRandom = .bay_rand_names(dd), progress = FALSE),
		iterations = max(4000L, as.integer(draws)), progress = FALSE))
	po <- as.matrix(po)
	eff <- po[, paste0("g-", lev), drop = FALSE]
	m <- sweep(eff, 1L, po[, "mu"], "+")
	# With replicates (D21), the unit is the spread of one measure of a subject in a condition:
	# residual plus subject x condition, sig2 (1 + sum of their g). This is what an
	# experimenter calls the variability of a measure once the subject is known.
	gx <- grep("^g_id[0-9]+$", colnames(po), value = TRUE)
	sig <- sqrt(po[, "sig2"] * (1 + if (length(gx)) rowSums(po[, gx, drop = FALSE]) else 0))
	hdi_of <- function(x) {
		h <- bayestestR::hdi(as.numeric(x), ci = 0.95)
		c(stats::median(x), h$CI_low, h$CI_high)
	}
	st <- t(vapply(seq_len(ncol(pairs)), function(j) {
		dif <- m[, pairs[2L, j]] - m[, pairs[1L, j]]
		c(hdi_of(dif), hdi_of(dif / sig))
	}, numeric(6)))
	res$diff <- sc * st[, 1L]; res$diff_low <- sc * st[, 2L]; res$diff_high <- sc * st[, 3L]
	res$delta <- st[, 4L]; res$delta_low <- st[, 5L]; res$delta_high <- st[, 6L]
	# The draws, for a mix with the Student law (D24).
	draws_diff <- sc * vapply(seq_len(ncol(pairs)), function(j) m[, pairs[2L, j]] - m[, pairs[1L, j]], numeric(nrow(m)))
	draws_delta <- vapply(seq_len(ncol(pairs)), function(j) (m[, pairs[2L, j]] - m[, pairs[1L, j]]) / sig, numeric(nrow(m)))
	list(pairs = res, draws_diff = if (keep_draws) draws_diff, draws_delta = if (keep_draws) draws_delta, partitions = parts_tab, k = k, alphas = alphas, draws = nrow(po),
		computed = TRUE, method = method, n_partitions = .bay_bell(k), sweeps = n_sweeps, mc_error = mc_error,
		max_exact = .bay_posthoc_max_groups(rm || !is.null(cov)), delta_unit = if (is.data.frame(id) && ncol(id) > 1L) "measure" else if (rm) "within" else if (!is.null(cov)) "adjusted" else "pooled",
		dg = dg0)
}

# Paired data, two levels (D15): size of the difference from the posterior of ttestBF on the
# differences [BayesFactor]: mean difference (second level minus first) and standardised
# difference dz (over the standard deviation of the differences), medians and 95 % HDI.
.bay_posthoc_paired <- function(y, g, id, rscale, draws, seed) {
	g <- droplevels(as.factor(g))
	lev <- levels(g)
	sc <- max(abs(y - stats::median(y)))
	if (!is.finite(sc) || sc <= 0) sc <- 1
	pr <- .bay_paired_xy((y - mean(y)) / sc, g, id)
	po <- .bay_with_seed(seed, BayesFactor::posterior(BayesFactor::ttestBF(x = pr$b, y = pr$a, paired = TRUE, rscale = rscale),
		iterations = max(4000L, as.integer(draws)), progress = FALSE))
	po <- as.matrix(po)
	hdi_of <- function(x) {
		h <- bayestestR::hdi(as.numeric(x), ci = 0.95)
		c(stats::median(x), h$CI_low, h$CI_high)
	}
	a <- hdi_of(po[, "mu"]); b <- hdi_of(po[, "delta"])
	res <- data.frame(group1 = lev[1L], group2 = lev[2L], p_differ = NA_real_, p_min = NA_real_, p_max = NA_real_,
		diff = sc * a[1L], diff_low = sc * a[2L], diff_high = sc * a[3L], delta = b[1L], delta_low = b[2L], delta_high = b[3L],
		stringsAsFactors = FALSE)
	list(pairs = res, partitions = NULL, k = 2L, draws = nrow(po), computed = TRUE, method = "paired",
		n_partitions = 2, sweeps = NA_integer_, mc_error = 0, delta_unit = "diff", n_pairs = length(pr$id),
		dg = .bay_dg(tapply(y, g, mean)))
}

# Slopes that differ between groups (D16): one adjusted gap no longer has a meaning, it
# depends on the covariate. From the posterior of lmBF(y ~ g + x + g:x) [BayesFactor], for
# each pair the gap at the quartiles of x, and the zone of the observed range where its 95 %
# interval excludes zero (Bayesian counterpart of Johnson and Neyman 1936). x: the covariate
# in interaction, centred at centre; other: the other covariates (centred), kept in the model.
.bay_posthoc_slopes <- function(y, g, x, centre, other, rscale, draws, seed) {
	g <- droplevels(as.factor(g))
	lev <- levels(g)
	k <- length(lev)
	sc <- max(abs(y - stats::median(y)))
	if (!is.finite(sc) || sc <= 0) sc <- 1
	d <- data.frame(y = (y - mean(y)) / sc, g = g, x = x)
	if (!is.null(other) && ncol(other)) for (i in seq_along(other)) d[[paste0("c", i)]] <- other[[i]]
	rhs <- paste(c("g", "x", "g:x", setdiff(names(d), c("y", "g", "x"))), collapse = " + ")
	po <- .bay_with_seed(seed, BayesFactor::posterior(BayesFactor::lmBF(stats::as.formula(paste("y ~", rhs)), data = d,
		rscaleFixed = rscale, rscaleCont = .bay_rcont(), progress = FALSE), iterations = max(4000L, as.integer(draws)), progress = FALSE))
	po <- as.matrix(po)
	cn <- colnames(po)
	a <- po[, paste0("g-", lev), drop = FALSE]
	b <- po[, grep("^g:x-", cn), drop = FALSE]
	if (ncol(b) != k) b <- matrix(0, nrow(po), k)
	hdi_of <- function(v) {
		h <- bayestestR::hdi(as.numeric(v), ci = 0.95)
		c(stats::median(v), h$CI_low, h$CI_high)
	}
	xo <- x + centre
	at <- stats::quantile(xo, c(0.25, 0.5, 0.75), names = FALSE)
	grid <- seq(min(xo), max(xo), length.out = 101L)
	pairs <- utils::combn(k, 2L)
	rows <- list(); zones <- character(ncol(pairs)); zone_list <- vector("list", ncol(pairs))
	for (j in seq_len(ncol(pairs))) {
		i1 <- pairs[1L, j]; i2 <- pairs[2L, j]
		gap <- function(x0) sc * ((a[, i2] - a[, i1]) + (b[, i2] - b[, i1]) * (x0 - centre))
		for (q in seq_along(at)) {
			h <- hdi_of(gap(at[q]))
			rows[[length(rows) + 1L]] <- data.frame(group1 = lev[i1], group2 = lev[i2], at = at[q], quartile = q,
				diff = h[1L], diff_low = h[2L], diff_high = h[3L], stringsAsFactors = FALSE)
		}
		ex <- vapply(grid, function(x0) { h <- hdi_of(gap(x0)); h[2L] > 0 || h[3L] < 0 }, logical(1))
		r <- rle(ex); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
		zone_list[[j]] <- data.frame(from = grid[starts[r$values]], to = grid[ends[r$values]])
	}
	list(method = "slopes", k = k, at = at, gaps = do.call(rbind, rows), zones = zone_list,
		pairs = data.frame(group1 = lev[pairs[1L, ]], group2 = lev[pairs[2L, ]], stringsAsFactors = FALSE),
		range = range(xo), draws = nrow(po), computed = TRUE)
}

# Several factors (PERIMETRE section 7): which levels differ, term by term. Rule:
#   - an interaction whose probability is at least 1/2, and that is not inside a more probable
#     higher-order one: its cells are compared (interaction of its factors);
#   - a main effect at least 1/2 that is in no such interaction: its levels are compared;
#   - the other factors, the covariates and the subjects are nuisance terms present in every
#     grouping, so that each gap is read at equal other factors.
# Each family is computed by .bay_posthoc_oneway (partitions of the levels, Ewens prior, same
# rules as with one factor). Terms with a covariate (slopes) are not compared here.
# Gaps of each factor of a term at each level of its other factors, from the posterior of the
# model of the cells (lmBF, posterior, hdi), the nuisance terms in it. Returns, per factor,
# the pairs (group1, group2, diff, diff_low, diff_high) and the level of the others ("at").
.bay_cell_gaps <- function(y, fl, nuis, id, rscale, draws, seed) {
	cell <- interaction(fl, sep = " x ", drop = TRUE, lex.order = TRUE)
	dd <- .bay_nuisance_data(y, id, nuis)
	dd$g <- cell
	rhs <- .bay_nuisance_rhs(dd)
	fo <- stats::as.formula(paste("y ~ g", if (nzchar(rhs)) paste("+", rhs) else ""))
	po <- .bay_with_seed(seed, BayesFactor::posterior(BayesFactor::lmBF(fo, data = dd, rscaleFixed = rscale,
			rscaleCont = .bay_rcont(), whichRandom = .bay_rand_names(dd), progress = FALSE),
		iterations = max(4000L, as.integer(draws)), progress = FALSE))
	po <- as.matrix(po)
	lev <- levels(cell)
	m <- po[, paste0("g-", lev), drop = FALSE]
	# Levels of each factor in each cell.
	first <- match(lev, as.character(cell))
	pairs <- list()
	at <- list()
	for (v in names(fl)) {
		others <- setdiff(names(fl), v)
		key <- do.call(paste, c(lapply(fl[others], function(f) as.character(f[first])), sep = " x "))
		lv <- as.character(fl[[v]][first])
		rows <- NULL
		ats <- character(0)
		for (k in unique(key)) {
			ix <- which(key == k)
			if (length(ix) < 2L) next
			cb <- utils::combn(ix, 2L)
			for (j in seq_len(ncol(cb))) {
				d <- m[, cb[2L, j]] - m[, cb[1L, j]]
				h <- bayestestR::hdi(as.numeric(d), ci = 0.95)
				rows <- rbind(rows, data.frame(group1 = lv[cb[1L, j]], group2 = lv[cb[2L, j]], diff = stats::median(d),
					diff_low = h$CI_low, diff_high = h$CI_high, stringsAsFactors = FALSE))
				ats <- c(ats, k)
			}
		}
		pairs[[v]] <- rows
		at[[v]] <- ats
	}
	list(pairs = pairs, at = at)
}

.bay_posthoc_factorial <- function(y, factors, effects, rscale, prior_h1, draws, seed, id = NULL, cov = NULL,
                                   undecided = character(0), t_weight = 0) {
	keys <- .bay_term_key(effects$term)
	fs <- strsplit(keys, ":", fixed = TRUE)
	is_fac <- vapply(fs, function(f) all(f %in% names(factors)), logical(1))
	p <- effects$p_h1
	prob <- is_fac & !is.na(p) & p >= 0.5
	ord <- lengths(fs)
	chosen <- integer(0)
	for (i in order(-ord, -p)) {
		if (!prob[i]) next
		# Inside an interaction already chosen: covered by its cells.
		if (any(vapply(chosen, function(j) all(fs[[i]] %in% fs[[j]]), logical(1)))) next
		chosen <- c(chosen, i)
	}
	fams <- lapply(chosen, function(i) {
		f <- fs[[i]]
		# Cells named in the order of the formula (the key sorts the factors).
		fo <- strsplit(gsub("`", "", effects$term[i], fixed = TRUE), ":", fixed = TRUE)[[1]]
		if (setequal(fo, f)) f <- fo
		g <- if (length(f) == 1L) factors[[f]] else interaction(factors[f], sep = " x ", drop = TRUE, lex.order = TRUE)
		other <- setdiff(names(factors), f)
		nuis <- c(if (length(other)) as.list(factors[other]), if (!is.null(cov)) as.list(cov))
		nuis <- if (length(nuis)) as.data.frame(nuis, check.names = FALSE) else NULL
		ph <- .bay_posthoc_oneway(y, g, rscale, prior_h1, draws, seed + i, id = id, cov = nuis, keep_draws = TRUE)
		if (!is.null(nuis)) ph$delta_unit <- if (is.data.frame(id) && ncol(id) > 1L) "measure" else if (!is.null(id)) "within" else "others"
		# A factor between subjects: its gaps are not measured on the residual within subject,
		# so no standardised gap is given.
		if (!is.null(id) && all(tapply(as.integer(g), .bay_id1(id), function(v) length(unique(v)) == 1L))) {
			ph$pairs$delta <- ph$pairs$delta_low <- ph$pairs$delta_high <- NA_real_
			ph$between <- TRUE
			ph$obs_means <- as.numeric(tapply(y, g, mean))
		}
		# Heavy tails (D24, audit of a206ba5): the gaps come from the same mix of laws as the
		# answer, the normal draws and the Student draws in proportion of P(Student).
		if (isTRUE(t_weight > 0.01) && is.null(id) && is.null(cov)) {
			tg <- .bay_ft_gaps(y, factors[f], factors[other], rscale, max(4000L, as.integer(draws)), seed + 31L + i)
			nn <- nrow(ph$draws_diff)
			n_t <- round(t_weight * nn)
			pick_n <- if (nn - n_t > 0L) sample.int(nn, nn - n_t) else integer(0)
			pick_t <- if (n_t > 0L) sample.int(nrow(tg$diff), n_t, replace = TRUE) else integer(0)
			dd <- rbind(ph$draws_diff[pick_n, , drop = FALSE], tg$diff[pick_t, , drop = FALSE])
			de <- rbind(ph$draws_delta[pick_n, , drop = FALSE], tg$delta[pick_t, , drop = FALSE])
			hd <- function(x) { h <- bayestestR::hdi(as.numeric(x), ci = 0.95); c(stats::median(x), h$CI_low, h$CI_high) }
			a <- t(apply(dd, 2L, hd)); b <- t(apply(de, 2L, hd))
			ph$pairs$diff <- a[, 1L]; ph$pairs$diff_low <- a[, 2L]; ph$pairs$diff_high <- a[, 3L]
			ph$pairs$delta <- b[, 1L]; ph$pairs$delta_low <- b[, 2L]; ph$pairs$delta_high <- b[, 3L]
			ph$t_weight <- t_weight
			if (t_weight >= 0.5) ph$delta_unit <- "tscale"
		}
		ph$draws_diff <- NULL
		ph$draws_delta <- NULL
		ph$term <- effects$term[i]
		ph$p_term <- p[i]
		ph$cells <- length(f) > 1L
		ph$adjusted_for <- c(other, names(cov))
		ph
	})
	skipped <- which(is_fac & !(seq_along(keys) %in% chosen) &
		!vapply(seq_along(keys), function(i) any(vapply(chosen, function(j) all(fs[[i]] %in% fs[[j]]), logical(1))), logical(1)))
	# An undecided interaction over a factor whose levels were compared as a main effect: the
	# gaps of that factor at each level of the others are given too (audit of 23c941e, point 2),
	# from the model of the cells, so that a single gap does not hide an interaction.
	by <- NULL
	for (i in which(is_fac & ord > 1L & keys %in% undecided & !(seq_along(keys) %in% chosen))) {
		fo <- strsplit(gsub("`", "", effects$term[i], fixed = TRUE), ":", fixed = TRUE)[[1]]
		if (!setequal(fo, fs[[i]])) fo <- fs[[i]]
		mains <- intersect(fo, vapply(fams, function(ph) if (ph$cells) "" else ph$term, character(1)))
		if (!length(mains)) next
		other <- setdiff(names(factors), fo)
		nuis <- c(if (length(other)) as.list(factors[other]), if (!is.null(cov)) as.list(cov))
		nuis <- if (length(nuis)) as.data.frame(nuis, check.names = FALSE) else NULL
		gp <- .bay_cell_gaps(y, factors[fo], nuis, id, rscale, draws, seed + 100L + i)
		for (v in mains) {
			at <- setdiff(fo, v)
			for (lv in unique(gp$at[[v]])) {
				sel <- gp$at[[v]] == lv
				if (!any(sel)) next
				by <- rbind(by, data.frame(interaction = effects$term[i], p_interaction = p[i], factor = v,
					at = paste0(paste(at, collapse = " x "), " = ", lv), gp$pairs[[v]][sel, , drop = FALSE],
					stringsAsFactors = FALSE))
			}
		}
	}
	rownames(by) <- NULL
	list(method = "factorial", families = fams, by = by, k = if (length(fams)) max(vapply(fams, `[[`, numeric(1), "k")) else 0L,
		skipped = data.frame(term = effects$term[skipped], p_h1 = p[skipped], stringsAsFactors = FALSE),
		pairs = if (length(fams)) do.call(rbind, lapply(fams, function(ph) data.frame(term = ph$term, ph$pairs, stringsAsFactors = FALSE))) else NULL)
}
