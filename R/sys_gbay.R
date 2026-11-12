# sys_gbay.R
# Engine, story and bilan of g.bay(), the Bayesian counterpart of the G test (DescTools::GTest).
# Decisions and references: dev/gbay/DECISIONS.md (G1 to G12). The layout functions of m.bay
# (sys_bay_story.R, sys_bay_layout.R, sys_bay_bilan.R) are called, never modified.
#
#   Goodness of fit (G2): multinomial, point H0 (the proportions given) against a Dirichlet(a)
#     alternative, closed form (Good 1967; Sarafoglou et al. 2023, multibridge).
#   Independence (G12, replaces G3 and G11): Dirichlet-multinomial with neutral margins,
#     .bay_dm_neutral_bf() (sys_bay_dm_neutral.R, shared with m.bay), a = 2, sensitivity
#     a = 1, 2, 4; the same number whether the total, the rows or the columns were fixed.
#   Estimation (G4, G6, G7): draws of the cell probabilities, uniform Dirichlet posterior
#     (contingencyTableBF(posterior = TRUE), joint multinomial: estimation only, its Bayes factor
#     is not used); odds ratio (2 x 2), ratio of each cell to its value under independence,
#     Cramer's V.

# ---------------------------------------------------------------------------------------------
# Bayes factors

# Log Dirichlet-multinomial marginal likelihood without the multinomial coefficient.
.gbay_ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))

# Goodness of fit: log BF of "proportions free, Dirichlet(a)" against "proportions = p" (G2).
.gbay_lbf_gof <- function(n, p, a = 1) {
	keep <- p > 0
	.gbay_ldm(n[keep], rep(a, sum(keep))) - sum(n[keep] * log(p[keep]))
}

# Independence (G12): Dirichlet-multinomial with neutral margins (sys_bay_dm_neutral.R).
.gbay_lbf_ct <- function(tab, a = .bay_dm_neutral_a()) .bay_dm_neutral_bf(tab, a)
# G statistic and its degrees of freedom, given for reference (the Bayes factor does not use it).
.gbay_G <- function(tab) {
	E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
	list(G = 2 * sum(ifelse(tab > 0, tab * log(tab / E), 0)), df = (nrow(tab) - 1L) * (ncol(tab) - 1L))
}

.gbay_prob <- function(lb, prior_h1) 1 / (1 + exp(-lb) * (1 - prior_h1) / prior_h1)

# ---------------------------------------------------------------------------------------------
# Posterior draws and summaries

.gbay_hdi <- function(v) {
	h <- bayestestR::hdi(as.numeric(v), ci = 0.95)
	c(stats::median(v), h$CI_low, h$CI_high)
}

# Draws of the cell probabilities (S x r x c): uniform Dirichlet posterior, Dirichlet(n + 1), from
# BayesFactor (joint multinomial; estimation only, G4, G6, G7, G11).
.gbay_draws_ct <- function(tab, draws, seed) {
	tab <- unclass(as.matrix(tab))
	r <- nrow(tab); cc <- ncol(tab)
	ps <- .bay_with_seed(seed, BayesFactor::contingencyTableBF(tab, sampleType = "jointMulti", posterior = TRUE, iterations = draws, progress = FALSE))
	m <- as.matrix(ps)[, seq_len(r * cc), drop = FALSE]
	m <- m / rowSums(m)
	array(m, c(nrow(m), r, cc))
}

# Summaries of the draws: cells, Cramer's V, odds ratio.
.gbay_summ_ct <- function(P, tab) {
	S <- dim(P)[1L]; r <- dim(P)[2L]; cc <- dim(P)[3L]
	pr <- apply(P, c(1L, 2L), sum)
	pc <- apply(P, c(1L, 3L), sum)
	E <- array(0, dim(P))
	for (j in seq_len(cc)) E[, , j] <- pr * pc[, j]
	th <- P / E
	cells <- expand.grid(row = rownames(tab), col = colnames(tab), stringsAsFactors = FALSE)
	ex <- outer(rowSums(tab), colSums(tab)) / sum(tab)
	cells$observed <- as.numeric(tab)
	cells$expected <- as.numeric(ex)
	st <- t(vapply(seq_len(r * cc), function(k) {
		i <- (k - 1L) %% r + 1L; j <- (k - 1L) %/% r + 1L
		lt <- log(th[, i, j])
		h <- .gbay_hdi(lt)
		c(exp(h), mean(lt > 0))
	}, numeric(4)))
	cells$ratio <- st[, 1L]; cells$ratio_low <- st[, 2L]; cells$ratio_high <- st[, 3L]; cells$p_more <- st[, 4L]
	cells$reading <- ifelse(cells$p_more >= 0.975, "more", ifelse(cells$p_more <= 0.025, "less", "-"))
	phi2 <- apply((P - E)^2 / E, 1L, sum)
	v <- sqrt(phi2 / (min(r, cc) - 1))
	out <- list(cells = cells, cramer_v = stats::setNames(.gbay_hdi(v), c("median", "low", "high")))
	if (r == 2L && cc == 2L) {
		lor <- log(P[, 1L, 1L]) + log(P[, 2L, 2L]) - log(P[, 1L, 2L]) - log(P[, 2L, 1L])
		h <- .gbay_hdi(lor)
		out$odds_ratio <- c(median = exp(h[1L]), low = exp(h[2L]), high = exp(h[3L]), p_gt1 = mean(lor > 0))
	}
	out
}

# Goodness of fit: Dirichlet(n + a) draws (definition through gamma variables, G2).
.gbay_summ_gof <- function(n, p, a, draws, seed) {
	k <- length(n)
	P <- .bay_with_seed(seed, {
		g <- matrix(stats::rgamma(draws * k, shape = rep(n + a, each = draws)), draws, k)
		g / rowSums(g)
	})
	keep <- p > 0
	cats <- data.frame(category = names(n), observed = as.numeric(n), expected = sum(n) * p, p = p, stringsAsFactors = FALSE)
	st <- t(vapply(seq_len(k), function(i) {
		pi <- P[, i]
		h <- .gbay_hdi(pi)
		lr <- log(pi / p[i])
		hr <- .gbay_hdi(lr)
		c(h, exp(hr), mean(lr > 0))
	}, numeric(7)))
	cats$prop <- st[, 1L]; cats$prop_low <- st[, 2L]; cats$prop_high <- st[, 3L]
	cats$ratio <- st[, 4L]; cats$ratio_low <- st[, 5L]; cats$ratio_high <- st[, 6L]; cats$p_more <- st[, 7L]
	cats$reading <- ifelse(cats$p_more >= 0.975, "more", ifelse(cats$p_more <= 0.025, "less", "-"))
	w <- sqrt(rowSums(sweep(P[, keep, drop = FALSE], 2L, p[keep])^2 / rep(p[keep], each = draws)))
	list(cells = cats, w = stats::setNames(.gbay_hdi(w), c("median", "low", "high")))
}

# ---------------------------------------------------------------------------------------------
# Input (G1; PERIMETRE)

.gbay_err <- function(en, fr) stop(.msg(en, fr), call. = FALSE)

.gbay_whole <- function(v) is.numeric(v) && all(is.finite(v)) && all(v >= 0) && all(abs(v - round(v)) < 1e-8)

.gbay_input <- function(x, y, p, data, rescale.p, x_name, y_name, env) {
	res <- list(changes = list(), n_removed = 0L, form = NULL, fixed = NULL, resp = FALSE)
	change <- function(en, fr) res$changes[[length(res$changes) + 1L]] <<- list(en = en, fr = fr)
	col <- function(v, nm) {
		if (is.character(v) && length(v) == 1L && !is.null(data) && v %in% names(data)) return(list(v = data[[v]], name = v))
		if (is.numeric(v) && length(v) == 1L && !is.null(data) && v >= 1 && v <= ncol(data) && v == round(v)) return(list(v = data[[v]], name = names(data)[v]))
		list(v = v, name = nm)
	}
	if (inherits(x, "formula")) {
		f <- x
		if (length(f) == 3L) { lhs <- f[[2L]]; rhs <- f[[3L]] } else { lhs <- NULL; rhs <- f[[2L]] }
		tl <- attr(stats::terms(stats::as.formula(call("~", rhs))), "term.labels")
		if (any(grepl(":", tl, fixed = TRUE)) || !all(tl %in% all.vars(rhs)) || !length(tl) || length(tl) > 2L) {
			.gbay_err("The formula takes one or two variables joined by +: ~ A, ~ A + B, y ~ g or n ~ A + B.",
				"La formule prend une ou deux variables jointes par + : ~ A, ~ A + B, y ~ g ou n ~ A + B.")
		}
		src <- if (is.null(data)) environment(f) else data
		get_v <- function(e) {
			v <- tryCatch(eval(e, if (is.data.frame(src)) src else NULL, if (is.environment(src)) src else env), error = function(err) NULL)
			if (is.null(v)) .gbay_err(paste0("Variable not found: ", deparse1(e), "."), paste0("Variable introuvable : ", deparse1(e), "."))
			v
		}
		rv <- lapply(tl, function(t) get_v(str2lang(t)))
		lv <- if (!is.null(lhs)) get_v(lhs) else NULL
		vars <- c(if (!is.null(lv)) list(lv), rv)
		if (length(unique(lengths(vars))) != 1L) .gbay_err("The variables of the formula have different lengths.", "Les variables de la formule n'ont pas la m\u00eame longueur.")
		ok <- stats::complete.cases(as.data.frame(lapply(vars, function(v) if (is.factor(v)) as.character(v) else v), stringsAsFactors = FALSE))
		res$n_removed <- sum(!ok)
		rv <- lapply(rv, function(v) v[ok]); if (!is.null(lv)) lv <- lv[ok]
		rf <- lapply(rv, function(v) if (is.factor(v)) droplevels(v) else factor(v))
		if (is.null(lv)) {
			res$form <- if (length(tl) == 1L) "formula_one" else "formula_two"
			res$names <- tl
			tab <- if (length(tl) == 1L) table(rf[[1L]]) else table(rf[[1L]], rf[[2L]])
		} else {
			key <- interaction(rf, drop = TRUE)
			weights <- .gbay_whole(lv) && !anyDuplicated(key) && length(lv) >= 2L
			if (weights) {
				res$form <- if (length(tl) == 1L) "formula_one" else "formula_two"
				res$weights <- deparse1(lhs)
				res$names <- tl
				tab <- if (length(tl) == 1L) tapply(lv, rf[[1L]], sum) else tapply(lv, list(rf[[1L]], rf[[2L]]), sum)
				tab[is.na(tab)] <- 0
				tab <- as.table(tab)
			} else {
				if (is.numeric(lv) && length(unique(lv)) > 2L) {
					.gbay_err(paste0("The left side (", deparse1(lhs), ") is neither categories nor counts (one row per combination): g.bay needs counts. For a numeric response, use m.bay()."),
						paste0("Le membre de gauche (", deparse1(lhs), ") n'est ni des modalit\u00e9s ni des effectifs (une ligne par combinaison) : g.bay demande des effectifs. Pour une r\u00e9ponse num\u00e9rique, utiliser m.bay()."))
				}
				if (length(tl) != 1L) {
					.gbay_err("A response with two factors makes three variables: g.bay handles two. Use m.bay(y ~ A * B).",
						"Une r\u00e9ponse avec deux facteurs fait trois variables : g.bay en traite deux. Utiliser m.bay(y ~ A * B).")
				}
				res$form <- "formula_resp"
				res$resp <- TRUE
				res$names <- c(tl, deparse1(lhs))
				yf <- if (is.factor(lv)) droplevels(lv) else factor(lv)
				tab <- table(rf[[1L]], yf)
				res$fixed <- "rows"
			}
		}
	} else {
		xc <- col(x, x_name); x <- xc$v
		if (!is.null(y)) {
			yc <- col(y, y_name); y <- yc$v
			if (is.matrix(x) || is.data.frame(x) || length(dim(x)) >= 2L) {
				.gbay_err("y is given with a table x: give either a table, or two vectors x and y.",
					"y est donn\u00e9 avec une table x : donner soit une table, soit deux vecteurs x et y.")
			}
			if (length(x) != length(y)) .gbay_err("x and y must have the same length.", "x et y doivent avoir la m\u00eame longueur.")
			ok <- !is.na(x) & !is.na(y)
			res$n_removed <- sum(!ok)
			xf <- if (is.factor(x)) droplevels(x[ok]) else factor(x[ok])
			yf <- if (is.factor(y)) droplevels(y[ok]) else factor(y[ok])
			if (nlevels(xf) < 2L || nlevels(yf) < 2L) .gbay_err("x and y must each have at least 2 levels.", "x et y doivent avoir chacun au moins 2 modalit\u00e9s.")
			tab <- table(xf, yf)
			res$names <- c(xc$name, yc$name)
			res$form <- "xy"
		} else {
			if (is.data.frame(x)) x <- as.matrix(x)
			if (length(dim(x)) > 2L) {
				.gbay_err("Tables with three dimensions or more are not handled (nor by GTest).",
					"Les tables \u00e0 trois dimensions ou plus ne sont pas trait\u00e9es (ni par GTest).")
			}
			if (is.factor(x) || is.character(x) || is.logical(x)) {
				ok <- !is.na(x); res$n_removed <- sum(!ok)
				tab <- table(if (is.factor(x)) droplevels(x[ok]) else factor(x[ok]))
				res$names <- xc$name
				res$form <- "raw"
			} else if (length(dim(x)) == 2L) {
				if (!is.numeric(x)) .gbay_err("x must hold counts.", "x doit contenir des effectifs.")
				if (min(dim(x)) == 1L) {
					change(paste0("a table with a single ", if (nrow(x) == 1L) "row" else "column", ", treated as one variable (goodness of fit), as GTest does."),
						paste0("table \u00e0 une seule ", if (nrow(x) == 1L) "ligne" else "colonne", ", trait\u00e9e comme une seule variable (ajustement), comme le fait GTest."))
					nm <- if (nrow(x) == 1L) colnames(x) else rownames(x)
					x <- as.vector(x); names(x) <- nm
					tab <- x
					res$names <- xc$name
					res$form <- "vector"
				} else {
					tab <- as.table(x)
					dn <- names(dimnames(tab))
					res$names <- if (length(dn) == 2L && all(nzchar(dn))) dn else NULL
					res$form <- "table"
				}
			} else {
				if (!is.numeric(x)) .gbay_err("x must hold counts.", "x doit contenir des effectifs.")
				tab <- x
				if (!is.null(dim(x))) { nm <- names(x); tab <- as.vector(x); names(tab) <- nm }
				res$names <- xc$name
				res$form <- "vector"
			}
		}
	}
	vals <- as.numeric(tab)
	if (anyNA(vals) || !.gbay_whole(vals)) {
		.gbay_err("Counts must be whole numbers, positive or zero, without missing values.",
			"Les effectifs doivent \u00eatre des nombres entiers, positifs ou nuls, sans valeur manquante.")
	}
	if (sum(vals) == 0) .gbay_err("At least one count must be positive.", "Au moins un effectif doit \u00eatre positif.")
	if (length(dim(tab)) == 2L) {
		tab <- unclass(as.matrix(tab)); storage.mode(tab) <- "double"
		if (is.null(rownames(tab))) rownames(tab) <- paste0("r", seq_len(nrow(tab)))
		if (is.null(colnames(tab))) colnames(tab) <- paste0("c", seq_len(ncol(tab)))
		zr <- rowSums(tab) == 0; zc <- colSums(tab) == 0
		if (any(zr) || any(zc)) {
			wh <- c(rownames(tab)[zr], colnames(tab)[zc])
			change(paste0("empty row or column removed (", paste(wh, collapse = ", "), "): a level never observed says nothing about the link and would change the answer (DECISIONS, G3)."),
				paste0("ligne ou colonne vide retir\u00e9e (", paste(wh, collapse = ", "), ") : une modalit\u00e9 jamais observ\u00e9e ne dit rien du lien et changerait la r\u00e9ponse (DECISIONS, G3)."))
			tab <- tab[!zr, !zc, drop = FALSE]
		}
		if (nrow(tab) < 2L || ncol(tab) < 2L) {
			.gbay_err("Fewer than two rows or columns hold counts: no link can be studied.",
				"Moins de deux lignes ou colonnes portent des effectifs : aucun lien n'est \u00e9tudiable.")
		}
		names(dimnames(tab)) <- NULL
		res$kind <- "independence"
		if (!is.null(p)) {
			change("p is ignored for a two-way table (independence), as in GTest.", "p est ignor\u00e9 pour une table \u00e0 deux entr\u00e9es (ind\u00e9pendance), comme dans GTest.")
		}
	} else {
		n <- as.numeric(tab); nm <- names(tab)
		if (is.null(nm) || !all(nzchar(nm))) nm <- paste0("c", seq_along(n))
		names(n) <- nm
		if (length(n) < 2L) .gbay_err("x must have at least 2 categories.", "x doit avoir au moins 2 modalit\u00e9s.")
		if (is.null(p)) {
			p <- rep(1 / length(n), length(n))
			res$p_given <- FALSE
		} else {
			res$p_given <- TRUE
			if (!is.numeric(p) || length(p) != length(n) || anyNA(p)) {
				.gbay_err("x and p must have the same number of elements.", "x et p doivent avoir le m\u00eame nombre d'\u00e9l\u00e9ments.")
			}
			if (any(p < 0)) .gbay_err("Probabilities must be positive or zero.", "Les probabilit\u00e9s doivent \u00eatre positives ou nulles.")
			if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
				if (isTRUE(rescale.p)) p <- p / sum(p) else {
					.gbay_err("Probabilities must sum to 1 (or rescale.p = TRUE).", "Les probabilit\u00e9s doivent faire 1 au total (ou rescale.p = TRUE).")
				}
			}
			if (any(p == 0 & n > 0)) {
				.gbay_err("A category of theoretical proportion 0 was observed: the proportions given are impossible.",
					"Une modalit\u00e9 de proportion th\u00e9orique nulle a \u00e9t\u00e9 observ\u00e9e : les proportions donn\u00e9es sont impossibles.")
			}
			if (any(p == 0)) {
				change("category of theoretical proportion 0, never observed, removed.", "modalit\u00e9 de proportion th\u00e9orique nulle, jamais observ\u00e9e, retir\u00e9e.")
				keep <- p > 0; n <- n[keep]; p <- p[keep]
				if (length(n) < 2L) .gbay_err("Fewer than 2 categories left.", "Il reste moins de 2 modalit\u00e9s.")
			}
		}
		names(p) <- names(n)
		tab <- n
		res$p <- p
		res$kind <- "gof"
	}
	res$tab <- tab
	res
}

# ---------------------------------------------------------------------------------------------
# Words

.gbay_fmt <- function(x, dc, dg = 2L) gsub(".", dc, formatC(x, format = "f", digits = dg), fixed = TRUE)
# A ratio with 3 significant digits at most, 2 decimals under 10.
.gbay_rt <- function(x, dc) {
	if (!is.finite(x)) return("NA")
	dg <- if (x >= 100) 0L else if (x >= 10) 1L else 2L
	gsub(".", dc, formatC(x, format = "f", digits = dg, big.mark = ""), fixed = TRUE)
}
.gbay_ci <- function(v, dc, f = .gbay_rt) {
	sep <- if (dc == ",") " ; " else "; "
	paste0(f(v[1L], dc), " [", f(v[2L], dc), sep, f(v[3L], dc), "]")
}

# Labels of the hypotheses (G8).
.gbay_labels <- function(obj) {
	if (identical(obj$test, "gof")) {
		a <- obj$names[1L]
		return(list(h1_en = paste0("The proportions of ", a, " depart from those expected"),
			h1_fr = paste0("Les proportions de ", a, " s'\u00e9cartent de celles attendues"),
			h0_en = paste0("The proportions of ", a, " follow those expected"),
			h0_fr = paste0("Les proportions de ", a, " suivent celles attendues")))
	}
	if (isTRUE(obj$resp)) {
		g <- obj$names[1L]; y <- obj$names[2L]
		return(list(h1_en = paste0("The distribution of ", y, " depends on ", g), h1_fr = paste0("La r\u00e9partition de ", y, " d\u00e9pend de ", g),
			h0_en = paste0("The distribution of ", y, " does not depend on ", g), h0_fr = paste0("La r\u00e9partition de ", y, " ne d\u00e9pend pas de ", g)))
	}
	a <- obj$names[1L]; b <- obj$names[2L]
	list(h1_en = paste0("Link between ", a, " and ", b), h1_fr = paste0("Lien entre ", a, " et ", b),
		h0_en = paste0(a, " and ", b, " independent"), h0_fr = paste0(a, " et ", b, " ind\u00e9pendants"))
}

.gbay_verdict <- function(obj, undecided) {
	h <- obj$hypotheses
	p <- h$posterior[1L]
	gof <- identical(obj$test, "gof")
	if (undecided) {
		return(if (gof) list(en = "Verdict: undecided. The data allow no conclusion, neither a departure nor a fit.",
			fr = "Verdict : ind\u00e9cis. Les donn\u00e9es ne permettent de conclure ni \u00e0 un \u00e9cart, ni \u00e0 l'ajustement.") else
			list(en = "Verdict: undecided. The data allow no conclusion, neither a link nor independence.",
				fr = "Verdict : ind\u00e9cis. Les donn\u00e9es ne permettent de conclure ni \u00e0 un lien, ni \u00e0 l'ind\u00e9pendance."))
	}
	if (h$null[1L]) {
		if (p >= 0.95) {
			return(if (gof) list(en = "Verdict: very probably no departure of the size expected. This is not a proof that the proportions are exactly those given.",
				fr = "Verdict : tr\u00e8s probablement pas d'\u00e9cart de l'ordre attendu. Ce n'est pas une preuve que les proportions sont exactement celles donn\u00e9es.") else
				list(en = "Verdict: very probably no link of the size expected. This is not a proof of strict independence.",
					fr = "Verdict : tr\u00e8s probablement pas de lien de l'ordre attendu. Ce n'est pas une preuve d'ind\u00e9pendance stricte."))
		}
		return(if (gof) list(en = "Verdict: rather no departure, without concluding.", fr = "Verdict : plut\u00f4t pas d'\u00e9cart, sans conclure.") else
			list(en = "Verdict: rather no link, without concluding.", fr = "Verdict : plut\u00f4t pas de lien, sans conclure."))
	}
	what_en <- if (gof) "a departure" else "a link"; what_fr <- if (gof) "un \u00e9cart" else "un lien"
	tail_en <- if (gof) " This says neither that it matters nor why." else " This says neither that it is strong nor why (not a cause)."
	tail_fr <- if (gof) " Cela ne dit ni qu'il est important, ni pourquoi." else " Cela ne dit ni qu'il est fort, ni pourquoi (pas une cause)."
	if (p >= 0.95) return(list(en = paste0("Verdict: ", what_en, " very probably exists.", tail_en), fr = paste0("Verdict : ", what_fr, " existe tr\u00e8s probablement.", tail_fr)))
	if (p >= 0.75) return(list(en = paste0("Verdict: ", what_en, " probably exists.", tail_en), fr = paste0("Verdict : ", what_fr, " existe probablement.", tail_fr)))
	list(en = paste0("Verdict: rather ", what_en, ", without concluding."), fr = paste0("Verdict : plut\u00f4t ", what_fr, ", sans conclure."))
}

# Lines on the cells (or categories) that carry the gap, the strongest first, at most `top`.
.gbay_cell_lines <- function(obj, emit, top = 4L, short = FALSE) {
	ce <- obj$cells
	if (is.null(ce)) return(invisible(FALSE))
	gof <- identical(obj$test, "gof")
	lab <- if (gof) ce$category else paste0("(", ce$row, ", ", ce$col, ")")
	o <- order(-abs(log(ce$ratio)))
	shown <- FALSE
	for (side in c("more", "less")) {
		sel <- o[ce$reading[o] == side]
		if (!length(sel)) next
		shown <- TRUE
		more <- length(sel) > top
		sel <- utils::head(sel, top)
		if (short) {
			txt_en <- paste0(lab[sel], " x ", vapply(ce$ratio[sel], .gbay_rt, character(1), dc = "."), collapse = ", ")
			txt_fr <- paste0(lab[sel], " x ", vapply(ce$ratio[sel], .gbay_rt, character(1), dc = ","), collapse = ", ")
		} else {
			txt_en <- paste0(lab[sel], " ", vapply(sel, function(i) .gbay_ci(c(ce$ratio[i], ce$ratio_low[i], ce$ratio_high[i]), "."), character(1)), " times", collapse = "; ")
			txt_fr <- paste0(lab[sel], " ", vapply(sel, function(i) .gbay_ci(c(ce$ratio[i], ce$ratio_low[i], ce$ratio_high[i]), ","), character(1)), " fois", collapse = " ; ")
		}
		if (more) { txt_en <- paste0(txt_en, "; others in $cells"); txt_fr <- paste0(txt_fr, " ; autres dans $cells") }
		if (side == "more") {
			emit("fact", paste0(if (gof) "More frequent than expected: " else "More frequent than without a link: ", txt_en, "."),
				paste0(if (gof) "Plus fr\u00e9quentes que pr\u00e9vu : " else "Plus fr\u00e9quentes que sans lien : ", txt_fr, "."))
		} else {
			emit("fact", paste0(if (gof) "Less frequent than expected: " else "Less frequent than without a link: ", txt_en, "."),
				paste0(if (gof) "Moins fr\u00e9quentes que pr\u00e9vu : " else "Moins fr\u00e9quentes que sans lien : ", txt_fr, "."))
		}
	}
	if (!shown) {
		emit("fact", if (gof) "No category departs clearly from its expected proportion (every 95 % interval holds 1)." else "No cell departs clearly from independence (every 95 % interval holds 1).",
			if (gof) "Aucune modalit\u00e9 ne s'\u00e9carte nettement de sa proportion attendue (chaque intervalle \u00e0 95 % contient 1)." else "Aucune case ne s'\u00e9carte nettement de l'ind\u00e9pendance (chaque intervalle \u00e0 95 % contient 1).")
	}
	invisible(shown)
}

# Size of the effect, one line.
.gbay_size_line <- function(obj, emit, short = FALSE) {
	if (identical(obj$test, "gof")) {
		if (is.null(obj$w)) return(invisible(NULL))
		emit("fact", paste0("Size of the departure (Cohen's w): ", .gbay_ci(obj$w, ".", function(x, dc) .gbay_fmt(x, dc)), "."),
			paste0("Taille de l'\u00e9cart (w de Cohen) : ", .gbay_ci(obj$w, ",", function(x, dc) .gbay_fmt(x, dc)), "."))
		return(invisible(NULL))
	}
	if (!is.null(obj$odds_ratio)) {
		t <- obj$table; o <- obj$odds_ratio
		if (short) {
			emit("fact", paste0("Odds ratio (", rownames(t)[1L], " against ", rownames(t)[2L], "): ", .gbay_ci(o, "."), "."),
				paste0("Rapport des cotes (", rownames(t)[1L], " contre ", rownames(t)[2L], ") : ", .gbay_ci(o, ","), "."))
		} else {
			emit("fact", paste0("Odds ratio: the odds of ", colnames(t)[1L], " against ", colnames(t)[2L], " are multiplied by ", .gbay_ci(o, "."), " in ", rownames(t)[1L], " compared with ", rownames(t)[2L],
					"; probability that the ratio exceeds 1: ", .bay_pct(o[["p_gt1"]], "."), "."),
				paste0("Rapport des cotes : la cote de ", colnames(t)[1L], " contre ", colnames(t)[2L], " est multipli\u00e9e par ", .gbay_ci(o, ","), " en ", rownames(t)[1L], " par rapport \u00e0 ", rownames(t)[2L],
					" ; probabilit\u00e9 que ce rapport d\u00e9passe 1 : ", .bay_pct(o[["p_gt1"]], ","), "."))
		}
	}
	if (!is.null(obj$cramer_v)) {
		emit("fact", paste0("Size of the link (Cramer's V): ", .gbay_ci(obj$cramer_v, ".", function(x, dc) .gbay_fmt(x, dc)), "."),
			paste0("Taille du lien (V de Cram\u00e9r) : ", .gbay_ci(obj$cramer_v, ",", function(x, dc) .gbay_fmt(x, dc)), "."))
	}
	invisible(NULL)
}

.gbay_cite <- function(obj) {
	if (identical(obj$test, "gof")) return(list(en = "Dirichlet against the proportions given, Good 1967", fr = "Dirichlet contre les proportions donn\u00e9es, Good 1967"))
	list(en = "Dirichlet-multinomial with neutral margins (BDe), a = 2, Gunel and Dickey 1974", fr = "Dirichlet-multinomiale \u00e0 marges neutres (BDe), a = 2, Gunel et Dickey 1974")
}

# ---------------------------------------------------------------------------------------------
# Bilan

# Coherence guard (G12): independence (or the fit) strongly supported while cells depart from it.
.gbay_guard <- function(obj, li = .bay_lead_info(obj)) {
	n_flag <- if (is.null(obj$cells)) 0L else sum(obj$cells$reading != "-")
	list(on = isTRUE(li$h$null[1L]) && li$sup$level >= 4L && n_flag > 0L, n_flag = n_flag)
}

.gbay_print_bilan <- function(obj, detail = FALSE, narrated = FALSE) {
	emit <- .bay_block_emitter()
	li <- .bay_lead_info(obj)
	h <- li$h
	gof <- identical(obj$test, "gof")
	# Coherence guard (G12): no "strong" support for independence (or for the proportions given)
	# while the estimation finds cells whose 95 % interval excludes 1, i.e. a gap that it does not
	# attribute to chance. The support is brought back to moderate, and the bilan says why.
	gd <- .gbay_guard(obj, li)
	n_flag <- gd$n_flag
	guard <- gd$on
	if (guard) li$sup <- .bay_support(0.8)
	undecided <- li$open || li$sup$level <= 1L || h$posterior[1L] < 0.6
	lead_en <- if (undecided) paste0("undecided (", h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), ")") else
		paste0(h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), " (", li$sup$en, ")")
	lead_fr <- if (undecided) paste0("ind\u00e9cis (", h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), ")") else
		paste0(h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), " (", li$sup$fr, ")")
	if (isTRUE(narrated)) {
		cat("\n")
		k <- if (!is.null(obj$story$k)) obj$story$k + 1L else 7L
		emit("title", paste0(k, ") Bilan - ", lead_en, "."), paste0(k, ") Bilan - ", lead_fr, "."))
	} else {
		emit("title", paste0("Bilan of g.bay (", obj$label, ", N = ", obj$n, ") - ", lead_en, "."),
			paste0("Bilan de g.bay (", obj$label, ", N = ", obj$n, ") - ", lead_fr, "."))
	}
	v <- .gbay_verdict(obj, undecided)
	if (guard && !undecided) {
		v <- if (gof) list(en = "Verdict: rather no departure, without concluding.", fr = "Verdict : plut\u00f4t pas d'\u00e9cart, sans conclure.") else
			list(en = "Verdict: rather no link, without concluding.", fr = "Verdict : plut\u00f4t pas de lien, sans conclure.")
	}
	emit("fact", v$en, v$fr)
	if (guard) {
		emit("decide", paste0("Caution: ", if (gof) "the fit" else "independence", " leads, but ", n_flag, if (n_flag == 1L) (if (gof) " category departs" else " cell departs") else (if (gof) " categories depart" else " cells depart"),
				" from it with a 95 % interval that excludes 1; its support is brought back to moderate. The Bayes factor prefers the simpler model; the data do not rule out ",
				if (gof) "a departure." else "a link."),
			paste0("Prudence : ", if (gof) "l'ajustement" else "l'ind\u00e9pendance", " est en t\u00eate, mais ", n_flag, if (n_flag == 1L) (if (gof) " modalit\u00e9 s'en \u00e9carte" else " case s'en \u00e9carte") else (if (gof) " modalit\u00e9s s'en \u00e9cartent" else " cases s'en \u00e9cartent"),
				" avec un intervalle \u00e0 95 % qui exclut 1 ; son soutien est ramen\u00e9 \u00e0 mod\u00e9r\u00e9. Le facteur de Bayes pr\u00e9f\u00e8re le mod\u00e8le le plus simple ; les donn\u00e9es n'excluent pas ",
				if (gof) "un \u00e9cart." else "un lien."))
	}
	w_en <- max(nchar(h$label_en)); w_fr <- max(nchar(h$label_fr))
	fit <- max(w_en, w_fr) + 32L <= .bay_width() - 12L
	pad <- function(x, w) if (fit) paste0(x, strrep("\u00a0", w - nchar(x))) else x
	emit("fact", "Probability of each hypothesis:", "Probabilit\u00e9 de chaque hypoth\u00e8se :")
	for (i in seq_len(nrow(h))) {
		emit("sub", paste0(i, ". ", pad(h$label_en[i], w_en), " : ", .bay_pct(h$posterior[i], "."), if (i == 1L) paste0(" (", li$sup$en, ")") else ""),
			paste0(i, ". ", pad(h$label_fr[i], w_fr), " : ", .bay_pct(h$posterior[i], ","), if (i == 1L) paste0(" (", li$sup$fr, ")") else ""))
	}
	bf_en <- .bay_bf_txt(abs(li$lb), "."); bf_fr <- .bay_bf_txt(abs(li$lb), ",")
	emit("note", paste0("the counts are ", bf_en, " times more probable under hypothesis 1 than under 2 (Bayes factor; 1 to 3 weak, 3 to 20 moderate, 20 to 150 strong, beyond very strong)."),
		paste0("les effectifs sont ", bf_fr, " fois plus probables sous l'hypoth\u00e8se 1 que sous la 2 (facteur de Bayes ; 1 \u00e0 3 faible, 3 \u00e0 20 mod\u00e9r\u00e9, 20 \u00e0 150 fort, au-del\u00e0 tr\u00e8s fort)."))
	sens <- obj$sensitivity
	if (!is.null(sens$range_lead)) {
		rg <- sens$range_lead
		why_en <- why_fr <- ""
		if (li$open) {
			emit("decide", paste0("Depending on the starting point", why_en, ": ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order can reverse, not settled."),
				paste0("Selon le point de d\u00e9part", why_fr, " : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre peut s'inverser, non tranch\u00e9."))
		} else if (identical(.bay_pct(rg[1L], "."), .bay_pct(rg[2L], ".")) || rg[1L] > 0.999) {
			emit("fact", paste0("The answer ", .bay_range_txt(rg)$en, "."), paste0("La r\u00e9ponse ", .bay_range_txt(rg)$fr, "."))
		} else {
			emit("fact", paste0("Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order does not change."),
				paste0("Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre ne change pas."))
		}
	} else if (!isTRUE(narrated)) {
		emit("note", "sensitivity to the priors not computed (sensitivity = FALSE).", "sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}
	emit("block", "", "")
	.gbay_size_line(obj, emit, short = isTRUE(narrated))
	if (!is.null(obj$cells)) {
		if (guard) {
			.gbay_cell_lines(obj, emit, top = 3L, short = isTRUE(narrated))
		} else if (h$null[1L] || undecided) {
			emit("fact", if (gof) "No departure established: the categories are not read one by one (they are in $cells)." else "No link established: the cells are not read one by one (they are in $cells).",
				if (gof) "Aucun \u00e9cart \u00e9tabli : les modalit\u00e9s ne sont pas lues une \u00e0 une (elles sont dans $cells)." else "Aucun lien \u00e9tabli : les cases ne sont pas lues une \u00e0 une (elles sont dans $cells).")
		} else {
			.gbay_cell_lines(obj, emit, top = 3L, short = isTRUE(narrated))
		}
		if (isTRUE(narrated)) {
			emit("note", if (gof) "intervals at step 6, every category in $cells." else "intervals at step 6, every cell in $cells.", if (gof) "intervalles \u00e0 l'\u00e9tape 6, toutes les modalit\u00e9s dans $cells." else "intervalles \u00e0 l'\u00e9tape 6, toutes les cases dans $cells.")
		} else {
			emit("note", if (gof) "ratio = proportion estimated / proportion given; [ ; ] = 95 % credible interval." else "ratio = cell / its value without a link; [ ; ] = 95 % credible interval.",
				if (gof) "rapport = proportion estim\u00e9e / proportion donn\u00e9e ; [ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %." else "rapport = case / sa valeur sans lien ; [ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %.")
		}
	}
	emit("block", "", "")
	.gbay_observed(obj, emit)
	ct <- .gbay_cite(obj)
	h1 <- which(!h$null)[1L]
	size_en <- size_fr <- ""
	if (!gof && !is.null(obj$cramer_v)) {
		size_en <- paste0("; Cramer's V ", .gbay_ci(obj$cramer_v, ".", function(x, dc) .gbay_fmt(x, dc)))
		size_fr <- paste0(" ; V de Cram\u00e9r ", .gbay_ci(obj$cramer_v, ",", function(x, dc) .gbay_fmt(x, dc)))
	}
	if (gof && !is.null(obj$w)) {
		size_en <- paste0("; Cohen's w ", .gbay_ci(obj$w, ".", function(x, dc) .gbay_fmt(x, dc)))
		size_fr <- paste0(" ; w de Cohen ", .gbay_ci(obj$w, ",", function(x, dc) .gbay_fmt(x, dc)))
	}
	if (undecided) {
		what_en <- if (gof) "a departure" else "a link"; what_fr <- if (gof) "d'un \u00e9cart" else "d'un lien"
		against_en <- if (gof) " in favour of the proportions given" else " in favour of independence"
		against_fr <- if (gof) " en faveur des proportions donn\u00e9es" else " en faveur de l'ind\u00e9pendance"
		emit("note", paste0("for a report: \"The data do not settle it (probability of ", what_en, " ", .bay_pct(h$posterior[h1], "."), ", Bayes factor ", bf_en,
				if (h$null[1L]) against_en else "", "); g.bay {KefiR}, ", ct$en, ".\""),
			paste0("pour un compte rendu : \u00ab Les donn\u00e9es ne permettent pas de conclure (probabilit\u00e9 ", what_fr, " ", .bay_pct(h$posterior[h1], ","), ", facteur de Bayes ", bf_fr,
				if (h$null[1L]) against_fr else "", ") ; g.bay {KefiR}, ", ct$fr, ". \u00bb"))
	} else {
		emit("note", paste0("for a report: \"", h$label_en[1L], " (posterior probability ", .bay_pct(h$posterior[1L], "."), "), Bayes factor ", bf_en,
				if (!h$null[1L]) size_en else "", if (guard) "; support brought back to moderate, cells depart" else "", " (g.bay {KefiR}, ", ct$en, ").\""),
			paste0("pour un compte rendu : \u00ab ", h$label_fr[1L], " (probabilit\u00e9 a posteriori ", .bay_pct(h$posterior[1L], ","), "), facteur de Bayes ", bf_fr,
				if (!h$null[1L]) size_fr else "", if (guard) " ; soutien ramen\u00e9 \u00e0 mod\u00e9r\u00e9, des cases s'en \u00e9cartent" else "", " (g.bay {KefiR}, ", ct$fr, "). \u00bb"))
	}
	if (length(obj$notes$en) && !isTRUE(narrated)) {
		for (i in seq_along(obj$notes$en)) emit("fact", obj$notes$en[i], obj$notes$fr[i])
	}
	emit("note", "probabilities of hypotheses, not p-values; the method with print(result, detail = TRUE).",
		"probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values ; la m\u00e9thode avec print(resultat, detail = TRUE).")
	emit("flush", "", "")
	if (isTRUE(detail)) .gbay_print_detail(obj)
	invisible(obj)
}

# What was observed: the largest gap between counts and expected counts.
.gbay_observed <- function(obj, emit) {
	if (identical(obj$test, "gof")) {
		n <- obj$table; e <- obj$expected
		i <- which.max(abs(n - e))
		emit("fact", paste0("Observed: ", paste0(n, " ", names(n), collapse = ", "), " (of ", sum(n), "); largest gap: ", names(n)[i], ", ", n[i], " against ", .gbay_fmt(e[i], ".", 1L), " expected."),
			paste0("Observ\u00e9 : ", paste0(n, " ", names(n), collapse = ", "), " (sur ", sum(n), ") ; \u00e9cart le plus marqu\u00e9 : ", names(n)[i], ", ", n[i], " contre ", .gbay_fmt(e[i], ",", 1L), " attendus."))
		return(invisible(NULL))
	}
	t <- obj$table; e <- obj$expected
	k <- which.max(abs(t - e) / sqrt(e))
	i <- (k - 1L) %% nrow(t) + 1L; j <- (k - 1L) %/% nrow(t) + 1L
	emit("fact", paste0("Largest gap observed: cell (", rownames(t)[i], ", ", colnames(t)[j], "), ", t[i, j], " against ", .gbay_fmt(e[i, j], ".", 1L), " expected without a link."),
		paste0("\u00c9cart observ\u00e9 le plus marqu\u00e9 : case (", rownames(t)[i], ", ", colnames(t)[j], "), ", t[i, j], " contre ", .gbay_fmt(e[i, j], ",", 1L), " attendus sans lien."))
	invisible(NULL)
}

.gbay_print_detail <- function(obj) {
	say <- .bay_detail_head(obj)
	if (identical(obj$test, "gof")) {
		say("Prior: under a departure, Dirichlet(1, ..., 1) on the proportions (all equally possible); under none, the proportions given exactly. Closed-form Bayes factor (Good 1967; Sarafoglou et al. 2023).",
			"A priori : s'il y a un \u00e9cart, Dirichlet(1, ..., 1) sur les proportions (toutes \u00e9galement possibles) ; sinon, exactement les proportions donn\u00e9es. Facteur de Bayes en forme close (Good 1967 ; Sarafoglou et al. 2023).")
		say("Sensitivity: concentration a of the Dirichlet 1/2, 1 and 2, crossed with the starting chances 1:2, 1:1, 2:1.",
			"Sensibilit\u00e9 : concentration a de la Dirichlet 1/2, 1 et 2, crois\u00e9e avec les chances de d\u00e9part 1:2, 1:1, 2:1.")
		say("Estimation: posterior Dirichlet(n + 1); per category, ratio of the proportion to the one given, 95 % HDI (bayestestR); Cohen's w.",
			"Estimation : Dirichlet(n + 1) a posteriori ; par modalit\u00e9, rapport de la proportion \u00e0 celle donn\u00e9e, HDI \u00e0 95 % (bayestestR) ; w de Cohen.")
	} else {
		say("Prior: under a link, Dirichlet(2) on the cell probabilities; under independence, rows Dirichlet(2c) and columns Dirichlet(2r), the laws that a link gives to the margins (neutral margins; the BDe score of Heckerman et al. 1995). Closed form (Gunel and Dickey 1974; DECISIONS, G12).",
			"A priori : s'il y a un lien, Dirichlet(2) sur les probabilit\u00e9s des cases ; sous l'ind\u00e9pendance, lignes Dirichlet(2c) et colonnes Dirichlet(2r), les lois qu'un lien donne aux marges (marges neutres ; score BDe de Heckerman et al. 1995). Forme close (Gunel et Dickey 1974 ; DECISIONS, G12).")
		say("The same Bayes factor whether the total, the rows or the columns were fixed (the row margin cancels).",
			"Le m\u00eame facteur de Bayes que le total, les lignes ou les colonnes aient \u00e9t\u00e9 fix\u00e9s (la marge des lignes se simplifie).")
		say("a = 2 chosen in advance: the prior median of Cramer's V under a link is 0.2 to 0.4 depending on the table, a medium effect. Sensitivity: a = 1, 2, 4 crossed with the starting chances. Measures: dev/gbay/VALIDATION.md, sections 10 and 12.",
			"a = 2 choisi d'avance : sous un lien, le V de Cram\u00e9r m\u00e9dian a priori vaut 0,2 \u00e0 0,4 selon la table, un effet moyen. Sensibilit\u00e9 : a = 1, 2, 4 crois\u00e9 avec les chances de d\u00e9part. Mesures : dev/gbay/VALIDATION.md, sections 10 et 12.")
		say("Estimation (cells, V, odds ratio): uniform Dirichlet posterior on the cell probabilities, Dirichlet(n + 1).",
			"Estimation (cases, V, rapport des cotes) : loi a posteriori de Dirichlet uniforme sur les probabilit\u00e9s des cases, Dirichlet(n + 1).")
		say("Cells: ratio pi_ij / (pi_i+ pi_+j) on the posterior draws, 95 % HDI on its log; a cell is read when the interval excludes 1 (Albert 1997; Kruschke 2015, ch. 24).",
			"Cases : rapport pi_ij / (pi_i+ pi_+j) sur les tirages a posteriori, HDI \u00e0 95 % sur son log ; une case est lue quand l'intervalle exclut 1 (Albert 1997 ; Kruschke 2015, chap. 24).")
	}
	invisible(obj)
}

# ---------------------------------------------------------------------------------------------
# Whole analysis

.gbay_run <- function(inp, prior_h1, sensitivity, draws, seed, detail, verbose, return, story, note, get_notes, call_expr, correct_given) {
	gof <- identical(inp$kind, "gof")
	tab <- inp$tab
	N <- sum(tab)
	nm <- inp$names
	if (gof) {
		if (is.null(nm) || !length(nm)) nm <- "x"
		nm <- nm[1L]
	} else if (is.null(nm) || length(nm) < 2L) {
		nm <- c("rows", "columns")
	}
	nm_fr <- if (identical(nm, c("rows", "columns"))) c("lignes", "colonnes") else nm
	obj <- list(test = if (gof) "gof" else "independence", names = nm, resp = isTRUE(inp$resp))
	lab <- .gbay_labels(obj)
	if (!gof && identical(nm, c("rows", "columns"))) {
		lab <- list(h1_en = "Link between rows and columns", h1_fr = "Lien entre lignes et colonnes",
			h0_en = "Rows and columns independent", h0_fr = "Lignes et colonnes ind\u00e9pendantes")
	}

	# Step 1: the plan.
	if (gof) {
		.bay_tell(story, paste0("Plan recognised: goodness of fit of ", nm, " (", length(tab), " categories) to given proportions; N = ", N, "."),
			paste0("Plan reconnu : ajustement de ", nm, " (", length(tab), " modalit\u00e9s) \u00e0 des proportions donn\u00e9es ; N = ", N, "."), "step")
		pp_en <- paste0(names(tab), " ", vapply(inp$p, .bay_pct, character(1), dec = "."), collapse = ", ")
		pp_fr <- paste0(names(tab), " ", vapply(inp$p, .bay_pct, character(1), dec = ","), collapse = ", ")
		.bay_tell(story, paste0("Proportions expected", if (!isTRUE(inp$p_given)) " (equal, p not given)" else "", ": ", pp_en, "."),
			paste0("Proportions attendues", if (!isTRUE(inp$p_given)) " (\u00e9gales, p non donn\u00e9)" else "", " : ", pp_fr, "."))
		.bay_tell(story, paste0("One question: do the proportions of ", nm, " depart from those expected? (counterpart of the G test of goodness of fit)"),
			paste0("Une question : les proportions de ", nm, " s'\u00e9cartent-elles de celles attendues ? (pendant du test G d'ajustement)"), "decide")
	} else {
		dims <- paste0(nrow(tab), " x ", ncol(tab))
		unnamed <- identical(nm, c("rows", "columns"))
		.bay_tell(story, paste0("Plan recognised: independence between ", nm[1L], if (!unnamed) " (rows)", " and ", nm[2L], if (!unnamed) " (columns)", ", table ", dims, "; N = ", N, "."),
			paste0("Plan reconnu : ind\u00e9pendance entre ", nm_fr[1L], if (!unnamed) " (lignes)", " et ", nm_fr[2L], if (!unnamed) " (colonnes)", ", table ", dims, " ; N = ", N, "."), "step")
		form_en <- switch(inp$form, table = "a contingency table", xy = "two variables x and y", formula_two = "a formula with two variables",
			formula_resp = paste0("a response (", nm[2L], ") and a factor (", nm[1L], ")"), "a table")
		form_fr <- switch(inp$form, table = "une table de contingence", xy = "deux variables x et y", formula_two = "une formule \u00e0 deux variables",
			formula_resp = paste0("une r\u00e9ponse (", nm[2L], ") et un facteur (", nm[1L], ")"), "une table")
		.bay_tell(story, paste0("Input: ", form_en, if (!is.null(inp$weights)) paste0(", counts in ", inp$weights) else "", "."),
			paste0("Entr\u00e9e : ", form_fr, if (!is.null(inp$weights)) paste0(", effectifs dans ", inp$weights) else "", "."))
		.bay_tell(story, paste0("One question: are ", nm[1L], " and ", nm[2L], " linked? (counterpart of the G test of independence)"),
			paste0("Une question : ", nm_fr[1L], " et ", nm_fr[2L], " sont-ils li\u00e9s ? (pendant du test G d'ind\u00e9pendance)"), "decide")
	}
	for (ch in inp$changes) {
		note(paste0("Change of strategy: ", ch$en), paste0("Changement de strat\u00e9gie : ", ch$fr), tell = FALSE)
		.bay_tell(story, paste0("Change of strategy: ", ch$en), paste0("Changement de strat\u00e9gie : ", ch$fr), "decide")
	}

	# Step 2: checks.
	.bay_tell(story, "Checks on the data", "Contr\u00f4les des donn\u00e9es", "step")
	.bay_step_result(story, "exact computation, no minimum count", "calcul exact, sans effectif minimum")
	if (inp$n_removed > 0L) {
		note(paste0(inp$n_removed, " observation(s) with missing values removed."), paste0(inp$n_removed, " observation(s) avec valeurs manquantes retir\u00e9e(s)."))
	}
	if (gof) {
		E <- N * inp$p
		.bay_tell(story, paste0("Counts: ", paste0(tab, " ", names(tab), collapse = ", "), " (of ", N, ")."),
			paste0("Effectifs : ", paste0(tab, " ", names(tab), collapse = ", "), " (sur ", N, ")."))
	} else {
		E <- outer(rowSums(tab), colSums(tab)) / N
		show <- seq_len(min(nrow(tab), 8L))
		for (i in show) {
			.bay_tell(story, paste0(rownames(tab)[i], ": ", paste0(tab[i, ], " ", colnames(tab), collapse = ", "), " (of ", sum(tab[i, ]), ")."),
				paste0(rownames(tab)[i], " : ", paste0(tab[i, ], " ", colnames(tab), collapse = ", "), " (sur ", sum(tab[i, ]), ")."))
		}
		if (nrow(tab) > 8L) .bay_tell(story, paste0("... ", nrow(tab) - 8L, " more rows in $table."), paste0("... ", nrow(tab) - 8L, " lignes de plus dans $table."))
	}
	small <- sum(E < 5)
	if (small > 0L) {
		.bay_tell(story, paste0(small, " of ", length(E), if (gof) " categories" else " cells", " with an expected count under 5: the G test and the chi-squared would be approximate there (Williams' correction)."),
			paste0(small, " ", if (gof) "modalit\u00e9" else "case", if (small > 1L) "s" else "", " sur ", length(E), " d'effectif attendu inf\u00e9rieur \u00e0 5 : le test G et le khi-deux y seraient approximatifs (correction de Williams)."))
	}
	if (gof || small == 0L) {
		.bay_tell(story, "small counts do not make the answer wrong, only less decisive.", "les petits effectifs ne faussent pas la r\u00e9ponse, ils la rendent moins tranch\u00e9e.", "note")
	} else {
		.bay_tell(story, "small counts make the answer less decisive: the Bayes factor stays near 1 or leans towards independence.",
			"les petits effectifs rendent la r\u00e9ponse moins tranch\u00e9e : le facteur de Bayes reste pr\u00e8s de 1 ou penche vers l'ind\u00e9pendance.", "note")
	}
	if (correct_given) {
		.bay_tell(story, "correct (Williams, Yates) is ignored: these corrections bring the p-value of G closer to its law; the Bayes factor is read on G itself.",
			"correct (Williams, Yates) est ignor\u00e9 : ces corrections rapprochent la p-value de G de sa loi ; le facteur de Bayes se lit sur G lui-m\u00eame.", "note_first")
		note("correct is ignored.", "correct est ignor\u00e9.", tell = FALSE)
	}

	# Step 3: law of the data.
	if (gof) {
		.bay_tell(story, "Law of the data", "Loi des donn\u00e9es", "step")
		.bay_step_result(story, "multinomial, total N fixed", "multinomiale, total N fix\u00e9")
		.bay_tell(story, paste0("Departure: the ", length(tab), " proportions are free, every set equally possible in advance."),
			paste0("\u00c9cart : les ", length(tab), " proportions sont libres, toutes \u00e9galement possibles d'avance."))
		.bay_tell(story, "No departure: the proportions are exactly those expected.", "Pas d'\u00e9cart : les proportions sont exactement celles attendues.")
	} else {
		.bay_tell(story, "Law of the data", "Loi des donn\u00e9es", "step")
		.bay_step_result(story, "multinomial, margins neutral", "multinomiale, marges neutres")
		unnamed_t <- identical(nm, c("rows", "columns"))
		.bay_tell(story, if (unnamed_t) "Link: each row its own proportions; independence: the same proportions in every row." else
				paste0("Link: each level of ", nm[1L], " its own proportions of ", nm[2L], "; independence: the same proportions for every level."),
			if (unnamed_t) "Lien : chaque ligne a ses propres proportions ; ind\u00e9pendance : les m\u00eames proportions dans toutes les lignes." else
				paste0("Lien : chaque modalit\u00e9 de ", nm_fr[1L], " a ses propres proportions de ", nm_fr[2L], " ; ind\u00e9pendance : les m\u00eames proportions pour toutes."))
		.bay_tell(story, "The same answer whether the total, the rows or the columns were fixed when the data were collected.",
			"M\u00eame r\u00e9ponse que le total, les lignes ou les colonnes aient \u00e9t\u00e9 fix\u00e9s au recueil des donn\u00e9es.", "decide")
		.bay_tell(story, "the margins follow the same law under both hypotheses: balanced margins say nothing of a link.",
			"les marges suivent la m\u00eame loi sous les deux hypoth\u00e8ses : des marges \u00e9quilibr\u00e9es ne disent rien d'un lien.", "note")
	}

	# Step 4: priors.
	.bay_tell(story, "Starting point (priors)", "Point de d\u00e9part (a priori)", "step")
	.bay_step_result(story, paste0(.bay_pct(prior_h1, "."), " chance of ", if (gof) "a departure" else "a link"),
		paste0(.bay_pct(prior_h1, ","), " de chances d'", if (gof) "un \u00e9cart" else "un lien"))
	if (gof) {
		.bay_tell(story, "If there is a departure, every set of proportions equally possible in advance.",
			"S'il y a un \u00e9cart, toutes les proportions sont \u00e9galement possibles d'avance.")
		.bay_tell(story, "Prior: Dirichlet(1, ..., 1) under a departure, the proportions given under none (Good 1967; Sarafoglou et al. 2023).",
			"A priori : Dirichlet(1, ..., 1) s'il y a un \u00e9cart, les proportions donn\u00e9es sinon (Good 1967 ; Sarafoglou et al. 2023).", "detail")
	} else {
		.bay_tell(story, "If there is a link, a link of medium size is expected (prior median of Cramer's V about 0.3).",
			"S'il y a un lien, on attend un lien de taille moyenne (V de Cram\u00e9r m\u00e9dian a priori vers 0,3).")
		.bay_tell(story, "Prior: Dirichlet(2) per row under a link, Dirichlet(2r) common under independence (neutral margins, BDe score; DECISIONS, G12).",
			"A priori : Dirichlet(2) par ligne s'il y a un lien, Dirichlet(2r) commune sous l'ind\u00e9pendance (marges neutres, score BDe ; DECISIONS, G12).", "detail")
	}
	.bay_tell(story, "these values are not tuned on the data; step 5 varies them.",
		"ces valeurs ne sont pas ajust\u00e9es sur les donn\u00e9es ; l'\u00e9tape 5 les fait varier.", "note")

	# Step 5: Bayes factor and sensitivity.
	ll <- NULL
	if (gof) {
		.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = c("g.bay() {KefiR}", "g.bay() {KefiR}"))
		lb <- .gbay_lbf_gof(tab, inp$p, 1)
	} else {
		.bay_tell(story, "Solidity of the answer", "Solidit\u00e9 de la r\u00e9ponse", "step", fun = c("g.bay() {KefiR}", "g.bay() {KefiR}"))
		ll <- .gbay_G(tab)
		lb <- .gbay_lbf_ct(tab)
	}
	p_h1 <- .gbay_prob(lb, prior_h1)
	inv_en <- if (lb < 0) paste0(", i.e. ", .bay_bf_txt(-lb, "."), " in favour of ", if (gof) "the proportions given" else "independence") else ""
	inv_fr <- if (lb < 0) paste0(", soit ", .bay_bf_txt(-lb, ","), " en faveur ", if (gof) "des proportions donn\u00e9es" else "de l'ind\u00e9pendance") else ""
	.bay_tell(story, paste0("Bayes factor: ", .bay_bf_txt(lb, "."), inv_en, "."), paste0("Facteur de Bayes : ", .bay_bf_txt(lb, ","), inv_fr, "."))
	.bay_tell(story, paste0("The Bayes factor says how many times more probable the counts are with ", if (gof) "a departure" else "a link", " than without."),
		paste0("Le facteur de Bayes dit combien de fois les effectifs sont plus probables avec ", if (gof) "un \u00e9cart" else "un lien", " que sans."), "note")
	.bay_tell(story, "Precision: exact computation, without draws.", "Pr\u00e9cision : calcul exact, sans tirage.")
	sens <- NULL
	if (isTRUE(sensitivity)) {
		concs <- if (gof) c(1 / 2, 1, 2) else .bay_dm_neutral_grid()
		phs <- c(1 / 3, 1 / 2, 2 / 3)
		lbs <- vapply(concs, function(a) if (gof) .gbay_lbf_gof(tab, inp$p, a) else .gbay_lbf_ct(tab, a), numeric(1))
		grid <- outer(lbs, phs, function(l, p) .gbay_prob(l, p))
		lead_h1 <- p_h1 >= 0.5
		rng <- range(if (lead_h1) grid else 1 - grid)
		verdict <- .bay_sensitivity_verdict(max(p_h1, 1 - p_h1), rng)
		sens <- list(conc = concs, prior_h1 = phs, p = grid, log_bf = lbs, range_lead = rng, verdict = verdict)
		stb <- .bay_stab_txt(verdict)
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", ", stb$en), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", ", stb$fr))
		.bay_tell(story, paste0("Sensitivity: the leading answer ", .bay_range_txt(rng)$en, "."),
			paste0("Sensibilit\u00e9 : la r\u00e9ponse en t\u00eate ", .bay_range_txt(rng)$fr, "."))
		.bay_tell(story, if (gof) "what varies is the expected spread of the proportions (a = 1/2, 1, 2) and the starting chances (1 against 2 to 2 against 1)." else
				"what varies is the expected size of a link (a = 1, 2, 4) and the starting chances (1 against 2 to 2 against 1).",
			if (gof) "on fait varier l'\u00e9talement attendu des proportions (a = 1/2, 1, 2) et les chances de d\u00e9part (1 contre 2 \u00e0 2 contre 1)." else
				"on fait varier la taille attendue d'un lien (a = 1, 2, 4) et les chances de d\u00e9part (1 contre 2 \u00e0 2 contre 1).", "note")
		if (identical(verdict, "changes")) {
			.bay_tell(story, "Caution added: the answer can reverse; it will be called not settled.",
				"Prudence ajout\u00e9e : la r\u00e9ponse peut s'inverser ; elle sera dite non tranch\u00e9e.", "decide")
		} else {
			.bay_tell(story, "The answer does not change side.", "La r\u00e9ponse ne change pas de c\u00f4t\u00e9.", "decide")
		}
	} else {
		.bay_step_result(story, paste0("Bayes factor ", .bay_bf_txt(lb, "."), ", sensitivity not computed"), paste0("facteur de Bayes ", .bay_bf_txt(lb, ","), ", sensibilit\u00e9 non calcul\u00e9e"))
		.bay_tell(story, "Sensitivity to the priors not computed (sensitivity = FALSE).", "Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	hyp <- data.frame(label_en = c(lab$h1_en, lab$h0_en), label_fr = c(lab$h1_fr, lab$h0_fr),
		posterior = c(p_h1, 1 - p_h1), prior = c(prior_h1, 1 - prior_h1), null = c(FALSE, TRUE), stringsAsFactors = FALSE)
	hyp <- hyp[order(-hyp$posterior), , drop = FALSE]
	rownames(hyp) <- NULL
	if (isFALSE(return)) return(p_h1)

	# Step 6: where the gap lies, and its size.
	est <- if (gof) .gbay_summ_gof(tab, inp$p, 1, max(1000L, as.integer(draws)), seed) else
		.gbay_summ_ct(.gbay_draws_ct(tab, max(1000L, as.integer(draws)), seed), tab)
	out <- list(hypotheses = hyp, prob_effect = p_h1, bf10 = exp(lb), log_bf = lb,
		effects = data.frame(term = if (gof) nm else paste0(nm[1L], ":", nm[2L]), bf10 = exp(lb), log_bf = lb, p_h1 = p_h1, p_h0 = 1 - p_h1,
			prior_h1 = prior_h1, analogue = if (gof) "G test of goodness of fit" else "G test of independence", stringsAsFactors = FALSE),
		test = obj$test, names = nm, names_fr = nm_fr, resp = obj$resp, G = if (gof) NULL else ll$G, df = if (gof) length(tab) - 1L else ll$df,
		table = tab, expected = E, p = if (gof) inp$p else NULL, cells = est$cells,
		cramer_v = est$cramer_v, odds_ratio = est$odds_ratio, w = est$w,
		sensitivity = sens, contradiction = FALSE, detail = isTRUE(detail),
		priors = list(concentration = if (gof) 1 else .bay_dm_neutral_a(), prior_h1 = prior_h1, law = if (gof) "multinomial" else "Dirichlet-multinomial, neutral margins"),
		method = if (gof) "Goodness of fit, Dirichlet-multinomial against a point (Good 1967)" else
			"Independence, Dirichlet-multinomial with neutral margins (BDe), a = 2 (G12)",
		label = if (gof) nm else paste(nm, collapse = " x "), n = N, call = call_expr)
	out$label <- if (gof) nm else paste0(nm[1L], " x ", nm[2L])
	pk <- if (gof) c("g.bay() {KefiR}", "hdi() {bayestestR}") else c("posterior() {BayesFactor}", "hdi() {bayestestR}")
	.bay_tell(story, if (gof) "Where the departure lies, and its size" else "Where the link lies, and its size",
		if (gof) "O\u00f9 se trouve l'\u00e9cart, et sa taille" else "O\u00f9 se trouve le lien, et sa taille", "step", fun = pk[1L])
	.bay_step_fun(story, pk[2L])
	if (gof) {
		.bay_step_result(story, paste0("Cohen's w ", .gbay_ci(est$w, ".", function(x, dc) .gbay_fmt(x, dc))), paste0("w de Cohen ", .gbay_ci(est$w, ",", function(x, dc) .gbay_fmt(x, dc))))
	} else {
		.bay_step_result(story, paste0("Cramer's V ", .gbay_ci(est$cramer_v, ".", function(x, dc) .gbay_fmt(x, dc))), paste0("V de Cram\u00e9r ", .gbay_ci(est$cramer_v, ",", function(x, dc) .gbay_fmt(x, dc))))
	}
	tell_emit <- function(kind, en, fr) .bay_tell(story, en, fr, kind)
	if (!is.null(est$odds_ratio)) .gbay_size_line(out, tell_emit)
	if (p_h1 < 0.5) {
		.bay_tell(story, if (gof) "No departure is in the lead: the lines below are read only as an indication." else "Independence is in the lead: the lines below are read only as an indication.",
			if (gof) "L'absence d'\u00e9cart est en t\u00eate : les lignes ci-dessous ne se lisent qu'\u00e0 titre indicatif." else "L'ind\u00e9pendance est en t\u00eate : les lignes ci-dessous ne se lisent qu'\u00e0 titre indicatif.", "decide")
	}
	if (gof && length(tab) == 2L) {
		ce <- est$cells
		pc <- function(v, dc) gsub(".", dc, formatC(100 * v, format = "f", digits = 1), fixed = TRUE)
		.bay_tell(story, paste0("Share of ", ce$category[1L], ": ", pc(ce$prop[1L], "."), " % [", pc(ce$prop_low[1L], "."), "; ", pc(ce$prop_high[1L], "."), "], expected ", pc(ce$p[1L], "."), " %."),
			paste0("Part de ", ce$category[1L], " : ", pc(ce$prop[1L], ","), " % [", pc(ce$prop_low[1L], ","), " ; ", pc(ce$prop_high[1L], ","), "], attendue ", pc(ce$p[1L], ","), " %."))
	}
	.gbay_cell_lines(out, tell_emit, top = 4L)
	.bay_tell(story, if (gof) "ratio = proportion estimated / proportion given, read when its 95 % interval excludes 1." else "ratio = cell / its value without a link, read when its 95 % interval excludes 1.",
		if (gof) "rapport = proportion estim\u00e9e / proportion donn\u00e9e, lu quand son intervalle \u00e0 95 % exclut 1." else "rapport = case / sa valeur sans lien, lu quand son intervalle \u00e0 95 % exclut 1.", "note_first")
	.bay_tell(story, if (gof) "w is positive even without a departure: a size, not a proof (the proof is the Bayes factor)." else "V is positive even without a link: a size, not a proof (the proof is the Bayes factor).",
		if (gof) "w est positif m\u00eame sans \u00e9cart : une taille, pas une preuve (la preuve est le facteur de Bayes)." else "V est positif m\u00eame sans lien : une taille, pas une preuve (la preuve est le facteur de Bayes).", "note")

	out$notes <- get_notes()
	out$guard <- .gbay_guard(out)$on
	class(out) <- "gbay"
	out$story <- .bay_story_lines(story)
	if (isTRUE(verbose)) .gbay_print_bilan(out, detail = detail, narrated = TRUE)
	if (isTRUE(verbose)) invisible(out) else out
}
