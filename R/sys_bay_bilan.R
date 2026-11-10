# sys_bay_bilan.R
# The bilan of the Bayesian family: hypotheses ranked with their posterior
# probability, a plain qualifier, the answers per term, the priors, the sensitivity,
# and the difference between a probability, a Bayes factor and a p-value. Reusable
# by cor.bay and g.bay. Nothing is exported except the print method of m.bay.
#
# References (comments only: dev/bp.log, 7.4 and BP-008).
#   Kass, R. E., & Raftery, A. E. (1995). Bayes factors. JASA, 90(430), 773-795.
#     https://doi.org/10.1080/01621459.1995.10476572 (thresholds 3, 20, 150 on the
#     BF; at even odds they read 0.75, 0.95, 0.993 on the posterior probability.)

# Qualifier of the support for the leading hypothesis. Never called "confidence", so
# that it cannot be mistaken for a confidence interval or a frequentist threshold.
.bay_support_breaks <- function() c(0.6, 0.75, 0.95, 0.993)

# Qualifier of the support when the data bring little: with a Bayes factor between 1/3 and 3
# ("not worth more than a bare mention", Kass and Raftery 1995), the data hardly move the odds
# and the answer is "not settled by the data", whatever probability the prior gives; the
# probability itself is still shown (DECISIONS, D14).
.bay_support_bf <- function(p, log_bf) {
	if (is.na(log_bf) || abs(log_bf) >= log(3) || p < 0.6) return(.bay_support(p))
	list(en = "not settled by the data", fr = "non tranch\u00e9 par les donn\u00e9es", level = 1L)
}

.bay_support <- function(p) {
	breaks <- .bay_support_breaks()
	idx <- findInterval(p, breaks) + 1L
	en <- c("no support, close to a coin toss", "weak support", "moderate support", "strong support", "very strong support")
	fr <- c("aucun soutien, proche de pile ou face", "soutien faible", "soutien mod\u00e9r\u00e9", "soutien fort", "soutien tr\u00e8s fort")
	list(en = en[idx], fr = fr[idx], level = idx)
}

# Percentages: never "100 %" or "0 %" for a probability that is not exactly so.
.bay_pct <- function(p, dec = ".") {
	if (is.na(p)) return("NA")
	txt <- if (p > 0.999) {
		"> 99.9"
	} else if (p < 0.001) {
		"< 0.1"
	} else {
		format(round(100 * p, 1), nsmall = 1, trim = TRUE)
	}
	paste0(if (dec == ",") gsub(".", ",", txt, fixed = TRUE) else txt, " %")
}

# Bayes factor as text, from its log.
.bay_bf_txt <- function(log_bf, dec = ".") {
	if (is.na(log_bf)) return("NA")
	txt <- if (abs(log_bf) > log(1e4)) {
		paste0("10^", round(log_bf / log(10)))
	} else {
		format(signif(exp(log_bf), 3), trim = TRUE, scientific = FALSE)
	}
	if (dec == ",") gsub(".", ",", txt, fixed = TRUE) else txt
}

# What the sensitivity says about the conclusion, given the central probability of
# the hypothesis in the lead and the range of that same probability over the priors.
#   "changes"  : the lead changes somewhere: "the data do not settle it".
#   "crosses"  : the lead stays, a threshold of the qualifier is crossed.
#   "stable"   : neither.
.bay_sensitivity_verdict <- function(p_lead, range_lead) {
	if (anyNA(range_lead)) return("changes")
	if (min(range_lead) <= 1/2) return("changes")
	a <- .bay_support(min(range_lead))
	b <- .bay_support(max(range_lead))
	if (a$level != b$level) "crosses" else "stable"
}

# The sentence a reader looks for first. Never an equivalence, a size or a cause.
.bay_verdict <- function(obj) {
	v <- .bay_verdict_raw(obj)
	if (!is.null(obj$posthoc) && (isTRUE(obj$posthoc$k > 2L) || identical(obj$posthoc$method, "slopes"))) {
		v$en <- sub("neither that it matters, nor which groups differ, nor why.", "neither that it matters nor why; the pairs are below.", v$en, fixed = TRUE)
		v$fr <- sub("ni qu'elle est importante, ni quels groupes diff\u00e8rent, ni pourquoi.", "ni qu'elle est importante, ni pourquoi ; les paires sont plus bas.", v$fr, fixed = TRUE)
	}
	v
}

.bay_verdict_raw <- function(obj) {
	h <- obj$hypotheses
	p <- h$posterior[1L]
	# A probability that comes from the prior rather than the data does not make a firm
	# verdict (Bayes factor of the leading answer between 1/3 and 3; DECISIONS, D14).
	lb <- log(h$posterior[1L] / h$posterior[2L]) - log(h$prior[1L] / h$prior[2L])
	if (isTRUE(abs(lb) < log(3))) p <- min(p, 0.59)
	lead_h1 <- isTRUE(obj$prob_effect >= 0.5)
	open <- identical(obj$sensitivity$verdict, "changes")
	if (isTRUE(obj$contradiction)) {
		return(list(en = "Verdict: undecided, the answer depends on the shape of the noise.",
			fr = "Verdict : ind\u00e9cis, la r\u00e9ponse d\u00e9pend de la forme du bruit."))
	}
	if (p < 0.6) {
		return(list(en = "Verdict: undecided. The data allow no conclusion, neither a difference nor its absence.",
			fr = "Verdict : ind\u00e9cis. Les donn\u00e9es ne permettent de conclure ni \u00e0 une diff\u00e9rence, ni \u00e0 son absence."))
	}
	if (!lead_h1) {
		if (p >= 0.95 && !open) {
			return(list(en = "Verdict: very probably no difference of the size expected. This is not a proof of equivalence.",
				fr = "Verdict : tr\u00e8s probablement pas de diff\u00e9rence de l'ordre attendu. Ce n'est pas une preuve d'\u00e9quivalence."))
		}
		return(list(en = "Verdict: rather no difference, without concluding. This does not prove that the groups are equivalent.",
			fr = "Verdict : plut\u00f4t pas de diff\u00e9rence, sans conclure. Cela ne prouve pas que les groupes sont \u00e9quivalents."))
	}
	if (p >= 0.95 && !open) {
		return(list(en = "Verdict: a difference very probably exists. This says neither that it matters, nor which groups differ, nor why.",
			fr = "Verdict : une diff\u00e9rence existe tr\u00e8s probablement. Cela ne dit ni qu'elle est importante, ni quels groupes diff\u00e8rent, ni pourquoi."))
	}
	if (p >= 0.75 && !open) {
		return(list(en = "Verdict: a difference probably exists. This says neither that it matters, nor which groups differ, nor why.",
			fr = "Verdict : une diff\u00e9rence existe probablement. Cela ne dit ni qu'elle est importante, ni quels groupes diff\u00e8rent, ni pourquoi."))
	}
	list(en = "Verdict: rather a difference, without concluding.",
		fr = "Verdict : plut\u00f4t une diff\u00e9rence, sans conclure.")
}

.bay_law_labels <- function() {
	list(
		en = c(normal_eq = "normal, equal variances", normal_het = "normal, unequal variances",
			t_eq = "heavy tails, equal variances", t_het = "heavy tails, unequal variances"),
		fr = c(normal_eq = "normale, variances \u00e9gales", normal_het = "normale, variances in\u00e9gales",
			t_eq = "queues lourdes, variances \u00e9gales", t_het = "queues lourdes, variances in\u00e9gales")
	)
}

# Prints the bilan. obj is the list returned by m.bay(). Text goes through .msg: one
# string in each language, numbers with the decimal mark of that language. The default is
# short enough for a student or an engineer; detail = TRUE adds the shapes of the noise,
# the varied priors and how to read the numbers.
.bay_print_bilan <- function(obj, detail = FALSE, narrated = FALSE) {
	sens <- obj$sensitivity
	p <- obj$priors
	dec <- if (.bay_lang_fr()) "," else "."
	num <- function(x) gsub(".", dec, format(signif(x, 3), trim = TRUE), fixed = TRUE)
	body <- function(emit, li, narrated) {
		d <- obj$descriptive_stats
		dg <- .bay_dg(d$mean)
		ph <- obj$posthoc
		undecided <- li$h$posterior[1L] < 0.6 || li$open
		if (!is.null(ph) && ph$k > 2L) .bay_pairs_short(ph, emit, undecided, narrated)
		if (!is.null(ph) && ph$k == 2L) {
			pr <- ph$pairs
			ci <- function(dc) paste0(.bay_fx(pr$diff, dg, dc), " [", .bay_fx(pr$diff_low, dg, dc), if (dc == ",") " ; " else "; ", .bay_fx(pr$diff_high, dg, dc), "]")
			emit("fact", paste0("Estimated gap ", pr$group2, " - ", pr$group1, ": ", ci("."), if (narrated) " (step 6)." else "."),
				paste0("\u00c9cart estim\u00e9 ", pr$group2, " - ", pr$group1, " : ", ci(","), if (narrated) " (\u00e9tape 6)." else "."))
		}
		centres <- isTRUE(obj$heavy) || isTRUE(obj$contradiction) || isTRUE(obj$p_tails > 0.5)
		val <- if (centres) d$median else d$mean
		same_n <- length(unique(d$n)) == 1L
		obs <- function(dc) paste0(paste0(d$group, " = ", .bay_fx(val, dg, dc), if (!same_n) paste0(" (n = ", d$n, ")") else "", collapse = if (dc == ",") " ; " else "; "),
			if (same_n) paste0(" (n = ", d$n[1L], if (dc == ",") " par groupe)" else " per group)") else "")
		emit("fact", paste0(if (centres) "Observed medians: " else "Observed means: ", obs("."), "."),
			paste0(if (centres) "M\u00e9dianes observ\u00e9es : " else "Moyennes observ\u00e9es : ", obs(","), "."))
		if (isTRUE(obj$contradiction)) {
			emit("note", "centres and means do not tell the same story (step 5): support capped at weak.",
				"centres et moyennes ne racontent pas la m\u00eame chose (\u00e9tape 5) : soutien plafonn\u00e9 \u00e0 faible.")
		}
		if (!is.null(ph) && isTRUE(ph$caution)) {
			emit("note", "the pairs assume a normal noise with a common variance, which step 3 did not find: read them with caution.",
				"les paires supposent un bruit normal \u00e0 variance commune, ce que l'\u00e9tape 3 n'a pas trouv\u00e9 : les lire avec r\u00e9serve.")
		}
	}
	.bay_bilan_common(obj, narrated, body)
	if (!isTRUE(detail)) return(invisible(obj))
	say <- .bay_detail_head(obj)

	# Detail.
	if (!is.null(sens)) {
		tip <- log(sens$odds_to_reverse)
		say(paste0("To level the two hypotheses, one would have to believe in the other one in advance at ",
				.bay_bf_txt(tip, "."), " to 1."),
			paste0("Pour mettre les deux hypoth\u00e8ses \u00e0 \u00e9galit\u00e9, il faudrait croire d'avance \u00e0 l'autre \u00e0 ",
				.bay_bf_txt(tip, ","), " contre 1."))
	}
	if (!is.null(obj$laws)) {
		lb <- .bay_law_labels()
		laws <- obj$laws[order(-obj$laws$posterior), , drop = FALSE]
		if (nrow(laws) > 1L) {
			say("Shape of the noise, compared and not chosen (posterior probabilities):",
				"Forme du bruit, compar\u00e9e et non choisie (probabilit\u00e9s a posteriori) :")
		} else {
			say("Shape of the noise: one law only, not compared:",
				"Forme du bruit : une seule loi retenue, non compar\u00e9e :")
		}
		say(paste0("  ", paste0(lb$en[laws$law], " ", vapply(laws$posterior, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
			paste0("  ", paste0(lb$fr[laws$law], " ", vapply(laws$posterior, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))
		if (!is.null(sens) && nrow(obj$laws) > 1L) {
			bl <- sens$by_law[laws$law]
			say(paste0("  If one shape were imposed, the probability of a difference would be: ",
					paste0(lb$en[laws$law], " ", vapply(bl, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
				paste0("  Si l'on imposait une seule forme, la probabilit\u00e9 d'une diff\u00e9rence serait : ",
					paste0(lb$fr[laws$law], " ", vapply(bl, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))
		}
	}
	say("Priors in detail:", "A priori en d\u00e9tail :")
	if (isTRUE(nrow(obj$descriptive_stats) > 2L)) {
		say("  Forms of a difference, at even odds: all groups may differ; or one group alone stands apart (scale r x sqrt(2)).",
			"  Formes d'une diff\u00e9rence, \u00e0 poids \u00e9gaux : tous les groupes peuvent diff\u00e9rer ; ou un seul groupe s'\u00e9carte (\u00e9chelle r x racine de 2).")
	}
	if (!is.null(obj$laws)) {
		kept <- obj$laws$law
		extra_en <- c(if (any(kept %in% c("normal_het", "t_het"))) paste0("proportions of precision Dirichlet(", format(p$u, trim = TRUE), ")"),
			if (any(kept %in% c("t_eq", "t_het"))) "heavy tails nu - 2 ~ Exp(1)")
		extra_fr <- c(if (any(kept %in% c("normal_het", "t_het"))) paste0("proportions de pr\u00e9cision Dirichlet(", num(p$u), ")"),
			if (any(kept %in% c("t_eq", "t_het"))) "queues lourdes nu - 2 ~ Exp(1)")
		head_en <- if (length(kept) > 1L) paste0("each of the ", length(kept), " shapes starts equal") else "one shape only"
		head_fr <- if (length(kept) > 1L) paste0("chacune des ", length(kept), " formes part \u00e0 \u00e9galit\u00e9") else "une seule forme"
		say(paste0("  Noise: ", paste(c(head_en, extra_en), collapse = "; "), "."),
			paste0("  Bruit : ", paste(c(head_fr, extra_fr), collapse = " ; "), "."))
	}
	if (!is.null(sens)) {
		say(paste0("  Varied in the sensitivity: r from ", format(min(sens$r), trim = TRUE), " to ",
				format(max(sens$r), trim = TRUE), ", odds H1:H0 ", paste(sens$odds[c(1L, length(sens$odds))], collapse = " to "),
				if (any(obj$laws$law %in% c("normal_het", "t_het"))) paste0(", Dirichlet ",
					format(min(sens$u), trim = TRUE), " to ", format(max(sens$u), trim = TRUE)), "."),
			paste0("  Vari\u00e9s dans la sensibilit\u00e9 : r de ", num(min(sens$r)), " \u00e0 ",
				num(max(sens$r)), ", cotes H1:H0 ", paste(sens$odds[c(1L, length(sens$odds))], collapse = " \u00e0 "),
				if (any(obj$laws$law %in% c("normal_het", "t_het"))) paste0(", Dirichlet ",
					num(min(sens$u)), " \u00e0 ", num(max(sens$u))), "."))
		if (isTRUE(nrow(obj$descriptive_stats) > 2L)) {
			say(paste0("  With the first form alone, the probability of a difference would be ", .bay_pct(sens$all_groups, "."), "."),
				paste0("  Avec la premi\u00e8re forme seule, la probabilit\u00e9 d'une diff\u00e9rence serait ", .bay_pct(sens$all_groups, ","), "."))
		}
	}
	breaks <- .bay_support_breaks()
	say("How to read it:", "Comment lire :")
	say("  The Bayes factor says how much the data move the odds. It does not depend on the prior odds.",
		"  Le facteur de Bayes dit de combien les donn\u00e9es d\u00e9placent les cotes. Il ne d\u00e9pend pas des cotes a priori.")
	say("  A p-value is the probability of data at least this extreme if there were no difference: neither of the above.",
		"  Une p-value est la probabilit\u00e9 de donn\u00e9es au moins aussi extr\u00eames s'il n'y avait pas de diff\u00e9rence : ni l'un ni l'autre ci-dessus.")
	say("  A difference between groups observed without randomisation does not show a cause.",
		"  Une diff\u00e9rence entre groupes observ\u00e9s sans tirage au sort ne montre pas une cause.")
	say(paste0("  Support: none under ", 100 * breaks[1L], " %, weak under ", 100 * breaks[2L], " %, moderate under ", 100 * breaks[3L],
			" %, strong under ", 100 * breaks[4L], " %, then very strong."),
		paste0("  Soutien : aucun sous ", 100 * breaks[1L], " %, faible sous ", 100 * breaks[2L], " %, mod\u00e9r\u00e9 sous ", 100 * breaks[3L],
			" %, fort sous ", gsub(".", ",", 100 * breaks[4L], fixed = TRUE), " %, puis tr\u00e8s fort."))
	invisible(obj)
}

# What each hypothesis rests on, one factor: the largest observed gap, in units of the
# within-group spread, and how much more probable it is under one hypothesis than under
# the other (the Bayes factor, against a difference of the size expected, r). Rows in the
# order of obj$hypotheses. The gap is given in full on the first row only.
.bay_hyp_args_oneway <- function(obj) {
	h <- obj$hypotheses
	if (is.null(h$null)) return(NULL)
	d <- obj$descriptive_stats
	e <- obj$effects
	lb <- e$log_bf[1L]
	if (is.na(lb)) return(NULL)
	hi <- which.max(d$mean)
	lo <- which.min(d$mean)
	gap <- d$mean[hi] - d$mean[lo]
	dg <- max(0L, 2L - floor(log10(max(abs(d$mean), 1e-300))))
	g_en <- formatC(gap, format = "f", digits = dg)
	sd_en <- formatC(obj$largest_gap, format = "f", digits = 2)
	g_fr <- gsub(".", ",", g_en, fixed = TRUE)
	sd_fr <- gsub(".", ",", sd_en, fixed = TRUE)
	r_en <- format(signif(obj$priors$rscale, 3), trim = TRUE)
	r_fr <- gsub(".", ",", r_en, fixed = TRUE)
	centres <- isTRUE(obj$heavy) || isTRUE(obj$contradiction) || isTRUE(obj$p_tails > 0.5)
	# Under heavy tails the groups are compared through their centres: the gap is then shown
	# on the medians, closer to them, with the gap of the means beside it.
	if (centres) {
		hi <- which.max(d$median)
		lo <- which.min(d$median)
		g_en <- formatC(d$median[hi] - d$median[lo], format = "f", digits = dg)
		g_fr <- gsub(".", ",", g_en, fixed = TRUE)
		gm_en <- formatC(diff(range(d$mean)), format = "f", digits = dg)
		gm_fr <- gsub(".", ",", gm_en, fixed = TRUE)
	}
	k <- nrow(d)
	what_en <- if (centres) "medians, close to the centres compared" else "means"
	small_en <- "the smallness of "
	small_fr <- "la petitesse de "
	what_fr <- if (centres) "m\u00e9dianes, proches des centres compar\u00e9s" else "moyennes"
	full_en <- paste0(if (k > 2L) "all the gaps between groups, the largest between " else "the gap between ", what_en, " (observed: ",
		d$group[hi], " - ", d$group[lo], " = ", g_en, if (!centres) paste0(", i.e. ", sd_en, " pooled within-group SD") else paste0("; gap of the means ", gm_en), ")")
	full_fr <- paste0(if (k > 2L) "l'ensemble des \u00e9carts entre groupes, le plus grand entre " else "l'\u00e9cart entre ", what_fr, " (observ\u00e9 : ",
		d$group[hi], " - ", d$group[lo], " = ", g_fr, if (!centres) paste0(", soit ", sd_fr, " \u00e9cart-type commun au sein des groupes") else paste0(" ; \u00e9cart des moyennes ", gm_fr), ")")
	bf_en <- .bay_bf_txt(abs(lb), ".")
	bf_fr <- .bay_bf_txt(abs(lb), ",")
	weak <- abs(lb) < log(3)
	more_en <- if (weak) paste0("hardly more probable (", bf_en, " times)") else paste0(bf_en, " times more probable")
	more_fr <- if (weak) paste0("\u00e0 peine plus probables (", bf_fr, " fois)") else paste0(bf_fr, " fois plus probables")
	en <- fr <- character(nrow(h))
	for (i in seq_len(nrow(h))) {
		gap_en <- if (i == 1L) full_en else "these data"
		gap_fr <- if (i == 1L) full_fr else "ces donn\u00e9es"
		for_it <- xor(isTRUE(h$null[i]), lb >= 0)
		if (!isTRUE(h$null[i])) {
			if (for_it) {
				en[i] <- paste0("Rests on ", gap_en, ": the data are ", more_en, " with a difference than without.")
				fr[i] <- paste0("Repose sur ", gap_fr, " : les donn\u00e9es sont ", more_fr, " avec une diff\u00e9rence que sans.")
			} else {
				en[i] <- paste0("Against it: ", gap_en, " are ", more_en, " without a difference; it rests mainly on its prior (", .bay_pct(h$prior[i], "."), "), lowered by the data.")
				fr[i] <- paste0("Contre elle : ", gap_fr, " sont ", more_fr, " sans diff\u00e9rence ; elle tient surtout \u00e0 son a priori (", .bay_pct(h$prior[i], ","), "), abaiss\u00e9 par les donn\u00e9es.")
			}
		} else {
			if (for_it) {
				en[i] <- paste0("Rests on ", if (i == 1L) small_en else "", gap_en, ": the data are ", more_en, " without a difference than with one of the size expected (r = ", r_en, ").")
				fr[i] <- paste0("Repose sur ", if (i == 1L) small_fr else "", gap_fr, " : les donn\u00e9es sont ", more_fr, " sans diff\u00e9rence qu'avec une diff\u00e9rence de l'ordre attendu (r = ", r_fr, ").")
			} else {
				en[i] <- paste0("Against it: ", gap_en, " are ", more_en, " with a difference.")
				fr[i] <- paste0("Contre elle : ", gap_fr, " sont ", more_fr, " avec une diff\u00e9rence.")
			}
		}
	}
	list(en = en, fr = fr)
}

# Which groups differ, in the conclusion: one line per pair, qualifier of the leading side.
.bay_print_posthoc <- function(ph, say) {
	if (identical(ph$method, "slopes")) {
		say(paste0("Which groups differ along ", ph$xname, " (slopes that differ; gap and 95 % credible interval at the quartiles, zone where they differ):"),
			paste0("Quels groupes diff\u00e8rent selon ", ph$xname, " (pentes diff\u00e9rentes ; \u00e9cart et intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 % aux quartiles, zone o\u00f9 ils diff\u00e8rent) :"))
		for (L in .bay_slopes_lines(ph, ph$xname)) say(paste0("  ", L$en), paste0("  ", L$fr))
		return(invisible(NULL))
	}
	pr <- ph$pairs
	f <- function(v, dc) gsub(".", dc, format(signif(v, 2), trim = TRUE), fixed = TRUE)
	ci <- function(i, dc, a, l, h) paste0(f(pr[[a]][i], dc), " [", f(pr[[l]][i], dc), if (dc == ",") " ; " else "; ", f(pr[[h]][i], dc), "]")
	if (ph$k > 2L) {
		if (isTRUE(ph$delta_unit %in% c("logor", "tvd"))) {
			say(paste0("Which groups differ in their proportions (probability, then estimated gap", if (ph$delta_unit == "logor") paste0(" of the proportion of ", ph$category, ", log odds ratio") else " over the categories, 0 to 1", ", 95 % credible intervals):"),
				paste0("Quels groupes diff\u00e8rent par leurs proportions (probabilit\u00e9, puis \u00e9cart estim\u00e9", if (ph$delta_unit == "logor") paste0(" de la proportion de ", ph$category, ", log du rapport de cotes") else " sur les modalit\u00e9s, de 0 \u00e0 1", ", intervalles de cr\u00e9dibilit\u00e9 \u00e0 95 %) :"))
		} else if (identical(ph$delta_unit, "adjusted")) {
			say("Which groups differ at the mean of the covariates (normal law; probability, then adjusted gap and 95 % credible interval):",
				"Quels groupes diff\u00e8rent \u00e0 la moyenne des covariables (loi normale ; probabilit\u00e9, puis \u00e9cart ajust\u00e9 et intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %) :")
		} else if (identical(ph$delta_unit, "within")) {
			say("Which conditions differ (normal law, random subject effect; probability, then estimated gap and 95 % credible interval):",
				"Quelles conditions diff\u00e8rent (loi normale, effet sujet al\u00e9atoire ; probabilit\u00e9, puis \u00e9cart estim\u00e9 et intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %) :")
		} else say("Which groups differ (normal law with a common variance; probability, then estimated gap and 95 % credible interval):",
			"Quels groupes diff\u00e8rent (loi normale \u00e0 variance commune ; probabilit\u00e9, puis \u00e9cart estim\u00e9 et intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %) :")
		shown <- .bay_posthoc_shown(pr)
		for (i in shown) {
			p <- pr$p_differ[i]
			open <- isTRUE(pr$p_min[i] < 0.5 && pr$p_max[i] > 0.5)
			s <- .bay_support(max(p, 1 - p))
			side_en <- if (open || s$level == 1L) "not settled" else if (p >= 0.5) paste0("different, ", s$en) else paste0("alike, ", s$en)
			side_fr <- if (open || s$level == 1L) "non tranch\u00e9" else if (p >= 0.5) paste0("diff\u00e9rents, ", s$fr) else paste0("semblables, ", s$fr)
			say(paste0("  ", pr$group1[i], " / ", pr$group2[i], ": ", side_en, " (different at ", .bay_pct(p, "."), "); estimated gap ",
					pr$group2[i], " - ", pr$group1[i], " = ", ci(i, ".", "diff", "diff_low", "diff_high"), if (is.na(pr$delta[i])) "" else paste0(", delta = ", ci(i, ".", "delta", "delta_low", "delta_high")), "."),
				paste0("  ", pr$group1[i], " / ", pr$group2[i], " : ", side_fr, " (diff\u00e9rents \u00e0 ", .bay_pct(p, ","), ") ; \u00e9cart estim\u00e9 ",
					pr$group2[i], " - ", pr$group1[i], " = ", ci(i, ",", "diff", "diff_low", "diff_high"), if (is.na(pr$delta[i])) "" else paste0(", delta = ", ci(i, ",", "delta", "delta_low", "delta_high")), "."))
		}
		if (length(shown) < nrow(pr)) say(paste0("  The other ", nrow(pr) - length(shown), " pairs lean towards alike; all pairs: $posthoc$pairs."),
			paste0("  Les ", nrow(pr) - length(shown), " autres paires penchent vers semblables ; toutes les paires : $posthoc$pairs."))
		if (identical(ph$method, "gibbs")) say(paste0("  Groupings explored by draws (more than ", ph$max_exact, " groups): about ", max(1L, ceiling(200 * ph$mc_error)), " point(s) of sampling error."),
			paste0("  Regroupements explor\u00e9s par tirage (plus de ", ph$max_exact, " groupes) : environ ", max(1L, ceiling(200 * ph$mc_error)), " point(s) d'erreur de tirage."))
		say(paste0("  Most probable grouping: ", ph$partitions$partition[1L], " (", .bay_pct(ph$partitions$posterior[1L], "."), ")."),
			paste0("  Regroupement le plus probable : ", ph$partitions$partition[1L], " (", .bay_pct(ph$partitions$posterior[1L], ","), ")."))
	} else {
		unit <- .bay_delta_unit(ph)
		if (isTRUE(ph$delta_unit %in% c("logor", "tvd"))) say(paste0("Size of the difference (estimated gap", if (ph$delta_unit == "logor") paste0(" of the proportion of ", ph$category, ", then log odds ratio") else " over the categories, 0 to 1", "; 95 % credible intervals):"),
			paste0("Taille de l'\u00e9cart (\u00e9cart estim\u00e9", if (ph$delta_unit == "logor") paste0(" de la proportion de ", ph$category, ", puis log du rapport de cotes") else " sur les modalit\u00e9s, de 0 \u00e0 1", " ; intervalles de cr\u00e9dibilit\u00e9 \u00e0 95 %) :")) else
		say(paste0("Size of the difference, if it exists (estimated gap: posterior median, and 95 % credible interval; delta in ", unit$en, "):"),
			paste0("Taille de l'\u00e9cart, s'il existe (\u00e9cart estim\u00e9 : m\u00e9diane a posteriori, et intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 % ; delta en ", unit$fr, ") :"))
		for (i in seq_len(nrow(pr))) {
			say(paste0("  ", pr$group2[i], " - ", pr$group1[i], " = ", ci(i, ".", "diff", "diff_low", "diff_high"), if (is.na(pr$delta[i])) "" else paste0(", delta = ", ci(i, ".", "delta", "delta_low", "delta_high")), "."),
				paste0("  ", pr$group2[i], " - ", pr$group1[i], " = ", ci(i, ",", "diff", "diff_low", "diff_high"), if (is.na(pr$delta[i])) "" else paste0(", delta = ", ci(i, ",", "delta", "delta_low", "delta_high")), "."))
		}
	}
	if (isTRUE(ph$caution_overlap)) say("  To read with caution: these gaps compare the groups where some have few or no values of the covariate.",
		"  \u00c0 lire avec r\u00e9serve : ces \u00e9carts comparent les groupes l\u00e0 o\u00f9 certains ont peu ou pas de valeurs de la covariable.")
	if (isTRUE(ph$caution)) say("  To read with caution: these pairs assume a normal noise with a common variance, which the comparison of noise forms did not find.",
		"  \u00c0 lire avec r\u00e9serve : ces paires supposent un bruit normal \u00e0 variance commune, ce que la comparaison des formes du bruit n'a pas trouv\u00e9.")
	invisible(NULL)
}

# Precision of a computation, one rule: below 0.1, "< 0.1".
.bay_prec_txt <- function(x, unit = "%") {
	if (!is.finite(x)) return(list(en = "unknown", fr = "inconnue"))
	if (x < 0.1) return(list(en = paste0("< 0.1 ", unit), fr = paste0("< 0,1 ", unit)))
	v <- format(round(x, 1), nsmall = 1, trim = TRUE)
	list(en = paste0("about ", v, " ", unit), fr = paste0("environ ", gsub(".", ",", v, fixed = TRUE), " ", unit))
}
