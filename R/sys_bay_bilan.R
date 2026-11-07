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
.bay_support_breaks <- function() c(0.75, 0.95, 0.993)

.bay_support <- function(p) {
	breaks <- .bay_support_breaks()
	idx <- findInterval(p, breaks) + 1L
	en <- c("weak support", "moderate support", "strong support", "very strong support")
	fr <- c("soutien faible", "soutien mod\u00e9r\u00e9", "soutien fort", "soutien tr\u00e8s fort")
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

.bay_law_labels <- function() {
	list(
		en = c(normal_eq = "normal, equal variances", normal_het = "normal, unequal variances",
			t_eq = "heavy tails, equal variances", t_het = "heavy tails, unequal variances"),
		fr = c(normal_eq = "normale, variances \u00e9gales", normal_het = "normale, variances in\u00e9gales",
			t_eq = "queues lourdes, variances \u00e9gales", t_het = "queues lourdes, variances in\u00e9gales")
	)
}

# Prints the bilan. obj is the list returned by m.bay(). Text goes through .msg: one
# string in each language, numbers with the decimal mark of that language.
.bay_print_bilan <- function(obj) {
	lang_fr <- identical(.msg("en", "fr"), "fr")
	dec <- if (lang_fr) "," else "."
	say <- function(en, fr) cat(.msg(en, fr), "\n", sep = "")
	num <- function(x) gsub(".", dec, format(signif(x, 3), trim = TRUE), fixed = TRUE)
	say(paste0("Bilan of m.bay: ", obj$label, " (N = ", obj$n, ")"),
		paste0("Bilan de m.bay : ", obj$label, " (N = ", obj$n, ")"))

	# Hypotheses.
	say("Hypotheses, from the most to the least probable:",
		"Hypoth\u00e8ses, de la plus probable \u00e0 la moins probable :")
	h <- obj$hypotheses
	lead <- h$posterior[1L]
	sup <- .bay_support(lead)
	width_en <- max(nchar(h$label_en))
	width_fr <- max(nchar(h$label_fr))
	for (i in seq_len(nrow(h))) {
		tail_en <- if (i == 1L) paste0(" (", sup$en, ")") else ""
		tail_fr <- if (i == 1L) paste0(" (", sup$fr, ")") else ""
		say(paste0("  ", i, ". ", formatC(h$label_en[i], width = -width_en), " : ",
				.bay_pct(h$posterior[i], "."), tail_en),
			paste0("  ", i, ". ", formatC(h$label_fr[i], width = -width_fr), " : ",
				.bay_pct(h$posterior[i], ","), tail_fr))
	}

	# Verdict of the sensitivity, next to what it qualifies.
	sens <- obj$sensitivity
	if (!is.null(sens)) {
		rng <- sens$range_lead
		rng_en <- paste0(.bay_pct(rng[1L], "."), " to ", .bay_pct(rng[2L], "."))
		rng_fr <- paste0(.bay_pct(rng[1L], ","), " \u00e0 ", .bay_pct(rng[2L], ","))
		if (identical(sens$verdict, "changes")) {
			say(paste0("  The data do not settle it: over the priors, this probability goes from ",
					rng_en, " and the leading hypothesis changes."),
				paste0("  Les donn\u00e9es ne tranchent pas : selon les a priori, cette probabilit\u00e9 va de ",
					rng_fr, " et l'hypoth\u00e8se en t\u00eate change."))
		} else if (identical(sens$verdict, "crosses")) {
			a <- .bay_support(rng[1L])
			b <- .bay_support(rng[2L])
			say(paste0("  The leading hypothesis does not change over the priors, but the support does: ",
					rng_en, " (", a$en, " to ", b$en, ")."),
				paste0("  L'hypoth\u00e8se en t\u00eate ne change pas selon les a priori, mais le soutien si : ",
					rng_fr, " (", a$fr, " \u00e0 ", b$fr, ")."))
		} else {
			say(paste0("  Stable over the priors: ", rng_en, "."),
				paste0("  Stable selon les a priori : ", rng_fr, "."))
		}
	} else {
		say("  Sensitivity to the priors not computed (sensitivity = FALSE).",
			"  Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}

	# Per term.
	say("Answers per term:", "R\u00e9ponses par terme :")
	e <- obj$effects
	for (i in seq_len(nrow(e))) {
		say(paste0("  Does ", e$term[i], " have an effect? Yes ", .bay_pct(e$p_h1[i], "."),
				", no ", .bay_pct(e$p_h0[i], "."),
				" (prior ", .bay_pct(e$prior_h1[i], "."), "; Bayes factor ",
				.bay_bf_txt(e$log_bf[i], "."), ")."),
			paste0("  ", e$term[i], " a-t-il un effet ? Oui ", .bay_pct(e$p_h1[i], ","),
				", non ", .bay_pct(e$p_h0[i], ","),
				" (a priori ", .bay_pct(e$prior_h1[i], ","), " ; facteur de Bayes ",
				.bay_bf_txt(e$log_bf[i], ","), ")."))
	}

	# Size, descriptive.
	d <- obj$descriptive_stats
	fmt <- function(x, dc) gsub(".", dc, format(signif(x, 4), trim = TRUE), fixed = TRUE)
	say("Observed means (description of the sample, not a test):",
		"Moyennes observ\u00e9es (description de l'\u00e9chantillon, pas un test) :")
	say(paste0("  ", paste0(d$group, " = ", vapply(d$mean, fmt, character(1), dc = "."), " (n = ", d$n, ")", collapse = "; "),
			". Largest gap: ", fmt(obj$largest_gap, "."), " within-group standard deviations."),
		paste0("  ", paste0(d$group, " = ", vapply(d$mean, fmt, character(1), dc = ","), " (n = ", d$n, ")", collapse = " ; "),
			". Plus grand \u00e9cart : ", fmt(obj$largest_gap, ","), " \u00e9cart-type intra-groupe."))

	# Noise.
	if (!is.null(obj$laws)) {
		lb <- .bay_law_labels()
		laws <- obj$laws[order(-obj$laws$posterior), , drop = FALSE]
		txt_en <- paste0(lb$en[laws$law], " ", vapply(laws$posterior, .bay_pct, character(1), dec = "."))
		txt_fr <- paste0(lb$fr[laws$law], " ", vapply(laws$posterior, .bay_pct, character(1), dec = ","))
		say("Shape of the noise, compared and not chosen (posterior probabilities):",
			"Forme du bruit, compar\u00e9e et non choisie (probabilit\u00e9s a posteriori) :")
		say(paste0("  ", paste(txt_en, collapse = "; "), "."),
			paste0("  ", paste(txt_fr, collapse = " ; "), "."))
		if (!is.null(sens) && nrow(obj$laws) > 1L) {
			bl <- sens$by_law[laws$law]
			say(paste0("  If one shape were imposed, the probability of an effect would be: ",
					paste0(lb$en[laws$law], " ", vapply(bl, .bay_pct, character(1), dec = "."), collapse = "; "), "."),
				paste0("  Si l'on imposait une seule forme, la probabilit\u00e9 d'un effet serait : ",
					paste0(lb$fr[laws$law], " ", vapply(bl, .bay_pct, character(1), dec = ","), collapse = " ; "), "."))
		}
	}

	# Priors.
	p <- obj$priors
	say("Priors, fixed in advance:", "A priori, fix\u00e9s d'avance :")
	say(paste0("  Effect: ", .bay_pct(p$prior_h1, "."), " before seeing the data. Size of an effect: ",
			"scale r = ", format(signif(p$rscale, 3), trim = TRUE), " on the standardised difference."),
		paste0("  Effet : ", .bay_pct(p$prior_h1, ","), " avant de voir les donn\u00e9es. Taille d'un effet : ",
			"\u00e9chelle r = ", num(p$rscale), " sur l'\u00e9cart standardis\u00e9."))
	if (isTRUE(nrow(obj$descriptive_stats) > 2L)) {
		say("  Forms of an effect, at even odds: all groups may differ; or one group alone stands apart (scale r x sqrt(2)).",
			"  Formes d'un effet, \u00e0 poids \u00e9gaux : tous les groupes peuvent diff\u00e9rer ; ou un seul groupe s'\u00e9carte (\u00e9chelle r x racine de 2).")
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
				format(max(sens$r), trim = TRUE), ", odds H1:H0 ", paste(sens$odds[c(1L, length(sens$odds))], collapse = " to "), ", Dirichlet ",
				format(min(sens$u), trim = TRUE), " to ", format(max(sens$u), trim = TRUE), "."),
			paste0("  Vari\u00e9s dans la sensibilit\u00e9 : r de ", num(min(sens$r)), " \u00e0 ",
				num(max(sens$r)), ", cotes H1:H0 ", paste(sens$odds[c(1L, length(sens$odds))], collapse = " \u00e0 "), ", Dirichlet ",
				num(min(sens$u)), " \u00e0 ", num(max(sens$u)), "."))
	}

	if (!is.null(sens) && isTRUE(nrow(obj$descriptive_stats) > 2L)) {
		say(paste0("  With the first form of an effect alone, the probability of an effect would be ",
				.bay_pct(sens$all_groups, "."), "."),
			paste0("  Avec la premi\u00e8re forme d'effet seule, la probabilit\u00e9 d'un effet serait ",
				.bay_pct(sens$all_groups, ","), "."))
	}

	# How to read it.
	breaks <- .bay_support_breaks()
	say("How to read it:", "Comment lire :")
	say("  A probability here is the probability of a hypothesis given the data and the priors above.",
		"  Une probabilit\u00e9 est ici celle d'une hypoth\u00e8se, sachant les donn\u00e9es et les a priori ci-dessus.")
	say("  The Bayes factor says how much the data move the odds. It does not depend on the prior odds.",
		"  Le facteur de Bayes dit de combien les donn\u00e9es d\u00e9placent les cotes. Il ne d\u00e9pend pas des cotes a priori.")
	say("  Neither is a p-value: a p-value is the probability of data at least this extreme if there were no effect.",
		"  Ni l'un ni l'autre n'est une p-value : une p-value est la probabilit\u00e9 de donn\u00e9es au moins aussi extr\u00eames s'il n'y avait pas d'effet.")
	say(paste0("  Support: weak under ", 100 * breaks[1L], " %, moderate under ", 100 * breaks[2L],
			" %, strong under ", 100 * breaks[3L], " %, then very strong."),
		paste0("  Soutien : faible sous ", 100 * breaks[1L], " %, mod\u00e9r\u00e9 sous ", 100 * breaks[2L],
			" %, fort sous ", gsub(".", ",", 100 * breaks[3L], fixed = TRUE), " %, puis tr\u00e8s fort."))
	if (length(obj$notes$en)) {
		say("Warnings:", "Avertissements :")
		for (i in seq_along(obj$notes$en)) say(paste0("  ", obj$notes$en[i]), paste0("  ", obj$notes$fr[i]))
	}
	invisible(obj)
}
