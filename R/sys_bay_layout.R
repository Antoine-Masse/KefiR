# sys_bay_layout.R
# The bilan of m.bay, one skeleton for every plan (relecture/lisibilite_2026-11-09): the
# last step of the story, visible at once ("7) Bilan - <leading hypothesis>: xx % (support)."),
# then the verdict, the hypotheses numbered and aligned, what the answer rests on, its
# stability, the answers of the plan (terms, pairs, slopes), what was observed, and short
# notes. The pairs are given in detail at step 6 only; here, their short list. Every line
# is folded within 80 columns by the emitter of sys_bay_story.R.

.bay_lang_fr <- function() identical(.msg("en", "fr"), "fr")

# Number of decimals of the response: that of the means (3 significant digits).
.bay_dg <- function(x) {
	m <- max(abs(x[is.finite(x)]), 1e-300)
	min(6L, max(0L, 2L - floor(log10(m))))
}
# Three significant digits (a standard deviation in the bilan, audit of a206ba5).
.bay_sig3 <- function(x, dc) gsub(".", dc, format(signif(x, 3), trim = TRUE, scientific = FALSE), fixed = TRUE)

# A value that rounds to zero is written without its sign ("-0.00" read as negative).
.bay_fx <- function(x, dg, dc) {
	x[!is.na(x) & round(x, dg) == 0] <- 0
	gsub(".", dc, formatC(x, format = "f", digits = dg), fixed = TRUE)
}

# Leading hypothesis: its label, probability, support and the Bayes factor in its favour.
.bay_lead_info <- function(obj) {
	h <- obj$hypotheses
	lb <- log(h$posterior[1L] / h$posterior[2L]) - log(h$prior[1L] / h$prior[2L])
	sup <- .bay_support_bf(h$posterior[1L], lb)
	if (isTRUE(obj$contradiction)) sup <- .bay_support(min(h$posterior[1L], 0.6))
	if (isTRUE(any(obj$overlap$poor)) && sup$level > 2L) sup <- .bay_support(0.7)
	if (isTRUE(obj$het_cells) && !isTRUE(h$null[1L]) && sup$level > 2L) sup <- .bay_support(0.7)
	open <- identical(obj$sensitivity$verdict, "changes")
	if (open) sup <- list(en = "not settled", fr = "non tranch\u00e9", level = 1L)
	list(h = h, lb = lb, sup = sup, open = open)
}

# How to cite the method, per plan.
.bay_cite <- function(obj) {
	switch(obj$plan,
		oneway = list(en = "comparison of noise laws, Rouder et al. 2012", fr = "comparaison de lois du bruit, Rouder et al. 2012"),
		paired = list(en = "Bayesian paired t test, Rouder et al. 2009", fr = "test t appari\u00e9 bay\u00e9sien, Rouder et al. 2009"),
		contingency = list(en = "Gunel and Dickey 1974", fr = "Gunel et Dickey 1974"),
		ordinal = list(en = "proportional odds, MASS::polr, BIC approximation", fr = "cotes proportionnelles, MASS::polr, approximation BIC"),
		counts = list(en = "Poisson and negative binomial, BIC approximation, Raftery 1995", fr = "Poisson et binomiale n\u00e9gative, approximation BIC, Raftery 1995"),
		discrete_fact = if (identical(obj$kind, "binary")) list(en = "logistic regression, BIC approximation, HOP prior, Raftery 1995", fr = "r\u00e9gression logistique, approximation BIC, a priori HOP, Raftery 1995") else
			list(en = "Poisson and negative binomial, BIC approximation, HOP prior, Raftery 1995", fr = "Poisson et binomiale n\u00e9gative, approximation BIC, a priori HOP, Raftery 1995"),
		list(en = "BayesFactor, Rouder et al. 2012", fr = "BayesFactor, Rouder et al. 2012"))
}

.bay_bilan_common <- function(obj, narrated = FALSE, body = NULL) {
	emit <- .bay_block_emitter()
	fr <- .bay_lang_fr(); dc <- if (fr) "," else "."
	li <- .bay_lead_info(obj)
	h <- li$h
	pct <- function(p) .bay_pct(p, dc)
	k <- if (!is.null(obj$story$k)) obj$story$k + 1L else 7L
	# Undecided: the title does not assert the leading hypothesis (relecture apres, point 1).
	undecided <- li$open || li$sup$level <= 1L || h$posterior[1L] < 0.6
	lead_en <- if (undecided) paste0("undecided (", h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), ")") else
		paste0(h$label_en[1L], ": ", .bay_pct(h$posterior[1L], "."), " (", li$sup$en, ")")
	lead_fr <- if (undecided) paste0("ind\u00e9cis (", h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), ")") else
		paste0(h$label_fr[1L], " : ", .bay_pct(h$posterior[1L], ","), " (", li$sup$fr, ")")
	# Several terms with at least one "yes" (audit of 9c44239, point 1): the title, the verdict
	# and the report sentence lead with the answers per term, the global "at least one
	# difference" comes second, and "undecided" is never said.
	e <- obj$effects
	term_lead <- FALSE
	if (!is.null(e) && nrow(e) > 1L && !is.null(e$p_h0)) {
		ta <- .bay_term_answers(obj)
		term_lead <- any(!ta$open & ta$yes)
		if (term_lead) {
			tl_en <- paste0(e$term, " ", ta$short_en, " (", vapply(e$p_h1, .bay_pct, character(1), dec = "."), ")", collapse = ", ")
			tl_fr <- paste0(e$term, " ", ta$short_fr, " (", vapply(e$p_h1, .bay_pct, character(1), dec = ","), ")", collapse = ", ")
			lead_en <- tl_en
			lead_fr <- tl_fr
			undecided <- FALSE
			# The global hypothesis keeps the support of its own Bayes factor, on the scale of
			# the bilan, never "not settled" under a "yes" (audit of aed2319).
			bfa <- exp(abs(li$lb))
			lv <- findInterval(bfa, c(3, 20, 150)) + 1L
			li$sup <- list(en = c("weak support", "moderate support", "strong support", "very strong support")[lv],
				fr = c("soutien faible", "soutien mod\u00e9r\u00e9", "soutien fort", "soutien tr\u00e8s fort")[lv], level = lv + 1L)
		}
	}
	if (isTRUE(narrated)) {
		cat("\n")
		emit("title", paste0(k, ") Bilan - ", lead_en, "."), paste0(k, ") Bilan - ", lead_fr, "."))
	} else {
		ns <- obj$n_subjects
		sub_en <- if (!is.null(ns) && !is.na(ns)) paste0(", ", ns, " subjects") else ""
		sub_fr <- if (!is.null(ns) && !is.na(ns)) paste0(", ", ns, " sujets") else ""
		emit("title", paste0("Bilan of m.bay (", obj$label, ", N = ", obj$n, sub_en, ") - ", lead_en, "."),
			paste0("Bilan de m.bay (", obj$label, ", N = ", obj$n, sub_fr, ") - ", lead_fr, "."))
	}
	v <- .bay_verdict(obj)
	first <- function(x) sub("^([^.]*\\.).*$", "\\1", x)
	if (term_lead) {
		yes_t <- e$term[!ta$open & ta$yes]
		emit("fact", paste0("Verdict: a difference according to ", paste(yes_t, collapse = ", "), "; the other terms are below."),
			paste0("Verdict : une diff\u00e9rence selon ", paste(yes_t, collapse = ", "), " ; les autres termes sont plus bas."))
	} else {
		emit("fact", first(v$en), first(v$fr))
	}
	w_en <- max(nchar(h$label_en)); w_fr <- max(nchar(h$label_fr))
	# Padding in non-breaking spaces, which the folding keeps.
	# Aligned only if the line still fits; otherwise, no padding.
	fit <- max(w_en, w_fr) + 32L <= .bay_width() - 12L
	pad <- function(x, w) if (fit) paste0(x, strrep("\u00a0", w - nchar(x))) else x
	emit("fact", "Probability of each hypothesis:", "Probabilit\u00e9 de chaque hypoth\u00e8se :")
	for (i in seq_len(nrow(h))) {
		emit("sub", paste0(i, ". ", pad(h$label_en[i], w_en), " : ", .bay_pct(h$posterior[i], "."), if (i == 1L) paste0(" (", li$sup$en, ")") else ""),
			paste0(i, ". ", pad(h$label_fr[i], w_fr), " : ", .bay_pct(h$posterior[i], ","), if (i == 1L) paste0(" (", li$sup$fr, ")") else ""))
	}
	bf_en <- .bay_bf_txt(abs(li$lb), "."); bf_fr <- .bay_bf_txt(abs(li$lb), ",")
	emit("note", paste0("the data are ", bf_en, " times more probable under hypothesis 1 than under 2 (Bayes factor; 1 to 3 weak, 3 to 20 moderate, 20 to 150 strong, beyond very strong)."),
		paste0("les donn\u00e9es sont ", bf_fr, " fois plus probables sous l'hypoth\u00e8se 1 que sous la 2 (facteur de Bayes ; 1 \u00e0 3 faible, 3 \u00e0 20 mod\u00e9r\u00e9, 20 \u00e0 150 fort, au-del\u00e0 tr\u00e8s fort)."))
	sens <- obj$sensitivity
	if (!is.null(sens$range_lead)) {
		rg <- sens$range_lead
		if (li$open) {
			emit("decide", paste0("Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order can reverse, not settled."),
				paste0("Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre peut s'inverser, non tranch\u00e9."))
		} else if (identical(.bay_pct(rg[1L], "."), .bay_pct(rg[2L], ".")) || rg[1L] > 0.999) {
			emit("fact", paste0("The answer ", .bay_range_txt(rg)$en, "."), paste0("La r\u00e9ponse ", .bay_range_txt(rg)$fr, "."))
		} else {
			emit("fact", paste0("Depending on the starting point: ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "; the order does not change."),
				paste0("Selon le point de d\u00e9part : ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), " ; l'ordre ne change pas."))
		}
	}
	if (is.null(sens) && !isTRUE(narrated)) {
		emit("note", "sensitivity to the priors not computed (sensitivity = FALSE).",
			"sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
	}
	emit("block", "", "")
	if (is.function(body)) body(emit, li, narrated)
	ct <- .bay_cite(obj)
	h1 <- which(!h$null)[1L]
	terms_en <- terms_fr <- ""
	e <- obj$effects
	if (!is.null(e) && nrow(e) > 1L) {
		# Same words as the lines of the bilan (.bay_term_answers).
		ta <- .bay_term_answers(obj)
		st_en <- ta$short_en
		st_fr <- ta$short_fr
		# A main effect inside a probable interaction: its "yes" may come from the interaction
		# alone (relecture D19-D20, point 2).
		fs <- strsplit(.bay_term_key(e$term), ":", fixed = TRUE)
		in_int <- vapply(seq_along(fs), function(i) length(fs[[i]]) == 1L && any(vapply(seq_along(fs), function(j) {
			length(fs[[j]]) > 1L && fs[[i]] %in% fs[[j]] && isTRUE(e$p_h1[j] >= 0.5)
		}, logical(1))), logical(1))
		yes_in <- in_int & !ta$open & ta$yes
		st_en[yes_in] <- "yes, alone or in interaction,"
		st_fr[yes_in] <- "oui, seul ou en interaction,"
		terms_en <- paste0("; ", paste0(e$term, " ", st_en, " (", vapply(e$p_h1, .bay_pct, character(1), dec = "."), ")", collapse = ", "))
		terms_fr <- paste0(" ; ", paste0(e$term, " ", st_fr, " (", vapply(e$p_h1, .bay_pct, character(1), dec = ","), ")", collapse = ", "))
	}
	if (term_lead) {
		# Terms first, the global answer second (audit of 9c44239, point 1).
		lb_any <- if (h$null[1L]) -li$lb else li$lb
		emit("note", paste0("for a report: \"", sub("^; ", "", terms_en), "; at least one difference: ", .bay_pct(h$posterior[h1], "."),
				", Bayes factor ", .bay_bf_txt(lb_any, "."), " (m.bay {KefiR}, ", ct$en, ").\""),
			paste0("pour un compte rendu : \u00ab ", sub("^ ; ", "", terms_fr), " ; au moins une diff\u00e9rence : ", .bay_pct(h$posterior[h1], ","),
				", facteur de Bayes ", .bay_bf_txt(lb_any, ","), " (m.bay {KefiR}, ", ct$fr, "). \u00bb"))
	} else if (undecided) {
		emit("note", paste0("for a report: \"The data do not settle it (probability of a difference ", .bay_pct(h$posterior[h1], "."), ", Bayes factor ", bf_en,
				if (h$null[1L]) " in favour of no difference" else "", terms_en, "); m.bay {KefiR}, ", ct$en, ".\""),
			paste0("pour un compte rendu : \u00ab Les donn\u00e9es ne permettent pas de conclure (probabilit\u00e9 d'une diff\u00e9rence ", .bay_pct(h$posterior[h1], ","), ", facteur de Bayes ", bf_fr,
				if (h$null[1L]) " en faveur de l'absence de diff\u00e9rence" else "", terms_fr, ") ; m.bay {KefiR}, ", ct$fr, ". \u00bb"))
	} else {
		emit("note", paste0("for a report: \"", h$label_en[1L], " (posterior probability ", .bay_pct(h$posterior[1L], "."), "), Bayes factor ", bf_en, terms_en,
				" (m.bay {KefiR}, ", ct$en, ").\""),
			paste0("pour un compte rendu : \u00ab ", h$label_fr[1L], " (probabilit\u00e9 a posteriori ", .bay_pct(h$posterior[1L], ","), "), facteur de Bayes ", bf_fr, terms_fr,
				" (m.bay {KefiR}, ", ct$fr, "). \u00bb"))
	}
	if (length(obj$notes$en) && !isTRUE(narrated)) {
		for (i in seq_along(obj$notes$en)) emit("fact", obj$notes$en[i], obj$notes$fr[i])
	}
	cause_en <- if (!h$null[1L] && !undecided) "; neither the size nor the cause of a difference" else ""
	cause_fr <- if (!h$null[1L] && !undecided) " ; ni l'importance, ni la cause d'une diff\u00e9rence" else ""
	emit("note", paste0("probabilities of hypotheses, not p-values", cause_en, "; the method with print(result, detail = TRUE)."),
		paste0("probabilit\u00e9s d'hypoth\u00e8ses, pas des p-values", cause_fr, " ; la m\u00e9thode avec print(resultat, detail = TRUE)."))
	emit("flush", "", "")
	invisible(emit)
}

# The details kept by the story (method, references), then those of the printer.
.bay_detail_head <- function(obj) {
	emit <- .bay_emitter()
	cat("\n")
	emit("title", "Method in detail:", "M\u00e9thode en d\u00e9tail :")
	d_en <- obj$story$detail_en; d_fr <- obj$story$detail_fr
	for (i in seq_along(d_en)) emit("note", d_en[i], d_fr[i])
	function(en, fr) emit("note", sub("^ +", "", en), sub("^ +", "", fr))
}

# Short list of the pairs, grouped by answer and support; "alike" only when the bilan
# leans towards a difference somewhere, otherwise "not settled" (relecture D17, point 3).
.bay_pairs_short <- function(ph, emit, undecided, narrated, note = TRUE) {
	pr <- ph$pairs
	if (is.null(pr) || anyNA(pr$p_differ)) return(invisible(NULL))
	lab <- paste0(pr$group1, "/", pr$group2)
	cls <- character(nrow(pr)); s_en <- s_fr <- character(nrow(pr))
	for (i in seq_len(nrow(pr))) {
		p <- pr$p_differ[i]
		s <- .bay_support(max(p, 1 - p))
		open <- isTRUE(pr$p_min[i] < 0.5 && pr$p_max[i] > 0.5) || s$level == 1L
		# One fixed rule (relecture apres, point 6): different from 75 %, alike up to 25 %.
		cls[i] <- if (isTRUE(pr$p_min[i] < 0.5 && pr$p_max[i] > 0.5)) "open" else if (p >= 0.75) "diff" else if (p <= 0.25) "same" else "open"
		s_en[i] <- sub(" support", "", s$en); s_fr[i] <- sub("soutien ", "", s$fr)
	}
	grp <- function(sel, by_support) {
		if (!any(sel)) return(NULL)
		key <- if (by_support) s_en[sel] else rep("", sum(sel))
		pe <- pf <- character(0)
		for (kk in unique(key)) {
			lb <- paste(lab[sel][key == kk], collapse = ", ")
			kf <- unique(s_fr[sel][key == kk])[1L]
			pe <- c(pe, paste0(lb, if (nzchar(kk)) paste0(" (", kk, ")") else ""))
			pf <- c(pf, paste0(lb, if (nzchar(kk)) paste0(" (", kf, ")") else ""))
		}
		list(en = paste(pe, collapse = "; "), fr = paste(pf, collapse = " ; "))
	}
	a <- grp(cls == "diff", TRUE); b <- grp(cls == "same", TRUE); c0 <- grp(cls == "open", FALSE)
	if (!is.null(a)) emit("fact", paste0("Different: ", a$en, "."), paste0("Diff\u00e9rents : ", a$fr, "."))
	if (!is.null(b)) emit("fact", paste0("Alike: ", b$en, "."), paste0("Semblables : ", b$fr, "."))
	if (!is.null(c0)) emit("fact", paste0("Not settled: ", c0$en, "."), paste0("Non tranch\u00e9s : ", c0$fr, "."))
	top <- .bay_top_grouping(ph)
	if (!is.null(top)) {
		emit("fact", paste0("If a difference exists, most probable grouping: ", top$label, " (", .bay_pct(top$p, "."), ")."),
			paste0("Si une diff\u00e9rence existe, regroupement le plus probable : ", top$label, " (", .bay_pct(top$p, ","), ")."))
	}
	if (note) emit("note", paste0("different from 75 % of probability, alike up to 25 %, otherwise not settled; gaps in figures ", if (narrated) "at step 6." else "in $posthoc$pairs."),
		paste0("diff\u00e9rents \u00e0 partir de 75 % de probabilit\u00e9, semblables jusqu'\u00e0 25 %, sinon non tranch\u00e9s ; \u00e9carts chiffr\u00e9s ", if (narrated) "\u00e0 l'\u00e9tape 6." else "dans $posthoc$pairs."))
	invisible(NULL)
}

# Leading grouping among those with a difference, its probability given a difference.
.bay_top_grouping <- function(ph) {
	pt <- ph$partitions
	if (is.null(pt) || !nrow(pt)) return(NULL)
	null <- !grepl("[<>]", pt$partition)
	if (all(null)) return(NULL)
	nn <- pt[!null, , drop = FALSE]
	list(label = nn$partition[1L], p = nn$posterior[1L] / sum(nn$posterior))
}
