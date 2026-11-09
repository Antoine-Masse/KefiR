# sys_bay_story.R
# The story told by m.bay(verbose = TRUE), in the manner of m.test(): numbered steps
# printed as the computation goes, "==>" for what was found, "-->" for what the function
# decides from it, the function and its package in brackets. Every line is also kept in
# the returned object (story), so print(result, story = TRUE) replays it. Lines are
# printed through .vbse(); the language is chosen first by .msg(), because .vbse() keeps
# only its French argument when given one string (it serves m.test() as it is).

.bay_story <- function(verbose) {
	st <- new.env(parent = emptyenv())
	st$k <- 0L
	st$en <- character(0)
	st$fr <- character(0)
	st$verbose <- isTRUE(verbose)
	st$pending <- list()
	st
}

# kind: "step" (numbered title), "fact" (==>), "decide" (-->), "info" (indented).
.bay_tell <- function(st, en, fr, kind = "fact", fun = NULL) {
	if (is.null(st)) return(invisible(NULL))
	if (!is.null(fun)) {
		en <- paste0(en, " [", fun[1L], "]")
		fr <- paste0(fr, " [", fun[length(fun)], "]")
	}
	if (identical(kind, "step")) {
		st$k <- st$k + 1L
		lead <- paste0(st$k, ") ")
	} else {
		lead <- c(fact = "\t==> ", decide = "\t--> ", info = "\t")[[kind]]
	}
	st$en <- c(st$en, paste0(lead, en))
	st$fr <- c(st$fr, paste0(lead, fr))
	if (st$verbose) {
		txt <- .msg(en, fr)
		if (identical(kind, "step")) {
			.vbse(txt, txt, k = st$k - 1L, cpt = "on", verbose = TRUE)
		} else {
			body <- paste0(sub("^\t", "", lead), txt)
			.vbse(body, body, cpt = "off", verbose = TRUE)
		}
	}
	if (identical(kind, "step") && length(st$pending)) {
		p <- st$pending
		st$pending <- list()
		for (x in p) .bay_tell(st, x$en, x$fr, x$kind)
	}
	invisible(NULL)
}

# A note said before the first step waits for it.
.bay_tell_note <- function(st, en, fr, kind = "fact") {
	if (is.null(st)) return(invisible(NULL))
	if (st$k == 0L) {
		st$pending[[length(st$pending) + 1L]] <- list(en = en, fr = fr, kind = kind)
	} else {
		.bay_tell(st, en, fr, kind)
	}
	invisible(NULL)
}

.bay_story_lines <- function(st) list(en = st$en, fr = st$fr)

# Replays a stored story.
.bay_story_replay <- function(story) {
	if (is.null(story) || !length(story$en)) return(invisible(NULL))
	for (i in seq_along(story$en)) cat(.msg(story$en[i], story$fr[i]), "\n", sep = "")
	invisible(NULL)
}

# Numbers in the language of the session.
.bay_num_txt <- function(x, digits = 3) {
	en <- format(signif(x, digits), trim = TRUE)
	list(en = en, fr = gsub(".", ",", en, fixed = TRUE))
}

# Values far from their group (Tukey's rule, per group of at least 4 values):
# beyond 1.5 and 3 interquartile ranges from the quartiles. Descriptive only.
.bay_extremes <- function(y, g) {
	g <- droplevels(as.factor(g))
	out <- c(mild = 0L, far = 0L)
	for (v in split(y, g)) {
		if (length(v) < 4L) next
		q <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
		iqr <- q[2L] - q[1L]
		if (!isTRUE(iqr > 0)) next
		out["mild"] <- out["mild"] + sum(v < q[1L] - 1.5 * iqr | v > q[2L] + 1.5 * iqr)
		out["far"] <- out["far"] + sum(v < q[1L] - 3 * iqr | v > q[2L] + 3 * iqr)
	}
	out
}

# Step 2, common to both paths: sizes, ties, constant groups, values far from their group.
# Returns the counts so that the caller can decide on them.
.bay_story_checks <- function(st, y, g, unit_en = "group", unit_fr = "groupe") {
	g <- droplevels(as.factor(g))
	n <- as.integer(table(g))
	N <- length(y)
	.bay_tell(st, paste0("Sizes: ", min(n), if (max(n) > min(n)) paste0(" to ", max(n)), " per ", unit_en,
			" (", nlevels(g), " ", unit_en, "s, N = ", N, ")", if (max(n) > min(n)) ", unbalanced." else "."),
		paste0("Effectifs : ", min(n), if (max(n) > min(n)) paste0(" \u00e0 ", max(n)), " par ", unit_fr,
			" (", nlevels(g), " ", unit_fr, "s, N = ", N, ")", if (max(n) > min(n)) ", d\u00e9s\u00e9quilibr\u00e9." else "."))
	tied <- 1 - length(unique(y)) / N
	if (tied > 0) {
		pc <- .bay_num_txt(100 * tied, 2)
		.bay_tell(st, paste0("Ties: ", pc$en, " % of the values repeat an earlier one."),
			paste0("Ex-aequo : ", pc$fr, " % des valeurs r\u00e9p\u00e8tent une valeur d\u00e9j\u00e0 vue."))
	} else {
		.bay_tell(st, "Ties: none.", "Ex-aequo : aucun.")
	}
	const <- levels(g)[vapply(split(y, g), function(v) length(v) > 1L && isTRUE(diff(range(v)) == 0), logical(1))]
	if (length(const)) {
		.bay_tell(st, paste0("Constant ", unit_en, "(s): ", paste(const, collapse = ", "), "."),
			paste0(if (unit_fr == "case") "Case(s) constante(s) : " else "Groupe(s) constant(s) : ", paste(const, collapse = ", "), "."))
	}
	if (length(const)) {
		.bay_tell(st, paste0("To check: a value repeated identically through a whole ", unit_en, " often signals a measurement ceiling or an entry error; the analysis goes on without removing it."),
			paste0("\u00c0 v\u00e9rifier : une valeur r\u00e9p\u00e9t\u00e9e \u00e0 l'identique dans tout un", if (unit_fr == "case") "e " else " ", unit_fr, " signale souvent un plafond de mesure ou une erreur de saisie ; l'analyse continue sans le retirer."), "decide")
	}
	ex <- .bay_extremes(y, g)
	.bay_tell(st, paste0("Values far from their ", unit_en, ": ", ex[["mild"]], " of ", N, " beyond 1.5 interquartile ranges, ",
			ex[["far"]], " beyond 3."),
		paste0("Valeurs \u00e9loign\u00e9es de leur ", unit_fr, " : ", ex[["mild"]], " sur ", N, " au-del\u00e0 de 1,5 \u00e9cart interquartile, ",
			ex[["far"]], " au-del\u00e0 de 3."), fun = "Tukey, quantile() {stats}")
	list(n = n, tied = tied, constant = const, extremes = ex)
}

# Step 1: the plan recognised. A plan that is not handled is named before the refusal.
.bay_story_plan <- function(st, input) {
	if (is.null(st)) return(invisible(NULL))
	N <- length(input$response)
	fac <- input$factors
	nm <- names(fac)
	nm[nm == ".g"] <- "g"
	if (identical(input$plan, "oneway")) {
		f <- fac[[1L]]
		.bay_tell(st, paste0("Plan recognised: one factor, ", nm[1L], ", with ", nlevels(f), " groups (",
				paste(levels(f), collapse = ", "), "); N = ", N, "."),
			paste0("Plan reconnu : un facteur, ", nm[1L], ", \u00e0 ", nlevels(f), " groupes (",
				paste(levels(f), collapse = ", "), ") ; N = ", N, "."), "step")
		.bay_tell(st, paste0("One question: is there a difference in ", input$response_name, " according to ", nm[1L],
				"? No preliminary test chooses the method: the shapes of the noise will be compared, not chosen."),
			paste0("Une seule question : y a-t-il une diff\u00e9rence de ", input$response_name, " selon ", nm[1L],
				" ? Aucun test pr\u00e9alable ne choisit la m\u00e9thode : les formes du bruit seront compar\u00e9es, pas choisies."), "decide")
	} else if (identical(input$plan, "factorial")) {
		lv <- vapply(fac, nlevels, integer(1))
		cells <- prod(lv)
		terms <- gsub("`", "", input$terms, fixed = TRUE)
		.bay_tell(st, paste0("Plan recognised: ", length(fac), " crossed factors, ",
				paste0(nm, " (", lv, " levels)", collapse = ", "), "; ", cells, " cells, N = ", N, "."),
			paste0("Plan reconnu : ", length(fac), " facteurs crois\u00e9s, ",
				paste0(nm, " (", lv, " niveaux)", collapse = ", "), " ; ", cells, " cases, N = ", N, "."), "step")
		sp <- .bay_fact_space(terms)
		.bay_tell(st, paste0("Terms: ", paste(terms, collapse = ", "), ". ", length(sp$models),
				" combinations of terms are possible (an interaction enters only with its factors: heredity)."),
			paste0("Termes : ", paste(terms, collapse = ", "), ". ", length(sp$models),
				" combinaisons de termes sont possibles (une interaction n'entre qu'avec ses facteurs : h\u00e9r\u00e9dit\u00e9)."))
		.bay_tell(st, paste0("One answer per term, obtained by weighing these ", length(sp$models), " combinations together."),
			paste0("Une r\u00e9ponse par terme, obtenue en pesant ensemble ces ", length(sp$models), " combinaisons."), "decide")
	} else if (identical(input$plan, "ancova")) {
		lv <- vapply(fac, nlevels, integer(1))
		cv <- names(input$covariates)
		.bay_tell(st, paste0("Plan recognised: analysis of covariance, ", paste0(nm, " (", lv, " groups)", collapse = ", "),
				", covariate(s) ", paste(cv, collapse = ", "), "; N = ", N, "."),
			paste0("Plan reconnu : analyse de covariance, ", paste0(nm, " (", lv, " groupes)", collapse = ", "),
				", covariable(s) ", paste(cv, collapse = ", "), " ; N = ", N, "."), "step")
		terms <- .bay_space_terms(gsub("`", "", input$terms, fixed = TRUE), cv)
		invisible(.bay_fact_space(terms, cv))
		slopes <- terms[grepl(":", terms, fixed = TRUE)]
		.bay_tell(st, paste0("Question: does ", input$response_name, " differ according to ", paste(nm, collapse = ", "), " once ",
				paste(cv, collapse = ", "), " is accounted for?", if (length(slopes)) paste0(" And do the slopes differ (", paste(slopes, collapse = ", "), ")?") else ""),
			paste0("Question : ", input$response_name, " diff\u00e8re-t-il selon ", paste(nm, collapse = ", "), " une fois ",
				paste(cv, collapse = ", "), " pris en compte ?", if (length(slopes)) paste0(" Et les pentes diff\u00e8rent-elles (", paste(slopes, collapse = ", "), ") ?") else ""), "decide")
	} else if (input$plan %in% c("paired", "repeated")) {
		ns <- nlevels(input$id)
		lv <- vapply(fac, nlevels, integer(1))
		wi <- intersect(input$within, names(fac)); bw <- setdiff(names(fac), wi)
		desc_en <- paste0(nm, " (", lv, " levels", ifelse(names(fac) %in% wi, ", within subject", ", between subjects"), ")", collapse = ", ")
		desc_fr <- paste0(nm, " (", lv, " niveaux", ifelse(names(fac) %in% wi, ", intra-sujet", ", entre sujets"), ")", collapse = ", ")
		if (identical(input$plan, "paired")) {
			.bay_tell(st, paste0("Plan recognised: paired data, ", nm[1L], " (", paste(levels(fac[[1L]]), collapse = ", "), ") measured on the same ", ns, " subjects; N = ", N, "."),
				paste0("Plan reconnu : donn\u00e9es appari\u00e9es, ", nm[1L], " (", paste(levels(fac[[1L]]), collapse = ", "), ") mesur\u00e9 sur les m\u00eames ", ns, " sujets ; N = ", N, "."), "step")
		} else {
			.bay_tell(st, paste0("Plan recognised: repeated measures on ", ns, " subjects, ", desc_en, "; N = ", N, "."),
				paste0("Plan reconnu : mesures r\u00e9p\u00e9t\u00e9es sur ", ns, " sujets, ", desc_fr, " ; N = ", N, "."), "step")
		}
		terms <- gsub("`", "", input$terms, fixed = TRUE)
		if (length(terms) == 1L) {
			.bay_tell(st, paste0("One question: is there a difference in ", input$response_name, " according to ", nm[1L], ", each subject compared with itself?"),
				paste0("Une seule question : y a-t-il une diff\u00e9rence de ", input$response_name, " selon ", nm[1L], ", chaque sujet compar\u00e9 \u00e0 lui-m\u00eame ?"), "decide")
		} else {
			sp <- .bay_fact_space(terms)
			.bay_tell(st, paste0("Terms: ", paste(terms, collapse = ", "), "; one answer per term, obtained by weighing the ", length(sp$models), " combinations of terms together."),
				paste0("Termes : ", paste(terms, collapse = ", "), " ; une r\u00e9ponse par terme, obtenue en pesant ensemble les ", length(sp$models), " combinaisons de termes."), "decide")
		}
	} else {
		lab <- list(
			en = c(ancova = "covariates (ANCOVA)", paired = "paired data", repeated = "repeated measures",
				mixed = "mixed model", manova = "several responses (MANOVA)"),
			fr = c(ancova = "covariables (ANCOVA)", paired = "donn\u00e9es appari\u00e9es", repeated = "mesures r\u00e9p\u00e9t\u00e9es",
				mixed = "mod\u00e8le mixte", manova = "plusieurs r\u00e9ponses (MANOVA)"))
		.bay_tell(st, paste0("Plan recognised: ", lab$en[[input$plan]], "; N = ", N, "."),
			paste0("Plan reconnu : ", lab$fr[[input$plan]], " ; N = ", N, "."), "step")
		.bay_tell(st, "This plan is not validated in m.bay() yet: refused rather than answered wrongly.",
			"Ce plan n'est pas encore valid\u00e9 dans m.bay() : refus\u00e9 plut\u00f4t que mal trait\u00e9.", "decide")
	}
	invisible(NULL)
}

# Step 5, one factor: what the sensitivity does to the leading hypothesis.
.bay_story_sensitivity <- function(st, sens, hyp, het) {
	if (is.null(st)) return(invisible(NULL))
	if (is.null(sens)) {
		.bay_tell(st, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
		.bay_tell(st, "The conclusion is not tested against other starting points.",
			"La conclusion n'est pas \u00e9prouv\u00e9e contre d'autres points de d\u00e9part.", "decide")
		return(invisible(NULL))
	}
	r1 <- .bay_num_txt(min(sens$r)); r2 <- .bay_num_txt(max(sens$r))
	u1 <- .bay_num_txt(min(sens$u)); u2 <- .bay_num_txt(max(sens$u))
	o <- sens$odds[c(1L, length(sens$odds))]
	rg <- sens$range_lead
	.bay_tell(st, paste0("Sensitivity: expected size r from ", r1$en, " to ", r2$en, ", starting chances of a difference from 1 against 2 to 2 against 1",
			if (het) ", spread of the variances between groups", ". \"", hyp$label_en[1L], "\" then goes from ",
			.bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "."),
		paste0("Sensibilit\u00e9 : taille attendue r de ", r1$fr, " \u00e0 ", r2$fr, ", chances de d\u00e9part d'une diff\u00e9rence de 1 contre 2 \u00e0 2 contre 1",
			if (het) ", r\u00e9partition des variances entre groupes", ". \u00ab ", hyp$label_fr[1L], " \u00bb va alors de ",
			.bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), "."), fun = "m.bay() {KefiR}")
	if (identical(sens$verdict, "changes")) {
		.bay_tell(st, "Caution added: the leading hypothesis can reverse with the starting point; the conclusion will not assert it.",
			"Prudence ajout\u00e9e : l'hypoth\u00e8se en t\u00eate peut s'inverser selon le point de d\u00e9part ; la conclusion ne l'affirmera pas.", "decide")
	} else if (identical(sens$verdict, "crosses")) {
		.bay_tell(st, "The order holds; only the degree of support (weak, moderate, strong) can change.",
			"L'ordre tient ; seul le degr\u00e9 de soutien (faible, mod\u00e9r\u00e9, fort) peut changer.", "decide")
	} else {
		.bay_tell(st, "The answer does not depend on these choices.",
			"La r\u00e9ponse ne d\u00e9pend pas de ces choix.", "decide")
	}
	invisible(NULL)
}

# Step 6, one factor: which groups differ, and by how much.
.bay_story_posthoc <- function(st, ph, caution = FALSE, subjects = FALSE) {
	if (is.null(st) || is.null(ph)) return(invisible(NULL))
	if (identical(ph$method, "paired")) {
		.bay_tell(st, "By how much the two levels differ.", "De combien les deux niveaux diff\u00e8rent.", "step",
			fun = "posterior(ttestBF()) {BayesFactor}, hdi() {bayestestR}")
	} else if (subjects) {
		.bay_tell(st, "Which conditions differ, and by how much.", "Quelles conditions diff\u00e8rent, et de combien.", "step",
			fun = "lmBF(whichRandom = id), posterior() {BayesFactor}, hdi() {bayestestR}")
	} else {
		.bay_tell(st, "Which groups differ, and by how much.", "Quels groupes diff\u00e8rent, et de combien.", "step",
			fun = "lmBF(), posterior() {BayesFactor}, hdi() {bayestestR}")
	}
	pr <- ph$pairs
	ci_txt <- function(x, lo, hi, dc) {
		f <- function(v) gsub(".", dc, format(signif(v, 2), trim = TRUE), fixed = TRUE)
		paste0(f(x), " [", f(lo), if (dc == ",") " ; " else "; ", f(hi), "]")
	}
	big <- function(n) formatC(n, format = "d", big.mark = " ")
	if (ph$k > 2L && identical(ph$method, "gibbs")) {
		mc <- max(1L, ceiling(200 * ph$mc_error))
		.bay_tell(st, paste0("Change of strategy: ", big(ph$n_partitions), " possible groupings, too many to compute them all (more than ",
				ph$max_exact, " groups); they are explored by draws (Gibbs sampler, Gopalan and Berry 1998), ",
				big(ph$sweeps), " rounds; the probabilities carry a sampling error of about ", mc, " percentage point(s)."),
			paste0("Changement de strat\u00e9gie : ", big(ph$n_partitions), " regroupements possibles, trop pour les calculer tous (plus de ",
				ph$max_exact, " groupes) ; ils sont explor\u00e9s par tirage (\u00e9chantillonneur de Gibbs, Gopalan et Berry 1998), ",
				big(ph$sweeps), " tours ; les probabilit\u00e9s portent une erreur de tirage d'environ ", mc, " point(s)."), "decide")
	}
	if (ph$k > 2L) {
		top <- ph$partitions[seq_len(min(3L, nrow(ph$partitions))), , drop = FALSE]
		n_en <- if (identical(ph$method, "gibbs")) paste0(big(nrow(ph$partitions)), " visited") else paste0(big(nrow(ph$partitions)), " compared")
		n_fr <- if (identical(ph$method, "gibbs")) paste0(big(nrow(ph$partitions)), " visit\u00e9s") else paste0(big(nrow(ph$partitions)), " compar\u00e9s")
		.bay_tell(st, paste0("Most probable groupings (", n_en, "): ",
				paste0(top$partition, " (", vapply(top$posterior, .bay_pct, character(1), dec = "."), ")", collapse = "; "), "."),
			paste0("Regroupements les plus probables (", n_fr, ") : ",
				paste0(top$partition, " (", vapply(top$posterior, .bay_pct, character(1), dec = ","), ")", collapse = " ; "), "."))
	}
	unit <- .bay_delta_unit(ph)
	shown <- .bay_posthoc_shown(pr)
	for (i in shown) {
		lead_en <- if (ph$k > 2L) paste0(pr$group1[i], " and ", pr$group2[i], ": different at ", .bay_pct(pr$p_differ[i], "."), "; estimated gap ") else "Estimated gap "
		lead_fr <- if (ph$k > 2L) paste0(pr$group1[i], " et ", pr$group2[i], " : diff\u00e9rents \u00e0 ", .bay_pct(pr$p_differ[i], ","), " ; \u00e9cart estim\u00e9 ") else "\u00c9cart estim\u00e9 "
		.bay_tell(st, paste0(lead_en, pr$group2[i], " - ", pr$group1[i], " = ", ci_txt(pr$diff[i], pr$diff_low[i], pr$diff_high[i], "."),
				", i.e. ", ci_txt(pr$delta[i], pr$delta_low[i], pr$delta_high[i], "."), " ", unit$en, "."),
			paste0(lead_fr, pr$group2[i], " - ", pr$group1[i], " = ", ci_txt(pr$diff[i], pr$diff_low[i], pr$diff_high[i], ","),
				", soit ", ci_txt(pr$delta[i], pr$delta_low[i], pr$delta_high[i], ","), " ", unit$fr, "."))
	}
	if (length(shown) < nrow(pr)) {
		.bay_tell(st, paste0("The other ", nrow(pr) - length(shown), " pairs lean towards alike (different below 50 %); all pairs: $posthoc$pairs."),
			paste0("Les ", nrow(pr) - length(shown), " autres paires penchent vers semblables (diff\u00e9rentes \u00e0 moins de 50 %) ; toutes les paires : $posthoc$pairs."), "info")
	}
	if (identical(ph$delta_unit %in% c("within", "adjusted", "diff"), TRUE)) {
		.bay_tell(st, "Estimated gap: posterior median, pulled slightly towards zero by the prior. Brackets: 95 % credible interval, if a difference exists (the value lies there with 95 % probability; not a confidence interval).",
			"\u00c9cart estim\u00e9 : m\u00e9diane a posteriori, l\u00e9g\u00e8rement ramen\u00e9e vers z\u00e9ro par l'a priori. Crochets : intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %, si une diff\u00e9rence existe (la valeur s'y trouve avec 95 % de probabilit\u00e9 ; ce n'est pas un intervalle de confiance).", "info")
	} else .bay_tell(st, "Estimated gap: posterior median, pulled slightly towards zero by the prior, hence a little smaller than the observed gap of step 7. Brackets: 95 % credible interval, if a difference exists (the value lies there with 95 % probability; not a confidence interval).",
		"\u00c9cart estim\u00e9 : m\u00e9diane a posteriori, l\u00e9g\u00e8rement ramen\u00e9e vers z\u00e9ro par l'a priori, donc un peu plus petite que l'\u00e9cart observ\u00e9 de l'\u00e9tape 7. Crochets : intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %, si une diff\u00e9rence existe (la valeur s'y trouve avec 95 % de probabilit\u00e9 ; ce n'est pas un intervalle de confiance).", "info")
	if (ph$k > 2L) {
		.bay_tell(st, "All pairs come from one single law over the groupings: no correction pair by pair is needed, and the answers cannot contradict each other.",
			"Toutes les paires viennent d'une seule loi sur les regroupements : aucune correction paire par paire n'est n\u00e9cessaire, et les r\u00e9ponses ne peuvent pas se contredire.", "decide")
	}
	if (isTRUE(ph$caution_overlap)) {
		.bay_tell(st, "These gaps compare the groups at the mean of the covariate, where some have few or no data (step 2): read them with caution.",
			"Ces \u00e9carts comparent les groupes \u00e0 la moyenne de la covariable, l\u00e0 o\u00f9 certains ont peu ou pas de donn\u00e9es (\u00e9tape 2) : les lire avec r\u00e9serve.", "decide")
	}
	if (caution) {
		.bay_tell(st, "These pairs assume a normal noise with a common variance, which step 4 did not find: read them with caution.",
			"Ces paires supposent un bruit normal \u00e0 variance commune, ce que l'\u00e9tape 4 n'a pas trouv\u00e9 : les lire avec r\u00e9serve.", "decide")
	}
	invisible(NULL)
}

# Pairs shown in the story and the bilan: all of them up to 15; beyond, those that lean
# towards different, most probable first (the others are summed up in one line).
.bay_posthoc_shown <- function(pr) {
	if (nrow(pr) <= 15L || anyNA(pr$p_differ)) return(seq_len(nrow(pr)))
	ix <- which(pr$p_differ >= 0.5)
	ix[order(-pr$p_differ[ix])]
}

# Unit of the standardised difference: pooled SD between groups, residual SD within subject
# (subject's own level removed), or SD of the paired differences (dz).
.bay_delta_unit <- function(ph) {
	switch(if (is.null(ph$delta_unit)) "pooled" else ph$delta_unit,
		within = list(en = "residual SD within subject", fr = "\u00e9cart-type r\u00e9siduel au sein des sujets"),
		adjusted = list(en = "residual SD (covariates accounted for)", fr = "\u00e9cart-type r\u00e9siduel (covariables prises en compte)"),
		diff = list(en = "SD of the differences (dz)", fr = "\u00e9cart-type des diff\u00e9rences (dz)"),
		list(en = "SD", fr = "\u00e9cart-type"))
}

# Steps 2 to 4 for paired data and repeated measures (D15).
.bay_story_rm <- function(st, input, y, paired, n_dropped, dropped_ids = character(0)) {
	if (is.null(st)) return(invisible(NULL))
	id <- input$id
	fac <- input$factors
	ns <- nlevels(id)
	if (n_dropped > 0L) {
		lst <- paste(utils::head(dropped_ids, 10L), collapse = ", ")
		if (length(dropped_ids) > 10L) lst <- paste0(lst, ", ...")
		.bay_tell(st, paste0(n_dropped, " subject(s) without a measure at both levels removed (", lst, "): a pair needs both."),
			paste0(n_dropped, " sujet(s) sans mesure aux deux niveaux retir\u00e9(s) (", lst, ") : une paire demande les deux."), "decide")
	}
	cells <- interaction(fac, drop = TRUE)
	wi <- intersect(input$within, names(fac))
	wc <- if (length(wi)) interaction(fac[wi], drop = TRUE) else cells
	full <- sum(tapply(wc, id, function(v) length(unique(v))) == nlevels(wc))
	.bay_tell(st, paste0("Subjects: ", ns, "; measured in every condition within subject: ", full, "."),
		paste0("Sujets : ", ns, " ; mesur\u00e9s dans toutes les conditions intra-sujet : ", full, "."))
	if (full < ns) {
		.bay_tell(st, "Incomplete subjects are kept: the model with a subject effect uses all their measures.",
			"Les sujets incomplets sont gard\u00e9s : le mod\u00e8le \u00e0 effet sujet utilise toutes leurs mesures.", "decide")
	}
	dev <- y - stats::ave(y, id) + mean(y)
	.bay_tell(st, "Values far from their condition are sought once each subject's own level is removed.",
		"Les valeurs \u00e9loign\u00e9es sont cherch\u00e9es une fois retir\u00e9 le niveau propre de chaque sujet.", "info")
	chk <- .bay_story_checks(st, dev, cells, "condition", "condition")
	.bay_tell(st, "No value is removed.", "Aucune valeur n'est retir\u00e9e.", "decide")
	.bay_tell(st, "Shape of the noise.", "Forme du bruit.", "step")
	if (paired) {
		.bay_tell(st, "Paired data: the comparison rests on the difference within each subject, which removes what varies from one subject to another.",
			"Donn\u00e9es appari\u00e9es : la comparaison porte sur la diff\u00e9rence au sein de chaque sujet, ce qui retire ce qui varie d'un sujet \u00e0 l'autre.")
	} else {
		.bay_tell(st, "Repeated measures: each subject keeps its own level (random subject effect, present in every model compared, \"no difference\" included).",
			"Mesures r\u00e9p\u00e9t\u00e9es : chaque sujet garde son niveau propre (effet sujet al\u00e9atoire, pr\u00e9sent dans tous les mod\u00e8les compar\u00e9s, \u00ab aucune diff\u00e9rence \u00bb compris).")
	}
	if (!paired) {
		.bay_tell(st, "A random subject effect assumes that all conditions are equally correlated within a subject (compound symmetry); this is not tested here (m.test() checks sphericity with Mauchly's test).",
			"Un effet sujet al\u00e9atoire suppose que toutes les conditions sont \u00e9galement corr\u00e9l\u00e9es au sein d'un sujet (sym\u00e9trie compos\u00e9e) ; ce n'est pas test\u00e9 ici (m.test() contr\u00f4le la sph\u00e9ricit\u00e9 par le test de Mauchly).")
	}
	.bay_tell(st, "Only the normal law is computed: heavy tails and unequal variances between conditions are not validated for subjects yet, so they are not compared.",
		"Seule la loi normale est calcul\u00e9e : queues lourdes et variances in\u00e9gales entre conditions ne sont pas encore valid\u00e9es avec des sujets, elles ne sont donc pas compar\u00e9es.")
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(st, paste0(far, " value(s) lie beyond 3 interquartile ranges of their condition: the answer assumes a normal noise, read it with caution."),
			paste0(far, " valeur(s) au-del\u00e0 de 3 \u00e9carts interquartiles de leur condition : la r\u00e9ponse suppose un bruit normal, la lire avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(st, "No value lies far enough to put the normal law in doubt at first sight; it is not tested against the others.",
			"Aucune valeur n'est assez \u00e9loign\u00e9e pour mettre d'embl\u00e9e la loi normale en doute ; elle n'est pas mise \u00e0 l'\u00e9preuve des autres.", "decide")
	}
	invisible(chk)
}

# Steps 2 and 3 for the analysis of covariance (D16).
.bay_story_ancova <- function(st, input, y) {
	fac <- input$factors
	cov <- input$covariates
	cells <- interaction(fac, drop = TRUE)
	res <- stats::residuals(stats::lm(y ~ ., data = data.frame(y = y, cov, check.names = FALSE)))
	.bay_tell(st, "Values far from their group are sought on what the covariates do not explain.",
		"Les valeurs \u00e9loign\u00e9es sont cherch\u00e9es sur ce que les covariables n'expliquent pas.", "info")
	chk <- .bay_story_checks(st, res + mean(y), cells)
	ov <- .bay_cov_overlap(cov, cells)
	for (i in seq_len(nrow(ov))) {
		v <- ov$covariate[i]
		lo <- .bay_num_txt(ov$low[i]); hi <- .bay_num_txt(ov$high[i]); sh <- .bay_pct(ov$min_share[i], ","); she <- .bay_pct(ov$min_share[i], ".")
		if (ov$low[i] > ov$high[i]) {
			.bay_tell(st, paste0("Covariate ", v, ": the groups do not overlap at all."),
				paste0("Covariable ", v, " : les groupes ne se recouvrent pas du tout."))
		} else {
			.bay_tell(st, paste0("Covariate ", v, ": common range of the groups from ", lo$en, " to ", hi$en, "; it holds at least ", she, " of each group (", ov$worst[i], " the least)."),
				paste0("Covariable ", v, " : \u00e9tendue commune aux groupes de ", lo$fr, " \u00e0 ", hi$fr, " ; elle contient au moins ", sh, " de chaque groupe (", ov$worst[i], " le moins)."))
		}
	}
	if (any(ov$poor)) {
		.bay_tell(st, "Caution added: the groups hardly overlap on the covariate, so the adjustment compares them where some have no data (extrapolation); the support of the answers on the groups will be capped at weak.",
			"Prudence ajout\u00e9e : les groupes se recouvrent peu sur la covariable, l'ajustement les compare donc l\u00e0 o\u00f9 certains n'ont pas de donn\u00e9es (extrapolation) ; le soutien des r\u00e9ponses sur les groupes sera plafonn\u00e9 \u00e0 faible.", "decide")
	}
	ctr <- vapply(cov, mean, numeric(1))
	ct <- .bay_num_txt(ctr)
	.bay_tell(st, paste0("Covariate(s) centred on their mean (", paste0(names(cov), " = ", ct$en, collapse = ", "), "): the groups are compared at this value, inside the data, not at 0 (Aiken and West 1991)."),
		paste0("Covariable(s) centr\u00e9e(s) sur leur moyenne (", paste0(names(cov), " = ", ct$fr, collapse = ", "), ") : les groupes sont compar\u00e9s \u00e0 cette valeur, dans les donn\u00e9es, et non \u00e0 0 (Aiken et West 1991)."), "decide")
	.bay_tell(st, "No value is removed.", "Aucune valeur n'est retir\u00e9e.", "decide")
	.bay_tell(st, "Shape of the noise.", "Forme du bruit.", "step")
	.bay_tell(st, "Only the normal law with a common variance is computed, and the covariate acts in a straight line (linearity not tested): heavy tails, unequal variances and curvature are not validated with covariates yet.",
		"Seule la loi normale \u00e0 variance commune est calcul\u00e9e, et la covariable agit en ligne droite (lin\u00e9arit\u00e9 non test\u00e9e) : queues lourdes, variances in\u00e9gales et courbure ne sont pas encore valid\u00e9es avec des covariables.")
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(st, paste0(far, " value(s) lie beyond 3 interquartile ranges once the covariates are accounted for: read the answer with caution."),
			paste0(far, " valeur(s) au-del\u00e0 de 3 \u00e9carts interquartiles une fois les covariables prises en compte : lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(st, "No value lies far enough to put the normal law in doubt at first sight; it is not tested against the others.",
			"Aucune valeur n'est assez \u00e9loign\u00e9e pour mettre d'embl\u00e9e la loi normale en doute ; elle n'est pas mise \u00e0 l'\u00e9preuve des autres.", "decide")
	}
	invisible(ov)
}

# Step 6 when the slopes differ (D16): gaps at the quartiles of the covariate, and where
# the groups differ.
.bay_story_slopes <- function(st, ph, xname, term, p_int) {
	if (is.null(st)) return(invisible(NULL))
	.bay_tell(st, paste0("Which groups differ, and by how much, along ", xname, "."),
		paste0("Quels groupes diff\u00e8rent, et de combien, selon ", xname, "."), "step",
		fun = "lmBF(y ~ g * x), posterior() {BayesFactor}, hdi() {bayestestR}")
	.bay_tell(st, paste0("Change of strategy: slopes that differ are probable (", term, " ", .bay_pct(p_int, "."),
			"), so the gap between groups depends on ", xname, ": one adjusted gap would mean nothing. Gaps at the quartiles of ", xname,
			", and the zone where the groups differ (Bayesian counterpart of Johnson and Neyman)."),
		paste0("Changement de strat\u00e9gie : des pentes diff\u00e9rentes sont probables (", term, " ", .bay_pct(p_int, ","),
			"), donc l'\u00e9cart entre groupes d\u00e9pend de ", xname, " : un \u00e9cart ajust\u00e9 unique n'aurait pas de sens. \u00c9carts aux quartiles de ", xname,
			", et zone o\u00f9 les groupes diff\u00e8rent (\u00e9quivalent bay\u00e9sien de Johnson et Neyman)."), "decide")
	for (L in .bay_slopes_lines(ph, xname)) .bay_tell(st, L$en, L$fr)
	.bay_tell(st, "Brackets: 95 % credible interval. Zone: values of the observed range where this interval excludes 0.",
		"Crochets : intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %. Zone : valeurs de l'\u00e9tendue observ\u00e9e o\u00f9 cet intervalle exclut 0.", "info")
	invisible(NULL)
}

.bay_slopes_lines <- function(ph, xname) {
	f <- function(v, dc) gsub(".", dc, format(signif(v, 2), trim = TRUE), fixed = TRUE)
	sep <- function(dc) if (dc == ",") " ; " else "; "
	out <- list()
	for (j in seq_len(nrow(ph$pairs))) {
		a <- ph$pairs$group1[j]; b <- ph$pairs$group2[j]
		gp <- ph$gaps[ph$gaps$group1 == a & ph$gaps$group2 == b, , drop = FALSE]
		txt <- function(dc) paste0(vapply(seq_len(nrow(gp)), function(q) paste0(xname, " = ", f(gp$at[q], dc), " : ", f(gp$diff[q], dc),
			" [", f(gp$diff_low[q], dc), sep(dc), f(gp$diff_high[q], dc), "]"), character(1)), collapse = sep(dc))
		z <- ph$zones[[j]]
		zt <- function(dc, none, word, and) if (!nrow(z)) none else paste(paste0(word, " ", f(z$from, dc), " ", and, " ", f(z$to, dc)), collapse = if (dc == ",") " et " else " and ")
		f <- function(v, dc) vapply(v, function(w) gsub(".", dc, format(signif(w, 2), trim = TRUE), fixed = TRUE), character(1))
		out[[j]] <- list(
			en = paste0(b, " - ", a, ": ", txt("."), "; they differ ", zt(".", "nowhere in the observed range", "from", "to"), "."),
			fr = paste0(b, " - ", a, " : ", txt(","), " ; ils diff\u00e8rent ", zt(",", "nulle part dans l'\u00e9tendue observ\u00e9e", "de", "\u00e0"), "."))
	}
	out
}
