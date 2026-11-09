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
