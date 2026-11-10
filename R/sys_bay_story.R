# sys_bay_story.R
# The story told by m.bay(verbose = TRUE), in the manner of m.test(): numbered steps,
# "==>" for a finding, "-->" for a decision, "Note :" for anything else; the function and
# its package in brackets on the title line; the result of the step in its title
# ("4) Shapes of the noise [m.bay() {KefiR}] - normal, equal variances ahead (51.9 %).").
# Readability (relecture/lisibilite_2026-11-09): every line holds within
# min(getOption("width"), 80) columns, a tab counted as 8, long text folded with strwrap()
# and a hanging indent. The lines of a step are held until the next step, so that its
# result can go into its title; they are kept in the returned object (story), so that
# print(result, story = TRUE) replays them. Lines of kind "detail" (method, references,
# jargon) are kept apart and shown by print(result, detail = TRUE).

.bay_story <- function(verbose) {
	st <- new.env(parent = emptyenv())
	st$k <- 0L
	st$en <- character(0)
	st$fr <- character(0)
	st$detail_en <- character(0)
	st$detail_fr <- character(0)
	st$verbose <- isTRUE(verbose)
	st$pending <- list()
	st$cur <- NULL
	st
}

# Width of a line: the console, at most 80 columns (Antoine, 10 November 2026).
.bay_width <- function() max(40L, min(as.integer(getOption("width", 80L)), 80L))

# Visible width, a tab counted as 8 columns.
.bay_vis <- function(s) nchar(gsub("\t", "        ", s, fixed = TRUE), type = "width")

# Folds txt behind a lead ("\t==> "); the following lines get cont, of the same width.
.bay_wrap <- function(txt, lead, cont = NULL) {
	if (is.null(cont)) cont <- paste0(if (startsWith(lead, "\t")) "\t" else "", strrep(" ", .bay_vis(sub("^\t", "", lead))))
	avail <- max(20L, .bay_width() - .bay_vis(lead))
	# French spacing and intervals are not split: non-breaking space before % : ; ! ? and
	# inside guillemets and brackets, turned back into a space once folded.
	nb <- "\u00a0"
	txt <- gsub(" ([%:;!?])", paste0(nb, "\\1"), txt)
	txt <- gsub("\u00ab ", paste0("\u00ab", nb), txt)
	txt <- gsub(" \u00bb", paste0(nb, "\u00bb"), txt)
	txt <- gsub("(\\[[^] ]*) ([^] ]*\\])", paste0("\\1", nb, "\\2"), txt)
	txt <- gsub("(\\[[^] ]*) ([^] ]*) ([^] ]*\\])", paste0("\\1", nb, "\\2", nb, "\\3"), txt)
	txt <- gsub(" ([=<>x]) ", paste0(nb, "\\1", nb), txt)
	p <- strwrap(txt, width = avail + 1L)
	p <- gsub(nb, " ", p, fixed = TRUE)
	if (!length(p)) p <- ""
	c(paste0(lead, p[1L]), if (length(p) > 1L) paste0(cont, p[-1L]))
}

.bay_lead <- function(kind, lang) {
	switch(kind,
		fact = "\t==> ", decide = "\t--> ",
		note = , info = if (lang == "fr") "\tNote : " else "\tNote: ",
		sub = "\t    ", plain = "\t", "\t")
}

# kind: "step" (numbered title), "fact" (==>), "decide" (-->), "note" or "info" (Note :),
# "detail" (kept for detail = TRUE, not told). fun: the function and package, c(en, fr)
# possible; on a line that is not a title it goes to the title of the current step.
.bay_tell <- function(st, en, fr, kind = "fact", fun = NULL) {
	if (is.null(st)) return(invisible(NULL))
	if (identical(kind, "detail")) {
		st$detail_en <- c(st$detail_en, en)
		st$detail_fr <- c(st$detail_fr, fr)
		return(invisible(NULL))
	}
	if (identical(kind, "step")) {
		.bay_story_flush(st)
		st$k <- st$k + 1L
		st$cur <- list(k = st$k, en = en, fr = fr, fun_en = character(0), fun_fr = character(0),
			res_en = NULL, res_fr = NULL, lines = list())
		if (!is.null(fun)) .bay_step_fun(st, fun)
		if (length(st$pending)) {
			p <- st$pending
			st$pending <- list()
			for (x in p) .bay_tell(st, x$en, x$fr, x$kind)
		}
		return(invisible(NULL))
	}
	if (is.null(st$cur)) {
		st$pending[[length(st$pending) + 1L]] <- list(en = en, fr = fr, kind = kind)
		return(invisible(NULL))
	}
	if (!is.null(fun)) .bay_step_fun(st, fun)
	st$cur$lines[[length(st$cur$lines) + 1L]] <- list(kind = kind, en = en, fr = fr)
	invisible(NULL)
}

.bay_step_fun <- function(st, fun) {
	if (is.null(st$cur)) return(invisible(NULL))
	st$cur$fun_en <- unique(c(st$cur$fun_en, fun[1L]))
	st$cur$fun_fr <- unique(c(st$cur$fun_fr, fun[length(fun)]))
	invisible(NULL)
}

# The result of the current step, said in its title after " - ".
.bay_step_result <- function(st, en, fr) {
	if (is.null(st) || is.null(st$cur)) return(invisible(NULL))
	st$cur$res_en <- en
	st$cur$res_fr <- fr
	invisible(NULL)
}

# Renders a step in one language: its title, then its lines, folded.
.bay_render_step <- function(cur, lang) {
	title <- sub("[.:]$", "", if (lang == "fr") cur$fr else cur$en)
	fun <- if (lang == "fr") cur$fun_fr else cur$fun_en
	res <- if (lang == "fr") cur$res_fr else cur$res_en
	head <- paste0(title, if (length(fun)) paste0(" [", paste(fun, collapse = ", "), "]"),
		if (!is.null(res)) paste0(" - ", sub("[.]$", "", res)), ".")
	lead <- paste0(cur$k, ") ")
	out <- .bay_wrap(head, lead, strrep(" ", nchar(lead)))
	for (x in cur$lines) out <- c(out, .bay_wrap(.bay_note_case(if (lang == "fr") x$fr else x$en, x$kind), .bay_lead(x$kind, lang)))
	out
}

# After "Note :", a lowercase first letter (but not for a name such as "Gibbs" or "R").
.bay_note_case <- function(txt, kind) {
	if (!kind %in% c("note", "info")) return(txt)
	if (grepl("^[A-Z\u00c0-\u00dd][a-z\u00e0-\u00ff]", txt) && !grepl("^(Gibbs|Gopalan|Rouder|Tukey|Bayes)", txt)) {
		substr(txt, 1L, 1L) <- tolower(substr(txt, 1L, 1L))
	}
	txt
}

# Notes go to the end of their step, at most two; the others are kept for detail = TRUE
# (audit of D18 by the parent session, point 1).
.bay_notes_last <- function(lines, max_notes = 2L, spill = NULL) {
	is_note <- vapply(lines, function(x) x$kind %in% c("note", "info"), logical(1))
	notes <- lines[is_note]
	if (length(notes) > max_notes && is.function(spill)) for (x in notes[-seq_len(max_notes)]) spill(x)
	c(lines[!is_note], utils::head(notes, max_notes))
}

# Prints and keeps the current step.
.bay_story_flush <- function(st) {
	if (is.null(st) || is.null(st$cur)) return(invisible(NULL))
	st$cur$lines <- .bay_notes_last(st$cur$lines, 2L, function(x) {
		st$detail_en <- c(st$detail_en, x$en)
		st$detail_fr <- c(st$detail_fr, x$fr)
	})
	en <- .bay_render_step(st$cur, "en")
	fr <- .bay_render_step(st$cur, "fr")
	st$en <- c(st$en, en)
	st$fr <- c(st$fr, fr)
	if (st$verbose) cat(paste0(.msg(en, fr), "\n"), sep = "")
	st$cur <- NULL
	invisible(NULL)
}

# A note said before the first step waits for it.
.bay_tell_note <- function(st, en, fr, kind = "fact") .bay_tell(st, en, fr, kind)

.bay_story_lines <- function(st) {
	.bay_story_flush(st)
	list(en = st$en, fr = st$fr, k = st$k, detail_en = st$detail_en, detail_fr = st$detail_fr)
}

# Replays a stored story.
.bay_story_replay <- function(story) {
	if (is.null(story) || !length(story$en)) return(invisible(NULL))
	txt <- if (identical(.msg("en", "fr"), "fr")) story$fr else story$en
	cat(paste0(txt, "\n"), sep = "")
	invisible(NULL)
}

# The bilan in short blocks separated by an empty line (hypotheses; terms and pairs;
# observations and report), each with at most two notes, at its end. emit(kind, en, fr);
# emit("block", "", "") opens the next block; emit("flush", "", "") prints.
.bay_block_emitter <- function() {
	out <- .bay_emitter()
	blocks <- list(list())
	function(kind, en, fr) {
		if (identical(kind, "block")) {
			if (length(blocks[[length(blocks)]])) blocks[[length(blocks) + 1L]] <<- list()
			return(invisible(NULL))
		}
		if (identical(kind, "flush")) {
			first <- TRUE
			for (b in blocks) {
				if (!length(b)) next
				if (!first) cat("\n")
				first <- FALSE
				for (x in .bay_notes_last(b, 2L)) out(x$kind, x$en, x$fr)
			}
			blocks <<- list(list())
			return(invisible(NULL))
		}
		if (identical(kind, "title")) return(out(kind, en, fr))
		blocks[[length(blocks)]][[length(blocks[[length(blocks)]]) + 1L]] <<- list(kind = kind, en = en, fr = fr)
		invisible(NULL)
	}
}

# Lines of the bilan and of the detail section, folded the same way. emit(kind, en, fr).
.bay_emitter <- function() {
	lang <- if (identical(.msg("en", "fr"), "fr")) "fr" else "en"
	function(kind, en, fr) {
		txt <- if (lang == "fr") fr else en
		lines <- if (identical(kind, "title")) .bay_wrap(txt, "", "   ") else .bay_wrap(.bay_note_case(txt, kind), .bay_lead(kind, lang))
		cat(paste0(lines, "\n"), sep = "")
		invisible(NULL)
	}
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
	uf <- function(k) if (k == 1L) unit_fr else paste0(unit_fr, "s")
	ue <- function(k) if (k == 1L) unit_en else paste0(unit_en, "s")
	.bay_tell(st, paste0("Sizes: ", min(n), if (max(n) > min(n)) paste0(" to ", max(n)), " per ", unit_en, " (N = ", N, ")",
			if (max(n) > min(n)) ", unbalanced." else "."),
		paste0("Effectifs : ", min(n), if (max(n) > min(n)) paste0(" \u00e0 ", max(n)), " par ", unit_fr, " (N = ", N, ")",
			if (max(n) > min(n)) ", d\u00e9s\u00e9quilibr\u00e9." else "."))
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
		k <- length(const)
		.bay_tell(st, paste0("Constant ", ue(k), ": ", paste(const, collapse = ", "), "."),
			paste0(if (unit_fr == "case") (if (k == 1L) "Case constante : " else "Cases constantes : ") else
				(if (k == 1L) "Groupe constant : " else "Groupes constants : "), paste(const, collapse = ", "), "."))
		.bay_tell(st, "To check: a value repeated through a whole group often signals a ceiling or an entry error; kept.",
			"\u00c0 v\u00e9rifier : une valeur r\u00e9p\u00e9t\u00e9e dans tout un groupe signale souvent un plafond ou une erreur de saisie ; gard\u00e9e.", "decide")
	}
	ex <- .bay_extremes(y, g)
	.bay_tell(st, paste0("Values far from their ", unit_en, ": ", ex[["mild"]], " beyond 1.5 interquartile ranges, ",
			ex[["far"]], " beyond 3 (of ", N, ")."),
		paste0("Valeurs \u00e9loign\u00e9es de leur ", unit_fr, " : ", ex[["mild"]], " au-del\u00e0 de 1,5 \u00e9cart interquartile, ",
			ex[["far"]], " au-del\u00e0 de 3 (sur ", N, ")."), fun = "quantile() {stats}")
	.bay_tell(st, "Far values: Tukey's rule, 1.5 and 3 interquartile ranges beyond the quartiles of each group of at least 4 values.",
		"Valeurs \u00e9loign\u00e9es : r\u00e8gle de Tukey, 1,5 et 3 \u00e9carts interquartiles au-del\u00e0 des quartiles de chaque groupe d'au moins 4 valeurs.", "detail")
	.bay_step_result(st, if (ex[["far"]] > 0L) paste0(ex[["far"]], if (ex[["far"]] == 1L) " far value, kept" else " far values, kept") else "nothing to remove",
		if (ex[["far"]] > 0L) paste0(ex[["far"]], if (ex[["far"]] == 1L) " valeur \u00e9loign\u00e9e, gard\u00e9e" else " valeurs \u00e9loign\u00e9es, gard\u00e9es") else "rien \u00e0 retirer")
	list(n = n, tied = tied, constant = const, extremes = ex)
}

# Words of the stability of the answer, for the title of step 5.
.bay_stab_txt <- function(verdict) {
	switch(if (is.null(verdict)) "none" else verdict,
		changes = list(en = "the order can reverse", fr = "l'ordre peut s'inverser"),
		crosses = list(en = "the order does not change", fr = "l'ordre ne change pas"),
		stable = list(en = "the order does not change", fr = "l'ordre ne change pas"),
		list(en = "sensitivity not computed", fr = "sensibilit\u00e9 non calcul\u00e9e"))
}

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
		yv <- input$response
		if (isTRUE(input$categorical) || (is.numeric(yv) && length(unique(yv)) == 2L)) {
			cats <- if (is.factor(yv)) levels(yv) else sort(unique(yv))
			.bay_tell(st, paste0("Categorical response, ", length(cats), " categories: ", paste(cats, collapse = ", "), "."),
				paste0("R\u00e9ponse cat\u00e9gorielle, ", length(cats), " modalit\u00e9s : ", paste(cats, collapse = ", "), "."))
			.bay_tell(st, paste0("One question: do the proportions differ according to ", nm[1L], "? (counterpart of the chi-squared test)"),
				paste0("Une question : les proportions diff\u00e8rent-elles selon ", nm[1L], " ? (pendant du test du khi-deux)"), "decide")
		} else {
			.bay_tell(st, paste0("One question: does ", input$response_name, " differ according to ", nm[1L], "?"),
				paste0("Une question : ", input$response_name, " diff\u00e8re-t-il selon ", nm[1L], " ?"), "decide")
			.bay_tell(st, "No preliminary test chooses the method; the shapes of the noise are compared.",
				"Aucun test pr\u00e9alable ne choisit la m\u00e9thode ; les formes du bruit sont compar\u00e9es.", "note")
		}
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
			.bay_tell(st, paste0("One question: does ", input$response_name, " differ according to ", nm[1L], ", each subject compared with itself?"),
				paste0("Une question : ", input$response_name, " diff\u00e8re-t-il selon ", nm[1L], ", chaque sujet compar\u00e9 \u00e0 lui-m\u00eame ?"), "decide")
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
.bay_story_sensitivity <- function(st, sens, hyp, het, log_bf = NULL) {
	if (is.null(st)) return(invisible(NULL))
	if (!is.null(log_bf)) {
		stb <- .bay_stab_txt(sens$verdict)
		.bay_step_result(st, paste0("Bayes factor ", .bay_bf_txt(log_bf, "."), ", ", stb$en),
			paste0("facteur de Bayes ", .bay_bf_txt(log_bf, ","), ", ", stb$fr))
	}
	if (is.null(sens)) {
		.bay_tell(st, "Sensitivity to the priors not computed (sensitivity = FALSE).",
			"Sensibilit\u00e9 aux a priori non calcul\u00e9e (sensitivity = FALSE).")
		.bay_tell(st, "The conclusion is not tested against other starting points.",
			"La conclusion n'est pas \u00e9prouv\u00e9e contre d'autres points de d\u00e9part.", "decide")
		return(invisible(NULL))
	}
	r1 <- .bay_num_txt(min(sens$r)); r2 <- .bay_num_txt(max(sens$r))
	rg <- sens$range_lead
	if (rg[1L] > 0.999) {
		.bay_tell(st, paste0("Sensitivity: \"", hyp$label_en[1L], "\" stays above 99.9 %."),
			paste0("Sensibilit\u00e9 : \u00ab ", hyp$label_fr[1L], " \u00bb reste au-dessus de 99,9 %."))
	} else {
		.bay_tell(st, paste0("Sensitivity: \"", hyp$label_en[1L], "\" goes from ", .bay_pct(rg[1L], "."), " to ", .bay_pct(rg[2L], "."), "."),
			paste0("Sensibilit\u00e9 : \u00ab ", hyp$label_fr[1L], " \u00bb va de ", .bay_pct(rg[1L], ","), " \u00e0 ", .bay_pct(rg[2L], ","), "."))
	}
	.bay_tell(st, paste0("r goes from ", r1$en, " to ", r2$en, ", the starting chances from 1 against 2 to 2 against 1", if (het) ", the spread of the variances varies" else "", "."),
		paste0("r va de ", r1$fr, " \u00e0 ", r2$fr, ", les chances de d\u00e9part de 1 contre 2 \u00e0 2 contre 1", if (het) ", la r\u00e9partition des variances varie" else "", "."), "note")
	if (identical(sens$verdict, "changes")) {
		.bay_tell(st, "Caution added: the leading hypothesis can reverse; the bilan will not assert it.",
			"Prudence ajout\u00e9e : l'hypoth\u00e8se en t\u00eate peut s'inverser ; le bilan ne l'affirmera pas.", "decide")
	} else if (identical(sens$verdict, "crosses")) {
		.bay_tell(st, "The order does not change; only the degree of support can.",
			"L'ordre ne change pas ; seul le degr\u00e9 de soutien peut changer.", "decide")
	} else {
		.bay_tell(st, "The order does not change.", "L'ordre ne change pas.", "decide")
	}
	invisible(NULL)
}

# Step 6, one factor: which groups differ, and by how much.
.bay_story_posthoc <- function(st, ph, caution = FALSE, subjects = FALSE) {
	if (is.null(st) || is.null(ph)) return(invisible(NULL))
	cont <- isTRUE(ph$delta_unit %in% c("logor", "tvd"))
	if (cont) {
		.bay_tell(st, "Which groups differ, and by how much", "Quels groupes diff\u00e8rent, et de combien", "step",
			fun = c("Dirichlet-multinomial {KefiR}, hdi() {bayestestR}", "Dirichlet-multinomiale {KefiR}, hdi() {bayestestR}"))
	} else if (identical(ph$method, "paired")) {
		.bay_tell(st, "By how much the two levels differ", "De combien les deux niveaux diff\u00e8rent", "step",
			fun = "posterior(ttestBF()) {BayesFactor}, hdi() {bayestestR}")
	} else if (subjects) {
		.bay_tell(st, "Which conditions differ, and by how much", "Quelles conditions diff\u00e8rent, et de combien", "step",
			fun = "lmBF(), posterior() {BayesFactor}, hdi() {bayestestR}")
	} else {
		.bay_tell(st, "Which groups differ, and by how much", "Quels groupes diff\u00e8rent, et de combien", "step",
			fun = "lmBF(), posterior() {BayesFactor}, hdi() {bayestestR}")
	}
	pr <- ph$pairs
	dg <- if (!is.null(ph$dg)) ph$dg else .bay_dg(c(pr$diff_low, pr$diff_high))
	sep <- function(dc) if (dc == ",") " ; " else "; "
	gap <- function(i, dc) {
		if (cont) {
			f <- function(v) gsub(".", dc, formatC(100 * v, format = "f", digits = 1), fixed = TRUE)
			return(paste0(f(pr$diff[i]), " points [", f(pr$diff_low[i]), sep(dc), f(pr$diff_high[i]), "]"))
		}
		f <- function(v) .bay_fx(v, dg, dc)
		paste0(f(pr$diff[i]), " [", f(pr$diff_low[i]), sep(dc), f(pr$diff_high[i]), "]")
	}
	unit <- .bay_delta_unit(ph)
	dl <- function(i, dc) if (is.na(pr$delta[i]) || cont) "" else paste0(if (dc == ",") ", soit " else ", i.e. ",
		gsub(".", dc, formatC(pr$delta[i], format = "f", digits = 2), fixed = TRUE), " ", if (dc == ",") unit$fr else unit$en)
	big <- function(n) formatC(n, format = "d", big.mark = " ")
	if (ph$k > 2L && identical(ph$method, "gibbs")) {
		mc <- max(1L, ceiling(200 * ph$mc_error))
		.bay_tell(st, paste0("Change of strategy: ", big(ph$n_partitions), " groupings, too many to compute all; explored by draws, ", big(ph$sweeps), " rounds."),
			paste0("Changement de strat\u00e9gie : ", big(ph$n_partitions), " regroupements, trop pour tous les calculer ; explor\u00e9s par tirage, ", big(ph$sweeps), " tours."), "decide")
		.bay_tell(st, paste0("Sampling error of about ", mc, if (mc == 1L) " point" else " points", " on the probabilities."),
			paste0("Erreur de tirage d'environ ", mc, if (mc == 1L) " point" else " points", " sur les probabilit\u00e9s."), "note")
		.bay_tell(st, "Groupings explored by the Gibbs sampler of Gopalan and Berry (1998).",
			"Regroupements explor\u00e9s par l'\u00e9chantillonneur de Gibbs de Gopalan et Berry (1998).", "detail")
	}
	if (ph$k > 2L) {
		top <- .bay_top_grouping(ph)
		if (!is.null(top)) {
			.bay_step_result(st, paste0(top$label, " ahead (", .bay_pct(top$p, "."), ")"), paste0(top$label, " en t\u00eate (", .bay_pct(top$p, ","), ")"))
			pt <- ph$partitions[grepl("[<>]", ph$partitions$partition), , drop = FALSE]
			pt <- pt[seq_len(min(3L, nrow(pt))), , drop = FALSE]
			pc <- pt$posterior / sum(ph$partitions$posterior[grepl("[<>]", ph$partitions$partition)])
			.bay_tell(st, paste0("If a difference exists, leading groupings: ", paste0(pt$partition, " (", vapply(pc, .bay_pct, character(1), dec = "."), ")", collapse = "; "), "."),
				paste0("Si une diff\u00e9rence existe, regroupements en t\u00eate : ", paste0(pt$partition, " (", vapply(pc, .bay_pct, character(1), dec = ","), ")", collapse = " ; "), "."))
		}
	} else {
		.bay_step_result(st, paste0(pr$group2[1L], " - ", pr$group1[1L], " = ", gap(1L, ".")), paste0(pr$group2[1L], " - ", pr$group1[1L], " = ", gap(1L, ",")))
	}
	shown <- .bay_posthoc_shown(pr)
	for (i in shown) {
		if (ph$k > 2L) {
			.bay_tell(st, paste0(pr$group1[i], "/", pr$group2[i], ": different at ", .bay_pct(pr$p_differ[i], "."), "; ",
					pr$group2[i], " - ", pr$group1[i], " = ", gap(i, "."), dl(i, "."), "."),
				paste0(pr$group1[i], "/", pr$group2[i], " : diff\u00e9rents \u00e0 ", .bay_pct(pr$p_differ[i], ","), " ; ",
					pr$group2[i], " - ", pr$group1[i], " = ", gap(i, ","), dl(i, ","), "."))
		} else {
			.bay_tell(st, paste0("Estimated gap ", pr$group2[i], " - ", pr$group1[i], " = ", gap(i, "."), dl(i, "."), "."),
				paste0("\u00c9cart estim\u00e9 ", pr$group2[i], " - ", pr$group1[i], " = ", gap(i, ","), dl(i, ","), "."))
		}
	}
	if (length(shown) < nrow(pr)) {
		.bay_tell(st, paste0("The other ", nrow(pr) - length(shown), " pairs lean towards alike; all pairs in $posthoc$pairs."),
			paste0("Les ", nrow(pr) - length(shown), " autres paires penchent vers semblables ; toutes dans $posthoc$pairs."), "note")
	}
	if (cont) {
		.bay_tell(st, paste0(if (ph$delta_unit == "logor") paste0("gap of the share of ", ph$category) else "total gap of the proportions", ", in points; [ ; ] = 95 % credible interval."),
			paste0(if (ph$delta_unit == "logor") paste0("\u00e9cart de la part de ", ph$category) else "\u00e9cart total des proportions", ", en points ; [ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 %."), "note")
		.bay_tell(st, "Gaps: posterior medians under a uniform Dirichlet prior; the log odds ratio is in $posthoc$pairs$delta.",
			"\u00c9carts : m\u00e9dianes a posteriori sous une Dirichlet uniforme ; le log du rapport de cotes est dans $posthoc$pairs$delta.", "detail")
	} else {
		.bay_tell(st, "[ ; ] = 95 % credible interval (95 % chance that the gap lies there).",
			"[ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 % (95 % de chances d'y trouver l'\u00e9cart).", "note")
		.bay_tell(st, "the estimated gap is a little smaller than the observed one, by caution with small samples.",
			"l'\u00e9cart estim\u00e9 est un peu plus petit que l'\u00e9cart observ\u00e9, par prudence avec les petits \u00e9chantillons.", "note")
		.bay_tell(st, "Estimated gap: posterior median of the model \"all groups may differ\", pulled slightly towards zero by the prior, hence a little smaller than the observed gap.",
			"\u00c9cart estim\u00e9 : m\u00e9diane a posteriori du mod\u00e8le \u00ab tous les groupes peuvent diff\u00e9rer \u00bb, l\u00e9g\u00e8rement ramen\u00e9e vers z\u00e9ro par l'a priori, donc un peu plus petite que l'\u00e9cart observ\u00e9.", "detail")
	}
	if (ph$k > 2L) {
		.bay_tell(st, "all pairs are judged together, without a Tukey-like correction, and cannot contradict each other.",
			"toutes les paires sont jug\u00e9es ensemble, sans correction \u00e0 la Tukey, et ne peuvent pas se contredire.", "note")
	}
	if (isTRUE(ph$caution_overlap)) {
		.bay_tell(st, "These gaps compare the groups where some have few values of the covariate: with caution.",
			"Ces \u00e9carts comparent les groupes l\u00e0 o\u00f9 certains ont peu de valeurs de la covariable : avec r\u00e9serve.", "decide")
	}
	if (caution) {
		.bay_tell(st, "These pairs assume a normal noise with a common variance, which step 3 did not find: with caution.",
			"Ces paires supposent un bruit normal \u00e0 variance commune, ce que l'\u00e9tape 3 n'a pas trouv\u00e9 : avec r\u00e9serve.", "decide")
	}
	invisible(NULL)
}

# Step 6 when no pair is computed: the step is still shown (same skeleton for every plan).
.bay_story_no_posthoc <- function(st, why = c("factors", "off")) {
	why <- match.arg(why)
	if (why == "factors") {
		.bay_tell(st, "Which groups differ - not available with several factors", "Quels groupes diff\u00e8rent - non disponible \u00e0 plusieurs facteurs", "step")
	} else {
		.bay_tell(st, "Which groups differ - not asked (posthoc = FALSE)", "Quels groupes diff\u00e8rent - non demand\u00e9 (posthoc = FALSE)", "step")
	}
	invisible(NULL)
}

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
		logor = list(en = "in log odds ratio", fr = "en log du rapport de cotes"),
		tvd = list(en = "", fr = ""),
		pooled = list(en = "SD", fr = "\u00e9cart-type"),
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
		.bay_tell(st, paste0(if (n_dropped == 1L) "Subject removed, without a measure at both levels: " else paste0(n_dropped, " subjects removed, without a measure at both levels: "), lst, "."),
			paste0(if (n_dropped == 1L) "Sujet retir\u00e9, sans mesure aux deux niveaux : " else paste0(n_dropped, " sujets retir\u00e9s, sans mesure aux deux niveaux : "), lst, "."), "decide")
	}
	cells <- interaction(fac, drop = TRUE)
	wi <- intersect(input$within, names(fac))
	wc <- if (length(wi)) interaction(fac[wi], drop = TRUE) else cells
	full <- sum(tapply(wc, id, function(v) length(unique(v))) == nlevels(wc))
	.bay_tell(st, paste0("Subjects: ", ns, "; complete: ", full, "."), paste0("Sujets : ", ns, " ; complets : ", full, "."))
	if (full < ns) {
		.bay_tell(st, "Incomplete subjects are kept: the model uses all their measures.",
			"Les sujets incomplets sont gard\u00e9s : le mod\u00e8le utilise toutes leurs mesures.", "decide")
	}
	dev <- y - stats::ave(y, id) + mean(y)
	chk <- .bay_story_checks(st, dev, cells, "condition", "condition")
	.bay_tell(st, "far values sought once each subject's own level is removed.",
		"valeurs \u00e9loign\u00e9es cherch\u00e9es une fois retir\u00e9 le niveau propre de chaque sujet.", "note")
	.bay_tell(st, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")
	.bay_tell(st, "Shape of the noise", "Forme du bruit", "step")
	.bay_step_result(st, "normal, only", "normale, seule")
	if (paired) {
		.bay_tell(st, "Paired data: the difference within each subject is compared.",
			"Donn\u00e9es appari\u00e9es : on compare la diff\u00e9rence au sein de chaque sujet.")
	} else {
		.bay_tell(st, "Repeated measures: each subject keeps its own level.",
			"Mesures r\u00e9p\u00e9t\u00e9es : chaque sujet garde son niveau propre.")
		.bay_tell(st, "conditions assumed equally correlated within a subject (not tested; m.test() checks sphericity).",
			"conditions suppos\u00e9es \u00e9galement corr\u00e9l\u00e9es au sein d'un sujet (non test\u00e9 ; m.test() contr\u00f4le la sph\u00e9ricit\u00e9).", "note")
		.bay_tell(st, "Random subject effect: compound symmetry assumed, as in the uncorrected repeated-measures ANOVA (Mauchly's test in m.test()).",
			"Effet sujet al\u00e9atoire : sym\u00e9trie compos\u00e9e suppos\u00e9e, comme dans l'ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es non corrig\u00e9e (test de Mauchly dans m.test()).", "detail")
	}
	.bay_tell(st, "heavy tails and unequal variances are not validated with subjects yet.",
		"queues lourdes et variances in\u00e9gales ne sont pas encore valid\u00e9es avec des sujets.", "note")
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(st, paste0(far, if (far == 1L) " value lies" else " values lie", " beyond 3 interquartile ranges: read the answer with caution."),
			paste0(far, if (far == 1L) " valeur" else " valeurs", " au-del\u00e0 de 3 \u00e9carts interquartiles : lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(st, "No far value puts the normal law in doubt at first sight.",
			"Aucune valeur \u00e9loign\u00e9e ne met d'embl\u00e9e la loi normale en doute.", "decide")
	}
	invisible(chk)
}

.bay_story_ancova <- function(st, input, y) {
	fac <- input$factors
	cov <- input$covariates
	cells <- interaction(fac, drop = TRUE)
	res <- stats::residuals(stats::lm(y ~ ., data = data.frame(y = y, cov, check.names = FALSE)))
	chk <- .bay_story_checks(st, res + mean(y), cells)
	.bay_tell(st, "far values sought on what the covariates do not explain.",
		"valeurs \u00e9loign\u00e9es cherch\u00e9es sur ce que les covariables n'expliquent pas.", "note")
	ov <- .bay_cov_overlap(cov, cells)
	for (i in seq_len(nrow(ov))) {
		v <- ov$covariate[i]
		lo <- .bay_num_txt(ov$low[i]); hi <- .bay_num_txt(ov$high[i])
		if (ov$low[i] > ov$high[i]) {
			.bay_tell(st, paste0("Covariate ", v, ": the groups do not overlap at all."),
				paste0("Covariable ", v, " : les groupes ne se recouvrent pas du tout."))
		} else {
			.bay_tell(st, paste0("Covariate ", v, ": common range ", lo$en, " to ", hi$en, ", at least ", .bay_pct(ov$min_share[i], "."), " of each group."),
				paste0("Covariable ", v, " : \u00e9tendue commune ", lo$fr, " \u00e0 ", hi$fr, ", au moins ", .bay_pct(ov$min_share[i], ","), " de chaque groupe."))
		}
	}
	if (any(ov$poor)) {
		cg <- vapply(cov, function(v) tapply(v, cells, mean), numeric(nlevels(cells)))
		cg <- matrix(cg, nrow = nlevels(cells), dimnames = list(levels(cells), names(cov)))
		.bay_tell(st, paste0("Mean ", colnames(cg)[1L], " per group: ", paste0(rownames(cg), " ", format(signif(cg[, 1L], 3), trim = TRUE), collapse = ", "), ": group and covariate nearly confounded."),
			paste0(colnames(cg)[1L], " moyen par groupe : ", paste0(rownames(cg), " ", gsub(".", ",", format(signif(cg[, 1L], 3), trim = TRUE), fixed = TRUE), collapse = ", "), " : groupe et covariable presque confondus."))
		.bay_tell(st, "Caution added: the adjustment extrapolates; support capped at weak.",
			"Prudence ajout\u00e9e : l'ajustement extrapole ; soutien plafonn\u00e9 \u00e0 faible.", "decide")
	}
	ctr <- vapply(cov, mean, numeric(1))
	ct <- .bay_num_txt(ctr)
	.bay_tell(st, paste0("Covariate centred on its mean (", paste0(names(cov), " = ", ct$en, collapse = ", "), "): groups compared at this value."),
		paste0("Covariable centr\u00e9e sur sa moyenne (", paste0(names(cov), " = ", ct$fr, collapse = ", "), ") : groupes compar\u00e9s \u00e0 cette valeur."), "decide")
	.bay_tell(st, "Centring (Aiken and West 1991): the group effect is read at the mean of the covariate, inside the data, not at 0.",
		"Centrage (Aiken et West 1991) : l'effet du groupe se lit \u00e0 la moyenne de la covariable, dans les donn\u00e9es, et non \u00e0 0.", "detail")
	.bay_tell(st, "No value is removed.", "Aucune valeur retir\u00e9e.", "decide")
	.bay_tell(st, "Shape of the noise", "Forme du bruit", "step")
	.bay_step_result(st, "normal with a common variance, only", "normale \u00e0 variance commune, seule")
	.bay_tell(st, "The covariate acts in a straight line (linearity not tested).",
		"La covariable agit en ligne droite (lin\u00e9arit\u00e9 non test\u00e9e).")
	.bay_tell(st, "heavy tails, unequal variances and curvature are not validated with covariates yet.",
		"queues lourdes, variances in\u00e9gales et courbure ne sont pas encore valid\u00e9es avec des covariables.", "note")
	far <- chk$extremes[["far"]]
	if (far > 0L) {
		.bay_tell(st, paste0(far, if (far == 1L) " value lies" else " values lie", " beyond 3 interquartile ranges: read the answer with caution."),
			paste0(far, if (far == 1L) " valeur" else " valeurs", " au-del\u00e0 de 3 \u00e9carts interquartiles : lire la r\u00e9ponse avec r\u00e9serve."), "decide")
	} else {
		.bay_tell(st, "No far value puts the normal law in doubt at first sight.",
			"Aucune valeur \u00e9loign\u00e9e ne met d'embl\u00e9e la loi normale en doute.", "decide")
	}
	invisible(ov)
}

# Step 6 when the slopes differ (D16): gaps at the quartiles of the covariate, and where
# the groups differ.
.bay_story_slopes <- function(st, ph, xname, term, p_int) {
	if (is.null(st)) return(invisible(NULL))
	.bay_tell(st, paste0("Which groups differ, along ", xname), paste0("Quels groupes diff\u00e8rent, selon ", xname), "step",
		fun = "lmBF(y ~ g * x), posterior() {BayesFactor}, hdi() {bayestestR}")
	.bay_step_result(st, "the gap depends on the covariate", "l'\u00e9cart d\u00e9pend de la covariable")
	.bay_tell(st, paste0("Slopes that differ probable (", term, " ", .bay_pct(p_int, "."), ")."),
		paste0("Pentes diff\u00e9rentes probables (", term, " ", .bay_pct(p_int, ","), ")."))
	.bay_tell(st, paste0("Change of strategy: no single adjusted gap; gaps at the quartiles of ", xname, "."),
		paste0("Changement de strat\u00e9gie : pas d'\u00e9cart ajust\u00e9 unique ; \u00e9carts aux quartiles de ", xname, "."), "decide")
	for (L in .bay_slopes_lines(ph, xname)) .bay_tell(st, L$en, L$fr)
	.bay_tell(st, "[ ; ] = 95 % credible interval; the zone gathers the values where it excludes 0.",
		"[ ; ] = intervalle de cr\u00e9dibilit\u00e9 \u00e0 95 % ; la zone r\u00e9unit les valeurs o\u00f9 il exclut 0.", "note")
	.bay_tell(st, "Zone where the groups differ: Bayesian counterpart of the Johnson-Neyman region (Johnson and Neyman 1936; Preacher et al. 2006), on a grid of 101 values of the observed range.",
		"Zone o\u00f9 les groupes diff\u00e8rent : \u00e9quivalent bay\u00e9sien de la r\u00e9gion de Johnson-Neyman (Johnson et Neyman 1936 ; Preacher et al. 2006), sur une grille de 101 valeurs de l'\u00e9tendue observ\u00e9e.", "detail")
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
