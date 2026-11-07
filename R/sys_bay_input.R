# sys_bay_input.R
# Normalises the inputs of m.bay() into one description of the plan. Same entry
# forms as m.test(), by reusing its internal helpers without changing them
# (decision D9, option B, PROPOSED: dev/mbay/DECISIONS.md). The forms that fail in
# m.test() (formula = named, g as a data.frame without data) work here on purpose:
# dev/mbay/PERIMETRE.md, 2.3.
#
# Result (class "bay_input"):
#   response   numeric vector, or matrix for MANOVA
#   response_name, response_names
#   factors    data.frame of the factors, one column per factor
#   covariates data.frame of the numeric covariates (possibly 0 columns)
#   id         factor of the subjects, or NULL
#   terms      fixed terms, spelled as in the formula expansion ("A", "B", "A:B")
#   slopes     predictors with a random slope on id (character, possibly empty)
#   within, between  factor names
#   plan       "oneway", "factorial", "ancova", "paired", "repeated", "mixed",
#              "manova"
#   discrete   TRUE if the response looks discrete (discret.test)
#   label      text of the model for the announcements
#   n_removed  rows dropped because of missing values


# ---------------------------------------------------------------------------
# Formula helpers
# ---------------------------------------------------------------------------

# Splits a right-hand side into its additive pieces, keeping calls such as
# Error(id/T) and (1 | id) whole.
.bay_rhs_pieces <- function(rhs) {
	if (is.call(rhs) && identical(rhs[[1L]], as.name("+")) && length(rhs) == 3L) {
		return(c(.bay_rhs_pieces(rhs[[2L]]), .bay_rhs_pieces(rhs[[3L]])))
	}
	if (is.call(rhs) && identical(rhs[[1L]], as.name("(")) ) {
		return(list(rhs))
	}
	list(rhs)
}

.bay_refuse_id <- function(id) {
	.exit(paste0("Subject id must be one variable, not: ", id, "."),
		paste0("L'identifiant du sujet doit \u00eatre une seule variable, pas : ", id, "."))
}

# Reads what the random or stratum parts of a formula say. Returns the formula
# without them, with id, within and the random slopes.
.bay_split_formula <- function(f) {
	lhs <- f[[2L]]
	rhs <- f[[3L]]
	id <- NULL
	within <- character(0)
	slopes <- character(0)
	# y ~ F | id
	if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
		id <- deparse1(rhs[[3L]])
		rhs <- rhs[[2L]]
		within <- all.vars(rhs)
	}
	keep <- list()
	for (p in .bay_rhs_pieces(rhs)) {
		is_error <- is.call(p) && identical(p[[1L]], as.name("Error"))
		is_bar <- is.call(p) && identical(p[[1L]], as.name("(")) && is.call(p[[2L]]) &&
			identical(p[[2L]][[1L]], as.name("|"))
		if (is_error) {
			inner <- p[[2L]]
			if (is.call(inner) && identical(inner[[1L]], as.name("/"))) {
				id <- deparse1(inner[[2L]])
				if (!is.name(inner[[2L]])) .bay_refuse_id(id)
				within <- union(within, all.vars(inner[[3L]]))
			} else {
				vars <- all.vars(inner)
				id <- vars[1L]
				within <- union(within, vars[-1L])
			}
		} else if (is_bar) {
			bar <- p[[2L]]
			id <- deparse1(bar[[3L]])
			if (!is.name(bar[[3L]])) .bay_refuse_id(id)
			lhs_vars <- all.vars(bar[[2L]])
			slopes <- union(slopes, lhs_vars)
		} else {
			keep[[length(keep) + 1L]] <- p
		}
	}
	rhs_new <- if (length(keep) == 0L) 1 else Reduce(function(a, b) call("+", a, b), keep)
	fixed <- f
	fixed[[3L]] <- rhs_new
	list(fixed = fixed, lhs = lhs, id = id, within = within, slopes = slopes)
}

# Variables of a formula without data: looked up in its environment. Handles d$y ~ d$F.
.bay_data_from_env <- function(f, env) {
	txt <- deparse1(f)
	if (grepl("\\$", txt)) {
		nm <- unique(gsub("\\$.*$", "", regmatches(txt,
			gregexpr("[A-Za-z_.][A-Za-z0-9_.]*\\$", txt))[[1L]]))
		if (length(nm) == 1L && exists(nm, envir = env) &&
			is.data.frame(get(nm, envir = env))) {
			data <- get(nm, envir = env)
			f2 <- stats::as.formula(gsub(paste0(nm, "$"), "", txt, fixed = TRUE),
				env = environment(f))
			return(list(formula = f2, data = data))
		}
	}
	vars <- all.vars(f)
	missing_vars <- vars[!vapply(vars, exists, logical(1), envir = env)]
	if (length(missing_vars) > 0L) {
		.exit(paste0("Variable(s) not found: ", paste(missing_vars, collapse = ", "), "."),
			paste0("Variable(s) introuvable(s) : ", paste(missing_vars, collapse = ", "), "."))
	}
	vals <- lapply(vars, get, envir = env)
	names(vals) <- vars
	list(formula = f, data = as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE))
}

# Numeric column: covariate or factor? Same rule as m.test (m.test.R:793): a numeric
# variable with more than 50 distinct values, or more than half of its values
# distinct, is continuous. 0/1 stays a factor.
.bay_is_continuous <- function(v) {
	if (!is.numeric(v)) return(FALSE)
	u <- length(unique(v[!is.na(v)]))
	if (u <= 2L) return(FALSE)
	isTRUE(u > 50L) || isTRUE(u > 0.5 * sum(!is.na(v)))
}


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

.bay_parse_input <- function(x = NULL, g = NULL, data = NULL, formula = NULL,
                             paired = FALSE, id = NULL, wt = NULL, within = NULL,
                             between = NULL, maxcat = 50, env = parent.frame()) {
	if (!is.null(within) && is.null(wt)) wt <- within
	if (!(identical(paired, TRUE) || identical(paired, FALSE))) {
		.exit("paired must be TRUE or FALSE.", "paired doit valoir TRUE ou FALSE.")
	}
	if (!is.numeric(maxcat) || length(maxcat) != 1L || is.na(maxcat)) {
		.exit("maxcat must be one number.", "maxcat doit \u00eatre un seul nombre.")
	}
	if (!is.null(between) && !is.null(data) && !all(between %in% names(data))) {
		.exit(paste0("between not found in data: ", paste(setdiff(between, names(data)), collapse = ", "), "."),
			paste0("between introuvable dans les donn\u00e9es : ", paste(setdiff(between, names(data)), collapse = ", "), "."))
	}
	# id: a simple name, an index with data, or a vector.
	if (!is.null(id) && length(id) == 1L && (is.na(id) || (is.numeric(id) && is.null(data)))) {
		.exit("id is not valid.", "id n'est pas valide.")
	}
	# x as a formula; a data.frame in second position is data.
	if (is.null(formula) && inherits(x, "formula")) {
		formula <- x
		x <- NULL
	}
	if (inherits(g, "formula") && is.null(formula)) {
		formula <- g
		g <- NULL
	}
	if (!is.null(formula) && !inherits(formula, "formula")) {
		.exit("formula must be a formula.", "formula doit \u00eatre une formule.")
	}
	if (!is.null(formula) && is.data.frame(g) && is.null(data)) {
		data <- g
		g <- NULL
	}
	if (!is.null(data) && !is.data.frame(data)) {
		.exit("data must be a data.frame.", "data doit \u00eatre un data.frame.")
	}

	id_vec <- NULL
	if (!is.null(id) && !(is.character(id) && length(id) == 1L && !is.null(data) && id %in% names(data))) {
		if (is.numeric(id) && length(id) == 1L && !is.null(data)) {
			id <- names(data)[id]
		} else if (length(id) > 1L) {
			id_vec <- id
			id <- NULL
		} else if (is.character(id) && length(id) == 1L) {
			.exit(paste0("id not found in data: ", id, "."),
				paste0("id introuvable dans les donn\u00e9es : ", id, "."))
		}
	}
	if (is.character(wt) && !is.null(data) && !all(wt %in% names(data))) {
		.exit(paste0("wt not found in data: ", paste(setdiff(wt, names(data)), collapse = ", "), "."),
			paste0("wt introuvable dans les donn\u00e9es : ", paste(setdiff(wt, names(data)), collapse = ", "), "."))
	}

	slopes <- character(0)
	if (!is.null(formula)) {
		if (is.null(data)) {
			found <- .bay_data_from_env(formula, environment(formula))
			formula <- found$formula
			data <- found$data
		}
		stripped <- .strip_factor_calls(formula, data)
		formula <- stripped$formula
		data <- stripped$data
		parts <- .bay_split_formula(formula)
		missing_vars <- setdiff(all.vars(formula), names(data))
		missing_vars <- missing_vars[!vapply(missing_vars, exists, logical(1),
			envir = environment(formula))]
		if (length(missing_vars) > 0L) {
			.exit(paste0("Variable(s) not found in data: ", paste(missing_vars, collapse = ", "), "."),
				paste0("Variable(s) introuvable(s) dans les donn\u00e9es : ",
					paste(missing_vars, collapse = ", "), "."))
		}
		if (!is.null(parts$id)) id <- parts$id
		wt <- union(wt, parts$within)
		slopes <- parts$slopes
		# Response evaluated once (BP-033): log(y) is computed here and never again.
		resp <- eval(parts$lhs, data, environment(formula))
		resp_label <- deparse1(parts$lhs)
		fx <- stats::terms(parts$fixed, data = data)
		terms_fixed <- attr(fx, "term.labels")
		# A term must be a variable: I(x > 2) or ifelse(...) would otherwise be analysed
		# through its raw variables under a label that says something else.
		bad_terms <- terms_fixed[!vapply(strsplit(terms_fixed, ":", fixed = TRUE), function(v) {
			all(vapply(v, function(s) identical(make.names(s), s) || grepl("^`.*`$", s), logical(1)))
		}, logical(1))]
		if (length(bad_terms)) {
			.exit(paste0("Terms must be variables. Create a column for: ", paste(bad_terms, collapse = ", "), "."),
				paste0("Les termes doivent \u00eatre des variables. Cr\u00e9er une colonne pour : ",
					paste(bad_terms, collapse = ", "), "."))
		}
		# A variable found outside data joins data, or is refused.
		for (v in setdiff(all.vars(parts$fixed[[3L]]), names(data))) {
			val <- get(v, envir = environment(formula))
			if (!is.atomic(val) || length(val) != nrow(data)) {
				.exit(paste0("Variable not found in data, or of the wrong length: ", v, "."),
					paste0("Variable absente des donn\u00e9es, ou de mauvaise longueur : ", v, "."))
			}
			data[[v]] <- val
		}
		pred_names <- setdiff(all.vars(parts$fixed[[3L]]), character(0))
		if (is.numeric(pred_names) || length(pred_names) == 0L) pred_names <- character(0)
	} else {
		# x, g forms: vectors, names or indices.
		if (is.null(x)) {
			.exit("A response is required: x or a formula.",
				"Une r\u00e9ponse est requise : x ou une formule.")
		}
		if (is.null(g) && is.null(id_vec) && is.null(id) && is.null(wt)) {
			.exit("A grouping is required: g or a formula.",
				"Un groupement est requis : g ou une formule.")
		}
		if (!is.null(data) && (is.character(x) || (is.numeric(x) && length(x) <= ncol(data)))) {
			resp <- if (length(x) > 1L) as.matrix(data[, x, drop = FALSE]) else data[[x]]
			resp_label <- if (is.character(x)) {
				paste(x, collapse = ", ")
			} else {
				paste(names(data)[x], collapse = ", ")
			}
		} else {
			resp <- if (is.data.frame(x)) as.matrix(x) else x
			resp_label <- if (is.matrix(resp)) "Y" else "x"
		}
		if (!is.null(data) && !is.null(g) && (is.character(g) || is.numeric(g)) &&
			length(g) <= ncol(data) && length(g) != nrow(data)) {
			if (is.numeric(g)) g <- names(data)[g]
			pred_names <- g
		} else if (is.data.frame(g)) {
			pred_names <- names(g)
			if (nrow(g) != NROW(resp)) {
				.exit("x and g do not have the same number of rows.",
					"x et g n\u0027ont pas le m\u00eame nombre de lignes.")
			}
			data <- if (is.null(data)) g else cbind(data, g[setdiff(names(g), names(data))])
		} else if (!is.null(g)) {
			if (is.list(g) || is.matrix(g) || !is.atomic(g)) {
				.exit("g must be a vector, a factor or a data.frame.",
					"g doit \u00eatre un vecteur, un facteur ou un data.frame.")
			}
			if (length(g) != NROW(resp)) {
				.exit("x and g do not have the same length.",
					"x et g n'ont pas la m\u00eame longueur.")
			}
			if (is.numeric(g) && any(is.infinite(g))) {
				.exit("g contains infinite values.", "g contient des valeurs infinies.")
			}
			pred_names <- ".g"
			if (is.null(data)) data <- data.frame(row.names = seq_along(as.matrix(resp)[, 1L]))
			data[[".g"]] <- g
		} else {
			pred_names <- character(0)
		}
		if (is.null(data)) data <- data.frame(row.names = seq_len(NROW(resp)))
		terms_fixed <- pred_names
	}

	resp_is_matrix <- is.matrix(resp) && ncol(resp) > 1L
	if (is.matrix(resp) && ncol(resp) == 1L) resp <- drop(resp)
	if (!is.numeric(resp)) {
		.exit("The response must be numeric. For a binary or categorical response, code it numerically.",
			"La r\u00e9ponse doit \u00eatre num\u00e9rique. Pour une r\u00e9ponse binaire ou cat\u00e9gorielle, la coder en nombres.")
	}
	if (NROW(resp) != nrow(data)) {
		.exit("The response and the data do not have the same number of rows.",
			"La r\u00e9ponse et les donn\u00e9es n'ont pas le m\u00eame nombre de lignes.")
	}
	if (!is.null(id_vec)) {
		if (length(id_vec) != NROW(resp)) {
			.exit("id and the response do not have the same length.",
				"id et la r\u00e9ponse n'ont pas la m\u00eame longueur.")
		}
		data[[".id"]] <- id_vec
		id <- ".id"
	}
	if (is.null(id) && (isTRUE(paired) || length(wt) > 0L)) {
		# paired = TRUE, or a within factor, without id: as in m.test, the order of the
		# rows pairs the observations.
		f_name <- if (length(wt)) wt else pred_names
		if (length(f_name) != 1L) {
			.exit("paired = TRUE or wt without id needs exactly one factor.",
				"paired = TRUE ou wt sans id demande exactement un facteur.")
		}
		fac <- as.character(data[[f_name]])
		data[[".id"]] <- stats::ave(seq_along(fac), fac, FUN = seq_along)
		id <- ".id"
	}
	if (is.character(wt) && length(wt)) pred_names <- union(pred_names, wt)
	pred_names <- setdiff(pred_names, id)

	# Classify the predictors. A sole numeric predictor is a factor unless continuous,
	# as in m.test; beside a factor, a continuous one is a covariate.
	is_cont <- vapply(pred_names, function(p) .bay_is_continuous(data[[p]]), logical(1))
	cov_names <- pred_names[is_cont]
	fac_names <- pred_names[!is_cont]
	terms_all <- if (!is.null(formula)) terms_fixed else fac_names
	terms_all <- union(terms_all, wt)
	terms_all <- setdiff(terms_all, id)

	# Missing values: complete cases on what the model uses, announced by the caller.
	used <- unique(c(fac_names, cov_names, id))
	keep <- stats::complete.cases(resp) &
		(if (length(used)) stats::complete.cases(data[used]) else TRUE)
	n_removed <- sum(!keep)
	resp <- if (is.matrix(resp)) resp[keep, , drop = FALSE] else resp[keep]
	data <- data[keep, , drop = FALSE]

	factors <- data[fac_names]
	for (p in fac_names) factors[[p]] <- droplevels(factor(factors[[p]]))
	covariates <- data[cov_names]
	id_out <- if (!is.null(id)) droplevels(factor(data[[id]])) else NULL

	for (p in fac_names) {
		if (nlevels(factors[[p]]) > maxcat) {
			.exit(paste0("More than maxcat levels in ", p, "."),
				paste0("Plus de maxcat niveaux dans ", p, "."))
		}
		if (nlevels(factors[[p]]) < 2L) {
			.exit(paste0("Factor ", p, " has fewer than 2 levels."),
				paste0("Le facteur ", p, " a moins de 2 niveaux."))
		}
	}
	if (length(fac_names) == 0L && length(cov_names) == 0L) {
		.exit("No factor to test: at least one grouping variable is required.",
			"Aucun facteur \u00e0 tester : il faut au moins une variable de groupement.")
	}

	# Plan.
	replicated <- FALSE
	if (!is.null(id_out)) {
		cell <- interaction(c(list(id_out), as.list(factors)), drop = TRUE)
		replicated <- any(table(cell) > 1L)
	}
	within_names <- intersect(wt, fac_names)
	if (is.null(id_out)) {
		between_names <- fac_names
	} else {
		# Within a subject, a factor that varies is within, one that does not is between.
		vary <- vapply(fac_names, function(p) {
			any(tapply(as.integer(factors[[p]]), id_out, function(v) length(unique(v))) > 1L)
		}, logical(1))
		within_names <- union(within_names, fac_names[vary])
		between_names <- fac_names[!vary]
	}
	if (!is.null(between)) between_names <- union(intersect(between, fac_names), between_names)
	plan <- if (resp_is_matrix) {
		"manova"
	} else if (!is.null(id_out)) {
		only_one <- length(fac_names) == 1L && length(cov_names) == 0L && !replicated
		if (replicated || length(slopes)) {
			"mixed"
		} else if (only_one && nlevels(factors[[1L]]) == 2L) {
			"paired"
		} else {
			"repeated"
		}
	} else if (length(cov_names)) {
		"ancova"
	} else if (length(fac_names) >= 2L) {
		"factorial"
	} else {
		"oneway"
	}

	label_pred <- if (identical(pred_names, ".g")) "g" else paste(terms_all, collapse = " + ")
	y_for_disc <- if (is.matrix(resp)) resp[, 1L] else resp
	list_out <- list(
		response = resp, response_name = resp_label,
		response_names = if (is.matrix(resp)) colnames(resp) else resp_label,
		factors = factors, covariates = covariates, id = id_out,
		terms = terms_all, slopes = slopes,
		within = within_names, between = between_names, plan = plan,
		discrete = isTRUE(tryCatch(discret.test(y_for_disc), error = function(e) FALSE)),
		label = paste0(resp_label, " ~ ", label_pred),
		n_removed = n_removed
	)
	class(list_out) <- "bay_input"
	list_out
}
