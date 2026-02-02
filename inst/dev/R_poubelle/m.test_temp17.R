#' Automatic average/median/variance comparison function.
#'
#' @param x Quantitative response. Can be:
#'   \itemize{
#'     \item a numeric vector;
#'     \item a matrix/data frame (multivariate response for MANOVA);
#'     \item a character or integer vector of column names/indices if `data` is supplied.
#'   }
#'   Ignored when `formula` is used.
#' @param g Grouping factor (or vector). May also be a data frame of factors
#'   (several explanatory variables) when `formula` is not used. Ignored when
#'   `formula` is used.
#' @param data A data frame that contains the variables referenced by `x`, `g`,
#'   and/or `formula`. If omitted, symbols are looked up in the calling
#'   environment. Formulas like `data$Y ~ data$G` are also accepted provided
#'   the object `data` is visible.
#' @param formula A model formula such as `Y ~ G1 + G2` (or
#'   `cbind(Y1, Y2, ...) ~ ...` for MANOVA). Notation `Y ~ G | id` is internally
#'   converted to `Y ~ G + Error(id/...)` for repeated-measures designs.
#' @param paired Logical. Forces paired/repeated analysis. It is also enabled
#'   automatically when `id` is supplied, or in specific repeated-measures patterns.
#' @param id Subject identifier (a column name or index in `data`, or a vector
#'   the same length as `x`). Supplying `id` implies `paired = TRUE`.
#' @param wt Within-subject factor(s) for repeated measures. Can be a name,
#'   index, vector, or a data frame of factors present (or injected) in `data`.
#' @param within Alias of `wt`.
#' @param alpha Global p-value threshold (default `0.05`).
#' @param control Name of the category used as control in directed post-hoc
#'   procedures (e.g., Dunnett-type comparisons), when applicable.
#' @param verbose Logical. Print the step-by-step reasoning.
#' @param plot Logical. Draw distribution plots (boxplots, violins, etc.).
#' @param return Logical. If `TRUE`, returns a summary object (p-values, compact
#'   letter displays, etc.). If `FALSE`, returns raw p-values.
#' @param boot Logical. Enable bootstrap for means/medians where relevant.
#' @param iter Number of bootstrap iterations when `boot = TRUE`. If `0`,
#'   a default is chosen as `iter <- 1/alpha * 5`.
#' @param conf Confidence level for bootstrap intervals.
#' @param maxcat Maximum number of allowed groups; beyond this, some procedures
#'   may fail.
#' @param silent Logical. Suppress secondary warnings.
#' @param code Logical. Print a simplified R “recipe” that reproduces the analysis
#'   (disables `verbose`).
#' @param debug Logical. Emit detailed diagnostic messages (forces `code = FALSE`).
#'
#' @return m.test() runs a decision tree to choose the most appropriate test series for sample comparison.
#' @return She chooses the tests, justifies her choices.
#' @return It can output groups of means or a comparison to a control.
#' @return Finally, it will measure the robustness of the results by bootstrap.
#' @importFrom fda.usc fanova.hetero
#' @importFrom agricolae kurtosis
#' @importFrom agricolae skewness
#' @importFrom agricolae SNK.test
#' @importFrom lawstat levene.test
#' @importFrom WRS2 med1way
#' @importFrom WRS2 medpb
#' @importFrom WRS2 t1way
#' @importFrom WRS2 lincon
#' @importFrom FSA dunnTest
#' @importFrom onewaytests bf.test
#' @importFrom vioplot vioplot
#' @importFrom DescTools DunnettTest
#' @import methods
#' @export
#'
#' @examples
#' data(iris)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[1:100,1],iris[1:100,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,4],iris[,5],verbose=TRUE, plot=FALSE, return=FALSE, boot=FALSE)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE,control="virginica")
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,4],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
m.test <- function (x = NULL, g = NULL, data = NULL, formula = NULL,
			 paired = FALSE, id = NULL, wt=NULL, within=NULL,	 between = NULL,
			   alpha = 0.05, control = NULL, verbose = TRUE, plot = TRUE,
			   return = TRUE, boot = TRUE, iter = 0, conf = 0.95,
			   maxcat=50,silent=TRUE,
			   code=FALSE,debug=FALSE){
	call_expr <- match.call()
	k <- 0 # Compteur de messages
	if (code==TRUE) {verbose <- FALSE}
	if (debug==TRUE) {verbose <- FALSE ; code=FALSE}
	########################
	# 	iter
	########################
	if (iter==0) {iter <- 1/alpha*5}
	########################
	# 	Si data présente un problème
	########################
	if (!is.null(data) && !is.data.frame(data))
	  .exit("'data' must be a data.frame.",
	        "'data' doit être un data.frame.",
	        verbose = verbose, return = return)
	########################
	# 	Si formula n'est pas fournie, mais que x semble être une formula
	#   contrôle au passage de l'aspect multivarié.
	########################
	check_manova <- FALSE
	check_anova <- FALSE
	if (is.null(formula)) {
		if (inherits(x, "formula")) {
			.dbg("x appears to be a formula. Transferring to formula.",
			"`x` semble être une formule. Transfert vers `formula`.", debug=debug)
			formula <- x
			x <- NULL  # On supprime la valeur de `x` pour éviter des conflits
		}
	}
	call_expr <- match.call()
	if ("formula" %in% names(call_expr)) {
		formula <- call_expr$formula  # Capture l'expression brute AVANT évaluation
	}
	########################
	# 	 Contrôle au passage de l'aspect multivarié de x, si x fourni
	########################
	########################
	# Vérifier que x est multivarié - plusieurs variables dépendantes
	# A vérifier avant et après dans le code, but étant de passer vers anova <- TRUE (où un traitement manova est prévu dans .multi_factor_analysis.
	########################
	if (!is.null(x)) {
		if (is.data.frame(x) || is.matrix(x)) {
			num_cols <- sapply(x, is.numeric)
			if (!all(num_cols)) {
				.exit("'x' contains non-numeric columns.",
					"'x' présente des colonnes non numériques.",verbose = verbose, return = return)
			}
			if ((is.data.frame(g))&&(nrow(g)!=nrow(x))) {
				.exit("'x' and 'g' length differ.",
					"'x' et 'g' n'ont pas les mêmes dimensions.",verbose = verbose, return = return)
			}
			if (ncol(x)==1) {x <- x[,1] ; x <- as.vector(x)
				} else {check_manova <- TRUE} # x multivarié
		} else if (is.vector(x)&is.null(data)) {
			if (!is.numeric(x)) {
				.exit("'x' non-numeric.",
					"'x' non numérique.",verbose = verbose, return = return)
			} else if (is.vector(g)&&(length(g)!=length(x))) { # si g vecteur
				.exit("'x' and 'g' length differ.",
					"'x' et 'g' n'ont pas la même taille.",verbose = verbose, return = return)
				# comparer length de g
				# RAS
			} else if ((is.data.frame(g))&&(nrow(g)!=length(x))) { # si g data.frame
				.exit("'x' and 'g' length differ.",
					"'x' et 'g' n'ont pas les mêmes dimensions.",verbose = verbose, return = return)
			}
		} else if (!is.null(data)) {
		  # Cas 4a : x est numérique (on l'interprète comme une liste d'indices de colonnes dans data)
		  if (is.numeric(x)) {
			if (!all(x %in% seq_len(ncol(data)))) {
			  .exit("x contains invalid indices for data.",
				"x contient des indices invalides pour data.",verbose = verbose, return = return)
			}
			x <- data[, x, drop = FALSE]
			#message("x a été extrait de data à partir des indices fournis.")
		  } else if (is.character(x)) { # Cas 4b : x est de type caractère (on l'interprète comme une liste de noms de colonnes dans data)
			if (!all(x %in% colnames(data))) {
			  .exit("x contains column names that do not exist in data.",
				"x contient des noms de colonnes qui n'existent pas dans data.",verbose = verbose, return = return)
			}
			x <- data[, x, drop = FALSE]
			#message("x a été extrait de data à partir des noms de colonnes fournis.")
		  }  else {
			.exit("Unrecognized format of x. x must be a numeric vector, a data.frame, a matrix, or a list of indices/column names (if data is provided).",
			      "format de x non reconnu. x doit être un vecteur numérique, une data.frame, une matrice ou une liste d'indices/noms de colonnes (si data est fourni).",
				verbose = verbose, return = return)
		  }
		  if (ncol(x)==1) {x <- x[,1] ; x <- as.vector(x)
		  } else {check_manova <- TRUE} # x multivarié
	  } else {
		 .exit("Unrecognized format of x. x must be a numeric vector, a data.frame, a matrix, or a list of indices/column names (if data is provided).",
				verbose = verbose, return = return)
	  }
	}
	########################
	# 	Si data n'est pas fourni, mais que g semble être une data.frame
	########################
	.dbg(NULL,
	     "Si data n'est pas fourni, mais que g semble être une data.frame.",
	     debug=debug)
	if (is.null(data)) {
	  if (is.data.frame(g)) {
		# Vérifier si `g` contient au moins deux colonnes
		if (ncol(g) >= 2) {
		  # Vérifier si une colonne de `g` est identique à `x`
		  if (any(sapply(g, function(col) all(col == x)))) {
			.dbg("g contains a column identical to x. g is assigned to data without being combined with x.",
				"`g` contient une colonne identique à `x`. `g` est attribué à `data` sans combinaison avec `x`.",
				debug=debug)
			data <- g
		  } else {
			.dbg("g does not contain a column identical to x. g is combined with x and assigned to data.",
				"`g` ne contient pas de colonne identique à `x`. `g` est combiné à `x` et attribué à `data`.",
				debug=debug)
			# Combiner `x` et `g` dans un nouveau data.frame
			data <- cbind(x = x, g)
		  }
		  g <- NULL  # `g` est attribué à `data`, donc on le met à NULL
		} else {
		  .dbg("`g` is a data.frame with fewer than 2 columns. It is kept as a vector or a single column.",
			"`g` est un data.frame avec moins de 2 colonnes. Conservé en tant que vecteur ou unique colonne.",
			debug=debug)
		  # Si `g` est un data.frame avec une seule colonne, on le traite comme un vecteur
		  g <- as.vector(g[[1]]) ; g <- as.factor(g)
		}
	  } else if (is.atomic(g) && !is.null(g)) {
		.dbg("`g` is a vector. It is kept as is.","`g` est un vecteur. Conservé tel quel.", debug=debug)
		# Rien à faire, `g` reste inchangé
	  } else if (is.null(formula)) {
		# `g` n'est ni un vecteur ni un data.frame approprié
		.exit("g must be a vector or a data.frame if data is NULL.", "`g` doit être un vecteur ou un data.frame si `data` est NULL.",
			verbose=verbose,return=return)
	  }
	}
  if (debug) {
    print("g")
    print(head(g))
    if (is.character(id) && length(id) == 1L) {
      cat("id (name):", id, "\n")
    } else {
      print(head(id, 10))
    }
    print("formula")
    print(formula)
  }
	########################
	# Si formula est fournie, vérifier et garantir son format
	########################
	# balises à suppr une fois la fonction validée.
	if (!is.null(formula)) {
		#print(class(formula))
		#formula_expr <- eval(call_expr$formula, parent.frame())
	  # Ci-dessous, ligne 224, le passage de .formulator_save peut bousiller la formule...
	  formula <- if (!is.null(data)) .strip_data_dollar_safe(formula, data) else formula
	}
	########################
  # Si un identifiant est précisé (id pour mesures répétées) - key ididid
  # Vérifier que cet id est une colonne du jeu de données, sinon l'ajouter.
	# Sauvegarde du nom réel de `id` dans `data` (ou extraction si nécessaire)
	########################
	pair_vector <- c()
	.dbg("Checking for the presence of a pairing/repeated measures variable.",
		"Contrôle de la présence d'une variable d'appariement/répétitions.", debug=debug)
	if (!is.null(data) && !is.null(id)) {
		.dbg(NULL,"Tentative d'extraction de id de data si data fourni.", debug=debug)
	  if (is.character(id) && id %in% colnames(data)) {
		.dbg(NULL,"id est bien le nom d'une colonne de data.", debug=debug)
		# `id` est déjà le nom de la colonne
		#pair_name <- id
	  } else if (is.numeric(id) && id <= ncol(data) && length(id) == 1) {
		.dbg(NULL,"id est à considérer comme l'index d'une colonne de data.",debug=debug)
		# `id` est un index numérique, convertir en nom de colonne
		id <- colnames(data)[id]
	  } else if (length(id) == nrow(data.frame(data))) {
		if (debug==TRUE) {cat(.msg("\n","\tid est un vecteur externe à data dont il faut extraire nom et contenu.\n"))}
		# `id` est un vecteur externe, sauvegarder son contenu et son nom
		pair_vector <- id
		pair_expr <- substitute(id)  # Capturer l'expression originale
		id <- tryCatch(
		  deparse(pair_expr),  # Récupérer le nom ou une expression
		  error = function(e) "id"  # Utiliser un nom générique en cas d'erreur
		)
		if (is.null(id) || id == "") {
		  id <- "id"  # Défaut si aucun nom exploitable
		}
	  } else {
		if (verbose==TRUE) {
			warning(.msg("Warning: id is not corrected, and ignored.\n","Warning : 'id' doit être un nom ou un index valide d'une colonne de 'data'.\nid est ignoré.\n"))
		} else {return(1)}
	  }
	} else if (!is.null(id)) {
		  .dbg("Attempt to extract id while data is not provided.",
		       "Tentative d'extraction de id alors que data non fourni.", debug=debug)
		  # `id` est un vecteur mais `data` est absent
		  pair_vector <- id
		  pair_expr <- substitute(id)
		  id <- tryCatch(
			deparse(pair_expr),
			error = function(e) "id"
		  )
		  if (is.null(id) || id == "") {
			id <- "id"
		  }
		  #print("id");print(id)
	}
	.dbg(NULL,"Identifiant unique des individus (mesures répétées)",debug=debug)
	if (debug) {
	  if (is.character(id) && length(id) == 1L) {
	    cat("id (name):", id, "\n")
	  } else {
	    print(head(id, 10))
	  }
	}
	########################
	# Sauvegarde du nom réel de `wt` dans `data` (ou extraction si nécessaire)
	# Attention, si `within` est fourni, l'attribuer à `wt`
	# `wt` et `within` ne peuvent être fournis les 2, il s'agit de la même variable.
	########################
	########################
	# Advanced Handling of Within-Subject Variable (`wt`)
	########################
	wt_vector <- NULL
	wt_names <- NULL

	.dbg("Advanced within-subject variable processing", "Traitement avancé de la variable intra-sujet", debug=debug)

	# Handle 'within' parameter if provided
	if (!is.null(within)) {
	  wt <- within
	  .dbg("'within' parameter used for wt", "'within' utilisé pour wt", debug=debug)
	}

	if (!is.null(wt)) {
	  # Multiple handling strategies
	  if (!is.null(data)) {
		# Case 1: Multiple column names in data
		if (is.character(wt) && all(wt %in% colnames(data))) {
		  wt_names <- wt
		  .dbg("wt as multiple column names in data", "wt comme noms de colonnes multiples dans data", debug=debug)
		}
		# Case 2: Multiple column indexes
		else if (is.numeric(wt) && all(wt <= ncol(data))) {
		  wt_names <- colnames(data)[wt]
		  .dbg("wt as multiple column indexes", "wt comme index de colonnes multiples", debug=debug)
		}
		# Case 3: External vector matching data rows
		else if (length(wt) == nrow(data)) {
		  wt_vector <- wt
		  wt_expr <- substitute(wt)
		  wt_names <- tryCatch(
			deparse(wt_expr),
			error = function(e) "wt"
		  )

		  # Ajouter le vecteur à data
		  if (is.null(wt_names) || wt_names == "") wt_names <- "wt"
		  data[[wt_names]] <- wt_vector

		  .dbg("wt as external vector", "wt comme vecteur externe", debug=debug)
		}
		# Case 4: External dataframe or matrix
		else if (is.data.frame(wt) || is.matrix(wt)) {
		  .dbg("Complex wt: external dataframe/matrix", "wt complexe : dataframe/matrice externe", debug=debug)

		  # Validate cross-referencing with id
		  if (!is.null(id)) {


			id_values <- if (!is.null(data) && !is.null(id)) {
			  data[[id]]
			} else {
			  id
			}

			# Check completeness of data for each unique id
			unique_ids <- unique(id_values)
			complete_check <- sapply(unique_ids, function(unique_id) {
			  nrow(wt[wt[id] == unique_id, ]) > 0
			})

			if (!all(complete_check)) {
			  warning(.msg(
				"Warning: Some ids lack complete within-subject data.\n",
				"Attention : Certains identifiants manquent de données intra-sujet complètes.\n"
			  ))
			  return(1)
			}
		  }

		  # Extract column names
		  wt_names <- if (is.data.frame(wt)) colnames(wt) else colnames(wt)

		  # Ajouter les données à data si nécessaire
		  if (is.data.frame(wt)) {
			for (col in wt_names) {
			  data[[col]] <- wt[[col]]
			}
		  }
		}
		else {
		  warning(.msg(
			"Warning: 'wt' specification is invalid.\n",
			"Warning : La spécification de 'wt' est invalide.\n"
		  ))
		  return(1)
		}
	  }
	  # If data is null, handle external wt
	  else {
		wt_vector <- wt
		wt_expr <- substitute(wt)
		wt_names <- tryCatch(
		  deparse(wt_expr),
		  error = function(e) "wt"
		)

		# Créer data avec le vecteur
		data <- data.frame(wt_vector)
		names(data) <- if (is.null(wt_names) || wt_names == "") "wt" else wt_names
	  }
	}


	########################
	#	Si formula est précisée ainsi que data, il faut nettoyer data pour ne conserver que
	# ... les variables dépendantes (x simple ou multivarié) & les variables explicatives (g)
	# Exemple : si la formule est Y ~ I(log(X1))*X2, data = dt
	# et que dt fait en réalité 10 colonnes
	#	Il faut à la fin que le data soit nettoyé des colonnes inutiles (tout sauf X1, X2 et Y)
	# => Et que g corresponde à X1 et X2
	########################
  if (debug) {
    print("g")
    print(head(g))
    print("x")
    print(x[1:10])
    print("data")
    print(head(data))
    print("formula")
    print(formula)
    print("check_manova approx ligne 401")
    print(check_manova)
  }
  #* Note perso : de nombreux blocs ne sont pas titrés, encadrés, commentés.
	if (!is.null(formula)) {
	  check_multivariate_response <- function(formula, data = NULL, env = parent.frame()) {
	    stopifnot(inherits(formula, "formula"))
	    lhs <- formula[[2L]]

	    data_sym <- if (!is.null(data)) as.name(deparse(substitute(data))) else NULL

	    # --- PATCH ICI : éviter d’indexer un symbole ---
	    strip_dollar <- function(expr) {
	      if (!is.language(expr)) return(expr)
	      if (is.name(expr)) return(expr)  # <-- AJOUT essentiel

	      if (is.call(expr) && identical(expr[[1L]], as.name("$"))) {
	        if (!is.null(data_sym) && identical(expr[[2L]], data_sym)) return(expr[[3L]])
	        return(expr)
	      }
	      for (i in seq_along(expr)) expr[[i]] <- strip_dollar(expr[[i]])
	      expr
	    }
	    lhs <- strip_dollar(lhs)

	    if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) return(TRUE)

	    rho <- new.env(parent = env)
	    if (!is.null(data)) assign(deparse(substitute(data)), data, envir = rho)
	    val <- try(eval(lhs, envir = rho), silent = TRUE)

	    if (inherits(val, "try-error") && !is.null(data)) {
	      rho2 <- list2env(as.list(data), parent = env)
	      val <- try(eval(lhs, envir = rho2), silent = TRUE)
	    }

	    if (inherits(val, "try-error")) return(FALSE)
	    if (is.matrix(val) || is.data.frame(val)) return(ncol(val) > 1L)
	    FALSE
	  }

		# Usage
		check_manova <- check_multivariate_response(formula, data)
	}

	############################################
	# Mode formula ou mode classique
	############################################
	.dbg("Data extraction based on the presence of a formula.",
		"Extraction des données en fonction de la présence d'une formule.", debug=debug)

		########################
		# Si la formule est fournie et que l'on n'est pas dans une MANOVA
		#########################
	if (!check_manova && !is.null(formula) && inherits(formula, "formula")) {


	  formula_vars <- all.vars(formula)
		if (!is.null(data)){
		  # Vérifier et ajouter les variables manquantes à data
		  for (var_seule in formula_vars) {
			# Si la variable n'est pas dans data
			if (!var_seule %in% colnames(data)) {
			  # Chercher la variable dans l'environnement global
			  if (exists(var_seule, envir = .GlobalEnv)) {
				# Ajouter la variable à data
				data[[var_seule]] <- get(var_seule, envir = .GlobalEnv)
				.dbg(paste0("External variable added to data: ", var_seule), paste0("Variable externe ajoutée à data : ", var_seule), debug = debug)
			  } else {
				.exit(
				  paste0("Variable '", var_seule, "' not found in data or global environment."),
				  paste0("Variable '", var_seule, "' non trouvée dans data ou l'environnement global."),
				  return = return, verbose = verbose
				)
			  }
			}
		  }
		}
		# Gestion des formules avec `|`
		if (grepl("\\|", deparse(formula))) {
			#print("|")
			split_formula <- strsplit(deparse(formula), "\\|")[[1]]
			if (length(split_formula) == 2) {
				#print(split_formula)
				left_side <- strtrim(split_formula[1])
				right_side <- strtrim(split_formula[2])
				# Reconstruire une formule compatible
				formula <- as.formula(paste(left_side, "+ Error(", right_side, ")"))
				.dbg(paste0("\nFormula adjusted with 'Error()': ", deparse(formula)),
					paste0("\nFormule ajustée avec 'Error()' : ", deparse(formula)),
					debug=debug)
			} else {
				.exit("The formula contains a malformed '|'.",
					"La formule contient un '|' mal formé.",return=return,verbose=verbose)
			}
		}
		# Évaluer la formule avec l'environnement mis à jour
#print("On applique data et formula")
#print("formula") ; print(formula)
  if (debug) {
    print("===========================")
    print("TITRE A PERSONNALISER AVEC .dbg() selon la position dans le code")
    print("===========================")
  }



  # 6) extractions
  norm <- .normalize_from_formula(formula, data, id = id, wt_names = wt_names, na_strategy = "preserve")
  x  <- norm$x
  g  <- norm$g     # g reste data.frame (même à 1 colonne) ; tu factoriseras plus tard si tu pars en 1 facteur
  mf <- norm$mf

  #if (!is.factor(g)) g <- factor(g)


		########################
		# Si data est fournie
		#########################
		# Créer un environnement combinant le data.frame et l'environnement global
		if (!is.null(data) && is.data.frame(data)) {
		  env <- new.env(parent = parent.frame())
		  list2env(data, env)  # Ajouter les colonnes du data.frame à l'environnement
		  vars <- all.vars(formula)
		  # Ajouter les variables globales référencées dans la formule
			# Ne pas utiliser 'var' comme nom de variable
			for (variable in vars) {
			  if (!exists(variable, envir = env, inherits = FALSE)) {
				assign(variable, get(variable, envir = parent.frame(), inherits = TRUE), envir = env)
			  }
			}
		  # Associer l'environnement fusionné à la formule
		  environment(formula) <- env
		########################
		# Si data est fournie et n'est pas au bon format
		#########################
		} else if(!is.null(data)&&!is.data.frame(data)) {
			.exit("The 'data' parameter must be a valid data frame or NULL if the formula includes the data.",
				"Le paramètre 'data' doit être un data.frame valide ou NULL si la formule inclut les données.",
				verbose=verbose, return=return)
		}
		########################
		# Si id est fournie
		#########################
		# Modification de la formule pour inclure les variables d'identification et contextuelles
		# Initialiser la formule temporaire
		formula_temp <- formula

		# Gestion de la variable d'identification (id)
		if (!is.null(id)) {
		  paired <- TRUE  # Forcer `paired = TRUE` si `id` est spécifié

		  # Vérification de l'existence de la colonne id dans `data`
		  if (!id %in% colnames(data) && length(pair_vector) == 0) {
			.exit("The specified id column does not exist in data.",
				  "La colonne `id` spécifiée n'existe pas dans `data`.",
				  verbose = verbose, return = return)
		  }

		# Gestion des variables contextuelles (wt)
		if (!is.null(wt_names)) {
			###################
			#	Si wt spécifié
			###################
		  # Vérification de l'existence des colonnes spécifiées par `wt_names` dans `data`
		  for (wt_var in wt_names) {
			if (!wt_var %in% colnames(data)) {
			  .exit(paste0("The specified variable `", wt_var, "` does not exist in data."),
					paste0("La variable spécifiée `", wt_var, "` n'existe pas dans `data`."),
					verbose = verbose, return = return)
			}
		  }

		  # Ajouter chaque variable wt si absente de la formule
		  for (wt_var in wt_names) {
			if (!wt_var %in% all.vars(formula_temp)) {
			  formula_temp <- update(formula_temp, paste0(". ~ . + ", wt_var))
			}
		  }
		  # Construire la structure Error() selon le nombre de variables de répétition
		  if (length(wt_names) > 1) {
			interaction_terms <- paste(wt_names, collapse = " * ")
			formula_temp <- update(formula_temp, paste0(". ~ . + Error(", id, "/(", interaction_terms, "))"))
		  } else {
			formula_temp <- update(formula_temp, paste0(". ~ . + Error(", id, "/", wt_names[1], ")"))
		  }

		} else {
			###################
			#	Si wt non spécifié
			###################)
		  # Aucune variable de répétition spécifiée : gérer Error(id) ou Error(id/g)

		  # Vérifier le nombre de prédicteurs dans la formule
		  pred_vars <- all.vars(formula[[3]])  # Extraction des prédicteurs à droite de ~

		  if (length(pred_vars) == 1) {
			# Un seul prédicteur explicatif (g)
			formula_temp <- update(formula_temp, paste0(". ~ . + Error(", id, "/", pred_vars[1], ")"))
		  } else {
			# Plusieurs prédicteurs ou pas de prédicteur spécifique
			formula_temp <- update(formula_temp, paste0(". ~ . + Error(", id, ")"))
		  }
		}
		}
  if (debug) {
    print("===========================")
    print("TITRE A PERSONNALISER AVEC .dbg() selon la position dans le code")
    print("===========================")
    print("x")
    print(x[1:10])
    print("g")
    print(head(g))
    print("formula_temp")
    print(formula_temp)
    print("data")
    print(head(data))
  }
    ###################################
    # Extraire les données nécessaires à l'application de la formule
    ###################################
    #* Note perso : il faudra peut-être revoir ce code pour autoriser les NA dans le multifacteurs.
    # Conserver la formule complète (avec Error) pour les modèles ANOVA/aov ultérieurs
    f_with_error <- formula_temp

    # Version nettoyée pour construire model.frame()
    f_tmp <- if (!is.null(data)) .strip_data_dollar_safe(f_with_error, data) else f_with_error
    f_mf  <- .drop_error_term(f_tmp)
    if (debug) {
      print("f_tmp"); print(f_tmp)
      cat("Formule AVEC Error():  ", deparse(f_with_error), "\n")
      cat("Formule SANS Error():  ", deparse(f_mf), "\n")
    }
    # Évaluation robuste (on conserve tryCatch)
    mf <- tryCatch(
      model.frame(
        f_mf,
        data = if (!is.null(data)) data else parent.frame(),
        na.action = na.omit,
        drop.unused.levels = TRUE
      ),
      error = function(e) {
        .exit("Error evaluating the formula. Check your data and formula.",
              "Erreur lors de l'évaluation de la formule. Vérifiez vos données et votre formule.",
              return = return, verbose = verbose)
      }
    )
		###################
		#
		# Trier les données en fonction des variables explicatives et de l'appariement (par catégories croisées)
		#
		####################
		.dbg("Sorting data based on explanatory variables and pairing if provided.\n",
			 "Tri des données en fonction des variables explicatives et de l'appariement si fourni.",
			 debug = debug)
		# Récupérer les noms des colonnes explicatives, en excluant la colonne quantitative
		explicatives <- colnames(mf)[-1]  # Supposons que la première colonne est quantitative
		# Vérification de la présence de la variable d'appariement
		if (!is.null(id) && id %in% explicatives) {
		  # Déplacer la colonne d'appariement à la fin de la liste pour tri final
		  explicatives <- c(setdiff(explicatives, id), id)
		}
		# Tri des données par toutes les colonnes explicatives
		#mf <- mf[do.call(order, mf[explicatives]), ]
		# Tri robuste colonne par colonne (évite order() direct sur data.frame)
		ord <- do.call(order, lapply(mf[explicatives], function(col) {
		  if (is.factor(col)) as.integer(col)              # facteurs -> codes entiers
		  else if (is.character(col)) match(col, unique(col)) # chars -> ordre d’apparition
		  else col                                          # numériques -> OK direct
		}))
		mf <- mf[ord, , drop=FALSE]



		if (debug == TRUE) {
		  .dbg("Sorted data:", "Données triées :",debug=debug)
		  print(head(mf))
		}

		# Mettre à jour les données triées dans `data` et `cat`
		#x <- mf_sorted[, 1]  # Colonne quantitative
		#cat <- mf_sorted[, -1]  # Autres colonnes explicatives
		# Gestion des formules avec plusieurs variables explicatives
		##########################
		#	n variables explicatives
		##########################
		if (ncol(mf) > 2) {
		  # On a data et formula, pas besoin de x et g
		  .dbg("Multiple explanatory variables detected.",
				"Plusieurs variables explicatives détectées.", debug=debug)
		  #data <- mf
		  # Supprimer explicitement Error(...) avant model.frame
		  f_lin_multi <- .formulator_safe(formula_temp, mode = "linear", debug=debug)
		  f_lin_multi_noErr <- .drop_error_term(f_lin_multi)

		  # (debug utile)
		  if (debug) {
		    cat("formula BIS (lin): ", deparse(f_lin_multi), "\n")
		    cat("formula BIS (lin, no Error): ", deparse(f_lin_multi_noErr), "\n")
		  }


		  mf <- tryCatch(
		    model.frame(
		      f_lin_multi_noErr,
		      data = data,
		      na.action = na.omit,
		      drop.unused.levels = TRUE
		    ),
		    error = function(e) {
		      .exit("While evaluating the formula. Please check your data and formula.",
		            "Lors de l'évaluation de la formule. Vérifiez vos données et votre formule.",
		            verbose=verbose, return=return)
		    }
		  )
		  # formula <- as.formula(f_lin_multi_noErr)

		  # Ne pas écraser formula ! Garder f_with_error pour .multi_factor_analysis()
		  # La version sans Error() (f_lin_multi_noErr) est déjà utilisée pour model.frame() ci-dessus
		  formula <- f_with_error  # Conserver la formule avec Error()
		  # DEBUG CRITIQUE
		  if (debug) {
		    cat("\n=== DEBUG après affectation formula <- f_with_error ===\n")
		    cat("ncol(mf) = ", ncol(mf), "\n")
		    cat("paired = ", paired, "\n")
		    cat("formula après affectation = ", deparse(formula), "\n")
		    cat("f_with_error contient = ", deparse(f_with_error), "\n")
		    cat("=========================================================\n\n")
		  }
		  #data <- cbind(data,cat)
		  if ((ncol(mf)==3)&(paired==TRUE)) {
			.dbg("One-way ANOVA mode on paired data.",
				"Mode ANOVA à 1 facteur sur données appariées.",debug=debug)

			# Remplacer la partie problématique par :
			if (is.call(formula[[2]])) {
			  x <- mf[[1]]  # La transformation est déjà appliquée dans mf
			} else {
			  x <- mf[[1]]
			}
			# Un raisonnement sur x et g suffit
			# En théorie, on pourrait effacer data et formula ou définir une formula simple ici.
			x <- mf[[1]]  # Variable dépendante (quantitative)
			g <- data[,2]

		  } else {
			.dbg("Two-way or higher ANOVA mode activated.",
				"Mode ANOVA à 2 facteurs ou plus activé.",
					debug=debug)
			check_anova <- TRUE # plus de 2 colonnes, plus d'une variable explicative ! ==> ANOVA à 2 facteurs ou +
			# CORRECTION : Réassigner formula_temp à formula pour préserver Error()
			formula <- formula_temp
		  }

		##########################
		#	1 variable explicative
		##########################
		} else {
			# Réimporter mf avec les transformations de base (log et autre)
		  # --- SECOND BLOC (cas « 1 variable explicative ») ---
		  # Réimporter mf avec les transformations de base (log et autre)
		  f_tmp <- if (!is.null(data)) .strip_data_dollar_safe(formula_temp, data) else formula_temp
		  f_mf  <- .drop_error_term(f_tmp)

		  # Lineariser puis supprimer explicitement Error(...) AVANT model.frame
		  f_lin <- .formulator_safe(formula_temp, mode = "linear", debug=debug)
		  f_lin_noErr <- .drop_error_term(f_lin)

		  # (debug utile)
		  if (debug) {
		    cat("formula_temp (lin): ", deparse(f_lin), "\n")
		    cat("formula_temp (lin, no Error): ", deparse(f_lin_noErr), "\n")
		  }


		  mf <- tryCatch(
		    model.frame(
		      f_lin_noErr,
		      data = data,
		      na.action = na.omit,
		      drop.unused.levels = TRUE
		    ),
		    error = function(e) {
		      .exit("While evaluating the formula. Please check your data and formula.",
		            "Lors de l'évaluation de la formule. Vérifiez vos données et votre formule.",
		            verbose=verbose, return=return)
		    }
		  )
			# Remplacer la partie problématique par :
			if (is.call(formula[[2]])) {
			  x <- mf[[1]]  # La transformation est déjà appliquée dans mf
			} else {
			  x <- mf[[1]]
			}
			# Un raisonnement sur x et g suffit
			# En théorie, on pourrait effacer data et formula ou définir une formula simple ici.
			#x <- mf[[1]]  # Variable dépendante (quantitative)
			#g <- mf[-1]    # Variable explicative (qualitatives)
			#g <- as.vector(unlist(g))
		  # Ce code risque de bloquer sur l'introduction de between ? Peut-être...
		  expl <- colnames(mf)[-1]
		  to_drop <- c()
		  if (!is.null(id))        to_drop <- c(to_drop, id)
		  if (!is.null(wt_names))  to_drop <- c(to_drop, wt_names)
		  g_cols <- setdiff(expl, to_drop)

		  if (length(g_cols) != 1L) {
		    .exit("Ambiguous grouping: none or multiple grouping columns after excluding id/wt.",
		          "Groupement ambigu : aucune ou plusieurs colonnes de groupement après exclusion de id/wt.",
		          verbose=verbose, return=return)
		  }

		  g <- droplevels(factor(mf[[ g_cols[1] ]]))
		  # CORRECTION : Réassigner formula_temp à formula pour préserver Error()
		  formula <- formula_temp
		}
	} else if (!check_manova) {
	  # DEBUG CRITIQUE
	  if (debug) {
	    cat("\n=== DEBUG: Entrée dans bloc !check_manova ===\n")
	    cat("check_manova = ", check_manova, "\n")
	    cat("formula avant modifications = ", deparse(formula), "\n")
	    cat("===========================================\n\n")
	  }
		########################
		#	Si pas de formule fournie
		########################
		.dbg("Data analysis without a specified formula.",
			"Analyse des données sans formule spécifiée.", debug=debug)

		# Si pas de formule, utiliser les arguments explicites 'x' et 'g' sous la forme de simples  vecteurs
		if (is.null(x) || is.null(g)) {
			.exit("If no 'formula' is provided, the arguments 'x' and 'g' must be specified.","Si aucun 'formula' n'est fourni, les arguments 'x' et 'g' doivent être spécifiés.",verbose=verbose, return=return)
		}
		if (!is.numeric(x)) {
			.exit("'x' must be a numeric vector representing the quantitative variable.","'x' doit être un vecteur numérique représentant la variable quantitative.",verbose=verbose, return=return)
		}
		if (!is.factor(g)) {
		  g <- as.factor(g)  # Conversion en facteur si nécessaire
		}
		if (nrow(data.frame(x)) != nrow(data.frame(g))) {
			.exit("x and g do not have the same dimensions.","x et g ne présentent pas les mêmes dimensions.",verbose=verbose, return=return)
		}
		#
		if (paired == TRUE) {
			########################
			#	Si pas de formule fournie et données appariées
			########################
			if (is.null(id)) {
			#	paired est TRUE, il faut obligatoirement que id soit fourni
			#			(sauf si cat ne contient que 2 catégories) - sinon message d'erreur.
				# Aucun `id` n'est spécifié
				if (length(unique(g)) == 2) {
					# Deux catégories : autoriser un test apparié direct
					if (verbose == TRUE) {cat("Mode apparié activé pour deux catégories.\n")}
					# Trier selon les catégories et maintenir la correspondance
					combined_data <- data.frame(x, g)

					combined_data <- combined_data[order(combined_data$g), ]
					x <- combined_data$x
					g <- combined_data$g
					data <- combined_data
					formula <- formula(x~g)
				} else {
					# Plus de deux catégories, `id` obligatoire
					.exit("paired=TRUE requires an identification variable id if g contains more than two categories.\nOtherwise, use a formula to specify the consideration of identifiers via Error().",
						"`paired=TRUE` nécessite une variable d'identification des individidus`id`si `g` contient plus de deux catégories.\nSinon utiliser formula pour indiquer la prise en compte de identifiants via Error().",
						verbose=verbose, return=return)
				}
			} else {
				# `id` est spécifié, il doit être ajouté à un tableau global si plus de 2 catégories.
				# ou sinon servir à trier data et g si 2 catégories
				if (!is.null(data) && is.data.frame(data)) {
					# Vérifier si `id` est valide
					if (is.character(id) && id %in% colnames(data)) {
						# Nom de colonne
						id <- data[[id]]
					} else if (is.numeric(id) && id <= ncol(data) && length(id) == 1) {
						# Position de colonne
						id <- data[[id]]
					} else if (length(id) == nrow(data.frame(data))) {
						# Vecteur externe valide
						id <- id
					} else {
						.exit("id must be a valid name, column position, or vector.",
						      "`id` doit être un nom, une position de colonne, ou un vecteur valide.")
					}
				} else {
					# `data` non spécifié, `id` doit être un vecteur externe
					if (length(id) != length(g)) {
						.exit("id must be the same length as g.",
							"`id` doit avoir la même longueur que `g`.", verbose=verbose, return=return)
					}
				}
				#
				if (!is.null(data) && is.data.frame(data)) {
				  # Vérifier si wt est une seule variable (nom, position ou vecteur)
				  if (is.character(wt) && all(wt %in% colnames(data))) {
					# Un ou plusieurs noms de colonnes
					wt <- data[ , wt, drop = FALSE]
				  } else if (is.numeric(wt) && all(wt <= ncol(data))) {
					# Une ou plusieurs positions de colonnes
					wt <- data[ , wt, drop = FALSE]
				  } else if (is.data.frame(wt) || is.matrix(wt)) {
					# wt déjà fourni sous forme de data.frame ou de matrice
					if (nrow(wt) != nrow(data)) {
					  .exit("wt as a data.frame or matrix must have the same number of rows as data.",
							"`wt` sous forme de data.frame ou de matrice doit avoir le même nombre de lignes que `data`.",
							verbose=verbose, return=return)
					}
				  } else if (length(wt) == nrow(data)) {
					# wt est un vecteur externe valide
					wt <- data.frame(wt)
				  } else {
					.exit("wt must be a valid name, column position, vector, data.frame, or matrix.",
						  "`wt` doit être un nom, une position de colonne, un vecteur, un data.frame ou une matrice valide.")
				  }
				} else {
				  # `data` non spécifié, `wt` doit être un vecteur externe ou une structure cohérente
				  if (is.vector(wt) && length(wt) == length(g)) {
					wt <- data.frame(wt)
				  } else if ((is.data.frame(wt) || is.matrix(wt)) && nrow(wt) == length(g)) {
					wt <- wt
				  } else {
					.exit("wt must be a vector or structure with the same length as g.",
						  "`wt` doit être un vecteur ou une structure avec la même longueur que `g`.")
				  }
				}

				# Trier les données selon `id`

				########
				# A rectifier en indiquant si on a un jeu global (>2 catégories) ou si comme ici on reste à date et cat...
				########




#		CE CODE NE SEMBLE PAS GERER UNE G MULTIFACTORIELLE
######

				combined_data <- data.frame(x, g, id)
				# Ajout des variables de répétition (si multivariées, ajout de toutes les colonnes)
				if (!is.null(wt)) {
				  if (is.vector(wt)) {
					combined_data$wt <- wt  # Ajout direct si wt est un vecteur simple
				  } else if (is.data.frame(wt) || is.matrix(wt)) {
					combined_data <- cbind(combined_data, wt)  # Ajout de plusieurs colonnes si wt est multivarié
				  }
				}

				# Tri des données selon id et g
				#combined_data <- combined_data[order(combined_data$id, combined_data$g), ]
				ord <- if (.is_df_g(combined_data$g)) {
				  # Tri par id puis par chaque colonne de g (convertie en codes)
				  do.call(order, c(list(combined_data$id),
				                   lapply(combined_data$g, function(col) as.integer(factor(col)))))
				} else {
				  order(combined_data$id, combined_data$g)
				}
				combined_data <- combined_data[ord, , drop = FALSE]

				# Réassignation des variables triées
				x <- combined_data$x
				g <- combined_data$g
				id <- combined_data$id

				# Construction de la formule selon la nature de wt
				if (is.null(wt)) {
				  # Pas de variable de répétition
				  data <- combined_data
				  formula <- formula(x ~ g + Error(id/g))
				} else {
				  # Variable(s) de répétition présente(s)
				  if (is.vector(wt)) {
					# wt est un vecteur simple
					data <- combined_data
					formula <- formula(x ~ g + wt + Error(id/wt))
				  } else if (is.data.frame(wt) || is.matrix(wt)) {
					# wt est multivarié, combiner les colonnes sous une interaction (wt1 * wt2 * ...)
					wt_vars <- paste0(names(combined_data)[-(1:3)], collapse = " * ")
					formula <- as.formula(paste("x ~ g +", wt_vars, "+ Error(id/(", wt_vars, "))"))
					data <- combined_data

				  }
				}
			}
		} else {
			# Pas d'appariements des données, la formule de base est
			formula <- formula(x~g)
		}
	}
    if (paired==TRUE) {.dbg("Paired analysis activated.",
		                        "Analyse appariée activée.",debug=debug)}
    ###########################################
    #	Fonctions internes
    ###########################################
    .dbg("Compilation of base functions.","Compilation des fonctions de base.", debug=debug)
  	# -------------------------------------------------------------------
  	# Helpers robustes pour traiter g (facteur simple OU data.frame)
  	# -------------------------------------------------------------------

  	# g est-il un data.frame (cas multi-facteurs) ?
  	.is_df_g <- function(g) is.data.frame(g)

  	# Convertit g en une clé d'interaction sûre (ex: F:G),
  	# en conservant les niveaux utiles seulement (drop=TRUE).
  	.grp_interact <- function(g, drop=TRUE) {
  	  if (.is_df_g(g)) interaction(g, drop=drop) else droplevels(factor(g))
  	}

  	# Nombre de niveaux (groupes) effectifs dans g (mono ou multi-facteurs)
  	.n_levels_g <- function(g) nlevels(.grp_interact(g))

  	# tapply(x, groupe, FUN, ...) qui marche même si g est data.frame
  	.tapply_g <- function(x, g, FUN, ...) {
  	  tapply(x, .grp_interact(g), FUN, ...)
  	}

  	# Test de NA dans g, qu’il soit vecteur ou data.frame
  	.anyNA_g <- function(g) {
  	  if (.is_df_g(g)) anyNA(as.data.frame(g)) else anyNA(g)
  	}




	 ###################
	# paired doit être développé dans cette fonction : mis en argument mais tout à faire.

	#

  if (anyNA(x) || .anyNA_g(g)) {
  	if (verbose==TRUE) {cat("Warning! Missing values.\n")}
  	temp <- data.frame(x,g)
  	temp <- na.omit(temp)
  	x <- temp[,1]
  	g <- temp[,-1]
  }
  .dbg("Analysis of the presence of infinite values.",
     "Analyse de la présence de valeurs infinies.",
     debug = debug)
  if(any(!is.finite(x))) {
	.exit("Infinite values in data. Analysis impossible.",
	      "Valeurs infinies dans les données. Analyse impossible.")
  }
	.dbg("Control of repetitions, number of values per category.",
		 "Contrôle des répétitions, du nombre de valeurs par catégories.",
		 debug = debug)
  ###################
  # Si les variables explicatives ne sont pas numériques et qu'il y a qu'une seule valeur par variables explicatives croisées,
  # alors, l'on va autoriser paired=TRUE et dans ce cas :
	# id reçoit le nom de la dernière colonne de data
	# message d'alerte pour prévenir que les catégories croisées donnent toute 1 valeur
	# et que la variable x est désignée automatiquement comme variable d'appariement.
	# invitation de l'utilisateur à utiliser les arguments paired et id pour paramétrer selon son goût
	# si il ne reste qu'une variable explicative en plus de id (soit 3 colonnes en tout dans data), ANOVA<-FALSE
	# message if(debug==TRUE)&ANOVA==FALSE {cat(.msg("","Retour vers un scénario d'ANOVA à 1 facteur sur données appariées."))}
	# message if(debug==TRUE)&ANOVA==FALSE {cat(.msg("","Retour vers un scénario d'ANOVA à 2 facteurs sur données appariées."))}
  ###################
	.dbg(NULL,"=======================================",debug=debug)
	.dbg("Analyzing explanatory variables to determine if paired=TRUE can be automatically enabled.",
		"Analyse des variables explicatives pour déterminer si paired=TRUE peut être activé automatiquement.",
			debug=debug)
	.dbg(NULL,"=======================================",debug=debug)
  if (debug) {
    print(head(data))
    print("x");print(x[1:10])
    print("g");print(head(g,3))
    print("id")
    if (is.character(id) && length(id) == 1L) {
      cat("id (name):", id, "\n")
    } else {
      print(head(id, 10))
    }
    print("formula - done 3");print(formula)
  }
# Vérifier si les variables explicatives sont non numériques
if (!all(sapply(g, is.numeric))) {
	# Calcul des combinaisons croisées des variables explicatives
	combinaisons <- table(g)

	if (all(combinaisons == 1)) { # Vérifier si toutes les combinaisons donnent une seule valeur
		# Activer `paired=TRUE`
		paired <- TRUE

		# Attribuer la dernière colonne comme variable d'appariement
		id <- colnames(data)[ncol(data)]

		# Messages d'alerte
		cat(.msg("","Attention : Les catégories croisées des variables explicatives donnent toutes une seule valeur.\n"))
		cat(.msg("",paste("La variable", id, "est désignée automatiquement comme variable d'appariement.\n")))
		cat(.msg("","Veuillez utiliser les arguments `paired` et `id` pour ajuster les paramètres selon votre goût.\n"))

		# Vérifier le nombre de colonnes restantes dans `data`
		if (ncol(data) == 3) {
			# Une seule variable explicative restante
			check_anova <- FALSE
			x <- data[,1]
			g <- data[,2]
			if (!check_anova) {
				.dbg(NULL, "Retour vers un scénario d'ANOVA à 1 facteur sur données appariées.", debug=debug)
			}
		} else if (ncol(data) > 3) {
			# Plus d'une variable explicative restante
			#check_anova <- TRUE

			if (check_anova) {
				.dbg("Returning to a two-factor ANOVA scenario on paired data.",
					"Retour vers un scénario d'ANOVA à 2 facteurs sur données appariées.",
					debug=debug)
			}
		}
	}
}

    # --- Cas 1 facteur (pas MANOVA, pas ANOVA à ≥2 facteurs) ---
    #if (check_anova == FALSE && check_manova == FALSE) {
    if (!isTRUE(paired) && check_anova == FALSE && check_manova == FALSE) {

      # Comptes d'effectifs par groupe (clé d'interaction si g est multi-facteurs)
      if (length(x) != length(.grp_interact(g))) {
        .exit("Length mismatch between x and g after preprocessing.",
              "Longueurs incohérentes entre x et g après prétraitements.",
              verbose=verbose, return=return)
      }
      n_per_grp <- .tapply_g(x, g, length)

      # Alerte "échantillons trop petits" : max < 3
      if (max(n_per_grp, na.rm = TRUE) < 3) {
        if (return == FALSE) {
          .exit("No enough values in the samples.", "Pas assez de valeur dans l'échantillon",
                verbose = verbose, return = return)
        } else {
          .vbse("Warning! No enough values in the samples.",
                "Attention ! Pas assez de valeurs dans les groupes.",
                verbose = verbose)
        }
      }

      # Groupes sans données
      if (any(is.na(n_per_grp))) {
        warning("Warning! Some levels do not have corresponding data.")
        # Si g est un facteur simple, on peut encore droplevels
        if (!.is_df_g(g)) {
          g <- factor(g)
          if (length(unique(g)) < 2) stop("Not enough levels presenting g.")
        }
      }

      # Supprimer les groupes avec n<3
      if (min(n_per_grp, na.rm = TRUE) < 3) {
        if (verbose == TRUE) warning("Warning! No enough values for some samples. The categories concerned are ignored.")
        ind_drop <- names(n_per_grp)[n_per_grp < 3]
        `%notin%` <- Negate(`%in%`)
        keep <- .grp_interact(g) %notin% ind_drop
        x <- x[keep]
        if (.is_df_g(g)) g <- g[keep, , drop = FALSE] else g <- droplevels(factor(g[keep]))
      }

      # Variances par groupe
      var_per_grp <- .tapply_g(x, g, var, na.rm = TRUE)

      # Variabilité nulle partout
      if (max(var_per_grp, na.rm = TRUE) == 0) {
        .exit("No variability in samples.", "Absence de variabilité dans les échantillons.",
              verbose = verbose, return = return)
      }

      # Supprimer les groupes de variance nulle
      if (min(var_per_grp, na.rm = TRUE) == 0) {
        if (verbose == TRUE) warning("Warning! Some samples do not vary. Non-variable categories are ignored.")
        ind_drop <- names(var_per_grp)[var_per_grp == 0]
        `%notin%` <- Negate(`%in%`)
        keep <- .grp_interact(g) %notin% ind_drop
        x <- x[keep]
        if (.is_df_g(g)) g <- g[keep, , drop = FALSE] else g <- droplevels(factor(g[keep]))
      }

      # Nombre réel de groupes (mono ou multi-facteurs)
      nlev <- .n_levels_g(g)
      if (nlev <= 1) {
        .exit("Only one category in 'g'.", "Une seule catégorie dans 'g'.",
              verbose = verbose, return = return)
      }
      if (nlev > maxcat) {
        .exit("Too much categories. See maxcat.", "Trop de catégories. Voir l'argument maxcat.",
              verbose = verbose, return = return)
      }
    } else if (check_anova==TRUE && check_manova==FALSE) {
	# Construire la table des combinaisons des catégories
	# Extraire les données de la formule
  if (debug) {
    print('.formulator_safe(formula,mode="linear")')
    print(.formulator_safe(formula,mode="linear"))
    print(head(data))
  }

	data_model <- model.frame(.formulator_safe(formula,mode="linear",debug=debug), data=data)

	# Identifier les variables explicatives (sans la réponse)
	f_lin2 <- .formulator_safe(formula, mode = "linear", debug=debug)
	f_lin2_noErr <- .drop_error_term(f_lin2)
	variables_explicatives <- attr(terms(f_lin2_noErr), "term.labels")


	# Vérifier si toutes les variables explicatives sont numériques
	non_numeriques <- variables_explicatives[!sapply(variables_explicatives, function(var) is.numeric(data_model[[var]]))]
	if (length(non_numeriques) > 0) {
	  # Ne traiter que les variables non numériques comme catégoriques
	  variables_categoriques <- non_numeriques

	  # Construire la table des combinaisons pour les variables catégoriques uniquement
	  combinaisons <- table(data_model[variables_categoriques])

	  # Vérification des combinaisons croisées
	  if (min(combinaisons) == 0) {
		.dbg("Some crossed categories of the categorical variables are not specified. Check your data.",
			"Certaines catégories croisées des variables catégoriques ne sont pas renseignées. Vérifiez vos données.", debug=debug)
	  } else {
		.dbg("The crossed categories of the categorical variables are correctly specified.","Les catégories croisées des variables catégoriques sont bien renseignées.", debug=debug)
	  }
	} else {
	  .dbg("All explanatory variables are numeric. No interaction check is required.",
			"Toutes les variables explicatives sont numériques. Aucun contrôle de croisements n'est nécessaire.", debug=debug)
	}
	# paired est par défaut en 'auto'
	# si paired == TRUE, un contrôle est effectué pour vérifier qu'il y a autant de valeurs pour chaque catégories de cat.
	# Analyser tous types de formules Conditions ~ Conditions | Individus ou Conditions + Error (Individus/Conditions) ... etc.
	# Analyser dès qu'il y a des croisements des facteurs avec table() toutes à 1 (tout croisement renseigné) dans un Individus que appariés (signaler message). Dans ce cas on passe à paired = TRUE
	# ==> un message averti qu'on passe en apparié, et le nom du vecteur d'individus (ex : "Individus") est attribué à id.
	# Enfin pour chaque catégorie de cat, data est trié dans le même ordre des individus (afin de pouvoir comparar id par id en post-hoc).
	# 2 catégories ! Autoriser le paired en Wilcoxon et Student
	# >2 catégories : => Passer au test de friedman ou une aov sinon de ce type aov(Conditions + Error (Individus/Conditions))
	# envisager selon le test une mise en forme appropriée de la formule donnée en entrée, ou de l'intégration de la variable en id.
	# Remarque : que faire quand variances non égales ? aov(), friadman, autres ?
	# Prévoir aussi le contrôle des assomptions 1)  Normalité des résidus, 2) Homogénéité des variances et des co-variances (puisque les observations ne sont plus indépendantes puisqu'on mesure plusieurs fois les individus).
	# Nemiyi en post-hoc à prévoir si friedman à ajouter à la fonction posthoc

  }
  if (debug) {
    #* A mettre au propre
    print("========================")
    print("   Ligne 1241")
    print("========================")
    print("formula") ; print(formula)
    print("head(data)") ; print(head(data))
    print("x") ; print(x[1:10])
    print("g") ; print(head(g))
  }



  # --- Flag "répété à 1 facteur avec >2 niveaux" comme multi-facteur (répété)
  #lev_g <- if (is.factor(g)) nlevels(g) else nlevels(droplevels(factor(g)))
  lev_g <- .n_levels_g(g)
  has_error <- (!is.null(formula) && grepl("Error\\(", deparse(formula)))
  if (isTRUE(paired) && (has_error || (!is.null(id) && lev_g > 2))) {
    check_anova <- TRUE
  }

  #############################################################
  #
  #
  #			ANOVA à n facteurs
  #
  #
  #############################################################
  # --- ROUTAGE PAIRED / NON PAIRED -------------------------------------------

  .dbg(NULL,"========================", debug=debug)
  .dbg(NULL,"   ANOVA à n facteurs", debug=debug)
  .dbg(NULL,"========================", debug=debug)
  route <- NULL
  if (isTRUE(paired)) {
    if (isTRUE(check_manova) || isTRUE(check_anova)) {
      route <- if (isTRUE(check_manova)) "MANOVA(repeated)" else "ANOVA multi-facteurs (repeated)"
      .vbse(paste0("Routing: ", route), paste0("Routage : ", route), verbose=verbose, k=k, cpt="off")->k
      if (isTRUE(check_manova)) {
        bilan <- .manova_analysis(x=x, g=g, formula=formula, data=data,
                                  alpha=alpha, paired=TRUE, id=id,
                                  k=k, code=code, debug=debug, verbose=verbose)
      } else {
        bilan <- .multi_factor_analysis(x=x, g=g, formula=formula, data=data,
                                        alpha=alpha, paired=TRUE, id=id,
                                        k=k, code=code, debug=debug, verbose=verbose)
      }
    } else {
      lev <- nlevels(droplevels(factor(g)))
      if (lev > 2 && is.null(id)) {
        .exit("For paired designs with >2 levels, 'id' is required.",
              "Pour un plan apparié avec >2 modalités, 'id' est requis.",
              verbose=verbose, return=return)
      }
      route <- if (lev==2) "Paired (2 levels)" else "Paired (k>2 levels)"
      .vbse(paste0("Routing: ", route), paste0("Routage : ", route), verbose=verbose, k=k, cpt="off")->k
      bilan <- .one_factor_analysis(x=x, g=g, formula=formula, data=data,
                                    alpha=alpha, paired=TRUE, id=id,
                                    k=k, code=code, debug=debug, verbose=verbose)
    }
  } else {
    if (isTRUE(check_manova)) {
      .vbse("Routing: MANOVA", "Routage : MANOVA", verbose=verbose, k=k, cpt="off")->k
      bilan <- .manova_analysis(x=x, g=g, formula=formula, data=data,
                                alpha=alpha, paired=FALSE, id=id,
                                k=k, code=code, debug=debug, verbose=verbose)
    } else if (isTRUE(check_anova)) {
      .vbse("Routing: multi-factor ANOVA", "Routage : ANOVA multi-facteurs", verbose=verbose, k=k, cpt="off")->k
      bilan <- .multi_factor_analysis(x=x, g=g, formula=formula, data=data,
                                      alpha=alpha, paired=FALSE, id=id,
                                      k=k, code=code, debug=debug, verbose=verbose)
    } else {
      .vbse("Routing: one-factor (unpaired)", "Routage : 1 facteur (non apparié)", verbose=verbose, k=k, cpt="off")->k
      bilan <- .one_factor_analysis(x=x, g=g, formula=formula, data=data,
                                    alpha=alpha, paired=FALSE, id=id,
                                    k=k, code=code, debug=debug, verbose=verbose)
    }
  }
  # --- FIN ROUTAGE ------------------------------------------------------------

  # ---
	################################################################
	#		bilan
	################################################################
	x <- bilan[[1]]
	g <- bilan[[2]]
	check_normality <- bilan[[3]]
	check_variance_equal <- bilan[[4]]
	k <- bilan$k
	#################
	# post-hoc and bootstrap
	#################

	  if ((plot==TRUE)&(length(unique(g))<maxcat)) {
		.dbg(NULL,"Représentations graphiques.",debug=debug)
		boxplot(x~g,col="cyan")
		vioplot(x~g,col="#00AA0077",add=TRUE)
		stripchart(x~g,col="#FF000088",pch=16,vertical=TRUE,add=T,method="jitter",jitter=1/(length(unique(g))+2))
	  }
	.dbg(NULL,"Tests posts-hocs.",debug=debug)
	if (return==TRUE) {
		synth <- .posthoc(x,g,alpha=alpha,normal=check_normality,
		var.equal=check_variance_equal,
		control=control, code=code,debug=debug,verbose=verbose, paired=paired, boot=boot,iter=iter, conf=conf, k=k)
		return(synth)
	} else {return(pvals)}
}
