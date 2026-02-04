#' Automation function for pairwise calculations.
#'
#' @param x numerical vector
#' @param g category vector
#' @param type Character string indicating the type of test to perform. Options are:
#'   - "mean" for pairwise.t.test(p.adjust.method = "holm"),
#'   - "median" for pairwise.wilcox.test(p.adjust.method = "BH"),
#'   - "ks" for ks.test(),
#'   - "lincon" for lincon() from the \pkg{WRS2} package,
#'   - "dunn" for dunnTest() from the \pkg{FSA} package,
#'   - "BM" for brunner.munzel.test() from the \pkg{lawstat} package (independent samples),
#'   - "neminye" for frdAllPairsNemenyiTest() from the \pkg{PMCMRplus} package,
#'   - "boot": Establishes pairwise differences through a bootstrap approach applied to a descriptor (e.g., mean, median). The p-value is derived from the proportion of bootstrap samples where the descriptor difference is consistently positive or negative. Currently, `paired = TRUE` is not implemented for this option.
#' @param alpha threshold value of p-value to establish the groups.
#' @param control name of the category that will be used as a control to establish differences with '*', '**' and '***'.
#' @param pool.sd switch to allow/disallow the use of a pooled SD.
#' @param silent for displaying or not warnings.
#' @param boot to activate the bootstrap on 'mean', 'median', 'lincon', 'BM' or 'dunn'.
#' @param boot_type Character string specifying the descriptor for `pairwise.boot`. Options are:
#'   - "meanbp" for the moving average,
#'   - "mean" for the mean,
#'   - "median" for the median,
#'   - "sd" for the standard deviation.
#' @param tr trimming parameter for lincon test (default: 0.2).
#' @param iter number of bootstrap iterations when boot=TRUE.
#' @param conf confidence level of bootstrap.
#' @param paired Logical. If TRUE, performs paired tests (e.g., paired t-tests or Wilcoxon signed-rank tests). Not compatible with 'ks', 'lincon', 'dunn', 'BM', or 'boot' types.
#' @param debug logical, if TRUE enables debug mode for troubleshooting.
#'
#' @return This function automates the work of the ks.test(), lincon() functions of \pkg{WRS2}, brunner.munzel.test() from \pkg{lawstat}, pairwise.t.test(), pairwise.wilcox.test(), dunnTest() from \pkg{FSA}, and frdAllPairsNemenyiTest() from \pkg{PMCMRplus}, and extracts groups of means or comparisons to a control with the catego() function.
#' @return It pre-sorts the means/medians to ensure that the groups are identified in ascending order.
#' @return It also identifies the robustness of these groups by establishing a bootstrap.
#' @importFrom WRS2 lincon
#' @importFrom PMCMRplus frdAllPairsNemenyiTest
#' @importFrom FSA dunnTest
#' @export
#'
#' @examples
#' data(iris)
#' pairwise(iris[,1], iris[,5], type = "mean") # t.test
#' pairwise(iris[,1], iris[,5], type = "median", alpha = 0.01, boot = TRUE) # wilcox
#' pairwise(iris[,1], iris[,5], type = "ks")
#' pairwise(iris[,1], iris[,5], type = "lincon")
#' pairwise(iris[,1], iris[,5], type = "dunn")
#' pairwise(iris[,1], iris[,5], type = "BM") # Brunner-Munzel
#' pairwise(iris[,1], iris[,5], type = "neminye")
pairwise <- function(x, g, type = "mean", alpha = 0.05, control = c(),
                    pool.sd = FALSE, silent = TRUE, boot = FALSE,
                    boot_type = mean, tr = 0.2, iter = 500, conf = 0.95,
                    paired = FALSE, debug = FALSE) {
  .dbg(NULL,"Exécution de pairwise().",debug=debug)
  g <- factor(g)
  init_order <- levels(g)
  which(levels(g)%in%control)-> ind_control
  # Vérification de compatibilité
  valid_types <- c("mean", "median", "ks", "lincon", "dunn","boot","neminye","BM")
  if (!type %in% valid_types) stop("Error! arg type must be one of: ", paste(valid_types, collapse = ", "))
  if (paired && !type %in% c("mean", "median","neminye")) stop("'paired' is only compatible with 'mean', 'median' et 'neminye' types.")
  # Gestion des contrôles et des catégories : permet de mettre la catégorie de contrôle en première.
  reorder_groups <- function(x, g, ind_control,stat_fn) {
    if (length(ind_control) == 1) {
		categories <- c(levels(g)[ind_control], levels(g)[-ind_control])
		g <- ordered(g, levels = categories)
    } else {
		stat_values <- by(x, g, stat_fn)
		categories <- levels(g)
		indices <- order(stat_values)
		g <- ordered(g, levels = categories[indices])
    }
  }
  		# A mettre dans les importations PMCMRplus::frdAllPairsNemenyiTest
		# Calculer automatiquement un  id pour chaque groupe
		NemenyTest <- function(x,g) {
			id <- rep(1:table(g)[1],time = length(unique(g)))
			result <- frdAllPairsNemenyiTest(x~g|id)
			return(result)
		}
  		lincon_to_pairwise <- function(x,g,tr=tr,KB=TRUE) {
			#t1 <- suppressWarnings(suppressMessages(lincon(x ~ g, tr = tr)))
			t1 <- lincon(x~g, tr=tr)

		#	t1 <- suppressWarnings(suppressMessages({
		#		sink(tempfile()) # Redirige la sortie vers un fichier temporaire
		#		on.exit(sink())  # S'assure de restaurer la sortie normale
		#		lincon(x ~ g, tr = tr)
		#	}))

			mymat <- matrix(rep(NA, (length(t1$fnames) - 1)^2),
			                ncol = length(t1$fnames) - 1,
			                nrow = length(t1$fnames) - 1)
			rownames(mymat) <- t1$fnames[2:length(t1$fnames)] ; colnames(mymat) <- t1$fnames[1:length(t1$fnames)-1]
			for (i in 1:nrow(t1$comp)) {
				mymat[(t1$comp[i,2]-1),t1$comp[i,1]] <- t1$comp[i,6]
			}
			result <- list() ; result$p.value <- mymat
			return(result)
		}
	  # Faire une fonction de sample bootstrap
	bootstrap <- function(x, g, mat,test_fn=pairwise.t.test,p.adjust.method="holm",
			alpha=alpha,control=control, paired=FALSE,debug=FALSE,pool.sd=FALSE, tr=0.2) {
		.dbg(NULL,"bootstrap() - Exécution de bootstrap. - Récupération du format de sortie des p-values.",debug=debug)
		mydim <- dim(mat)
		mymat <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
		# sample
		if (paired) {
			.dbg("pairwise() - bootstrap() - Bootstrap for paired data, collecting indices of pairs.","pairwise() - bootstrap() - Bootstrap pour données appariées, collecte des indices des paires.",debug=debug)


			# Bootstrap pour données appariées
			paired_indices <- split(seq_along(x), g) # Indices par groupe
			n_pairs <- min(sapply(paired_indices, length)) # Nombre minimal d'appariements
			for (i in seq_len(iter)) {
			  x_boot <- x
			  for (j in seq_along(paired_indices)) {
				sampled_indices <- sample(paired_indices[[j]], n_pairs, replace = TRUE)
				x_boot[paired_indices[[j]]] <- x[sampled_indices]
			  }
			  if (identical(test_fn, pairwise.t.test)) {
				temp <- try(test_fn(x_boot, g, paired = paired, p.adjust.method = p.adjust.method,pool.sd=pool.sd,exact=FALSE)$p.value,silent=silent)
			  } else if (identical(test_fn, NemenyTest)) {
				temp <- try(test_fn(x_boot, g))$p.value
			  } else {
				temp <- try(test_fn(x_boot, g, paired = paired, p.adjust.method = p.adjust.method,exact=FALSE)$p.value,silent=silent)
			  }
			  # Vérification de la taille de temp
			  #if (!is.null(temp) && is.matrix(temp) && dim(temp)[1] == mydim[1] && dim(temp)[2] == mydim[2]) {
				#mymat[i, , ] <- temp
			  #} else {
				#.dbg(NULL, paste("Erreur de taille à l'itération", i, ": temp n'a pas la bonne taille."), debug = debug)
			  #}
			  mymat[i, , ] <- temp
			}
		} else {
			.dbg(NULL,"bootstrap() - Bootstrap pour données non appariées. Tirage des indices.",debug=debug)
			# Bootstrap classique
			for (i in seq_len(iter)) {
			  x_boot <- x
			  # sampling par catégories
			  for (j in levels(g)) {
				x_boot[g == j] <- sample(x[g == j], replace = TRUE)
			  }
			  if (identical(test_fn, lincon_to_pairwise)) {
				  # Appel spécifique pour lincon_to_pairwise()
				  temp <- try(test_fn(x_boot, g, tr=tr), silent = silent)
			  } else if (identical(test_fn, pairwise.t.test)) {
				  # Appel générique pour d'autres fonctions
				  temp <- try(test_fn(x_boot, g, paired = paired, p.adjust.method = p.adjust.method, pool.sd = pool.sd), silent = silent)
			  } else {
				  # Appel générique pour d'autres fonctions
				  temp <- try(test_fn(x_boot, g, paired = paired, p.adjust.method = p.adjust.method), silent = silent)
			  }
			   # Vérification du résultat de temp
			   if (inherits(temp, "try-error")) {
				  .dbg(NULL, paste("Erreur capturée à l'itération", i, ":", temp), debug = debug)
				  next  # Passer à l'itération suivante si une erreur est capturée
			   } else if (is.null(temp)) {
				  .dbg(NULL, paste("Résultat NULL à l'itération", i), debug = debug)
				  next
			   } else {
				  # Extraction des p-values si le test retourne une liste avec $p.value
				  if (is.list(temp) && !is.null(temp$p.value)) {
					temp <- temp$p.value
				  }
				  mymat[i, , ] <- temp
			   }
			  #temp <- test_fn(x_boot, g, paired = paired, p.adjust.method = p.adjust.method)$p.value
			  #mymat[i, , ] <- temp
			}
		}
		.dbg("pairwise() - bootstrap() - Confidence interval estimation.","pairwise() - bootstrap() - Estimation des intervalles de confiance.",debug=debug)
		# Calcul des intervalles de confiance
		output <- list()
		apply(mymat, c(2, 3), quantile, probs = conf, na.rm = TRUE) -> output$p.value
		colnames(output$p.value) <- colnames(mat)
		rownames(output$p.value) <- rownames(mat)
		.dbg(NULL,"Matrice de 'p-values' correspondants aux différences.",debug=debug)
		output <- catego(output,alpha=alpha,control=control, debug=debug)
		return(output)
	}
	#----------------------------------------
	#     Moyennes
	#----------------------------------------
	if (type=="mean") {
		.dbg(NULL,"pairwise() - Bootstrap de pairwise.t.test().",debug=debug)
		g <- reorder_groups(x, g, ind_control,mean)
		pairwise.t.test(x,g,pool.sd=pool.sd,p.adjust.method="holm",paired=paired)-> result
		groups <- catego(result,control=control, debug=debug)
		if (boot==TRUE) {
				.dbg(NULL,"Bootstrap.",debug=debug)
				groups$bootstrap <- bootstrap(x,g,groups$p.value,pairwise.t.test,alpha=alpha,
					control=control,paired=paired,pool.sd=pool.sd)
		}
		.dbg(NULL,paste("Résultats du bootstrap.", groups$bootstrap),debug=debug)
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		.dbg(NULL,"pairwise() - Fin de pairwise().",debug=debug)
		return(groups)
	#----------------------------------------
	#     Médianes//Rangs
	#----------------------------------------
	} else if (type == "median" ){
		.dbg(NULL,"pairwise.wilcox.test()",debug=debug)
		p.w.t <- function(x,g, paired = FALSE, exact = FALSE, silent = TRUE,p.adjust.method="holm"){#,pool.sd=pool.sd,p.adjust.method="holm",paired=paired){
			result <- try(pairwise.wilcox.test(x,g,
						p.adjust.method=p.adjust.method,
						paired=paired,
						exact = exact),silent=silent)
			return(result)
		}
		g <- reorder_groups(x, g, ind_control,median)
		result <- p.w.t(x,g)
		groups <- catego(result,control=control, debug=debug)
		if (boot==TRUE) {
			groups$bootstrap <- bootstrap(x,g,groups$p.value,p.w.t,
					alpha=alpha,control=control,
					paired=paired, debug=debug)
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		.dbg(NULL,"pairwise() - Fin de pairwise().",debug=debug)
		return(groups)
		#----------------------------------------
		#     Brunner-Munzel (rangs solides)
		#----------------------------------------
	} else if (type == "BM" ){ # Brunner–Munzel test
	  .dbg(NULL,"lawstat::brunner.munzel.test",debug=debug)
	  # Fonction wrapper pour Brunner-Munzel en format pairwise
	  BM_to_pairwise <- function(x, g, p.adjust.method = "holm", paired = FALSE, silent = TRUE) {
	    # BM_to_pairwise ne gère pas paired=TRUE (test indépendant par nature)
	    # arg paired conservé pour normaliser l'appel en bootstrap()
	    if (is.ordered(g)) {
	      g <- factor(as.character(g), levels = levels(g))
	    }
	    unique_g <- levels(g)
	    n_groups <- length(unique_g)
	    # Créer une matrice triangulaire inférieure pour les p-values
	    mymat <- matrix(NA, nrow = n_groups - 1, ncol = n_groups - 1)
	    rownames(mymat) <- unique_g[-1]
	    colnames(mymat) <- unique_g[-n_groups]
	    # Stocker toutes les p-values brutes pour ajustement
	    all_pvals <- c()
	    comparisons <- list()
	    # Effectuer tous les tests pairwise
	    for (i in 1:(n_groups - 1)) {
	      for (j in (i + 1):n_groups) {
	        x1 <- x[g == unique_g[i]]
	        x2 <- x[g == unique_g[j]]
	        # Test de Brunner-Munzel
	        bm_result <- try(lawstat::brunner.munzel.test(x1, x2), silent = silent)
	        if (!inherits(bm_result, "try-error")) {
	          all_pvals <- c(all_pvals, bm_result$p.value)
	          comparisons[[length(comparisons) + 1]] <- c(i, j, bm_result$p.value)
	        } else {
	          all_pvals <- c(all_pvals, NA)
	          comparisons[[length(comparisons) + 1]] <- c(i, j, NA)
	        }
	      }
	    }
	    # Ajustement des p-values
	    adjusted_pvals <- p.adjust(all_pvals, method = p.adjust.method)
	    # Remplir la matrice avec les p-values ajustées
	    idx <- 1
	    for (i in 1:(n_groups - 1)) {
	      for (j in (i + 1):n_groups) {
	        mymat[j - 1, i] <- adjusted_pvals[idx]
	        idx <- idx + 1
	      }
	    }
	    result <- list()
	    result$p.value <- mymat
	    return(result)
	  }
	  # Réordonner les groupes
	  g <- reorder_groups(x, g, ind_control, median)

	  # Effectuer le test
	  result <- BM_to_pairwise(x, g, p.adjust.method = "holm", silent = silent)
	  # Catégoriser les résultats
	  groups <- catego(result, control = control, alpha = alpha, debug = debug)
	  # Bootstrap si demandé
	  if (boot == TRUE) {
	    groups$bootstrap <- bootstrap(x, g, groups$p.value, BM_to_pairwise,
	                                  alpha = alpha, control = control,
	                                  paired = FALSE, debug = debug)
	  }
	  # Réorganiser selon l'ordre initial
	  match(init_order, levels(g)) -> indices
	  groups$groups <- groups$groups[indices, ]
	  if (boot == TRUE) {
	    groups$bootstrap$groups <- groups$bootstrap$groups[indices, ]
	  }

	  .dbg(NULL, "pairwise() - Fin de pairwise() sur Brunner-Munzel.", debug = debug)
	  return(groups)
	} else if (type == "ks" ){
		unique_g <- unique(g)
		# centrer-réduire
		for (i in unique_g) {
			x[g==i] <- (x[g==i]-median(x[g==i]))/sd(x[g==i])
		}
		#return(data.frame(x,g))
		mymat <- matrix(rep(NA, (length(unique_g) - 1)^2),
		                ncol = (length(unique_g) - 1),
		                nrow = (length(unique_g) - 1))
		rownames(mymat) <- unique_g[2:length(unique_g)] ; colnames(mymat) <- unique_g[1:(length(unique_g)-1)]
		#print(mymat)
		ks_func <- function(x,g,mymat,unique_g) {
			#print(unique_g)
			#print(unique_g[1])
			for (i in 1:(length(unique_g)-1)) {
				for (j in (i+1):length(unique_g)) {
					#cat("i",unique_g[i], " et ")
					#cat("j",unique_g[j]," pour ")
					#print(x[g==unique_g[i]])
					#print(x[g==unique_g[j]])

					# Vérifier qu'il y a assez de données pour ks.test (min 2 obs par groupe)
					x_i <- x[g==unique_g[i]]
					x_j <- x[g==unique_g[j]]

					if (length(x_i) >= 2 && length(x_j) >= 2) {
						pv <- ks.test(x_i, x_j)$p.value
					} else {
						# Pas assez de données pour ks.test, retourner NA
						pv <- NA
					}
					#cat("pv",pv,"\n")
					mymat[(j-1),i] <- pv
				}
			}
			return(mymat)
		}
		if (boot == FALSE) {
			mymat <- suppressWarnings(ks_func(x,g,mymat,unique_g)) ; output <- list() ; output$p.value <- mymat
		} else if (boot==TRUE) {
		  mydim <- dim(mymat)
		  myarray <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
		  for (k in 1: iter) {
			x_temp <- x
			for (j in levels(g)) {
			  x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
			}
			temp <- suppressWarnings(ks_func(x_temp,g,mymat,unique_g))
			myarray[k,,] <- temp
		  }
		  output <- list()
		  apply(myarray,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
		  colnames(output$p.value) <- colnames(mymat)
		  rownames(output$p.value) <- rownames(mymat)
		}
		.dbg("pairwise() - End","pairwise() - Fin de pairwise().",debug=debug)
		return(output)
	} else if (type == "lincon" ){
		.dbg(NULL,"pairwise() - lincon()",debug=debug)
		g <- reorder_groups(x, g, ind_control,median)
		result <- lincon_to_pairwise(x,g,tr=tr)
		groups <- catego(result,alpha=alpha, control=control, debug=debug)
		if (boot==TRUE) {
				groups$bootstrap <- bootstrap(x,g,groups$p.value,lincon_to_pairwise,alpha=alpha,
					tr=tr,
					control=control,paired=FALSE, debug=debug)
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		# pas de boot pour lincon...
		.dbg(NULL,"pairwise() - lincon() - Fin de pairwise().",debug=debug)
		return(groups)
	} else if (type == "dunn") {
		if (!requireNamespace("FSA", quietly = TRUE)) {
		  .exit("The 'FSA' package is required for 'dunn' type. Please install it.",NULL)
		}
		dunn_to_pairwise <- function(x,g,p.adjust.method=c(),paired=FALSE) {
			#.dbg(NULL,"Exécution de dunn_to_pairwise().",debug=debug)
			# dunn_to_pairwise ne gère ni les ajustements de p ni les paired=TRUE
			# arg conservés pour normalisé l'appel en bootstrap()
			if (is.ordered(g)) {
				g <- factor(as.character(g), levels = levels(g))
			}
			result <- FSA::dunnTest(x ~ g, method = "holm")
			# Extraire les comparaisons et les p-valeurs ajustées
			comparisons <- result$res[, "Comparison"]
			p_values <- result$res[, "P.adj"]

			# Créer une matrice pairwise pour les p-valeurs (triangulaire inférieure)
			unique_groups <- levels(g)
			#print(unique_groups)
			#print(unique_groups[-1])
			#dimnames <-  list(unique_groups[-1],unique_groups[-length(unique_groups)])
			#print(dimnames)
			result$p.value <- matrix(NA, nrow = length(unique_groups) - 1, ncol = length(unique_groups) - 1)
			rownames(result$p.value) <- unique_groups[-1]
			colnames(result$p.value) <- unique_groups[-length(unique_groups)]
			#print(result$p.value)
			# Remplir la matrice
			for (i in seq_along(comparisons)) {
			  comparison <- strsplit(comparisons[i], " - ")[[1]]
			  #print(comparison)
			  if (length(comparison) == 2) {
				row_group <- comparison[2]
				col_group <- comparison[1]
				# Trouver les indices corrects
				row_index <- match(row_group, rownames(result$p.value))
				col_index <- match(col_group, colnames(result$p.value))

				#print("row");print(row_index)
				#print("col");print(col_index)

				if (!is.na(row_index) && !is.na(col_index)) {
					#print(p_values[i])
				  result$p.value[row_index, col_index] <- p_values[i]
				}
				row_group <- comparison[1]
				col_group <- comparison[2]
				# Trouver les indices corrects
				row_index <- match(row_group, rownames(result$p.value))
				col_index <- match(col_group, colnames(result$p.value))
				if (!is.na(row_index) && !is.na(col_index)) {
				  result$p.value[row_index, col_index] <- p_values[i]
				}
			  } else {
				cat("Invalid comparison:", comparisons[i], "\n")
			  }
			}
			return(result)
		}
		g <- reorder_groups(x, g, ind_control,median)
		dunn_to_pairwise(x,g)-> result
		groups <- catego(result,control=control,alpha=alpha,debug=debug)
		if (boot==TRUE) {
				# paired est mis à FALSE
				groups$bootstrap <- bootstrap(x,g,groups$p.value,dunn_to_pairwise,alpha=alpha,control=control,paired=FALSE)
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		.dbg(NULL,"pairwise() - Fin de pairwise().",debug=debug)
		return(groups)
	#----------------------------------------
	#     Neminye
	#----------------------------------------
	} else if (type == "neminye") {

		nb_id <- table(g)
		if (any(nb_id[1]!=nb_id)) {.exit("Unbalanced data for a Nemenyi test.","Données non équilibrés pour un Neminye.")
		} else {
			id <- rep(1:nb_id[1],time = length(unique(g)))
			result <- frdAllPairsNemenyiTest(x~g|id)
		}
		g <- reorder_groups(x, g, ind_control,median)
print("catego de base")
		groups <- catego(result,control=control,alpha=alpha,debug=debug)
print("av bootstrap")
print(groups$p.value)
		if (boot==TRUE) {
				groups$bootstrap <- bootstrap(x,g,groups$p.value,NemenyTest,alpha=alpha,control=control,paired=TRUE, debug=debug)
		}
print("après bootstrap")
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		.dbg(NULL,"pairwise() - Fin de pairwise() sur Neminye.",debug=debug)
		return(groups)
	#----------------------------------------
	#     Boot
	#----------------------------------------
	} else if (type == "boot"){
		# Note: pairwise.boot() réordonne les groupes en interne par statistique
		# Le réordonnement était auparavant fait ici, mais cela créait une incohérence
		# entre pairwise(type="boot") et pairwise.boot() appelé directement.
		pairwise.boot(x,g,iter=iter,mu=boot_type,debug=debug)-> result
		groups <- catego(result,alpha=alpha,control=control, debug=debug)
		# Restaurer l'ordre initial des groupes pour l'affichage
		reordered_levels <- c(colnames(result$p.value)[1], rownames(result$p.value))
		indices <- match(init_order, reordered_levels)
		if (!any(is.na(indices))) {
			groups$groups <- groups$groups[indices,]
		}
		return(groups)
	} else {.exit("Error! arg type is 'mean','median', 'ks', 'lincon','dunn' or 'boot'.\n",NULL)}
}

