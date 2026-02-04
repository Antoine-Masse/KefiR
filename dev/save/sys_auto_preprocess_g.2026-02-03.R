auto_preprocess_g <- function(
  x,
  g,
  gamma = 5,
  start_bins = 10,
  cor_threshold = 0.8,
  k = 0,
  debug = FALSE,
  verbose = FALSE
) {
  #---------------------------------------------------------
  # 0) Vérifications basiques
  #---------------------------------------------------------
  if (!is.numeric(x)) {
    .exit("`x` must be numeric (response variable).","`x` doit être numérique (variable réponse).")
  }
  if (!is.data.frame(g)) {
    .exit("`g` must be a data.frame of explanatory variables.","`g` doit être un data.frame de variables explicatives.")
  }
  if (length(x) != nrow(g)) {
    .exit("`x` and `g` must have the same length (same number of observations).","`x` et `g` doivent avoir la même longueur (même nombre d'observations).")
  }

  if (debug) {
    cat(.msg("\n[DEBUG] Start of auto_preprocess_g\n","\n[DEBUG] Début de auto_preprocess_g\n"))
    cat("[DEBUG] length(x) =", length(x), "\n")
    cat("[DEBUG] dim(g)    =", dim(g), "\n")
    cat("[DEBUG] colnames(g) =", colnames(g), "\n\n")
  }

  #---------------------------------------------------------
  # 1) Identifier types de variables (facteurs, binaires, continues)
  #---------------------------------------------------------
  numeric_vars <- character()
  factor_vars  <- character()

  for (colname in colnames(g)) {
    vec <- g[[colname]]
    if (is.numeric(vec)) {
      # Vérif binaire
      unique_vals <- unique(na.omit(vec))
      if (length(unique_vals) == 2 && all(sort(unique_vals) == c(0,1))) {
        g[[colname]] <- factor(vec, levels = c(0,1))
        factor_vars  <- c(factor_vars, colname)
        if (debug)   cat(.msg(paste0("[DEBUG] ", colname, " => binary => transformed into factor.\n"),paste0("[DEBUG] ", colname, " => binaire => transformé en facteur.\n")))
      } else {
        numeric_vars <- c(numeric_vars, colname)
        if (debug) cat(.msg(paste0("[DEBUG] ", colname, " => numeric variable.\n"),
                 paste0("[DEBUG] ", colname, " => variable numérique.\n")))
      }
    } else {
      # Déjà facteur ou character
      if (!is.factor(vec)) {
        g[[colname]] <- factor(vec)
      }
      factor_vars <- c(factor_vars, colname)
      if (debug) cat(.msg(paste0("[DEBUG] ", colname, " => factor variable.\n"),
                 paste0("[DEBUG] ", colname, " => variable factor.\n")))
    }
  }

	if (debug) {
	  cat(
		.msg(
		  paste0(
			"\n[DEBUG] numeric_vars = ", paste(numeric_vars, collapse = ", "), "\n",
			"[DEBUG] factor_vars  = ", paste(factor_vars, collapse = ", "), "\n\n"
		  ),
		  paste0(
			"\n[DEBUG] numeric_vars = ", paste(numeric_vars, collapse = ", "), "\n",
			"[DEBUG] factor_vars  = ", paste(factor_vars, collapse = ", "), "\n\n"
		  )
		)
	  )
	}


# Trier les variables numériques par corrélation absolue avec x
if (length(numeric_vars) > 1) {
    cor_x <- sapply(numeric_vars, function(vn) {
        val <- suppressWarnings(cor(x, g[[vn]], use="pairwise.complete.obs"))
        if (is.na(val)) 0 else abs(val)  # On prend la valeur absolue
    })
    # Trier numeric_vars par corrélation décroissante
    numeric_vars <- numeric_vars[order(cor_x, decreasing = TRUE)]
	if (debug) cat(
	  .msg(
		paste0("[DEBUG] Order of numeric variables for binning (sorted by absolute correlation): ", 
			   paste(numeric_vars, collapse = ", "), "\n"),
		paste0("[DEBUG] Ordre des variables numériques pour le binning (par corrélation absolue) : ", 
			   paste(numeric_vars, collapse = ", "), "\n")
	  )
	)
}


  #---------------------------------------------------------
  # 2) Exclure les variables numériques trop corrélées entre elles
  #---------------------------------------------------------
  if (length(numeric_vars) > 1) {
    df_num <- g[, numeric_vars, drop=FALSE]
    colnames(df_num) <- numeric_vars

    mat_cor <- suppressWarnings(cor(df_num, use="pairwise.complete.obs"))
    dimnames(mat_cor) <- list(numeric_vars, numeric_vars)

	if (debug) {
	  cat(
		.msg(
		  paste0(
			"[DEBUG] dim(mat_cor) = ", paste(dim(mat_cor), collapse = " x "), "\n",
			"[DEBUG] rownames(mat_cor) = ", paste(rownames(mat_cor), collapse = ", "), "\n",
			"[DEBUG] colnames(mat_cor) = ", paste(colnames(mat_cor), collapse = ", "), "\n",
			"[DEBUG] mat_cor :\n"
		  ),
		  paste0(
			"[DEBUG] dim(mat_cor) = ", paste(dim(mat_cor), collapse = " x "), "\n",
			"[DEBUG] rownames(mat_cor) = ", paste(rownames(mat_cor), collapse = ", "), "\n",
			"[DEBUG] colnames(mat_cor) = ", paste(colnames(mat_cor), collapse = ", "), "\n",
			"[DEBUG] mat_cor :\n"
		  )
		)
	  )
	  print(mat_cor)
	  cat("\n")
	}

    cor_x <- sapply(numeric_vars, function(vn) {
      val <- suppressWarnings(cor(x, g[[vn]], use="pairwise.complete.obs"))
      if (is.na(val)) 0 else val
    })
	if (debug) {
	  cat(.msg("[DEBUG] cor_x :\n", "[DEBUG] cor_x :\n"))
	  print(cor_x)
	  cat("\n")
	}

    to_remove <- character()

    # *** BOUCLE STANDARD ***

    for (i in 1:(length(numeric_vars)-1))  {
      for (j in seq(i+1, length(numeric_vars))) {
        var_i <- numeric_vars[i]
        var_j <- numeric_vars[j]

		if (debug) cat(
		  .msg(
			paste0("[DEBUG] Compare pair (i=", i, ", j=", j, ") : ", var_i, " - ", var_j, "\n"),
			paste0("[DEBUG] Comparer la paire (i=", i, ", j=", j, ") : ", var_i, " - ", var_j, "\n")
		  )
		)
        if (!(var_i %in% to_remove) && !(var_j %in% to_remove)) {
          rho_ij <- mat_cor[var_i, var_j]
			if (debug) cat(
			  .msg(
				paste0("   [DEBUG]  rho_ij = ", rho_ij, "\n"),
				paste0("   [DEBUG]  rho_ij = ", rho_ij, "\n")
			  )
			)

          if (abs(rho_ij) >= cor_threshold) {
            # exclure la moins corrélée à x
            if (abs(cor_x[var_i]) < abs(cor_x[var_j])) {
              to_remove <- c(to_remove, var_i)
              if (debug) cat("   [DEBUG]  => remove var_i =", var_i, "\n")
            } else {
              to_remove <- c(to_remove, var_j)
              if (debug) cat("   [DEBUG]  => remove var_j =", var_j, "\n")
            }
          }
        }
      }
    }

    numeric_vars <- setdiff(numeric_vars, to_remove)
	if (debug && length(to_remove) > 0) {
	  cat(
		.msg(
		  paste0(
			"[DEBUG] Exclusion for corr >= ", cor_threshold, ": ", paste(to_remove, collapse = ", "), "\n",
			"[DEBUG] numeric_vars final = ", paste(numeric_vars, collapse = ", "), "\n"
		  ),
		  paste0(
			"[DEBUG] Exclusion pour corr >= ", cor_threshold, " : ", paste(to_remove, collapse = ", "), "\n",
			"[DEBUG] numeric_vars final = ", paste(numeric_vars, collapse = ", "), "\n"
		  )
		)
	  )
	}
  } else if (debug) {
		 cat(
		  .msg(
			"[DEBUG] Less than 2 numeric variables => no comparison.\n",
			"[DEBUG] Moins de 2 variables numériques => pas de comparaison.\n"
		  )
		)
  }

  #---------------------------------------------------------
  # 3) Binning
  #---------------------------------------------------------
  cross_nobs <- function(...) {
    tab_ <- table(...)
    min(tab_)
  }

  create_bins <- function(vec, nbins) {
    probs <- seq(0, 1, length.out=nbins+1)
    brks  <- unique(quantile(vec, probs, na.rm=TRUE))
    if (length(brks)<2) {
      return(rep(NA, length(vec)))
    }
    cut(vec, breaks=brks, include.lowest=TRUE, right=TRUE)
  }

  for (v in numeric_vars) {
    vec <- g[[v]]
    valid_idx <- which(!is.na(x) & !is.na(vec))
    if (length(valid_idx) < (2*gamma)) {
      if (debug) cat(
		  .msg(
			paste0("[DEBUG] => No binning on ", v, " (n < 2*gamma)\n"),
			paste0("[DEBUG] => Pas de binning sur ", v, " (n < 2*gamma)\n")
		  )
		)
      next
    }

    bins_candidate <- start_bins
    success <- FALSE
    binned_factor <- NULL

    repeat {
      if (bins_candidate<2) break
      tmp_cut <- create_bins(vec, bins_candidate)
      if (all(is.na(tmp_cut))) {
        if (debug) cat(
		  .msg(
			paste0("[DEBUG] ", v, " bins=", bins_candidate, " => all NA => decreasing\n"),
			paste0("[DEBUG] ", v, " bins=", bins_candidate, " => tout NA => on descend\n")
		  )
		)
        bins_candidate <- bins_candidate - 1
        next
      }
      df_tmp <- data.frame(tmp_cut, g[factor_vars], stringsAsFactors=TRUE)
      min_n <- cross_nobs(df_tmp)
      if (debug) cat(
		  .msg(
			paste0("[DEBUG] Trying binning ", v, " -> bins=", bins_candidate, " => min freq = ", min_n, "\n"),
			paste0("[DEBUG] Essai binning ", v, " -> bins=", bins_candidate, " => min freq = ", min_n, "\n")
		  )
		)

      if (min_n>=gamma) {
        binned_factor <- tmp_cut
        success <- TRUE
		if (debug) cat(
		  .msg(
			paste0("[DEBUG] => ", v, " transformed to factor (", bins_candidate, " classes )\n"),
			paste0("[DEBUG] => ", v, " transformée en facteur (", bins_candidate, " classes )\n")
		  )
		)
        break
      } else {
        bins_candidate <- bins_candidate - 1
      }
    }
    if (success) {
      g[[v]] <- binned_factor
      if (!(v %in% factor_vars)) factor_vars <- c(factor_vars, v)
    } else {
      if (debug) cat(
		  .msg(
			paste0("[DEBUG] => Binning failed on ", v, " => remains numeric.\n"),
			paste0("[DEBUG] => Échec binning sur ", v, " => reste numérique.\n")
		  )
		)
    }
  }
  # Affichage du niveau de binning atteint si verbose est activé
  if (isTRUE(verbose)) {
    # Calcul du nombre de variables factorielles et du total des niveaux
    total_vars <- length(factor_vars)
    total_levels <- sum(sapply(factor_vars, function(v) {
      # Vérifie si la colonne existe dans g et est un facteur, puis compte les niveaux
      if (v %in% colnames(g) && is.factor(g[[v]])) length(levels(g[[v]])) else 0
    }))

    cat(
      paste0(
        k,
        .msg(
          paste0(") Binning level reached: ", total_vars, 
                 " variables with a total of ", total_levels, " levels.\n"),
          paste0(") Niveau de binning atteint : ", total_vars, 
                 " variables avec un total de ", total_levels, " niveaux.\n")
        )
      )
    )
  }
	if (debug) cat(
	  .msg(
		"\n[DEBUG] End of auto_preprocess_g\n",
		"\n[DEBUG] Fin de auto_preprocess_g\n"
	  )
	)
  return(g)
}
