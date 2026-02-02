################################################################################
#
#                           .pairwise_medpb2()
#                     Wrapper for WRS2::medpb2 pairwise comparisons
#
################################################################################

#' Pairwise median comparisons using WRS2::medpb2 with Holm correction
#'
#' @description
#' Performs pairwise comparisons of medians using percentile bootstrap
#' (medpb2 from WRS2 package) with Holm correction for multiple comparisons.
#' Designed for k>2 groups with strong asymmetry or heavy tails (med1way context).
#'
#' @param x Numeric vector of observations
#' @param g Factor of grouping variable
#' @param alpha Significance level (default 0.05)
#' @param control Optional control group name for Dunnett-style comparisons
#' @param nboot Number of bootstrap iterations (default 2000, as per WRS2::medpb2)
#' @param debug Logical for debug messages
#'
#' @return List with:
#'   \itemize{
#'     \item groups: data.frame with categories and compact letter display
#'     \item p.value: matrix of pairwise p-values (raw, before Holm adjustment)
#'     \item p.adjusted: matrix of Holm-adjusted p-values
#'     \item method: "medpb2_Holm"
#'   }
#'
#' @details
#' medpb2() from WRS2 compares two independent groups using medians with
#' percentile bootstrap. It is the only known method that works well in
#' simulations when tied values are likely to occur (Wilcox, 2017).
#'
#' This wrapper extends medpb2() to k>2 groups by:
#' 1. Calling medpb2() on all pairs (i,j) where i<j
#' 2. Applying Holm correction for multiple comparisons
#' 3. Generating compact letter display
#'
#' @references
#' Wilcox, R. R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.).
#' Academic Press. Chapter 7: Comparing two groups.
#' DOI: 10.1016/C2010-0-67044-1
#'
#' @keywords internal
#' @export
.pairwise_medpb2 <- function(x, g, alpha = 0.05, control = NULL, nboot = 2000, debug = FALSE) {

  # ============================================================================
  # VALIDATION ET PRÉPARATION
  # ============================================================================

  # Coercitions
  if (is.data.frame(x)) x <- x[[1]]
  if (is.matrix(x)) x <- x[, 1, drop = TRUE]
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.vector(x)

  if (is.data.frame(g)) g <- g[[1]]
  if (is.matrix(g)) g <- g[, 1, drop = TRUE]
  if (is.list(g)) g <- unlist(g, use.names = FALSE)

  g <- droplevels(factor(g))
  lev <- levels(g)
  ng <- length(lev)

  if (ng < 2) {
    stop(".pairwise_medpb2(): Need at least 2 groups")
  }

  if (ng == 2) {
    warning(".pairwise_medpb2(): Only 2 groups detected. Consider using medpb2() directly.")
  }

  # ============================================================================
  # MATRICES DE P-VALUES
  # ============================================================================

  # Matrice de p-values brutes
  p_matrix <- matrix(1, nrow = ng, ncol = ng)
  rownames(p_matrix) <- colnames(p_matrix) <- lev
  diag(p_matrix) <- 1

  # Matrice de p-values ajustées
  p_adjusted_matrix <- matrix(1, nrow = ng, ncol = ng)
  rownames(p_adjusted_matrix) <- colnames(p_adjusted_matrix) <- lev
  diag(p_adjusted_matrix) <- 1

  # Vecteur pour stocker toutes les p-values (pour Holm)
  all_pvals <- c()
  pair_indices <- list()

  # ============================================================================
  # GESTION DU CONTRÔLE (Dunnett-style)
  # ============================================================================

  if (!is.null(control)) {
    control_chr <- as.character(control)[1]
    if (!control_chr %in% lev) {
      warning(paste0(".pairwise_medpb2(): Control '", control_chr, "' not found in groups. Using all pairwise comparisons."))
      control <- NULL
    }
  }

  # ============================================================================
  # COMPARAISONS PAR PAIRES
  # ============================================================================

  if (!is.null(control)) {
    # ------------------------------------------------------------------
    # Mode Dunnett: seulement vs contrôle
    # ------------------------------------------------------------------
    control_chr <- as.character(control)[1]
    ind_ctrl <- match(control_chr, lev)

    for (i in seq_along(lev)) {
      if (i == ind_ctrl) next

      # Extraire données des 2 groupes
      x1 <- x[g == lev[ind_ctrl]]
      x2 <- x[g == lev[i]]

      # Vérifier que les groupes ont des données
      if (length(x1) < 2 || length(x2) < 2) {
        .dbg(paste0("Skipping ", lev[ind_ctrl], " vs ", lev[i], ": insufficient data"),
             paste0("Ignorer ", lev[ind_ctrl], " vs ", lev[i], " : données insuffisantes"),
             debug = debug)
        next
      }

      # Appeler medpb2
      test_result <- tryCatch({
        WRS2::medpb2(formula = value ~ group,
                     data = data.frame(value = c(x1, x2),
                                      group = factor(c(rep(lev[ind_ctrl], length(x1)),
                                                      rep(lev[i], length(x2))))),
                     nboot = nboot)
      }, error = function(e) {
        .dbg(paste0("medpb2 failed for ", lev[ind_ctrl], " vs ", lev[i], ": ", e$message),
             paste0("medpb2 échoué pour ", lev[ind_ctrl], " vs ", lev[i], " : ", e$message),
             debug = debug)
        return(NULL)
      })

      if (!is.null(test_result) && !is.null(test_result$p.value)) {
        p_matrix[ind_ctrl, i] <- p_matrix[i, ind_ctrl] <- test_result$p.value
        all_pvals <- c(all_pvals, test_result$p.value)
        pair_indices[[length(pair_indices) + 1]] <- c(ind_ctrl, i)
      }
    }

  } else {
    # ------------------------------------------------------------------
    # Mode all pairwise
    # ------------------------------------------------------------------
    for (i in 1:(ng-1)) {
      for (j in (i+1):ng) {
        # Extraire données des 2 groupes
        x1 <- x[g == lev[i]]
        x2 <- x[g == lev[j]]

        # Vérifier que les groupes ont des données
        if (length(x1) < 2 || length(x2) < 2) {
          .dbg(paste0("Skipping ", lev[i], " vs ", lev[j], ": insufficient data"),
               paste0("Ignorer ", lev[i], " vs ", lev[j], " : données insuffisantes"),
               debug = debug)
          next
        }

        # Appeler medpb2
        test_result <- tryCatch({
          WRS2::medpb2(formula = value ~ group,
                       data = data.frame(value = c(x1, x2),
                                        group = factor(c(rep(lev[i], length(x1)),
                                                        rep(lev[j], length(x2))))),
                       nboot = nboot)
        }, error = function(e) {
          .dbg(paste0("medpb2 failed for ", lev[i], " vs ", lev[j], ": ", e$message),
               paste0("medpb2 échoué pour ", lev[i], " vs ", lev[j], " : ", e$message),
               debug = debug)
          return(NULL)
        })

        if (!is.null(test_result) && !is.null(test_result$p.value)) {
          p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
          all_pvals <- c(all_pvals, test_result$p.value)
          pair_indices[[length(pair_indices) + 1]] <- c(i, j)
        }
      }
    }
  }

  # ============================================================================
  # CORRECTION DE HOLM
  # ============================================================================

  if (length(all_pvals) > 0) {
    adjusted_pvals <- p.adjust(all_pvals, method = "holm")

    # Remplir matrice ajustée
    for (idx in seq_along(pair_indices)) {
      i <- pair_indices[[idx]][1]
      j <- pair_indices[[idx]][2]
      p_adjusted_matrix[i, j] <- p_adjusted_matrix[j, i] <- adjusted_pvals[idx]
    }

    # ============================================================================
    # COMPACT LETTER DISPLAY
    # ============================================================================

    # Créer matrice de significativité
    signif_matrix <- (p_adjusted_matrix <= alpha)
    diag(signif_matrix) <- FALSE  # Groupe avec lui-même : non significatif

    # Algorithme de génération de lettres
    # Approche : groupes qui ne diffèrent pas significativement partagent lettres
    letters_vec <- rep("", ng)

    # Initialiser avec "a" pour tous
    letters_vec <- rep("a", ng)

    # Pour chaque groupe, vérifier s'il diffère des précédents
    for (i in 2:ng) {
      # Trouver groupes différents de i
      differs_from <- which(signif_matrix[i, 1:(i-1)])

      if (length(differs_from) > 0) {
        # i diffère d'au moins un groupe précédent
        # Trouver la lettre maximale utilisée jusqu'ici
        max_letter_code <- max(sapply(letters_vec[1:(i-1)], function(l) {
          if (nchar(l) == 0) return(0)
          utf8ToInt(substr(l, 1, 1)) - utf8ToInt("a") + 1
        }))

        # Assigner lettre suivante
        letters_vec[i] <- intToUtf8(utf8ToInt("a") + max_letter_code)
      }
      # Sinon, garder "a" (ne diffère de personne)
    }

    groups_df <- data.frame(
      categories = lev,
      medpb2_Holm = letters_vec,
      stringsAsFactors = FALSE
    )

  } else {
    # Aucune comparaison réussie
    warning(".pairwise_medpb2(): No successful pairwise comparisons.")
    groups_df <- data.frame(
      categories = lev,
      medpb2_Holm = rep("", ng),
      stringsAsFactors = FALSE
    )
  }

  # ============================================================================
  # RETOUR
  # ============================================================================

  return(list(
    groups = groups_df,
    p.value = p_matrix,
    p.adjusted = p_adjusted_matrix,
    method = "medpb2_Holm"
  ))
}
