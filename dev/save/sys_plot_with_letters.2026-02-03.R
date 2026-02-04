#' Graphique avec lettres de significativite
#'
#' Fonction interne pour creer boxplot + violin + stripchart avec lettres
#' de significativite superposees selon les resultats de .posthoc()
#'
#' @param x Vecteur numerique - Variable dependante
#' @param g Facteur - Variable groupante
#' @param posthoc_result Liste retournee par .posthoc() ou .posthoc_ANCOVA()
#'   contenant $groups avec colonnes: categories, letters
#' @param main Titre du graphique (defaut: "")
#' @param ylab Etiquette axe Y (defaut: "")
#' @param xlab Etiquette axe X (defaut: "")
#' @param col_box Couleur boxplot (defaut: "cyan")
#' @param col_violin Couleur violinplot (defaut: "#00AA0077")
#' @param col_points Couleur points (defaut: "#FF000088")
#' @param cex_letters Taille lettres (defaut: 1.2)
#' @param offset_letters Decalage vertical lettres en % de range (defaut: 0.05)
#' @param verbose Logique - Afficher messages informatifs
#' @param debug Logique - Afficher messages de debogage
#' @param code Chaine de caracteres - Code de langue (ignore, pour compatibilite)
#'
#' @return NULL (affiche graphique)
#'
#' @details
#' Cette fonction cree un graphique composite:
#' - Boxplot de base (avec mediane, quartiles, whiskers)
#' - Violinplot superpose (si package vioplot disponible)
#' - Stripchart (points individuels avec jitter)
#' - Lettres de significativite au-dessus de chaque groupe
#'
#' L'ordre des groupes sur le graphique respecte l'ordre de posthoc_result$groups,
#' garantissant la coherence avec les tests statistiques.
#'
#' Les lettres sont positionnees automatiquement au-dessus du max de chaque groupe.
#'
#' @note
#' Fonction interne - Non exportee
#' Utilisee par m.test() pour Priority 5 (Mode graphique avance)
#'
#' @keywords internal
.plot_with_letters <- function(x, g,
                               posthoc_result = NULL,
                               main = "",
                               ylab = "",
                               xlab = "",
                               col_box = "cyan",
                               col_violin = "#00AA0077",
                               col_points = "#FF000088",
                               cex_letters = 1.2,
                               offset_letters = 0.05,
                               verbose = FALSE,
                               debug = FALSE,
                               code = NULL) {

  # Reinitialiser le layout graphique (evite problemes si layout precedent actif)
  # Sauvegarde des parametres actuels pour restauration ulterieure
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # Reset a un layout 1x1 standard
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

  # Validation
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (!is.factor(g)) {
    g <- factor(g)
  }

  # Creer data.frame pour compatibilite avec formule
  plot_data <- data.frame(x = x, g = g)

  # ETAPE 1: Boxplot de base
  ylim <- c(min(x),max(x)+0.05*(max(x)-min(x)))
  bp <- boxplot(x ~ g, data = plot_data,
                col = col_box,ylim=ylim,
                main = main,
                ylab = ylab,
                xlab = xlab,
                las = 1)  # Etiquettes horizontales

  # ETAPE 2: Violinplot (si disponible)
  if (requireNamespace("vioplot", quietly = TRUE)) {
    vioplot::vioplot(x ~ g, data = plot_data,
                     col = col_violin,
                     add = TRUE,
                     wex = 0.8)  # Largeur violons
  }

  # ETAPE 3: Stripchart (points individuels)
  n_groups <- length(levels(g))
  jitter_amount <- 1 / (n_groups + 2)

  stripchart(x ~ g, data = plot_data,
             col = col_points,
             pch = 16,
             vertical = TRUE,
             add = TRUE,
             method = "jitter",
             jitter = jitter_amount)

  # ETAPE 4: Lettres de significativite
  if (!is.null(posthoc_result) && !is.null(posthoc_result$groups)) {

    groups_df <- posthoc_result$groups

    # Detecter colonne contenant les lettres
    # Chercher colonne avec nom contenant "group" ou qui a des lettres (a,b,c...)
    letter_col_idx <- NULL

    for (i in seq_along(names(groups_df))) {
      col_name <- names(groups_df)[i]
      # Skip colonne "categories" ou colonnes numeriques
      if (grepl("categor|mean|std|se", col_name, ignore.case = TRUE)) {
        next
      }
      # Verifier si contient des lettres (a, b, c, ab, etc.)
      col_values <- as.character(groups_df[[i]])
      if (any(grepl("^[a-z]+$", col_values, ignore.case = TRUE))) {
        letter_col_idx <- i
        break
      }
    }

    if (is.null(letter_col_idx)) {
      # Fallback: deuxieme colonne si pas trouve
      if (ncol(groups_df) >= 2) {
        letter_col_idx <- 2
      }
    }

    if (!is.null(letter_col_idx)) {

      # Extraire categories et lettres
      categories <- as.character(groups_df[[1]])  # Premiere colonne = categories
      letters <- as.character(groups_df[[letter_col_idx]])

      # S'assurer que l'ordre correspond aux niveaux de g
      g_levels <- levels(g)

      # Calculer position verticale pour chaque lettre
      # Positionner au-dessus du max de chaque groupe + offset
      y_range <- diff(range(x, na.rm = TRUE))
      offset <- offset_letters * y_range

      for (i in seq_along(g_levels)) {
        group_level <- g_levels[i]

        # Trouver la lettre correspondante
        idx_letter <- which(categories == group_level)

        if (length(idx_letter) > 0) {
          current_letter <- letters[idx_letter[1]]

          # Calculer y position: max du groupe + offset
          group_values <- x[g == group_level]
          if (length(group_values) > 0) {
            y_pos <- max(group_values, na.rm = TRUE) + offset

            # Ajouter la lettre
            text(x = i, y = y_pos,
                 labels = current_letter,
                 cex = cex_letters,
                 font = 2,  # Bold
                 col = "black")
          }
        }
      }

      if (debug) {
        cat(sprintf("[.plot_with_letters] Lettres ajout\u00e9es: %s\n",
                    paste(letters, collapse = ", ")))
      }
    } else {
      if (verbose) {
        cat("Note: Pas de colonne 'lettres' trouv\u00e9e dans posthoc$groups\n")
      }
    }
  } else {
    if (debug) {
      cat("[.plot_with_letters] Pas de r\u00e9sultats post-hocs fournis - graphique sans lettres\n")
    }
  }

  invisible(NULL)
}
