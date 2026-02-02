#' Graphique avec lettres de significativité
#'
#' Fonction interne pour créer boxplot + violin + stripchart avec lettres
#' de significativité superposées selon les résultats de .posthoc()
#'
#' @param x Vecteur numérique - Variable dépendante
#' @param g Facteur - Variable groupante
#' @param posthoc_result Liste retournée par .posthoc() ou .posthoc_ANCOVA()
#'   contenant $groups avec colonnes: categories, letters
#' @param main Titre du graphique (défaut: "")
#' @param ylab Étiquette axe Y (défaut: "")
#' @param xlab Étiquette axe X (défaut: "")
#' @param col_box Couleur boxplot (défaut: "cyan")
#' @param col_violin Couleur violinplot (défaut: "#00AA0077")
#' @param col_points Couleur points (défaut: "#FF000088")
#' @param cex_letters Taille lettres (défaut: 1.2)
#' @param offset_letters Décalage vertical lettres en % de range (défaut: 0.05)
#' @param verbose Logique - Afficher messages informatifs
#' @param debug Logique - Afficher messages de débogage
#' @param code Chaîne de caractères - Code de langue (ignoré, pour compatibilité)
#'
#' @return NULL (affiche graphique)
#'
#' @details
#' Cette fonction crée un graphique composite:
#' - Boxplot de base (avec médiane, quartiles, whiskers)
#' - Violinplot superposé (si package vioplot disponible)
#' - Stripchart (points individuels avec jitter)
#' - Lettres de significativité au-dessus de chaque groupe
#'
#' L'ordre des groupes sur le graphique respecte l'ordre de posthoc_result$groups,
#' garantissant la cohérence avec les tests statistiques.
#'
#' Les lettres sont positionnées automatiquement au-dessus du max de chaque groupe.
#'
#' @note
#' Fonction interne - Non exportée
#' Utilisée par m.test() pour Priority 5 (Mode graphique avancé)
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

  # Réinitialiser le layout graphique (évite problèmes si layout précédent actif)
  # Sauvegarde des paramètres actuels pour restauration ultérieure
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # Reset à un layout 1x1 standard
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

  # Validation
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (!is.factor(g)) {
    g <- factor(g)
  }

  # Créer data.frame pour compatibilité avec formule
  plot_data <- data.frame(x = x, g = g)

  # ÉTAPE 1: Boxplot de base
  ylim <- c(min(x),max(x)+0.05*(max(x)-min(x)))
  bp <- boxplot(x ~ g, data = plot_data,
                col = col_box,ylim=ylim,
                main = main,
                ylab = ylab,
                xlab = xlab,
                las = 1)  # Étiquettes horizontales

  # ÉTAPE 2: Violinplot (si disponible)
  if (requireNamespace("vioplot", quietly = TRUE)) {
    vioplot::vioplot(x ~ g, data = plot_data,
                     col = col_violin,
                     add = TRUE,
                     wex = 0.8)  # Largeur violons
  }

  # ÉTAPE 3: Stripchart (points individuels)
  n_groups <- length(levels(g))
  jitter_amount <- 1 / (n_groups + 2)

  stripchart(x ~ g, data = plot_data,
             col = col_points,
             pch = 16,
             vertical = TRUE,
             add = TRUE,
             method = "jitter",
             jitter = jitter_amount)

  # ÉTAPE 4: Lettres de significativité
  if (!is.null(posthoc_result) && !is.null(posthoc_result$groups)) {

    groups_df <- posthoc_result$groups

    # Détecter colonne contenant les lettres
    # Chercher colonne avec nom contenant "group" ou qui a des lettres (a,b,c...)
    letter_col_idx <- NULL

    for (i in seq_along(names(groups_df))) {
      col_name <- names(groups_df)[i]
      # Skip colonne "categories" ou colonnes numériques
      if (grepl("categor|mean|std|se", col_name, ignore.case = TRUE)) {
        next
      }
      # Vérifier si contient des lettres (a, b, c, ab, etc.)
      col_values <- as.character(groups_df[[i]])
      if (any(grepl("^[a-z]+$", col_values, ignore.case = TRUE))) {
        letter_col_idx <- i
        break
      }
    }

    if (is.null(letter_col_idx)) {
      # Fallback: deuxième colonne si pas trouvé
      if (ncol(groups_df) >= 2) {
        letter_col_idx <- 2
      }
    }

    if (!is.null(letter_col_idx)) {

      # Extraire catégories et lettres
      categories <- as.character(groups_df[[1]])  # Première colonne = catégories
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
        cat(sprintf("[.plot_with_letters] Lettres ajoutées: %s\n",
                    paste(letters, collapse = ", ")))
      }
    } else {
      if (verbose) {
        cat("Note: Pas de colonne 'lettres' trouvée dans posthoc$groups\n")
      }
    }
  } else {
    if (debug) {
      cat("[.plot_with_letters] Pas de résultats post-hocs fournis - graphique sans lettres\n")
    }
  }

  invisible(NULL)
}
