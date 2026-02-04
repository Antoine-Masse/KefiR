#' Function to compute Cook's distance for mixed models
#'
#' @param model A mixed-effects model of class lmerMod (fitted with lme4).
#' @param obs Logical, if TRUE, computes Cook's distance for individual observations as well.
#'
#' @return A named vector with Cook's distance values for each group of random effects and optionally for each observation.
#' @examples
#' # Example 1
#' library(lme4)
#' data("sleepstudy")
#' model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' cooks_distances <- cooks.distance_lmer(model)
#' 
#' # Example 2
#' data(mtcars)
#' model <- lmer(cyl ~ disp + (1 | hp) + (1 | gear), data = mtcars)
#' cooks_distances <- cooks.distance_lmer(model)
#' @export

# Internal function to calculate Cook's distance for mixed models
cooks.distance_lmer <- function(model, obs = FALSE) {
  # Fonction pour détecter la langue et afficher un message en conséquence
  detect_lang <- function(message_fr, message_en) {
    locale <- Sys.getlocale("LC_TIME")
    # Détecter la langue à partir de la locale
    if (grepl("fr_", locale)) {
      message <- message_fr
    } else {
      message <- message_en
    }
    return(message)
  }
  
  # Utilisation de la fonction
  message_fr <- "Le modèle doit être de type lmerMod (ajusté avec lme4)."
  message_en <- "The model must be of type lmerMod (adjusted with lme4)."
  if (!inherits(model, "lmerMod")) {
    message <- detect_lang(message_fr, message_en)
    stop(message)
  }
  
  # Extraire les groupes d'effets aléatoires du modèle
  random_effects <- names(ranef(model))
  
  # Initialiser une liste pour stocker les distances de Cook pour chaque groupe d'effets aléatoires
  cooks_dist_all <- list()
  
  # Calculer la distance de Cook pour chaque groupe d'effets aléatoires
  for (group in random_effects) {
    influ <- influence(model = model, group = group)
    cooks_dist_group <- cooks.distance(influ)
    # Ajouter les noms des sous-groupes
    names(cooks_dist_group) <- rownames(ranef(model)[[group]])
    # Ajouter les distances de Cook à la liste avec le nom du groupe
    cooks_dist_all[[group]] <- cooks_dist_group
  }
  
  # Calculer la distance de Cook pour les observations si spécifié
  if (obs) {
    influ_obs <- influence(model = model, obs = TRUE)
    cooks_dist_obs <- cooks.distance(influ_obs)
    names(cooks_dist_obs) <- paste("Observation", seq_along(cooks_dist_obs), sep = ":")
    # Ajouter les distances de Cook des observations à la liste avec le nom "Observations"
    cooks_dist_all[["Observations"]] <- cooks_dist_obs
  }
  
  return(cooks_dist_all)
}	


