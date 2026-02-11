#' Génère une nouvelle variable corrélée à la première colonne d'un dataframe
#'
#' Cette fonction génère une nouvelle variable numérique qui est corrélée à la première colonne d'un dataframe fourni. Elle permet de spécifier le coefficient de corrélation souhaité, la moyenne et l'écart-type de la nouvelle variable, ainsi que des options pour contrôler la précision de la corrélation et des paramètres statistiques.
#'
#' @param dataframe Un dataframe contenant au moins une colonne numérique.
#' @param cor.f Un nombre entre -1 et 1 spécifiant le coefficient de corrélation souhaité entre la nouvelle variable et la première colonne du dataframe.
#' @param meani (Optionnel) La moyenne souhaitée pour la nouvelle variable. Par défaut à 0.
#' @param sdi (Optionnel) L'écart-type souhaité pour la nouvelle variable. Par défaut à 1.
#' @param strict (Optionnel) Booléen indiquant si la nouvelle variable doit avoir exactement la moyenne et l'écart-type spécifiés (\code{TRUE}), ou si elle peut fluctuer naturellement (\code{FALSE}). Par défaut à \code{FALSE}.
#' @param strict_corr (Optionnel) Booléen indiquant si la corrélation entre la nouvelle variable et la première colonne doit être exactement égale à \code{cor.f} (\code{TRUE}), ou si elle peut fluctuer légèrement (\code{FALSE}). Par défaut à \code{FALSE}.
#'
#' @return Un vecteur numérique de même longueur que le nombre de lignes du dataframe, représentant la nouvelle variable générée.
#'
#' @details
#' La fonction suit les étapes suivantes :
#' \enumerate{
#'   \item Convertit le dataframe en matrice.
#'   \item Génère une variable aléatoire \code{z} selon une distribution normale avec la moyenne \code{meani} et l'écart-type \code{sdi}.
#'   \item Centre et standardise la première colonne du dataframe.
#'   \item Selon la valeur de \code{strict_corr} :
#'     \begin{itemize}
#'       \item Si \code{TRUE}, orthogonalise \code{z} par rapport à la première colonne pour obtenir une corrélation exacte.
#'       \item Si \code{FALSE}, utilise \code{z} standardisé pour générer la nouvelle variable avec des fluctuations naturelles de la corrélation.
#'     \end{itemize}
#'   \item Calcule la nouvelle variable en combinant la première colonne standardisée et \code{z} (orthogonalisé ou standardisé), pondérés par \code{cor.f} et \code{sqrt(1 - cor.f^2)} respectivement.
#'   \item Ajuste la moyenne et l'écart-type de la nouvelle variable selon la valeur de \code{strict} :
#'     \begin{itemize}
#'       \item Si \code{TRUE}, impose exactement \code{meani} et \code{sdi}.
#'       \item Si \code{FALSE}, conserve les fluctuations naturelles autour de \code{meani} et \code{sdi}.
#'     \end{itemize}
#' }
#'
#' @examples
#' # Exemple d'utilisation
#' df <- data.frame(x = rnorm(100))
#' nouvelle_var <- cor.e(df, cor.f = 0.5, meani = 10, sdi = 2)
#' # Vérifier la corrélation
#' cor(df$x, nouvelle_var)
#' # Vérifier la moyenne et l'écart-type
#' mean(nouvelle_var)
#' sd(nouvelle_var)
#'
#' @export
cor.e <- function(dataframe, cor.f, meani=0, sdi=1, strict=FALSE, strict_corr=FALSE){
  if (is.vector(dataframe)) dataframe <- data.frame(V1=dataframe)
  if (!is.data.frame(dataframe)) stop("L'entrée doit être un dataframe ou un vecteur.")
  if (!is.numeric(dataframe[[1]])) stop("La première colonne doit être numérique.")
  if (cor.f < -1 || cor.f > 1) stop("cor.f dans [-1,1].")
  if (sdi <= 0) stop("sdi > 0.")
  mat.ini <- as.matrix(dataframe)

  # helpers sûrs
  safe_sd <- function(v){ s <- sd(v, na.rm=TRUE); if(!is.finite(s) || s==0) 1e-8 else s }
  safe_var <- function(v){ vv <- var(v, na.rm=TRUE); if(!is.finite(vv) || vv==0) 1e-8 else vv }

  z <- rnorm(nrow(mat.ini), mean=meani, sd=sdi)
  z_cent <- z - mean(z)

  x1 <- mat.ini[,1]
  x1_cent <- x1 - mean(x1)
  x1_std  <- x1_cent / safe_sd(x1_cent)

  if (strict_corr) {
    denom <- sum(x1_cent^2); if(!is.finite(denom) || denom==0) denom <- 1e-8
    proj_coef <- sum(z_cent * x1_cent) / denom
    z_orth <- z_cent - proj_coef * x1_cent
    z_orth_std <- z_orth / safe_sd(z_orth)
    new_var <- cor.f * x1_std + sqrt(1 - cor.f^2) * z_orth_std
  } else {
    z_std <- z_cent / safe_sd(z_cent)        # <- sécurisée (avant: sd(z_cent))
    new_var <- cor.f * x1_std + sqrt(1 - cor.f^2) * z_std
  }

  if (strict) {
    new_var <- (new_var - mean(new_var)) / safe_sd(new_var)  # <- sécurisé
    new_var <- new_var * sdi + meani
  } else {
    new_var <- new_var * safe_sd(z) + mean(z)                # <- sécurisé
  }
  new_var
}