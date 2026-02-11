.exit <- function(ang = NULL, fr = NULL, verbose = TRUE, return = TRUE) {
  if (is.null(ang) || length(ang) == 1) {
    ang <- fr
  }
  if (verbose || return) {
    # Stop avec un message bilingue, sans répétition
    stop(.msg(ang, fr), call. = FALSE) 
  } else if (!verbose && !return) {
	signalCondition(simpleCondition("forced_exit"))
    return(1) # Mode forcé pour les boucles automatisées
  }
}