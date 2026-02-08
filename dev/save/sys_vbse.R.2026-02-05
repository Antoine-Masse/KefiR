.vbse <- function(ang=NULL,fr=NULL, k=NULL,cpt="on",verbose=FALSE, code=FALSE) {
	# Protection contre code=NULL
	if (is.null(code)) code <- FALSE

	# En mode code=TRUE, on n'affiche rien mais on incrémente k si nécessaire
	if (isTRUE(code)) {
		if (cpt=="on") {
			k <- k + 1
		}
		return(k)
	}

	# Code normal : afficher les messages si verbose==TRUE
	if (isTRUE(verbose)) {
		#if (cpt=="on") {k <- k + 1}
		if ((is.null(ang)) | (length(ang)==1)) {ang <- fr}
		if (cpt=="on") {
			k <- k + 1
			ang <- paste0(k,") ",ang,"\n")
			fr <- paste0(k,") ",fr,"\n")
			cat(.msg(ang,fr))
			return(k)
		} else {
			ang <- paste0("\t",ang,"\n")
			fr <- paste0("\t",fr,"\n")
			cat(.msg(ang,fr))
		}
	}
	return(k)
}
