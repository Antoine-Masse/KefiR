.dbg <- function(ang=NULL,fr=NULL,debug=FALSE, ...) {
	if ((is.null(ang)) | (length(ang)==1)) {ang <- fr}
	ang <- paste0("DEBUG: ",ang,"\n")
	fr <- paste0("DEBUG : ",fr,"\n")	
	if (debug==TRUE) {
		cat(.msg(ang,fr))
	} 
}
