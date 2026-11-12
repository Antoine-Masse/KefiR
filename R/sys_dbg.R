.dbg <- function(ang=NULL,fr=NULL,debug=FALSE, ...) {
	if (.ang_missing(ang)) {ang <- fr}
	ang <- paste0("DEBUG: ",ang,"\n")
	fr <- paste0("DEBUG : ",fr,"\n")	
	if (debug==TRUE) {
		cat(.msg(ang,fr))
	} 
}
