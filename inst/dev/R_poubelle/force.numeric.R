#' Force numerical variable to be numeric in a dataframe
#'
#' @param data a dataframe.
#'
#' @return The dataframe with numerical variables where it is possible.
#' @export
#'
#' @examples
#' tablo <- data.frame("A"=as.character(rnorm(20)),"B"=c(rnorm(5),"C3",rnorm(4)),"C"=rep("A",10))
#' apply(tablo,2,is.numeric)
#' force.numeric(tablo)-> tablo2 ; apply(tablo2,2,is.numeric)
#' 
force.numeric <- function(data) {
	temp1 <- as.numeric(as.character(data[,1]))
	if (all(is.na(temp1))) {temp1 <- data[,1]}
	temp2 <- as.numeric(as.character(data[,2]))
	if (all(is.na(temp2))) {temp2 <- data[,2]}
	dt <- data.frame(temp1,temp2)
	colnames(dt) <- colnames(data)[1:2]
	if (ncol(data)>2) {
  		for (i in 3:ncol(data)){
    			#if (length(unique(data[,i]))>1) {
    			temp <- as.numeric(as.character(data[,i]))
    			if (all(is.na(temp))) {temp <- data[,i]}
    			dt <- data.frame(dt,temp)
    			colnames(dt)[ncol(dt)] <- colnames(data)[i]
    		#}
 		 }
	}
  return(dt)
}