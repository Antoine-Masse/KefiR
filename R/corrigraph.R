#' igraph of correlated variables
#'
#' @param data a data.frame
#' @param pval the maximum permissible p-value for the display
#' @param exclude the minimum threshold of displayed correlations
#' @param ampli coefficient of amplification of vertices
#' @param return if return=T, returns the correlation matrix of significant correlations
#' @param wash automatically eliminates variables using the signal-to-noise ratio when there are too many variables
#'
#' @return Correlation graph network (igraph) of the variables of a data.frame. Pay attention to the possible presence of non-numeric variables or missing data. Grouping of correlated variables: the vertices (circles) correspond to the variables. The more a variable is connected, the larger it appears. The color of the lines reflects the nature of the correlation (positive or negative).  The size of the lines is the value of the correlation from 0 to 1. All these correlations are significant (pval < 0.01). The coloured groupings reflect families of inter-correlated variables. BLUE: positive correlation - RED: negative correlation
#' @import igraph
#' @importFrom stats na.omit
#' @importFrom stats cor
#' @export
#'
#' @examples
#' # Example 1
#' data(swiss)
#' corrigraph(swiss)
#' # Example 2
#' data(airquality)
#' corrigraph(airquality)
corrigraph <- function(data,pval=0.01,exclude=0.3, ampli=4,return=FALSE,wash=TRUE) {
  # Fonction réalisée par Antoine Massé
  # Ctrl Alt Shift R
  # Version 01
  # Janvier 2021
  # Control 1 - is.numeric ?
  sapply(data, is.numeric) -> temp_id
  data <- data[,temp_id]
  if (any(temp_id==FALSE)) {cat("Warning ! Presence of non-numeric variables that cannot be taken into account.\n")}
  # Control 2 - var is NULL ?
  sapply(data, var, na.rm=T)->temp_var
  which(temp_var!=0) -> temp_id
  if (length(temp_id) != ncol(data)) {cat("Warning ! Some variables have a null variance and cannot be taken into account.\n")}
  data <- data[,temp_id]
  # Control 3 - is.na ?
  if (any(is.na(data))==TRUE) {cat("Warning ! Presence of missing values.\n")}
  if (ncol(data)>50) {cat("warning : The calculation time increases exponentially with the number of variables (do not exceed 50).\n")}
  # Matrice de correlation
  matrix(rep(0,ncol(data)^2),ncol(data),ncol(data))->mymat
  rownames(mymat) <- colnames(data)
  colnames(mymat) <- colnames(data)
  warning<-0
  for (i in 1: ncol(data)) {
    for (j in 1: ncol(data)) {
      temp_tab <- data.frame(data[,i],data[,j])
      na.omit(temp_tab) -> temp_tab
      if ((var(temp_tab[,1])!=0)&(var(temp_tab[,2])!=0)&(nrow(temp_tab)>2)){
        temp <- cor.test(data[,i],data[,j],na.rm=T)
        temp <- ifelse(temp$p.value<= pval,temp$estimate,0)
      }else{
        temp<-0
        if(warning==0){
          warning <- 1
          cat("Warning : Failure to account for missing data generated zero variances on some variables that had to be ignored.\n")}
      }
      mymat[i,j] <- temp
    }
  }
  pas = (ncol(mymat)-50)/20
  if (wash==TRUE) {
    while(ncol(mymat)>50){
      #mymat <- ifelse(abs(mymat)< exclude,0,mymat)
      indices <- apply(abs(mymat),2,sum)
      indices <- which(indices>1)
      mymat <- mymat[indices,indices]
      if (ncol(mymat) > 50) {
        moyennes <- apply(abs(mymat),2,mean)
        deviation <- apply(abs(mymat),2,sd)
        snp <- moyennes/deviation
        indices <- 1:ncol(mymat)
        indices [order(snp,decreasing=T)]
        indices <- indices[1:(ncol(mymat)-pas)]
        mymat <- mymat[indices,indices]
        indices <- apply(abs(mymat),2,sum)
        indices <- which(indices>1)
        mymat <- mymat[indices,indices]
      }}
  }
  #return(mymat)
  net <- graph_from_adjacency_matrix(mymat, weighted=T,mode="lower")
  net <- simplify(net, remove.multiple = T, remove.loops = TRUE) # élaguer les liens redondants
  net <- delete.edges(net, E(net)[ abs(weight) < exclude ])
  E(net)$colour <- ifelse(E(net)$weight<0,"red","blue")
  E(net)$weight <- abs(E(net)$weight)
  if (ncol(mymat)<50) {
    clp <- cluster_optimal(net)
    class(clp)
    l <- layout_with_fr(net)
    plot(clp, net, layout = l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color="yellow",
         edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
  } else {
    plot(net,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color="yellow",
         edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
  }
  if (return==TRUE){return(mymat)}
}

