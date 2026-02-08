#' Determining groups based on the results of a pairwise function
#'
#' @param result output results of the pairwise(), pairwise.t.test() or pairwise.wilcox.test() function
#' @param control name of the category that will be used as a control to establish differences with '*', '**' and '***'.
#' @param alpha threshold of p-value to establish groups of averages or medians of type "a", "ab", "b" ...
#' @param debug logical, if TRUE enables debug mode for troubleshooting.
#'
#' @importFrom stringr str_sort
#'
#' @return catego() can be used to establish groups with significantly different mean/median values relative to each other or to a reference control. It exploits the results of the functions pairwise.t.test(), pairwise.wilcox.test() and pairwise(type = "median" or "mean").
#' @export
#'
#' @examples
#' data(iris)
#' output <- pairwise.t.test(iris[,1],iris$Species)
#' catego(output)
#' catego(output,control="setosa")
catego <- function(result,control=c(),alpha=0.05,debug=FALSE) {
  .dbg("catego() - start.","Execution de catego().",debug=debug)

  # Fallback pour str_sort si stringr n'est pas disponible
  .str_sort <- if (requireNamespace("stringr", quietly = TRUE)) {
    stringr::str_sort
  } else {
    function(x, ...) sort(x)
  }

  result$p.value -> mymat
  #.dbg(NULL, paste0("Matrice des p-valeurs :\n", paste(capture.output(print(signif(mymat,2))), collapse = "\n")), debug = debug)
  categories <- c(colnames(mymat),rownames(mymat)[nrow(mymat)])
  #print("a") ; print(categories)
  #print("b") ; print(control)
  # ERREUR A CORRIGER : no match.
  which(categories%in%control)-> ind_control
  #which(control%in%levels(categories))-> ind_control
  if (length(ind_control)==1) {
    c(1,mymat[,1]) -> cats
    ifelse(cats<=0.001,"***",ifelse(cats<=0.01,"**",ifelse(cats<=0.05,"*","")))->cats
  } else {
  if (length(categories)<=26) {
	possibility <- letters[seq( from = 1, to = length(categories))]
  } else {possibility <- letters[seq( from = 1, to = 26)]}
  
    start = 0 ; pos=1 ; cats <- c("a") ; change <- 0
    for (i in 1:ncol(mymat)) {
      #cat("Passage",start,"\n")
      for (j in i:nrow(mymat)) {
        pvals <- mymat[j,i]
        # Gérer les NA : traiter comme non-significatif (conserver le groupe actuel)
        if (is.na(pvals)) {
          if (start == 0) {
            cats <- c(cats, possibility[pos])
          }
          next
        }
        if (start==0) {

          if (pvals > alpha){
            #print("non-signif")
            cats <- c(cats,possibility[pos])
          } else {
            #print("signif")
            cats <- c(cats,possibility[pos+1])
            change <- 1
          }
        } else {
          if (pvals <= alpha){
            #print(length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]])))
            if (length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]]))>0){
              #print("signif avec changement de cat.")
              # a on avait b vs b, l'autre b devient c
              cats[j+1]<- possibility[pos+1]
              change <- 1

            } else {
              #print("signif")
              next
            }
          } else {
            if (length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]]))>0) {
              #print("non-signif")
              next
            } else {
			  #.dbg(NULL, paste0("Appel a .str_sort avec: ", c(cats[i], cats[j + 1])), debug = debug)
              #print("signif mais sensé etre de cat differente (situation intermediaire)")
			  if (.str_sort(c(cats[i],cats[j+1]))[1] == cats[i]) {
				cats[i] <- paste0(.str_sort(c(cats[i],cats[j+1])),collapse="")
			  } else {
				cats[i+1] <- paste0(.str_sort(c(cats[i],cats[j+1])),collapse="")
			  }
            }
          }
        }
#.dbg(NULL, paste("Groupes intermédiaires en j :", paste(cats, collapse = ",")), debug = debug)
      }
#.dbg(NULL, paste("Groupes intermédiaires en i :", paste(cats, collapse = ",")), debug = debug)	  
      #print(cats)
      #print("########")
      start <- start+1
      if (change == 1) {pos <- pos+1 ; change <- 0}
    }
  }
  synth <- list()
  groups <- cbind(categories,groups=cats) ; rownames(groups) <- rep("",nrow(groups))
  synth$groups <- groups
  synth$p.value <- result$p.value
  .dbg("catego() - End","Fin de catego().",debug=debug)
  return(synth)
}
