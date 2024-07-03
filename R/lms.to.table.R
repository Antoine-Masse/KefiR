#' Make a table to summarize coefficients of a list of linear models
#' 
#' @author Julien Bousquet (2021)
#' @param L A list of linear model objects.
#'
#' @return A dataframe with all the terms of the different linear models and their coefficients.
#' The minus sign '-' is set for absence of the term in the model.
#'  
#' @details The different linear models must be computed separately, and then 
#' passed to the argument in a list format, using the function list().
#'
#' @examples
#' # Simple example
#' # Assuming lm1, lm2, and lm3 are linear models
#' # lm_list <- list(lm1, lm2, lm3)
#' # lms.to.table(L=lm_list)
#'
#' @export
lms.to.table <- function(L){
  DF <- data.frame(x=0)
  foreach(reg=L)%do%{
    # find terms of different lm, and place them as name in DF first column

  }
  foreach(reg=L, packages='broom')%do%{
    # write coeffs of different lm, and place them in DF on the row of their terms 
    GLANCE <- broom::glance(reg) # get main informations from reg
    adj.r.squared <- GLANCE$`adj.r.squared`
    pval <- GLANCE$`p.value`
  }
  return(DF)
}
