#' Plot data.frame and titles
#'
#' To help print graphics in M.test()
#' @author Julien Bousquet (2021)
#' @param x abscisse position of print
#' @param y ordinate position of print
#' @param info the information to print, a data.frame or a character string
#' @param title a character string title printed in bald font. If only title is given, print in bald and italic
device.print <- function(x, y, info=NULL, title=NULL){
  police <- police.title <- "monospace"
  linespace <- 2.25
  font <- 1 # 1 =plain, 2=bold, 4=bold+italic
  if(is.character(title)) {
    text(x,y, title, family=police.title,   adj=c(0.5-0.5*is.null(info),1), font=2+2*is.null(info))
    y <- y+linespace 
  }

  if(is.vector(info)){
     text(x, y, info,  family=police,   adj=c(0.5,1))
     y <- y+linespace 
  }else if(is.data.frame(info)){
        #Print dataframe info  into the current graphic device
         text(x,y, paste(capture.output(print(info, digits=4, row.names=FALSE)), collapse='\n'),  family=police,   adj=c(0.5,1))
         y <- y+linespace*(nrow(info)+1)
  }
 return(c(x,y))
}