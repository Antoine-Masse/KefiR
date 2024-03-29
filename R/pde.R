#' Vary two parameters in a recipe and model the ideal recipe.
#'
#' @param x x values of first parameter
#' @param y y values of second parameter
#' @param z evaluation of the 6 recipes
#' @param xlim manual x limitation
#' @param ylim manual y limitation
#' @param zlim manual z limitation
#' @param dim resolution of the modelisation
#' @param xlab x title
#' @param ylab y title
#' @param zlab z title
#' @param main main
#' @param col color of th modelisation
#' @param alpha transparency
#' @param pch cf. plot function
#' @param col.pt color of points
#' @param cex.pt length of points
#' @param mode Two graphical modes (1 or 2)
#'
#' @return : By varying 2 parameters, we can obtain different recipes that we will evaluate. A 3D modeling will thus allow to anticipate what would have been the ideal parameterization to obtain the best recipe.
#' @import plot3D
#' @import rgl
#' @rawNamespace import(plotly, except = groups)
#' @import plot3Drgl
#' @import plot3D
#' @export
#'
#' @examples
#' #Example 1
#' x1   <- c(0,2,3,4,5,6)   # sugar
#' y1   <- c(0,2,1,0,4,4)   # salt
#' hedo <- c(0,3,3,2,0,0)   # hedonique taste
#' pde(x1,y1,hedo,xlab="Sucre",ylab="Sel",zlab="Note hedonique",main="Receipe",dim=30,mode=1)
#' #Example 2
#' x1   <- c(0,2,3,4,4,0,1,3,2,0,2)  ; # sucre
#' y1   <- c(0,2,1,0,4,4,3,3,0,2,4)  ;  # sel
#' hedo <- c(0,3,3,0,0,0,3,5,0,0,0)  ; # note hédonique obtenue en goûtant le produit
#' pde(x1,y1,hedo,xlab="Sucre",mode=2,pch=c(10:16),alpha=0.5,
#'     ylab="Sel",zlab="Note hédonique",main="Receipe",dim=50)
pde <- function(x,y,z,xlim=c(),ylim=c(),zlim=c(),dim=45,xlab="",ylab="",
                zlab="",main="3D plot",col=c("blue","cyan","yellow","red","#990000"),alpha=0.1,pch=16,
                col.pt="blue",cex.pt=6,mode=1) {
  # Made by Antoine Masse & Julien Bousquet
  if (length(xlim) == 0) {xlim <- range(x)}
  if (length(ylim) == 0) {ylim <- range(y)}
  if (length(zlim) == 0) {zlim <- range(z)}
  # Debogage
  BUG<-0
  if(cor.test(x,y)$p.value<0.05) {BUG<-1;print("x and y should not be correlated!\n") }
  if((length(x)!=length(y)) | (length(x)!=length(z)) | (length(unique(paste(x,y)))<6) | (length(unique(x)) < 3) | (length(unique(y)) < 3)){
    cat("It can't work.\n")
    if((length(x)!=length(y)) | (length(x)!=length(z))) {print("Vector length are not identical.\n") }
    if (length(unique(paste(x,y)))<6) {print("There is no minimum of 6 different x and y combinations.") }
    if ((length(unique(x))) < 3) {print("Beware x does not vary enough!")}
    if ((length(unique(y))) < 3) {print("Beware y does not vary enough!")}
  }else if (BUG==0) {
    lm(formula = z ~ I(x^2) + I(y^2) +I(x*y) + x + y ) -> lmpoly
    cat( "Equation de la nappe dans lmpoly : \n" )
    print(lmpoly)
    if(length(lmpoly$coefficients[is.na(lmpoly$coefficients)==TRUE])>0){
      cat("WARNING: The x and y parameters given are incompatible with a 3D regression.\nPDE() is stopped.")
      plot(x,y,main="Distribution of your values",sub="These values are not dispersed irregularly enough to allow a 3D regression.",cex=(z-min(z))/max(z-min(z))*4+1,col="blue",pch=16)
    } else {
      # Fabrication de la grille de la nappe
      x.grid <- seq(xlim[1],xlim[2],length.out=dim)
      y.grid <- seq(ylim[1],ylim[2],length.out=dim)
      xy<-expand.grid(x=x.grid,y=y.grid)
      # Calcul des ordonnées des points de la nappe
      z_nappe<-with(xy,
                    lmpoly$coefficients[2]*x^2+
                      lmpoly$coefficients[3]*y^2+
                      lmpoly$coefficients[4]*x*y+
                      lmpoly$coefficients[5]*x+
                      lmpoly$coefficients[6]*y+
                      lmpoly$coefficients[1]
      )
      z.mat<- matrix(z_nappe,nrow=dim, ncol=dim,byrow=F)
      max <- max(z_nappe)
      xy.max  =c()
      xy.max$x <- xy[,1][z_nappe==max]
      xy.max$y <- xy[,2][z_nappe==max]
      cat("Coordinates of the top :\n")
      print(xy.max)
      colfunc<-colorRampPalette(col)
      colors <- (colfunc(dim))
      #col2hex <- function(col, alpha=1) rgb(t(col2rgb(col)), alpha=alpha*255, maxColorValue=255)
      #colors <- col2hex(colors,alpha=alpha)
      colors_tp <- colfunc(length(z)*2)
      pas = (1/(length(z)*2-1))
      colorscale <- cbind(seq(0, 1, by=pas),colors_tp)
      if (mode==1) {
        print("mode 1")
        font <- list(    family = "Courier New, monospace",    size = 12,    color = "#7f7f7f" )
        xlabel <- list(    title = paste(xlab," (x)"),    titlefont = font )
        ylabel <- list(    title = paste(ylab," (y)"),    titlefont = font )
        zlabel <- list(    title = paste(zlab," (z)"),    titlefont = font )
        scene = list(    xaxis = xlabel ,    yaxis = ylabel ,    zaxis = zlabel )
        plot.pde <- plot_ly(x = x.grid, y = y.grid, z = t(z.mat), colorscale = colorscale) %>%
          layout(title = main, scene = scene)  %>% add_surface() %>%
          add_trace(x = x, y = y, z = z,name="Data brutes",
                    marker = list(
                      size = cex.pt,
                      color = col.pt,
                      line = list( color = 'rgb(231, 99, 250)', width = 2 )
                    )
          )
        #plot_ly(x = x, y = y, z = z, size = I(3))
        return(plot.pde)
        # Autres sources : https://plot.ly/r/line-and-scatter/
      }
      if (mode == 2) {
        scatter3D(x=xy.max$x,y=xy.max$y,z=max,
                  col='red',    pch=16,
                  cex=1.5,    surf=NULL,
                  xlim=xlim,    ylim=ylim,
                  zlim=zlim,    main =main,
                  xlab=xlab,     ylab=ylab,
                  zlab=zlab    )
        segments3D(x0=xy.max$x,y0=xy.max$y,z0=max,x1=xy.max$x,y1=xy.max$y,z1=0,lty='dashed',col='red',add=TRUE)
        segments3D(x0=xy.max$x,y0=xy.max$y,z0=0,x1=xy.max$x,y1=0,z1=0,lty='dashed',col='red',add=TRUE)
        segments3D(x0=xy.max$x,y0=xy.max$y,z0=0, x1=0,y1=xy.max$y,z1=0, lty='dashed',col='red',add=TRUE)
        text3D(0,xy.max$y,0,paste('ymax=',round(xy.max$y,2)),col='red',add=TRUE)
        text3D(xy.max$x,0,0,paste('ymax=',round(xy.max$x,2)),col='red',add=TRUE)
        scatter3D(x, y, z,pch=pch ,cex=1,
                  #sphere.     #size=5     #surface.alpha=0.5,    #revolution=300,
                  # colvar = NULL,     #col='blue',
                  col=colors,
                  surf=list(
                    x=x.grid,y=y.grid,z=z.mat,
                    facets = NA, col=colors,
                    #       colvar=z.mat,
                    alpha=0.1
                  ),
                  xlab="", ylab="",zlab="",add=TRUE)
        scatter3D(x, y, z,
                  pch=pch ,cex=1,col="black",
                  #xlab="",ylab="",zlab="",
                  add=TRUE
        )
        Unaccent <- function(text) {
          text <- gsub("Ç","C",text); text <- gsub("ç","c",text); text <- gsub("é","e",text); text <- gsub("È","e",text);
          text <- gsub("î","i",text); text <- gsub("ï","i",text); text <- gsub("è","e",text); text <- gsub("ê","e",text);
          text <- gsub(" ","_",text); text <- gsub("-","_",text);   text <- gsub("['`^~\"]", " ", text)
          text <- iconv(text, to="ASCII//TRANSLIT//IGNORE");   text <- gsub("['`^~\"]", "", text)
          return(text)
        }
        main <- Unaccent(main);xlab <- Unaccent(xlab);ylab <- Unaccent(ylab);zlab <- Unaccent(zlab);
        scatter3Drgl(x=xy.max$x,y=xy.max$y,z=max,
                     col='red',    pch=16,
                     cex=1.5,    surf=NULL,
                     xlim=xlim,    ylim=ylim,
                     zlim=zlim,    main =main,
                     xlab=xlab,     ylab=ylab,
                     zlab=zlab    )
        segments3Drgl(x0=xy.max$x,y0=xy.max$y,z0=max,x1=xy.max$x,y1=xy.max$y,z1=0,lty='dashed',col='red',add=TRUE)
        segments3Drgl(x0=xy.max$x,y0=xy.max$y,z0=0,x1=xy.max$x,y1=0,z1=0,lty='dashed',col='red',add=TRUE)
        segments3Drgl(x0=xy.max$x,y0=xy.max$y,z0=0, x1=0,y1=xy.max$y,z1=0, lty='dashed',col='red',add=TRUE)
        text3Drgl(0,xy.max$y,0, paste('ymax=',round(xy.max$y,2)),col='red',add=TRUE)
        text3Drgl(xy.max$x,0,0,paste('ymax=',round(xy.max$x,2)),col='red',add=TRUE)
        persp3Drgl(x=x.grid,y=y.grid,z=z.mat,col=colors,alpha=alpha ,add=T)
        scatter3Drgl(x, y, z,
                     pch=pch ,cex=1,col="black",
                     #xlab="",ylab="",zlab="",
                     add=TRUE
        )

      }
    }}}
