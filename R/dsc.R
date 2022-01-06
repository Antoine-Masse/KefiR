#' Obtain the best configuration to meet the objectives determined by one or more linear models.
#'
#' @param data A data.frame with X(s) and Y(s).
#' @param reg A linear model or a list of linear models.
#' @param Y Values that we want to achieve for the different Y predicted using the model (s).
#' @param ymin List of minimum values tolerated for the different Y.
#' @param ymax List of maximum values tolerated for the different Y.
#' @param pop Population of parameters which will cross randomly to generate better parameters.
#' @param iter Number of iterations in the scalable approach (should ideally be much greater than the popupulation (pop) of settings.
#' @param wash The maximum number of desired settings.
#' @param plot If TRUE, displays interactive parallel coordinates (plot_ly) to identify the best possible settings.
#' @param verbose If TRUE, gives information about the analysis.
#' @param save For saving the graph (html format)
#' @param file Name of the html page "xxx.html"
#'
#' @return A dataframe containing all the selected settings sorted from best (top) to worst (bottom).
#' @export
#'
#' @examples
#' data(mtcars)
#' colnames(mtcars)
#' myreg1 <- evolreg(mtcars,"mpg")
#' myreg2 <- evolreg(mtcars,"cyl")
#' reg <- list()
#' reg[[1]] <- myreg1
#' reg[[2]] <- myreg2
#' output <- dsc(mtcars,reg,Y=c(23.4,5.4),pop=400,iter=200)
#' # Aggregation of several trials
#' for (i in 1:10) {
#' 	output <- rbind(output,dsc(mtcars,reg,Y=c(23.4,5.4),plot=F))
#' } ; parco(output,"Distance")
#' # With filtration of min and max y.
#' output <- dsc(mtcars,reg,Y=c(15,5),ymin=c(14,4),ymax=c(15,6),pop=5000,iter=10000)
dsc <- function(data,reg,Y=c(),ymin=c(),ymax=c(),pop=iter/20,iter=4000,wash=pop/2,plot=T,verbose=T,save=F,file="file.html") {
	# Mise au format list()
	if (is(reg)[1]=="lm") {
		temp <- reg ; reg = list()
		reg[[1]] <- temp ; rm(temp)
	} else if (is(reg)[1]!="list") {
		stop("Error! reg is not a list of linear model or a linear models.")
	}
  if (length(data)==0){stop("data is null.")}
  if (is(data)[1]!="data.frame"){stop("data is not data.frame.")}
  if (length(Y)==0){stop("Y must have the same number of values as models in reg.")}
	# Recuperation des variables utiles
	my_colnames <- c() ; my_Y <- c()
	for (i in 1:length(reg)) {
		temp <- reg[[i]]
		my_names <- names(get_all_vars(formula(temp$terms),data))
		# model.frame(formula, data = NULL, …)
		my_Y <- c(my_Y,my_names[1])
		my_colnames <- union(my_colnames,my_names)
	}
	if (verbose == TRUE) {
		print("Names of useful variables : ")
		print(setdiff(my_names,my_Y))
		print("Names of useful Y : ")
		print(my_Y)
	}
	data <- data[,which(colnames(data)%in%my_colnames)]
	# Captation des proprietes des Y pour centration-reduction
	if (length(my_Y)==1) {
		crY <- data.frame(c(mean(data[,which(colnames(data)%in%my_Y)]),sd(data[,which(colnames(data)%in%my_Y)])))
	} else {
		crY <- apply(data[,which(colnames(data)%in%my_Y)],2,mean)
		crY <- rbind(crY,apply(data[,which(colnames(data)%in%my_Y)],2,sd))
	}
	# Generation des x parentaux
	n_x <- length(my_colnames) - length(my_Y)
	x <- apply(data,2,function(x){rbeta(pop,0.5/n_x,0.5/n_x)*(max(x)-min(x))+min(x)}) # min(x) max(x)
	x <- data.frame(x)
	# Generation of Y parentaux
	for (i in 1:length(reg)) {
		temp <- reg[[i]]
		prediction<-predict(temp,x)
		x[,which(colnames(x)%in%my_Y[i])] <- prediction
	}
	# Centration-Reduction des Y (Objectifs)
	for (colonne in 1:ncol(crY)) {
		Y[colonne]<- ((Y[colonne]-crY[1,colonne])/crY[2,colonne])
	}
	if (verbose == TRUE) {
		print("center-reduced Y : ")
		print(Y)}
	# Centration-Reduction des Y issues des x puis calcul d'une distance quadratique par rapport à chaque objectif Y
	xY <- data.frame(x[,which(colnames(x)%in%my_Y)] )
	for (colonne in 1:ncol(xY)) {
		xY[,colonne] <- (xY[,colonne]-crY[1,colonne])/crY[2,colonne]
	}
	# Distance quadratique des objectifs
	diff_parents <- xY
	for (colonne in 1:ncol(diff_parents)) {
		diff_parents[,colonne] <- (diff_parents[,colonne]-Y[colonne])^2
	}
	diff_parents <- apply(diff_parents,1,mean)
	if (verbose==TRUE) {
		print("Before - best combination")
		print(x[which(diff_parents==min(diff_parents)),])
	}
	# Progress bar
	pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
	# Evolutiv approach
	for (i in 1:iter) {
		ind <- sample(1:pop,2)
		x_enfant <- apply(x[ind,],2,mean)
		for (j in 1:length(reg)) {
			temp <- reg[[j]]
			prediction<-predict(temp,data.frame(t(unlist(x_enfant))))
			x_enfant[my_Y[j]] <- prediction
		}
		# Distance quadratiques des objectifs des enfants
		# Centration-Reduction puis calcul d'une distance quadratique pour chaque Y
		xY_enfant <- x_enfant[which(colnames(x)%in%my_Y)]
		for (colonne in 1:length(xY_enfant)) {
			xY_enfant[colonne] <- (xY_enfant[colonne]-crY[1,colonne])/crY[2,colonne]
		}
		diff_enfant <- xY_enfant
		for (colonne in 1:length(diff_enfant)) {
			diff_enfant[colonne] <- ( diff_enfant[colonne]-Y[colonne] )^2
		}
		diff_enfant <- mean(diff_enfant)
		if ((diff_enfant < max(diff_parents[ind]))&(diff_enfant > min(diff_parents[ind]))) {
			ind_p <- which(diff_parents[ind]==max(diff_parents[ind]))
			x[ind[ind_p],] <- x_enfant
			diff_parents[ind[ind_p]]  <- diff_enfant
		} else if (diff_enfant < min(diff_parents[ind])) {
			ind_p <- which(diff_parents==max(diff_parents))
			x[ind_p,] <- x_enfant
			diff_parents[ind_p]  <- diff_enfant
		}
		setTxtProgressBar(pb, i)
	}
	Distance <- diff_parents
	x <- cbind(x,Distance)
	tri <- order(diff_parents)
	tri <- x[tri,]
	if ((length(ymin)==length(my_Y))&(length(ymax)==length(my_Y))) {
		indice <- 1
		for (colonne in which(colnames(tri)%in%my_Y)) {
			filtre <- ((tri[,colonne]>ymin[indice])&(tri[,colonne]<ymax[indice]))
			if (any(filtre) == FALSE) {
				stop("ymin and ymax too restrictive.")
			} else {
				tri <- tri[filtre,]
			}
			indice <- indice+1
		}
	} else if ((length(ymin)>0)|(length(ymax)>0)){
		stop("Error ! ymin or ymax do not correspond to the number of Y")
	}
	if (wash<nrow(tri)){
		tri <- tri[1:wash,]
	}
	if (plot==TRUE) {
		parco(tri, Y = "Distance", X = c(),save=save,file=file)
	}
	if (verbose==TRUE) {
		print("")
		print("After")
		print(tri[1,])
	}
	return(tri)
}

