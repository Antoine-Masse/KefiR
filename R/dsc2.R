#' Check the robustness of a suggested dsc() prediction of X for objectives in Y by bootstrapping
#'
#' @param data A data.frame with X(s) and Y(s).
#' @param reg A linear model or a list of linear models.
#' @param dsc A dataframe line containing the X's to be predicted (ideally a dsc() output).
#' @param iter Number of iterations in the scalable approach (should ideally be much greater than the popupulation (pop) of settings.
#' @param plot If TRUE, displays interactive parallel coordinates (plot_ly) to identify the best possible settings.
#' @param retun If TRUE, return the data.frame of values predicted by bootstrapping.
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
#' # Bootstrap on the best line with dsc2()
#' output2 <- dsc2(mtcars,reg,dsc=output[1,],plot=TRUE)
dsc2 <- function(data,reg,dsc,iter=500,plot=TRUE,return=FALSE) {
	# Mise au format list()
	if (is(reg)[1]=="lm") {
		temp <- reg ; reg = list()
		reg[[1]] <- temp ; rm(temp)
	} else if (is(reg)[1]!="list") {
		stop("Error! reg is not a list of linear model or a linear models.")
	}
  if (length(data)==0){stop("data is null.")}
  if (is(data)[1]!="data.frame"){stop("data is not data.frame.")}
	if (is(dsc)[1]!="data.frame"){stop("dsc is not data.frame. Try dsc() function for making it.")}
	reg2 <- list() ; j<-1 
	for (i in 1:length(reg)) {
		if (!is.null(reg[[i]])) {
			reg2[[j]] <- reg[[i]]
			j <- j+1
		}
	}
	reg <- reg2
	# Progress bar
	#pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
    #                 max = iter, # Maximum value of the progress bar
    ##                 style = 3,    # Progress bar style (also available style = 1 and style = 2)
     #                width = 50,   # Progress bar width. Defaults to getOption("width")
     #                char = "=")   # Character used to create the bar
#print("là")
	for (i in 1:length(reg)) {
		Y_predict_temp <-c()
		temp <- reg[[i]]
		temp_formula <- formula(temp)
		my_names <- names(get_all_vars(formula(temp$terms),data))
		my_Y <-my_names[1]
		ind_without_NA <- which(!is.na(data[,which(colnames(data)%in%my_Y)]))
		data2 <- data[ind_without_NA,]
		#print(my_Y)
		#print(nrow(data2))
		#print(temp_formula)
		for (j in 1:iter) {
			temp_reg <- lm(temp_formula,data=data2[sample(1:nrow(data2),nrow(data2),replace=TRUE),])
			prediction <- try(predict(temp_reg,dsc[1,]))		
			if (is(Y_predict_temp[length(Y_predict_temp)])[1]=="try-error"){
				Y_predict_temp <- c(Y_predict_temp,NA)	
			} else{Y_predict_temp <- c(Y_predict_temp,prediction)	
			}
		}
		if (i==1) {
			output <- data.frame(Y_predict_temp)
			colnames(output) <- my_Y
		} else {
			output <- cbind(output,Y_predict_temp)
			colnames(output)[i] <- my_Y
		}
		#setTxtProgressBar(pb, i)
	}
	#print("boucle")
	#return(output)
	for (i in 1:length(reg)) {
		print(c(min(data[,colnames(data)%in%colnames(output)[i]],na.rm=T),max(data[,colnames(data)%in%colnames(output)[i]],na.rm=T)))
	}
	if (plot==TRUE) {
		layout(matrix(1:ncol(output),1,ncol(output)))
		for (i in 1:length(reg)) {
			boxplot(output[,i],main=colnames(output)[i],ylim=c(min(data[,colnames(data)%in%colnames(output)[i]],na.rm=T),max(data[,colnames(data)%in%colnames(output)[i]],na.rm=T)))
		}
	}
	if (return==TRUE) {return(output)}
}