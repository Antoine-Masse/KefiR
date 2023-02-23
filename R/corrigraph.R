#' igraph of correlated variables global or in relation to y
#'
#' @param data a data.frame
#' @param colX a vector of indices or variables to follow. We will only keep the variables that are connected to them on 1 or more levels (level parameter).
#' @param colY a vector of indices or variables to predict. To force the correlogram to display only the variables correlated to a selection of Y.
#' @param type "x" or "y". To force the display in correlogram mode (colX, type = "x") or in prediction mode (colY, type = "y").
#' @param alpha the maximum permissible p-value for the display
#' @param exclude the minimum threshold of displayed correlations - or a vector of threshold in this order : c(cor,mu,prop)
#' @param ampli coefficient of amplification of vertices
#' @param return if return=T, returns the correlation matrix of significant correlation.
#' @param wash automatically eliminates variables using differnts methods when there are too many variables (method = NA, stn (signal-to-noise ratio), sum, length).
#' @param multi to ignore multiple regressions and control only single regressions.
#' @param mu to display the effect on median/mean identified by m.test().
#' @param prop to display the dependencies between categorical variables identified by GTest().
#' @param layout to choose the network organization method - choose "fr", "circle", "kk" or "3d".
#' @param cluster to make automatic clustering of variables or not.
#' @param verbose to see the comments.
#' @param NAfreq from 0 to 1. NA part allowed in the variables. 1 by default (100% of NA tolerate).
#' @param NAcat TRUE or FALSE. Requires recognition of missing data as categories.
#' @param level to be used with colY. Number of variable layers allowed (minimum 2, default 5).
#' @param evolreg TRUE or FALSE. Not yet available. Allows you to use the evolreg function to improve the predictive ability (R squared) for the variables specified in colY.
#'
#' @return Correlation graph network (igraph) of the variables of a data.frame. Pay attention to the possible presence of non-numeric variables or missing data. Grouping of correlated variables: the vertices (circles) correspond to the variables. The more a variable is connected, the larger it appears. The color of the lines reflects the nature of the correlation (positive or negative).  The size of the lines is the value of the correlation from 0 to 1. All these correlations are significant (pval < 0.01). The coloured groupings reflect families of inter-correlated variables. BLUE: positive correlation - RED: negative correlation
#' @return When mu is TRUE or prop : we see the connexion with mean effect (orange) and G (~chisq) effect (pink)
#' @return The size of orange edge and pink edge depend of p-values (-1*log10(p-value)/10) of kruskal.test() and GTest().
#' @return
#' @return  When indicating Y's in colY, the correlogram will identify the correlated X's, then the remaining X's correlated to these X's, and so on.
#' @return X's not related to these Y's are excluded.
#' @return The blue always displays the positive correlations and the red, negative correlations. When the display is green, it means that the predictive (~correlation) capacity of the variable can be reinforced by adding a 2nd variable in a multiple regression model (interaction X1+X2, X1*X2 or X1+X1:X2) better than X1 or X2 alone.
#' @return Correlations between X or Y of the same level are neglected.
#' @return The color of the vertices makes it possible to identify the correlated variables alone in a significant way (blue: positive, red: negative, purple: positive or negative depending on the Y).
#' @return The values displayed to the right of the Ys (colY) correspond to the maximum predictive capacity of these Ys by one or two variables.
#'
#' @import utils
#' @rawNamespace import(igraph, except = c(decompose,spectrum))
#' @importFrom stats na.omit
#' @importFrom stats cor
#' @importFrom DescTools GTest
#'
#' @export
#'
#' @examples
#' # Example 1
#' data(swiss)
#' corrigraph(swiss)
#' # Example 2
#' data(airquality)
#' corrigraph(airquality,layout="3d")
#' # Example 3
#' data(airquality)
#' corrigraph(airquality,c("Ozone","Wind"),type="y")
#' # Example 4
#' data(iris)
#' corrigraph(iris,mu=TRUE)
#' # Example 5
#' require(MASS) ; data(Aids2)
#' corrigraph(Aids2 ,prop=TRUE,mu=TRUE,exclude=c(0.3,0.3,0))
#' # Example 6
#' data(airquality)
#' corrigraph(airquality,c("Ozone","Wind"),type="x")
corrigraph <- function(data,colY=c(),colX=c(),type="x",alpha=0.05,exclude=c(0,0,0), ampli=4,return=FALSE,wash="stn",multi=TRUE,
					   mu=FALSE,prop=FALSE,layout="fr",cluster=TRUE,verbose=FALSE,NAfreq=1,NAcat=TRUE,level=2,evolreg=FALSE) {
  # Fonction réalisée par Antoine Massé
  # Ctrl Alt Shift R
  # Version 04
  # December 2021
  #igraph::graph_from_adjacency_matrix
  databrut<-data # saving
  # Control 1 - is.numeric ?
  which(sapply(data, is.numeric)) -> temp_id_num ; temp_id_factor <- c()
  if ((mu==FALSE)&(prop==FALSE)){
	if (length(temp_id_num)<2) {stop("Error! Not enough numerical variable. Try mu=TRUE or prop=TRUE ?\n")}
	data <- data[,temp_id_num]
	which(sapply(data, is.numeric)) -> temp_id_num
	if (length(temp_id_num)!=length(1:ncol(data))) {warning("Warning! Presence of non-numeric variables that cannot be taken into account.\n")}
  } else {
	temp_id_factor <- setdiff(1:ncol(data),temp_id_num)
  }
  # temp_id_factor = categorical values identified
  # Control 2 - var is NULL ?
  which(sapply(data.frame(data[,temp_id_num]), var, na.rm=T)!=0)-> temp_var ; temp_var	<- temp_id_num[temp_var]
  if (length(temp_var) != length(temp_id_num)) {warning("Warning! Some variables have a null variance and cannot be taken into account.\n")}
  data <- data[,union(temp_var,temp_id_factor)] ; temp_id_num <- which(sapply(data, is.numeric)) ; temp_id_factor  <- setdiff(1:ncol(data),temp_id_num)
  # Control 3 - is.na ?
  if (any(is.na(data))==TRUE) {if (NAfreq == 1) {warning("Warning ! Presence of missing values. Please check NAfreq argument.\n")}}
  if (ncol(data)>50) {warning("Warning! The calculation time increases exponentially with the number of variables (do not exceed 50).\n")}
  ##################################################
  # Prediction mode (colY)
  ##################################################
  if (length(colY)>0) {
	#print("mode 2")
	# If modelization of Y
	if (is.numeric(colY)) {
	  colnames(databrut)[colY] -> colY
	}
	which(colnames(data)%in%colY) -> indY
	if (length(colY)!=length(indY)){stop("ERROR! some Y's provided cannot be analyzed.\n")}
	if (length(indY)==0) {stop()}
	if (type=="x"){
		indX <- indY
		colX <- colY
	} else {
		type<-"y"
		if((mu==TRUE)|(prop==TRUE)){stop("Error! mu and prop can't be used width colY.\n")}
	}
}
if ((type=="y")&(length(colY)>0)){
	indices_all <- 1:ncol(data)
	# Matrice de correlation
	matrix(rep(0,ncol(data)^2),ncol(data),ncol(data))->mymat
	rownames(mymat) <- colnames(data)
	colnames(mymat) <- colnames(data)
	#mymat2 <- mymat
	sum_mymat2 <- ncol(data)^2
	levels <- rep(0,ncol(mymat))
	borders <- rep("black",ncol(mymat))
	vertices <- rep("white",ncol(mymat))
	level_count <- 1
	levels[indY] <- level_count
	warning<-0
	##############################
	# Function for calculation of different BUC for different interaction X1 & X2 for Y
	##############################
	check_mdl <- function(x) {
	  data_temp <- data.frame("data[,i]"=data[,i],"data[,j]"=data[,j],"x"=x)
	  total <- nrow(data_temp)
	  data_temp <- na.omit(data_temp)
	  total_NA <- total - nrow(data_temp)
	  seuil_NA <- total_NA/total
	  nb_subset <- nrow(data_temp)
	  #cat("seuil_NA",NAfreq,"\n")
	  #print(seuil_NA)
	  if ((min(apply(data_temp,2,var))!=0)&(nb_subset>3)&(seuil_NA<=NAfreq)) {
		#cat("seuil_NA",seuil_NA,"\n")
		regA <- lm(data_temp[,1]~data_temp[,2]+data_temp[,3])
		regB <- lm(data_temp[,1]~data_temp[,2]+data_temp[,2]:data_temp[,3])
		regC <- lm(data_temp[,1]~data_temp[,2]*data_temp[,3])
		synthese <- data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA)
		oldw <- getOption("warn")
		options(warn = -1)
		pvalsA <- c(
		  pf(summary(regA)$fstatistic[1],summary(regA)$fstatistic[2],summary(regA)$fstatistic[3],lower.tail=FALSE),
		  summary(regA)[[4]][,4]	)
		pvalsB <- c(
		  pf(summary(regB)$fstatistic[1],summary(regB)$fstatistic[2],summary(regB)$fstatistic[3],lower.tail=FALSE),
		  summary(regB)[[4]][,4]	)
		pvalsC <- c(
		  pf(summary(regC)$fstatistic[1],summary(regC)$fstatistic[2],summary(regC)$fstatistic[3],lower.tail=FALSE),
		  summary(regC)[[4]][,4]	)
		options(warn = oldw)
		# Comparison with k alone
		data_temp2 <- data.frame(data[,i],x)
		data_temp2 <- na.omit(data_temp2)
		nb_subset2 <- nrow(data_temp2)
		if ((min(apply(data_temp2,2,var))!=0)&(nb_subset2>3)) {
		  temp2 <- cor.test(temp_tab[,1],temp_tab[,2],na.rm=T)
		  temp2 <- ifelse(temp2$p.value<= alpha,abs(temp2$estimate),0)
		  options(warn=-1)
		  if (temp2 != 0 ) {
			temp2 <- mean(c(summary(lm(temp_tab[,1]~temp_tab[,2]))$adj.r.squared,summary(lm(temp_tab[,2]~temp_tab[,1]))$adj.r.squared))
		  }
		  if (nrow(na.omit(data.frame(data[,i],data[,j]))) < length(data[,i])&temp>0) {
				temp <- mean(c(summary(lm(data[,i]~data[,j]))$adj.r.squared,summary(lm(data[,j]~data[,i]))$adj.r.squared))
			}
		  if (temp2>0) {
			reg2 <- lm(data_temp2[,1]~data_temp2[,2],data_temp2)
			oldw <- getOption("warn")
			options(warn = -1)
			pvals2 <- try(c(
			  pf(summary(reg2)$fstatistic[1],summary(reg2)$fstatistic[2],summary(reg2)$fstatistic[3],lower.tail=FALSE),
			  summary(reg2)[[4]][,4]	))
			options(warn = oldw)
			if ((any(is.na(pvals2)))|(length(pvals2)<=1)){synthese2 <- data.frame(R2_2 = 0,BIC_2 = myBIC)
			} else {synthese2 <- data.frame(R2_2 = abs(summary(reg2)$adj.r.squared),BIC_2 = BIC(reg2))}
		  } else {synthese2 <- data.frame(R2_2 = 0,BIC_2 = myBIC)}
		}else {synthese2 <- data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA)}
		# Agglomeration of datas
		if (any(is.na(pvalsA))){	rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else if (max(pvalsA)>alpha){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else {
		  tablo <- data.frame(R2 = abs(summary(regA)$adj.r.squared),BIC = BIC(regA))
		  synthese <- rbind(synthese,cbind(tablo,synthese2))}
		if (any(is.na(pvalsB))){	rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else if (max(pvalsB)>alpha){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else {
		  tablo <- data.frame(R2 = abs(summary(regB)$adj.r.squared),BIC = BIC(regB))
		  synthese <- rbind(synthese,cbind(tablo,synthese2))
		}
		if (any(is.na(pvalsC))){	rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else if (max(pvalsC)>alpha){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
		} else {
		  tablo <- data.frame(R2 = abs(summary(regC)$adj.r.squared),BIC = BIC(regC))
		  synthese <- rbind(synthese,cbind(tablo,synthese2))}
	  } else {return(data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))}
	  return(synthese)
	}
	mymat_interconnect <- mymat
	while ((sum_mymat2  != sum(mymat))|(level_count <= level)) { # Tant que la matrice evolue à chaque cycle...
	  sum_mymat2 <- sum(mymat)
	  indices <- c()
	  for (i in indY) { # Modelization of Y (i = Y)
		if (evolreg==TRUE) {
			level <- 2
			X_retenus <- evolreg(data,colnames(data)[i],NAfreq=NAfreq, interaction=FALSE, multix=FALSE,multidiv=FALSE)
			# ATTENTION DE PROTEGER SI PAS DE R2 EN SORTIE
			if ((length(X_retenus$R2)>0)&(X_retenus$R2!=0)){ # Coloration selon lien neg ou pos
				temp <- X_retenus$R2
				for (j in X_retenus$indices) {
						indices <- c(indices,j) ; mymat[i,j] <- temp ; mymat[j,i] <- temp
						#if (abs(temp) > abs(temp_save)){
						mymat_interconnect[i,j]<-2;mymat_interconnect[j,i]<-2
						#} else if (temp_save<0) {
						#  mymat_interconnect[i,j]<- -1;mymat_interconnect[j,i]<- -1;
						#} else if (temp_save>0){
						#  mymat_interconnect[i,j]<-1;mymat_interconnect[j,i]<- 1;
						#}
				}
			}
		} else {
		indices_X <- setdiff(indices_all,indY)
		if (length(indices_X)==0){break}
		for (j in indices_X) { # En fonction d'un X
			temp_tab <- data.frame(data[,i],data[,j])
		  	total <- nrow(temp_tab)
			na.omit(temp_tab) -> temp_tab
			total_NA <- total - nrow(temp_tab)
			seuil_NA <- total_NA/total
			nb_subset <- nrow(temp_tab)
			if ((min(apply(temp_tab,2,var)) != 0)&(nb_subset>3)&(seuil_NA<=NAfreq)) {
				temp <- cor.test(data[,i],data[,j],na.rm=T)
				temp_save <- ifelse(temp$p.value<= alpha,temp$estimate,0) # saving cor neg or pos
				temp <- ifelse(temp$p.value<= alpha,abs(temp$estimate),0) # saving abs(cor)
				options(warn=-1)
				if (temp != 0) {
					temp <- mean(c(summary(lm(data[,i]~data[,j]))$adj.r.squared,summary(lm(data[,j]~data[,i]))$adj.r.squared))
				}
				options(warn=0)
				if (temp_save > 0) {  # coloring vertices if significant correlation neg or pos // Y
				  if (borders[j]=="purple") {
				  } else if (borders[j]=="black") {
					borders[j]<-"blue" ; vertices[j] <- "cyan"
				  } else if (borders[j]=="red") {
					borders[j]<-"purple"  ;  vertices[j] <- "#FBC5F8"
				  }
				} else if (temp_save < 0) {
				  if (borders[j]=="purple") { # do nothing
				  } else if (borders[j]=="black") {
					borders[j]<-"red" ; vertices[j] <- "#F9A09E"
				  } else if (borders[j]=="blue") {
					borders[j]<-"purple";  vertices[j] <- "#FBC5F8"
				  }
				}
				# Regresiv analysis (only for cor < 0.9) , for acceleration
				k <- setdiff(indices_all,c(indY,j))
				if ((temp<0.90)&(length(k)>0)&(multi==TRUE)) {
				  reg <- lm(data[,i]~data[,j])
				  myBIC <- BIC(reg)
				  temps <- data.frame(R2=temp,BIC=myBIC,R2_2=0,BIC_2=myBIC)
				  dataij <- data.frame(data[,i],data[,j])
				  datak <- data.frame(data[,k])
				  temps <- rbind(temps,do.call("rbind", apply(datak,2,check_mdl)))
				  if (all(is.na(temps))) {# do nothing
				  } else {
					temps <- na.omit(temps)
					temps <- temps[union(1,which(temps$R2>temps$R2_2)),]
					temp_tp <- temps$R2[which((temps$BIC == min(temps$BIC))&(temps$BIC <= min(temps$BIC_2)))]
					if (length(temp_tp )>=1){temp <- max(temp,temp_tp)}
				  }
				}
			}else{
				temp<-0
				if(warning==0){
				  warning <- 1
				  warning("Warning! Failure to account for missing data generated zero variances on some variables that had to be ignored.\n")}
			}
			if (temp!=0){ # Coloration selon lien neg ou pos
				indices <- c(indices,j) ; mymat[i,j] <- temp ; mymat[j,i] <- temp
				if (abs(temp) > abs(temp_save)){
				  mymat_interconnect[i,j]<-2;mymat_interconnect[j,i]<-2
				} else if (temp_save<0) {
				  mymat_interconnect[i,j]<- -1;mymat_interconnect[j,i]<- -1;
				} else if (temp_save>0){
				  mymat_interconnect[i,j]<-1;mymat_interconnect[j,i]<- 1;
				}
			}
		}
		}
	  }
	  indices_all <- setdiff(indices_all,indY)
	  indY <- unique(indices)
	  level_count <- level_count+1
	  levels[indY] <- level_count
	  if (level_count==level) {break}
	}
	indices <- which(levels!=0)
	mymat <- mymat[indices,indices]
	if (sum(mymat)==0) {stop("Error! No x reliebable to Y\n")}
	mymat_interconnect <- mymat_interconnect[indices,indices]
	levels <- levels[indices]
	borders <- borders[indices]
	vertices <- vertices[indices]
	layout <- matrix(rep(0,length(levels)*2),length(levels),2)
	maximum <- max(table(levels))
	largeur = 4*maximum
	for (i in unique(levels)) {
	  rev(which(levels==i)) -> pos
	  layout[pos,1] <- (max(levels)-i)*(largeur/(max(levels)))+1 # Positionnement horiz
	  step <- maximum/(length(pos)+1)
	  layout[pos,2] <- seq(maximum/(length(pos)+1),to = step*length(pos),by=step)  # Positionnement vertic
	}
	net_interconnect <- graph_from_adjacency_matrix(mymat_interconnect, weighted=TRUE,mode="directed")
	net <- graph_from_adjacency_matrix(mymat, weighted=TRUE,mode="directed")
	net <- simplify(net, remove.multiple = T, remove.loops = TRUE) # elaguer les liens redondants
	net_interconnect  <- delete.edges(net_interconnect , E(net_interconnect)[ abs(E(net)$weight) < exclude[1] ])
	net <-               delete.edges(net, E(net)[ abs(weight) < exclude[1] ])
	E(net)$colour <- ifelse(E(net_interconnect)$weight<0,"red",ifelse(E(net_interconnect)$weight==2,"green","blue"))
	E(net)$weight <- abs(E(net)$weight)
	plot.igraph(net,layout=layout,vertex.size=8*ampli,vertex.color=vertices,vertex.frame.color=borders,
				edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
	# Display max correlation
	tabulation <-  table(levels)
	seq_max <- layout[which(levels==as.numeric(names(tabulation)[which(tabulation==max(tabulation))])),2]
	seq <- layout[layout[,1]==max(layout[,1]),2]
	seq <- seq -min(seq_max)
	seq <- seq/(max(seq_max)-min(seq_max))
	seq <- (seq*2)-1
	correl  <- apply(mymat,2,max)
	correl <- round(correl[which(levels==1)],2)
	ifelse(correl>0.9,">0.9",correl)-> correl
	text(rep(1.3,length(seq_max)),seq,
		 correl)
	if (return==TRUE) {return(mymat)}
  } else {
##################################################
# Correlation matrice for basal corrigraph()
##################################################
	#print("mode1")
	#	 Matrice de correlation - Initialization
	matrix(rep(0,ncol(data)^2),ncol(data),ncol(data))->mymat
	rownames(mymat) <- colnames(data)
	colnames(mymat) <- colnames(data)
	matrix(rep(1,ncol(data)^2),ncol(data),ncol(data))->mymat2
	rownames(mymat2) <- colnames(data)
	colnames(mymat2) <- colnames(data)
	warning<-0
	mymat_interconnect <- mymat
	vertices <- rep("white",ncol(mymat)) ; vertices[temp_id_factor]<-"grey"
	if ((prop==TRUE)|(mu==TRUE)) {
	  # Factor interactions
	  if (verbose==TRUE) {
		cat("Non-numeric variables.\n")
		print(colnames(data)[temp_id_factor])
	  }
		if (length(temp_id_factor)>5){
			barre <- txtProgressBar(min=0,max=length(temp_id_factor),width=50)
			cat("\nprop calculation\n")
			cat(paste(rep("=",50),collapse=""),"\n")
		}
	  # For each categorical
	  for (i in temp_id_factor) {
		if ((length(temp_id_factor)>5)&(prop==TRUE)){setTxtProgressBar(barre,i)}
		if (min(table(data[,i]))<3) {
		  temp_id_factor <- setdiff(temp_id_factor,i)
		  next
		}
		#########################
		#	prop=T
		#########################
		if ((prop==TRUE)&(length(temp_id_factor)>1)) {
		  if (verbose==TRUE) {cat("\nAnalysis of the cat:") ; cat(colnames(data)[i])}
		  for (j in setdiff(temp_id_factor,i)) {
			if (verbose==TRUE) {
				cat("\n\t vs categorical : ") ; cat(colnames(data)[j])
			}

			data_temp <- data.frame(data[,i],data[,j])
			#print(colnames(table(data_temp)))
			#print(head(data_temp))
			#print(which(is.na(data_temp[,1])))
			#print(which(is.na(data_temp[,2])))


			if (identical(which(is.na(data[,i])),which(is.na(data[,j])))==TRUE) {
				#cat("na omit")
				data_temp <- na.omit(data_temp)
				#print(table(data_temp))
			}
			if (NAcat==TRUE) {
				for (k in 1:2){
					data_temp[which(is.na(data_temp[,k])),k]<-"MANQUANTES"
				}
			}
			# Pourquoi plante qd mu activé seul sur big_fusion

			# NA are automatically eliminated by table()
			tablo <- table(data_temp)

			if ((ncol(tablo)>=2) & (nrow(tablo)>=2)) {
			  #cat("Test entre",i,"&",j,"\n")
			  #print(tablo)
			  if (min(tablo)<5) {
				#cat("Effectif insuffisant\n")
				if (quantile(tablo,probs=c(0.2))<5) {
				  if (verbose==TRUE) {
				    cat("\nNot enough individuals per category.\n")
				  }
				  next
				}
			  }
			  oldw <- getOption("warn")
			  options(warn = -1)
			  pvals <- try(GTest(tablo)$p.value,silent=TRUE)
			  options(warn = oldw)
			  #print("GTEST")
			  #print(pvals)
			  #cat("Test entre",colnames(data)[i],"&",colnames(data)[j],"pval",pvals,"\n")
			  if (pvals <= alpha) {
				#print(pvals)
				log10(1/pvals)/10 -> temp ; ifelse(temp>1 ,1 ,temp)-> temp
				mymat[i,j] <- temp ; mymat[j,i] <- temp
				mymat2[i,j] <- pvals ; mymat2[j,i] <- pvals
				mymat_interconnect[i,j] <- 2 ; mymat_interconnect[j,i] <- 2
			  }
			} else {next}
		  }
		}
	  }
	  if ((length(temp_id_factor)>5)&(prop==TRUE)){close(barre)}
	  #print(colnames(data)[temp_id_factor])
	  #########################
	  #	mu=TRUE
	  #########################
	  if ((mu == TRUE)&(length(temp_id_factor)>0)) {
		if (length(temp_id_factor)>5){
			barre <- txtProgressBar(min=0,max=length(temp_id_factor),width=50)
			cat("\nmu calculation\n")
			cat(paste(rep("=",50),collapse=""),"\n")
		}
		for (i in temp_id_factor) {
			if ((length(temp_id_factor)>5)&(mu==TRUE)){
				#print("Barre de progression")
				setTxtProgressBar(barre,i)
			}
		  if (verbose==TRUE) {cat("\nAnalysis of the cat:") ; print(colnames(data)[i])}
		  for (j in temp_id_num) {
			if (verbose==TRUE) {cat("\n\tvs numerical: ") ; print(colnames(data)[j])}
			data_temp <- data.frame(data[,i],data[,j])
			if (NAcat==TRUE) {
				data_temp[which(is.na(data_temp[,1])),1]<-"MANQUANTES"
			}
			pvals <- m.test(data_temp[,2],data_temp[,1],alpha=alpha,return=FALSE,verbose=FALSE,plot=FALSE,boot=FALSE)
			if (pvals < alpha) {
			  log10(1/pvals)/10 -> temp ; ifelse(temp>1 ,1 ,temp)-> temp
			  mymat[i,j] <- temp ; mymat[j,i] <- temp
			  mymat2[i,j] <- pvals ; mymat2[j,i] <- pvals
			  mymat_interconnect[i,j] <- 3 ; mymat_interconnect[j,i] <- 3
			}
		  }
		}
		if (length(temp_id_factor)>5){close(barre)}
	  }
	}

	###############################
	# Numerical correlations
	###############################
	if (length(temp_id_num)>1) {
	  barre <- txtProgressBar(min=0,max=length(temp_id_num),width=50)
	  cat("\ncor calculation\n")
	  cat(paste(rep("=",50),collapse=""),"\n")
	  for (i in temp_id_num) {
		setTxtProgressBar(barre,i)
		for (j in setdiff(temp_id_num,i)) {
			temp_tab <- data.frame(data[,i],data[,j])
			total <- nrow(temp_tab)
			na.omit(temp_tab) -> temp_tab
			total_NA <- total - nrow(temp_tab)
			seuil_NA <- total_NA/total
			# L'ajout de NAfreq fait planter prop = T et mu = T?
			# Ajoutons que si mu=T (prop = T) apparait d'office si NAfreq < 0
			# Etablir le lien entre NAfreq et prop et mu...






		  if ((var(temp_tab[,1])!=0)&(var(temp_tab[,2])!=0)&(nrow(temp_tab)>2)&(seuil_NA<=NAfreq)){
			tempA <- cor.test(data[,i],data[,j],na.rm=T)
			#temp <- ifelse(temp$p.value<= alpha,temp$estimate,0) # Working with correlation
			pvals  <- ifelse(tempA$p.value<= alpha,tempA$p.value,1) # Working with p-value
			log10(1/pvals)/10 -> temp ; ifelse(temp>1 ,1 ,temp)-> temp
			ifelse(tempA$estimate<0,-temp,temp)->temp
		  }else{
			temp<-0
			pvals <- 1
			if(warning==0){warning <- 1}
		  }
		  mymat[i,j] <- temp ; mymat[j,i] <- temp
		  mymat2[i,j] <- pvals ; mymat2[j,i] <- pvals
		  mymat_interconnect[i,j] <- temp ; mymat_interconnect[j,i] <- temp
		}
	  }
	  close(barre)
	}
	cat("\n")
	if (warning==1) {warning("Warning! Failure to account for missing data generated zero variances on some variables that had to be ignored.\n")}
	#print(mymat)
	#print(mymat_interconnect)
	#########################
	#	colX
	#########################
	# Modelization arround X
	if ((type=="x")&(length(colY)>0)){
		#if (length(colX)!=length(indX)){stop("ERROR! some X's provided cannot be analyzed.\n")}
		#if (level<=1){stop("ERROR! Level must exceed 1.")
		#} else {
			if (level<2) {level<-2}
			for (lev in 1:(level-1)) {
				#print(mymat[indX,])
				indX_temp <- union(indX, which(apply(data.frame(mymat[indX,]),1,sum)!=0))
				if (length(indX_temp)>length(indX)){
					indX <- unique(sort(indX_temp))
				} else {
					break
				}
			}
		#}
		mymat <- mymat[indX,indX]
		mymat2 <- mymat2[indX,indX]
		mymat_interconnect <- mymat_interconnect[indX,indX]

	}
	#########################
	#	Wash
	#########################
	if (verbose ==TRUE) {cat("	Reduction of the number of variables (wash).\n")}
	pas = (ncol(mymat)-50)/20
	if ((is.na(wash)!=TRUE)&(wash!=FALSE)) {
	  indices <- apply(abs(mymat),2,sum)
	  indices <- which(indices>0)
	  mymat <- mymat[indices,indices]
	  mymat2 <- mymat2[indices,indices]
	  mymat_interconnect <- mymat_interconnect[indices,indices]
	  vertices <- vertices[indices]
	  while(ncol(mymat)>50){
		#mymat <- ifelse(abs(mymat)< exclude,0,mymat)
		if (wash=="length"){
			if(verbose==T) {cat("Cleaning of the variables (wash parameter) in length mode.\n")}
		  taille <- apply(abs(mymat),2,function(x){return(length(x[x>0]))})
		  indices <- 1:ncol(mymat)
		  indices [order(taille ,decreasing=T)]
		  indices <- indices[1:(ncol(mymat)-pas)]
		  mymat <- mymat[indices,indices]
		  mymat2 <- mymat2[indices,indices]
		  mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
		  indices <- apply(abs(mymat),2,sum)
		  indices <- which(indices>1)
		  mymat <- mymat[indices,indices]
		  mymat2 <- mymat2[indices,indices]
		  mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
		} else if (wash=="sum"){
			if(verbose==T) {cat("Cleaning of the variables (wash parameter) in sum mode.\n")}
		  somme <- apply(abs(mymat),2,sum)
		  indices <- 1:ncol(mymat)
		  indices [order(somme,decreasing=T)]
		  indices <- indices[1:(ncol(mymat)-pas)]
		  mymat <- mymat[indices,indices]
		  mymat2 <- mymat2[indices,indices]
		  mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
		  indices <- apply(abs(mymat),2,sum)
		  indices <- which(indices>1)
		  mymat <- mymat[indices,indices]
		  mymat2 <- mymat2[indices,indices]
		  mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
		} else { #(wash=="stn"){
			if(verbose==T) {cat("Cleaning of the variables (wash parameter) in stn mode.\n")}
		  if (ncol(mymat) > 50) {
			moyennes <- apply(abs(mymat),2,mean)
			deviation <- apply(abs(mymat),2,sd)
			snp <- moyennes/deviation
			indices <- 1:ncol(mymat)
			indices [order(snp,decreasing=T)]
			indices <- indices[1:(ncol(mymat)-pas)]
			mymat <- mymat[indices,indices]
			mymat2 <- mymat2[indices,indices]
			mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
			indices <- apply(abs(mymat),2,sum)
			indices <- which(indices>1)
			mymat <- mymat[indices,indices]
			mymat2 <- mymat2[indices,indices]
			mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
		  }
		}
	  }
	}
	if (ncol(mymat)==0) {stop("No interaction identified.\n")}
	if ((mu==FALSE)&(prop==FALSE)) {
	  net <- graph_from_adjacency_matrix(mymat, weighted=T,mode="lower")
	  net <- simplify(net, remove.multiple = T, remove.loops = TRUE) # élaguer les liens redondants
	  net <- delete.edges(net, E(net)[ abs(weight) < exclude[1] ])
	  E(net)$colour <- ifelse(E(net)$weight<0,"red","blue")
	  E(net)$weight <- abs(E(net)$weight)
	  if (layout=="circle") {l <- layout_in_circle(net)
	  } else if (layout=="kk"){l <- layout_with_kk(net)
	  } else  {l <- layout_with_fr(net)}
	  if ((layout=="3d") |(layout=="3D")) {
		l <- layout_with_fr(net,dim=3)
		rglplot( net, layout = l, vertex.size = sqrt(betweenness(net))*ampli/4+10,
				 vertex.color = vertices, edge.arrow.size = 0,
				 arrow.mode = 0, edge.width = (abs(E(net)$weight) * 2)^2,
				 edge.color = E(net)$colour,label.dist=sqrt(betweenness(net))*ampli/4+10)
	  }else if (ncol(mymat)<50) {
		clp <- cluster_optimal(net)
		class(clp)
		plot(clp, net, layout = l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color="yellow",
			 edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^2, edge.color =E(net)$colour)
	  } else {
		plot(net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
			 edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^2, edge.color =E(net)$colour)
	  }
	} else {
	  ######################
	  # mu = TRUE or prop = TRUE
	  ######################
	  net <- graph_from_adjacency_matrix(mymat, weighted=T,mode="directed")
	  #net <- simplify(net, remove.multiple = TRUE, remove.loops = TRUE) # elaguer les liens redondants
	  net_interconnect <- graph_from_adjacency_matrix(mymat_interconnect, weighted=T,mode="directed")
	  #net_interconnect <- simplify(net_interconnect, remove.multiple = TRUE, remove.loops = TRUE) # elaguer les liens redondants
	  ##############################
	  #		Nettoyage différentielle des liens en fonction de la p-value
	  ##############################
	  if (length(exclude)==3) {
		net2 <- net ; net_interconnect2 <- net_interconnect
		net_interconnect  <- delete.edges(net_interconnect , E(net_interconnect)[ abs(E(net2)$weight) < exclude[1] & abs(E(net_interconnect2)$weight) < 2 ])
		net <- 	delete.edges(net, E(net)[ abs(E(net2)$weight) < exclude[1] & abs(E(net_interconnect2)$weight) < 2])
		net2 <- net ; net_interconnect2 <- net_interconnect
		net_interconnect  <- delete.edges(net_interconnect , E(net_interconnect)[ abs(E(net2)$weight) < exclude[2] & abs(E(net_interconnect2)$weight) == 3 ])
		net <- 	delete.edges(net, E(net)[ abs(E(net2)$weight) < exclude[2] & abs(E(net_interconnect2)$weight) ==3])
		net2 <- net ; net_interconnect2 <- net_interconnect
		net_interconnect  <- delete.edges(net_interconnect , E(net_interconnect)[ abs(E(net2)$weight) < exclude[3] & abs(E(net_interconnect2)$weight) == 2 ])
		net <- 	delete.edges(net, E(net)[ (abs(E(net2)$weight) < exclude[3]) & (abs(E(net_interconnect2)$weight) ==2)])
	  } else {
		net_interconnect  <- delete.edges(net_interconnect , E(net_interconnect)[ abs(E(net)$weight) < exclude[1]])
		net <- 	delete.edges(net, E(net)[ abs(E(net)$weight) < exclude[1] ])
	  }
	  # Coloration des liens
	  E(net)$weight <- abs(E(net)$weight)
	  #print(E(net_interconnect)$weight)
	  E(net)$colour <- ifelse(E(net_interconnect)$weight<0,"red",
							  ifelse(E(net_interconnect)$weight==2,"#FF99BB",
									 ifelse(E(net_interconnect)$weight==3,"orange","blue")))
	  # Mise en place de la layout
	  if (layout=="circle") {l <- layout_in_circle(net)
	  } else if (layout=="kk"){l <- layout_with_kk(net)
	  } else  {l <- layout_with_fr(net)}
	  # Basculer en affichage 3D ou non
	  if ((layout=="3d") |(layout=="3D")) {
		l <- layout_with_fr(net,dim=3)
		rglplot( net, layout = l, vertex.size = sqrt(betweenness(net))*ampli/4+10,
				 vertex.color = vertices, edge.arrow.size = 0,
				 arrow.mode = 0, edge.width = (abs(E(net)$weight) * 2)^2,
				 edge.color = E(net)$colour,label.dist=sqrt(betweenness(net))*ampli/4+10)
	  } else {
		# Si AFFICHAGE Non-3D : clustering sur les p-values
		#plot(net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
		#	 edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^2, edge.color =E(net)$colour)
		# Matrice de distances
		as.dist(mymat2) -> mydi
		if ((length(mydi)>1)&(cluster==TRUE)) {
			# Clustering
			myclust <- hclust(mydi, method="ward.D2")
			inertie <- sort(myclust$height, decreasing = TRUE)
			inertie <- inertie[2:(length(inertie)-1)]-inertie[3:length(inertie)]
			which(inertie==max(inertie,na.rm=T))+1->n_cluster
			#print(n_cluster)
			mytree = cutree(myclust, k=n_cluster)
			# ATTENTION NE MARCHE QUE SI WASH = TRUE, faire en sorte que mymat2 suiva mymat
			# Mettre les communautés en option pour ne pas masquer la couleur des vertices catego vs num
			# INACTIVE OU PAS LE CLUSTERING POUR QU'ON PUISSE FAIRE LA DISTINCTION ENTRE VARIABLE CATEGORIELLES ET QUANTITATIVES.
			##method1: communities object type clone
			members <- as.double(mytree)
			nam <- as.character(names(mytree))
			comms <- list(membership=members, vcount=vcount(net), names=nam, algorithm="by.hand")
			class(comms) <- "communities"
			#try(modularity(net, membership(comms) ),silent=TRUE)
			plot(comms,net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
				 edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^2, edge.color =E(net)$colour)
		} else {
			plot(net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
				 edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^2, edge.color =E(net)$colour)
		}
	  }
	}
	if (return==TRUE){return(mymat2)}
  }
}
