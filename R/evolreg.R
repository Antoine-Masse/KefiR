#' Identify the best linear, logistic or mixed regression model using an evolutionary approach.
#'
#' @param data a dataframe.
#' @param Y the y to predict.
#' @param X vector of variables whose presence we want to force in the model.
#' @param pval 0 to 1. If there are too many variables and the argument wash=TRUE, use this p-value threshold to eliminate the variables whose effect is too insignificant (Risk of eliminating the variables that will have an effect once transformed or in interaction).
#' @param nvar Maximum number of variables in the model. A default value is proposed according to the number of individuals.
#' @param iter Number of iterations.
#' @param NAfreq from 0 to 1. NA part allowed in the variables. 1 by default (100% of NA tolerate).
#' @param interaction FALSE or TRUE. To allow interactions between variables.
#' @param multix FALSE or TRUE. To allow variable variants (log, exp, polynomial, ^2).
#' @param multidiv FALSE or TRUE. To allow the synthesis of variables combining the ratio of one variable divided by another.
#' @param nbind Number of simulated individuals in the population of models to be crossed in an evolutionary approach. A default value is proposed according to the number of variables.
#' @param family "lm", "logical" or "lmer". Type of regression
#' @param wash TRUE or FALSE.To select the best variables when there are too many.
#' @param plot To visualize the evolution of the R2 of the models obtained after each crossing.
#' @param verbose To display a summary of the intermediate models.
#' @param fast Paramètre qui stoppe les itérations lorsque le R carré des modèles n'évolue plus de façon significative.
#'
#' @importFrom e1071 naiveBayes
#'
#' @return A strongest possible regression model chosen from the available variables.
#' @export
#'
#' @examples
#' data(mtcars)
#' evolreg(mtcars,"mpg",plot=TRUE)
evolreg <- function(data,Y, X=c(),pval=0.05, nvar = 0,
                    iter = 4000,multix=TRUE, interaction = TRUE, multidiv=TRUE,nbind=c(),wash=TRUE,NAfreq=1,
                    family = "lm",plot=FALSE,verbose=TRUE,fast=TRUE) {
  # Identification des variables d'intérêt
  #	print("Id var interest")
  data -> data_save
  # Substraction of Y and washing NA
  if ((is.numeric(Y)==TRUE)&(Y>0)&(Y<=ncol(data))){ind_Y <- Y ; Y <- colnames(data)[Y]
  } else {ind_Y <- which(colnames(data)%in%Y)}
  if (length(ind_Y)==0){stop("Impossible to find Y.")} 
  data <- data[which(!is.na(data[,ind_Y])),]
  #
  data_glob <- dvar(data,Y,pval=0.05,X=X,NAfreq=NAfreq,wash=wash,interaction=interaction,verbose=FALSE)
  data <- data_glob$data
  dt <- data_glob$variables
  if (verbose==TRUE) {print(head(dt))}
  # Une détection devrait être proposée par défaut de la famille de Y
  #
  # Ainsi qu'une élimination optionnelle des var quantitatives
  ###########
  #	Changer la valeur de pval proportionnellement à l'augmentation des variables par interaction
  ###########
  #print(head(data))
  #return(data)
  #	print("End Interaction")
  # Si interaction==TRUE, remettre en place un nettoyage dvar pour ne pas se noyer...
  # Corriger la fonction dvar pour qu'elle sorte automatiquement le jeu de données nettoyées avec un verbose sur le tableau des variables.
  # Inclure le wash directement dans dans le dvar : on peut ainsi faire un dvar avant et après interaction.
  # Si interaction =TRUE , autorise un nvar de plus.
  # Supprimer les parents en double si apparait un doublon
  # Calcule du nombre de variables optimal
  if (nvar==0) {
    reste <- 5 ; individus <- nrow(data)
    #			print(reste)
    #			print(individus)
    #			print(nvar)
    while (reste >= 5) {
      nvar <- nvar+1
      reste <- individus/2^nvar
      #			print(reste)
      #			print(individus)
      #			print(nvar)
    }
    nvar <- nvar-1
  }
  if (nvar > (ncol(data)-1)) {
    nvar <- ncol(data)-1
  }
  if (verbose==TRUE) {
    print("Nombre de variable maximal conseillé :")
    print(nvar)
    #print(colnames(data))
  }
  ###########################################
  #
  ###########################################
  #
  bornage_global <- 0 ; bornage2 <- 1 ; bornage3 <- 0 ; bornage4 <- 1000
  globalBIC <- 100000 ; resultat_global <- 0 ; distance_global <- 1000000
  global_save <- c() ; BIC_save <- c() ; resultat_min <- 0 ; super_reg <- NA
  my_i_model <- 0 ; formul <- NA ; formule <- NA ; global <- 0
  '%notin%' <- Negate('%in%')
  individuals <- nrow(dt)*(ceiling(40/log(nrow(dt))))^2;
  individuals <- ifelse(individuals>20000,20000,individuals)
  if(length(nbind)>0){if (nbind>individuals){individuals <- nbind}}
  nb <- c()
  #print("la")
  for (i in 1:nvar) {
    nb <- c(nb,  rep(i,choose(ncol(data),i)))
  }
  ###################################################
  #			FAIRE LES PARENTS
  ###################################################
  my_i_model <- c() ; 
  if (verbose==TRUE) {print(paste("Making the ",individuals," parents."))}
  for (indiv in 1:individuals) {
    nb_temp <- sample(nb,1) ; nvar_temp <- 0 ; formul <- c() ; my_i <- c()
    #print("new indiv")
    if (length(X)>0) {
      ind <- which(dt$variables%in%X)
      nvar_temp <- nvar_temp + dt$weight[ind]
      formul <- c(formul,dt$variables[ind])
      my_i <- c(my_i,ind)
    }
	breaker <- 0
    while(nvar_temp <= nb_temp) { # nb_temp : nombre de var à atteindre, nvar_temp : nombre réel
      #print(nb_temp)
      ind <- sample(nrow(dt),1)
      nvar_temp <- nvar_temp + dt$weight[ind]
      #print(nvar_temp)
      if (nvar_temp<= nb_temp){
        formul <- c(formul,dt$variables[ind])
        my_i <- c(my_i,ind)
		if (nvar_temp == nb_temp) {break}
      } else if ((nvar_temp>nb_temp)&(length(formul)==0)){
        nvar_temp <- 0
		breaker <- breaker+1
		if ((breaker>3)&(nvar_temp<=(nb_temp+1))) {
			formul <- c(formul,dt$variables[ind])
			my_i <- c(my_i,ind)		
			break
		}
      }
	  if (breaker>5) {
		formul <- c(formul,dt$variables[ind])
        my_i <- c(my_i,ind)  
		break}
    }
	ssave <- formul
	formul_temp <- paste0(Y,"~",paste0(formul,collapse='+'))
    formul <- try(formula(formul_temp))
#if (is(formul)[1]=="try-error"){
#	print("La formule a planté")
#			print(formul_temp)
#			print(ssave)
#			return(ssave)
#}
	var_a_croiser <- names(get_all_vars(formul,data))
	data_a_croiser <- data[,colnames(data)%in%var_a_croiser]
	if ((nrow(na.omit(data_a_croiser))==0)|(min(apply(data_a_croiser,2,function(x){length(unique(x,na.rm=T))}))<2)) {
		next
	}
    if (family == "lm") {
      reg <- try(lm(formula=formul,data=data))
	  if (is(reg)[1]=="try-error"){
			print(dt$variables[my_i])
			print(formul)
			print(colnames(data))
			#return(data)
		}
	  if (is(reg)[1]=="lm") {
		  global <- summary(reg)$adj.r.squared
		  bic_temp <- BIC(reg)
	  } else {
		  global <- 0
		  bic_temp <- 100000
		}
		#global_save <- c(global_save,global)
		#BIC_save <- c(BIC_save,bic_temp)
    } else if (family=="logit") {
      # DES COMPILATIONS DE DATA3 à REMPLACER PAR DATA...
      data3[,1] <- as.factor(data3[,1]) # Mettre en facteur
      model <- naiveBayes(formul,data=data3)
      prediction <- try(predict(object=model,newdata=data3))
      tempdf <- data.frame("test" = data3[,which(colnames(data3)%in%Y)],prediction)
      tempdf <- table(tempdf)
      tempdf <- round(prop.table(tempdf), 2)
      global <- sum(tempdf[1,1],tempdf[2,2])
    }
	if (length(global)==0) {global<-0}
	if (length(bic_temp)==0) {bic_temp<-100000}
    parent <- list()
    parent[[1]] <- my_i ; parent[[2]] <- global ; parent[[3]] <- bic_temp
    if (indiv == 1) {parents <- list()}
    parents[[indiv]] <- parent
    if ((length(X) == 0) | (length(intersect(which(dt$variables%in%X),my_i))>0)) {
      if ((global >= bornage_global)&(bic_temp<globalBIC)) {
          if (family == "lm") {
            pvals <- summary(reg)[[4]][,4]
          } else if (family=="logit") {
            # Encore du data3 à nettoyer
            reg <- glm(formul,data=data3,family=binomial(logit))
            reg2 <- drop1(reg, test="F")
            pvals <- unlist(reg2[5])[-1]
          }
		 if (any(is.na(reg$coefficients))==FALSE) {
		  if (any(is.na(pvals))==FALSE) {
            if (any(pvals>0.05)==FALSE) {
              if (identical(sort(my_i),sort(my_i_model))==FALSE){
			    resultat <- c()
				for (i in 1:250) {
              		apprentissage <- sample(1:nrow(data),replace=T)
              		test <- c()
              		test <- setdiff(1:nrow(data),apprentissage)
	var_a_croiser <- names(get_all_vars(formul,data))
	data_a_croiser <- data[apprentissage,colnames(data)%in%var_a_croiser]
	if ((nrow(na.omit(data_a_croiser))==0)|(min(apply(data_a_croiser,2,function(x){length(unique(x,na.rm=T))}))<2)) {
		next
	}
              		if (family == "lm") {		
						reg <- try(lm(formul,data=data[apprentissage,]))
						#prediction <- try(predict(reg,newdata=data[test,]))
						if (is(reg)[1]=="lm") { 
							prediction <- summary(reg)$adj.r.squared
#print(paste0("prediction ",prediction))
							resultat <- c(resultat,prediction)
							if (length(resultat)>100) {
								if (quantile(resultat,probs=0.05,na.rm=T)<resultat_min){break}
							}
						} else {next}
              		} else if ((family=="logit")&(length(test)>0)&(length(unique(data3[,which(colnames(data3)%in%Y)][apprentissage]))==2)&(min(table(data3[,which(colnames(data3)%in%Y)][test]))>2)&(min(table(data3[,which(colnames(data3)%in%Y)][apprentissage]))>2)) {
						ctrl <- c()
              			for (i in 1:ncol(data)){
              				if (is.numeric(data[,i])==FALSE){
              					part_common <- sort(intersect(unique(data[apprentissage,i]),unique(data[test,i])))
              					part_common <- as.factor(part_common)
              					ctrl <- c(ctrl,identical(part_common,as.factor(sort(unique(data[test,i])))))
              				}
              			}
						if (all(ctrl==TRUE)==TRUE) {
              				reg <- glm(formul,data=data3[apprentissage,],family=binomial(logit))
              				prediction <- try(predict(reg,newdata=data3[test,]))
              				prediction <- ifelse(prediction>0.5,"oui","non")
              				tempdf <- data.frame("test" = data3[,which(colnames(data3)%in%Y)][test],prediction)
              				tempdf <- table(tempdf)
              				tempdf <- round(prop.table(tempdf ), 2)
              				if (ncol(tempdf)==2) {resultat <- c(resultat,sum(tempdf[1,1],tempdf[2,2]))
              				} else {resultat <- c(resultat,tempdf[1,1])}
              			}
					}
				}
				if (length(resultat)==0) {resultat<-0}
				if (quantile(resultat,probs=0.05,na.rm=T)>=resultat_min){
					resultat_min<-quantile(resultat,probs=0.05,na.rm=T)
					resultat_global<-median(resultat)
					# Global function after bootstrap
					if (family == "lm") {
						reg <- lm(formul,data=data)
					} else if (family=="logit") {
						reg <- glm(formul,data=data3,family=binomial(logit))
					}
					reg2 <- drop1(reg, test="F")
					super_reg <- reg ; super_reg2 <- reg2
					my_i_model <- my_i
					formule = formul
					formule <- formula(formule)
					globalBIC <- BIC(reg)
					bornage_global <- global
					if (family == "lm") {
					} else if (family=="logit") {
						modelG <- naiveBayes(form,data=data3)
						# print(formule)
						# cat("Prévision médiane : ",median(resultat),"\n")
						# cat("Prévision inférieure 95% : ",quantile(resultat,probs=0.05),"\n")
						prediction <- predict(model,data)
						tempdf <- data.frame("test" = data[,which(colnames(data)%in%Y)],prediction)
						tempdf <- table(tempdf)
						tempdf <- round(prop.table(tempdf ), 2) ;
						tempdf_save <- tempdf
						#print("Prédiction : ")
						#		    print(tempdf)
					}
					if ((verbose==TRUE)&(bornage_global == 0)&(bornage_global!=global)) {print("At least one significant model identified.")}
					if (verbose == TRUE) {
						cat("Modèle sur parents avec R2 : ",round(global,2),"\n")
						print(paste("BIC ",BIC(reg)))
						print(paste("Performance inf 95%",round(resultat_min,2)))
						print(formule)
					}
				}
              }
            }
          }
		 #}
		}
	  }
    }
  }
  ###########################################
  ###########################################
  tosupp <- which(sapply(parents, is.null))
  if (length(tosupp)>0) {parents <- parents[-tosupp]}
  ############################################################################################
  ############################################################################################
  # Making evolutiv approach
  ############################################################################################
  ############################################################################################
  if (bornage_global>0) { # Ne pas exec les croisements parentaux si pas de R² > 0
  super_reg_global <- bornage_global
  if (verbose==TRUE) {print(paste("Making evolutiv approach on ",iter," iterations."))}
  stocking <- c() ; stocking2<- c() ; pval_stop2 <- 1
  #rognage <- 400/round(log2(length(parents)),0) # Nombre de part de prélèvements...
  rognage <- 400
  prelevement <- 2 + 1 + round(length(parents)/rognage,0)
  prelevement_save <- round(((prelevement-3)/2),0) + 2
  if (prelevement_save>prelevement) {prelevement_save = 2 ; prelevement = 3}
  for (gene in 1:iter) {
	########################
	#	Prélever un échantillon parental
	########################
    #if (verbose == TRUE) {cat("Prélèvement",prelevement," et save ",prelevement_save," sur ",length(parents)," parents.\n")}
    if (length(parents)< (prelevement)) {
      break
    }
    sample(1:length(parents),prelevement) -> parents_temp
	#
    # Indices parents_temp des individus à éliminer
    ind_to_kill <- which((sapply(parents,"[[",2)[parents_temp])%in%(sort(sapply(parents,"[[",2)[parents_temp])[1:(prelevement-prelevement_save)]))
	if ((length(parents_temp)-length(ind_to_kill))<2){
		#print("Nettoyage")
		if (length(ind_to_kill) > length(parents_temp)) {ind_to_kill <- c()
		} else if ((length(parents_temp) - length(ind_to_kill))==0) {ind_to_kill <- ind_to_kill[-c(1:2)]
		} else {ind_to_kill <- ind_to_kill[-1]
		}
		#print(ind_to_kill)
	}
    ind_to_kill2 <- parents_temp[ind_to_kill]
    parents_temp <- parents_temp[-ind_to_kill]
	#print("parents_temp2")
	#print(parents_temp)
    ind_to_parents <- which(sapply(parents,"[[",2)[parents_temp]%in%sort(sapply(parents,"[[",2)[parents_temp],decreasing = TRUE)[1:2])
    #	print(ind_to_parents)
    parents_temp <- parents_temp[ind_to_parents]
	#print("parents_temp3")
	#print(parents_temp)
    dad <- sort(parents[[parents_temp[1]]][[1]])
    mom <- sort(parents[[parents_temp[2]]][[1]])
    #print(mom)
    #print(dad)
    if (identical(mom,dad)==TRUE) {
      parents <- parents[-parents_temp[1]]
      stocking <- c(stocking,0)
      if (length(stocking)>10) {
        stocking2 <- c(stocking2,median(stocking[(length(stocking)-9):length(stocking)]))
      }
      next
    }
    nb_temp <- sample(nb,1) # Material of son/daughter
    demi_nb_temp <- round(nb_temp/2,0) # Length of gamets
    #print(demi_nb_temp)
    if (sum(dt$weight[dad])>=demi_nb_temp) {
      if (demi_nb_temp>=length(dad)) {spermatozoid <- dad
      } else {sample(dad,demi_nb_temp)->spermatozoid}
    }else{spermatozoid<-dad}
    demi_nb_temp <- nb_temp - length(spermatozoid)
    if (sum(dt$weight[mom])>=demi_nb_temp) {
      if (demi_nb_temp>=length(mom)) {ovul <- mom
      } else {sample(mom,demi_nb_temp)->ovul}
    }else{ovul <- mom}
    my_i <- union(spermatozoid,ovul)
    poids <- sum(dt$weight[my_i])
    if (poids < nb_temp) {
      my_i_max <- union(mom,dad)
      my_i_max <- setdiff(my_i_max,my_i)
      if (length(my_i_max)>0) {
        if (length(my_i_max)>=(nb_temp-poids)){
          my_i <- my_i_max[1:(nb_temp-round(poids,0))]
        } else {
          my_i <- union(mom,dad)
        }
      }
    }
#print("B")
	breaker <- 0
    while((poids>nb_temp)|(poids>nvar)){
      my_i_temp<-sample(my_i)[-1]
      poids <- sum(dt$weight[my_i_temp])
      if (length(my_i_temp)>=1) {
        my_i <- my_i_temp
      } else { poids <- 0 ; breaker <- breaker+1
      }
	  if (breaker > 3) {break}
    }
    formul <- formula(paste0(Y,"~",paste0(dt$variables[my_i],collapse='+')))
	var_a_croiser <- names(get_all_vars(formul,data))
	data_a_croiser <- data[,colnames(data)%in%var_a_croiser]
	if ((nrow(na.omit(data_a_croiser))==0)|(min(apply(data_a_croiser,2,function(x){length(unique(x,na.rm=T))}))<2)) {
		next
	}
    if (family == "lm") {
      reg <- try(lm(formula=formul,data=data))
	  if (is(reg)[1]=="lm") {
		global <- summary(reg)$adj.r.squared
		bic_temp <- BIC(reg)
		if (length(bic_temp)==0) {bic_temp<- 1E6}
	  } else {global <- 0 ; bic_temp<- 1E6}
    } else if (family=="logit") {
      # Là aussi du data3 à remplacer par data
      data3[,1] <- as.factor(data3[,1]) # Mettre en facteur
      model <- naiveBayes(formul,data=data3)
      prediction <- predict(object=model,newdata=data3) # Mettre try et condition is...
      tempdf <- data.frame("test" = data3[,which(colnames(data3)%in%Y)],prediction)
      tempdf <- table(tempdf)
      tempdf <- round(prop.table(tempdf), 2)
      global <- sum(tempdf[1,1],tempdf[2,2])
    }
    stocking <- c(stocking,global)
    if (plot==TRUE) {
		#plot(global_save,BIC_save)
      if (length(stocking)>10) {
        stocking2 <- c(stocking2,median(stocking[(length(stocking)-9):length(stocking)]))
      }
      plot(stocking,type="l",lwd=3,col="red",ylim=c(0,1))
      points(stocking2,type="l",lwd=3,col="blue",ylim=c(0,1))

    }
    # Considerer aussi l'apparition de nouveaux Rmax et le t.test par rapport au premier lot en n//vs N+1 à N+5
    win <- 200
    #cat("gene",gene,"avec stocking",length(stocking),"\n")
    if (((gene-(gene%/%win)*win)==0)&(gene>win)&(fast==TRUE)) {
      if (verbose==TRUE) {cat("De : ",(gene-((1+pval_stop2)*win-1)),"-",(gene-win*(pval_stop2)),"à :",(gene-(win-1)),"-",gene)}
      t.test(stocking[(gene-((1+pval_stop2)*win-1)):(gene-win*(pval_stop2))],stocking[(gene-(win-1)):gene])$p.value -> pval_stop
      if (verbose==TRUE) {cat(" - ",pval_stop)
		cat(" pour ",length(sapply(parents,"[[",2))," parents en prélèvements ",prelevement,".\n")}
      if (pval_stop>0.05) {
        pval_stop2 <- pval_stop2 + 1
      } else {
        pval_stop2 <- 1
        #rognage <- 300/round(log2(length(parents)),0) # Nombre de part de prélèvements...
        prelevement <- 3 + round(length(parents)/rognage,0)
        prelevement_save <- round(((prelevement-2)/2),0)+2
      }
    }
    if ((pval_stop2 == 5)&(verbose==TRUE)) {cat("Fin accélérée.\n") ; break}
    #############
    #	Si le modèle est meilleur que le meilleur modèle disponible
    #############
    if ((length(X) == 0) | (length(intersect(which(dt$variables%in%X),my_i))>0)) {
		# Centrer-Réduire BIC et R² pour trouver un compromis
	#	if (length(global_save)>20) {
#			distance_g <- (global-mean(global_save))/sd(global_save)#
			#distance_b <- (bic_temp-mean(BIC_save))/sd(BIC_save)
			#BIC_save_cr <- (BIC_save-mean(BIC_save))/sd(BIC_save)
			#distance <- (distance_g-((1-mean(global_save))/sd(global_save)))^2+(ifelse((distance_b-min(distance_g))<0,0,(distance_b-min(distance_b))))^2
#			distance <- (distance_g-((1-mean(global_save))/sd(global_save)))^2+(distance_b-min(BIC_save_cr))^2
#		} else {
#			distance <- 10000
#		}
      if ((global > bornage_global)|(bic_temp < globalBIC))  { # | (distance<distance_global))
#print("Un modèle enfant ressort.") ; print(paste0(formula(reg)," vs ",formula(super_reg)))
        if (family == "lm") {
          reg <- try(lm(formul,data=data))
          pvals <- summary(reg)[[4]][,4]
        } else if (family=="logit") {
          # Encore du data3 à nettoyer
          reg <- glm(formul,data=data3,family=binomial(logit))
          reg2 <- drop1(reg, test="F")
          pvals <- unlist(reg2[5])[-1]
        }
		if (any(is.na(reg$coefficients))==FALSE) {
        if (any(is.na(pvals))==FALSE) {
          if (any(pvals>0.05)==FALSE) {
			#if (verbose==TRUE){print("Un modèle enfant ressort.") ; print(formula(reg));print(formula(super_reg))}
            # Si BON BIC identifié, faire bootstrap, sauf si modèle déjà sorti !
            if (identical(sort(my_i),sort(my_i_model))==FALSE) {
			#if (verbose==TRUE){print("my_i") ; print(sort(my_i));print(sort(my_i_model))}
			if (formula(reg) != formula(super_reg)){
              #if (verbose==TRUE){print(formula(reg));print(formula(super_reg))}
              #print("On se retrouve dans un modèle à bootstraper")
              resultat <- c()
              #if (verbose == TRUE) {print("Bootstrap d'un modèle enfant potentiellement meilleur...")}
              for (i in 1:500) {
                apprentissage <- sample(1:nrow(data),replace=T)
                test <- c()
                test <- setdiff(1:nrow(data),apprentissage)
				var_a_croiser <- names(get_all_vars(formul,data))
				data_a_croiser <- data[apprentissage,colnames(data)%in%var_a_croiser]
				if ((nrow(na.omit(data_a_croiser))==0)|(min(apply(data_a_croiser,2,function(x){length(unique(x,na.rm=T))}))<2)) {
					next
				}
                if (family == "lm") {
                  reg <- try(lm(formul,data=data[apprentissage,]))
                  #prediction <- predict(reg,newdata=data3[test,])
				  if (is(reg)[1]=="lm"){
					prediction <- summary(reg)$adj.r.squared
					resultat <- c(resultat,prediction)
				  } 
                } else if ((family=="logit")&(length(test)>0)&(length(unique(data3[,which(colnames(data3)%in%Y)][apprentissage]))==2)&(min(table(data3[,which(colnames(data3)%in%Y)][test]))>2)&(min(table(data3[,which(colnames(data3)%in%Y)][apprentissage]))>2)) {
					  # Données catégorielles dans reg logistique :
					# vérifier que les variantes de Y (0 et 1) apparaissent dans le test et l'apprentissage
					ctrl <- c()#
					for (i in 1:ncol(data)){
					  if (is.numeric(data[,i])==FALSE){
						part_common <- sort(intersect(unique(data[apprentissage,i]),unique(data[test,i])))
						part_common <- as.factor(part_common)
						ctrl <- c(ctrl,identical(part_common,as.factor(sort(unique(data[test,i])))))
					  }
					}
					if (all(ctrl==TRUE)==TRUE){ # Intégrer si NA dans test... avec try()
						reg <- glm(formul,data=data3[apprentissage,],family=binomial(logit))
						 prediction <- predict(reg,newdata=data3[test,])
						 prediction <- ifelse(prediction>0.5,"oui","non")
						 tempdf <- data.frame("test" = data3[,which(colnames(data3)%in%Y)][test],prediction)
						 tempdf <- table(tempdf)
						 tempdf <- round(prop.table(tempdf ), 2)
						 if (ncol(tempdf)==2) {resultat <- c(resultat,sum(tempdf[1,1],tempdf[2,2]))
						 } else {resultat <- c(resultat,tempdf[1,1])}
					}
                }
              }
              if (round(quantile(resultat,probs=0.05,na.rm=T),2)>=round(resultat_min,2)){
				if (verbose==TRUE){print("A new model identified for its better minimum capacities.")}
				  if (family == "lm") {
					reg <- try(lm(formul,data=data))
				  } else if (family=="logit") {
					reg <- glm(formul,data=data3,family=binomial(logit))
				  }
				  reg2 <- drop1(reg, test="F")
				  resultat_global<-median(resultat,na.rm=T)
				  resultat_min<-quantile(resultat,probs=0.05,na.rm=T)
				  super_reg <- reg ; super_reg2 <- reg2
				  super_reg_global <- global
				  if (BIC(reg)<globalBIC) {globalBIC <- BIC(reg)}
				  if (bornage_global < global) {bornage_global <- global}
				  my_i_model <- my_i
				  formule = formul
				  #if (!is.na(formule)) {formule <- formula(formule)}				  
				  if (family == "lm") {
				  } else if (family=="logit") {
					modelG <- naiveBayes(form,data=data3)
					prediction <- predict(modelG,data)
					tempdf <- data.frame("test" = data[,which(colnames(data)%in%Y)],prediction)
					tempdf <- table(tempdf)
					tempdf <- round(prop.table(tempdf ), 2) ;
					tempdf_save <- tempdf
				  }
				  if ((verbose==TRUE)&(bornage_global == 0)&(bornage_global!=global)) {print("At least one significant model identified.")}
				  if (verbose==TRUE) {
						cat("Modèle sur enfants avec R2 : ",round(global,2),"\n")
						print(paste("BIC ",BIC(reg)))
						print(paste("Performance inf 95%",round(resultat_min,2)))
						print(formule)
				   }
				}
			}}
          }}}}
    }
    #############
    #	Si l'enfant est meilleur que 1 des parents
    #############
    #print(global)
    #print(dad)
    #print(mom)
    dad_global <- sort(parents[[parents_temp[1]]][[2]])
	if (length(dad_global)==0) {dad_global<-0}
    mom_global <- sort(parents[[parents_temp[2]]][[2]])
	if (length(mom_global)==0) {mom_global<-0}
    bic_dad <- 	parents[[parents_temp[1]]][[3]]
	if (length(bic_dad)==0) {bic_dad<-1E6}	
    bic_mom <- 	parents[[parents_temp[2]]][[3]]
	if (length(bic_mom)==0) {bic_mom<-1E6}	
	#print("global") ; print(global)
	#print("mom_global") ; print(mom_global)
	#print("bic_temp") ; print(bic_temp)
	#print("bic_mom") ; print(bic_mom)
    which(c(dad_global,mom_global)==min(c(dad_global,mom_global),na.rm=T)) -> ind
    #print("B1")
    if (ind[1] == 1) { # Si papa est le plus faible
      #print("B2")
      if ((global >= mom_global)&(bic_temp < bic_mom)) { # Si l'enfant est meilleur que maman, il demeure avec ses parents
        #print("C1")
        parent <- list() ; parent[[1]] <- my_i ; parent[[2]] <- global ; parent[[3]] <- bic_temp
        #print("parent")
        #print(parent)
        parents[length(parents)+1][[1]] <- parent
        #print("integration")
      } else if ((global >= dad_global)&(bic_temp < bic_dad)) { # Si l'enfant est meilleur que le père, il remplace son père
        #print("C2")
        parents[[parents_temp[1]]][[1]] <- my_i
        parents[[parents_temp[1]]][[2]] <- global[1]
        parents[[parents_temp[1]]][[3]] <- bic_temp
      }
    } else { # Si maman est la plus faible
      #print("B3")
      #print(global)
      #print(dad_global)
      #print(mom_global)
      #print(bic_temp)
      #print(bic_dad)
      if ((global >= dad_global)&(bic_temp < bic_dad)) { # Si l'enfant est meilleur que le père, il demeure avec ses parents
        parent <- list() ; parent[[1]] <- my_i ; parent[[2]] <- global ; parent[[3]] <- bic_temp
        #parents[length(parents)+1][[1]] <- parent
        parents[length(parents)+1][[1]] <- parent
      } else if ((global >= mom_global)&(bic_temp < bic_mom)) { # Sinon, il remplace sa mère
		#print(parents[[parents_temp[2]]][[1]])
	    #print(my_i)
		#print(parents_temp[2])
		#if (is.null(parents[[parents_temp[2]]][[1]])==TRUE) {print("Encore nul.") ;return(parents)}
		#print("ABC")
        parents[[parents_temp[2]]][[1]] <- my_i
        parents[[parents_temp[2]]][[2]] <- global[1]
        parents[[parents_temp[2]]][[3]] <- bic_temp[1]
      }
    }
    ########################
    # Les faibles sont tués
    ########################
    #			print("finito")
    #			print(paste("Le faible est tué ",length(parents)))
    parents <- parents[-ind_to_kill2]
    #			print(paste("Il reste ",length(parents)))
    #print(paste("il reste parents : ",length(sapply(parents,"[[",2))))
    #		return(parents)
    # Ajouter les mutations ? si pas de différence, créer une mutation ? mais pas toujours ?
    #print("F")
  }
  if (verbose==TRUE) {
    print(paste("Evolutive approach really made on ",gene,"iterations."))
    print(paste("il reste parents : ",length(sapply(parents,"[[",2))))
  }
  # Vérifier car la formule qui ressort est parfois celle d'un objet formule calculé hors de la fonction...
  # Scénario en cas d'absence de global ?

  # Problème lorsque combiné avec corrigraph var(x) = NULL sur COAT de Jamet.
  if (verbose==TRUE) {
    print("Meilleur modèle sélectionné :")
    print(formule)
    cat("Capacité de prédiction brute : ",round(super_reg_global,2),"\n")
    cat("Capacité de prédiction médiane (bootstrap): ",round(resultat_global,2),"\n")
    cat("Capacité de prédiction inférieure 95% : ",round(resultat_min,2),"\n")
    #cat("Prédiction obtenue (détails)\n")
    #print(tempdf_save
  }
  # CORRIGER LES INDICES i pour qu'ils correspondent bien au tableau...
  #		print("A")
  #		print(colnames(data_save))
  #		print("B")
  #		colnames(data)[my_i_model]
  #		print("C")
  my_i_model <- which(colnames(data_save)%in%colnames(data)[my_i_model])
  #	print(my_i_model)
  #sortie <- list()
  #sortie$R2 <- resultat_global
  #sortie$indices <- my_i_model
  #sortie$BIC <- globalBIC
  #return(sortie)
  } else {print("No model")}
  return(super_reg)
}
