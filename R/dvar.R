#' Diversifies the variables of a dataframe by testing interactions, polynomials, logs... so that evolreg can draw a larger number of model combinations.
#'
#' @param data a dataframe.
#' @param Y the y to predict.
#' @param X variables whose presence we want to force in the model.
#' @param pval 0 to 1. If there are too many variables and the argument wash=TRUE, use this p-value threshold to eliminate the variables whose effect is too insignificant (Risk of eliminating the variables that will have an effect once transformed or in interaction).
#' @param family "lm", "logical" or "lmer". Type of regression
#' @param wash TRUE or FALSE.To select the best variables when there are too many.
#' @param NAfreq from 0 to 1. NA part allowed in the variables. 1 by default (100% of NA tolerate).
#' @param interaction FALSE or TRUE. To allow interactions between variables.
#' @param multix FALSE or TRUE. To allow variable variants (log, exp, polynomial, ^2).
#' @param multidiv FALSE or TRUE. To allow the synthesis of variables combining the ratio of one variable divided by another.
#' @param verbose With an operating report.
#'
#' @return The dataframe of the data and a list of interaction formulas and transformations with associated p-values in their ability to predict Y.
#' @export
#'
#' @examples
#' data(iris)
#' dvar(iris,"Sepal.Length")
dvar <- function(data, Y, X=c(), pval=0.05, family="lm", wash=TRUE, NAfreq=1, interaction=TRUE,
                 multix=TRUE,multidiv=FALSE, verbose=FALSE){
  # Compilation of numerical var
  X_i_num <- c()
  for (i in 1:ncol(data)){
    if ((is.numeric(data[,i])==TRUE)&(length(unique(data[,i]))>1)) {
      X_i_num <- c(X_i_num,i)
    }
  }
  # Compilation of non-numerical var
  X_i_char <- 1:ncol(data)
  X_i_char <- setdiff(X_i_char,X_i_num)
  # Substraction of Y
  if ((is.numeric(Y)==TRUE)&(Y>0)&(Y<=ncol(data))){ind_Y <- Y ; Y <- colnames(data)[Y]
  } else {ind_Y <- which(colnames(data)%in%Y)}
  if (length(ind_Y)==0){stop("Impossible to find Y.")}
  X_i_num <- setdiff(X_i_num,ind_Y)
  X_i_char <- setdiff(X_i_char,ind_Y)
 if (verbose == TRUE){print(paste("Indice of Y: ",ind_Y))}
  # Type of Y
  if (length(unique(data[,which(colnames(data)%in%Y)])) == 2) {
    Y_type <- "binary"
    data[,which(colnames(data)==Y)] <-ifelse(data[,which(colnames(data)%in%Y)]==unique(data[,which(colnames(data)%in%Y)])[1],0,1)
  }
  # NA like category in categorical
  for (i in X_i_char){
    data[which(is.na(data[,i])),i]<-"MANQUANTES"
  }
  #
  variables <- c()
  type <- c()
  var_i <- c()
  p.value <- c()
  weight <- c()
  if (is.numeric(data[,which(colnames(data)%in%Y)])==TRUE) {
	cor.test.NA <- function(x,y) {
		na.omit(data.frame(x,y)) -> temp
		nrow(temp) -> nligne
		if (any(apply(temp,2,function(x){any(is.infinite(x))}))==FALSE) {
			if ((nligne >= 2)&(sd(temp[,1])>0)&(sd(temp[,2])>0)) {
				return(cor.test(x,y)$p.value)
			} else {
				return(NA)
			}
		} else {return(NA)}
	}
    pvals <- apply(data[,X_i_num],2,function(x){cor.test.NA(x,data[,ind_Y])})
    ind_noNA <- which(!is.na(pvals))
    var_i <- c(var_i, X_i_num[ind_noNA])
    variables <- c(variables, colnames(data)[X_i_num[ind_noNA]])
    type <- c(type,rep("Ynum_vs_num",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    #print("If Y presents some categories")
    if (length(unique(data[,ind_Y]))/length(data[,which(colnames(data)%in%Y)]) < 0.1){
      #			print("m.test")
      for (i in X_i_num){
        if (is.numeric(data[,i])==TRUE) {
          #					print(colnames(data)[i])
          temp <- m.test(data[,i],data[,ind_Y],return=FALSE,plot=FALSE,boot=FALSE,verbose=FALSE)
          #print("temp") ; print(temp)
          if (length(temp)>1) {
            var_i <- c(var_i, i)
            variables <- c(variables, colnames(data)[i])
            type <- c(type,"Ycal_vs_num")
            p.value <- c(p.value,temp$p.value)
          }}}
      if (length(X_i_char)>0){
        for (i in X_i_char){
          chisq.test(table(data[,i],data[,ind_Y]))$p.value -> pvals
          if (!is.na(pvals)==TRUE) {
            var_i <- c(var_i, i)
            variables <- c(variables, colnames(data)[i])
            type <- c(type,"Qualitativ")
            p.value <- c(p.value,pvals)
          }
        }
      }
    }
  }
 #return(variables)
  weight <- rep(1,length(var_i))
  dt <- data.frame(var_i,variables,type,p.value,weight)
  #	print(head(dt))
  #############################
  #		Conserve les valeurs qui correspondent au type de Y selon type reg
  #############################
  if (family=="lm") {
    #synthese <-synthese
    dt[dt$type=="Ynum_vs_num",]
  } else if (family=="logit") {
    #		print("logit")
    # AJOUTER LES DONNEES MANQUANTES DANS LES CATEGORIELS
    # A REVOIR +++
    # Remplacer dt par synthese
    # Compilation of numerical var
    dt2 <- data[,1:2]
    for (i in 3:ncol(data)){
      if (length(unique(data[,i]))>1) {
        temp <- as.numeric(as.character(data[,i]))
        if (all(is.na(temp))) {temp <- data[,i]}
        dt2 <- data.frame(dt2,temp)
        colnames(dt2)[ncol(dt2)] <- colnames(data)[i]
      }
    }
    dt2 -> data
    X_i_num <- c()
    for (i in 1:ncol(data)){
      if (is.numeric(data[,i])==TRUE) {
        X_i_num <- c(X_i_num,i)
      }
    }
	if (verbose==TRUE) {print(paste("X_i_num",X_i_num))}
    # Compilation of non-numerical var
    X_i_char <- 1:ncol(data)
    X_i_char <- setdiff(X_i_char,X_i_num)
    # Substraction of Y
    ind_Y <- which(colnames(data)%in%Y)
    X_i_char <- setdiff(X_i_char,ind_Y)
    # Rebinariser Y si necessaire
    for (k in X_i_char){
      data[which(is.na(data[,i])),k]<-"MANQUANTES"
    }
    # Ajouter les variables retenues
    #return(data)

  } else {
    # Autre scenario
    stop()
  }
  #
#return(dt)
  # A controler
  if (multix == TRUE) {
    variables <- c()
    type <- c()
    var_i <- c()
    p.value <- c()
    weight <- c()
	if (verbose==TRUE) {print("log")}
    X_i_num_log <- X_i_num[apply(data[,X_i_num],2,function(x){any(!is.finite(log(x)))})]
    #print(head(data[,X_i_num]))
	if (verbose==TRUE) {print(head(data[,X_i_num_log]))}
    pvals <- apply(data.frame(data[,X_i_num_log]),2,function(x){cor.test.NA(log(x),data[,ind_Y])})
    ind_noNA <- which(!is.na(pvals))
    formule0 <- paste("I(log(",colnames(data)[X_i_num_log],"))")
    var_i <- c(var_i, X_i_num_log[ind_noNA])
    #			print("var_i")
    #			print(var_i)
    variables <- c(variables, formule0[ind_noNA])
    #			print("variables")
    #			print(variables)
    type <- c(type,rep("Ynum_transfo",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    weight	<- weight<-c(rep(1,length(ind_noNA)))
	if (verbose==TRUE) {print("x^2")}
    pvals <- apply(data[,X_i_num],2,function(x){cor.test.NA((x^2),data[,ind_Y])})
    #		print(pvals)
    ind_noNA <- which(!is.na(pvals))
    formule0 <- paste("I(",colnames(data)[X_i_num],"^2)")
    var_i <- c(var_i, X_i_num[ind_noNA])
    #		print(var_i)
    variables <- c(variables, formule0[ind_noNA])
    #		print(variables)
    type <- c(type,rep("Ynum_transfo",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    weight	<- weight<-c(weight,c(rep(1,length(ind_noNA))))
    #	cat("1",length(p.value),"2",length(var_i),"3",length(variables),"4",length(weight))
	if (verbose==TRUE) {print("exp")}
    pvals <- apply(data[,X_i_num],2,function(x){cor.test.NA((exp(x)),data[,ind_Y])})
    #		print(pvals)
    ind_noNA <- which(!is.na(pvals))
    formule0 <- paste("I(exp(",colnames(data)[X_i_num],"))")
    var_i <- c(var_i, X_i_num[ind_noNA])
    #		print(var_i)
    variables <- c(variables, formule0[ind_noNA])
    #		print(variables)
    type <- c(type,rep("Ynum_transfo",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    weight	<- weight<-c(weight,c(rep(1,length(ind_noNA))))
    dt_inter <- data.frame(var_i,variables,type,p.value,weight)
    dt <- rbind(dt,dt_inter)
  }
  ##########
  if (multidiv == TRUE) {
    variables <- c()
    type <- c()
    var_i <- c()
    p.value <- c()
	if (verbose==TRUE) {print("Interaction j/i & i/j")}
    for (i in which(dt$type=="Ynum_vs_num"|dt$type=="Ynum_transfo")) {
      for (j in setdiff(which(dt$type=="Ynum_vs_num"|dt$type=="Ynum_transfo"),i)) {
        # INTERACTION i/j
        #print("Interaction i/j")
        #print(dt$variables[j])
        if ((any(within(data, result <- eval(parse(text =dt$variables[j])))$result==0)==FALSE)&(any(is.na(within(data, result <- eval(parse(text =dt$variables[j])))$result))==FALSE)){
          #print("go")
          formule0 <- paste("I(",dt$variables[i],"/",dt$variables[j],")")
          formule1 <- paste(Y,'~',formule0)
          formule2 <- as.formula(formule1)
          reg <- lm(formule2,data=data)
          #print(summary(reg))
          preg <- pf(summary(reg)$fstatistic[1], summary(reg)$fstatistic[2],
                     summary(reg)$fstatistic[3], lower.tail = FALSE)
          if (!is.na(preg)==TRUE) {
            var_i <- c(var_i, dt$var_i[i])
            variables <- c(variables, formule0)
            type <- c(type,"Interaction")
            p.value <- c(p.value,preg)
          }
        }
        #if (any(eval(parse(text = dt$variables[i])))==FALSE){
        #print(dt$variables[i])
        if ((any(within(data, result <- eval(parse(text =dt$variables[i])))$result==0)==FALSE)&(any(is.na(within(data, result <- eval(parse(text =dt$variables[i])))$result))==FALSE)){
          #print("go1")
          #print("Interaction j/i")
          # INTERACTION j/i
          formule0 <- paste("I(",dt$variables[j],"/",dt$variables[i],")")
          formule1 <- paste(Y,'~',formule0)
          formule2 <- as.formula(formule1)
          reg <- lm(formule2,data=data)
          preg <- pf(summary(reg)$fstatistic[1], summary(reg)$fstatistic[2],
                     summary(reg)$fstatistic[3], lower.tail = FALSE)
          if (!is.na(preg)==TRUE) {
            var_i <- c(var_i, dt$var_i[i])
            variables <- c(variables, formule0)
            type <- c(type,"x1/x2")
            p.value <- c(p.value,preg)
          }
        }
      }
    }
	weight <- rep(1,length(var_i))
	dt_inter <- data.frame(var_i,variables,type,p.value,weight)
	dt <- rbind(dt,dt_inter)
  }
  if (multix == TRUE) {
    variables <- c()
    type <- c()
    var_i <- c()
    p.value <- c()
    weight <- c()
    ##########
    #	cat("1",length(p.value),"2",length(var_i),"3",length(variables),"4",length(weight))
	if (verbose==TRUE) {print("1/x")}
    X_i_num_inv <- apply(data[,X_i_num],2,function(x){return(any((x==0)))})
    ind_noNA <- which(!is.na(X_i_num_inv))
    #			print(X_i_num_inv)
    X_i_num_temp <- X_i_num[ind_noNA]
    X_i_num_inv <- X_i_num_inv[ind_noNA]
    X_i_num_inv <- X_i_num_temp[-(c(1:length(X_i_num_temp))[X_i_num_inv])]
    #			print(X_i_num_inv )
    #print(head(data[,X_i_num_log]))
    pvals <- apply(data[,X_i_num_inv],2,function(x){cor.test.NA((1/x),data[,ind_Y])})
    ind_noNA <- which(!is.na(pvals))
    formule0 <- paste("I(1/",colnames(data)[X_i_num_inv],")")
    var_i <- c(var_i, X_i_num_inv[ind_noNA])
    #			print("var_i")
    #			print(var_i)
    variables <- c(variables, formule0[ind_noNA])
    #			print("variables")
    #			print(variables)
    type <- c(type,rep("Ynum_transfo",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    weight	<- weight<-c(weight,rep(1,length(ind_noNA)))
    #
    #	cat("1",length(p.value),"2",length(var_i),"3",length(variables),"4",length(weight))
	if (verbose==TRUE) {print("poly")}
    X_i_num_poly <- apply(data[,X_i_num],2,function(x){return(any(is.na(x)))})
    X_i_num_poly <- which(X_i_num_poly==FALSE)
    X_i_num_temp <- X_i_num[X_i_num_poly]
    #print(X_i_num_temp)*
    X_i_num_poly <- apply(data.frame(data[,X_i_num_temp]),2,function(x){return(length(unique(x)))})
    X_i_num_poly <- which(X_i_num_poly>2)
	#print(X_i_num_poly)
	if (length(X_i_num_poly)>0){
		X_i_num_temp <- X_i_num_temp[X_i_num_poly]
		#print(colnames(data)[X_i_num_temp])
		X_i_num_temp_vari <- apply(data.frame(data[,X_i_num_temp]),2,function(x){
			return(cor(x, data[,ind_Y]))})
		#print(X_i_num_temp_vari)
		X_i_num_temp <- X_i_num_temp[which(!is.na(X_i_num_temp_vari))]
		#print(colnames(data)[X_i_num_temp])
		if(length(X_i_num_temp)==0){
			pvals <- NA
		}else {
			#print('pb')
			pvals <- apply(data.frame(data[,X_i_num_temp]),2,function(x){
			  reg_temp <- lm( data[,ind_Y]~poly(x,2))
			  pval_temp <- pf(summary(reg_temp)$fstatistic[1], summary(reg_temp)$fstatistic[2],summary(reg_temp)$fstatistic[3], lower.tail = FALSE)
			  return(pval_temp)})
		}
	}else {pvals <- NA}
#	print('nop')
    #		print(pvals)
    ind_noNA <- which(!is.na(pvals))
    formule0 <- paste("poly(",colnames(data)[X_i_num_temp],",2)")
    var_i <- c(var_i, X_i_num_temp[ind_noNA])
    #		print(var_i)
    variables <- c(variables, formule0[ind_noNA])
    #		print(variables)
    type <- c(type,rep("Ynum_vs_poly",length(ind_noNA)))
    p.value <- c(p.value,pvals[ind_noNA])
    weight	<- weight<-c(weight,c(rep(1.5,length(ind_noNA))))
    dt_inter <- data.frame(var_i,variables,type,p.value,weight)
    dt <- rbind(dt,dt_inter)
  }
  #print("B")
  #print("multix")
  if (wash==TRUE) {
    if (verbose == TRUE) {print(paste("Il y a pour le moment ",nrow(dt)," variables."))}
    ########################
    #	Prenettoyage
    ########################
    #if (length(wash)>0) {
    if (nrow(dt) > 40) {
      var_alea <- qbinom(0.95,nrow(dt),pval)
      dt0 <- dt[dt$p.value<pval,]
      dt1 <- dt[dt$p.value<pval/nrow(dt),]
      ind <- order(dt0$p.value)
      if (var_alea < length(ind)) {
        dt2 <- dt[ind[1:(length(ind)-var_alea)],]
        n1 <- nrow(dt1)-40
        n2 <- nrow(dt2)-40
        n12 <- c(n1 , n2)
        n12[n12>0]->n12
		if (length(n12)>0){n12 <- min(n12)
		}else{n12 <- 0}
        if (n12 == n1) {
          dt <- dt1
        } else if (n12 == n2) {
          dt <- dt2
        } else {
          # nothing
        }
      }
    }
    #}
	if (verbose==TRUE) {
    print("Apres nettoyage")
    print(paste("Il reste ",nrow(dt)," variables."))
	}
    #		print(dt)
  }
  #print(head(dt))
  ########################
  #	Interactions
  ########################
  if (interaction == TRUE){
  if (verbose==TRUE) {print("Interaction")}
    variables <- c()
    type <- c()
    var_i <- c()
    p.value <- c()
    #			for (i in dt$var_i[dt$type=="Ynum_vs_num"|dt$type=="Ynum_vs_log"][1:(length(dt$var_i[dt$type=="Ynum_vs_num"|dt$type=="Ynum_vs_log"])-1)]) {
    #				for (j in dt$var_i[dt$type=="Ynum_vs_num"|dt$type=="Ynum_vs_log"][2:length(dt$var_i[dt$type=="Ynum_vs_num"|dt$type=="Ynum_vs_log"])]) {
    X_i_num_for_inter <- which((dt$type=="Ynum_vs_num")|(dt$type=="Ynum_transfo"))
    #print(dt$variables[X_i_num_for_inter])
    #print(length(X_i_num_for_inter))
    #print((length(X_i_num_for_inter)-1))
    #print(dt$variables[X_i_num_for_inter[1:(length(X_i_num_for_inter)-1)]])
    for (i in 1:(length(X_i_num_for_inter)-1)) {
      #			for (i in X_i_num_for_inter[listeri]) {
      #print("D")
      #print(paste("i",X_i_num_for_inter[i]))
      #print("C")
      for (j in (i+1):length(X_i_num_for_inter)) {
        # INTERACTION
        #	print(1:(length(X_i_num_for_inter)-1))
        #print(i)
        #	print((i+1):length(X_i_num_for_inter))
        #print("A")
        #print(paste("i ",X_i_num_for_inter[i]," et j ",X_i_num_for_inter[j]))
        #print("B")
        formule0 <- paste(dt$variables[X_i_num_for_inter[i]],'*',dt$variables[X_i_num_for_inter[j]])
		formule1 <- paste(Y,'~',formule0)
		formule2 <- as.formula(formule1)
		# Check NA values
		ind_temp <- which(colnames(data)%in%names(get_all_vars(formule2, data = data)))
		if (nrow(na.omit(data[,ind_temp]))>2){
			#print(formule0)
			reg <- lm(formule2,data=data)
			preg <- pf(summary(reg)$fstatistic[1], summary(reg)$fstatistic[2],
					   summary(reg)$fstatistic[3], lower.tail = FALSE)
		} else {
			preg <- NA
		}
        if (!is.na(preg)==TRUE) {
          var_i <- c(var_i, dt$var_i[X_i_num_for_inter[i]])
          variables <- c(variables, formule0)
          type <- c(type,"Interaction")
          p.value <- c(p.value,preg)
        }
      }
    }
    weight <- rep(2,length(var_i))
    dt_inter <- data.frame(var_i,variables,type,p.value,weight)
    dt <- rbind(dt,dt_inter)
  }
  if (wash==TRUE) {
  if (verbose==TRUE) {
    print("Nettoyage")
    print(paste("Il y a pour le moment ",nrow(dt)," variables."))
	}
    #######################
    #	Nettoyage
    #######################
    seuil <- 100
    if (nrow(dt) > seuil) {
	if (verbose==TRUE) {
      print("Nettoyage")
      print(nrow(dt))
	  }
      var_alea <- qbinom(0.95,nrow(dt),pval)
      dt0 <- dt[dt$p.value<pval,]
      dt1 <- dt[dt$p.value<pval/nrow(dt),]
      ind <- order(dt0$p.value)
      if (var_alea < length(ind)) {
        dt2 <- dt[ind[1:(length(ind)-var_alea)],]
        n1 <- nrow(dt1)-seuil
		if (verbose==TRUE) {
        print(paste("reste_",n1))}
        n2 <- nrow(dt2)-seuil
		if (verbose==TRUE) {
        print(paste("reste_",n2))}
        n12 <- c(n1 , n2)
        n12[n12>0]->n12
        if (length(n12)>0){n12 <- min(n12)
        }else{n12 <- 0}
        if (n12 == n1) {
          dt <- dt1
        } else if (n12 == n2) {
          dt <- dt2
        } else {
		if (verbose==TRUE) {
          print(n12)}
          # nothing
        }
      }
    }
    #	if (length(wash) > 0) {
    #		ind <- order(dt$p.value)
    #		if (wash > nrow(dt)) {wash <- nrow(dt)}
    #		ind <- ind[1:wash]
    #		dt <- dt[ind,]
    #	}
    #}
	if (verbose==TRUE) {
    print("Apres nettoyage")
    print(paste("Nombre de variables retenues ",nrow(dt)))
	}
  }
  if (length(X)>0) {
    #ind_X <- which(colnames(data)%in%X)
    ind_X_bis <- which(dt$variables%in%X)
    if (length(X) != length(ind_X_bis)) {
      for (i in X) {
        ind_X_bis <- which(dt$variables%in%i)
        if (length(ind_X_bis)==0) {
          dtligne <- data.frame(var_i=which(colnames(data)%in%i),
                                variables=i,
                                type = "Ynum_vs_num",
                                p.value = 0,
                                weight = 1)
          dt <- rbind(dt,dtligne)
        }
      }
    }
  }
  sortie <- list()
  sortie$variables <- dt
  sortie$data <- data
  return(sortie)
}
