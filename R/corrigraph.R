#' igraph of correlated variables global or in relation to y
#'
#' @param data a data.frame
#' @param colY a vector of indices or variables to predict. To force the correlogram to display only the variables correlated to a selection of Y.
#' @param pval the maximum permissible p-value for the display
#' @param exclude the minimum threshold of displayed correlations - or a vector of threshold in this order : c(cor,mu,prop)
#' @param ampli coefficient of amplification of vertices
#' @param return if return=T, returns the correlation matrix of significant correlation.
#' @param wash automatically eliminates variables using differnts methods when there are too many variables (method = NA, stn (signal-to-noise ratio), sum, length).
#' @param multi to ignore multiple regressions and control only single regressions.
#' @param mu to display the effect on median/mean identified by m.test().
#' @param prop to display the dependencies between categorical variables identified by chisq.test().
#' @param layout to choose the network organization method - choose "fr", "circle", "kk" or "3d".
#' @param verbose to see the comments.
#'
#' @return Correlation graph network (igraph) of the variables of a data.frame. Pay attention to the possible presence of non-numeric variables or missing data. Grouping of correlated variables: the vertices (circles) correspond to the variables. The more a variable is connected, the larger it appears. The color of the lines reflects the nature of the correlation (positive or negative).  The size of the lines is the value of the correlation from 0 to 1. All these correlations are significant (pval < 0.01). The coloured groupings reflect families of inter-correlated variables. BLUE: positive correlation - RED: negative correlation
#' @return When mu is TRUE or prop : we see the connexion with mean effect (orange) and chisq effect (pink)
#' @return The size of orange edge and pink edge depend of p-values (-1*log10(p-value)/10) of kruskal.test() and chisq.test().
#' @return
#' @return  When indicating Y's in colY, the correlogram will identify the correlated X's, then the remaining X's correlated to these X's, and so on.
#' @return X's not related to these Y's are excluded.
#' @return The blue always displays the positive correlations and the red, negative correlations. When the display is green, it means that the predictive (~correlation) capacity of the variable can be reinforced by adding a 2nd variable in a multiple regression model (interaction X1+X2, X1*X2 or X1+X1:X2) better than X1 or X2 alone.
#' @return Correlations between X or Y of the same level are neglected.
#' @return The color of the vertices makes it possible to identify the correlated variables alone in a significant way (blue: positive, red: negative, purple: positive or negative depending on the Y).
#' @return The values displayed to the right of the Ys (colY) correspond to the maximum predictive capacity of these Ys by one or two variables.
#'
#' @import utils
#' @rawNamespace import(igraph, except = decompose)
#' @importFrom igraph graph_from_adjacency_matrix
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
#' corrigraph(airquality,layout="3d")
#' # Example 3
#' data(airquality)
#' corrigraph(airquality,c("Ozone","Wind"))
#' # Example 4
#' data(iris)
#' corrigraph(iris,mu=TRUE)
#' # Example 5
#' require(MASS) ; data(Aids2)
#' corrigraph(Aids2 ,prop=TRUE,mu=TRUE,exclude=c(0.3,0.3,0))
corrigraph <- function(data,colY=c(),pval=0.05,exclude=c(0,0,0), ampli=4,return=FALSE,wash="stn",multi=TRUE,
                       mu=FALSE,prop=FALSE,layout="fr",verbose=FALSE) {
  # Fonction réalisée par Antoine Massé
  # Ctrl Alt Shift R
  # Version 03
  # Février 2021
  databrut<-data # saving
  # Control 1 - is.numeric ?
  which(sapply(data, is.numeric)) -> temp_id_num ; temp_id_factor <- c()
  if ((mu==FALSE)&(prop==FALSE)){
    if (length(temp_id_num)<2) {stop("Error! Not enough numerical variable. Try mu=TRUE or prop=TRUE ?\n")}
    data <- data[,temp_id_num]
    which(sapply(data, is.numeric)) -> temp_id_num
    if (length(temp_id_num)!=length(1:ncol(data))) {cat("Warning! Presence of non-numeric variables that cannot be taken into account.\n")}
  } else {
    temp_id_factor <- setdiff(1:ncol(data),temp_id_num)
  }
  # Control 2 - var is NULL ?
  which(sapply(data.frame(data[,temp_id_num]), var, na.rm=T)!=0)-> temp_var ; temp_var	<- temp_id_num[temp_var]
  if (length(temp_var) != length(temp_id_num)) {cat("Warning! Some variables have a null variance and cannot be taken into account.\n")}
  data <- data[,union(temp_var,temp_id_factor)] ; temp_id_num <- which(sapply(data, is.numeric)) ; temp_id_factor  <- setdiff(1:ncol(data),temp_id_num)
  # Control 3 - is.na ?
  if (any(is.na(data))==TRUE) {cat("Warning ! Presence of missing values.\n")}
  if (ncol(data)>50) {cat("Warning! The calculation time increases exponentially with the number of variables (do not exceed 50).\n")}
  if (length(colY)>0) {
    if((mu==TRUE)|(prop==TRUE)){stop("Error! mu and prop can't be used width colY.\n")}
    #print("mode 2")
    # If modelization of Y
    if (is.numeric(colY)) {
      colnames(databrut)[colY] -> colY
    }
    which(colnames(data)%in%colY) -> indY
    if (length(colY)!=length(indY)){stop("ERROR! some Y's provided cannot be analyzed.\n")}
    if (length(indY)==0) {stop()}
    indices_all <- 1:ncol(data)
    # Matrice de correlation
    matrix(rep(0,ncol(data)^2),ncol(data),ncol(data))->mymat
    rownames(mymat) <- colnames(data)
    colnames(mymat) <- colnames(data)
    sum_mymat2 <- ncol(data)^2
    levels <- rep(0,ncol(mymat))
    borders <- rep("black",ncol(mymat))
    vertices <- rep("white",ncol(mymat))
    level <- 1
    levels[indY] <- level
    warning<-0
    check_mdl <- function(x) {
      data_temp <- data.frame("data[,i]"=data[,i],"data[,j]"=data[,j],"x"=x)
      data_temp <- na.omit(data_temp)
      nb_subset <- nrow(data_temp)
      if ((min(apply(data_temp,2,var))!=0)&(nb_subset>3)) {
        regA <- lm(data_temp[,1]~data_temp[,2]+data_temp[,3])
        regB <- lm(data_temp[,1]~data_temp[,2]+data_temp[,2]:data_temp[,3])
        regC <- lm(data_temp[,1]~data_temp[,2]+data_temp[,2]*data_temp[,3])
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
          temp2 <- ifelse(temp2$p.value<= pval,abs(temp2$estimate),0)
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
        } else if (max(pvalsA)>pval){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
        } else {
          tablo <- data.frame(R2 = abs(summary(regA)$adj.r.squared),BIC = BIC(regA))
          synthese <- rbind(synthese,cbind(tablo,synthese2))}
        if (any(is.na(pvalsB))){	rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
        } else if (max(pvalsB)>pval){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
        } else {
          tablo <- data.frame(R2 = abs(summary(regB)$adj.r.squared),BIC = BIC(regB))
          synthese <- rbind(synthese,cbind(tablo,synthese2))
        }
        if (any(is.na(pvalsC))){	rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
        } else if (max(pvalsC)>pval){rbind(synthese,data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))
        } else {
          tablo <- data.frame(R2 = abs(summary(regC)$adj.r.squared),BIC = BIC(regC))
          synthese <- rbind(synthese,cbind(tablo,synthese2))}
      } else {return(data.frame(R2 = NA,BIC = NA, R2_2 = NA, BIC_2 = NA))}
      return(synthese)
    }
    mymat_interconnect <- mymat
    while (sum_mymat2  != sum(mymat)) { # Tant que la matrice evolue à chaque cycle...
      sum_mymat2 <- sum(mymat)
      indices <- c()
      for (i in indY) { # Modelization of Y
        indices_X <- setdiff(indices_all,indY)
        if (length(indices_X)==0){break}
        for (j in indices_X) { # En fonction d'un X
          temp_tab <- data.frame(data[,i],data[,j])
          na.omit(temp_tab) -> temp_tab
          nb_subset <- nrow(temp_tab)
          if ((min(apply(temp_tab,2,var)) != 0)&(nb_subset>3)) {
            temp <- cor.test(data[,i],data[,j],na.rm=T)
            temp_save <- ifelse(temp$p.value<= pval,temp$estimate,0) # saving cor neg or pos
            temp <- ifelse(temp$p.value<= pval,abs(temp$estimate),0) # saving abs(cor)
            if (temp_save>0) {  # coloring vertices if significant correlation neg or pos // Y
              if (borders[j]=="purple") {
              } else if (borders[j]=="black") {
                borders[j]<-"blue" ; vertices[j] <- "cyan"
              } else if (borders[j]=="red") {
                borders[j]<-"purple"  ;  vertices[j] <- "#FBC5F8"
              }
            } else if (temp_save < 0) {
              if (borders[j]=="purple") {
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
              cat("Warning! Failure to account for missing data generated zero variances on some variables that had to be ignored.\n")}
          }
          if (temp!=0){
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
      indices_all <- setdiff(indices_all,indY)
      indY <- unique(indices)
      level <- level+1
      levels[indY] <- level
      if (level==5) {break}
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
    #print("mode1")
    #	 Matrice de correlation
    matrix(rep(0,ncol(data)^2),ncol(data),ncol(data))->mymat
    rownames(mymat) <- colnames(data)
    colnames(mymat) <- colnames(data)
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
      for (i in temp_id_factor) {
        if (length(temp_id_factor)>5){setTxtProgressBar(barre,i)}
        if (min(table(data[,i]))<3) {
          temp_id_factor <- setdiff(temp_id_factor,i)
          next
        }
        if ((prop==TRUE)&(length(temp_id_factor)>1)) {
          for (j in setdiff(temp_id_factor,i)) {
            tablo <- table(data[,i],data[,j])
            #print(tablo)
            if ((ncol(tablo)>=2) & (nrow(tablo)>=2)) {
              #cat("Test entre",i,"&",j,"\n")
              #print(tablo)
              if (min(tablo)<5) {
                #cat("Effectif insuffisant\n")
                if (quantile(tablo,probs=c(0.2))<5) {
                  #round(length(as.vector(tablo))/100*20,0) -> temp
                  #temp <- sort(as.vector(tablo))[temp]
                  #if (temp==0) {print(tablo) ; next}
                  #if (temp<5) {
                  #cat("Même sur 80%\n")
                  next
                }
              }
              oldw <- getOption("warn")
              options(warn = -1)
              pvals <- try(chisq.test(tablo)$p.value,silent=TRUE)
              options(warn = oldw)
              #cat("Test entre",colnames(data)[i],"&",colnames(data)[j],"pval",pvals,"\n")
              if (pvals <= pval) {
                #print(pvals)
                -1*log10(pvals)/10 -> temp ; ifelse(temp>1 ,1 ,temp)-> temp
                #print(temp)
                mymat[i,j] <- temp ; mymat[j,i] <- temp
                mymat_interconnect[i,j] <- 2 ; mymat_interconnect[j,i] <- 2
              }
            } else {next}
          }
        }
      }
      if (length(temp_id_factor)>5){close(barre)}
      #print(colnames(data)[temp_id_factor])
      if ((mu == TRUE)&(length(temp_id_factor)>0)) {
        if (length(temp_id_factor)>5){
          barre <- txtProgressBar(min=0,max=length(temp_id_factor),width=50)
          cat("\nmu calculation\n")
          cat(paste(rep("=",50),collapse=""),"\n")
        }
        for (i in temp_id_factor) {
          if (length(temp_id_factor)>5){setTxtProgressBar(barre,i)}
          if (verbose==TRUE) {cat("\nAnalysis of the cat:") ; print(colnames(data)[i])}
          for (j in temp_id_num) {
            if (verbose==TRUE) {cat("\n\tvs numerical: ") ; print(colnames(data)[j])}
            if (m.test(data[,j],data[,i],pval=pval,return=FALSE,verbose=FALSE,plot=FALSE,boot=FALSE)) {
              #by(data[,j],data[,i],median,na.rm=T)-> temp
              #temp <- mean(abs(temp-median(temp)))/mean(by(data[,j],data[,i],function(x){return(hdi(x)[2]-hdi(x)[1])}))
              #temp <- sd(temp-median(temp,na.rm=T))/mean(by(data[,j],data[,i],function(x){return(hdi(x)[2]-hdi(x)[1])}))
              #temp <- mean(temp,na.rm=TRUE)/mean(by(data[,j],data[,i],function(x){return(hdi(x)[2]-hdi(x)[1])}))
              #print(temp)
              #temp <- ifelse(temp>1,1,ifelse(temp<0.1,0.1,temp))
              #print(temp)
              kruskal.test(data[,j],data[,i])$p.value -> pvals
              #print("pval") ;  print(pvals)
              -1*log10(pvals)/10 -> temp ; ifelse(temp>1 ,1 ,temp)-> temp
              #print(temp)
              mymat[i,j] <- temp ; mymat[j,i] <- temp
              mymat_interconnect[i,j] <- 3 ; mymat_interconnect[j,i] <- 3
            }
          }
        }
        if (length(temp_id_factor)>5){close(barre)}
      }
    }
    # Numerical correlations
    if (length(temp_id_num)>1) {
      barre <- txtProgressBar(min=0,max=length(temp_id_num),width=50)
      cat("\ncor calculation\n")
      cat(paste(rep("=",50),collapse=""),"\n")
      for (i in temp_id_num) {
        setTxtProgressBar(barre,i)
        for (j in setdiff(temp_id_num,i)) {
          temp_tab <- data.frame(data[,i],data[,j])
          na.omit(temp_tab) -> temp_tab
          if ((var(temp_tab[,1])!=0)&(var(temp_tab[,2])!=0)&(nrow(temp_tab)>2)){
            temp <- cor.test(data[,i],data[,j],na.rm=T)
            temp <- ifelse(temp$p.value<= pval,temp$estimate,0)
          }else{
            temp<-0
            if(warning==0){warning <- 1}
          }
          mymat[i,j] <- temp ; mymat[j,i] <- temp
          mymat_interconnect[i,j] <- temp ; mymat_interconnect[j,i] <- temp
        }
      }
      close(barre)
    }
    cat("\n")
    if (warning==1) {cat("Warning! Failure to account for missing data generated zero variances on some variables that had to be ignored.\n")}
    #print(mymat)
    #print(mymat_interconnect)
    pas = (ncol(mymat)-50)/20
    if ((is.na(wash)!=TRUE)&(wash!=FALSE)) {
      indices <- apply(abs(mymat),2,sum)
      indices <- which(indices>0)
      mymat <- mymat[indices,indices]
      mymat_interconnect <- mymat_interconnect[indices,indices]
      vertices <- vertices[indices]
      while(ncol(mymat)>50){
        #mymat <- ifelse(abs(mymat)< exclude,0,mymat)
        if (wash=="stn"){
          if (ncol(mymat) > 50) {
            moyennes <- apply(abs(mymat),2,mean)
            deviation <- apply(abs(mymat),2,sd)
            snp <- moyennes/deviation
            indices <- 1:ncol(mymat)
            indices [order(snp,decreasing=T)]
            indices <- indices[1:(ncol(mymat)-pas)]
            mymat <- mymat[indices,indices]
            mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
            indices <- apply(abs(mymat),2,sum)
            indices <- which(indices>1)
            mymat <- mymat[indices,indices]
            mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
          }
        } else if (wash=="sum"){
          somme <- apply(abs(mymat),2,sum)
          indices <- 1:ncol(mymat)
          indices [order(somme,decreasing=T)]
          indices <- indices[1:(ncol(mymat)-pas)]
          mymat <- mymat[indices,indices]
          mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
          indices <- apply(abs(mymat),2,sum)
          indices <- which(indices>1)
          mymat <- mymat[indices,indices]
          mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
        } else if (wash=="length"){
          taille <- apply(abs(mymat),2,function(x){return(length(x[x>0]))})
          indices <- 1:ncol(mymat)
          indices [order(taille ,decreasing=T)]
          indices <- indices[1:(ncol(mymat)-pas)]
          mymat <- mymat[indices,indices]
          mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
          indices <- apply(abs(mymat),2,sum)
          indices <- which(indices>1)
          mymat <- mymat[indices,indices]
          mymat_interconnect <- mymat_interconnect[indices,indices] ; vertices <- vertices[indices]
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
                 arrow.mode = 0, edge.width = (abs(E(net)$weight) * 2)^3,
                 edge.color = E(net)$colour,label.dist=sqrt(betweenness(net))*ampli/4+10)
      }else if (ncol(mymat)<50) {
        clp <- cluster_optimal(net)
        class(clp)
        plot(clp, net, layout = l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color="yellow",
             edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
      } else {
        plot(net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
             edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
      }
    } else {
      #print(mymat)
      net <- graph_from_adjacency_matrix(mymat, weighted=T,mode="directed")
      #net <- simplify(net, remove.multiple = TRUE, remove.loops = TRUE) # elaguer les liens redondants
      net_interconnect <- graph_from_adjacency_matrix(mymat_interconnect, weighted=T,mode="directed")
      #net_interconnect <- simplify(net_interconnect, remove.multiple = TRUE, remove.loops = TRUE) # elaguer les liens redondants
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
      E(net)$weight <- abs(E(net)$weight)
      #print(E(net_interconnect)$weight)
      E(net)$colour <- ifelse(E(net_interconnect)$weight<0,"red",
                              ifelse(E(net_interconnect)$weight==2,"pink",
                                     ifelse(E(net_interconnect)$weight==3,"orange","blue")))
      if (layout=="circle") {l <- layout_in_circle(net)
      } else if (layout=="kk"){l <- layout_with_kk(net)
      } else  {l <- layout_with_fr(net)}
      if ((layout=="3d") |(layout=="3D")) {
        l <- layout_with_fr(net,dim=3)
        rglplot( net, layout = l, vertex.size = sqrt(betweenness(net))*ampli/4+10,
                 vertex.color = vertices, edge.arrow.size = 0,
                 arrow.mode = 0, edge.width = (abs(E(net)$weight) * 2)^3,
                 edge.color = E(net)$colour,label.dist=sqrt(betweenness(net))*ampli/4+10)
      } else {
        plot(net,layout=l,vertex.size=sqrt(betweenness(net))*ampli+10,vertex.color=vertices,
             edge.arrow.size =0,arrow.mode=0,edge.width=(abs(E(net)$weight)*2)^3, edge.color =E(net)$colour)
      }
    }
    if (return==TRUE){return(mymat)}
  }
}
