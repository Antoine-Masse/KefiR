#' Identify the sample size needed to detect a difference between two samples
#'
#' @param mean1 Sample Average 1
#' @param mean2 Sample Average 2
#' @param sd1 Sample standard deviation 1
#' @param sd2 Sample standard deviation 2
#' @param iter Number of iterations per bootstrap
#' @param conf.level Expected confidence on Student's test when comparing 2 samples (0 to 1)
#' @param conf Percentage of times I expect to have a p-value that responds to my test if there is an effect when comparing 2 samples
#' @param nmin Minimum number of individuals accepted for the sample.
#' @param nmax Maximum number of individuals attainable for the sample.
#'
#' @return The required sample size to detect a difference in means by a Student test with the desired significance (\code{conf.level}) and in a minimum number of cases (\code{conf}). If the required sample size exceeds \code{nmax}, a warning message is printed.
#' 
#' @export
#'
#' @examples
#' identify_ech(mean1=100, mean2=101, sd1=2, sd2=3, nmax=10000)
identify_ech <- function(mean1,mean2,sd1,sd2,iter=500,conf.level = 0.99,conf=0.95,nmin=2,nmax=1000) {
  # Version 02 - 17/11/2020
  # par Antoine Masse
  conftest <- 0 ;  intervalle <- 1 ; test <- 0 ; nprec <- 0 ;
  nmint = nmin ; nmaxt = nmax
  while ((intervalle >= 1) & (nprec < nmax) & (test==0)) {
    n <- nmint
    intervalle <- round((nmaxt-nmint)/9,0)
    if (intervalle < 1) {intervalle <-1}
    if(intervalle == 1) {test = test+1}
    cat("Pas testé : ",intervalle,"\n")
    while ((conftest < conf) & (n <= nmax)) {
      pval <- c()
      for (i in 1:iter) {
        pop_ref <- rnorm(n, mean = mean1, sd = sd1)
        pop_test <- rnorm(n, mean = mean2, sd = sd2)
        pval <- c(pval,t.test(pop_ref,pop_test)$p.value)
      }
      conftest <- length(pval[pval<(1-conf.level)])/iter
      nprec_prec = n-intervalle
      nprec <- n
      n <- n+intervalle
    }
    cat("Next loop: from ",nprec_prec," to ",nprec,"\n")
    nmint = nprec_prec
    nmaxt = nprec
    confdef <- conftest
    conftest <- 0
  }
  if (nprec == nmax) {
    cat("Warning : The required sample size requires exceeding the maximum population indicated.\n")
  } else {cat("To detect an effect under experimental conditions,\nthe population must have a minimum size of : ",nprec,"\n")}
  cat("Confidence expected : ",confdef,"\n")
}




#' Determine the differences in mean that I can detect for a known sample size
#'
#' @param n Size of my sample.
#' @param mean Average of my reference sample.
#' @param sd Standard deviation of my sample.
#' @param conf.level Expected confidence on Student's test when comparing 2 samples (0 to 1)
#' @param conf Frequency of times we wish to have the desired p-value defined
#' @param iter Number of iterations per bootstrap
#' @param cut Resolution (number of steps)
#'
#' @return The minimum detectable effect size (mean difference) that can be identified with the given sample size, confidence level, and other parameters. If it is not possible to find a significant variation with the provided sample size, a message will be printed.
#' 
#' @export
#'
#' @examples
#' # What minimum variation could I detect for a sample of 20 individuals (mean 170, sd 20)?
#' check_ech(n=20, mean = 170, sd = 20)
check_ech <- function(n, mean,sd,conf.level=0.99,conf=0.95, iter=200,cut=100){
  # Version 02 - 17 novembre 2020
  # Author : Antoine Masse
  seq <- (10:-10)
  seq <- 10^(seq)
  # conf
  conftest <- 1
  i <- 1
  res_min <- res_max <- 0
  while ((conftest > conf) & (i < length(seq))) {
    moy_sup = mean+mean*seq[i]
    moy_inf = mean-mean*seq[i]
    pval_sup <- c()
    pval_inf <- c()
    for (j in 1:iter) {
      pop_ref <- rnorm(n, mean = mean, sd = sd)
      pop_test_sup <- rnorm(n, mean = moy_sup, sd = sd)
      pop_test_inf <- rnorm(n, mean = moy_inf, sd = sd)
      pval_sup <- c(pval_sup, t.test(pop_ref,pop_test_sup)$p.value)
      pval_inf <- c(pval_inf, t.test(pop_ref,pop_test_inf)$p.value)
    }
    conftest_sup <- length(pval_sup[pval_sup<(1-conf.level)])/iter
    conftest_inf <- length(pval_inf[pval_inf<(1-conf.level)])/iter
    conftest <- min(conftest_sup,conftest_inf)
    if ((i > 1)&(conftest < conf)) {
      res_min <- seq[i-1]
      res_max <- seq[i]
    }
    i <- i+1
  }
  if (conftest > 1) {cat("It was not possible to find significant variation with a sample size:",n,"\n")
  } else {
    if (res_min != res_max) {
      pas <- abs(mean*res_min-mean*res_max)/cut
      i <- cut; conftest <- 1
      while ((i > -1) & (conftest > conf) ) {
        pvalsup <- c() ; pvalinf <- c()
        moysup = (mean+i*pas)
        moyinf = (mean-i*pas)
        for (j in 1:iter) {
          pop_ref <- rnorm(n, mean = mean, sd = sd)
          pop_test_sup <- rnorm(n, mean = moysup, sd = sd)
          pop_test_inf <- rnorm(n, mean = moyinf, sd = sd)
          pvalsup <- c(pvalsup, t.test(pop_ref,pop_test_sup)$p.value)
          pvalinf <- c(pvalinf, t.test(pop_ref,pop_test_inf)$p.value)
        }
        conftest_sup <- length(pvalsup[pvalsup<(1-conf.level)])/iter
        conftest_inf <- length(pvalinf[pvalinf<(1-conf.level)])/iter
        conftest <- min(conftest_sup,conftest_inf)
        if (conftest < conf) {
          cat("This sample of ",n," individuals is discriminant with a p-value < ",1-conf,"\n\tfor ",conf.level*100,"% of the cases\n\tfor mean differences of +/- :\n")
          return((i-1)*pas)
        }
        i <- i-1
      }
    }}
}



#' Anticipate the sample size needed to detect a correlation (if there is one)
#'
#' @param x The values of my x variable
#' @param y The values of my y variable
#' @param iter Number of iterations of bootstrapping
#'
#' @return The minimum sample size needed to detect a significant correlation with the given data, if a significant correlation exists. The function returns this sample size.
#' 
#' @export
#' @importFrom stats cor.test
#'
#' @examples
#' x <- c(1,5,6,9,10)
#' y <- c(0,4.5,10,9,30)
#' plot(x,y,cex=2,pch=16)
#' cor.test(x,y)$p.value # No significant
#' cor_ech(x,y)
cor_ech <- function(x,y,iter=100) {
  nb_return <- c()
  for (i in 1:40 ) {
    nb_final <- c()
    for (i in iter) {
      pvalue = 1
      nb = 0
      while (pvalue > 0.05) {
        nb <- nb+1
        indices <- 1:length(x)
        indices_temp <- sample(indices ,nb,replace=T)
        x1 <- c(x,x[indices_temp])
        y1 <- c(y,y[indices_temp])
        pval <- c() ; j=0
        for (i in 1:iter) {
          indices <- 1:length(x1)
          indices_temp <- sample(indices ,length(indices), replace=T)
          if (length(unique(x1[indices_temp])) > 2) {
            j=j+1
            x_temp <- x1[indices_temp]
            y_temp  <- y1[indices_temp]
            #print(x_temp)
            pval[j] <- cor.test(x_temp,y_temp)$p.value
          }
        }
        pvalue <- quantile(pval,probs=0.95)
        if (nb > 1000) {pvalue<-0}
      }
      nb_final <- c(nb_final,nb)
    }
    max(nb_return,(quantile(nb_final,probs=0.95,names=F)+length(x)))-> nb_return
  }
  cat("To have a significant correlation,\nit would have required at least a number of values of : \n")
  return(nb_return)
}
