#' Skewness test for normality
#'
#' @param x a numeric vector of data values.
#' @param nrepl the number of replications in Monte Carlo simulation.
#'
#' @return Performs skewness test for the composite hypothesis of normality, see, e.g., Shapiro, Wilk and Chen (1968).
#'
#' @examples
#' # This function is extracted from the deprecated library normtest
#' #skewness.norm.test(rnorm(100))
skewness.norm.test <- function(x, nrepl = 2000) {
    DNAME <- deparse(substitute(x))
    l <- 0
    n <- length(x)
    x1 <- sum(x)/n
    t <- sqrt(n) * (sum((x - x1)^3))/((sum((x - x1)^2))^(3/2))
    for (i in 1:nrepl) {
        z <- rnorm(n)
        z1 <- sum(z)/n
        T <- sqrt(n) * (sum((z - z1)^3))/((sum((z - z1)^2))^(3/2))
        if (abs(T) > abs(t))
            l = l + 1
    }
    p.value <- l/nrepl
    RVAL <- list(statistic = c(T = t), p.value = p.value, method = "Skewness test for normality",
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)

}


#' Kurtosis test for normality
#'
#' @param x a numeric vector of data values.
#' @param nrepl the number of replications in Monte Carlo simulation.
#'
#' @return Performs kurtosis test for the composite hypothesis of normality, see, e.g., Shapiro, Wilk and Chen (1968).
#'
#' @examples
#' # This function is extracted from the deprecated library normtest
#' #kurtosis.norm.test(rnorm(100))
kurtosis.norm.test <- function(x, nrepl = 2000) {
    DNAME <- deparse(substitute(x))
    l <- 0
    n <- length(x)
    x1 <- sum(x)/n
    t <- n * (sum((x - x1)^4))/((sum((x - x1)^2))^2)
    for (i in 1:nrepl) {
        z <- rnorm(n)
        z1 <- sum(z)/n
        T <- n * (sum((z - z1)^4))/((sum((z - z1)^2))^2)
        if (abs(T - 3) > abs(t - 3))
            l = l + 1
    }
    p.value <- l/nrepl
    RVAL <- list(statistic = c(T = t), p.value = p.value, method = "Kurtosis test for normality",
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)}


#' Jarque-Bera test for normality
#'
#' @param x a numeric vector of data values.
#' @param nrepl the number of replications in Monte Carlo simulation.
#'
#' @return Performs Jarque--Bera test for the composite hypothesis of normality, see Jarque and Bera (1987)..
#'
#' @examples
#' # This function is extracted from the deprecated library normtest
#' # jb.norm.test(rnorm(100))
#' # jb.norm.test(abs(runif(100,-2,5)))
jb.norm.test<-function(x, nrepl=2000) {
	DNAME <- deparse(substitute(x))
	l<-0
	n<-length(x)
	x1<-sum(x)/n
	b1<-sqrt(n)*(sum((x-x1)^3))/((sum((x-x1)^2))^(3/2))
	b2<-n*(sum((x-x1)^4))/((sum((x-x1)^2))^2)
	t<-(n/6)*((b1)^2 + ((b2-3)^2)/4)
	for(i in 1:nrepl)
	{
		z<-rnorm(n)
		z1<-sum(z)/n
		a1<-sqrt(n)*(sum((z-z1)^3))/((sum((z-z1)^2))^(3/2))
		a2<-n*(sum((z-z1)^4))/((sum((z-z1)^2))^2)
		T<-(n/6)*((a1)^2 + ((a2-3)^2)/4)
		if (T>t) l=l+1
	}
	p.value<-l/nrepl
	RVAL<-list(statistic=c(JB=t), p.value=p.value, method="Jarque-Bera test for normality",data.name = DNAME)
	class(RVAL)<-"htest"
	return(RVAL)
}

