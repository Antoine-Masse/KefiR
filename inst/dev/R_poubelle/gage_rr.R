#' Gage R&R
#'
#' @param Var measures.
#' @param optr identifiers of the three operators.
#' @param app identifiers of the 10 pieces.
#' @param sigma tolerated variations (multiplication factor of the standard deviation).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{repeatability}{Variability due to the measuring equipment.}
#'   \item{reproducibility}{Variability due to operators.}
#'   \item{RR}{Combined effect of repeatability and reproducibility.}
#'   \item{Vp}{Part - variability due to the measured parts.}
#'   \item{Vt}{Total variability.}
#'   \item{part_Equipement}{Percentage of total variability due to measuring equipment.}
#'   \item{part_Operators}{Percentage of total variability due to operators.}
#'   \item{part_R_R}{Percentage of total variability due to the combination of repeatability and reproducibility.}
#'   \item{part_Pieces}{Percentage of total variability due to the measured parts.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(2)
#' data <- data.frame("Measures"=rep(rnorm(10,20,0.007),9)*abs(rnorm(90,1,0.025)),
#' "Operator"=rep(c("A","B","C"),30),
#' "Equipement"=rep(1:10,9))
#' boxplot(data$Measures~data$Equipement)
#' rr(data$Measures,data$Operator,data$Equipement)
rr = function(Var=c(),optr=c(),app=c(),sigma=5.15) {
  means = function(data=c(),parametres=c()) {
    temp = split(data,parametres,drop=TRUE) #drop permet de ne faire des calculs que sur les paramètres présentant des valeurs associées
    x = c() ; y = unique(parametres) ; z = c() ; w = c() ; resultat = list()
    for (i in (1:(length(y)))) {
      x = c(x,mean(temp[[i]],na.rm = TRUE))
      if (length(temp[[i]]) > 1) {
        z = c(z,sd(temp[[i]],na.rm = TRUE)) ;  w = c(w,sd(temp[[i]],na.rm = TRUE)*1.96/sqrt(length(temp[[i]])))}
      else {z = c(z,0);w = c(w,0)}
    }
    resultat$moyennes = x[order(y)] ; resultat$sd = z[order(y)] ; resultat$ic = w[order(y)] ; resultat$parametres = y[order(y)] #Intervalles de confiance
    return(resultat)
  }
  d2_constants = c(
    1.41,1.91,2.24,2.48,2.67,2.83,2.96,3.08,3.18,3.27,3.35,3.42,3.49,3.55 ,1.28,1.81,2.15,2.40,2.60,2.77,2.91,3.02,3.13,3.22,3.30,3.38,
    3.45,3.51, 1.23,1.77,2.12,2.38,2.58,2.75,2.89,3.01,3.11,3.21,3.29,3.37,3.43,3.5, 1.21,1.75,2.11,2.37,2.57,2.74,2.88,3.00,3.1,3.2,3.28,3.36,
    3.43,3.49, 1.19,1.74,2.1,2.36,2.56,2.78,2.87,2.99,3.1,3.19,3.28,3.36,3.42,3.49, 1.18,1.73,2.09,2.35,2.56,2.73,2.87,2.99,3.10,3.19,3.27,3.35,
    3.42,3.49, 1.17,1.73,2.09,2.35,2.55,2.72,2.87,2.99,3.1,3.19,3.27,3.35,3.42,3.48, 1.17,1.72,2.08,2.35,2.55,2.72,2.87,2.98,3.09,3.19,3.27,3.35,
    3.42,3.48, 1.16, 1.72, 2.08, 2.34, 2.55, 2.72, 2.86, 2.98, 3.09, 3.19, 3.27, 3.35, 3.42, 3.48, 1.16, 1.72, 2.08, 2.34, 2.55, 2.72, 2.86, 2.98, 3.09, 3.18, 3.27, 3.34, 3.42, 3.48, 1.15, 1.71, 2.08, 2.34, 2.55, 2.72, 2.86, 2.98, 3.09, 3.18, 3.27, 3.34, 3.41, 3.48, 1.15, 1.71, 2.07, 2.34, 2.55, 2.72, 2.85, 2.98, 3.09, 3.18, 3.27, 3.34, 3.41, 3.48, 1.15, 1.71, 2.07, 2.34, 2.55, 2.71, 2.85, 2.98, 3.09, 3.18, 3.27, 3.34, 3.41, 3.48, 1.15, 1.71,
    2.07, 2.34, 2.54, 2.71, 2.85, 2.98, 3.09, 3.18, 3.27, 3.34, 3.41, 3.48, 1.15, 1.71, 2.07, 2.34, 2.54, 2.71, 2.85, 2.98, 3.08, 3.18, 3.26, 3.34, 3.41, 3.48, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.97, 3.078, 3.173, 3.258, 3.336, 3.407, 3.472)
  d2_constants = matrix(d2_constants,ncol=14, nrow=16, byrow=T)
  # d2 is deduced from a grid with z and w in relation :
  # z : product of the number of operators times the number of pieces
  # w: number of repetitions per test
  colnames(d2_constants) = c(2:15) #W
  rownames(d2_constants) = c(1:15,">15") #z
  if ((length(Var)==0)|(length(optr)==0)|(length(app)==0)) {stop("Un des paramètres est un vecteur vide.\n")}
  operateurs = unique(optr)
  appareils = unique(app)
  n = length(appareils)
  r = length(Var[optr==operateurs[1]&app==appareils[1]]) #repetition
  o = length(operateurs)
  W = r-1 ; z = n*o
  if (z > 15) {z=16}
  d2 = d2_constants[z,W] #d2 = ss.cc.getd2(repetition) avec le package SixSigma
  # Calculation of repeatability
  moyenne_par_operateur <- c()
  etendue = c()
  for (op in operateurs) {
    for (i in appareils) {
      etendue = c(etendue,(max(Var[optr==op&app==i])-min(Var[optr==op&app==i])))
    }
    moyenne_par_operateur = c(moyenne_par_operateur,mean(Var[optr==op]))
  }
  #cat("expanses: ",etendue[1:10],"\n");cat("expanses: ",etendue [11:20],"\n"); cat("expanses: ",etendue [21:30],"\n")
  R = mean(etendue)
  repeatability = sigma*R/d2
  W = o-1 # W number of operators
  z = 1
  d2 = d2_constants[z,W] #d2 = ss.cc.getd2(repetition) avec le package SixSigma
  # Calculation of reproducibility
  X = max(moyenne_par_operateur)-min(moyenne_par_operateur) # étendue entre opérateurs
  #reproducibility = sqrt((sigma*X/d2)^2-((repetability^2)/(n*r)))
  A = (sigma*X/d2)^2
  B = ((repeatability^2)/(n*r))
  if (B>A) {stop("too high repeatability does not allow the calculation of reproducibility.\n")}
  reproductibility = sqrt(A-B)
  #cat("R: ",R,"\n");cat("Repeatability: ",repeatability,"\n");cat("d2: ",d2,"\n");cat("average_per_operator: ",average_per_operator,"\n");cat("X: ",X,"\n");cat("Reproducibility: ",reproductibility,"\n")
  R_R = sqrt(repeatability^2+reproductibility^2)
  moyenne_par_appareil = means(Var,app)
  Rp = max(moyenne_par_appareil$moyennes) - min(moyenne_par_appareil$moyennes) # Rp extended between the device with a max average measurement and the one with a min average measurement
  W = n-1 # W nombre d'appareils
  z = 1
  d2 = d2_constants[z,W] #d2 = ss.cc.getd2(repetition) avec le package SixSigma
  Vp = sigma*Rp/d2
  Vt = sqrt(R_R^2+Vp^2)
  RR=c()
  RR$repeatability = repeatability
  RR$reproductibility = reproductibility
  RR$RR = R_R
  RR$Vp = Vp
  RR$Vt = Vt
  RR$part_Equipement = repeatability/Vt*100
  RR$part_Operators = reproductibility/Vt*100
  RR$part_R_R = R_R/Vt*100
  RR$part_Pieces = Vp/Vt*100
  cat("\nWith Vp : part variability and Vt : total variability :\n")
  # Add total variability --> capability
  return(RR)
}
