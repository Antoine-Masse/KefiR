#' Calculate an pairwise difference between samples by bootstrap without encountering a central-limit problem
#'
#' @param x numerical vector
#' @param g categorical vector
#' @param iter number of iteration
#' @param mu comparison criterion, by default 'mean' (arithmetic mean), otherwise 'meanbp' (moving average per iteration: an automatic compromise between mean and median), 'median' or 'sd'.
#' @param paired Logical. If TRUE, performs paired bootstrap tests.
#' @param debug logical, if TRUE enables debug mode for troubleshooting.
#' @param verbose Logical. If TRUE, prints verbose output.
#' @param return Logical. If TRUE, returns the results.
#'
#' @return This function returns an array of p-values like the other pairwise functions.
#' @return It calculates the number of differences between bootstrap samples that do not include 0 (all higher or all lower).
#' @return This number of differences divided by the number of iterations (iter) gives the maximum percentage of times that a convergent difference (all higher or all lower) is found: a confidence.
#' @return This value subtracted from 1 gives an equivalent of the p-value whose precision depends on the number of iterations.
#' @return Note: using the meanbp criterion is more relevant, it allows a compromise between mean and median by avoiding leverage effects.
#' @export
#'
#' @examples
#' # Example 1
#' data(iris)
#' pairwise.boot(iris[,2],iris$Species)
#' # Example 2 by using pairwise(type=="boot")
#' data(mtcars)
#' pairwise(mtcars$mpg[mtcars$carb<=4],mtcars$carb[mtcars$carb<=4],type="boot")
pairwise.boot <- function(x, g, iter = 500, mu = "meanbp", debug = FALSE,
                         paired = FALSE, verbose = F, return = T) {
	.dbg(NULL,"Exécution de pairwise.boot().",debug=debug)
  g <- as.factor(g)

  # Normaliser mu : si c'est une fonction, la convertir en string
  if (is.function(mu)) {
    mu_name <- deparse(substitute(mu))
    # Si deparse ne donne pas un nom simple, essayer de déduire
    if (identical(mu, mean)) mu <- "mean"
    else if (identical(mu, median)) mu <- "median"
    else if (identical(mu, sd)) mu <- "sd"
    else if (exists("meanbp") && identical(mu, meanbp)) mu <- "meanbp"
    else mu <- "mean"  # défaut
  }

  # Réordonner les groupes par statistique (cohérence avec pairwise(type="boot"))
  # Utilise la même statistique que celle utilisée pour le bootstrap
  stat_fn <- switch(mu,
    "mean" = mean,
    "median" = median,
    "meanbp" = meanbp,
    "sd" = sd,
    mean  # défaut
  )
  stat_values <- by(x, g, stat_fn)
  ordered_levels <- levels(g)[order(stat_values)]
  g <- factor(g, levels = ordered_levels)

	unique_g <- levels(g)
  # Vérification des données pour `paired`
  if (paired) {
    group_lengths <- tapply(x, g, length)
    if (length(unique(group_lengths)) != 1) {
      .exit("For paired data, all groups must have the same number of observations.",
	  NULL, verbose=verbose, return=return)
    }
  }
	#.dbg(NULL,"bim.",debug=debug)
	mymat <- matrix(rep(NA, (length(unique_g) - 1)^2),
	                ncol = (length(unique_g) - 1),
	                nrow = (length(unique_g) - 1))
	rownames(mymat) <- unique_g[2:length(unique_g)] ; colnames(mymat) <- unique_g[1:(length(unique_g)-1)]
	#.dbg(NULL,"bam.",debug=debug)	
	for (i in 1:(length(unique_g)-1)) {
		for (j in (i+1):length(unique_g)) {
			ech1 <- x[g==unique_g[i]]
			ech2 <- x[g==unique_g[j]]
			diff <- numeric(iter)
			for (k in 1:iter) {
				if (paired) {
				  # Rééchantillonnage apparié
				  sampled_indices <- sample(seq_along(ech1), replace = TRUE)
				  ech1_temp <- ech1[sampled_indices]
				  ech2_temp <- ech2[sampled_indices]
				} else {
				  # Rééchantillonnage indépendant
				  ech1_temp <- sample(ech1, replace = TRUE)
				  ech2_temp <- sample(ech2, replace = TRUE)
				}
				# Calcul du descripteur
				if (mu == "meanbp") {
				  diff[k] <- meanbp(ech1_temp) - meanbp(ech2_temp)
				} else if (mu == "mean") {
				  diff[k] <- mean(ech1_temp) - mean(ech2_temp)
				} else if (mu == "median") {
				  diff[k] <- median(ech1_temp) - median(ech2_temp)
				} else if (mu == "sd") {
				  diff[k] <- sd(ech1_temp) - sd(ech2_temp)
				} else {
				  .exit("Invalid value for 'mu'. Choose from 'meanbp', 'mean', 'median', or 'sd'.",
				  NULL,verbose=verbose, return=return)
				}
			}
			# Calcul de la p-value bilatérale
			length(diff[diff>0]) -> sup
			length(diff[diff<0]) -> inf
			if (sup>inf) {
				pv <- 1-(sup/iter)
			}else {pv <- 1-(inf/iter)}
			mymat[(j-1),i] <- pv*2 # *2 because bilateral
		}
	}
	#.dbg(NULL,"boum.",debug=debug)	
	result <- list() ; result$p.value <- mymat
	return(result)
}


