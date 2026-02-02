#' Perform One-Factor Analysis
#'
#' This function performs statistical analysis for one grouping factor.
#' It handles both normal and non-normal data, adjusting the statistical tests accordingly.
#'
#' @param x Numeric vector. Dependent variable to be analyzed.
#' @param g Factor. Grouping variable for the analysis.
#' @param formula Optional. A formula specifying the model (e.g., x ~ g).
#' @param data Optional. A data frame containing variables specified in the formula.
#' @param paired Logical. Indicates whether data are paired (default is FALSE).
#' @param id Optional. Identifier for paired data.
#' @param alpha Numeric. Significance level for hypothesis testing (default is 0.05).
#' @param k Optional. Counter for verbose messages.
#' @param debug Logical. If TRUE, detailed debugging messages are displayed (default is FALSE).
#' @param verbose Logical. If TRUE, provides detailed output during the analysis (default is FALSE).
#'
#' @details
#' The function performs the following:
#' - Adjusts the significance level using Sidak's correction.
#' - Tests normality of data using Shapiro-Wilk or Jarque-Bera tests, based on sample size.
#' - Applies appropriate statistical tests: t-test, Mann-Whitney, ANOVA, Kruskal-Wallis.
#' - Handles cases with equal or unequal variances and non-normal data.
#'
#' @importFrom fda.usc fanova.hetero
#' @importFrom agricolae kurtosis skewness
#' @importFrom WRS2 med1way t1way
#' @importFrom onewaytests bf.test
#'
#' @examples
#' # Example with normal data
#' set.seed(123)
#' x <- c(rnorm(50, mean = 5), rnorm(50, mean = 6))
#' g <- factor(rep(1:2, each = 50))
#' .one_factor_analysis(x, g, verbose = TRUE)
#'
#' # Example with non-normal data
#' x <- c(rpois(50, lambda = 5), rpois(50, lambda = 7))
#' g <- factor(rep(1:2, each = 50))
#' .one_factor_analysis(x, g, verbose = TRUE)
#'
#' @return
#' Returns results of the selected statistical tests, including:
#' - p-values for hypothesis testing.
#' - Results of normality checks and variance assessments.
#'
#' @seealso
#' - [fda.usc::fanova.hetero()]
#' - [agricolae::kurtosis()]
#' - [lawstat::levene.test()]
#' - [WRS2::med1way()]
#' - [onewaytests::bf.test()]
#'
#' @export
.one_factor_analysis <- function(x=NULL,g=NULL,formula=NULL,data=NULL,
	paired = FALSE, id = NULL,alpha = 0.05,return=TRUE,
	k=NULL,code=FALSE,debug=FALSE,verbose=FALSE, boot = TRUE, silent=TRUE) {

  ################################################################
  #
  #
  #		Si paired = TRUE et 2 groupes (k  = 2)
  #
  #
  ################################################################
  ## --- Vérifs d’intégrité pour le mode apparié ------------------------------
  if (isTRUE(paired)) {
    # Si identifiants fourni,le contrôle et réaligner les données dessus
    if(!is.null(id)) {
      print("id")
      print(id)
      g2 <- droplevels(factor(g))
      #if (nlevels(g2) != 2L) {
      #  .exit("Paired mode expects exactly two levels in 'g'.",
      #        "Le mode apparié requiert exactement deux modalités dans 'g'.",
      #        verbose=verbose, return=return)
      #}
      dfc <- data.frame(id = data[id], g = g2, x = as.numeric(x))
      #dfc <- dfc[complete.cases(dfc), ]
      print(head(dfc))
      tab <- table(dfc$id, dfc$g)              # matrice (#id x 2)
      # 1) chaque id doit être présent dans les 2 modalités
      miss <- rownames(tab)[rowSums(tab > 0) < 2L]
      if (length(miss)) {
        .exit(
          "Some ids are incomplete (no match for each modality of 'g').",
          "Certains identifiants de id sont incomplets (pas de correspondance pour chaque modalité de 'g').",
          verbose = verbose, return = return
        )
      }
      # 2) (option strict) exactement une obs par id et par modalité
      if (any(tab != 1L)) {
        .exit(
          "Each id must have exactly one observation per level of 'g' (no duplicates).",
          "Chaque identifiant doit avoir exactement une observation par modalité de 'g' (pas de doublons).",
          verbose = verbose, return = return
        )
      }
      # Réaligner
      k <- .vbse("Alignment by id for paired (2 levels).",
            "Alignement par id pour les données appariées.",
            verbose = verbose, k = k, cpt = "on") -> k
      temp <- .align_pairs(x,data[id],g)
      x <- temp$x
      g <- temp$g
    }
  #
    #		Si paired = TRUE et identifiants fournis...
    if (nlevels(g)==2) {
      diff <- x[g==levels(g)[1]]-x[g==levels(g)[2]]
      check_normality <- .normality(diff, paired=paired, k=k, verbose=verbose, alpha=alpha)
      k <- check_normality[[2]] ; check_normality <- check_normality[[1]]
      check_variance_equal <- TRUE
      return(check <- list(x,g,check_normality,check_variance_equal,k))
    }
  } # Fin du scénario paired sachant que l'alignement si n_g > 2 (k) passe en .multi_factor_analysis...

  #==================================================
  #         Initiation des variables de contrôle
  #==================================================
  check_discret <- discret.test(x)
  check_normality <- c()
  check_number <- c()
  check_variance_equal <-  c()

  print("c")
  #==================================================
  #         Contrôle du nombre de groupes
  #==================================================
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Code
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (code==TRUE){cat("length(unique(g))#2)\n")}

  check_number <- length(unique(g))
  if (check_number==2) { 							# 2 categories
    if (paired==TRUE) {
      k <- .vbse("Two conditions.",
                 "Deux conditions.",
                 verbose = verbose, k = k, cpt="on")
    } else{
      k <- .vbse("Two groups.",
                 "Deux groupes.",
                 verbose = verbose, k = k, cpt="on")
    }
  } else if (check_number>2) {
    if (paired==TRUE) {
      k <- .vbse("More than two conditions.",
        "Plus de deux conditions.",
        verbose = verbose, k = k)
    } else{
      k <- .vbse("More than two groups.",
                 "Plus de deux groupes.",
                 verbose = verbose, k = k)
    }
  } else {
    .exit("Not enough analyzable groups.","Pas assez de groupes analysables.")
  }

  #==================================================
  # Contrôle de l'indépendance des observations.
  #==================================================
  .dbg("Check of observation independence.",
       "Contrôle de l'indépendance des observations.",debug=debug)
  if (paired == FALSE) {
    # Données non appariées
    k <- .vbse("Control for independence: Ensure that observations between groups are independent...\n\t(no repeated measures, no clustering, and no carry-over effects).",
               "Contrôle de l'indépendance : vérifiez que les observations entre les groupes sont indépendantes\n\t(pas de mesures répétées, pas d'effet cluster, et pas d'effets report [influence sur l'ordre des mesures]).",
               verbose = verbose, k = k, cpt = "on")
  } else {
    # Données appariées
    k <- .vbse("Control for paired independence: Ensure that the dependence structure is only due to the intended pairing (e.g., same subject before/after).\n\tCheck for no additional clustering or carry-over effects between paired observations.",
               "Contrôle de l'indépendance pour données appariées :\n\tVérifiez que la structure de dépendance provienne uniquement de l'appariement prévu (ex. : même sujet avant/après).\n\tAssurez-vous qu'il n'y a pas de clustering ou d'effets report supplémentaires entre les observations appariées.",
               verbose = verbose, k = k, cpt = "on")
  }
  #==================================================
  # Détection des Outlier en X
  #==================================================
  if (paired == FALSE) {
    #---------------------------------
    #  Gestion des données extrêmes si non appariés
    #---------------------------------
    outliers_base <- as.vector(unlist(
      by(x, g, function(z) {
        length(which(identify_outliers(data.frame(z))$is.outlier == TRUE)) / length(z)
      })
    ))
    outliers_extrem <- as.vector(unlist(
      by(x, g, function(z) {
        length(which(identify_outliers(data.frame(z))$is.extreme == TRUE)) / length(z)
      })
    ))
    max_outliers_base  <- max(outliers_base)
    max_outliers_extrem <- max(outliers_extrem)
    # Messages gradués selon la sévérité
    if (max_outliers_extrem > 0) {
      k <- .vbse(
        paste0("Extreme outliers detected [identify_outliers() from {rstatix}].\n\t",
               "Maximum proportion of extreme outliers in a group: ",
               round(max_outliers_extrem * 100, 1), "%.\n\t",
               "Caution: these values may strongly influence parametric tests."),
        paste0("Outliers extrêmes détectés [identify_outliers() de {rstatix}].\n\t",
               "Proportion maximale d'outliers extrêmes dans un groupe : ",
               round(max_outliers_extrem * 100, 1), "%.\n\t",
               "Attention : ces valeurs peuvent fortement influencer les tests paramétriques."),
        verbose = verbose, k = k, cpt = "on"
      )
    } else if (max_outliers_base > 0) {
      k <- .vbse(
        paste0("Outliers detected [identify_outliers() from {rstatix}].\n\t",
               "Maximum proportion of outliers in a group: ",
               round(max_outliers_base * 100, 1), "%.\n\t",
               "These values are monitored but not excessive."),
        paste0("Outliers détectés [identify_outliers() de {rstatix}].\n\t",
               "Proportion maximale d'outliers dans un groupe : ",
               round(max_outliers_base * 100, 1), "%.\n\t",
               "Ces valeurs sont à surveiller mais non excessives."),
        verbose = verbose, k = k, cpt = "on"
      )
    }
  } else {
    #---------------------------------
    # Données appariées — contrôle des outliers sur les différences intra-sujet
    lv <- levels(droplevels(factor(g)))
    #---------------------------------
    if (length(lv) == 2L) {
      dif <- x[g == lv[2]] - x[g == lv[1]]
      io  <- suppressWarnings(identify_outliers(data.frame(dif)))

      prop_ext <- if ("is.extreme" %in% names(io)) mean(io$is.extreme, na.rm = TRUE) else 0
      prop_out <- if ("is.outlier" %in% names(io)) mean(io$is.outlier, na.rm = TRUE) else 0

      if (prop_ext > 0) {
        k <- .vbse(
          paste0("Extreme outliers in paired differences [identify_outliers() from {rstatix}].\n\t",
                 "Proportion of extreme outliers: ", round(prop_ext * 100, 1), "%.\n\t",
                 "These pairs may strongly influence paired parametric tests."),
          paste0("Outliers extrêmes dans les différences appariées [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d’outliers extrêmes : ", round(prop_ext * 100, 1), "%.\n\t",
                 "Ces paires peuvent fortement influencer les tests paramétriques appariés."),
          verbose = verbose, k = k, cpt = "on"
        )
      } else if (prop_out > 0) {
        k <- .vbse(
          paste0("Outliers in paired differences [identify_outliers() from {rstatix}].\n\t",
                 "Proportion of outliers: ", round(prop_out * 100, 1), "%.\n\t",
                 "Monitored but not extreme."),
          paste0("Outliers dans les différences appariées [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d’outliers : ", round(prop_out * 100, 1), "%.\n\t",
                 "À surveiller mais non extrêmes."),
          verbose = verbose, k = k, cpt = "on"
        )
      }
    } else {
      k <- .vbse(
        "Paired data with more than two levels: outlier check on differences is not defined here.",
        "Données appariées avec plus de deux modalités : contrôle des outliers sur les différences non défini ici.",
        verbose = verbose, k = k, cpt = "off"
      )
    }
  }
	  #==================================================
	  #
    #   Contrôle de la normalité
    #
	  #		On notera que formula est ignorée par cette fonction (à dev si nécessaire)
	  #
	  #==================================================
    .dbg("Normality check with .normality().",
         "Contrôle de normalité avec .normality().",debug=debug)
    skew <- function(vector) {return(abs(skewness(vector)))}
    skew2 <- function(vector) {return(skewness.norm.test(vector)$p.value)}
    kurto <- function(vector) {if (is.na(abs(kurtosis(vector)))){return(10)} ; return(abs(kurtosis(vector)))}
    kurto2 <- function(vector) {return(kurtosis.norm.test(vector)$p.value)}
    ##########################
    # Correction de Sidak
    ##########################
    # Sidak devient trop conservatif au-delà de 10 groupes...
    # Bender, R., & Lange, S. (2001). Adjusting for multiple testing—when and how?  DOI: 10.1016/S0895-4356(00)00314-0
    # Il est justifié de passer ensuite au bootstrap
    # Méthode FDR la correction par FDR (False Discovery Rate) – notamment via la méthode Benjamini-Hochberg (BH)
    # Murray, E. J., Berrie, L., & Matthews, J. N. S. (2021). Understanding multiplicity in clinical trials: the impact of multiple endpoints in the interpretation of treatment effects.
    # Attention : FDR est dans une démarche exploratoire ou descriptive, pas une validation réglementaire.
    .dbg("Sidak's correction rather than Bonferroni for independent tests (conservative).",
         "Correction de Sidak plutôt que Bonferroni pour les tests indépendants (conservateur).",debug=debug)
    pval <- 1-(1-alpha)^(1/length(unique(g)))
    #-*************************
    #-    code==TRUE
    #-*************************
    if (code==TRUE){
      cat(paste0("alpha1 <- 1-(1-",alpha,")^(1/length(unique(g))) # Sidak correction of alpha for test repetitions (like Shapiro, Skewness...)\n"))
      if (min(by(x,g,length))<=100) {
        cat("by(x,g,shapiro.test)#1) Control of the normality of small samples (<100)\n")
      }
      if (any((by(x,g,length)<=1000)&(by(x,g,length)>100))) {
        cat("#1A)\nby(x,g,jb.norm.test)#1B) Control of the normality of big samples (<1000)\n")
      }
      if (max(by(x,g,length))>1000) {
        cat("#1) Central limit theory: enough values for some samples to not have to check normality.\n")
      }
    }
    #-*************************
    #-    code==TRUE
    #-*************************
    if (code==TRUE){cat("library(agricolae)\nby(x,g,skewness)\nby(x,g,skewness.norm.test)\nby(x,g,kurtosis)\nby(x,g,kurtosis.norm.test)\nby(x,g,length)#3)\n")}
    if (check_number==2) { 							# 2 categories
      # Contrôle de la normalité des groupes
      check_normality <- .normality(x,g, alpha=alpha, k=k, verbose=verbose)
    } else {
      # Contrôle de la normalité des résidus
      myaov <- aov(x~g)
      n_residu <- length(x)
      check_normality <- .normality(myaov$residuals, alpha=alpha, k=k, verbose=verbose) #pvals <- .normality(x,g)
    }
    k <- check_normality[[2]] ; check_normality <- check_normality[[1]]
    #==================================================
    #
    #   Contrôle de la variance
    #
    #==================================================
    .dbg("Variance check with .variance().",
         "Contrôle de la variance avec .variance().",debug=debug)
    variance_temp <- .variance(x, g, check_normality=check_normality, alpha=alpha, paired=paired,
                  debug = debug, verbose=verbose, k=k, code=code)
    check_variance_equal <- variance_temp[[1]]
    k <- variance_temp[[2]]
    #=============================
    # Contrôle de l'équilibrage
    #=============================
    if (check_number>2) {
      .dbg("Balance check.",
           "Contrôle de l'équilibrage.",debug=debug)
      if (check_variance_equal == TRUE) {
        fm <- formula(x~g)
        mya <- suppressWarnings(aov(fm, data = data.frame(x, g)))
        #=============================
        # Contrôle de l'équilibrage
        #=============================
        # McKissick, S. K. (2017). Effectiveness of Pearson’s SuccessMaker Mathematics for Students with Disabilities. Journal of the American Academy of Special Education Professionals (JAASEP), 21(4), 103-135. Retrieved from https://files.eric.ed.gov/fulltext/EJ1129643.pdf
        k <- .vbse(
          "Procedure: assess data balance with n_max/n_min (size imbalance) and F_max = Var_max/Var_min (variance heterogeneity).",
          "Contrôle de l’équilibre par n_max/n_min (déséquilibre des tailles) et F_max = Var_max/Var_min (ratio des variances).",
          verbose = verbose, k = k, cpt = "on")
        tab_n    <- table(g)
        ratio_ng <- max(tab_n) / min(tab_n) # nmax/nmix # Vérifier que < 4:1.
        # Variances par groupe avec by()
        var_groupes <- by(x, g, var, na.rm = TRUE)
        Fmax <- max(var_groupes) / min(var_groupes) # VARmax/VARmix # Vérifier que < 4:1.
        if (ratio_ng <= 4 && Fmax <= 10) {
          # On reste vers une ANOVA classique
          k <- .vbse(
            "Decision: keep classical ANOVA (sizes balanced and F_max within tolerance).",
            "Décision : on conserve l’ANOVA classique (tailles équilibrées et F_max dans la tolérance).",
            verbose = verbose, k = k, cpt="off")
        } else {
          # On repart vers Welch
          check_variance_equal <- FALSE
          k <- .vbse(
            "Decision: classical ANOVA not robust (large size imbalance and/or high F_max) → switch to Welch/robust.",
            "Décision : ANOVA classique non robuste (déséquilibre fort et/ou F_max élevé) → bascule vers Welch/robuste.",
            verbose = verbose, k = k, cpt="off"
          )
        }
      }
    }
    #==================================================
    #
    #   Contrôle d'un éventuel retour vers paramétrique
    #
    #==================================================
    .dbg("Check for a possible return to parametric methods.",
         "Contrôle d'un éventuel retour vers paramétrique.",debug=debug)
  	##########################
  	#
  	#		Auto-analyse pour envisager un retour vers la situation paramétrique (Kurtosis & Skweness...)
  	#
  	##########################
	  # |Skewness| < 1 et |Kurtosis| < 1.5 sont les seuils recommandés pour considérer une distribution comme "approximativement normale", selon Kline (2011).
	  # Kline, R. B. (2011). Principles and practice of structural equation modeling (4th ed.). New York, NY: Guilford Press.
	  # Normalité variables indicates acceptable normality for most measures, based on Kline's (2023) guidelines, where skewness values should ideally be within ±3 and kurtosis values within ±10.
	  # Permissivité extrême : skewness > 2 et kurtosis > 7 Blanca et al. (2018)
	  # Blanca, M. J., Alarcón, R., Arnau, J., Bono, R., & Bendayan, R. (2017). Non-normal data: Is ANOVA still a valid option? Psicothema, 29(4), 552-557. https://diposit.ub.edu/dspace/bitstream/2445/122126/1/671797.pdf#:~\:text=sought%20to%20answer%20the%20following,In%20addition
    # Une violation de la normalité a peu d'effet sur l'ANOVA si variances homogènes.Caldwell (2022)
    # Caldwell, A. (2022, March 31). Chapter 12 Violations of Assumptions. In Power Analysis with Superpower. https://aaroncaldwell.us/SuperpowerBook/violations-of-assumptions.html

    auto_ku_sk <- function(x, g, check_normality = TRUE, alpha = 0.05, verbose = TRUE, k = NULL) {
      # Coercitions
      if (is.data.frame(x)) x <- x[[1]]
      if (is.matrix(x))     x <- x[,1,drop=TRUE]
      if (is.list(x))       x <- unlist(x, use.names = FALSE)
      x <- as.vector(x)
      if (!is.numeric(x)) stop("auto_ku_sk() attend un vecteur numérique 'x'.")

      if (is.data.frame(g)) g <- g[[1]]
      if (is.matrix(g))     g <- g[,1,drop=TRUE]
      if (is.list(g))       g <- unlist(g, use.names = FALSE)
      g <- droplevels(factor(g))
      lev <- levels(g)
      ng  <- nlevels(g)

      # Normalité (tolérante)
      if (ng == 2L) {
        tmp <- .normality(x, g, alpha = alpha, tolerance = "extrem", k = k, verbose = verbose)
        check_normality <- tmp[[1]]; k <- tmp[[2]]
      } else {
        fm <- x ~ g
        myaov <- suppressWarnings(aov(fm, data = data.frame(x = x, g = g)))
        tmp <- .normality(myaov$residuals, alpha = alpha, tolerance = "extrem", k = k, verbose = verbose)
        check_normality <- tmp[[1]]; k <- tmp[[2]]
      }
      # Outliers (par groupe) si normalité non rejetée
      if (isTRUE(check_normality)) {
        if (ng == 2L) { seuil_base <- 0.05; seuil_extrem <- 0.01 } else { seuil_base <- 0.10; seuil_extrem <- 0.03 }

        by_out <- by(x, g, function(z){
          tab <- rstatix::identify_outliers(data.frame(z))
          n   <- length(z)
          if (nrow(tab) == 0) return(c(out = 0, ext = 0))
          prop_out <- sum(tab$is.outlier %in% TRUE) / n
          prop_ext <- sum(tab$is.extreme %in% TRUE) / n
          c(out = prop_out, ext = prop_ext)
        })
        mat <- t(sapply(by_out, identity))
        outliers_base   <- max(mat[,"out"], na.rm = TRUE)
        outliers_extrem <- max(mat[,"ext"], na.rm = TRUE)

        viol_grp <- names(which(mat[,"out"] >= seuil_base | mat[,"ext"] >= seuil_extrem))
        n_viol   <- length(viol_grp)

        # --- Affichage conditionnel demandé ---
        if ((outliers_base >= seuil_base) || (outliers_extrem >= seuil_extrem)) {
          # Problème d'outliers : on signale (dans tous les cas)
          if (ng == 2L) {
            # Deux groupes : on affiche SEULEMENT s'il y a des groupes à citer
            if (n_viol > 0) {
              liste <- paste(viol_grp, collapse = ", ")
              ang <- paste0(
                "Outlier issue in group(s): ", liste, ". ",
                "\n\tMax proportions — outliers = ", round(100*outliers_base,1), "%, ",
                "extreme = ", round(100*outliers_extrem,1), "%. ",
                "\n\tPrefer a robust or non-parametric approach."
              )
              fr  <- paste0(
                "Problème d’outliers dans le(s) groupe(s) : ", liste, ". ",
                "\n\tProportions maximales — outliers = ", round(100*outliers_base,1), "%, ",
                "extrêmes = ", round(100*outliers_extrem,1), "%. ",
                "\n\tPréférer une approche robuste ou non paramétrique."
              )
              k <- .vbse(ang, fr, verbose = verbose, k = k)
            }
          } else {
            # ng > 2 : message synthétique
            ang <- paste0(
              "Outlier issues in ", n_viol, " group(s). ",
              "Max proportions — outliers = ", round(100*outliers_base,1), "%, ",
              "extreme = ", round(100*outliers_extrem,1), "%. ",
              "Prefer a robust or non-parametric approach."
            )
            fr  <- paste0(
              "Problèmes d’outliers dans ", n_viol, " groupe(s). ",
              "Proportions maximales — outliers = ", round(100*outliers_base,1), "%, ",
              "extrêmes = ", round(100*outliers_extrem,1), "%. ",
              "Préférer une approche robuste ou non paramétrique."
            )
            k <- .vbse(ang, fr, verbose = verbose, k = k)
          }
          check_normality <- FALSE
        } else {
          # Niveaux acceptables : on N'AFFICHE RIEN si ng == 2 (demande utilisateur).
          if (ng > 2L) {
            ang <- paste0(
              "Outlier check: low proportions — outliers = ", round(100*outliers_base,1),
              "%, extreme = ", round(100*outliers_extrem,1), "%. ",
              "Parametric approach can be maintained."
            )
            fr  <- paste0(
              "Contrôle des outliers : proportions faibles — outliers = ", round(100*outliers_base,1),
              "%, extrêmes = ", round(100*outliers_extrem,1), "%. ",
              "L’approche paramétrique peut être conservée."
            )
            k <- .vbse(ang, fr, verbose = verbose, k = k)
          }
        }
      }

      # Explication complémentaire si rejet
      if (!isTRUE(check_normality) && ng > 2L) {
        tmp2 <- .normality(x, g, alpha = alpha, tolerance = "extrem", k = k, verbose = verbose, cpt = "off")
        k <- tmp2[[2]]
      }

      # Sortie
      list(check_normality, k)
    }
    if ((check_normality==FALSE)&(check_variance_equal == TRUE)) {
      temp <- auto_ku_sk(x, g, check_normality=check_normality,
                             alpha=alpha, verbose = verbose, k = k)
      check_normality <- temp[[1]]
      k <- temp[[2]]
    }

	  ## --------------------------------------------------------------------------
	  ################################################################
	  #
	  #
	  #		NON-NORMAL
	  #
	  #
	  ################################################################
	  if (!check_normality) { # NON-NORMAL
      #========================================================
      # Contrôle de la distribution pour évaluer la fiabilité du test de Wilcoxon-Mann-Witney...
	    #========================================================
			sd_cr <- by(x,g,sd,na.rm=T) ; median_cr <- by(x,g,median,na.rm=T)
			data_cr <- x
			for (i in names(median_cr)) {
				data_cr[g==i] <- (data_cr[g==i]-median_cr[which(names(median_cr)==i)])/sd_cr[which(names(sd_cr)==i)]
			}
			temp <- pairwise(data_cr,g,type="ks",silent=silent,boot=FALSE)
			ks_result <- min(unlist(temp$p.value),na.rm=T)

			#========================================================
			#			NON-NORMAL		2 categories
			#========================================================
			if (check_number==2) { 							# 2 categories

			    if  (check_discret == TRUE)  {
					#========================================================
					#			NON-NORMAL		2 categories		Non acceptable for t.test()		Discret
					#========================================================
						ang <- paste0("The number of unique values [length() & unique()] suggests the presence of discrete data :\n\t",round(length(unique(x))/length(x),3)*100,"% of unique values.")
						fr <- paste0("Le nombre de valeurs uniques [length() & unique()] suggère la présence de données discrètes :\n\t", round(length(unique(x))/length(x),3)*100, "% de valeurs uniques.")
						k <- .vbse(ang,fr,verbose = verbose, k = k, cpt="on")
						pvals <- mood.test(x[g==unique(g)[1]],x[g==unique(g)[2]])$p.value
						if (pvals <= alpha) {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are different :\n\t==> Data need to be centered on the mediane for ansari.test(). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les médianes sont différentes :\n\t==> Les données doivent être centrées sur la médiane pour ansari.test(). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, k = k)
						  by(x,g,function(x){return(x-median(x))})->cent_med
						  pvals <- ansari.test(unlist(cent_med[1]),unlist(cent_med[2]))$p.value
						} else {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are the same (no need to center them for ansari.test). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les médianes sont identiques (pas besoin de les centrer pour ansari.test). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, k = k)
						  pvals <- ansari.test(x[g==unique(g)[1]],x[g==unique(g)[2]])$p.value
					  }
						if (pvals < alpha) {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data do not have the same variance. p-value: ",.format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les données n'ont pas la même variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, k = k)
						} else {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data have the same variance. p-value: ",.format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les données ont la même variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, k = k)
						}
  						#========================================================
  						#			NON-NORMAL		2 categories		Non acceptable for t.test() ==> Wilcoxon
  						#========================================================
					} else if ((check_variance_equal==FALSE) & (check_discret == FALSE)) {

  					#pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
  					#if (pvals <= alpha) {
  				#		ang <- paste0("Wilcoxon-Mann-Whitney test [wilcox.test()] -\n\t=> Significant differences between groups. p-value: ", .format_pval(pvals))
  			#			fr <- paste0("Test de Wilcoxon-Mann-Whitney [wilcox.test()] -\n\t=>  Différences significatives entre les groupes. p-value : ", .format_pval(pvals))
  			#			k <- .vbse(ang, fr, verbose = verbose, k = k)
  		#			} else {
  		#				ang <- paste0("Wilcoxon-Mann-Whitney test [wilcox.test()] -\n\t=> Non-significant differences between groups. p-value: ", .format_pval(pvals))
  		#				fr <- paste0("Test de Wilcoxon-Mann-Whitney [wilcox.test()] -\n\t=> Différences non significatives entre les groupes. p-value : ", .format_pval(pvals))
  		#				k <- .vbse(ang, fr, verbose = verbose, k = k)
  		#			}
					}
			} else { 											# > 2 categories
			  ###################################################
			  #			NON-NORMAL		k>2 categories
			  ###################################################
			    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			    #   CODE
			    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					if (code==TRUE){cat("kruskal.test(x,g) #6)\n")}


					##################
					#	Test de Fligner
					##################
					pvals3 <- fligner.test(x,g)$p.value
					if (is.na(pvals3)) {
					  .exit("Fligner-Killeen test [fligner.test()] - Error, return NA.",
						"Test de Fligner-Killeen [fligner.test()] - Erreur, retourne NA.",
						verbose=TRUE, return=TRUE)
					}
					##################
					#	Test de Kruskal
					##################
					pvals <- kruskal.test(x,g)$p.value
					if (verbose==TRUE) {
					  if (pvals3<=alpha) {
  						ang <- paste0("Fligner-Killeen test [fligner.test()] -\n\t=> Significant differences of variance between groups, p-value:", .format_pval(pvals3),"\n\t... go to med1way().")
  						fr <- paste0("Test de Fligner-Killeen [fligner.test()] -\n\t=> Différences significatives de variance entre les groupes, p-value : ", .format_pval(pvals3),"\n\t... aller vers med1way().")
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  } else {
  						ang <- paste0("Fligner-Killeen test [fligner.test()] -\n\t=> Non-significant differences of variance between groups, p-value:",.format_pval(pvals3),"\n\t\t...go to t1way().")
  						fr <- paste0("Test de Fligner-Killeen [fligner.test()] -\n\t=> Différences non significatives de variance entre les groupes, p-value : ",.format_pval(pvals3),"\n\t\t...aller vers t1way().")
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  }
					  temp <- pairwise(x,g,type="ks",silent=silent,boot=FALSE)$p.value
					  ks_result <- min(unlist(temp),na.rm=TRUE)
					  if (ks_result <= pval) {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\tWarning! the groups do not have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t... by comparing with Sidak corrected alpha ", .format_pval(pval), "\n\tThe Kruskal-Wallis test and Mann-Whitney-Wilcoxon test will be less reliable.\n\tPlease, check graphically the groups distributions.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\tAttention ! les groupes n'ont pas la même distribution. min(p-value) : ", .format_pval(ks_result), "\n\t... en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), "\n\tLe test de Kruskal-Wallis sera moins fiable.\n\tVeuillez vérifier graphiquement les distributions des groupes.")
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  } else {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\tThe groups have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t... by comparing with Sidak corrected alpha ", .format_pval(pval), "\n\tGood accuracy expected on the tests of Kruskal-Wallis and Mann-Whitney-Wilcoxon.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\tLes groupes ont la même distribution. min(p-value) : ", .format_pval(ks_result), "\n\t... en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), "\n\tBonne précision attendue pour le test de Kruskal-Wallis.")
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  }
  					  if (pvals <= alpha) {
  						ang <- paste0("Kruskal-Wallis test [kruskal.test()] - At least one group appears to show a difference. p-value:", .format_pval(pvals))
  						fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] - Au moins un groupe semble montrer une différence. p-value : ", .format_pval(pvals))
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  } else {
  						ang <- paste0("Kruskal-Wallis test [kruskal.test()] - No different group a priori. p-value: ", .format_pval(pvals))
  						fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] - Aucun groupe différent a priori. p-value : ", .format_pval(pvals))
  						k <- .vbse(ang, fr, verbose = verbose, k = k)
					  }
					}
					if ((return==TRUE) | (verbose==TRUE)) {
					  ############
					  # 	Allons jusqu'au bout du raisonnement
					  ############
					  if (pvals3<=alpha) {
						#################
						#	Grande hétéroscédasticité - var.equal=F ==> median med1way
						#################
						check_variance_equal <-  FALSE
						pvals3 <- med1way(x~g)$p.value
						if (is.na(pvals3)) {
						    ang <- paste0("Oneway ANOVA of medians [med1way()] - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.")
							fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] - Échec, retourne NA. Le test de Kruskal-Wallis doit être utilisé pour l'interprétation.")
							k <- .vbse(ang, fr, verbose = verbose, k = k)
						} else {
						  if (pvals3 <= alpha) {
							ang <- paste0("Oneway ANOVA of medians [med1way()] :\n\tSignificant differences between the median of groups. p-value:", .format_pval(pvals3))
							fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] :\n\tDifférences significatives entre les médianes des groupes. p-value :", .format_pval(pvals3))
							k <- .vbse(ang, fr, verbose = verbose, k = k)
							if (pvals > alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on medians give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance des médianes donnent des résultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
							}
						  } else if (pvals3 > alpha) {
							ang <- paste0("Oneway ANOVA of medians [med1way()] -\n\t=> Non-significant differences between the median of groups. p-value:", .format_pval(pvals3))
							fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] -\n\t=> Différences non significatives entre les médianes des groupes. p-value :", .format_pval(pvals3))
							k <- .vbse(ang, fr, verbose = verbose, k = k)
							if (pvals <= alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on medians give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance des médianes donnent des résultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
							}
						  }
						}
					  } else {
						#################
						#	Homoscédasticité tolérable var.equal=T ==> trimm t1way
						#################
						check_variance_equal <-  TRUE
						pvals3 <- t1way(x~g)$p.value
						if (is.na(pvals3)) {
							ang <- paste0("Oneway ANOVA on trimmed means [t1way()] - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.")
							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] - Échec, retourne NA. Le test de Kruskal-Wallis doit être utilisé pour l'interprétation.")
							k <- .vbse(ang, fr, verbose = verbose, k = k)
						} else {
						  if (pvals3 <= alpha) {
  							ang <- paste0("Oneway ANOVA on trimmed means [t1way()] :\n\t==> Significant differences between the trimmed groups. p-value:", .format_pval(pvals3))
  							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] :\n\t==> Différences significatives entre les groupes tronqués. p-value :", .format_pval(pvals3))
  							k <- .vbse(ang, fr, verbose = verbose, k = k)
							if (pvals > alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronquées donnent des résultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
							}
						  } else if (pvals3 > alpha) {
  							ang <- paste0("One-way ANOVA on trimmed means [t1way()] - Non-significant differences between the trimmed groups. p-value:", .format_pval(pvals3))
  							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] - Différences non significatives entre les groupes tronqués. p-value :", .format_pval(pvals3))
  							k <- .vbse(ang, fr, verbose = verbose, k = k)
  							if (pvals <= alpha) {
  								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
  								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronquées donnent des résultats contradictoires.")
  								k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
  							}
						  }
						}
					  } # fin var.equal
					} # test uniquement pour aller au bout du raisonnement
				} # fin d'ajustement à la normalité malgré tout
	  } else { # Si NORMAL
		################################################################
		#
		#
		#		NORMAL
		#
		#
		################################################################
		###################################################
		#			NORMAL
		###################################################
		if (check_number==2) { # 2 categories
		  if (code==TRUE){cat("length(unique(g)) #2)\n")}
		  #ang <- paste0("Two categories.")
			#fr <- paste0("Deux catégories.")
			#k <- .vbse(ang, fr, verbose = verbose, k = k)
		  fm <- formula(x~g)
		  pvals <- var.test(fm)$p.value
		  if (code==TRUE){cat("var.test(x~g) #3)\n")}
		  if (pvals>alpha) {
			###################################################
			#			NORMAL		2 categories	homogene variance
			###################################################
			pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=TRUE)$p.value
			if (verbose==TRUE) {
				if (pvals <= alpha) {
				  ##################################
				  #        NORMAL, Student significant (var.equal=T)
				  ##################################
				  ang <- paste0("Student test [t.test()] :\n\t==> Significant differences between groups. p-value:", .format_pval(pvals))
				  fr <- paste0("Test de Student [t.test()] :\n\t==> Différences significatives entre les groupes. p-value : ", .format_pval(pvals))
				  k <- .vbse(ang, fr, verbose = verbose, k = k)
				} else {
				  ##################################
				  #        NORMAL, Student unsignificant (var.equal=T)
				  ##################################
				  ang <- paste0("Student test [t.test()] - Non-significant differences between groups. p-value:", .format_pval(pvals))
				  fr <- paste0("Test de Student [t.test()] - Différences non significatives entre les groupes. p-value : ", .format_pval(pvals))
				  k <- .vbse(ang, fr, verbose = verbose, k = k)
				}
			}
		  } else {
			###################################################
			#			NORMAL		2 categories	non-homogene variance
			###################################################
			check_variance_equal <-  FALSE
			ang <- paste0("Fisher-Snedecor test [var.test()] - Non-identical group variances. p-value:", .format_pval(pvals))
			fr <- paste0("Test de Fisher-Snedecor [var.test()] - Variances des groupes non identiques. p-value : ", .format_pval(pvals))
			k <- .vbse(ang, fr, verbose = verbose, k = k)
			pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=FALSE)$p.value
			if (code==TRUE){cat("t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=TRUE) #4)\n")}
			if (pvals <= alpha) {
			  ##################################
			  #        NORMAL, Student significant (var.equal=F)
			  ##################################
				ang <- paste0("Student test [t.test()] :\n\t==> Significant differences between groups. p-value:", .format_pval(pvals))
				fr <- paste0("Test de Student [t.test()] :\n\t==> Différences significatives entre les groupes. p-value : ", .format_pval(pvals))
				k <- .vbse(ang, fr, verbose = verbose, k = k)

			} else {
			  ##################################
			  #        NORMAL, Student unsignificant (var.equal=F) 	Acceptable for t.test() de Welch
			  ##################################
				ang <- paste0("Student test (t.test()) - Non-significant differences between groups. p-value:", .format_pval(pvals))
				fr <- paste0("Test de Student (t.test()) - Différences non significatives entre les groupes. p-value :", .format_pval(pvals))
				k <- .vbse(ang, fr, verbose = verbose, k = k)

				if (code==TRUE){cat("t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=FALSE)   #4)\n")}
				pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=FALSE)$p.value
				if (pvals <= alpha) {
				  k <- .vbse("Student Test [t.test()]- Significant differences between groups.",
				             "Test de Student [t.test()] - Différences significatives entre les groupes.",
				             verbose = verbose, k = k)
				} else {
				  k <- .vbse("Student test [t.test()] -\n\t=> Non-significant differences between groups.",
				             "Test de Student [t.test()] -\n\t=> Différences non significatives entre les groupes.",
				             verbose = verbose, k = k)
				}


			}
		  }
		} else { 													# > 2 categories
		  ###################################################
		  #			NORMAL		>2 categories
		  ###################################################
		  if (check_variance_equal==TRUE) {											# Identical variances
			  ##################################
			  #        NORMAL, >2 categories var.equal=T => AOV
			  ##################################
			fm <- formula(x~g)
			mya <- suppressWarnings(aov(data.frame(x,g), formula=fm))
			if (code==TRUE){cat("mya <- aov(data.frame(x,g), formula=x~g) #4)\n")}
			pvals <- summary(mya)[[1]][["Pr(>F)"]][1]
				if (pvals<=alpha) {										# Significant AOV
				  k <- .vbse(paste0("One-way analysis of variance [aov()] :\n\t==> Significant differences between groups. p-value: ", .format_pval(pvals)),
					  paste0("Analyse de variance à un facteur [aov()] :\n\t==> Différences significatives entre les groupes. p-value : ", .format_pval(pvals)),
					  verbose = verbose, k = k)
				} else {											#	 Non-significant AOV
				  k <- .vbse(paste0("One-way analysis of variance [aov()] :\n\t==> Non-significant differences between groups. p-value: ", .format_pval(pvals)),
					paste0("Analyse de variance à un facteur [aov()] :\n\t==> Différences non significatives entre les groupes. p-value : ", .format_pval(pvals)),
					verbose = verbose, k = k)
				}
		  } else {												# Non-identical variances
			  ##################################
			  #        NORMAL, >2 categories var.equal=F => FANOVA.HETERO
			  ##################################
			pvals <- oneway.test(x~g,var.equal=FALSE)$p.value
			if (code==TRUE){cat("oneway.test(x~g,var.equal=FALSE) #4)\n")}
			myf <- try(fanova.hetero(data.frame(x,g = as.factor(g)),x~g),silent=silent)
			if (inherits(myf,"try-error")) {
			  #if (verbose==TRUE) {cat("Error on fanova.hetero()\n")}
			  pvals2 <- alpha
			} else {pvals2 <- myf$ans[4]}
			if (verbose==TRUE) {
				if (pvals<=alpha) {
					ang <- paste0("Welch’s heteroscedastic F test [oneway.test(var.equal=FALSE)] :\n\t==> Significant differences between groups. p-value:", .format_pval(pvals))
					fr <- paste0("Test F hétéroscédastique de Welch [oneway.test(var.equal=FALSE)] :\n\t==> Différences significatives entre les groupes. p-value :", .format_pval(pvals))
					k <- .vbse(ang, fr, verbose = verbose, k = k)
					if (pvals2 > alpha) {
						ang <- paste0("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:", .format_pval(pvals2))
						fr <- paste0("Attention ! fanova.hetero() ne donne pas le même résultat que oneway.test. p-value :", .format_pval(pvals2))
						k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
					}
				} else {
					ang <- paste0("Welch’s heteroscedastic F test [oneway.test(var.equal=FALSE)] - Non-significant differences between groups. p-value:", .format_pval(pvals))
					fr <- paste0("Test F hétéroscédastique de Welch [oneway.test(var.equal=FALSE)] - Différences non significatives entre les groupes. p-value :", .format_pval(pvals))
					k <- .vbse(ang, fr, verbose = verbose, k = k)
					if (pvals2 <= alpha) {
						ang <- paste0("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:", .format_pval(pvals))
						fr <- paste0("Attention ! fanova.hetero() ne donne pas le même résultat que oneway.test. p-value :", .format_pval(pvals))
						k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
					}
				}
			}


		  }
		}
	  }
	  return(check <- list(x,g,check_normality,check_variance_equal,k))
 }
