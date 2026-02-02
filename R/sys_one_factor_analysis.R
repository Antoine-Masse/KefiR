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
# Fonction auxiliaire pour afficher code R avec numérotation (ANOVA simple)
.code_anova <- function(step_num, title, code_lines) {
  cat(paste0("# ", step_num, ") ", title, "\n"))
  for (line in code_lines) {
    cat(paste0(line, "\n"))
  }
  cat("\n")
}

#' @export
.one_factor_analysis <- function(x=NULL,g=NULL,formula=NULL,data=NULL,
	paired = FALSE, id = NULL,alpha = 0.05,return=TRUE,
	k=NULL,code=FALSE,debug=FALSE,verbose=FALSE, boot = TRUE, silent=TRUE) {

  # Protection contre code=NULL (si utilisateur a utilisé code=F avec variable F dans données)
  if (is.null(code)) {
    warning("Le paramètre 'code' était NULL. Utilisez code=FALSE ou code=TRUE (pas F ou T). Réinitialisation à FALSE.")
    code <- FALSE
  }

  # Si code==TRUE, désactiver verbose MAIS garder k synchronisé
  if (isTRUE(code)) {
    verbose_original <- verbose
    verbose <- FALSE
    # NOTE: On utilisera k au lieu de k pour synchroniser la numérotation
  }

  ################################################################
  #
  #
  #		Si paired = TRUE et 2 groupes (k  = 2)
  #
  #
  ################################################################
  ## --- Vérifs d’intégrité pour le mode apparié ------------------------------
  if (isTRUE(paired)) {
    # Si identifiants fourni, le contrôler et réaligner les données dessus
    if(!is.null(id)) {
      g2 <- droplevels(factor(g))
      dfc <- data.frame(id = data[id], g = g2, x = as.numeric(x))
      tab <- table(dfc$id, dfc$g)              # matrice (#id x 2)
      # 1) chaque id doit être présent dans les 2 modalités
      miss <- rownames(tab)[rowSums(tab > 0) < 2L]
      if (length(miss)) {
        .exit(
          "Some ids are incomplete (no match for each modality of 'g').",
          "Certains identifiants de id sont incomplets (pas de correspondance pour chaque modalité de 'g').",
          verbose = verbose, code = code, return = return
        )
      }
      # 2) (option strict) exactement une obs par id et par modalité
      if (any(tab != 1L)) {
        .exit(
          "Each id must have exactly one observation per level of 'g' (no duplicates).",
          "Chaque identifiant doit avoir exactement une observation par modalité de 'g' (pas de doublons).",
          verbose = verbose, code = code, return = return
        )
      }
      # Réaligner
      k <- .vbse("Alignment by id for paired (2 levels).",
            "Alignement par id pour les données appariées.",
            verbose = verbose, code = code, k = k, cpt = "on") -> k
      temp <- .align_pairs(x,data[id],g)
      x <- temp$x
      g <- temp$g
    }
  #
    #		Si paired = TRUE et identifiants fournis	-->
    if (nlevels(g)==2) {
      diff <- x[g==levels(g)[1]]-x[g==levels(g)[2]]
      check_normality <- .normality(diff, paired=paired, k=k, verbose=verbose, alpha=alpha)
      k <- check_normality[[2]] ; check_normality <- check_normality[[1]]
      check_variance_equal <- TRUE

      # =============================================================================
      # TEST T APPARIÉ (ou Wilcoxon signé) - Ajouté pour afficher le résultat du test
      # =============================================================================
      lv <- levels(droplevels(factor(g)))
      if (isTRUE(check_normality)) {
        # Différences normales => t-test apparié
        pvals <- t.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
        if (isTRUE(code)) {
          .code_anova(8, "Test de Student apparié", c(
            "lv <- levels(droplevels(factor(g)))",
            "t.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
          ))
        }
        if (isTRUE(verbose)) {
          if (pvals <= alpha) {
            ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Student apparié [t.test(paired=TRUE)] :\n\t==> Différences significatives entre les conditions. p-value : ", .format_pval(pvals))
          } else {
            ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Student apparié [t.test(paired=TRUE)] :\n\t==> Différences non significatives entre les conditions. p-value : ", .format_pval(pvals))
          }
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
        }
      } else {
        # Différences non normales => Wilcoxon signé
        pvals <- wilcox.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
        if (isTRUE(code)) {
          .code_anova(8, "Test de Wilcoxon signé (données appariées non normales)", c(
            "lv <- levels(droplevels(factor(g)))",
            "wilcox.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
          ))
        }
        if (isTRUE(verbose)) {
          if (pvals <= alpha) {
            ang <- paste0("Wilcoxon signed-rank test [wilcox.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Wilcoxon signé [wilcox.test(paired=TRUE)] :\n\t==> Différences significatives entre les conditions. p-value : ", .format_pval(pvals))
          } else {
            ang <- paste0("Wilcoxon signed-rank test [wilcox.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Wilcoxon signé [wilcox.test(paired=TRUE)] :\n\t==> Différences non significatives entre les conditions. p-value : ", .format_pval(pvals))
          }
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
        }
      }

      return(check <- list(x = x, g = g, check_normality = check_normality,
                           check_variance_equal = check_variance_equal, k = k,
                           global_pvalue = pvals, chosen_test = NULL))
    }
  } # Fin du scénario paired sachant que l'alignement si n_g > 2 (k) passe en .multi_factor_analysis	-->

  #==================================================
  #         Initiation des variables de contrôle
  #==================================================
  check_discret <- discret.test(x)
  check_normality <- c()
  check_number <- c()
  check_variance_equal <-  c()
  pvals <- NA  # Initialize pvals to avoid "object not found" errors at return
  chosen_test <- NULL  # Initialize: NULL=parametric, "med1way", "t1way", or "kruskal" for non-parametric

  #==================================================
  #         Contrôle du nombre de groupes
  #==================================================

  check_number <- length(unique(g))
  if (check_number==2) { 							# 2 categories
    if (paired==TRUE) {
      k <- .vbse("Two conditions.",
                 "Deux conditions.",
                 verbose = verbose, code = code, k = k, cpt="on")
      if (isTRUE(code)) {
        .code_anova(1, "Deux conditions (données appariées)", c(
          "length(unique(g))  # Doit être 2"
        ))
      }
    } else{
      k <- .vbse("Two groups.",
                 "Deux groupes.",
                 verbose = verbose, code = code, k = k, cpt="on")
      if (isTRUE(code)) {
        .code_anova(1, "Deux groupes indépendants", c(
          "length(unique(g))  # Doit être 2"
        ))
      }
    }
  } else if (check_number>2) {
    if (paired==TRUE) {
      k <- .vbse("More than two conditions.",
        "Plus de deux conditions.",
        verbose = verbose, code = code, k = k)
    } else{
      k <- .vbse("More than two groups.",
                 "Plus de deux groupes.",
                 verbose = verbose, code = code, k = k)
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
    # Données non appariées - ASSOMPTION 1/3 pour ANOVA 1 facteur
    k <- .vbse(paste0("ASSUMPTION 1/3: Independence check - Ensure that observations between groups are independent.\n",
                      "\tVerify that:\n",
                      "\t  • No repeated measures\n",
                      "\t  • No clustering effects\n",
                      "\t  • No carry-over effects"),
               paste0("ASSOMPTION 1/3 : Contrôle de l'indépendance - vérifiez que les observations entre les groupes sont indépendantes.\n",
                      "\tVérifiez que :\n",
                      "\t  • Pas de mesures répétées\n",
                      "\t  • Pas d'effet cluster\n",
                      "\t  • Pas d'effets report [influence sur l'ordre des mesures]"),
               verbose = verbose, code = code, k = k, cpt = "on")
    if (isTRUE(code)) {
      .code_anova(2, "Contrôle indépendance (groupes indépendants)", c())
    }
  } else {
    # Données appariées - ASSOMPTION 1/3 pour ANOVA 1 facteur apparié
    k <- .vbse("ASSUMPTION 1/3: Paired independence check - Ensure that the dependence structure is only due to the intended pairing\n\t* e.g., same subject before/after\n\t* no additional clustering between paired observations\n\t* no carry-over effects between paired observations",
               "ASSOMPTION 1/3 : Contrôle de l'indépendance pour données appariées - Vérifiez que la structure de dépendance provienne uniquement de l'appariement prévu\n\t* ex. : même sujet avant/après\n\t* pas de clustering supplémentaire entre observations appariées\n\t* pas d'effets report supplémentaires entre observations appariées",
               verbose = verbose, code = code, k = k, cpt = "on")
    if (isTRUE(code)) {
      .code_anova(2, "Contrôle indépendance (données appariées)", c())
    }
  }
  #==================================================
  # Détection des Outlier en X
  #==================================================
  # Helper: safe outlier detection with fallback
  .safe_identify_outliers <- function(data) {
    if (requireNamespace("rstatix", quietly = TRUE)) {
      return(rstatix::identify_outliers(data))
    } else {
      # Fallback: IQR method
      x <- data[[1]]
      Q1 <- quantile(x, 0.25, na.rm = TRUE)
      Q3 <- quantile(x, 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      lower_extreme <- Q1 - 3 * IQR
      upper_extreme <- Q3 + 3 * IQR

      return(data.frame(
        data,
        is.outlier = (x < lower_bound | x > upper_bound),
        is.extreme = (x < lower_extreme | x > upper_extreme)
      ))
    }
  }

  if (paired == FALSE) {
    #---------------------------------
    #  Gestion des données extrêmes si non appariés
    #---------------------------------
    outliers_base <- as.vector(unlist(
      by(x, g, function(z) {
        base::length(which(.safe_identify_outliers(data.frame(z))$is.outlier == TRUE)) / base::length(z)
      })
    ))
    outliers_extrem <- as.vector(unlist(
      by(x, g, function(z) {
        base::length(which(.safe_identify_outliers(data.frame(z))$is.extreme == TRUE)) / base::length(z)
      })
    ))
    max_outliers_base  <- max(outliers_base)
    max_outliers_extrem <- max(outliers_extrem)
    # Messages gradués selon la sévérité
    if (max_outliers_extrem > 0) {
      k <- .vbse(
        paste0("Outlier detection [identify_outliers() {rstatix}]\n\t",
               "Result: ", round(max_outliers_extrem * 100, 1), "% EXTREME outliers (maximum in a group)\n\t",
               "==> Extreme values may strongly influence parametric tests"),
        paste0("Détection outliers [identify_outliers() {rstatix}]\n\t",
               "Résultat : ", round(max_outliers_extrem * 100, 1), "% outliers EXTRÊMES (maximum dans un groupe)\n\t",
               "==> Valeurs extrêmes peuvent fortement influencer tests paramétriques."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      if (isTRUE(code)) {
        .code_anova(3, "Détection outliers extrêmes (groupes indépendants)", c(
          "library(rstatix)",
          "outliers_by_group <- by(x, g, function(z) {",
          "  result <- identify_outliers(data.frame(z))",
          "  prop_extreme <- sum(result$is.extreme, na.rm = TRUE) / length(z)",
          "  return(prop_extreme)",
          "})",
          "max_outliers_extrem <- max(unlist(outliers_by_group))"
        ))
      }
    } else if (max_outliers_base > 0) {
      k <- .vbse(
        paste0("Outlier detection [identify_outliers() {rstatix}]\n\t",
               "Result: ", round(max_outliers_base * 100, 1), "% outliers (maximum in a group, not extreme)\n\t",
               "==> Values monitored but not excessive"),
        paste0("Détection outliers [identify_outliers() {rstatix}]\n\t",
               "Résultat : ", round(max_outliers_base * 100, 1), "% outliers (maximum dans un groupe, non extrêmes)\n\t",
               "==> Valeurs à surveiller, mais non excessives"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      if (isTRUE(code)) {
        .code_anova(3, "Détection outliers modérés (groupes indépendants)", c(
          "library(rstatix)",
          "outliers_by_group <- by(x, g, function(z) {",
          "  result <- identify_outliers(data.frame(z))",
          "  prop_outliers <- sum(result$is.outlier, na.rm = TRUE) / length(z)",
          "  return(prop_outliers)",
          "})",
          "max_outliers_base <- max(unlist(outliers_by_group))"
        ))
      }
    }
  } else {
    #---------------------------------
    # Données appariées — contrôle des outliers sur les différences intra-sujet
    lv <- levels(droplevels(factor(g)))
    #---------------------------------
    if (length(lv) == 2L) {
      dif <- x[g == lv[2]] - x[g == lv[1]]
      io  <- suppressWarnings(.safe_identify_outliers(data.frame(dif)))

      prop_ext <- if ("is.extreme" %in% names(io)) mean(io$is.extreme, na.rm = TRUE) else 0
      prop_out <- if ("is.outlier" %in% names(io)) mean(io$is.outlier, na.rm = TRUE) else 0

      if (prop_ext > 0) {
        k <- .vbse(
          paste0("Extreme outliers in paired differences [identify_outliers() from {rstatix}].\n\t",
                 "Proportion of extreme outliers: ", round(prop_ext * 100, 1), "%.\n\t",
                 "These pairs may strongly influence paired parametric tests."),
          paste0("Outliers extrêmes dans les différences appariées [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d'outliers extrêmes : ", round(prop_ext * 100, 1), "%.\n\t",
                 "Ces paires peuvent fortement influencer les tests paramétriques appariés."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        if (isTRUE(code)) {
          .code_anova(3, "Détection outliers extrêmes (différences appariées)", c(
            "lv <- levels(droplevels(factor(g)))",
            "dif <- x[g == lv[2]] - x[g == lv[1]]",
            "library(rstatix)",
            "io <- identify_outliers(data.frame(dif))",
            "prop_ext <- mean(io$is.extreme, na.rm = TRUE)"
          ))
        }
      } else if (prop_out > 0) {
        k <- .vbse(
          paste0("Outliers in paired differences [identify_outliers() from {rstatix}].\n\t",
                 "Proportion of outliers: ", round(prop_out * 100, 1), "%.\n\t",
                 "Monitored but not extreme."),
          paste0("Outliers dans les différences appariées [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d'outliers : ", round(prop_out * 100, 1), "%.\n\t",
                 "À surveiller mais non extrêmes."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        if (isTRUE(code)) {
          .code_anova(3, "Détection outliers modérés (différences appariées)", c(
            "lv <- levels(droplevels(factor(g)))",
            "dif <- x[g == lv[2]] - x[g == lv[1]]",
            "library(rstatix)",
            "io <- identify_outliers(data.frame(dif))",
            "prop_out <- mean(io$is.outlier, na.rm = TRUE)"
          ))
        }
      }
    } else {
      k <- .vbse(
        "Paired data with more than two levels: outlier check on differences is not defined here.",
        "Données appariées avec plus de deux modalités : contrôle des outliers sur les différences non défini ici.",
        verbose = verbose, code = code, k = k, cpt = "off"
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
    # Sidak devient trop conservatif au-delà de 10 groupes	-->
    # Bender, R., & Lange, S. (2001). Adjusting for multiple testing—when and how?  DOI: 10.1016/S0895-4356(00)00314-0
    # Il est justifié de passer ensuite au bootstrap
    # Méthode FDR la correction par FDR (False Discovery Rate) – notamment via la méthode Benjamini-Hochberg (BH)
    # Murray, E. J., Berrie, L., & Matthews, J. N. S. (2021). Understanding multiplicity in clinical trials: the impact of multiple endpoints in the interpretation of treatment effects.
    # Attention : FDR est dans une démarche exploratoire ou descriptive, pas une validation réglementaire.
    .dbg("Sidak's correction rather than Bonferroni for independent tests (conservative).",
         "Correction de Sidak plutôt que Bonferroni pour les tests indépendants (conservateur).",debug=debug)
    pval <- 1-(1-alpha)^(1/length(unique(g)))

    if (check_number==2) { 							# 2 categories
      # Contrôle de la normalité des groupes (ou des différences si apparié)
      check_normality <- .normality(x, g, alpha=alpha, k=k, verbose=verbose, paired=paired)
    } else {
      # Annonce de l'ajustement du modèle ANOVA
      k <- .vbse(
        "Fitting ANOVA model [aov()].",
        "Ajustement du modèle ANOVA [aov()].",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      if (isTRUE(code)) {
        .code_anova(4, "Ajustement du modèle ANOVA", c(
          "myaov <- aov(x ~ g)"
        ))
      }

      # Contrôle de la normalité des résidus
      myaov <- aov(x~g)
      n_residu <- length(x)

      # ASSOMPTION 2/3 pour ANOVA 1 facteur (1=indépendance, 2=normalité, 3=variance)
      k <- .vbse(
        "ASSUMPTION 2/3: Academic check of ANOVA model residuals normality.",
        "ASSOMPTION 2/3 : Contrôle ACADÉMIQUE de la normalité des résidus du modèle ANOVA.",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      check_normality <- .normality(myaov$residuals, alpha=alpha, k=k, verbose=verbose, cpt="off") #pvals <- .normality(x,g)
    }
    k <- check_normality[[2]] ; check_normality <- check_normality[[1]]

    # Code generation for normality check (after .normality() call)
    if (isTRUE(code)) {
      if (check_number == 2) {
        # Detect group sizes and heterogeneity (same logic as .normality())
        group_sizes <- table(g)
        min_group_size <- min(group_sizes)
        max_group_size <- max(group_sizes)

        # Detect heterogeneity (reajust flag in .normality())
        reajust <- (min_group_size <= 50 && max_group_size > 50) ||
                   (min_group_size <= 500 && max_group_size > 500)

        # Generate code following .normality() logic
        if (!reajust && min_group_size <= 50) {
          # Scenario 1: n ≤ 50 -> Shapiro-Wilk only
          code_lines <- c(
            paste0("alpha_sidak <- 1 - (1 - ", alpha, ")^(1/length(unique(g)))"),
            "by(x, g, shapiro.test)"
          )
        } else if (!reajust && min_group_size <= 500) {
          # Scenario 2: 50 < n ≤ 500 -> Jarque-Bera with Sidak correction
          code_lines <- c(
            "library(KefiR)  # Pour jb.norm.test()",
            paste0("alpha_sidak <- 1 - (1 - ", alpha, ")^(1/length(unique(g)))"),
            "by(x, g, function(z) jb.norm.test(z)$p.value)"
          )
        } else {
          # Scenario 3: n > 500 OR heterogeneity -> Skewness/Kurtosis
          code_lines <- c(
            "library(agricolae)",
            "by(x, g, function(z) {",
            "  sk <- abs(skewness(z))",
            "  ku <- abs(kurtosis(z))",
            "  list(skewness = sk, kurtosis = ku, normal = sk <= 1 && ku <= 1.5)",
            "})"
          )
        }
        .code_anova(5, "ASSOMPTION 2 : Contrôle de la normalité (2 groupes)", code_lines)
      } else {
        # Detect residuals size (same logic as .normality() for single group)
        n_residuals <- length(x)

        # Generate code following .normality() logic for residuals
        if (n_residuals <= 50) {
          # Scenario 1: n ≤ 50 -> Shapiro-Wilk only
          code_lines <- c(
            "residus <- residuals(myaov)",
            "shapiro.test(residus)"
          )
        } else if (n_residuals <= 500) {
          # Scenario 2: 50 < n ≤ 500 -> Jarque-Bera
          code_lines <- c(
            "library(KefiR)  # Pour jb.norm.test()",
            "residus <- residuals(myaov)",
            "jb.norm.test(residus)"
          )
        } else {
          # Scenario 3: n > 500 -> Skewness/Kurtosis
          code_lines <- c(
            "library(agricolae)",
            "residus <- residuals(myaov)",
            "sk <- abs(skewness(residus))",
            "ku <- abs(kurtosis(residus))",
            "# Normal si |skewness| <= 1 et |kurtosis| <= 1.5",
            "list(skewness = sk, kurtosis = ku, normal = sk <= 1 && ku <= 1.5)"
          )
        }
        .code_anova(5, "ASSOMPTION 2 : Contrôle de la normalité (résidus ANOVA)", code_lines)
      }
    }
    #==================================================
    #
    #   Contrôle de la variance
    #
    #==================================================
    .dbg("Variance check with .variance().",
         "Contrôle de la variance avec .variance().",debug=debug)
    # Pour ANOVA 1 facteur : ASSOMPTION 3/3 (1=indépendance, 2=normalité, 3=variance)
    variance_temp <- .variance(x, g, check_normality=check_normality, alpha=alpha, paired=paired,
                  debug = debug, verbose=verbose, k=k, code=FALSE, assumption_label="3/3")
    check_variance_equal <- variance_temp[[1]]
    k <- variance_temp[[2]]

    # Générer le code de variance ICI avec numérotation k
    if (isTRUE(code)) {
      if (check_number == 2) {
        .code_anova(6, "Test de variance (Fisher-Snedecor)", c(
          "var.test(x~g)"
        ))
      } else {
        if (check_normality == TRUE) {
          .code_anova(6, "Test d'homogénéité des variances (Bartlett)", c(
            "bartlett.test(x, g)"
          ))
        } else {
          .code_anova(6, "Test d'homogénéité des variances (Levene)", c(
            "library(car)",
            "leveneTest(x, g)"
          ))
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
  	#		Auto-analyse pour envisager un retour vers la situation paramétrique (Kurtosis & Skweness	-->)
  	#
  	##########################
	  # |Skewness| < 1 et |Kurtosis| < 1.5 sont les seuils recommandés pour considérer une distribution comme "approximativement normale", selon Kline (2011).
	  # Kline, R. B. (2011). Principles and practice of structural equation modeling (4th ed.). New York, NY: Guilford Press.
	  # Normalité variables indicates acceptable normality for most measures, based on Kline's (2023) guidelines, where skewness values should ideally be within ±3 and kurtosis values within ±10.
	  # Permissivité extrême : skewness > 2 et kurtosis > 7 Blanca et al. (2018)
	  # Blanca, M. J., Alarcón, R., Arnau, J., Bono, R., & Bendayan, R. (2017). Non-normal data: Is ANOVA still a valid option? Psicothema, 29(4), 552-557. https://diposit.ub.edu/dspace/bitstream/2445/122126/1/671797.pdf#:~\:text=sought%20to%20answer%20the%20following,In%20addition
    # Une violation de la normalité a peu d'effet sur l'ANOVA si variances homogènes.Caldwell (2022)
    # Caldwell, A. (2022, March 31). Chapter 12 Violations of Assumptions. In Power Analysis with Superpower. https://aaroncaldwell.us/SuperpowerBook/violations-of-assumptions.html

    auto_ku_sk <- function(x, g, check_normality = TRUE, alpha = 0.05, verbose = TRUE, k = NULL, code = NULL) {
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

      # Normalité résidus (tolérante) - ÉTAPE 7
      if (ng == 2L) {
        tmp <- .normality(x, g, alpha = alpha, tolerance = "extrem", k = k, verbose = verbose, code = code, cpt = "off")
        check_normality <- tmp[[1]]; k <- tmp[[2]]
      } else {
        fm <- x ~ g
        myaov <- suppressWarnings(aov(fm, data = data.frame(x = x, g = g)))
        tmp <- .normality(myaov$residuals, alpha = alpha, tolerance = "extrem", k = k, verbose = verbose, code = code, cpt = "off")
        check_normality <- tmp[[1]]; k <- tmp[[2]]
      }

      # Retour
      list(check_normality, k, if(ng > 2L && exists("myaov")) myaov else NULL)
    }

    # Fonction diagnostic ANOVA (ÉTAPE 8) - optionnelle pour garantir fiabilité
    diagnostic_anova <- function(x, g, myaov = NULL, alpha = 0.05, verbose = TRUE, k = NULL, code = NULL) {
      # Coercitions
      if (is.data.frame(x)) x <- x[[1]]
      if (is.matrix(x))     x <- x[,1,drop=TRUE]
      if (is.list(x))       x <- unlist(x, use.names = FALSE)
      x <- as.vector(x)

      if (is.data.frame(g)) g <- g[[1]]
      if (is.matrix(g))     g <- g[,1,drop=TRUE]
      if (is.list(g))       g <- unlist(g, use.names = FALSE)
      g <- droplevels(factor(g))
      ng  <- nlevels(g)

      # En-tête DIAGNOSTIC
      k <- .vbse(
        "Optional DIAGNOSTIC to ensure ANOVA reading reliability:",
        "DIAGNOSTIC optionnel pour garantir la fiabilité de la lecture de l'ANOVA :",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      # a) Outliers sur RÉSIDUS par groupe
      if (!is.null(myaov)) {
        residuals_data <- myaov$residuals

        # Détecter outliers sur résidus
        tab_res <- .safe_identify_outliers(data.frame(res = residuals_data))
        n_total <- length(residuals_data)
        n_outliers <- sum(tab_res$is.outlier %in% TRUE, na.rm = TRUE)
        n_extreme <- sum(tab_res$is.extreme %in% TRUE, na.rm = TRUE)

        if (n_outliers > 0 || n_extreme > 0) {
          prop_out <- round(100 * n_outliers / n_total, 1)
          prop_ext <- round(100 * n_extreme / n_total, 1)
          k <- .vbse(
            paste0("a) Outliers detected on residuals:\n",
                   "\t\t==> ", n_outliers, " outliers (", prop_out, "%), ",
                   n_extreme, " extreme (", prop_ext, "%)."),
            paste0("a) Outliers détectés sur les résidus :\n",
                   "\t\t==> ", n_outliers, " outliers (", prop_out, "%), ",
                   n_extreme, " extrêmes (", prop_ext, "%)."),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            "a) No outliers detected on residuals.",
            "a) Aucun outlier détecté sur les résidus.",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      }

      # b) Distance de Cook (diagnostic d'influence)
      if (!is.null(myaov) && ng > 2L) {
        cook_dist <- cooks.distance(myaov)
        n_total <- length(x)
        cook_threshold <- 4 / n_total
        influential <- which(cook_dist > cook_threshold)
        n_influential <- length(influential)

        if (n_influential > 0) {
          prop_influential <- round(100 * n_influential / n_total, 1)
          k <- .vbse(
            paste0("b) Diagnostic d'influence (distance de Cook) :\n",
                   "\t\t==> ", n_influential, " observations influentes détectées.\n",
                   "\t\t(", prop_influential, "% dépassent seuil 4/n = ", round(cook_threshold, 4), ")."),
            paste0("b) Diagnostic d'influence (distance de Cook) :\n",
                   "\t\t==> ", n_influential, " observations influentes détectées.\n",
                   "\t\t(", prop_influential, "% dépassent seuil 4/n = ", round(cook_threshold, 4), ")."),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            "b) Cook's distance: No influential observations detected (all < 4/n).",
            "b) Distance de Cook : Aucune observation influente détectée (toutes < 4/n).",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      }

      # c) Normalité diagnostique par groupe
      # Appeler .normality() silencieusement pour obtenir le résultat
      tmp2 <- .normality(x, g, alpha = alpha, tolerance = "extrem", k = k, verbose = FALSE, code = FALSE, cpt = "off")
      check_normality_diag <- tmp2[[1]]

      # Afficher avec indentation correcte
      k <- .vbse(
        paste0("c) Individual group normality diagnostic:\n",
               "\t\tNormality check with variable tolerance based on group size.\n",
               "\t\tSkewness/Kurtosis ≤1/4.5 (n≥20); ≤1.5/5 (n≥30); ≤2/6.5 (n≥50).\n",
               "\t\t",
               if (check_normality_diag==TRUE) {
                 "==> Groups show acceptable normality."
               } else {
                 "==> At least one group shows extreme non-normality."
               }),
        paste0("c) Normalité diagnostique par groupe :\n",
               "\t\tContrôle de la normalité avec tolérance variable selon la taille des groupes.\n",
               "\t\tSkewness/Kurtosis ≤1/4.5 (n≥20) ; ≤1.5/5 (n≥30) ; ≤2/6.5 (n≥50).\n",
               "\t\t",
               if (check_normality_diag==TRUE) {
                 "==> Les groupes présentent une normalité acceptable."
               } else {
                 "==> Au moins un des groupes présente une non-normalité trop extrême."
               }),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
      k <- tmp2[[2]]

      # Conclusion du diagnostic
      k <- .vbse(
        paste0("--> Recommendation: Consider these diagnostics when interpreting results.\n",
               "\t\tOutliers and influential points may affect the reliability of the analysis."),
        paste0("--> Recommandation : Tenir compte de ces diagnostics lors de l'interprétation.\n",
               "\t\tLes outliers et points influents peuvent affecter la fiabilité de l'analyse."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )

      # Retour
      list(k)
    }

    if ((check_normality==FALSE)&(check_variance_equal == TRUE)) {
      # Sauvegarder le résultat du premier contrôle
      check_normality_before <- check_normality

      # ÉTAPE 7 : Contrôle plus tolérant de la normalité
      k <- .vbse(
        "More tolerant normality check in case of homoscedasticity:",
        "Contrôle plus tolérant de la normalité en cas d'homoscédasticité :",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      temp <- auto_ku_sk(x, g, check_normality=check_normality,
                             alpha=alpha, verbose = verbose, code = code, k = k)
      check_normality <- temp[[1]]
      k <- temp[[2]]
      myaov_for_diagnostic <- temp[[3]]

      # Conclusion étape 7
      if (check_normality_before == FALSE && check_normality == TRUE) {
        k <- .vbse(
          "--> Maintaining classical ANOVA (parametric situation).",
          "--> Conservation de l'ANOVA classique (situation paramétrique).",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      } else if (check_normality_before == FALSE && check_normality == FALSE) {
        k <- .vbse(
          "--> Switching to a non-parametric comparison.",
          "--> On part vers une comparaison non paramétrique.",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }

      # ÉTAPE 8 : Diagnostic ANOVA (optionnel)
      if (isTRUE(check_normality)) {
        temp_diag <- diagnostic_anova(x, g, myaov = myaov_for_diagnostic,
                                      alpha = alpha, verbose = verbose, code = code, k = k)
        k <- temp_diag[[1]]
      }
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
      # Contrôle de la distribution pour évaluer la fiabilité du test de Wilcoxon-Mann-Witney	-->
	    #========================================================
			sd_cr <- by(x,g,stats::sd,na.rm=T) ; median_cr <- by(x,g,stats::median,na.rm=T)
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
						k <- .vbse(ang,fr,verbose = verbose, code = code, k = k, cpt="on")
						pvals <- mood.test(x[g==unique(g)[1]],x[g==unique(g)[2]])$p.value
						if (!is.na(pvals) && pvals <= alpha) {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are different :\n\t==> Data need to be centered on the mediane for ansari.test(). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les médianes sont différentes :\n\t==> Les données doivent être centrées sur la médiane pour ansari.test(). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						  by(x,g,function(x){return(x-stats::median(x))})->cent_med
						  # Supprimer le warning sur les ex-aequo (très fréquent avec données discrètes)
				  ansari_result <- withCallingHandlers(
				    ansari.test(unlist(cent_med[1]),unlist(cent_med[2])),
				    warning = function(w) invokeRestart("muffleWarning")
				  )
					  pvals <- ansari_result$p.value
					  has_ties <- length(unique(c(unlist(cent_med[1]), unlist(cent_med[2])))) < length(c(unlist(cent_med[1]), unlist(cent_med[2])))
						} else {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are the same (no need to center them for ansari.test). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les médianes sont identiques (pas besoin de les centrer pour ansari.test). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						  # Supprimer le warning sur les ex-aequo (très fréquent avec données discrètes)
				  ansari_result <- withCallingHandlers(
				    ansari.test(x[g==unique(g)[1]],x[g==unique(g)[2]]),
				    warning = function(w) invokeRestart("muffleWarning")
				  )
					  pvals <- ansari_result$p.value
					  has_ties <- length(unique(x)) < length(x)
					  }
						if (pvals < alpha) {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data do not have the same variance. p-value: ", .format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les données n'ont pas la même variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						} else {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data have the same variance. p-value: ", .format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les données ont la même variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						}
  						#========================================================
  						#			NON-NORMAL		2 categories		Variances hétérogènes + NON discret
  						#========================================================
					} else if ((check_variance_equal==FALSE) & (check_discret == FALSE)) {
						# Cas : non-normalité + variances hétérogènes + 2 groupes
						# Solution : Wilcoxon-Mann-Whitney (robuste aux violations)

						# Étape 1 : Évaluation de la fiabilité avec test KS
						pval_sidak <- 1-(1-alpha)^(1/length(unique(g)))
						if (ks_result <= pval_sidak) {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Warning! Groups do not have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
						               "Please check group distributions graphically.")
						  fr <- paste0("Estimation de la fiabilité du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\t",
						              "==> Attention ! Les groupes n'ont pas la même distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
						              "(Veuillez vérifier graphiquement les distributions des groupes.)")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						} else {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Groups have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test is expected to be reliable.")
						  fr <- paste0("Estimation de la fiabilité du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\t",
						              "==> Les groupes ont la même distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						}

						# Étape 2 : Tests non paramétriques
						pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
						if (pvals <= alpha) {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> Significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non paramétriques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des différences de médianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Différences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					} else {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> No significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non paramétriques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des différences de médianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Pas de différences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					}
					} else {
						# Cas : non-normalité + variances HOMOGÈNES + données continues + 2 groupes
						# Solution : Wilcoxon-Mann-Whitney (test non-paramétrique robuste)
						# Référence: Fay & Proschan (2010). Wilcoxon-Mann-Whitney or t-test? Statistics Surveys, 4, 1-39.

						# Étape 1 : Évaluation de la fiabilité avec test KS
						pval_sidak <- 1-(1-alpha)^(1/length(unique(g)))
						if (ks_result <= pval_sidak) {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Warning! Groups do not have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
						               "Please check group distributions graphically.")
						  fr <- paste0("Estimation de la fiabilité du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\t",
						              "==> Attention ! Les groupes n'ont pas la même distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
						              "(Veuillez vérifier graphiquement les distributions des groupes.)")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						} else {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Groups have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test is expected to be reliable.")
						  fr <- paste0("Estimation de la fiabilité du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\t",
						              "==> Les groupes ont la même distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						}

						# Étape 2 : Tests non paramétriques
						pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
						if (!is.na(pvals) && pvals <= alpha) {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> Significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non paramétriques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des différences de médianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Différences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					} else {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> No significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non paramétriques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des différences de médianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Pas de différences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					}
					}
			} else { 											# > 2 categories
			  ###################################################
			  #			NON-NORMAL		k>2 categories
			  ###################################################
					##########################################
					#   DIAGNOSTIC PRAGMATIQUE
					#   Décision du test non-paramétrique
					##########################################
					# Calculer les métriques nécessaires

					# 0. Test de Fligner (nécessaire pour le diagnostic)
					pvals3 <- fligner.test(x,g)$p.value
					if (is.na(pvals3)) {
					  .exit("Fligner-Killeen test [fligner.test()] - Error, return NA.",
						"Test de Fligner-Killeen [fligner.test()] - Erreur, retourne NA.",
						verbose=TRUE, return=TRUE)
					}

					# Calculer ks_result pour plus tard (test KS)
					temp <- pairwise(x,g,type="ks",silent=silent,boot=FALSE)$p.value
					ks_result <- min(unlist(temp),na.rm=TRUE)

					# 1. Asymétrie (skewness) GROUPE PAR GROUPE (prendre le max)
					# Utiliser agricolae::skewness qui est déjà dans les dépendances
					skew_by_group <- as.vector(unlist(
					  by(x, g, function(z) abs(agricolae::skewness(z)))
					))
					skew_val <- max(skew_by_group, na.rm = TRUE)

					# 2. Ratio des variances max/min
					var_groups <- by(x, g, var, na.rm = TRUE)
					var_ratio <- max(var_groups, na.rm = TRUE) / min(var_groups, na.rm = TRUE)

					# 3. Outliers extrêmes (déjà calculés plus haut)
					has_extreme_outliers <- exists("max_outliers_extrem") && max_outliers_extrem > 0.02  # > 2%

					# 4. Distance de Cook
					cook_dist <- cooks.distance(myaov)
					max_cook <- max(cook_dist, na.rm = TRUE)
					seuil_cook <- 4 / length(x)
					cook_influential <- max_cook > seuil_cook

					# Affichage du diagnostic pragmatique
					k <- .vbse(
					  "PRAGMATIC DIAGNOSTIC to decide non-parametric test:",
					  "DIAGNOSTIC PRAGMATIQUE pour décider du test non paramétrique :",
					  verbose = verbose, code = code, k = k, cpt = "on"
					)

					# Génération code pour mode code=TRUE
					if (isTRUE(code)) {
					  .code_anova(8, "DIAGNOSTIC PRAGMATIQUE", c(
					    "# Skewness max par groupe",
					    "skew_by_group <- by(x, g, function(z) abs(agricolae::skewness(z)))",
					    "skew_val <- max(unlist(skew_by_group), na.rm = TRUE)",
					    "",
					    "# Ratio variances",
					    "var_groups <- by(x, g, var, na.rm = TRUE)",
					    "var_ratio <- max(var_groups, na.rm = TRUE) / min(var_groups, na.rm = TRUE)",
					    "",
					    "# Fligner-Killeen",
					    "pvals3 <- fligner.test(x, g)$p.value",
					    "",
					    "# Cook's distance",
					    "cook_dist <- cooks.distance(myaov)",
					    "max_cook <- max(cook_dist, na.rm = TRUE)",
					    "seuil_cook <- 4 / length(x)",
					    "",
					    "# Outliers extrêmes",
					    "library(rstatix)",
					    "outliers_by_group <- by(x, g, function(z) {",
					    "  result <- identify_outliers(data.frame(z))",
					    "  sum(result$is.extreme, na.rm = TRUE) / length(z)",
					    "})",
					    "max_outliers_extrem <- max(unlist(outliers_by_group))"
					  ))
					}

					# Afficher les métriques
					k <- .vbse(
					  paste0("a) Skewness (asymmetry) - max |sk| across groups: ", round(skew_val, 2)),
					  paste0("a) Skewness (asymétrie) - max |sk| par groupe : ", round(skew_val, 2)),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					k <- .vbse(
					  paste0("b) Variance ratio (max/min): ", round(var_ratio, 2)),
					  paste0("b) Ratio de variances (max/min) : ", round(var_ratio, 2)),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					k <- .vbse(
					  paste0("c) Fligner-Killeen test [fligner.test()] - p-value: ", .format_pval(pvals3)),
					  paste0("c) Test de Fligner-Killeen [fligner.test()] - p-value : ", .format_pval(pvals3)),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					# d) Cook : expliquer franchissement du seuil
					cook_status_en <- if (cook_influential) {
					  paste0(" ==> INFLUENTIAL (exceeds threshold)")
					} else {
					  paste0(" ==> OK (below threshold)")
					}
					cook_status_fr <- if (cook_influential) {
					  paste0(" ==> INFLUENT (dépasse le seuil)")
					} else {
					  paste0(" ==> OK (sous le seuil)")
					}

					k <- .vbse(
					  paste0("d) Cook's distance [cooks.distance()] - max: ", round(max_cook, 4),
					         " (threshold: ", round(seuil_cook, 4), ")", cook_status_en),
					  paste0("d) Distance de Cook [cooks.distance()] - max : ", round(max_cook, 4),
					         " (seuil : ", round(seuil_cook, 4), ")", cook_status_fr),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					if (exists("max_outliers_extrem")) {
					  k <- .vbse(
					    paste0("e) Extreme outliers [rstatix::identify_outliers] - max across groups: ",
					           round(max_outliers_extrem * 100, 1), "%"),
					    paste0("e) Outliers extrêmes [rstatix::identify_outliers] - max par groupe : ",
					           round(max_outliers_extrem * 100, 1), "%"),
					    verbose = verbose, code = code, k = k, cpt = "off"
					  )
					}

					# Alpha corrigé de Šidák pour le test KS
					alpha_sidak <- 1 - (1 - alpha)^(1/length(unique(g)))

					# f) Test KS pour vérifier la comparabilité des formes
					k <- .vbse(
					  paste0("f) Kolmogorov-Smirnov test [ks.test()] - min p-value: ", .format_pval(ks_result),
					         " (α_Šidák = ", round(alpha_sidak, 3), ")"),
					  paste0("f) Test de Kolmogorov-Smirnov [ks.test()] - min p-value : ", .format_pval(ks_result),
					         " (α_Šidák = ", round(alpha_sidak, 3), ")"),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					# Décision séquentielle selon arbre de décision

					recommended_test <- NULL
					decision_reason <- NULL

					# ÉTAPE 1 : Outliers extrêmes > 5% ?
					if (has_extreme_outliers) {
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> EXTREME OUTLIERS detected (> 5%)\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> OUTLIERS EXTRÊMES détectés (> 5%)\n\t",
					               "==> Test recommandé : med1way() [ANOVA sur médianes]")
					  )
					# ÉTAPE 2 : Asymétrie forte (|skew| ≥ 2) ?
					} else if (skew_val >= 2) {
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> STRONG ASYMMETRY detected (|skew| = ", round(skew_val, 2), ")\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> ASYMÉTRIE FORTE détectée (|skew| = ", round(skew_val, 2), ")\n\t",
					               "==> Test recommandé : med1way() [ANOVA sur médianes]")
					  )
					# ÉTAPE 3 : Formes des distributions similaires (KS test avec Šidák) ?
					} else if (ks_result > alpha_sidak) {
					  # Formes similaires → Kruskal-Wallis (VALIDE même si variances différentes)
					  recommended_test <- "kruskal"
					  decision_reason <- list(
					    en = paste0("--> SIMILAR DISTRIBUTION SHAPES (KS p = ", .format_pval(ks_result), " > α_Šidák)\n\t",
					               "==> Recommended test: Kruskal-Wallis [rank-based test]\n\t",
					               "    (valid even with heteroscedasticity if shapes are similar)"),
					    fr = paste0("--> FORMES DE DISTRIBUTIONS SIMILAIRES (KS p = ", .format_pval(ks_result), " > α_Šidák)\n\t",
					               "==> Test recommandé : Kruskal-Wallis [test sur rangs]\n\t",
					               "    (valide même avec hétéroscédasticité si formes similaires)")
					  )
					# ÉTAPE 4 : Formes différentes → vérifier hétéroscédasticité
					} else if (pvals3 <= 0.05) {
					  # Formes différentes ET hétéroscédasticité → t1way
					  recommended_test <- "t1way"
					  decision_reason <- list(
					    en = paste0("--> DIFFERENT SHAPES (KS p ≤ α_Šidák) + HETEROSCEDASTICITY (Fligner p = ", .format_pval(pvals3), ")\n\t",
					               "==> Recommended test: t1way() [ANOVA on trimmed means]"),
					    fr = paste0("--> FORMES DIFFÉRENTES (KS p ≤ α_Šidák) + HÉTÉROSCÉDASTICITÉ (Fligner p = ", .format_pval(pvals3), ")\n\t",
					               "==> Test recommandé : t1way() [ANOVA sur moyennes tronquées]")
					  )
					} else {
					  # Formes différentes mais variances homogènes → med1way
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> DIFFERENT SHAPES (KS p ≤ α_Šidák) but homogeneous variances\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> FORMES DIFFÉRENTES (KS p ≤ α_Šidák) mais variances homogènes\n\t",
					               "==> Test recommandé : med1way() [ANOVA sur médianes]")
					  )
					}

					k <- .vbse(
					  decision_reason$en,
					  decision_reason$fr,
					  verbose = verbose, code = code, k = k, cpt = "off"
					)
					if ((return==TRUE) | (verbose==TRUE)) {
					  ############
					  # 	Allons jusqu'au bout du raisonnement
					  # NOUVELLE LOGIQUE : Utiliser la recommandation du diagnostic pragmatique
					  # (remplace l'ancien arbre de décision)
					  ############

					  # Variable pour stocker le test choisi (basée sur diagnostic pragmatique)
					  chosen_test <- recommended_test

					  # Exécution selon le test recommandé
					  if (chosen_test == "med1way") {
						#################
						#	Forte asymétrie / queues lourdes → med1way()
						#################
						check_variance_equal <- FALSE

						# Génération du code pour med1way
						if (isTRUE(code)) {
						  .code_anova(9, "ANOVA sur médianes (outliers/asymétrie)", c(
						    "library(WRS2)",
						    "med1way(x ~ g)"
						  ))
						}

						pvals3 <- med1way(x~g)$p.value
						if (is.na(pvals3)) {
						  ang <- paste0("Oneway ANOVA of medians [med1way()] - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.")
						  fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] - Échec, retourne NA. Le test de Kruskal-Wallis doit être utilisé pour l'interprétation.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  if (pvals3 <= alpha) {
						    ang <- paste0("Oneway ANOVA of medians [med1way()]:\n\t==> Significant differences between the median of groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
						    fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] :\n\t==> Différences significatives entre les médianes des groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  } else {
						    ang <- paste0("Oneway ANOVA of medians [med1way()]:\n\t==> Non-significant differences between the median of groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
						    fr <- paste0("Analyse de la variance à un facteur des médianes [med1way()] :\n\t==> Différences non significatives entre les médianes des groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  }
						}
					  } else if (chosen_test == "t1way") {
						#################
						#	Asymétrie modérée + problèmes variance/influence → t1way()
						#################
						check_variance_equal <- FALSE

						# Génération du code pour t1way
						if (isTRUE(code)) {
						  .code_anova(9, "ANOVA sur moyennes tronquées (variances hétérogènes)", c(
						    "library(WRS2)",
						    "t1way(x ~ g)"
						  ))
						}

						pvals3 <- t1way(x~g)$p.value
						if (is.na(pvals3)) {
						    ang <- paste0("Oneway ANOVA on trimmed means [t1way()] - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.")
							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] - Échec, retourne NA. Le test de Kruskal-Wallis doit être utilisé pour l'interprétation.")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  if (pvals3 <= alpha) {
							ang <- paste0("Oneway ANOVA on trimmed means [t1way()]:\n\t==> Significant differences between the trimmed groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] :\n\t==> Différences significatives entre les groupes tronqués.\n\t==> p-value : ", .format_pval(pvals3), ".")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
							if (!is.na(pvals) && pvals > alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronquées donnent des résultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "off")
							}
						  } else if (pvals3 > alpha) {
							ang <- paste0("One-way ANOVA on trimmed means [t1way()]:\n\t==> Non-significant differences between the trimmed groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
							fr <- paste0("Analyse de la variance à un facteur sur les moyennes tronquées [t1way()] :\n\t==> Différences non significatives entre les groupes tronqués.\n\t==> p-value : ", .format_pval(pvals3), ".")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
							if (!is.na(pvals) && pvals <= alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronquées donnent des résultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "off")
							}
						  }
						}
					  } else {
						#################
						#	Formes similaires → Kruskal-Wallis (test sur rangs)
						#################
						check_variance_equal <- TRUE  # Pour le routage post-hoc

						# Génération du code pour Kruskal-Wallis
						if (isTRUE(code)) {
						  .code_anova(9, "Test de Kruskal-Wallis (formes similaires)", c(
						    "kruskal.test(x, g)"
						  ))
						}

						# Exécuter Kruskal-Wallis
						pvals3 <- kruskal.test(x, g)$p.value
						if (is.na(pvals3)) {
						  ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()] - Failed, return NA.")
						  fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] - Échec, retourne NA.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  if (pvals3 <= alpha) {
						    ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Significant differences between groups.\n\t==> p-value: ", .format_pval(pvals3), ".")
						    fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Différences significatives entre les groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  } else {
						    ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Non-significant differences between groups.\n\t==> p-value: ", .format_pval(pvals3), ".")
						    fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Différences non significatives entre les groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  }
						}
					  } # fin choix test selon diagnostic pragmatique

					  # Mettre à jour pvals avec la p-value du test choisi (pour le retour)
					  pvals <- pvals3

					  # Affichage du test KS après med1way/t1way (avant les posthocs)
					  if (isTRUE(verbose)) {
					  if (ks_result <= pval) {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
  						              "Warning! the groups do not have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t",
  						              "(comparing with Sidak corrected alpha ", .format_pval(pval), ")\n\t",
  						              "==> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
  						              "\t(Please check graphically the groups distributions.)\n\t",
  						              "--> Adding Brunner-Munzel test as robust alternative.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\t",
  						             "Attention ! les groupes n'ont pas la même distribution. min(p-value) : ", .format_pval(ks_result), "\n\t",
  						             "(en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), ")\n\t",
  						             "==> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
  						             "\t(Veuillez vérifier graphiquement les distributions des groupes.)\n\t",
  						             "--> Ajout du test de Brunner-Munzel comme alternative robuste.")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
					  } else {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data.\n\t",
  						              "==> The groups have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t",
  						              "(comparing with Sidak corrected alpha ", .format_pval(pval), ")\n\t",
  						              "--> The Mann-Whitney-Wilcoxon test is reliable a priori.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites.\n\t",
  						             "==> Les groupes ont la même distribution. min(p-value) : ", .format_pval(ks_result), "\n\t",
  						             "(en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), ")\n\t",
  						             "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
					  }
					  } # Fin affichage KS
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
		  #ang <- paste0("Two categories.")
			#fr <- paste0("Deux catégories.")
			#k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)

		  if (isTRUE(paired)) {
			###################################################
			#			NORMAL		2 categories	APPARIÉ
			###################################################
			# Pour données appariées, pas de test de variance - t-test apparié directement
			lv <- levels(droplevels(factor(g)))
			pvals <- t.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
			if (isTRUE(code)) {
				.code_anova(8, "Test de Student apparié", c(
					"lv <- levels(droplevels(factor(g)))",
					"t.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
				))
			}
			if (isTRUE(verbose)) {
				if (pvals <= alpha) {
				  ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
				  fr <- paste0("Test de Student apparié [t.test(paired=TRUE)] :\n\t==> Différences significatives entre les conditions. p-value : ", .format_pval(pvals))
				} else {
				  ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
				  fr <- paste0("Test de Student apparié [t.test(paired=TRUE)] :\n\t==> Différences non significatives entre les conditions. p-value : ", .format_pval(pvals))
				}
				k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
			}

		  } else {
			###################################################
			#			NORMAL		2 categories	NON-APPARIÉ
			###################################################
			fm <- formula(x~g)
			pvals <- var.test(fm)$p.value
			if (pvals>alpha) {
			  ###################################################
			  #			NORMAL		2 categories	homogene variance
			  ###################################################
			  pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=TRUE)$p.value
			  if (isTRUE(code)) {
				  .code_anova(9, "Test de Student (variances homogènes)", c(
					  "levels_g <- levels(factor(g))",
					  "t.test(x[g==levels_g[1]], x[g==levels_g[2]], var.equal = TRUE, paired = FALSE)"
				  ))
			  }
			  if (isTRUE(verbose)) {
				  if (pvals <= alpha) {
					##################################
					#        NORMAL, Student significant (var.equal=T)
					##################################
					ang <- paste0("Student test [t.test()] :\n\t==> Significant differences between groups. p-value: ",  .format_pval(pvals))
					fr <- paste0("Test de Student [t.test()] :\n\t==> Différences significatives entre les groupes. p-value : ", .format_pval(pvals))
					k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
				  } else {
					##################################
					#        NORMAL, Student unsignificant (var.equal=T)
					##################################
					ang <- paste0("Student test [t.test()] - Non-significant differences between groups. p-value: ",  .format_pval(pvals))
					fr <- paste0("Test de Student [t.test()] - Différences non significatives entre les groupes. p-value : ", .format_pval(pvals))
					k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
				  }
			  }
			} else {
			  ###################################################
			  #			NORMAL		2 categories	non-homogene variance
			  ###################################################
			  check_variance_equal <-  FALSE
			  pvals_var <- pvals
			  pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=FALSE)$p.value
			  if (isTRUE(code)) {
				  .code_anova(9, "Test de Student-Welch (variances hétérogènes)", c(
					  "levels_g <- levels(factor(g))",
					  "t.test(x[g==levels_g[1]], x[g==levels_g[2]], var.equal = FALSE, paired = FALSE)"
				  ))
			  }
			  ang <- paste0("Tests for comparing two normal groups with heterogeneous variances.\n\t",
						   "a) Analysis of mean differences by bootstrap [pairwise.boot(mu='mean') from {KefiR}]\n\t",
						   "b) Welch test [t.test(var.equal=FALSE)].\n\t",
						   if (pvals <= alpha) "\t==> Significant differences between groups (p = " else "\t==> No significant differences between groups (p = ",
						   .format_pval(pvals), ").")
			  fr <- paste0("Tests de comparaison de deux groupes normaux à variances hétérogènes.\n\t",
						  "a) Analyse des différences de moyennes par bootstrap [pairwise.boot(mu='mean') de {KefiR}]\n\t",
						  "b) Test de Welch [t.test(var.equal=FALSE)].\n\t",
						  if (pvals <= alpha) "\t==> Différences de moyennes significatives entre groupes (p = " else "\t==> Aucune différence de moyennes significative entre groupes (p = ",
						  .format_pval(pvals), ").")
			  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
			  # NOTE: Le message ci-dessus contient déjà la conclusion (significatif ou non)
			  # Pas besoin de message supplémentaire - évite la redondance
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
			# Réutilisation de myaov créé à l'étape 4
			pvals <- summary(myaov)[[1]][["Pr(>F)"]][1]
			if (isTRUE(code)) {
				.code_anova(9, "ANOVA à un facteur (variances homogènes)", c(
					"summary(myaov)"
				))
			}
			ang <- paste0("One-way ANOVA [aov()].\n\t",
			             "Result: p = ", .format_pval(pvals),
			             if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
			             if (pvals <= alpha) "==> Significant differences between groups." else "==> No significant differences between groups.")
			fr <- paste0("ANOVA à un facteur [aov()].\n\t",
			            "Résultat : p = ", .format_pval(pvals),
			            if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
			            if (pvals <= alpha) "==> Différences de moyennes significatives entre groupes." else "==> Aucune différence de moyennes significative entre groupes.")
			k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
		  } else {												# Non-identical variances
			  ##################################
			  #        NORMAL, >2 categories var.equal=F => FANOVA.HETERO
			  ##################################
			pvals <- oneway.test(x~g,var.equal=FALSE)$p.value
			if (isTRUE(code)) {
				.code_anova(9, "Test de Welch (variances hétérogènes)", c(
					"oneway.test(x ~ g, var.equal = FALSE)"
				))
			}
			myf <- try(fanova.hetero(data.frame(x,g = as.factor(g)),x~g),silent=silent)
			if (inherits(myf,"try-error")) {
			  #if (isTRUE(verbose)) {cat("Error on fanova.hetero()\n")}
			  pvals2 <- alpha
			} else {pvals2 <- myf$ans[4]}
			if (isTRUE(verbose)) {
				ang <- paste0("Welch's heteroscedastic F-test [oneway.test()].\n\t",
				             "Result: p = ", .format_pval(pvals),
				             if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
				             if (pvals <= alpha) "==> Significant differences between groups." else "==> No significant differences between groups.")
				fr <- paste0("Test F hétéroscédastique de Welch [oneway.test()].\n\t",
				            "Résultat : p = ", .format_pval(pvals),
				            if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
				            if (pvals <= alpha) "==> Différences de moyennes significatives entre groupes." else "==> Aucune différence de moyennes significative entre groupes.")
				k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)

				# Check fanova.hetero consistency
				if ((pvals <= alpha && pvals2 > alpha) || (pvals > alpha && pvals2 <= alpha)) {
					ang <- paste0("\tNote: fanova.hetero() suggests distribution differences (p = ", .format_pval(pvals2), ").")
					fr <- paste0("\tNote : fanova.hetero() suggère des différences de distributions (p = ", .format_pval(pvals2), ").")
					k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "off")
				}
			}


		  }
		}
	  }

	  ################################################################################
	  # MODE CODE=TRUE
	  ################################################################################

	  # NOTE: Code generation for one-factor analyses is handled by individual
	  # functions (.normality.R, .variance.R, .posthoc.R) which output code
	  # synchronized with their verbose messages. This ensures correspondence
	  # between explanations and code. No centralized code_str needed here.

	  # NOTE PERSO: Ajouter p-value globale dans le bilan pour return=FALSE
	  # (cf. Cahier des charges priorité 2)
	  # NOTE: chosen_test indique quel test non-paramétrique a été utilisé (med1way, t1way, kruskal, ou NULL)
	  # chosen_test est initialisé à NULL au début de la fonction et mis à jour dans le diagnostic pragmatique
	  return(check <- list(x = x, g = g, check_normality = check_normality, check_variance_equal = check_variance_equal, k = k, global_pvalue = pvals, chosen_test = chosen_test))
 }
