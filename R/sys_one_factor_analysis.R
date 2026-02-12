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
#' @param return Logical. If TRUE, returns results; if FALSE, only displays output (default is TRUE).
#' @param k Optional. Counter for verbose messages.
#' @param code Logical. If TRUE, displays R code snippets (default is FALSE).
#' @param debug Logical. If TRUE, detailed debugging messages are displayed (default is FALSE).
#' @param verbose Logical. If TRUE, provides detailed output during the analysis (default is FALSE).
#' @param boot Logical. If TRUE, enables bootstrap-based procedures (default is TRUE).
#' @param silent Logical. If TRUE, suppresses warnings from underlying functions (default is TRUE).
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
#' \dontrun{
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
#' }
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
#' @keywords internal
.one_factor_analysis <- function(x=NULL,g=NULL,formula=NULL,data=NULL,
	paired = FALSE, id = NULL,alpha = 0.05,return=TRUE,
	k=NULL,code=FALSE,debug=FALSE,verbose=FALSE, boot = TRUE, silent=TRUE) {

  # Protection contre code=NULL (si utilisateur a utilis\u00e9 code=F avec variable F dans donn\u00e9es)
  if (is.null(code)) {
    warning("Le param\u00e8tre 'code' \u00e9tait NULL. Utilisez code=FALSE ou code=TRUE (pas F ou T). R\u00e9initialisation \u00e0 FALSE.")
    code <- FALSE
}

  # Fonction auxiliaire locale pour afficher code R avec num\u00e9rotation (ANOVA simple)
  .code_anova <- function(step_num, title, code_lines) {
    cat(paste0("# ", step_num, ") ", title, "\n"))
    for (line in code_lines) {
      cat(paste0(line, "\n"))
    }
    cat("\n")
  }

  # Si code==TRUE, d\u00e9sactiver verbose MAIS garder k synchronis\u00e9
  if (isTRUE(code)) {
    verbose_original <- verbose
    verbose <- FALSE
    # NOTE: On utilisera k au lieu de k pour synchroniser la num\u00e9rotation
  }

  ################################################################
  #
  #
  #		Si paired = TRUE et 2 groupes (k  = 2)
  #
  #
  ################################################################
  ## --- V\u00e9rifs d\u2019int\u00e9grit\u00e9 pour le mode appari\u00e9 ------------------------------
  if (isTRUE(paired)) {
    # Si identifiants fourni, le contr\u00f4ler et r\u00e9aligner les donn\u00e9es dessus
    if(!is.null(id)) {
      g2 <- droplevels(factor(g))
      dfc <- data.frame(id = data[id], g = g2, x = as.numeric(x))
      tab <- table(dfc$id, dfc$g)              # matrice (#id x 2)
      # 1) chaque id doit \u00eatre pr\u00e9sent dans les 2 modalit\u00e9s
      miss <- rownames(tab)[rowSums(tab > 0) < 2L]
      if (length(miss)) {
        .exit(
          "Some ids are incomplete (no match for each modality of 'g').",
          "Certains identifiants de id sont incomplets (pas de correspondance pour chaque modalit\u00e9 de 'g').",
          verbose = verbose, code = code, return = return
        )
      }
      # 2) (option strict) exactement une obs par id et par modalit\u00e9
      if (any(tab != 1L)) {
        .exit(
          "Each id must have exactly one observation per level of 'g' (no duplicates).",
          "Chaque identifiant doit avoir exactement une observation par modalit\u00e9 de 'g' (pas de doublons).",
          verbose = verbose, code = code, return = return
        )
      }
      # R\u00e9aligner
      k <- .vbse("Alignment by id for paired (2 levels).",
            "Alignement par id pour les donn\u00e9es appari\u00e9es.",
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
      # TEST T APPARI\u00c9 (ou Wilcoxon sign\u00e9) - Ajout\u00e9 pour afficher le r\u00e9sultat du test
      # =============================================================================
      lv <- levels(droplevels(factor(g)))
      if (isTRUE(check_normality)) {
        # Diff\u00e9rences normales => t-test appari\u00e9
        pvals <- t.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
        if (isTRUE(code)) {
          .code_anova(8, "Test de Student appari\u00e9", c(
            "lv <- levels(droplevels(factor(g)))",
            "t.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
          ))
        }
        if (isTRUE(verbose)) {
          if (pvals <= alpha) {
            ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Student appari\u00e9 [t.test(paired=TRUE)] :\n\t==> Diff\u00e9rences significatives entre les conditions. p-value : ", .format_pval(pvals))
          } else {
            ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Student appari\u00e9 [t.test(paired=TRUE)] :\n\t==> Diff\u00e9rences non significatives entre les conditions. p-value : ", .format_pval(pvals))
          }
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
        }
      } else {
        # Diff\u00e9rences non normales => Wilcoxon sign\u00e9
        pvals <- wilcox.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
        if (isTRUE(code)) {
          .code_anova(8, "Test de Wilcoxon sign\u00e9 (donn\u00e9es appari\u00e9es non normales)", c(
            "lv <- levels(droplevels(factor(g)))",
            "wilcox.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
          ))
        }
        if (isTRUE(verbose)) {
          if (pvals <= alpha) {
            ang <- paste0("Wilcoxon signed-rank test [wilcox.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Wilcoxon sign\u00e9 [wilcox.test(paired=TRUE)] :\n\t==> Diff\u00e9rences significatives entre les conditions. p-value : ", .format_pval(pvals))
          } else {
            ang <- paste0("Wilcoxon signed-rank test [wilcox.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
            fr <- paste0("Test de Wilcoxon sign\u00e9 [wilcox.test(paired=TRUE)] :\n\t==> Diff\u00e9rences non significatives entre les conditions. p-value : ", .format_pval(pvals))
          }
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
        }
      }

      return(check <- list(x = x, g = g, check_normality = check_normality,
                           check_variance_equal = check_variance_equal, k = k,
                           global_pvalue = pvals, chosen_test = NULL))
    }
  } # Fin du sc\u00e9nario paired sachant que l'alignement si n_g > 2 (k) passe en .multi_factor_analysis	-->

  #==================================================
  #         Initiation des variables de contr\u00f4le
  #==================================================
  check_discret <- discret.test(x)
  check_normality <- c()
  check_number <- c()
  check_variance_equal <-  c()
  pvals <- NA  # Initialize pvals to avoid "object not found" errors at return
  chosen_test <- NULL  # Initialize: NULL=parametric, "med1way", "t1way", or "kruskal" for non-parametric

  #==================================================
  #         Contr\u00f4le du nombre de groupes
  #==================================================

  check_number <- length(unique(g[!is.na(g)]))
  if (check_number==2) { 							# 2 categories
    if (paired==TRUE) {
      k <- .vbse("Two conditions.",
                 "Deux conditions.",
                 verbose = verbose, code = code, k = k, cpt="on")
      if (isTRUE(code)) {
        .code_anova(1, "Deux conditions (donn\u00e9es appari\u00e9es)", c(
          "length(unique(g))  # Doit \u00eatre 2"
        ))
      }
    } else{
      k <- .vbse("Two groups.",
                 "Deux groupes.",
                 verbose = verbose, code = code, k = k, cpt="on")
      if (isTRUE(code)) {
        .code_anova(1, "Deux groupes ind\u00e9pendants", c(
          "length(unique(g))  # Doit \u00eatre 2"
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
  # Contr\u00f4le de l'ind\u00e9pendance des observations.
  #==================================================
  .dbg("Check of observation independence.",
       "Contr\u00f4le de l'ind\u00e9pendance des observations.",debug=debug)
  if (paired == FALSE) {
    # Donn\u00e9es non appari\u00e9es - ASSOMPTION 1/3 pour ANOVA 1 facteur
    k <- .vbse(paste0("ASSUMPTION 1/3: Independence check - Ensure that observations between groups are independent.\n",
                      "\tVerify that:\n",
                      "\t  \u2022 No repeated measures\n",
                      "\t  \u2022 No clustering effects\n",
                      "\t  \u2022 No carry-over effects"),
               paste0("ASSOMPTION 1/3 : Contr\u00f4le de l'ind\u00e9pendance - v\u00e9rifiez que les observations entre les groupes sont ind\u00e9pendantes.\n",
                      "\tV\u00e9rifiez que :\n",
                      "\t  \u2022 Pas de mesures r\u00e9p\u00e9t\u00e9es\n",
                      "\t  \u2022 Pas d'effet cluster\n",
                      "\t  \u2022 Pas d'effets report [influence sur l'ordre des mesures]"),
               verbose = verbose, code = code, k = k, cpt = "on")
    if (isTRUE(code)) {
      .code_anova(2, "Contr\u00f4le ind\u00e9pendance (groupes ind\u00e9pendants)", c())
    }
  } else {
    # Donn\u00e9es appari\u00e9es - ASSOMPTION 1/3 pour ANOVA 1 facteur appari\u00e9
    k <- .vbse("ASSUMPTION 1/3: Paired independence check - Ensure that the dependence structure is only due to the intended pairing\n\t* e.g., same subject before/after\n\t* no additional clustering between paired observations\n\t* no carry-over effects between paired observations",
               "ASSOMPTION 1/3 : Contr\u00f4le de l'ind\u00e9pendance pour donn\u00e9es appari\u00e9es - V\u00e9rifiez que la structure de d\u00e9pendance provienne uniquement de l'appariement pr\u00e9vu\n\t* ex. : m\u00eame sujet avant/apr\u00e8s\n\t* pas de clustering suppl\u00e9mentaire entre observations appari\u00e9es\n\t* pas d'effets report suppl\u00e9mentaires entre observations appari\u00e9es",
               verbose = verbose, code = code, k = k, cpt = "on")
    if (isTRUE(code)) {
      .code_anova(2, "Contr\u00f4le ind\u00e9pendance (donn\u00e9es appari\u00e9es)", c())
    }
  }
  #==================================================
  # D\u00e9tection des Outlier en X
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
    #  Gestion des donn\u00e9es extr\u00eames si non appari\u00e9s
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
    # Messages gradu\u00e9s selon la s\u00e9v\u00e9rit\u00e9
    if (max_outliers_extrem > 0) {
      k <- .vbse(
        paste0("Outlier detection [identify_outliers() {rstatix}]\n\t",
               "Result: ", round(max_outliers_extrem * 100, 1), "% EXTREME outliers (maximum in a group)\n\t",
               "==> Extreme values may strongly influence parametric tests"),
        paste0("D\u00e9tection outliers [identify_outliers() {rstatix}]\n\t",
               "R\u00e9sultat : ", round(max_outliers_extrem * 100, 1), "% outliers EXTR\u00caMES (maximum dans un groupe)\n\t",
               "==> Valeurs extr\u00eames peuvent fortement influencer tests param\u00e9triques."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      if (isTRUE(code)) {
        .code_anova(3, "D\u00e9tection outliers extr\u00eames (groupes ind\u00e9pendants)", c(
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
        paste0("D\u00e9tection outliers [identify_outliers() {rstatix}]\n\t",
               "R\u00e9sultat : ", round(max_outliers_base * 100, 1), "% outliers (maximum dans un groupe, non extr\u00eames)\n\t",
               "==> Valeurs \u00e0 surveiller, mais non excessives"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      if (isTRUE(code)) {
        .code_anova(3, "D\u00e9tection outliers mod\u00e9r\u00e9s (groupes ind\u00e9pendants)", c(
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
    # Donn\u00e9es appari\u00e9es \u2014 contr\u00f4le des outliers sur les diff\u00e9rences intra-sujet
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
          paste0("Outliers extr\u00eames dans les diff\u00e9rences appari\u00e9es [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d'outliers extr\u00eames : ", round(prop_ext * 100, 1), "%.\n\t",
                 "Ces paires peuvent fortement influencer les tests param\u00e9triques appari\u00e9s."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        if (isTRUE(code)) {
          .code_anova(3, "D\u00e9tection outliers extr\u00eames (diff\u00e9rences appari\u00e9es)", c(
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
          paste0("Outliers dans les diff\u00e9rences appari\u00e9es [identify_outliers() de {rstatix}].\n\t",
                 "Proportion d'outliers : ", round(prop_out * 100, 1), "%.\n\t",
                 "\u00c0 surveiller mais non extr\u00eames."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        if (isTRUE(code)) {
          .code_anova(3, "D\u00e9tection outliers mod\u00e9r\u00e9s (diff\u00e9rences appari\u00e9es)", c(
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
        "Donn\u00e9es appari\u00e9es avec plus de deux modalit\u00e9s : contr\u00f4le des outliers sur les diff\u00e9rences non d\u00e9fini ici.",
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
  }
	  #==================================================
	  #
    #   Contr\u00f4le de la normalit\u00e9
    #
	  #		On notera que formula est ignor\u00e9e par cette fonction (\u00e0 dev si n\u00e9cessaire)
	  #
	  #==================================================
    .dbg("Normality check with .normality().",
         "Contr\u00f4le de normalit\u00e9 avec .normality().",debug=debug)
    skew <- function(vector) {return(abs(skewness(vector)))}
    skew2 <- function(vector) {return(skewness.norm.test(vector)$p.value)}
    kurto <- function(vector) {if (is.na(abs(kurtosis(vector)))){return(10)} ; return(abs(kurtosis(vector)))}
    kurto2 <- function(vector) {return(kurtosis.norm.test(vector)$p.value)}
    ##########################
    # Correction de Sidak
    ##########################
    # Sidak devient trop conservatif au-del\u00e0 de 10 groupes	-->
    # Bender, R., & Lange, S. (2001). Adjusting for multiple testing\u2014when and how?  DOI: 10.1016/S0895-4356(00)00314-0
    # Il est justifi\u00e9 de passer ensuite au bootstrap
    # M\u00e9thode FDR la correction par FDR (False Discovery Rate) \u2013 notamment via la m\u00e9thode Benjamini-Hochberg (BH)
    # Murray, E. J., Berrie, L., & Matthews, J. N. S. (2021). Understanding multiplicity in clinical trials: the impact of multiple endpoints in the interpretation of treatment effects.
    # Attention : FDR est dans une d\u00e9marche exploratoire ou descriptive, pas une validation r\u00e9glementaire.
    .dbg("Sidak's correction rather than Bonferroni for independent tests (conservative).",
         "Correction de Sidak plut\u00f4t que Bonferroni pour les tests ind\u00e9pendants (conservateur).",debug=debug)
    pval <- 1-(1-alpha)^(1/length(unique(g)))

    if (check_number==2) { 							# 2 categories
      # Contr\u00f4le de la normalit\u00e9 des groupes (ou des diff\u00e9rences si appari\u00e9)
      check_normality <- .normality(x, g, alpha=alpha, k=k, verbose=verbose, paired=paired)
    } else {
      # Annonce de l'ajustement du mod\u00e8le ANOVA
      k <- .vbse(
        "Fitting ANOVA model [aov()].",
        "Ajustement du mod\u00e8le ANOVA [aov()].",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      if (isTRUE(code)) {
        .code_anova(4, "Ajustement du mod\u00e8le ANOVA", c(
          "myaov <- aov(x ~ g)"
        ))
      }

      # Contr\u00f4le de la normalit\u00e9 des r\u00e9sidus
      myaov <- aov(x~g)
      n_residu <- length(x)

      # ASSOMPTION 2/3 pour ANOVA 1 facteur (1=ind\u00e9pendance, 2=normalit\u00e9, 3=variance)
      k <- .vbse(
        "ASSUMPTION 2/3: Academic check of ANOVA model residuals normality.",
        "ASSOMPTION 2/3 : Contr\u00f4le ACAD\u00c9MIQUE de la normalit\u00e9 des r\u00e9sidus du mod\u00e8le ANOVA.",
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
          # Scenario 1: n \u2264 50 -> Shapiro-Wilk only
          code_lines <- c(
            paste0("alpha_sidak <- 1 - (1 - ", alpha, ")^(1/length(unique(g)))"),
            "by(x, g, shapiro.test)"
          )
        } else if (!reajust && min_group_size <= 500) {
          # Scenario 2: 50 < n \u2264 500 -> Jarque-Bera with Sidak correction
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
        .code_anova(5, "ASSOMPTION 2 : Contr\u00f4le de la normalit\u00e9 (2 groupes)", code_lines)
      } else {
        # Detect residuals size (same logic as .normality() for single group)
        n_residuals <- length(x)

        # Generate code following .normality() logic for residuals
        if (n_residuals <= 50) {
          # Scenario 1: n \u2264 50 -> Shapiro-Wilk only
          code_lines <- c(
            "residus <- residuals(myaov)",
            "shapiro.test(residus)"
          )
        } else if (n_residuals <= 500) {
          # Scenario 2: 50 < n \u2264 500 -> Jarque-Bera
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
        .code_anova(5, "ASSOMPTION 2 : Contr\u00f4le de la normalit\u00e9 (r\u00e9sidus ANOVA)", code_lines)
      }
    }
    #==================================================
    #
    #   Contr\u00f4le de la variance
    #
    #==================================================
    .dbg("Variance check with .variance().",
         "Contr\u00f4le de la variance avec .variance().",debug=debug)
    # Pour ANOVA 1 facteur : ASSOMPTION 3/3 (1=ind\u00e9pendance, 2=normalit\u00e9, 3=variance)
    variance_temp <- .variance(x, g, check_normality=check_normality, alpha=alpha, paired=paired,
                  debug = debug, verbose=verbose, k=k, code=FALSE, assumption_label="3/3")
    check_variance_equal <- variance_temp[[1]]
    k <- variance_temp[[2]]

    # G\u00e9n\u00e9rer le code de variance ICI avec num\u00e9rotation k
    if (isTRUE(code)) {
      if (check_number == 2) {
        .code_anova(6, "Test de variance (Fisher-Snedecor)", c(
          "var.test(x~g)"
        ))
      } else {
        if (check_normality == TRUE) {
          .code_anova(6, "Test d'homog\u00e9n\u00e9it\u00e9 des variances (Bartlett)", c(
            "bartlett.test(x, g)"
          ))
        } else {
          .code_anova(6, "Test d'homog\u00e9n\u00e9it\u00e9 des variances (Levene)", c(
            "library(car)",
            "leveneTest(x, g)"
          ))
        }
      }
    }

    #==================================================
    #
    #   Contr\u00f4le d'un \u00e9ventuel retour vers param\u00e9trique
    #
    #==================================================
    .dbg("Check for a possible return to parametric methods.",
         "Contr\u00f4le d'un \u00e9ventuel retour vers param\u00e9trique.",debug=debug)
  	##########################
  	#
  	#		Auto-analyse pour envisager un retour vers la situation param\u00e9trique (Kurtosis & Skweness	-->)
  	#
  	##########################
	  # |Skewness| < 1 et |Kurtosis| < 1.5 sont les seuils recommand\u00e9s pour consid\u00e9rer une distribution comme "approximativement normale", selon Kline (2011).
	  # Kline, R. B. (2011). Principles and practice of structural equation modeling (4th ed.). New York, NY: Guilford Press.
	  # Normalit\u00e9 variables indicates acceptable normality for most measures, based on Kline's (2023) guidelines, where skewness values should ideally be within \u00b13 and kurtosis values within \u00b110.
	  # Permissivit\u00e9 extr\u00eame : skewness > 2 et kurtosis > 7 Blanca et al. (2018)
	  # Blanca, M. J., Alarc\u00f3n, R., Arnau, J., Bono, R., & Bendayan, R. (2017). Non-normal data: Is ANOVA still a valid option? Psicothema, 29(4), 552-557. https://diposit.ub.edu/dspace/bitstream/2445/122126/1/671797.pdf#:~\:text=sought%20to%20answer%20the%20following,In%20addition
    # Une violation de la normalit\u00e9 a peu d'effet sur l'ANOVA si variances homog\u00e8nes.Caldwell (2022)
    # Caldwell, A. (2022, March 31). Chapter 12 Violations of Assumptions. In Power Analysis with Superpower. https://aaroncaldwell.us/SuperpowerBook/violations-of-assumptions.html

    auto_ku_sk <- function(x, g, check_normality = TRUE, alpha = 0.05, verbose = TRUE, k = NULL, code = NULL) {
      # Coercitions
      if (is.data.frame(x)) x <- x[[1]]
      if (is.matrix(x))     x <- x[,1,drop=TRUE]
      if (is.list(x))       x <- unlist(x, use.names = FALSE)
      x <- as.vector(x)
      if (!is.numeric(x)) stop("auto_ku_sk() attend un vecteur num\u00e9rique 'x'.")

      if (is.data.frame(g)) g <- g[[1]]
      if (is.matrix(g))     g <- g[,1,drop=TRUE]
      if (is.list(g))       g <- unlist(g, use.names = FALSE)
      g <- droplevels(factor(g))
      lev <- levels(g)
      ng  <- nlevels(g)

      # Normalit\u00e9 r\u00e9sidus (tol\u00e9rante) - \u00c9TAPE 7
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

    # Fonction diagnostic ANOVA (\u00c9TAPE 8) - optionnelle pour garantir fiabilit\u00e9
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

      # En-t\u00eate DIAGNOSTIC
      k <- .vbse(
        "Optional DIAGNOSTIC to ensure ANOVA reading reliability:",
        "DIAGNOSTIC optionnel pour garantir la fiabilit\u00e9 de la lecture de l'ANOVA :",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      # a) Outliers sur R\u00c9SIDUS par groupe
      if (!is.null(myaov)) {
        residuals_data <- myaov$residuals

        # D\u00e9tecter outliers sur r\u00e9sidus
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
            paste0("a) Outliers d\u00e9tect\u00e9s sur les r\u00e9sidus :\n",
                   "\t\t==> ", n_outliers, " outliers (", prop_out, "%), ",
                   n_extreme, " extr\u00eames (", prop_ext, "%)."),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            "a) No outliers detected on residuals.",
            "a) Aucun outlier d\u00e9tect\u00e9 sur les r\u00e9sidus.",
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
                   "\t\t==> ", n_influential, " observations influentes d\u00e9tect\u00e9es.\n",
                   "\t\t(", prop_influential, "% d\u00e9passent seuil 4/n = ", round(cook_threshold, 4), ")."),
            paste0("b) Diagnostic d'influence (distance de Cook) :\n",
                   "\t\t==> ", n_influential, " observations influentes d\u00e9tect\u00e9es.\n",
                   "\t\t(", prop_influential, "% d\u00e9passent seuil 4/n = ", round(cook_threshold, 4), ")."),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            "b) Cook's distance: No influential observations detected (all < 4/n).",
            "b) Distance de Cook : Aucune observation influente d\u00e9tect\u00e9e (toutes < 4/n).",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      }

      # c) Normalit\u00e9 diagnostique par groupe
      # Appeler .normality() silencieusement pour obtenir le r\u00e9sultat
      tmp2 <- .normality(x, g, alpha = alpha, tolerance = "extrem", k = k, verbose = FALSE, code = FALSE, cpt = "off")
      check_normality_diag <- tmp2[[1]]

      # Afficher avec indentation correcte
      k <- .vbse(
        paste0("c) Individual group normality diagnostic:\n",
               "\t\tNormality check with variable tolerance based on group size.\n",
               "\t\tSkewness/Kurtosis \u22641/4.5 (n\u226520); \u22641.5/5 (n\u226530); \u22642/6.5 (n\u226550).\n",
               "\t\t",
               if (check_normality_diag==TRUE) {
                 "==> Groups show acceptable normality."
               } else {
                 "==> At least one group shows extreme non-normality."
               }),
        paste0("c) Normalit\u00e9 diagnostique par groupe :\n",
               "\t\tContr\u00f4le de la normalit\u00e9 avec tol\u00e9rance variable selon la taille des groupes.\n",
               "\t\tSkewness/Kurtosis \u22641/4.5 (n\u226520) ; \u22641.5/5 (n\u226530) ; \u22642/6.5 (n\u226550).\n",
               "\t\t",
               if (check_normality_diag==TRUE) {
                 "==> Les groupes pr\u00e9sentent une normalit\u00e9 acceptable."
               } else {
                 "==> Au moins un des groupes pr\u00e9sente une non-normalit\u00e9 trop extr\u00eame."
               }),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
      k <- tmp2[[2]]

      # Conclusion du diagnostic
      k <- .vbse(
        paste0("--> Recommendation: Consider these diagnostics when interpreting results.\n",
               "\t\tOutliers and influential points may affect the reliability of the analysis."),
        paste0("--> Recommandation : Tenir compte de ces diagnostics lors de l'interpr\u00e9tation.\n",
               "\t\tLes outliers et points influents peuvent affecter la fiabilit\u00e9 de l'analyse."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )

      # Retour
      list(k)
    }

    if ((check_normality==FALSE)&(check_variance_equal == TRUE)) {
      # Sauvegarder le r\u00e9sultat du premier contr\u00f4le
      check_normality_before <- check_normality

      # \u00c9TAPE 7 : Contr\u00f4le plus tol\u00e9rant de la normalit\u00e9
      k <- .vbse(
        "More tolerant normality check in case of homoscedasticity:",
        "Contr\u00f4le plus tol\u00e9rant de la normalit\u00e9 en cas d'homosc\u00e9dasticit\u00e9 :",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      temp <- auto_ku_sk(x, g, check_normality=check_normality,
                             alpha=alpha, verbose = verbose, code = code, k = k)
      check_normality <- temp[[1]]
      k <- temp[[2]]
      myaov_for_diagnostic <- temp[[3]]

      # Conclusion \u00e9tape 7
      if (check_normality_before == FALSE && check_normality == TRUE) {
        if (check_number == 2) {
          k <- .vbse(
            "--> Maintaining parametric comparison (Student's t-test).",
            "--> Conservation de la comparaison param\u00e9trique (test de Student).",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            "--> Maintaining classical ANOVA (parametric situation).",
            "--> Conservation de l'ANOVA classique (situation param\u00e9trique).",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      } else if (check_normality_before == FALSE && check_normality == FALSE) {
        k <- .vbse(
          "--> Switching to a non-parametric comparison.",
          "--> On part vers une comparaison non param\u00e9trique.",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }

      # \u00c9TAPE 8 : Diagnostic ANOVA (optionnel, uniquement pour n>2 groupes)
      if (isTRUE(check_normality) && check_number > 2) {
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
      # Contr\u00f4le de la distribution pour \u00e9valuer la fiabilit\u00e9 du test de Wilcoxon-Mann-Witney	-->
	    #========================================================
			sd_cr <- by(x,g,stats::sd,na.rm=T) ; median_cr <- by(x,g,stats::median,na.rm=T)
			data_cr <- x
			for (i in names(median_cr)) {
				sd_i <- sd_cr[which(names(sd_cr)==i)]
				if (is.na(sd_i) || sd_i == 0) sd_i <- 1  # Protection données constantes (ex: Poisson)
				data_cr[g==i] <- (data_cr[g==i]-median_cr[which(names(median_cr)==i)])/sd_i
			}
			ks_result <- tryCatch({
				temp <- pairwise(data_cr,g,type="ks",silent=silent,boot=FALSE)
				min(unlist(temp$p.value),na.rm=T)
			}, error = function(e) 1)  # 1 = assume same distribution (conservative)

			#========================================================
			#			NON-NORMAL		2 categories
			#========================================================
			if (check_number==2) { 							# 2 categories

			    if  (check_discret == TRUE)  {
					#========================================================
					#			NON-NORMAL		2 categories		Non acceptable for t.test()		Discret
					#========================================================
						ang <- paste0("The number of unique values [length() & unique()] suggests the presence of discrete data :\n\t",round(length(unique(x))/length(x),3)*100,"% of unique values.")
						fr <- paste0("Le nombre de valeurs uniques [length() & unique()] sugg\u00e8re la pr\u00e9sence de donn\u00e9es discr\u00e8tes :\n\t", round(length(unique(x))/length(x),3)*100, "% de valeurs uniques.")
						k <- .vbse(ang,fr,verbose = verbose, code = code, k = k, cpt="on")
						pvals <- mood.test(x[g==unique(g)[1]],x[g==unique(g)[2]])$p.value
						if (!is.na(pvals) && pvals <= alpha) {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are different :\n\t==> Data need to be centered on the mediane for ansari.test(). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les m\u00e9dianes sont diff\u00e9rentes :\n\t==> Les donn\u00e9es doivent \u00eatre centr\u00e9es sur la m\u00e9diane pour ansari.test(). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						  by(x,g,function(x){return(x-stats::median(x))})->cent_med
						  # Supprimer le warning sur les ex-aequo (tr\u00e8s fr\u00e9quent avec donn\u00e9es discr\u00e8tes)
				  ansari_result <- withCallingHandlers(
				    ansari.test(unlist(cent_med[1]),unlist(cent_med[2])),
				    warning = function(w) invokeRestart("muffleWarning")
				  )
					  pvals <- ansari_result$p.value
					  has_ties <- length(unique(c(unlist(cent_med[1]), unlist(cent_med[2])))) < length(c(unlist(cent_med[1]), unlist(cent_med[2])))
						} else {
						  ang <- paste0("Brown and Mood test [mood.test()] - The medianes are the same (no need to center them for ansari.test). p-value : ",.format_pval(pvals))
						  fr <- paste0("Test de Brown et Mood [mood.test()] - Les m\u00e9dianes sont identiques (pas besoin de les centrer pour ansari.test). p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						  # Supprimer le warning sur les ex-aequo (tr\u00e8s fr\u00e9quent avec donn\u00e9es discr\u00e8tes)
				  ansari_result <- withCallingHandlers(
				    ansari.test(x[g==unique(g)[1]],x[g==unique(g)[2]]),
				    warning = function(w) invokeRestart("muffleWarning")
				  )
					  pvals <- ansari_result$p.value
					  has_ties <- length(unique(x)) < length(x)
					  }
						if (pvals < alpha) {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data do not have the same variance. p-value: ", .format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les donn\u00e9es n'ont pas la m\u00eame variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						} else {
						  ang <- paste0("Ansari-Bradley test [ansari.test()] - Data have the same variance. p-value: ", .format_pval(pvals))
						  fr <- paste0("Test d'Ansari-Bradley [ansari.test()] - Les donn\u00e9es ont la m\u00eame variance. p-value : ", .format_pval(pvals))
						  k <- .vbse(ang,fr,verbose = verbose, code = code, k = k)
						}

						# Test de Wilcoxon pour la comparaison de position (compl\u00e8te Mood/Ansari)
						pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
						if (!is.na(pvals) && pvals <= alpha) {
						  ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
						               "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
						               "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
						               "\t==> Significant differences between groups (p = ", .format_pval(pvals), ").")
						  fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
						              "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
						              "b) Test de Wilcoxon-Mann-Whitney [wilcox.test()].\n\t",
						              "\t==> Diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
						               "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
						               "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
						               "\t==> No significant differences between groups (p = ", .format_pval(pvals), ").")
						  fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
						              "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
						              "b) Test de Wilcoxon-Mann-Whitney [wilcox.test()].\n\t",
						              "\t==> Pas de diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						}
  						#========================================================
  						#			NON-NORMAL		2 categories		Variances h\u00e9t\u00e9rog\u00e8nes + NON discret
  						#========================================================
					} else if ((check_variance_equal==FALSE) & (check_discret == FALSE)) {
						# Cas : non-normalit\u00e9 + variances h\u00e9t\u00e9rog\u00e8nes + 2 groupes
						# Solution : Wilcoxon-Mann-Whitney (robuste aux violations)

						# \u00c9tape 1 : \u00c9valuation de la fiabilit\u00e9 avec test KS
						pval_sidak <- 1-(1-alpha)^(1/length(unique(g)))
						if (ks_result <= pval_sidak) {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Warning! Groups do not have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
						               "Please check group distributions graphically.")
						  fr <- paste0("Estimation de la fiabilit\u00e9 du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites -\n\t",
						              "==> Attention ! Les groupes n'ont pas la m\u00eame distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
						              "(Veuillez v\u00e9rifier graphiquement les distributions des groupes.)")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						} else {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Groups have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test is expected to be reliable.")
						  fr <- paste0("Estimation de la fiabilit\u00e9 du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites -\n\t",
						              "==> Les groupes ont la m\u00eame distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						}

						# \u00c9tape 2 : Tests non param\u00e9triques
						pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
						if (pvals <= alpha) {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> Significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					} else {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> No significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Pas de diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					}
					} else {
						# Cas : non-normalit\u00e9 + variances HOMOG\u00c8NES + donn\u00e9es continues + 2 groupes
						# Solution : Wilcoxon-Mann-Whitney (test non-param\u00e9trique robuste)
						# R\u00e9f\u00e9rence: Fay & Proschan (2010). Wilcoxon-Mann-Whitney or t-test? Statistics Surveys, 4, 1-39.

						# \u00c9tape 1 : \u00c9valuation de la fiabilit\u00e9 avec test KS
						pval_sidak <- 1-(1-alpha)^(1/length(unique(g)))
						if (ks_result <= pval_sidak) {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Warning! Groups do not have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
						               "Please check group distributions graphically.")
						  fr <- paste0("Estimation de la fiabilit\u00e9 du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites -\n\t",
						              "==> Attention ! Les groupes n'ont pas la m\u00eame distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
						              "(Veuillez v\u00e9rifier graphiquement les distributions des groupes.)")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						} else {
						  ang <- paste0("Reliability assessment for Wilcoxon-Mann-Whitney test:\n\t",
						               "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
						               "==> Groups have the same distribution. p-value: ",  .format_pval(ks_result), "\n\t",
						               "	--> compared to Sidak corrected alpha ", .format_pval(pval_sidak), "\n\t",
						               "--> The Mann-Whitney-Wilcoxon test is expected to be reliable.")
						  fr <- paste0("Estimation de la fiabilit\u00e9 du test de Wilcoxon-Mann-Whitney :\n\t",
						              "Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites -\n\t",
						              "==> Les groupes ont la m\u00eame distribution. p-value : ", .format_pval(ks_result), "\n\t",
						              "	--> en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval_sidak), "\n\t",
						              "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "on")
						}

						# \u00c9tape 2 : Tests non param\u00e9triques
						pvals <- suppressWarnings(wilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]]))$p.value
						if (!is.na(pvals) && pvals <= alpha) {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> Significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					} else {
  						ang <- paste0("Non-parametric tests for comparing two groups.\n\t",
  						             "a) Bootstrap median differences analysis [pairwise.boot(mu='median') from {KefiR}]\n\t",
  						             "b) Wilcoxon-Mann-Whitney test [wilcox.test()].\n\t",
  						             "\t==> No significant differences between groups (p = ", .format_pval(pvals), ").")
  						fr <- paste0("Tests non param\u00e9triques de comparaison des deux groupes.\n\t",
  						            "a) Analyse des diff\u00e9rences de m\u00e9dianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]\n\t",
  						            "b) Test de Wilcoxon-Mann-Witney [wilcox.test()].\n\t",
  						            "\t==> Pas de diff\u00e9rences significatives entre les groupes (p = ", .format_pval(pvals), ").")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
  					}
					}
			} else { 											# > 2 categories
			  ###################################################
			  #			NON-NORMAL		k>2 categories
			  ###################################################
					##########################################
					#   DIAGNOSTIC PRAGMATIQUE
					#   D\u00e9cision du test non-param\u00e9trique
					##########################################
					# Calculer les m\u00e9triques n\u00e9cessaires

					# 0. Test de Fligner (n\u00e9cessaire pour le diagnostic)
					pvals3 <- fligner.test(x,g)$p.value
					if (is.na(pvals3)) {
					  .exit("Fligner-Killeen test [fligner.test()] - Error, return NA.",
						"Test de Fligner-Killeen [fligner.test()] - Erreur, retourne NA.",
						verbose=TRUE, return=TRUE)
					}

					# Calculer ks_result pour plus tard (test KS)
					temp <- pairwise(x,g,type="ks",silent=silent,boot=FALSE)$p.value
					ks_result <- min(unlist(temp),na.rm=TRUE)

					# 1. Asym\u00e9trie (skewness) GROUPE PAR GROUPE (prendre le max)
					# Utiliser agricolae::skewness qui est d\u00e9j\u00e0 dans les d\u00e9pendances
					skew_by_group <- as.vector(unlist(
					  by(x, g, function(z) abs(agricolae::skewness(z)))
					))
					skew_val <- max(skew_by_group, na.rm = TRUE)

					# 2. Ratio des variances max/min
					var_groups <- by(x, g, var, na.rm = TRUE)
					var_ratio <- max(var_groups, na.rm = TRUE) / min(var_groups, na.rm = TRUE)

					# 3. Outliers extr\u00eames (d\u00e9j\u00e0 calcul\u00e9s plus haut)
					has_extreme_outliers <- exists("max_outliers_extrem") && max_outliers_extrem > 0.02  # > 2%

					# 4. Distance de Cook
					cook_dist <- cooks.distance(myaov)
					max_cook <- max(cook_dist, na.rm = TRUE)
					seuil_cook <- 4 / length(x)
					cook_influential <- max_cook > seuil_cook

					# Affichage du diagnostic pragmatique
					k <- .vbse(
					  "PRAGMATIC DIAGNOSTIC to decide non-parametric test:",
					  "DIAGNOSTIC PRAGMATIQUE pour d\u00e9cider du test non param\u00e9trique :",
					  verbose = verbose, code = code, k = k, cpt = "on"
					)

					# G\u00e9n\u00e9ration code pour mode code=TRUE
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
					    "# Outliers extr\u00eames",
					    "library(rstatix)",
					    "outliers_by_group <- by(x, g, function(z) {",
					    "  result <- identify_outliers(data.frame(z))",
					    "  sum(result$is.extreme, na.rm = TRUE) / length(z)",
					    "})",
					    "max_outliers_extrem <- max(unlist(outliers_by_group))"
					  ))
					}

					# Afficher les m\u00e9triques
					k <- .vbse(
					  paste0("a) Skewness (asymmetry) - max |sk| across groups: ", round(skew_val, 2)),
					  paste0("a) Skewness (asym\u00e9trie) - max |sk| par groupe : ", round(skew_val, 2)),
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
					  paste0(" ==> INFLUENT (d\u00e9passe le seuil)")
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
					    paste0("e) Outliers extr\u00eames [rstatix::identify_outliers] - max par groupe : ",
					           round(max_outliers_extrem * 100, 1), "%"),
					    verbose = verbose, code = code, k = k, cpt = "off"
					  )
					}

					# Alpha corrig\u00e9 de \u0160id\u00e1k pour le test KS
					alpha_sidak <- 1 - (1 - alpha)^(1/length(unique(g)))

					# f) Test KS pour v\u00e9rifier la comparabilit\u00e9 des formes
					k <- .vbse(
					  paste0("f) Kolmogorov-Smirnov test [ks.test()] - min p-value: ", .format_pval(ks_result),
					         " (\u03b1_\u0160id\u00e1k = ", round(alpha_sidak, 3), ")"),
					  paste0("f) Test de Kolmogorov-Smirnov [ks.test()] - min p-value : ", .format_pval(ks_result),
					         " (\u03b1_\u0160id\u00e1k = ", round(alpha_sidak, 3), ")"),
					  verbose = verbose, code = code, k = k, cpt = "off"
					)

					# D\u00e9cision s\u00e9quentielle selon arbre de d\u00e9cision

					recommended_test <- NULL
					decision_reason <- NULL

					# \u00c9TAPE 1 : Outliers extr\u00eames > 5% ?
					if (has_extreme_outliers) {
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> EXTREME OUTLIERS detected (> 5%)\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> OUTLIERS EXTR\u00caMES d\u00e9tect\u00e9s (> 5%)\n\t",
					               "==> Test recommand\u00e9 : med1way() [ANOVA sur m\u00e9dianes]")
					  )
					# \u00c9TAPE 2 : Asym\u00e9trie forte (|skew| \u2265 2) ?
					} else if (skew_val >= 2) {
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> STRONG ASYMMETRY detected (|skew| = ", round(skew_val, 2), ")\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> ASYM\u00c9TRIE FORTE d\u00e9tect\u00e9e (|skew| = ", round(skew_val, 2), ")\n\t",
					               "==> Test recommand\u00e9 : med1way() [ANOVA sur m\u00e9dianes]")
					  )
					# \u00c9TAPE 3 : Formes des distributions similaires (KS test avec \u0160id\u00e1k) ?
					} else if (ks_result > alpha_sidak) {
					  # Formes similaires \u2192 Kruskal-Wallis (VALIDE m\u00eame si variances diff\u00e9rentes)
					  recommended_test <- "kruskal"
					  decision_reason <- list(
					    en = paste0("--> SIMILAR DISTRIBUTION SHAPES (KS p = ", .format_pval(ks_result), " > \u03b1_\u0160id\u00e1k)\n\t",
					               "==> Recommended test: Kruskal-Wallis [rank-based test]\n\t",
					               "    (valid even with heteroscedasticity if shapes are similar)"),
					    fr = paste0("--> FORMES DE DISTRIBUTIONS SIMILAIRES (KS p = ", .format_pval(ks_result), " > \u03b1_\u0160id\u00e1k)\n\t",
					               "==> Test recommand\u00e9 : Kruskal-Wallis [test sur rangs]\n\t",
					               "    (valide m\u00eame avec h\u00e9t\u00e9rosc\u00e9dasticit\u00e9 si formes similaires)")
					  )
					# \u00c9TAPE 4 : Formes diff\u00e9rentes \u2192 v\u00e9rifier h\u00e9t\u00e9rosc\u00e9dasticit\u00e9
					} else if (pvals3 <= 0.05) {
					  # Formes diff\u00e9rentes ET h\u00e9t\u00e9rosc\u00e9dasticit\u00e9 \u2192 t1way
					  recommended_test <- "t1way"
					  decision_reason <- list(
					    en = paste0("--> DIFFERENT SHAPES (KS p \u2264 \u03b1_\u0160id\u00e1k) + HETEROSCEDASTICITY (Fligner p = ", .format_pval(pvals3), ")\n\t",
					               "==> Recommended test: t1way() [ANOVA on trimmed means]"),
					    fr = paste0("--> FORMES DIFF\u00c9RENTES (KS p \u2264 \u03b1_\u0160id\u00e1k) + H\u00c9T\u00c9ROSC\u00c9DASTICIT\u00c9 (Fligner p = ", .format_pval(pvals3), ")\n\t",
					               "==> Test recommand\u00e9 : t1way() [ANOVA sur moyennes tronqu\u00e9es]")
					  )
					} else {
					  # Formes diff\u00e9rentes mais variances homog\u00e8nes \u2192 med1way
					  recommended_test <- "med1way"
					  decision_reason <- list(
					    en = paste0("--> DIFFERENT SHAPES (KS p \u2264 \u03b1_\u0160id\u00e1k) but homogeneous variances\n\t",
					               "==> Recommended test: med1way() [ANOVA on medians]"),
					    fr = paste0("--> FORMES DIFF\u00c9RENTES (KS p \u2264 \u03b1_\u0160id\u00e1k) mais variances homog\u00e8nes\n\t",
					               "==> Test recommand\u00e9 : med1way() [ANOVA sur m\u00e9dianes]")
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
					  # (remplace l'ancien arbre de d\u00e9cision)
					  ############

					  # Variable pour stocker le test choisi (bas\u00e9e sur diagnostic pragmatique)
					  chosen_test <- recommended_test

					  # Ex\u00e9cution selon le test recommand\u00e9
					  if (chosen_test == "med1way") {
						#################
						#	Forte asym\u00e9trie / queues lourdes \u2192 med1way()
						#################
						check_variance_equal <- FALSE

						# G\u00e9n\u00e9ration du code pour med1way
						if (isTRUE(code)) {
						  .code_anova(9, "ANOVA sur m\u00e9dianes (outliers/asym\u00e9trie)", c(
						    "library(WRS2)",
						    "med1way(x ~ g)"
						  ))
						}

						pvals3 <- tryCatch(med1way(x~g)$p.value, error = function(e) NA)
						if (is.na(pvals3)) {
						  k <- .vbse(
						    "Oneway ANOVA of medians [med1way()] - Failed, return NA.\n\t--> Fallback: Kruskal-Wallis test [kruskal.test()].",
						    "Analyse de la variance \u00e0 un facteur des m\u00e9dianes [med1way()] - \u00c9chec, retourne NA.\n\t--> Repli : test de Kruskal-Wallis [kruskal.test()].",
						    verbose = verbose, code = code, k = k
						  )
						  pvals3 <- kruskal.test(x, g)$p.value
						  chosen_test <- "kruskal"
						  if (!is.na(pvals3) && pvals3 <= alpha) {
						    k <- .vbse(
						      paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Significant differences between groups. p-value: ", .format_pval(pvals3), "."),
						      paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Diff\u00e9rences significatives entre les groupes. p-value : ", .format_pval(pvals3), "."),
						      verbose = verbose, code = code, k = k
						    )
						  } else if (!is.na(pvals3)) {
						    k <- .vbse(
						      paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Non-significant differences between groups. p-value: ", .format_pval(pvals3), "."),
						      paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Diff\u00e9rences non significatives entre les groupes. p-value : ", .format_pval(pvals3), "."),
						      verbose = verbose, code = code, k = k
						    )
						  }
						} else {
						  if (pvals3 <= alpha) {
						    ang <- paste0("Oneway ANOVA of medians [med1way()]:\n\t==> Significant differences between the median of groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
						    fr <- paste0("Analyse de la variance \u00e0 un facteur des m\u00e9dianes [med1way()] :\n\t==> Diff\u00e9rences significatives entre les m\u00e9dianes des groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  } else {
						    ang <- paste0("Oneway ANOVA of medians [med1way()]:\n\t==> Non-significant differences between the median of groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
						    fr <- paste0("Analyse de la variance \u00e0 un facteur des m\u00e9dianes [med1way()] :\n\t==> Diff\u00e9rences non significatives entre les m\u00e9dianes des groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  }
						}
					  } else if (chosen_test == "t1way") {
						#################
						#	Asym\u00e9trie mod\u00e9r\u00e9e + probl\u00e8mes variance/influence \u2192 t1way()
						#################
						check_variance_equal <- FALSE

						# G\u00e9n\u00e9ration du code pour t1way
						if (isTRUE(code)) {
						  .code_anova(9, "ANOVA sur moyennes tronqu\u00e9es (variances h\u00e9t\u00e9rog\u00e8nes)", c(
						    "library(WRS2)",
						    "t1way(x ~ g)"
						  ))
						}

						pvals3 <- t1way(x~g)$p.value
						if (is.na(pvals3)) {
						    ang <- paste0("Oneway ANOVA on trimmed means [t1way()] - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.")
							fr <- paste0("Analyse de la variance \u00e0 un facteur sur les moyennes tronqu\u00e9es [t1way()] - \u00c9chec, retourne NA. Le test de Kruskal-Wallis doit \u00eatre utilis\u00e9 pour l'interpr\u00e9tation.")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  if (pvals3 <= alpha) {
							ang <- paste0("Oneway ANOVA on trimmed means [t1way()]:\n\t==> Significant differences between the trimmed groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
							fr <- paste0("Analyse de la variance \u00e0 un facteur sur les moyennes tronqu\u00e9es [t1way()] :\n\t==> Diff\u00e9rences significatives entre les groupes tronqu\u00e9s.\n\t==> p-value : ", .format_pval(pvals3), ".")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
							if (!is.na(pvals) && pvals > alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronqu\u00e9es donnent des r\u00e9sultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "off")
							}
						  } else if (pvals3 > alpha) {
							ang <- paste0("One-way ANOVA on trimmed means [t1way()]:\n\t==> Non-significant differences between the trimmed groups.\n\t==> p-value: ",  .format_pval(pvals3), ".")
							fr <- paste0("Analyse de la variance \u00e0 un facteur sur les moyennes tronqu\u00e9es [t1way()] :\n\t==> Diff\u00e9rences non significatives entre les groupes tronqu\u00e9s.\n\t==> p-value : ", .format_pval(pvals3), ".")
							k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
							if (!is.na(pvals) && pvals <= alpha) {
								ang <- paste0("Warning! The Kruskal-Wallis test and ANOVA on trimmed means give contradictory results.")
								fr <- paste0("Attention ! Le test de Kruskal-Wallis et l'analyse de la variance sur les moyennes tronqu\u00e9es donnent des r\u00e9sultats contradictoires.")
								k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt = "off")
							}
						  }
						}
					  } else {
						#################
						#	Formes similaires \u2192 Kruskal-Wallis (test sur rangs)
						#################
						check_variance_equal <- TRUE  # Pour le routage post-hoc

						# G\u00e9n\u00e9ration du code pour Kruskal-Wallis
						if (isTRUE(code)) {
						  .code_anova(9, "Test de Kruskal-Wallis (formes similaires)", c(
						    "kruskal.test(x, g)"
						  ))
						}

						# Ex\u00e9cuter Kruskal-Wallis
						pvals3 <- kruskal.test(x, g)$p.value
						if (is.na(pvals3)) {
						  ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()] - Failed, return NA.")
						  fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] - \u00c9chec, retourne NA.")
						  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						} else {
						  if (pvals3 <= alpha) {
						    ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Significant differences between groups.\n\t==> p-value: ", .format_pval(pvals3), ".")
						    fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Diff\u00e9rences significatives entre les groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  } else {
						    ang <- paste0("Kruskal-Wallis rank sum test [kruskal.test()]:\n\t==> Non-significant differences between groups.\n\t==> p-value: ", .format_pval(pvals3), ".")
						    fr <- paste0("Test de Kruskal-Wallis [kruskal.test()] :\n\t==> Diff\u00e9rences non significatives entre les groupes.\n\t==> p-value : ", .format_pval(pvals3), ".")
						    k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
						  }
						}
					  } # fin choix test selon diagnostic pragmatique

					  # Mettre \u00e0 jour pvals avec la p-value du test choisi (pour le retour)
					  pvals <- pvals3

					  # Affichage du test KS apr\u00e8s med1way/t1way (avant les posthocs)
					  if (isTRUE(verbose)) {
					  if (ks_result <= pval) {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\t",
  						              "Warning! the groups do not have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t",
  						              "(comparing with Sidak corrected alpha ", .format_pval(pval), ")\n\t",
  						              "==> The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
  						              "\t(Please check graphically the groups distributions.)\n\t",
  						              "--> Adding Brunner-Munzel test as robust alternative.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites -\n\t",
  						             "Attention ! les groupes n'ont pas la m\u00eame distribution. min(p-value) : ", .format_pval(ks_result), "\n\t",
  						             "(en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval), ")\n\t",
  						             "==> Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
  						             "\t(Veuillez v\u00e9rifier graphiquement les distributions des groupes.)\n\t",
  						             "--> Ajout du test de Brunner-Munzel comme alternative robuste.")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
					  } else {
  						ang <- paste0("Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data.\n\t",
  						              "==> The groups have the same distribution. min(p-value): ", .format_pval(ks_result), "\n\t",
  						              "(comparing with Sidak corrected alpha ", .format_pval(pval), ")\n\t",
  						              "--> The Mann-Whitney-Wilcoxon test is reliable a priori.")
  						fr <- paste0("Test de Kolmogorov-Smirnov [ks.test()] sur les donn\u00e9es centr\u00e9es sur la m\u00e9diane et r\u00e9duites.\n\t",
  						             "==> Les groupes ont la m\u00eame distribution. min(p-value) : ", .format_pval(ks_result), "\n\t",
  						             "(en comparant avec l'alpha corrig\u00e9 de Sidak ", .format_pval(pval), ")\n\t",
  						             "--> Le test de Mann-Whitney-Wilcoxon est fiable a priori.")
  						k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
					  }
					  } # Fin affichage KS
					} else {
					  # return=FALSE & verbose=FALSE : fallback Kruskal-Wallis pour obtenir une p-value
					  pvals <- kruskal.test(x, g)$p.value
					} # test uniquement pour aller au bout du raisonnement
				} # fin d'ajustement \u00e0 la normalit\u00e9 malgr\u00e9 tout
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
			#fr <- paste0("Deux cat\u00e9gories.")
			#k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)

		  if (isTRUE(paired)) {
			###################################################
			#			NORMAL		2 categories	APPARI\u00c9
			###################################################
			# Pour donn\u00e9es appari\u00e9es, pas de test de variance - t-test appari\u00e9 directement
			lv <- levels(droplevels(factor(g)))
			pvals <- t.test(x[g==lv[1]], x[g==lv[2]], paired=TRUE)$p.value
			if (isTRUE(code)) {
				.code_anova(8, "Test de Student appari\u00e9", c(
					"lv <- levels(droplevels(factor(g)))",
					"t.test(x[g==lv[1]], x[g==lv[2]], paired = TRUE)"
				))
			}
			if (isTRUE(verbose)) {
				if (pvals <= alpha) {
				  ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Significant differences between conditions. p-value: ", .format_pval(pvals))
				  fr <- paste0("Test de Student appari\u00e9 [t.test(paired=TRUE)] :\n\t==> Diff\u00e9rences significatives entre les conditions. p-value : ", .format_pval(pvals))
				} else {
				  ang <- paste0("Paired Student test [t.test(paired=TRUE)] :\n\t==> Non-significant differences between conditions. p-value: ", .format_pval(pvals))
				  fr <- paste0("Test de Student appari\u00e9 [t.test(paired=TRUE)] :\n\t==> Diff\u00e9rences non significatives entre les conditions. p-value : ", .format_pval(pvals))
				}
				k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
			}

		  } else {
			###################################################
			#			NORMAL		2 categories	NON-APPARI\u00c9
			###################################################
			fm <- formula(x~g)
			pvals <- var.test(fm)$p.value
			if (pvals>alpha) {
			  ###################################################
			  #			NORMAL		2 categories	homogene variance
			  ###################################################
			  pvals <- t.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=TRUE)$p.value
			  if (isTRUE(code)) {
				  .code_anova(9, "Test de Student (variances homog\u00e8nes)", c(
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
					fr <- paste0("Test de Student [t.test()] :\n\t==> Diff\u00e9rences significatives entre les groupes. p-value : ", .format_pval(pvals))
					k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
				  } else {
					##################################
					#        NORMAL, Student unsignificant (var.equal=T)
					##################################
					ang <- paste0("Student test [t.test()] - Non-significant differences between groups. p-value: ",  .format_pval(pvals))
					fr <- paste0("Test de Student [t.test()] - Diff\u00e9rences non significatives entre les groupes. p-value : ", .format_pval(pvals))
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
				  .code_anova(9, "Test de Student-Welch (variances h\u00e9t\u00e9rog\u00e8nes)", c(
					  "levels_g <- levels(factor(g))",
					  "t.test(x[g==levels_g[1]], x[g==levels_g[2]], var.equal = FALSE, paired = FALSE)"
				  ))
			  }
			  ang <- paste0("Tests for comparing two normal groups with heterogeneous variances.\n\t",
						   "a) Analysis of mean differences by bootstrap [pairwise.boot(mu='mean') from {KefiR}]\n\t",
						   "b) Welch test [t.test(var.equal=FALSE)].\n\t",
						   if (pvals <= alpha) "\t==> Significant differences between groups (p = " else "\t==> No significant differences between groups (p = ",
						   .format_pval(pvals), ").")
			  fr <- paste0("Tests de comparaison de deux groupes normaux \u00e0 variances h\u00e9t\u00e9rog\u00e8nes.\n\t",
						  "a) Analyse des diff\u00e9rences de moyennes par bootstrap [pairwise.boot(mu='mean') de {KefiR}]\n\t",
						  "b) Test de Welch [t.test(var.equal=FALSE)].\n\t",
						  if (pvals <= alpha) "\t==> Diff\u00e9rences de moyennes significatives entre groupes (p = " else "\t==> Aucune diff\u00e9rence de moyennes significative entre groupes (p = ",
						  .format_pval(pvals), ").")
			  k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
			  # NOTE: Le message ci-dessus contient d\u00e9j\u00e0 la conclusion (significatif ou non)
			  # Pas besoin de message suppl\u00e9mentaire - \u00e9vite la redondance
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
			# R\u00e9utilisation de myaov cr\u00e9\u00e9 \u00e0 l'\u00e9tape 4
			pvals <- summary(myaov)[[1]][["Pr(>F)"]][1]
			if (isTRUE(code)) {
				.code_anova(9, "ANOVA \u00e0 un facteur (variances homog\u00e8nes)", c(
					"summary(myaov)"
				))
			}
			ang <- paste0("One-way ANOVA [aov()].\n\t",
			             "Result: p = ", .format_pval(pvals),
			             if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
			             if (pvals <= alpha) "==> Significant differences between groups." else "==> No significant differences between groups.")
			fr <- paste0("ANOVA \u00e0 un facteur [aov()].\n\t",
			            "R\u00e9sultat : p = ", .format_pval(pvals),
			            if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
			            if (pvals <= alpha) "==> Diff\u00e9rences de moyennes significatives entre groupes." else "==> Aucune diff\u00e9rence de moyennes significative entre groupes.")
			k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)
		  } else {												# Non-identical variances
			  ##################################
			  #        NORMAL, >2 categories var.equal=F => FANOVA.HETERO
			  ##################################
			pvals <- oneway.test(x~g,var.equal=FALSE)$p.value
			if (isTRUE(code)) {
				.code_anova(9, "Test de Welch (variances h\u00e9t\u00e9rog\u00e8nes)", c(
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
				fr <- paste0("Test F h\u00e9t\u00e9rosc\u00e9dastique de Welch [oneway.test()].\n\t",
				            "R\u00e9sultat : p = ", .format_pval(pvals),
				            if (pvals <= alpha) " (< alpha)." else " (>= alpha).", "\n\t",
				            if (pvals <= alpha) "==> Diff\u00e9rences de moyennes significatives entre groupes." else "==> Aucune diff\u00e9rence de moyennes significative entre groupes.")
				k <- .vbse(ang, fr, verbose = verbose, code = code, k = k)

				# Check fanova.hetero consistency
				if ((pvals <= alpha && pvals2 > alpha) || (pvals > alpha && pvals2 <= alpha)) {
					ang <- paste0("\tNote: fanova.hetero() suggests distribution differences (p = ", .format_pval(pvals2), ").")
					fr <- paste0("\tNote : fanova.hetero() sugg\u00e8re des diff\u00e9rences de distributions (p = ", .format_pval(pvals2), ").")
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
	  # (cf. Cahier des charges priorit\u00e9 2)
	  # NOTE: chosen_test indique quel test non-param\u00e9trique a \u00e9t\u00e9 utilis\u00e9 (med1way, t1way, kruskal, ou NULL)
	  # chosen_test est initialis\u00e9 \u00e0 NULL au d\u00e9but de la fonction et mis \u00e0 jour dans le diagnostic pragmatique
	  return(check <- list(x = x, g = g, check_normality = check_normality, check_variance_equal = check_variance_equal, k = k, global_pvalue = pvals, chosen_test = chosen_test))
 }
