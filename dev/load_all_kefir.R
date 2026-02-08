# Script pour charger tous les fichiers KefiR dans l'ordre correct
# Session 4 - Priorité 4

# Fonction helper pour charger packages avec gestion d'erreurs
load_package <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    warning(paste0("Package '", pkg_name, "' not available. Some features may be limited."))
    return(FALSE)
  }
  library(pkg_name, character.only = TRUE)
  return(TRUE)
}

# Charger packages (avec gestion erreurs)
load_package("rstatix")
load_package("lmerTest")
load_package("lmtest")
load_package("agricolae")
load_package("DescTools")
load_package("stringr")
load_package("PMCMRplus")
load_package("WRS2")
load_package("car")
load_package("onewaytests")
load_package("lmPerm")        # ANOVA par permutation
load_package("rcompanion")    # Scheirer-Ray-Hare test

# Supprimer le warning tidyselect de rstatix (externe, hors de notre contrôle)
suppressWarnings({
  options(lifecycle_verbosity = "quiet")
})

# Détecter OS et ajuster chemin
if (.Platform$OS.type == "windows") {
  setwd("C:/Users/masse/Desktop/KefiR/KefiR/R")
} else {
  # Linux/WSL
  setwd("/mnt/c/Users/masse/Desktop/KefiR/KefiR/R")
}

cat("\n=== CHARGEMENT DE TOUS LES FICHIERS KEFIR ===\n\n")

# Liste des fichiers à charger dans l'ordre de dépendance
fichiers <- c(
  "sys_msg.R",
  "sys_dbg.R",
  "sys_vbse.R",
  "sys_format_pval.R",
  "sys_drop_error_term.R",
  "sys_strip_data_dollar_safe.R",
  "sys_normalize_formula_dollar.R",
  "sys_normalize_from_formula.R",
  "sys_formulator.R",
  "sys_formulator_safe.R",
  "sys_detect_model_type.R",
  "sys_exit.R",
  "sys_filter_small_groups.R",       # Filtrage groupes trop petits
  "sys_select_ss_type.R",            # NOUVEAU: Sélection intelligente Type I/II/III
  "sys_analyze_ancova_structure.R",  # NOUVEAU: Détection ANCOVA vs pentes hétérogènes (après dépendances)
  "sys_diagnostic_influence.R",      # NOUVEAU: Diagnostics influence (Cook/Leverage/DFBETAS)
  "sys_ancova_analysis.R",
  "sys_align_pairs.R",
  "sys_auto_preprocess_g.R",
  "sys_control_independence.R",
  "sys_normality.R","catego.R",
  "sys_variance.R",
  "sys_boots.R",
    "pairwise.R", "pairwise.boot.R","sys_posthoc_mixed_model.R",
  "00_jb.norm.test.R",
  "discret.test.R",                # PRIORITÉ 7: Améliorations valreg()
  "sys_one_factor_analysis.R",
  "sys_manova_analysis.R",
  "print.manova_result.R",   # Méthode print() pour MANOVA (évite affichage verbeux)
  "sys_multi_factor_analysis.R",
  "sys_mixed_model_analysis.R",  # PRIORITÉ 6: Modèles mixtes
  "sys_posthoc.R",
  "sys_posthoc_MANOVA.R",
  "sys_posthoc_ANCOVA.R",
  "sys_plot_with_letters.R",    # PRIORITÉ 5: Graphics avec lettres
  "m.test.R",
  "valreg.R","bootreg.R","confidence.R",
  "sys_testeur_m.test.R",        # Suite de tests exhaustifs
  "sys_testeur_one_factor.R",    # Tests analyse à un facteur
  "sys_testeur_ANCOVA.R",         # Tests ANCOVA
  "sys_pairwise_medpb2.R"
)

errors <- c()
loaded <- c()

for (f in fichiers) {
  if (file.exists(f)) {
    result <- try(source(f, encoding = "UTF-8"), silent = TRUE)
    if (inherits(result, "try-error")) {
      errors <- c(errors, f)
      cat(sprintf("✗ %s (ERREUR)\n", f))
    } else {
      loaded <- c(loaded, f)
      cat(sprintf("✓ %s\n", f))
    }
  } else {
    cat(sprintf("- %s (non trouvé)\n", f))
  }
}

cat(sprintf("\n=== RÉSUMÉ ===\n"))
cat(sprintf("Chargés:  %d/%d\n", length(loaded), length(fichiers)))
cat(sprintf("Erreurs:  %d\n", length(errors)))

if (length(errors) > 0) {
  cat("\nFichiers avec erreurs:\n")
  for (e in errors) {
    cat(sprintf("  - %s\n", e))
  }
}



cat("\n")

neutralise <- 0
if (neutralise !=0) {

#' Génération de données aléatoires pour tests ANOVA mixtes
#'
#' @param n Nombre de sujets par cellule du design 2x3 (doit être divisible par 4)
#' @param k_H Nombre de modalités du facteur H déséquilibré (défaut: aléatoire entre 4 et 20)
#' @param desequilibre_H Degré de déséquilibre pour H :
#'   \itemize{
#'     \item "faible" : ratio max/min < 2
#'     \item "modere" : ratio max/min entre 2 et 5
#'     \item "fort" : ratio max/min > 5
#'     \item "aleatoire" : poids uniformément aléatoires (défaut)
#'   }
#' @param seed Graine aléatoire pour reproductibilité (défaut: NULL)
#'
#' @return Un data.frame avec :
#'   - Facteurs entre-sujets : F (2 niveaux), G (3 niveaux), H (k niveaux déséquilibrés)
#'   - Facteur intra-sujet : Temps (4 mesures répétées)
#'   - Variables numériques : A (dépend de F), B (dépend de G), C, D, E
#'   - Identifiants : id, idF (appariement sur F), idG (appariement sur G)
#'
#' @examples
#' # Cas 1 : Déséquilibre fort avec 5 groupes
#' data1 <- simul(n = 12, k_H = 5, desequilibre_H = "fort")
#' table(data1$H)  # Ratios typiques : 10:1 ou plus
#'
#' # Cas 2 : Déséquilibre modéré avec 10 groupes
#' data2 <- simul(n = 20, k_H = 10, desequilibre_H = "modere")
#' table(data2$H)
#'
#' # Cas 3 : Reproduction exacte (seed fixé)
#' data3 <- simul(n = 8, k_H = 7, desequilibre_H = "faible", seed = 123)
#'
#' @export
simul <- function(n = 12,
                  k_H = NULL,
                  desequilibre_H = c("aleatoire", "faible", "modere", "fort"),
                  seed = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1 : VALIDATION ET INITIALISATION
  # ═══════════════════════════════════════════════════════════════════════════
  if (is.null(seed)) {
    seed <- round(runif(1)*100)
  }
  set.seed(seed)

  if (n %% 4 != 0) {
    stop("n doit être divisible par 4 (design 2x3 équilibré)")
  }

  desequilibre_H <- match.arg(desequilibre_H)

  # Nombre de modalités H : aléatoire si non spécifié
  if (is.null(k_H)) {
    k_H <- sample(4:20, 1)
  } else {
    if (k_H < 2) stop("k_H doit être >= 2")
  }

  m <- n / 4  # Sujets par cellule F x G

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2 : DESIGN ENTRE-SUJETS (F x G)
  # ═══════════════════════════════════════════════════════════════════════════

  levels_F <- c("catF1", "catF2")
  levels_G <- c("catG1", "catG2", "catG3")

  design_sujet <- expand.grid(F = levels_F, G = levels_G, KEEP.OUT.ATTRS = FALSE)
  design_sujet <- design_sujet[rep(seq_len(nrow(design_sujet)), each = m), ]
  rownames(design_sujet) <- NULL
  design_sujet$id <- paste0("S", seq_len(nrow(design_sujet)))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3 : APPARIEMENTS (idF et idG)
  # ═══════════════════════════════════════════════════════════════════════════

  ## 3.1 : idF - Appariement entre niveaux de F (par G), désordonné
  design_sujet$idF <- NA_character_
  for (g in levels_G) {
    idx1 <- which(design_sujet$G == g & design_sujet$F == "catF1")
    idx2 <- which(design_sujet$G == g & design_sujet$F == "catF2")
    idx1 <- sample(idx1)  # Mélange aléatoire
    idx2 <- sample(idx2)
    for (i in seq_len(m)) {
      tag <- paste0("PF_", g, "_", i)
      design_sujet$idF[idx1[i]] <- tag
      design_sujet$idF[idx2[i]] <- tag
    }
  }

  ## 3.2 : idG - Appariement entre niveaux de G (par F), désordonné
  design_sujet$idG <- NA_character_
  for (f in levels_F) {
    i1 <- which(design_sujet$F == f & design_sujet$G == "catG1")
    i2 <- which(design_sujet$F == f & design_sujet$G == "catG2")
    i3 <- which(design_sujet$F == f & design_sujet$G == "catG3")
    i1 <- sample(i1); i2 <- sample(i2); i3 <- sample(i3)
    for (i in seq_len(m)) {
      tag <- paste0("PG_", f, "_", i)
      design_sujet$idG[i1[i]] <- tag
      design_sujet$idG[i2[i]] <- tag
      design_sujet$idG[i3[i]] <- tag
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4 : MESURES RÉPÉTÉES (Temps 1..4)
  # ═══════════════════════════════════════════════════════════════════════════

  data <- design_sujet[rep(seq_len(nrow(design_sujet)), each = 4), ]
  data$Temps <- rep(1:4, times = nrow(design_sujet))
  N <- nrow(data)  # = 6 * n (24 * m)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5 : VARIABLES NUMÉRIQUES (A, B, C, D, E)
  # ═══════════════════════════════════════════════════════════════════════════

  ## 5.1 : A dépend de F
  mu_A <- c("catF1" = sample(0:10, 1), "catF2" = sample(0:10, 1))
  data$A <- NA_real_
  data$A[data$F == "catF1"] <- rnorm(sum(data$F == "catF1"), mean = mu_A["catF1"], sd = 1)
  data$A[data$F == "catF2"] <- rnorm(sum(data$F == "catF2"), mean = mu_A["catF2"], sd = 1)

  ## 5.2 : B dépend de G
  mu_B <- c("catG1" = sample(0:6, 1),
            "catG2" = sample(0:6, 1),
            "catG3" = sample(0:6, 1))
  data$B <- NA_real_
  data$B[data$G == "catG1"] <- rnorm(sum(data$G == "catG1"), mean = mu_B["catG1"], sd = 2.5)
  data$B[data$G == "catG2"] <- rnorm(sum(data$G == "catG2"), mean = mu_B["catG2"], sd = 1)
  data$B[data$G == "catG3"] <- rnorm(sum(data$G == "catG3"), mean = mu_B["catG3"], sd = 2.5)

  ## 5.3 : C et E corrélées à A et B (r ≈ 0.3)
  noise_sd <- sqrt(1 - 0.3^2 - 0.3^2)  # ≈ 0.8246
  data$C <- 0.3 * data$A + 0.3 * data$B + rnorm(N, 0, noise_sd)
  data$E <- 0.3 * data$A + 0.3 * data$B + rnorm(N, 0, noise_sd)

  ## 5.4 : D aléatoire (Poisson ou Uniforme)
  if (runif(1) < 0.5) {
    data$D <- rpois(N, lambda = sample(1:10, 1))
  } else {
    data$D <- runif(N, 0, 10)
  }

  ## 5.5 : baseline - Covariable pour ANCOVA (corrélée avec A)
  data$baseline <- 70 + 0.5 * data$A + rnorm(N, 0, 10)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5bis : VARIABLES NON NORMALES I et J (NOUVEAU)
  # ═══════════════════════════════════════════════════════════════════════════

  ## 5bis.1 : Variable I - Distributions non normales selon F
  ## -----------------------------------------------------------------------
  ## Chaque niveau de F a une distribution différente aléatoire

  # Liste des distributions possibles
  distrib_types <- c("uniforme", "poisson", "exponentielle", "gamma_asymetrique",
                     "bimodale", "lognormale", "cauchy_tronque")

  # Sélectionner aléatoirement une distribution pour chaque niveau de F
  distrib_F <- sample(distrib_types, length(levels_F), replace = TRUE)
  names(distrib_F) <- levels_F

  data$I <- NA_real_

  for (f in levels_F) {
    idx <- data$F == f
    n_obs <- sum(idx)

    data$I[idx] <- switch(
      distrib_F[f],

      # Uniforme [0, 10]
      "uniforme" = {
        runif(n_obs, min = 0, max = 10)
      },

      # Poisson (lambda aléatoire)
      "poisson" = {
        rpois(n_obs, lambda = sample(2:8, 1))
      },

      # Exponentielle (asymétrie forte à droite)
      "exponentielle" = {
        rexp(n_obs, rate = 1/sample(2:5, 1))
      },

      # Gamma asymétrique (skewness positif)
      "gamma_asymetrique" = {
        rgamma(n_obs, shape = sample(c(1.5, 2, 2.5), 1), rate = 0.5)
      },

      # Bimodale (mélange de deux normales)
      "bimodale" = {
        n1 <- rbinom(1, n_obs, 0.5)
        n2 <- n_obs - n1
        c(rnorm(n1, mean = 2, sd = 1),
          rnorm(n2, mean = 7, sd = 1))
      },

      # Log-normale (asymétrie forte à droite)
      "lognormale" = {
        rlnorm(n_obs, meanlog = 1, sdlog = 0.5)
      },

      # Cauchy tronquée (queues lourdes + outliers)
      "cauchy_tronque" = {
        vals <- rcauchy(n_obs, location = 5, scale = 2)
        # Tronquer pour éviter valeurs extrêmes
        pmin(pmax(vals, 0), 15)
      }
    )
  }

  ## 5bis.2 : Variable J - Distributions non normales selon G
  ## -----------------------------------------------------------------------
  ## Chaque niveau de G a une distribution différente aléatoire

  # Sélectionner aléatoirement une distribution pour chaque niveau de G
  distrib_G <- sample(distrib_types, length(levels_G), replace = TRUE)
  names(distrib_G) <- levels_G

  data$J <- NA_real_

  for (g in levels_G) {
    idx <- data$G == g
    n_obs <- sum(idx)

    data$J[idx] <- switch(
      distrib_G[g],

      # Uniforme [0, 10]
      "uniforme" = {
        runif(n_obs, min = 0, max = 10)
      },

      # Poisson (lambda aléatoire)
      "poisson" = {
        rpois(n_obs, lambda = sample(2:8, 1))
      },

      # Exponentielle (asymétrie forte à droite)
      "exponentielle" = {
        rexp(n_obs, rate = 1/sample(2:5, 1))
      },

      # Gamma asymétrique (skewness positif)
      "gamma_asymetrique" = {
        rgamma(n_obs, shape = sample(c(1.5, 2, 2.5), 1), rate = 0.5)
      },

      # Bimodale (mélange de deux normales)
      "bimodale" = {
        n1 <- rbinom(1, n_obs, 0.5)
        n2 <- n_obs - n1
        c(rnorm(n1, mean = 2, sd = 1),
          rnorm(n2, mean = 7, sd = 1))
      },

      # Log-normale (asymétrie forte à droite)
      "lognormale" = {
        rlnorm(n_obs, meanlog = 1, sdlog = 0.5)
      },

      # Cauchy tronquée (queues lourdes + outliers)
      "cauchy_tronque" = {
        vals <- rcauchy(n_obs, location = 5, scale = 2)
        # Tronquer pour éviter valeurs extrêmes
        pmin(pmax(vals, 0), 15)
      }
    )
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6 : FACTEUR H DÉSÉQUILIBRÉ (VERSION AMÉLIORÉE)
  # ═══════════════════════════════════════════════════════════════════════════

  ## 6.1 : Génération des probabilités selon le degré de déséquilibre
  p_H <- switch(
    desequilibre_H,

    # -------------------------------------------------------------------------
    # Cas 1 : Aléatoire (comme version originale)
    # -------------------------------------------------------------------------
    "aleatoire" = {
      probs <- runif(k_H)
      probs / sum(probs)
    },

    # -------------------------------------------------------------------------
    # Cas 2 : Faible (ratio max/min < 2)
    # -------------------------------------------------------------------------
    "faible" = {
      # Poids proches (écart-type faible)
      probs <- rnorm(k_H, mean = 1, sd = 0.2)
      probs <- pmax(probs, 0.5)  # Éviter valeurs trop faibles
      probs / sum(probs)
    },

    # -------------------------------------------------------------------------
    # Cas 3 : Modéré (ratio max/min entre 2 et 5)
    # -------------------------------------------------------------------------
    "modere" = {
      # Distribution log-normale (asymétrie modérée)
      probs <- rlnorm(k_H, meanlog = 0, sdlog = 0.5)
      probs / sum(probs)
    },

    # -------------------------------------------------------------------------
    # Cas 4 : Fort (ratio max/min > 5)
    # -------------------------------------------------------------------------
    "fort" = {
      # 1-2 groupes majoritaires + plusieurs petits groupes
      probs <- numeric(k_H)
      n_major <- sample(1:2, 1)  # Nombre de groupes majoritaires

      # Groupes majoritaires : 50-70% des effectifs
      major_weight <- runif(1, 0.5, 0.7)
      probs[1:n_major] <- runif(n_major, 0.3, 0.7)
      probs[1:n_major] <- probs[1:n_major] / sum(probs[1:n_major]) * major_weight

      # Petits groupes : 30-50% restants (très déséquilibrés)
      if (k_H > n_major) {
        remaining <- 1 - major_weight
        small_weights <- rexp(k_H - n_major, rate = 3)  # Distribution exponentielle
        probs[(n_major + 1):k_H] <- small_weights / sum(small_weights) * remaining
      }

      probs
    }
  )

  ## 6.2 : Affectation du facteur H
  data$H <- factor(
    sample(paste0("catH", seq_len(k_H)),
           size = N,
           replace = TRUE,
           prob = p_H)
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7 : SUFFIXAGE DES IDENTIFIANTS (par Temps)
  # ═══════════================================================================

  ## 7.0 : Créer versions "paired-ready" AVANT suffixage
  ## Nécessaire pour m.test(..., paired=TRUE, id=data$idF_paired)
  ## Ces IDs permettent appariement 1:1 entre niveaux de F ou G
  data$idF_paired <- data$idF  # Conservation ID original pour tests appariés
  data$idG_paired <- data$idG  # Conservation ID original pour tests appariés

  ## 7.1 : Suffixer idF par Temps et F (pour mesures répétées)
  for (i in unique(data$idF)) {
    for (j in unique(data$F)) {
      mask <- data$idF == i & data$F == j
      if (sum(mask) > 0) {
        data$idF[mask] <- paste0(data$idF[mask], "_", 1:sum(mask))
      }
    }
  }

  ## 7.2 : Suffixer idG par Temps et G (pour mesures répétées)
  for (i in unique(data$idG)) {
    for (j in unique(data$G)) {
      mask <- data$idG == i & data$G == j
      if (sum(mask) > 0) {
        data$idG[mask] <- paste0(data$idG[mask], "_", 1:sum(mask))
      }
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8 : MÉLANGE FINAL ET RETOUR
  # ═══════════════════════════════════════════════════════════════════════════

  data <- data[sample.int(N, N), ]
  rownames(data) <- NULL

  # Ajouter attributs pour diagnostics
  attr(data, "design") <- list(
    n_sujets = nrow(design_sujet),
    k_H = k_H,
    desequilibre_H = desequilibre_H,
    effectifs_H = table(data$H),
    ratio_H = max(table(data$H)) / min(table(data$H)),
    seed = seed,
    # Informations sur les distributions non normales
    distributions_I_par_F = distrib_F,
    distributions_J_par_G = distrib_G
  )

  return(data)
}

###########################
# TESTS AUTOMATIQUES
###########################
# Pour exécuter les tests automatiques, décommentez les lignes suivantes:
# .testeur_one_factor()
# .testeur_ANCOVA()
dt <- simul(20)
m.test(B~F+G, data=dt)
m.test(B~F*G, data=dt)
m.test(B~G, data=dt)
m.test(B~G+A, data=dt)

n<-10 ; x <- c(rnorm(n,3,2),rnorm(n,4,2),rnorm(n,2,2))
g <- rep(c("A","B","C"),each=n)
m.test(x~g)
#
n<-10 ; x <- c(rnorm(n,3,2),rnorm(n,6,7),rnorm(n,2,2))
g <- rep(c("A","B","C"),each=n)
m.test(x~g)
#
n<-10 ; x <- c(runif(n,0,5),rnorm(n,2,8),rnorm(n,-2,4))
g <- rep(c("A","B","C"),each=n)
m.test(x~g)
# Distributions différentes
n<-50 ; x <- c(runif(n,0,5),rnorm(n,2,2), rpois(n,1)*runif(n,0,5))
g <- rep(c("A","B","C"),each=n)
m.test(x~g)


# Distributions différentes
n<-50 ; x <- c(runif(n,0,5),rnorm(n,2,2), rgamma(n, shape = 4, scale = 2)/3 )  #)
g <- rep(c("A","B","C"),each=n)
m.test(x~g)->res

#
m.test(D~H, data=dt, verbose=TRUE) -> res
#pairwise.boot(dt$D,dt$H, mu="median")
catego(pairwise.boot(dt$D,dt$H, mu="mean"))$groups
res$groups


}

