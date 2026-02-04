library(vioplot) ; library(agricolae) ; library(onewaytests) ; library(KefiR) ; library(lawstat) ; library(WRS2)
library(FSA)
library(car)
library(dae)
library(stringr)
library(DescTools)
library(rstatix)
library(PMCMRplus)

# Fonction interne : retourne un message en fonction de la langue système
.msg <- function(en, fr) {
  # Obtenir la locale système pour LC_CTYPE
  sys_locale <- Sys.getlocale("LC_CTYPE")

  # Extraire la langue principale (avant le premier "_" ou ".")
  lang <- sub("_.*", "", sys_locale)
  lang <- sub("\\..*", "", lang)

  # Retourner le message en fonction de la langue détectée
  if (grepl("French", lang)) {
    return(fr)
  } else {
    return(en)
  }
}
################################
#
################################
################################
################################
################################
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
  # ═══════════════════════════════════════════════════════════════════════════

  ## 7.1 : Suffixer idF par Temps et F
  for (i in unique(data$idF)) {
    for (j in unique(data$F)) {
      mask <- data$idF == i & data$F == j
      if (sum(mask) > 0) {
        data$idF[mask] <- paste0(data$idF[mask], "_", 1:sum(mask))
      }
    }
  }

  ## 7.2 : Suffixer idG par Temps et G
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



# Exemple d'utilisation
data <- simul(20)
data <- simul(60)
data <- simul(60) ; vioplot(data$A~data$H) ; .normality(data$A,data$H)
data <- simul(100)
data <- simul(200)
################################
#
################################
#* Note perso : changer l'adresse par défaut pour charger les scripts dans le dossier où est .testeur.R
source("C:/Users/masse/Desktop/KefiR/KefiR/R/m.test_temp17.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/.multi_factor_analysis02.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/discret.test.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/.normality.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/.variance.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/.boots.R")
source("C:/Users/masse/Desktop/KefiR/KefiR/R/.exit.R")
###############################
# Normalisation des formules
###############################
head(data)
.formulator(A~F)
.formulator(A~F+Error(G))
.formulator(cbind(A,B,C,D)~F+Error(G))
.formulator(data$B~data$F)
###############################
# Un seul facteur (m.test classique) - 2 niveaux
###############################
m.test(A~F, data=data, debug=TRUE)
m.test(A~F, data=data, debug=FALSE)
m.test(B~F, data=data)
m.test(B~F, data=data, debug=TRUE)
m.test(C~F, data=data)
m.test(C~F, data=data, debug=TRUE)
m.test(D~F, data=data)
m.test(E~F, data=data)
m.test(I~F, data=data)
m.test(J~F, data=data)
m.test(round(data$A/2),data$F)
m.test(data$B~data$F,debug=FALSE)
m.test(data$B~data$F,debug=TRUE)
.formulator(data$B~data$F,debug=TRUE) # Problèmes actuels liés à formulator()
m.test(C~data$F, data=data)
m.test(data$D~F, data=data)
m.test(data$E,data$F)
###############################
# Un seul facteur (m.test classique) - 3 niveaux
###############################
m.test(A~G, data=data)
m.test(A~G, data=data, debug=TRUE)
m.test(B~G, data=data)
m.test(C~G, data=data, debug=FALSE)
m.test(D~G, data=data)
m.test(E~G, data=data)
###############################
# Un seul facteur (m.test classique) - n niveaux
###############################
m.test(A~H, data=data)
m.test(A~H, data=data, debug=TRUE)
m.test(B~H, data=data)
m.test(C~H, data=data, debug=FALSE)
m.test(D~H, data=data)
m.test(E~H, data=data)
###############################
# Un seul facteur (m.test classique) avec control
###############################
m.test(A~F, data=data, control="catF1")
m.test(A~F, data=data, control="catF2")
m.test(B~F, data=data, control="catF1")
m.test(B~F, data=data, control="catF2")
m.test(C~F, data=data, control="catF1")
m.test(C~F, data=data, control="catF2")
m.test(D~F, data=data, control="catF1")
m.test(D~F, data=data, control="catF2")
m.test(E~F, data=data, control="catF1")
m.test(E~F, data=data, control="catF2")
# Forcer les bugs
m.test(A~G, data=data, control="catF1")->b
m.test(A~B, data=data, control="catF1")->b # Bug à finir.



# NE MARCHE PAS .boots(x, g, ctrl = TRUE, type = "mean", var.equal = TRUE, conf = conf,  : objet 'ind_control' introuvable
m.test(A~F, data=data, control="catF1")
m.test(A~F, data=data, control="catF2",boot=FALSE)
m.test(B~F, data=data, control="catF1")
m.test(C~F, data=data, control="catF1")
m.test(D~F, data=data, control="catF2")
m.test(E~F, data=data, control="catF1")
m.test(A~G, data=data, control="catG1")->b
###############################
# Données appariées
###############################
m.test(A~F, data=data,id="idF")
m.test(A~F, data=data,id="id") # Message  d'erreur exprès (test)
m.test(A~F, data=data,id="idG") # Message  d'erreur exprès (test)
m.test(B~F, data=data,id="idF")
m.test(C~F, data=data,id="idF")
m.test(D~F, data=data,id="idF")
m.test(E~F, data=data,id="idF")
m.test(A~G, data=data,id="idG")
# Vérifier que les données appariées se combinent bien au contrôle sur 2 catégories.
# Exemple où l'identifiant unique des valeurs ne se retrouve pas dans chaque catégorie.
m.test(A~F, data=data, control="catF1",id="id")



###############################
# Multifacteur
###############################
m.test(A~G, data=data)
m.test(A~G, data=data, id="idG")
m.test(A~G, data=data, id="idG", debug=TRUE)
m.test(A~F+G, data=data)
m.test(A~F+G, data=data, id="idG")
m.test(A~F+G, data=data, debug=TRUE)
m.test(A~F+G, data=data, id="idG")
m.test(A~F+G, data=data, id="idG", debug=TRUE)

###############################
# DIVERS
###############################
m.test(log(A)~G, data=data, control="catG3")->b
m.test(B~G, data=data)->b
m.test(B~G, data=data, control="catG3")->b
m.test(log(B)~G, data=data, control="catG3")->b
m.test(C~G, data=data)->b
m.test(C~G, data=data, control="catG3")->b
m.test(log(C)~G, data=data, control="catG3")->b
m.test(D~G, data=data)->b
m.test(D~G, data=data, control="catG3")->b
m.test(log(D)~G, data=data, control="catG3")->b
m.test(E~G, data=data)->b
m.test(E~G, data=data, control="catG3")->b
m.test(log(E)~G, data=data, control="catG3")->b
m.test(A~G, data=data, id="id", paired=TRUE, debug=F)->b
m.test(A~G, data=data, paired=TRUE, debug=TRUE)->b # BUG BLOC SUR BOOTSTRAP de neminye
m.test(D~G, data=data, id="id", paired=TRUE)->b # BUG BLOC SUR BOOTSTRAP de neminye
pairwise(data$A,data$G,type="neminye")
m.test(I(log(A))~F, data=data, debug=TRUE) # A n'est pas mis en log
m.test(I(log(A))~F+B, data=data, debug=TRUE)
m.test(I(log(A))~F+B*C, data=data, debug=TRUE)
m.test(I(log(A))~F+B*C, data=data, paired=TRUE, id="id", debug=TRUE) # id bug
m.test(I(log(A))~F+B*C+Error(id), data=data, paired=TRUE, id="id", debug=TRUE) # ERREUR
m.test(A~F, id = "id", wt="Temps", data=data, debug=TRUE)
m.test(A~F+Error(id), id = "id", wt="Temps", data=data, debug=TRUE)  # ERREUR
m.test(A~F+Error(G), id = "id", wt="Temps", data=data, debug=TRUE)  # ERREUR ==> lmer(A ~ F + Temps + (1|G) + (id/Temps), data = data)
m.test(cbind(A,B)~F+B, data=data)
m.test(data$A,data$F, debug=TRUE)->g
m.test(A~F,data, debug=TRUE)
m.test(A~F,data=data, debug=TRUE)
m.test(I(log(A))~F,data=data, debug=TRUE)
A_ext <- log(data$A)
m.test(A_ext~F,data=data, debug=TRUE)
F_ext <- data$F
m.test(I(log(A))~F+F_ext,data=data, debug=TRUE)
m.test(cbind(A,B,C)~F+F_ext,data=data, debug=TRUE)


###############################
# MANOVA
###############################
m.test(cbind(A,B)~F+G, data=data)
# NOTE: Les tests doivent passer par m.test() qui appelle .manova_analysis()
# Charger les fonctions
source("C:/Users/masse/Claude/m.test_temp18.R")
source("C:/Users/masse/Claude/.manova_analysis.R")
source("C:/Users/masse/Claude/.one_factor_analysis.R")
source("C:/Users/masse/Claude/.multi_factor_analysis02.R")

# Générer des données de test
data <- simul(60)
# Test 1: MANOVA classique avec 2 groupes et 2 variables dépendantes
m.test(cbind(A, B) ~ F, data=data, verbose=TRUE, debug=FALSE)
# Test 2: MANOVA classique avec 2 groupes et 3 variables dépendantes
m.test(cbind(A, B, C) ~ F, data=data, verbose=TRUE, debug=FALSE)
# Test 3: MANOVA avec 2 facteurs (interaction)
m.test(cbind(A, B, C) ~ F + G, data=data, verbose=TRUE, debug=FALSE)
# Test 4: MANOVA avec 3 variables dépendantes et 3 groupes
m.test(cbind(A, B, C) ~ G, data=data, verbose=TRUE, debug=FALSE)
# Test 5: MANOVA avec plusieurs groupes (H déséquilibré)
m.test(cbind(A, B, C) ~ H, data=data, verbose=TRUE, debug=FALSE)
# Test 6: MANOVA avec variables non normales (I et J) - devrait déclencher robuste
m.test(cbind(I, J, D) ~ F, data=data, verbose=TRUE, debug=FALSE)
# Test 7: MANOVA robuste avec 2 variables non normales
m.test(cbind(I, J) ~ G, data=data, verbose=TRUE, debug=FALSE)



##################################################
#   Testing à finir
##################################################
# ATTENTION, catego() doit se faire refaire dans la désignation des lettres
# Beaucoup de bugs...
# Wilcoxon est mis en situation d'erreur avec <NA>
data(TitanicSurvival)

m.test(age~survived+sex+passengerClass,data=TitanicSurvival,debug=TRUE)->b
m.test(age~survived+sex+passengerClass,data=TitanicSurvival,debug=FALSE)->b
# continuer ici car g_cat introuvable...
