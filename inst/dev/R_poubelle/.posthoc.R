#' Réaliser différents tests post-hoc en fonction du contexte statistique
#'
#' Cette fonction \code{.posthoc} réalise un ensemble de tests post-hoc
#' (tests de Student, Wilcoxon, SNK, Dunn, etc.) en fonction :
#' \itemize{
#'   \item du nombre de groupes à comparer (\code{number}) ;
#'   \item de la normalité présumée des données (\code{normal}) ;
#'   \item de l'égalité présumée des variances (\code{var.equal}) ;
#'   \item du caractère apparié (\code{paired}) ;
#'   \item de la présence éventuelle d'un groupe de contrôle (\code{control}).
#' }
#'
#' Lorsque plusieurs groupes sont comparés, la fonction applique les procédures
#' post-hoc adaptées (par exemple, \code{\link[agricolae]{SNK.test}} pour SNK,
#' \code{\link[FSA]{dunnTest}} pour Dunn, \code{\link[DescTools]{DunnettTest}}
#' pour Dunnett, etc.). En cas d'options supplémentaires (\code{boot}, \code{conf},
#' \code{iter}), la fonction peut également réaliser des tests bootstrap via
#' \code{.boots()} et des analyses de paires via \code{pairwise()} (fonctionnalités
#' internes ou fournies par des paquets tiers).
#'
#' @param x Numérique. Le vecteur de données à analyser.
#' @param g Facteur ou vecteur catégoriel. Indique à quel groupe chaque valeur
#'   de \code{x} est associée.
#' @param alpha Numérique. Seuil de signification (\emph{p-value}) pour tous
#'   les tests (par défaut \code{0.05}).
#' @param normal Logique. \code{TRUE} si les données sont supposées normales,
#'   \code{FALSE} sinon.
#' @param number Entier ou NULL. Nombre de groupes. Si NULL (défaut), calculé
#'   automatiquement via \code{length(unique(g))}.
#' @param var.equal Logique. \code{TRUE} si la variance est supposée homogène,
#'   \code{FALSE} sinon (par défaut \code{FALSE}).
#' @param control Chaîne de caractères ou vecteur ou NULL. Nom(s) du(des)
#'   groupe(s) utilisé(s) comme contrôle(s). Si un contrôle est spécifié,
#'   la fonction réalise également un test de Dunnett
#'   (\code{\link[DescTools]{DunnettTest}}).
#' @param verbose Logique. \code{TRUE} (défaut) pour afficher des messages de
#'   suivi, \code{FALSE} sinon.
#' @param debug Logique. \code{TRUE} pour afficher des messages de débogage
#'   détaillés, \code{FALSE} (défaut) sinon.
#' @param code Logique. \code{TRUE} pour afficher des exemples de codes R
#'   équivalents, \code{FALSE} (défaut) sinon.
#' @param paired Logique. \code{TRUE} pour des données appariées, \code{FALSE}
#'   (défaut) sinon. Certains tests (p. ex. Wilcoxon ou \code{t.test}) tiennent
#'   alors compte de cette information.
#' @param boot Logique. \code{TRUE} (défaut) pour réaliser des validations par
#'   bootstrap, \code{FALSE} sinon.
#' @param iter Entier. Nombre d'itérations pour les procédures bootstrap
#'   (par défaut \code{1/alpha * 10}).
#' @param conf Numérique. Niveau de confiance pour les intervalles bootstrap
#'   (par défaut \code{0.95}).
#' @param k Entier ou NULL. Compteur interne pour les messages (usage avancé).
#'
#' @details
#' \strong{Stratégie de sélection des tests :}
#' \enumerate{
#'   \item \strong{Deux groupes seulement} (\code{number = 2}):
#'   \itemize{
#'     \item Si \code{normal = TRUE} : \code{t.test} (Student ou Welch selon
#'       \code{var.equal})
#'     \item Si \code{normal = FALSE} : \code{wilcox.test} (Wilcoxon-Mann-Whitney
#'       ou Wilcoxon apparié)
#'     \item Pour données appariées non normales : ajout du test de signe si
#'       asymétrie détectée (test MGG)
#'   }
#'   \item \strong{Plus de deux groupes} (\code{number > 2}):
#'   \itemize{
#'     \item Si \code{normal = TRUE} et \code{var.equal = TRUE} :
#'       \code{pairwise.t.test} avec pool.sd=TRUE, plus tests SNK, Tukey,
#'       Scheffé, Waller-Duncan
#'     \item Si \code{normal = TRUE} et \code{var.equal = FALSE} :
#'       \code{pairwise.t.test} avec pool.sd=FALSE, plus Games-Howell et Dunnett T3
#'     \item Si \code{normal = FALSE} : \code{pairwise.wilcox.test}, test de Dunn,
#'       et test de Brunner-Munzel si distributions différentes détectées (test KS)
#'     \item Pour données appariées non normales : ajout du test de Nemenyi
#'     \item Si valeurs extrêmes détectées : ajout du test Lincon sur données
#'       tronquées
#'   }
#' }
#'
#' \strong{Validation de la fiabilité des tests :}
#' \itemize{
#'   \item Pour Wilcoxon non apparié : test de Kolmogorov-Smirnov (avec correction
#'     de Sidak) pour vérifier l'identité des distributions
#'   \item Pour Wilcoxon apparié : test MGG (Miao-Gel-Gastwirth) pour vérifier
#'     la symétrie des différences
#'   \item Si conditions non remplies : ajout automatique de tests alternatifs
#'     robustes (Brunner-Munzel, test de signe)
#' }
#'
#' Dans tous les cas, si un contrôle est spécifié (\code{control}), la fonction
#' effectue également un test de Dunnett (\code{\link[DescTools]{DunnettTest}})
#' afin de comparer ce contrôle à chacun des autres groupes.
#'
#' Par ailleurs, si l'option \code{boot = TRUE}, la fonction exécute des
#' comparaisons par bootstrap via \code{.boots()}, et peut ainsi vérifier la
#' robustesse des tests classiques (en affichant des warnings si le bootstrap
#' détecte des incohérences).
#'
#' @return
#' \code{.posthoc} renvoie un objet de classe \code{"posthoc"}, une liste
#' contenant généralement :
#' \itemize{
#'   \item \code{groups} : un data.frame indiquant la catégorisation de chaque
#'     groupe selon la comparaison (lettres pour les groupes équivalents, étoiles
#'     pour les comparaisons au contrôle, etc.) ;
#'   \item \code{p.value} : la ou les p-values associées aux tests effectués
#'     (valeur unique pour 2 groupes, matrice pour >2 groupes) ;
#'   \item \code{bootstrap} : (optionnel) les résultats des comparaisons par
#'     bootstrap, si \code{boot = TRUE}.
#' }
#' Des éléments supplémentaires peuvent être ajoutés selon le test post-hoc réalisé
#' (p. ex. \code{Dunnett}, \code{SNK}, \code{Wilcoxon_holm}, etc.).
#'
#' @references
#' \itemize{
#'   \item Lehmann, E. L., & Romano, J. P. (2005). \emph{Testing Statistical
#'     Hypotheses} (3ᵉ éd.). Springer. (test de signe, Wilcoxon)
#'   \item Zheng, T., & Gastwirth, J. L. (2010). On bootstrap tests of symmetry
#'     about an unknown median. \emph{Journal of Data Science}, 8(3), 397–412.
#'     (test MGG)
#'   \item Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher
#'     problem: Asymptotic theory and a small-sample approximation.
#'     \emph{Biometrical Journal}, 42(1), 17-25. (test de Brunner-Munzel)
#' }
#'
#' @importFrom DescTools DunnettTest
#' @importFrom agricolae SNK.test HSD.test scheffe.test waller.test
#' @importFrom FSA dunnTest
#' @importFrom WRS2 lincon
#' @importFrom stringr str_split_fixed
#' @importFrom rstatix games_howell_test identify_outliers
#' @importFrom PMCMRplus dunnettT3Test frdAllPairsNemenyiTest
#' @importFrom stats t.test wilcox.test aov pairwise.t.test pairwise.wilcox.test
#'   quantile median sd binom.test
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{t.test}}, \code{\link{wilcox.test}} pour les tests de base,
#'   \item \code{\link[agricolae]{SNK.test}} pour le test de Newman-Keuls,
#'   \item \code{\link[FSA]{dunnTest}} pour le test de Dunn,
#'   \item \code{\link[DescTools]{DunnettTest}} pour le test de Dunnett,
#'   \item \code{\link[WRS2]{lincon}} pour les méthodes robustes (linéaires),
#'   \item \code{\link[rstatix]{games_howell_test}} pour le test de Games-Howell,
#'   \item \code{.boots()}, \code{pairwise()}, \code{catego()} (fonctions internes)
#'     pour la gestion du bootstrap et des comparaisons multiples.
#' }
#'
#' @examples
#' \dontrun{
#' # Exemple minimal (2 groupes, données normales, variances égales)
#' set.seed(123)
#' x <- rnorm(30)
#' g <- factor(rep(c("GroupeA", "GroupeB"), each = 15))
#' res <- .posthoc(x, g, alpha = 0.05, normal = TRUE, var.equal = TRUE)
#' print(res)
#'
#' # Exemple avec plus de 2 groupes et données non normales
#' set.seed(123)
#' x2 <- runif(40)
#' g2 <- factor(rep(c("G1", "G2", "G3", "G4"), each = 10))
#' res2 <- .posthoc(x2, g2, alpha = 0.01, normal = FALSE, var.equal = FALSE)
#' print(res2)
#'
#' # Exemple avec groupe contrôle
#' set.seed(456)
#' x3 <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 0.5))
#' g3 <- factor(rep(c("Control", "TreatA", "TreatB"), each = 10))
#' res3 <- .posthoc(x3, g3, alpha = 0.05, normal = TRUE,
#'                  var.equal = TRUE, control = "Control")
#' print(res3)
#' }
#'
#' @keywords internal
#' @export
.posthoc <- function(
    x, g,
    alpha = 0.05,
    normal = TRUE,
    number = NULL,
    var.equal = FALSE,
    control = NULL,
    verbose = TRUE,
    debug = FALSE,
    code = FALSE,
    paired = FALSE,
    boot = TRUE,
    iter = 1/alpha*10,
    conf = 0.95,
    k = NULL
) {
  # -- Coercitions --
  if (is.data.frame(x)) x <- x[[1]]
  if (is.matrix(x)) x <- x[, 1, drop = TRUE]
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.vector(x)
  if (!is.numeric(x)) stop(".posthoc() attend un vecteur numérique 'x'.")

  if (is.data.frame(g)) g <- g[[1]]
  if (is.matrix(g)) g <- g[, 1, drop = TRUE]
  if (is.list(g)) g <- unlist(g, use.names = FALSE)
  g <- droplevels(factor(g))
  lev <- levels(g)
  ng <- length(lev)

  if (is.null(number)) number <- ng

  if (length(x) != length(g)) stop("Longueurs incompatibles entre x et g.")

  if (ng < 2L || number < 2L) {
    synth <- list(
      groups = data.frame(categories = lev),
      note = "Pas assez de catégories (number < 2) pour réaliser des post-hoc."
    )
    class(synth) <- "posthoc"
    .dbg(NULL, "Fin de .posthoc() - moins de 2 groupes.",
         debug = debug)
    return(synth)
  }

  g1 <- lev[1L]
  g2 <- lev[2L]

  if (normal == TRUE) {
    .dbg(".posthoc() - Post-hoc on normal data.",
         ".posthoc() - Post-hoc sur données normales.",
         debug = debug)

    if (number == 2) {
      if (var.equal == TRUE) {
        ########################
        #  STUDENT (var.equal = TRUE)
        ########################
        if (code == TRUE) {
          cat("# Test de Student\n",
              "t.test(x[g==g1], x[g==g2], var.equal = TRUE, paired = ", paired, ")\n", sep = "")
        }

        pvals <- t.test(x[g == g1], x[g == g2], var.equal = TRUE, paired = paired)$p.value
        control_chr <- if (is.null(control)) NULL else as.character(control)[1]
        ind_control <- if (is.null(control_chr)) integer(0L) else match(control_chr, lev)

        if (isTRUE(paired)) {
          # --- CAS APPARIÉ : PAS de bootstrap "indice moyenne"
          k <- .vbse(
            "Post-hoc comparison of the two paired levels using a paired Student's t-test [t.test(paired = TRUE)].",
            "Comparaison post-hoc des 2 niveaux appariés par un test de Student [t.test(paired == TRUE)].",
            verbose = verbose, k = k, cpt = "on"
          )

          # Résultat (significatif / non significatif)
          if (pvals <= alpha) {
            k <- .vbse(
              paste0("Significant difference between the paired levels (p = ", .format_pval(pvals), ")."),
              paste0("Différence significative entre les niveaux appariés (p = ", .format_pval(pvals), ")."),
              verbose = verbose, k = k, cpt = "off"
            )
          } else {
            k <- .vbse(
              paste0("No significant difference between the paired levels (p = ", .format_pval(pvals), ")."),
              paste0("Aucune différence significative entre les niveaux appariés (p = ", .format_pval(pvals), ")."),
              verbose = verbose, k = k, cpt = "off"
            )
          }

          synth <- list()
          if (length(ind_control) != 1L) {
            # Sans témoin : lettres a/b selon la significativité
            synth$groups <- data.frame(
              categories = lev,
              Student_Holm = if (pvals <= alpha) c("a", "b") else c("a", "a")
            )
          } else {
            # Avec témoin : étoiles sur la catégorie non témoin
            stars <- c("", "")
            stars[-ind_control] <- ifelse(pvals <= 0.001, "***",
                                          ifelse(pvals <= 0.01,  "**",
                                                 ifelse(pvals <= 0.05,   "*", "")))
            synth$groups <- data.frame(
              categories = lev,
              Student_Holm = stars
            )
          }
          synth$p.value <- pvals

          # On peut garder le bootstrap "robustesse du Student" (optionnel) :
          if (isTRUE(boot)) {
            synth$bootstrap <- .boots(
              x, g,
              ctrl = (length(ind_control) == 1L),
              type = "mean",
              var.equal = TRUE,
              conf = conf, iter = iter, alpha = alpha,
              paired = TRUE, control = control
            )
            colnames(synth$bootstrap$groups)[2] <- "Student_bootstrapped"
          }

        } else {

          k <- .vbse(
            "a) Posthoc - Analysis of mean differences by bootstrap [pairwise.boot() from {KefiR}].",
            "a) Posthoc - Analyse des différences de moyennes par bootstrap [pairwise.boot() de {KefiR}]",
            verbose = verbose, k = k, cpt = "off"
          )

          synth2 <- pairwise(x, g, type = "boot", alpha = alpha, control = control, boot = FALSE,
                             boot_type = "mean", conf = conf, iter = iter, debug = debug)

          if (length(ind_control) != 1) {
            ########################
            #	Pas de control
            ########################
            # b) Student
            k <- .vbse(
              "b) Posthoc - Student test [t.test()].",
              "b) Posthoc - Test de Student [t.test()]",
              verbose = verbose, k = k, cpt = "off"
            )

            synth <- list()
            if (pvals <= alpha) {
              synth$groups <- data.frame(categories = lev, Student = c("a", "b"))
            } else {
              synth$groups <- data.frame(categories = lev, Student = c("a", "a"))
            }

            synth$p.value <- pvals

          } else {
            #######################
            #	control
            ########################
            .dbg(".posthoc() - Student with control.",
                 ".posthoc() - Student avec présence d'un Témoin.",
                 debug = debug, verbose = verbose)

            k <- .vbse(
              "b) Posthoc - Student test [t.test()] with control.",
              "b) Posthoc - Test de Student [t.test()] avec témoin.",
              verbose = verbose, k = k, cpt = "off"
            )

            synth <- list()
            stars <- c("", "")
            stars[-ind_control] <- ifelse(pvals <= 0.001, "***",
                                          ifelse(pvals <= 0.01, "**",
                                                 ifelse(pvals <= 0.05, "*", "")))
            synth$groups <- data.frame(categories = lev, Student = stars)
            synth$p.value <- pvals
          }

          # Bootstrap sur le test de Student
          if (boot == TRUE) {
            synth$bootstrap <- .boots(
              x, g,
              ctrl = (length(ind_control) == 1),
              type = "mean", var.equal = TRUE,
              conf = conf, iter = iter,
              alpha = alpha, paired = paired, control = control
            )
            colnames(synth$bootstrap$groups)[2] <- "Student_bootstrapped"
          }

          if ((verbose == TRUE) & (boot == TRUE) & any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            ang <- "Warning! Bootstrap detects weaknesses in the significance of the results."
            fr <- "Attention ! Le bootstrap détecte des faiblesses dans la signification des résultats."
            k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
          }

          # Ajout du bootstrap sur la moyenne
          ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
          synth$groups <- data.frame(
            "categories" = synth$groups[, 1],
            "Bootstrap_Mean" = synth2$groups[ind_temp, 2],
            "Student_Holm" = synth$groups[, 2]
          )
        }

        ########################
        #	STUDENT var.equal=FALSE - Test t de Welch
        ########################
      } else if (var.equal == FALSE) {
        if (code == TRUE) {
          cat("# Test t de Welch\nt.test(x[g==unique(g)[1]],x[g==unique(g)[2]],var.equal=FALSE, paired=", paired, ")\n")
        }

        pvals <- t.test(x[g == g1], x[g == g2], var.equal = FALSE, paired = paired)$p.value
        control_chr <- if (is.null(control)) NULL else as.character(control)[1]
        ind_control <- if (is.null(control_chr)) integer(0L) else match(control_chr, lev)

        # a) Bootstrap sur l'indice moyen
        k <- .vbse(
          "a) Posthoc - Analysis of mean differences by bootstrap [pairwise.boot(mu='mean') from {KefiR}].",
          "a) Posthoc - Analyse des différences de moyennes par bootstrap [pairwise.boot(mu='mean') de {KefiR}]",
          verbose = verbose, k = k, cpt = "off"
        )

        synth2 <- pairwise(x, g, type = "boot", alpha = alpha, control = control, boot = FALSE,
                           boot_type = "mean", conf = conf, iter = iter, debug = debug)

        if (length(ind_control) != 1) {
          ########################
          #	Pas de control
          ########################
          # b) Student
          k <- .vbse(
            "b) Posthoc - Student test [t.test(var.equal=FALSE)].",
            "b) Posthoc - Test de Student à variances inégales [t.test(var.equal=FALSE)]",
            verbose = verbose, k = k, cpt = "off"
          )

          synth <- list()
          if (pvals <= alpha) {
            synth$groups <- data.frame(categories = lev, Student = c("a", "b"))
          } else {
            synth$groups <- data.frame(categories = lev, Student = c("a", "a"))
          }

          synth$p.value <- pvals

          if (boot == TRUE) {
            synth$bootstrap <- .boots(
              x, g, ctrl = FALSE,
              type = "mean", var.equal = FALSE, conf = conf, iter = iter, alpha = alpha,
              paired = paired, control = control
            )
            colnames(synth$bootstrap$groups)[2] <- "Student_bootstrapped"
          }

          if ((verbose == TRUE) & (boot == TRUE) & any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            ang <- "Warning! Bootstrap detects weaknesses in the significance of the results."
            fr <- "Attention ! Le bootstrap détecte des faiblesses dans la signification des résultats."
            k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
          }

          # Ajout du bootstrap sur la moyenne
          ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
          synth$groups <- data.frame(
            "categories" = synth$groups[, 1],
            "Bootstrap_Mean" = synth2$groups[ind_temp, 2],
            "Student_Holm" = synth$groups[, 2]
          )

        } else {
          ########################
          #	control
          ########################
          k <- .vbse(
            "b) Posthoc - Student test [t.test(var.equal=FALSE)] with control.",
            "b) Posthoc - Test de Student à variances inégales [t.test(var.equal=FALSE)] avec témoin.",
            verbose = verbose, k = k, cpt = "off"
          )

          synth <- list()
          stars <- c("", "")
          stars[-ind_control] <- ifelse(pvals <= 0.001, "***",
                                        ifelse(pvals <= 0.01, "**",
                                               ifelse(pvals <= 0.05, "*", "")))
          synth$groups <- data.frame(categories = lev, Student = stars)
          synth$p.value <- pvals

          ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
          synth$groups <- data.frame(
            "categories" = synth$groups[, 1],
            "Bootstrap_Mean" = synth2$groups[ind_temp, 2],
            "Student_Holm" = synth$groups[, 2]
          )

          if (boot == TRUE) {
            synth$bootstrap <- .boots(
              x, g, ctrl = TRUE,
              type = "mean", var.equal = FALSE, conf = conf, iter = iter, alpha = alpha,
              paired = paired, control = control
            )
            colnames(synth$bootstrap$groups)[2] <- "Student_bootstrapped"
          }

          if ((verbose == TRUE) & (boot == TRUE) & any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            ang <- "Warning! Bootstrap detects weaknesses in the significance of the results."
            fr <- "Attention ! Le bootstrap détecte des faiblesses dans la signification des résultats."
            k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
          }
        }
      }

    } else if (number > 2) {
      .dbg(NULL, ".posthoc() - '>2 catégories'.",
           debug = debug)

      if (paired) {
        k <- .vbse(
          "Posthoc - Post-hoc tests on paired normal groups.",
          "Posthoc - Tests post-hoc sur groupes normaux appariés.",
          verbose = verbose, k = k
        )
      } else {
        k <- .vbse(
          "Posthoc - Post-hoc tests on normal groups.",
          "Posthoc - Tests post-hoc sur groupes normaux.",
          verbose = verbose, k = k
        )
      }

      ################
      # a) Bootstrap sur l'indice moyen
      k <- .vbse(
        "a) Posthoc - Analysis of mean differences by bootstrap [pairwise.boot() from {KefiR}].",
        "a) Posthoc - Analyse des différences de moyennes par bootstrap [pairwise.boot() de {KefiR}].",
        verbose = verbose, k = k, cpt = "off"
      )

      synth2 <- pairwise(x, g, type = "boot", alpha = alpha, control = control, boot = FALSE,
                         boot_type = "mean", conf = conf, iter = iter, debug = debug)

      ################
      # b) Réaliser un test de Student en pairwise en tenant compte du sd.pool.
      if (paired == TRUE) {
        if (code == TRUE) {
          message("# Tests t appariés\nresult <- pairwise.t.test(x,g,paired=TRUE)\n# Identification des groupes\nlibrary(KefiR)\ncatego(result)\n")
        }

        k <- .vbse(
          "b) Posthoc - Student test [pairwise.t.test(paired=TRUE)].",
          "b) Posthoc - Test de Student [pairwise.t.test(paired=TRUE)].",
          verbose = verbose, k = k, cpt = "off"
        )

        synth <- pairwise(x, g, type = "mean", pool.sd = FALSE, alpha = alpha, control = control, boot = boot,
                          conf = conf, iter = iter, paired = paired)

      } else if (paired == FALSE) {
        if (var.equal == FALSE) {
          if (code == TRUE) {
            message("# Tests t appariés\nresult <- pairwise.t.test(x,g,pool.sd=FALSE)\n# Identification des groupes\nlibrary(KefiR)\ncatego(result)\n")
          }

          k <- .vbse(
            "b) Posthoc - Student test [pairwise.t.test(pool.sd=FALSE)].",
            "b) Posthoc - Test de Student [pairwise.t.test(pool.sd=FALSE)].",
            verbose = verbose, k = k, cpt = "off"
          )

          synth <- pairwise(x, g, type = "mean", pool.sd = FALSE, alpha = alpha, control = control, boot = boot,
                            conf = conf, iter = iter, paired = paired)

        } else if (var.equal == TRUE) {
          if (code == TRUE) {
            message("# Tests t appariés\nresult <- pairwise.t.test(x,g,pool.sd=TRUE)\n# Identification des groupes\nlibrary(KefiR)\ncatego(result)\n")
          }

          ang <- "b) Posthoc - Student test [pairwise.t.test(pool.sd=TRUE)]."
          fr <- "b) Posthoc - Test de Student [pairwise.t.test(pool.sd=TRUE)]."
          k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")

          synth <- pairwise(x, g, type = "mean", pool.sd = TRUE, alpha = alpha, control = control, boot = boot,
                            conf = conf, iter = iter, paired = paired)
        }
      }

      colnames(synth$groups)[2] <- "Student_Holm"

      # Agglomérat du bootstrap median bootstrap sur colonne 2, ...
      ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
      synth$groups <- data.frame(
        "categories" = synth$groups[, 1],
        "Bootstrap_Mean" = synth2$groups[ind_temp, 2],
        "Student_Holm" = synth$groups[, 2]
      )

      # Plus bootstrap intégré
      if (boot == TRUE) {
        if (any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
          k <- .vbse(
            "\tWarning! Bootstrap detects weaknesses in the significance of the results.",
            "\tAttention ! Le bootstrap détecte des faiblesses dans la signification des résultats.",
            verbose = verbose, k = k, cpt = "off"
          )
        }
        colnames(synth$bootstrap$groups)[2] <- "Student_Holm_bootstrapped"
      }

      if ((var.equal == TRUE & length(control) == 0 & paired == FALSE)) {
        ################
        # c) Test de Scheffé
        if (code == TRUE) {
          message("# ANOVA à 1 facteur\nmya <- aov(x~g, data=data)\n\n#package nécessaire aux tests post-hocs\nlibrary(agricolae)\n\n#Test de Scheffé\nscheffe.test(mya,'g',alpha=", alpha, ")")
        }

        k <- .vbse(
          "c) Posthoc - Conservative Scheffé test [scheffe.test() of {agricolae}].",
          "c) Posthoc - Test conservateur de Scheffé [scheffe.test() of {agricolae}].",
          verbose = verbose, k = k, cpt = "off"
        )

        formula <- formula(x ~ g)
        data <- data.frame(x, g)
        data$g <- factor(data$g)
        mya <- suppressWarnings(aov(formula = formula, data = data))
        myscheffe <- scheffe.test(mya, "g", alpha = alpha)

        ind_temp <- match(synth$groups[, 1], rownames(myscheffe$groups))
        synth$groups <- data.frame(synth$groups, "Scheffe" = myscheffe$groups$groups[ind_temp])

        ################
        # d) Test de Tukey
        if (code == TRUE) {
          message("#Test de Tukey\nHSD.test(mya,'g',alpha=", alpha, ")")
        }

        k <- .vbse(
          "d) Posthoc - Conservative Tukey test [HSD.test() of {agricolae}].",
          "d) Posthoc - Test conservateur de Tukey [HSD.test() of {agricolae}].",
          verbose = verbose, k = k, cpt = "off"
        )

        mytukey <- HSD.test(mya, "g", alpha = alpha)

        ind_temp <- match(synth$groups[, 1], rownames(mytukey$groups))
        synth$groups <- data.frame(synth$groups, "Tukey" = mytukey$groups$groups[ind_temp])

        ################
        # e) Test de Newman-Keuls
        if (code == TRUE) {
          message("#Test de Newman-Keuls\nSNK.test(mya,'g',alpha=", alpha, ")")
        }

        k <- .vbse(
          "e) Posthoc - Powerful Newman-Keuls test [SNK.test() of {agricolae}].",
          "e) Posthoc - Test puissant de Newman-Keuls [SNK.test() of {agricolae}].",
          verbose = verbose, k = k, cpt = "off"
        )

        mynk <- SNK.test(mya, "g", alpha = alpha)

        ind_temp <- match(synth$groups[, 1], rownames(mynk$groups))
        synth$groups <- data.frame(synth$groups, "SNK" = mynk$groups$groups[ind_temp])

        cat1 <- unique(unlist(strsplit(synth$groups[, 3], "")))
        cat2 <- unique(unlist(strsplit(mynk$groups$groups, "")))
        rownames(synth$groups) <- rep(c(), nrow(synth$groups))

        control_chr <- if (is.null(control)) NULL else as.character(control)[1]
        ind_control <- if (is.null(control_chr)) integer(0L) else match(control_chr, lev)

        if ((length(cat1) != length(cat2)) & (length(ind_control) != 1)) {
          k <- .vbse(
            "Warning! pairwise.t.test() and SNK.test() don't return the same number of groups.",
            "Attention ! pairwise.t.test() et SNK.test() ne renvoient pas le même nombre de groupes.",
            k = k, cpt = "off"
          )
        }

        ################
        # f) Test de Waller-Duncan
        mywd <- tryCatch(waller.test(mya, "g"), error = function(e) NULL)
        if (!is.null(mywd)) {
          if (code == TRUE) {
            message("#Test de Waller-Duncan\nwaller.test(mya,'g')\n")
          }

          k <- .vbse(
            "f) Posthoc - Powerful Waller-Duncan test [waller.test() of {agricolae}].",
            "f) Posthoc - Test puissant de Waller-Duncan [waller.test() of {agricolae}].",
            verbose = verbose, k = k, cpt = "off"
          )

          ind_temp <- match(synth$groups[, 1], rownames(mywd$groups))
          synth$groups <- data.frame(synth$groups, "Waller-Duncan" = mywd$groups$groups[ind_temp])
        }

      } else if ((var.equal == TRUE) & (length(control) > 0) & (paired == FALSE)) {
        if (code == TRUE) {
          message("# Test de Dunnett car présence de Témoin\nlibrary(DescTools)\nDunnettTest(x,g,control='", control, "')\n")
        }

        k <- .vbse(
          "c) Posthoc - See also post-hoc Dunnett test for the control [DunnettTest() from {DescTools}].",
          "c) Posthoc - Voir aussi le test post-hoc de Dunnett pour le contrôle [DunnettTest() de {DescTools}].",
          verbose = verbose, k = k, cpt = "off"
        )

        p_value_to_symbol <- function(p_value) {
          if (p_value < 0.001) {
            return("***")
          } else if (p_value < 0.01) {
            return("**")
          } else if (p_value < 0.05) {
            return("*")
          } else {
            return("")
          }
        }

        dunnett_results <- DunnettTest(x, g, control = control)
        dunnett_results$Significance <- sapply(dunnett_results[[1]][, 4], p_value_to_symbol)
        decoupage <- str_split_fixed(names(dunnett_results$Significance), "-", 2)

        if (length(unique(decoupage[, 2])) == 1) {
          categories <- c(decoupage[, 1], unique(decoupage[, 2]))
          lettre <- c(dunnett_results$Significance, "")
        } else {
          .exit(
            "Le traitement de la sortie DunnettTest() renvoie une erreur. SVP, lancez .posthoc() avec l'argument debug = TRUE.",
            verbose = verbose
          )
        }

        lettre <- lettre[match(synth$groups$categories, categories)]
        synth$groups <- cbind(synth$groups, "Dunnett" = lettre)
        rownames(synth$groups) <- NULL

      } else if ((var.equal == FALSE) & (paired == FALSE)) {
        ################
        # c) Test de Games-Howell
        if (code == TRUE) {
          message("#Test de Games-Howell\nlibrary(rstatix)\ngames_howell_test(x~g, data=data)\n")
        }

        k <- .vbse(
          "c) Posthoc - Games-Howell test on normally distributed data with non-homogeneous variances [games_howell_test() of {rstatix}].",
          "c) Posthoc - Test de Games-Howell sur donnnées normales à variance non-homogènes...\n\t\t... [games_howell_test() of {rstatix}].",
          verbose = verbose, k = k, cpt = "off"
        )

        dt <- data.frame("x" = x, "g" = g)
        resultGH <- games_howell_test(x ~ g, data = dt)

        # Extraire les noms des classes uniques
        group_names <- unique(c(resultGH$group1, resultGH$group2))

        # Initialiser une matrice de NA avec les noms des groupes comme lignes et colonnes
        p_value_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                                 dimnames = list(group_names, group_names))

        # Remplir la matrice avec les p-values
        for (i in seq_len(nrow(resultGH))) {
          group1 <- resultGH$group1[i]
          group2 <- resultGH$group2[i]
          p_value <- resultGH$p.adj[i]
          p_value_matrix[group1, group2] <- p_value
          p_value_matrix[group2, group1] <- p_value  # La matrice est symétrique
        }

        p_value_matrix <- p_value_matrix[-1, -ncol(p_value_matrix)]

        # Afficher la matrice de p-values
        a <- list()
        a$p.value <- p_value_matrix
        synthGamesHowell <- catego(a)

        lettre <- synthGamesHowell$groups[, 2]
        lettre <- lettre[match(synth$groups$categories, synthGamesHowell$groups[, 1])]
        synth$groups <- cbind(synth$groups, "Games-Howell" = lettre)
        rownames(synth$groups) <- NULL

        ################
        # d) Test de Dunnett T3
        if (code == TRUE) {
          message("#Test de Dunnett T3\nlibrary(PMCMRplus)\nresult <-dunnettT3Test(mya)\nKefiR\ncatego(result)\n")
        }

        k <- .vbse(
          "d) Posthoc - Dunnett T3 test on normally distributed data with non-homogeneous variances [dunnettT3Test() of {PMCMRplus}].",
          "d) Posthoc - Test de Dunnett T3 sur donnnées normales à variance non-homogènes [dunnettT3Test() of {PMCMRplus}].",
          verbose = verbose, k = k, cpt = "off"
        )

        formula <- formula(x ~ g)
        data <- data.frame(x, g)
        data$g <- factor(data$g)
        mya <- suppressWarnings(aov(formula = formula, data = data))
        DunnettT3 <- catego(dunnettT3Test(mya))

        lettre <- DunnettT3$groups[, 2]
        lettre <- lettre[match(synth$groups$categories, DunnettT3$groups[, 1])]
        synth$groups <- cbind(synth$groups, "Dunnett T3" = lettre)
        rownames(synth$groups) <- NULL
      }
    }

  } else if (normal == FALSE) {
    #================================================
    #   Non-paramétrique
    #================================================
    .dbg(NULL, "post-hoc sur données non normales - k==2 et k>2.",
         debug = debug)

    check_wilcox_fiability <- TRUE

    if (paired == FALSE) {
      #------------------------------------------
      # Contrôle de la distribution pour évaluer la fiabilité du test de Wilcoxon-Mann-Witney...
      #------------------------------------------
      # Centrer sur la médiane pour fiabilisé le contrôle des équivalences de distributions par KS
      sd_cr <- by(x, g, sd, na.rm = T)
      median_cr <- by(x, g, median, na.rm = T)
      data_cr <- x
      for (i in names(median_cr)) {
        data_cr[g == i] <- (data_cr[g == i] - median_cr[which(names(median_cr) == i)]) / sd_cr[which(names(sd_cr) == i)]
      }

      # Vérifier l'équivalence de distribution pour évaluer la fiabilité de WMW si non appariées.
      # Lehmann, E. L., & Romano, J. P. (2005). Testing Statistical Hypotheses (3ᵉ éd.). Springer.
      # KS vérifie l'identité de distribution...
      # Correction de Sidak
      pval <- 1 - (1 - alpha)^(1/length(unique(g)))
      temp <- pairwise(data_cr, g, type = "ks", boot = FALSE)
      ks_result <- min(unlist(temp$p.value), na.rm = T)

      if (ks_result < pval) {
        check_wilcox_fiability <- FALSE
        ang <- paste0(
          "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data:\n\t",
          "Warning! The data do not have the same distribution (p-value: ", .format_pval(ks_result), ")...\n\t",
          "...compared to Sidak-corrected alpha ", .format_pval(pval), ".\n\t",
          "→ The Mann-Whitney-Wilcoxon test will be less reliable.\n\t",
          "Warning! For wilcox.test(): Please verify graphically that groups have the same distribution...\n\t",
          "...Or rely more on the Brunner-Munzel test."
        )
        fr <- paste0(
          "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites :\n\t",
          "Attention ! Les données n'ont pas la même distribution (p-value : ", .format_pval(ks_result), ")...\n\t",
          "...en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), ".\n\t",
          "→ Le test de Mann-Whitney-Wilcoxon sera moins fiable.\n\t",
          "Attention ! Pour wilcox.test() : Veuillez vérifier graphiquement que les groupes ont la même distribution...\n\t",
          "...Ou accordez plus de confiance au test de Brunner-Munzel."
        )
        k <- .vbse(ang, fr, verbose = verbose, k = k)
      } else {
        ang <- paste0(
          "Kolmogorov-Smirnov test [ks.test()] on median-centered and reduced data -\n\tThe groups have the same distribution. p-value: ", .format_pval(ks_result), "...\n\t\t...by comparing with Sidak corrected alpha ", .format_pval(pval), "\n\tThe Mann-Whitney-Wilcoxon test will be reliable."
        )
        fr <- paste0(
          "Test de Kolmogorov-Smirnov [ks.test()] sur les données centrées sur la médiane et réduites -\n\tLes groupes ont la même distribution. p-value : ", .format_pval(ks_result), "...\n\t\t...en comparant avec l'alpha corrigé de Sidak ", .format_pval(pval), "\n\tLe test de Mann-Whitney-Wilcoxon sera fiable a priori."
        )
        k <- .vbse(ang, fr, verbose = verbose, k = k)
      }

      #------------------------------------------
      # Annoncer les tests post-hocs avec une formulation adapté au nombre de groupes
      #------------------------------------------
      if (number == 2) {
        ang <- "Non-parametric comparison tests for two groups."
        fr <- "Tests non paramétriques de comparaison des deux groupes."
        k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "on")
      } else {
        ang <- "Non-parametric post-hoc tests for group comparisons."
        fr <- "Tests post-hocs non paramétriques de comparaison des groupes."
        k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "on")
      }

    } else if (number == 2) {
      # Si 2 groupes appariés :
      # Vérifier la symétrie des différences pour évaluer la fiabilité de WMW si appariées - Miao–Gel–Gastwirth MGG
      # Zheng, T., & Gastwirth, J. L. (2010). On bootstrap tests of symmetry about an unknown median. Journal of Data Science, 8(3), 397–412.
      # Le test MGG, robuste et bootstrapable, permet de vérifier la symétrie requise pour que le test de Wilcoxon apparié (WMW) soit statistiquement valide.
      k <- .vbse(
        "Symmetry check of paired differences - Miao-Gel-Gastwirth (MGG) test.",
        "Contrôle de la symétrie des différences appariées - Test de Miao-Gel-Gastwirth (MGG).",
        verbose = verbose, k = k
      )

      # Calculer les différences
      unique_g <- levels(g)
      g1 <- unique_g[1]
      g2 <- unique_g[2]
      x1 <- x[g == g1]
      x2 <- x[g == g2]
      differences <- x1 - x2

      # Fonction pour le test MGG (bootstrap)
      mgg_test <- function(diff, iter = 1000, alpha = 0.05) {
        n <- length(diff)
        med_diff <- median(diff)
        # Centrer sur la médiane
        centered_diff <- diff - med_diff
        # Statistique observée : moyenne des différences centrées
        T_obs <- mean(centered_diff)
        # Bootstrap sous H0 de symétrie
        T_boot <- numeric(iter)
        for (i in 1:iter) {
          # Sous H0, les signes sont aléatoires
          signs <- sample(c(-1, 1), n, replace = TRUE)
          centered_boot <- abs(centered_diff) * signs
          T_boot[i] <- mean(centered_boot)
        }
        # P-value bilatérale
        p_value <- mean(abs(T_boot) >= abs(T_obs))
        return(list(
          statistic = T_obs,
          p.value = p_value,
          median_diff = med_diff
        ))
      }

      # Effectuer le test MGG
      mgg_result <- mgg_test(differences, iter = iter, alpha = alpha)

      # Interprétation
      if (mgg_result$p.value < alpha) {
        check_wilcox_fiability <- FALSE
        ang <- paste0(
          "MGG symmetry test - Asymmetric differences detected (p = ",
          .format_pval(mgg_result$p.value), ").\n\t",
          "Warning! The Wilcoxon signed-rank test may be less reliable.\n\t",
          "Consider using the Sign test as an alternative."
        )
        fr <- paste0(
          "Test de symétrie MGG - Différences asymétriques détectées (p = ",
          .format_pval(mgg_result$p.value), ").\n\t",
          "Attention ! Le test de Wilcoxon apparié peut être moins fiable.\n\t",
          "Envisagez d'utiliser le test de signe comme alternative."
        )
        k <- .vbse(ang, fr, verbose = verbose, k = k, cpt = "off")
      }
    }

    ##########
    # Échantillons non-normaux - WILCOXON
    ##########
    if (number == 2) {
      .dbg(NULL, "2 catégories.",
           debug = debug)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Code
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (code == TRUE) {
        message("# Test de Wilcoxon-Mann-Witney\nwilcox.test(x[g==unique(g)[1]],x[g==unique(g)[2]],paired=", paired, ")")
      }

      pvals <- suppressWarnings(wilcox.test(x[g == g1], x[g == g2], paired = paired))$p.value
      control_chr <- if (is.null(control)) NULL else as.character(control)[1]
      ind_control <- if (is.null(control_chr)) integer(0L) else match(control_chr, lev)

      ##########
      # Échantillons appariés à 2 niveaux
      ##########
      if (isTRUE(paired)) {
        #=====================================================
        # --- CAS APPARIÉ : PAS de bootstrap "indice médiane"
        #=====================================================
        k <- .vbse(
          "Post-hoc: paired Wilcoxon signed-rank test [wilcox.test(paired = TRUE)].",
          "Post-hoc : test de Wilcoxon apparié [wilcox.test(paired = TRUE)].",
          verbose = verbose, k = k, cpt = "on"
        )

        if (pvals <= alpha) {
          k <- .vbse(
            paste0("Significant difference between the paired levels (p = ", .format_pval(pvals), ")."),
            paste0("Différence significative entre les niveaux appariés (p = ", .format_pval(pvals), ")."),
            verbose = verbose, k = k, cpt = "off"
          )
        } else {
          k <- .vbse(
            paste0("No significant difference between the paired levels (p = ", .format_pval(pvals), ")."),
            paste0("Aucune différence significative entre les niveaux appariés (p = ", .format_pval(pvals), ")."),
            verbose = verbose, k = k, cpt = "off"
          )
        }

        # Attribution des lettres et des étoiles
        synth <- list()
        if (length(ind_control) != 1L) {
          # Sans témoin : lettres a/b
          synth$groups <- data.frame(
            categories = lev,
            Wilcoxon_Holm = if (pvals <= alpha) c("a", "b") else c("a", "a")
          )
        } else {
          # Avec témoin : étoiles sur la catégorie non témoin
          stars <- c("", "")
          stars[-ind_control] <- ifelse(pvals <= 0.001, "***",
                                        ifelse(pvals <= 0.01,  "**",
                                               ifelse(pvals <= 0.05,   "*", "")))
          synth$groups <- data.frame(
            categories = lev,
            Wilcoxon_Holm = stars
          )
        }
        synth$p.value <- pvals

        # (optionnel) robustesse bootstrap sur la médiane
        if (isTRUE(boot)) {
          synth$bootstrap <- .boots(
            x, g,
            ctrl = (length(ind_control) == 1L),
            type = "median",
            conf = conf, iter = iter, alpha = alpha,
            paired = TRUE, control = control
          )
          colnames(synth$bootstrap$groups)[2] <- "Wilcoxon_bootstrapped"

          if ((verbose == TRUE) && any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            k <- .vbse(
              "Warning! Bootstrap detects weaknesses in the significance of the results.",
              "Attention ! Le bootstrap détecte des faiblesses dans la signification des résultats.",
              verbose = verbose, k = k, cpt = "off"
            )
          }
        }

        #=====================================================
        # --- check_wilcox_fiability == FALSE, compléter par un test de signe avec binom.test()
        #=====================================================
        if (check_wilcox_fiability == FALSE) {
          #========================================
          # CAS APPARIÉ : Différences asymétriques détectées par MGG
          # Alternative : Test de signe EN COMPLÉMENT du Wilcoxon
          #========================================
          k <- .vbse(
            "Asymmetric differences detected - Adding Sign test as robust alternative.",
            "Différences asymétriques détectées - Ajout du test de signe comme alternative robuste.",
            verbose = verbose, k = k
          )

          # Calculer les différences
          x1 <- x[g == g1]
          x2 <- x[g == g2]
          differences <- x1 - x2

          # Compter les différences (exclure les ex-aequo)
          n_positive <- sum(differences > 0)
          n_negative <- sum(differences < 0)
          n_zero <- sum(differences == 0)
          n_total <- n_positive + n_negative  # Différences non nulles

          if (n_total == 0) {
            k <- .vbse(
              "All differences are zero. Sign test cannot be performed.",
              "Toutes les différences sont nulles. Le test de signe ne peut être effectué.",
              verbose = verbose, k = k, cpt = "off"
            )
            # Ajouter une colonne vide
            synth$groups$Sign_test <- c("a", "a")
            # p-value = 1 (aucune différence)
            pvals_sign <- 1

          } else {
            # Test de signe avec binom.test()
            if (code == TRUE) {
              message("# Test de signe (Sign test)\n",
                      "differences <- x[g==unique(g)[1]] - x[g==unique(g)[2]]\n",
                      "n_positive <- sum(differences > 0)\n",
                      "n_total <- sum(differences != 0)\n",
                      "binom.test(n_positive, n_total, p = 0.5, alternative = 'two.sided')\n")
            }

            k <- .vbse(
              paste0("Sign test [binom.test()] - Distribution-free test for paired data.\n\t",
                     "Positive differences: ", n_positive, ", Negative: ", n_negative,
                     ", Ties (excluded): ", n_zero),
              paste0("Test de signe [binom.test()] - Test sans hypothèse de distribution pour données appariées.\n\t",
                     "Différences positives : ", n_positive, ", Négatives : ", n_negative,
                     ", Ex-aequo (exclus) : ", n_zero),
              verbose = verbose, k = k, cpt = "off"
            )

            # Effectuer le test binomial
            binom_result <- binom.test(n_positive, n_total, p = 0.5, alternative = "two.sided")
            pvals_sign <- binom_result$p.value

            # Attribution des lettres ou étoiles pour le test de signe
            if (length(ind_control) != 1) {
              #========================================
              # Pas de contrôle - lettres a/b
              #========================================
              if (pvals_sign <= alpha) {
                sign_groups <- c("a", "b")
                k <- .vbse(
                  paste0("Sign test: Significant difference detected (p = ", .format_pval(pvals_sign), ")."),
                  paste0("Test de signe : Différence significative détectée (p = ", .format_pval(pvals_sign), ")."),
                  verbose = verbose, k = k, cpt = "off"
                )
              } else {
                sign_groups <- c("a", "a")
                k <- .vbse(
                  paste0("Sign test: No significant difference (p = ", .format_pval(pvals_sign), ")."),
                  paste0("Test de signe : Aucune différence significative (p = ", .format_pval(pvals_sign), ")."),
                  verbose = verbose, k = k, cpt = "off"
                )
              }

            } else {
              #========================================
              # Avec contrôle - étoiles
              #========================================
              sign_groups <- c("", "")
              sign_groups[-ind_control] <- ifelse(pvals_sign <= 0.001, "***",
                                                  ifelse(pvals_sign <= 0.01, "**",
                                                         ifelse(pvals_sign <= 0.05, "*", "")))

              if (pvals_sign <= alpha) {
                k <- .vbse(
                  paste0("Sign test: Significant difference vs control (p = ", .format_pval(pvals_sign), ")."),
                  paste0("Test de signe : Différence significative vs témoin (p = ", .format_pval(pvals_sign), ")."),
                  verbose = verbose, k = k, cpt = "off"
                )
              } else {
                k <- .vbse(
                  paste0("Sign test: No significant difference vs control (p = ", .format_pval(pvals_sign), ")."),
                  paste0("Test de signe : Aucune différence significative vs témoin (p = ", .format_pval(pvals_sign), ")."),
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }

            # AJOUTER la colonne Sign_test à synth$groups (qui contient déjà Wilcoxon)
            synth$groups$Sign_test <- sign_groups

            # ÉCRASER synth$p.value avec la p-value du Sign test (test de référence car asymétrie)
            synth$p.value <- pvals_sign

            # Comparaison entre Wilcoxon et Sign test si divergence
            wilcox_signif <- (pvals <= alpha)
            sign_signif <- (pvals_sign <= alpha)

            if (wilcox_signif != sign_signif) {
              k <- .vbse(
                paste0("Warning! Wilcoxon and Sign tests give different conclusions.\n\t",
                       "Wilcoxon p=", .format_pval(pvals), " vs Sign p=", .format_pval(pvals_sign), "\n\t",
                       "Sign test is more reliable with asymmetric differences (synth$p.value updated)."),
                paste0("Attention ! Les tests de Wilcoxon et de signe donnent des conclusions différentes.\n\t",
                       "Wilcoxon p=", .format_pval(pvals), " vs Signe p=", .format_pval(pvals_sign), "\n\t",
                       "Le test de signe est plus fiable avec des différences asymétriques (synth$p.value mise à jour)."),
                verbose = verbose, k = k, cpt = "off"
              )
            }

            # Bootstrap optionnel pour robustesse du test de signe
            if (boot == TRUE) {
              k <- .vbse(
                "Bootstrap validation of Sign test.",
                "Validation bootstrap du test de signe.",
                verbose = verbose, k = k, cpt = "off"
              )

              boot_pvals <- numeric(iter)
              for (i in 1:iter) {
                # Rééchantillonner les paires
                indices <- sample(1:length(differences), replace = TRUE)
                diff_boot <- differences[indices]
                non_zero_boot <- diff_boot[diff_boot != 0]

                if (length(non_zero_boot) > 0) {
                  n_pos_boot <- sum(non_zero_boot > 0)
                  boot_pvals[i] <- binom.test(n_pos_boot, length(non_zero_boot),
                                              p = 0.5, alternative = "two.sided")$p.value
                } else {
                  boot_pvals[i] <- 1
                }
              }

              # Proportion significative
              boot_significant <- sum(boot_pvals <= alpha) / iter

              # Ajouter à synth$bootstrap (qui existe déjà pour Wilcoxon)
              if (is.null(synth$bootstrap)) {
                synth$bootstrap <- list()
              }
              synth$bootstrap$sign_proportion_significant <- boot_significant
              synth$bootstrap$sign_median_pvalue <- median(boot_pvals)

              # Ajouter colonne bootstrap Sign_test
              if (length(ind_control) != 1) {
                # Lettres selon le bootstrap
                sign_boot_groups <- if (boot_significant >= 0.95) sign_groups else c("a", "a")
              } else {
                # Étoiles selon le bootstrap
                sign_boot_groups <- if (boot_significant >= 0.95) sign_groups else c("", "")
              }

              # Créer ou compléter synth$bootstrap$groups
              if (is.null(synth$bootstrap$groups)) {
                synth$bootstrap$groups <- data.frame(categories = lev)
              }
              synth$bootstrap$groups$Sign_test_bootstrapped <- sign_boot_groups

              if (boot_significant < 0.95 && pvals_sign <= alpha) {
                k <- .vbse(
                  paste0("Warning! Bootstrap shows instability for Sign test: only ",
                         round(boot_significant * 100, 1), "% of iterations significant."),
                  paste0("Attention ! Bootstrap montre une instabilité pour le test de signe : seulement ",
                         round(boot_significant * 100, 1), "% d'itérations significatives."),
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }
          }

          # Note académique
          k <- .vbse(
            "Note: Sign test is the only UMP (Uniformly Most Powerful) invariant test when symmetry is not guaranteed (Lehmann & Romano, 2005).",
            "Note : Le test de signe est le seul test UMP (Uniformly Most Powerful) invariant lorsque la symétrie n'est pas garantie (Lehmann & Romano, 2005).",
            verbose = verbose, k = k, cpt = "off"
          )
        }

      } else {
        # Si non apparié
        #=====================================================
        #       Si non apparié
        #=====================================================

        # a) Bootstrap sur l'indice moyen
        k <- .vbse(
          "a) Posthoc - Analysis of median differences by bootstrap [pairwise.boot(mu='median') from {KefiR}].",
          "a) Posthoc - Analyse des différences de médianes par bootstrap [pairwise.boot(mu='median') de {KefiR}]",
          verbose = verbose, k = k, cpt = "off"
        )

        synth2 <- pairwise(x, g, type = "boot", alpha = alpha, control = control, boot = FALSE,
                           boot_type = "median", conf = conf, iter = iter, debug = debug)

        # b) Test de Wilcoxon
        k <- .vbse(
          "b) Post-hoc - Test de Wilcoxon-Mann-Whitney [wilcox.test()].",
          "b) Posthoc - Test de Wilcoxon-Mann-Witney [wilcox.test()].",
          verbose = verbose, k = k, cpt = "off"
        )

        if (length(ind_control) != 1) {
          ##########
          #	Pas de contrôle
          ##########
          synth <- list()
          if (pvals <= alpha) {
            synth$groups <- data.frame(categories = lev, Wilcoxon = c("a", "b"))
          } else {
            synth$groups <- data.frame(categories = lev, Wilcoxon = c("a", "a"))
          }
          synth$p.value <- pvals

        } else {
          ##########
          #	Contrôle
          ##########
          synth <- list()
          stars <- c("", "")
          stars[-ind_control] <- ifelse(pvals <= 0.001, "***",
                                        ifelse(pvals <= 0.01, "**",
                                               ifelse(pvals <= 0.05, "*", "")))
          synth$groups <- data.frame(categories = lev, Wilcoxon = stars)
          synth$p.value <- pvals
        }

        if (boot == TRUE) {
          synth$bootstrap <- .boots(
            x, g,
            ctrl = (length(ind_control) == 1),
            type = "median", conf = conf, iter = iter, alpha = alpha,
            paired = paired, control = control
          )
          colnames(synth$bootstrap$groups)[2] <- "Wilcoxon_bootstrapped"
        }

        if ((boot == TRUE) & any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
          k <- .vbse(
            "Warning! Bootstrap detects weaknesses in the significance of the results.",
            "Attention ! Le bootstrap détecte des faiblesses dans la signification des résultats.",
            k = k, cpt = "off"
          )
        }

        # Ajout du bootstrap sur la médiane
        ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
        synth$groups <- data.frame(
          "categories" = synth$groups[, 1],
          "Bootstrap_Median" = synth2$groups[ind_temp, 2],
          "Wilcoxon_Holm" = synth$groups[, 2]
        )

        #=====================================================
        # --- check_wilcox_fiability == FALSE
        # Distributions différentes détectées par KS
        # Alternative : Test de Brunner-Munzel
        #=====================================================
        if (check_wilcox_fiability == FALSE) {
          k <- .vbse(
            "Different distributions detected by KS test - Brunner-Munzel test added as robust alternative.",
            "Distributions différentes détectées par le test KS - Test de Brunner-Munzel ajouté comme alternative robuste.",
            verbose = verbose, k = k
          )

          # c) Test de Brunner-Munzel
          if (code == TRUE) {
            message("# Test de Brunner-Munzel\n",
                    "library(KefiR)\n",
                    "pairwise(x, g, type='BM', alpha=", alpha,
                    if (length(ind_control) == 1) paste0(", control='", control, "'") else "", ")\n")
          }

          # Appel à pairwise() avec type="BM" (fonctionne pour 2 groupes ou plus)
          synth_BM <- pairwise(x, g, type = "BM", alpha = alpha, control = control,
                               boot = boot, conf = conf, iter = iter, debug = debug)

          # Vérifier que pairwise a réussi
          if (!is.null(synth_BM) && !is.null(synth_BM$groups)) {

            # AJOUTER la colonne Brunner-Munzel aux groupes existants
            ind_temp <- match(synth$groups[, 1], synth_BM$groups[, 1])

            if (number == 2) {
              synth$groups$`Brunner-Munzel` <- synth_BM$groups[ind_temp, 2]
            } else {
              synth$groups$`Brunner-Munzel_Holm` <- synth_BM$groups[ind_temp, 2]
            }

            # ÉCRASER synth$p.value avec la p-value/matrice de Brunner-Munzel
            synth$p.value <- synth_BM$p.value

            # Comparaison Wilcoxon vs Brunner-Munzel (pour 2 groupes)
            if (number == 2) {
              wilcox_signif <- (pvals <= alpha)
              BM_signif <- (synth_BM$p.value <= alpha)

              if (wilcox_signif != BM_signif) {
                k <- .vbse(
                  paste0("Warning! Wilcoxon and Brunner-Munzel tests give different conclusions.\n\t",
                         "Wilcoxon p = ", .format_pval(pvals), " vs Brunner-Munzel p = ", .format_pval(synth_BM$p.value), "\n\t",
                         "→ Brunner-Munzel is more reliable with different distributions.\n\t",
                         "→ synth$p.value now contains Brunner-Munzel p-value (reference test)."),
                  paste0("Attention ! Les tests de Wilcoxon et de Brunner-Munzel donnent des conclusions différentes.\n\t",
                         "Wilcoxon p = ", .format_pval(pvals), " vs Brunner-Munzel p = ", .format_pval(synth_BM$p.value), "\n\t",
                         "→ Brunner-Munzel est plus fiable avec des distributions différentes.\n\t",
                         "→ synth$p.value contient maintenant la p-value de Brunner-Munzel (test de référence)."),
                  verbose = verbose, k = k, cpt = "off"
                )
              } else {
                k <- .vbse(
                  paste0("Wilcoxon and Brunner-Munzel tests agree (both ",
                         if (BM_signif) "significant" else "non-significant",
                         ", p-value (BM) = ", .format_pval(synth_BM$p.value), ")."),
                  paste0("Les tests de Wilcoxon et de Brunner-Munzel concordent (les deux ",
                         if (BM_signif) "significatifs" else "non significatifs",
                         ", p-value (BM) = ", .format_pval(synth_BM$p.value), ")."),
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }

            # Ajouter le bootstrap de BM si disponible
            if (boot == TRUE && !is.null(synth_BM$bootstrap)) {

              # Compléter synth$bootstrap$groups avec les résultats BM
              if (is.null(synth$bootstrap$groups)) {
                synth$bootstrap$groups <- data.frame(categories = lev)
              }

              ind_temp_boot <- match(synth$bootstrap$groups[, 1], synth_BM$bootstrap$groups[, 1])

              if (number == 2) {
                synth$bootstrap$groups$`Brunner-Munzel_bootstrapped` <- synth_BM$bootstrap$groups[ind_temp_boot, 2]
              } else {
                synth$bootstrap$groups$`Brunner-Munzel_Holm_bootstrapped` <- synth_BM$bootstrap$groups[ind_temp_boot, 2]
              }

              # Vérifier si divergence entre Wilcoxon bootstrap et BM bootstrap
              if (any(synth$bootstrap$groups[, 2] != synth_BM$bootstrap$groups[ind_temp_boot, 2], na.rm = TRUE)) {
                k <- .vbse(
                  "\tWarning! Bootstrap on Brunner-Munzel shows differences with Wilcoxon bootstrap.",
                  "\tAttention ! Bootstrap sur Brunner-Munzel montre des différences avec le bootstrap de Wilcoxon.",
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }

          } else {
            k <- .vbse(
              "Error! pairwise(type='BM') failed. Keeping only Wilcoxon results.",
              "Erreur ! pairwise(type='BM') a échoué. Conservation uniquement des résultats de Wilcoxon.",
              verbose = verbose, k = k, cpt = "off"
            )
          }

          # Note académique
          k <- .vbse(
            "Note: Brunner-Munzel test does not assume identical distributions, unlike Wilcoxon-Mann-Whitney.",
            "Note : Le test de Brunner-Munzel ne suppose pas de distributions identiques, contrairement à Wilcoxon-Mann-Whitney.",
            verbose = verbose, k = k, cpt = "off"
          )
        }
      }

      #=======================================================
      #         Plus de 2 catégories, non paramétriques
      #=======================================================
    } else if (number > 2) {
      .dbg(NULL, ".posthoc() - '>2 catégories'.",
           debug = debug)
      # Détection des outliers
      outlier <- function(z) {
        tablo_outlier <- identify_outliers(data.frame(z))
        if (nrow(tablo_outlier) > 0) {
          small_outliers <- tablo_outlier[tablo_outlier[, 1] < median(z), 1]
          big_outliers <- tablo_outlier[tablo_outlier[, 1] > median(z), 1]
          small_outliers <- length(small_outliers) / length(z)
          big_outliers <- length(big_outliers) / length(z)
          return(max(c(small_outliers, big_outliers)))
        } else {
          return(0)
        }
      }

      #------------------------------------------
      # 0) Identifier les outliers
      #------------------------------------------
      trimmage <- max(by(x, g, outlier))
      if (trimmage > 0) {
        k <- .vbse(
          "Posthoc - Warning! Groups with extreme values:\n\t\trisk of leverage (identify_outliers() of {rstatix}).",
          "Posthoc - Attention ! Groupes avec des valeurs extrêmes :\n\t\trisque d'effet de levier (identify_outliers() de {rstatix}).",
          verbose = verbose, k = k, cpt = "off"
        )
      }

      #------------------------------------------
      # 0) Croisements - Voir le nombre de paires à analyser...
      #------------------------------------------
      croisement <- ng * (ng - 1) / 2
      if (croisement > 28) {
        k <- .vbse(
          paste0("Posthoc - Warning! Many cross-tests to perform.\n\t\t", ng, " categories for ", croisement, " pairwise comparisons."),
          paste0("Posthoc - Attention ! Beaucoup de tests croisés à réaliser.\n\t\t", ng, " catégories pour ", croisement, " croisements"),
          verbose = verbose, k = k, cpt = "off"
        )
      }

      #------------------------------------------
      # 0) Bootstrap sur l'indice median (Wilcoxon)
      #------------------------------------------
      if (code == TRUE) {
        message("#Analyse des différences de médianes par bootstrap\nlibrary(KefiR)pairwise(x,g,type='boot')")
      }

      k <- .vbse(
        "Posthoc - a) Analysis of median differences using bootstrap [pairwise.boot() from {KefiR}].",
        "Posthoc - a) Analyse des différences de médianes par bootstrap [pairwise.boot() de {KefiR}].",
        verbose = verbose, k = k, cpt = "off"
      )

      synth2 <- pairwise(x, g, type = "boot", alpha = alpha, control = control, boot = FALSE,
                         boot_type = "median", conf = conf, iter = iter, debug = debug)

      # 1) Faire un Wilcoxon (position 1) en pairwise avec Holm
      if (code == TRUE) {
        message("# Wilcoxon par paires#pairwise.wilcox.test(x,g,p.adjust.method='BH')\n")
      }

      k <- .vbse(
        "Posthoc - b) Wilcoxon-Mann-Whitney test [pairwise.wilcox.test()].",
        "Posthoc - b) Test de Wilcoxon-Mann-Whitney [pairwise.wilcox.test()].",
        verbose = verbose, k = k, cpt = "off"
      )

      synth <- pairwise(x, g, type = "median", alpha = alpha, control = control, boot = boot, conf = conf, iter = iter, paired = paired, debug = debug)
      colnames(synth$groups)[2] <- "Wilcoxon_Holm"

      # Aggloméra du bootstrap median bootstrap sur colonne 2, Wilcoxon sur colonne 3
      ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
      synth$groups <- data.frame(
        "categories" = synth$groups[, 1],
        "Bootstrap_Median" = synth2$groups[ind_temp, 2],
        "Wilcoxon_Holm" = synth$groups[, 2]
      )

      # Plus bootstrap intégré
      if (boot == TRUE) {
        if (any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
          k <- .vbse(
            "\tWarning! Bootstrap detects weaknesses in the significance of the results.",
            "\tAttention ! Le bootstrap détecte des faiblesses dans la signification des résultats.",
            verbose = verbose, k = k, cpt = "off"
          )
        }
        colnames(synth$bootstrap$groups)[2] <- "Wilcoxon_Holm_bootstrapped"
      }

      if (paired == TRUE) {
        if (code == TRUE) {
          message("# Nemenyi - post-hoc sur données appariées\nlibrary(PMCMRplus)\nfrdAllPairsNemenyiTest(x~g)\n")
        }

        k <- .vbse(
          "Posthoc - c) Paired groups [frdAllPairsNemenyiTest() of {PMCMRplus}].",
          "c) Posthoc - Test post-hoc sur groupes non normaux appariés [frdAllPairsNemenyiTest() of {PMCMRplus}].",
          verbose = verbose, k = k, cpt = "off"
        )

        # Données appariées, faire un Neminye
        synth2 <- pairwise(x, g, type = "neminye", alpha = alpha, control = control, boot = boot, conf = conf, iter = iter, paired = paired, debug = debug)

        ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
        colnames(synth2$groups)[2] <- "Neminyi"
        synth$groups <- cbind(synth$groups, "Neminyi" = synth2$groups[ind_temp, 2])

        if (boot == TRUE) {
          if (any(synth2$bootstrap$groups[, 2] != synth2$groups[, 2])) {
            k <- .vbse(
              "\tWarning! Bootstrap on Neminyi- weaknesses.",
              "\tAttention ! Bootstrap sur Neminyi - faiblesses.",
              verbose = verbose, k = k, cpt = "off"
            )
          }
          colnames(synth2$bootstrap$groups)[2] <- "Neminyi_bootstrapped"
          synth$bootstrap$groups <- cbind(synth$bootstrap$groups, "Neminyi_bootstrapped" = synth2$bootstrap$groups[ind_temp, 2])
        }

      } else if ((trimmage > 0) & (croisement < 28)) {
        # Si détection d'outliers avec rstatix::identify_outliers()
        if (code == TRUE) {
          message("# Lincon\nlibrary(WRS2)\nlincon(x~g,tr)\n")
        }

        k <- .vbse(
          paste0("Posthoc - c) Post-hoc on trimmed data [lincon() from {WRS2}].\n\t\tLower and upper trimming level of ", round(trimmage * 100, 1), "%."),
          paste0("Posthoc - c) Post-hoc sur données tronquées [lincon() de {WRS2}].\n\t\tNiveau de troncature inférieure et supérieure de ", round(trimmage * 100, 1), "%."),
          verbose = verbose, k = k, cpt = "off"
        )

        # Faire un lincon() sur données trimmées
        synth2 <- pairwise(x, g, type = "lincon", tr = trimmage,
                           alpha = alpha, control = control, boot = boot,
                           conf = conf, iter = iter, debug = debug)

        ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
        colnames(synth2$groups)[2] <- "Lincon"
        synth$groups <- cbind(synth$groups, "Lincon" = synth2$groups[ind_temp, 2])

        if (boot == TRUE) {
          if (any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            k <- .vbse(
              "\tWarning! Bootstrap on Lincon - weaknesses.",
              "\tAttention ! Bootstrap sur Lincon - faiblesses.",
              verbose = verbose, k = k, cpt = "off"
            )
          }
          colnames(synth2$bootstrap$groups)[2] <- "Lincon_bootstrapped"
          synth$bootstrap$groups <- cbind(synth$bootstrap$groups, synth2$bootstrap$groups[ind_temp, 2])
        }

      } else if (check_wilcox_fiability == FALSE) {
        #=====================================================
        # K GROUPES - Distributions différentes détectées
        # Wilcoxon et Dunn non fiables
        # Alternative : Brunner-Munzel
        #=====================================================
        k <- .vbse(
          paste0("Warning! Different distributions detected by KS test for ", ng, " groups.\n\t",
                 "Wilcoxon and Dunn tests may be unreliable.\n\t",
                 "→ Adding Brunner-Munzel test as robust alternative."),
          paste0("Attention ! Distributions différentes détectées par le test KS pour ", ng, " groupes.\n\t",
                 "Les tests de Wilcoxon et Dunn peuvent être peu fiables.\n\t",
                 "→ Ajout du test de Brunner-Munzel comme alternative robuste."),
          verbose = verbose, k = k
        )

        if (code == TRUE) {
          message("# Test de Brunner-Munzel par paires avec correction de Holm\n",
                  "# Pour chaque paire de groupes :\n",
                  "library(lawstat)\n",
                  "# Exemple pour groupes 1 et 2 :\n",
                  "brunner.munzel.test(x[g==levels(g)[1]], x[g==levels(g)[2]])\n",
                  "# Répéter pour toutes les paires (i,j) où i < j\n",
                  "# Puis appliquer la correction de Holm sur toutes les p-values :\n",
                  "p.adjust(p_values_vector, method='holm')\n")
        }

        k <- .vbse(
          paste0("Brunner-Munzel pairwise test [brunner.munzel.test() from {lawstat} with Holm correction] - ",
                 ng, " groups, robust to different distributions."),
          paste0("Test de Brunner-Munzel par paires [brunner.munzel.test() de {lawstat} avec correction de Holm] - ",
                 ng, " groupes, robuste aux distributions différentes."),
          verbose = verbose, k = k, cpt = "off"
        )

        # Appel à pairwise() avec type="BM" (en interne, utilise lawstat::brunner.munzel.test)
        synth_BM <- pairwise(x, g, type = "BM", alpha = alpha, control = control,
                             boot = boot, conf = conf, iter = iter, debug = debug)

        # Vérifier que pairwise a réussi
        if (!is.null(synth_BM) && !is.null(synth_BM$groups)) {

          # AJOUTER la colonne Brunner-Munzel aux groupes existants (après Wilcoxon/Dunn/Lincon)
          ind_temp <- match(synth$groups[, 1], synth_BM$groups[, 1])
          synth$groups$`Brunner-Munzel_Holm` <- synth_BM$groups[ind_temp, 2]

          # ÉCRASER synth$p.value avec la matrice de p-values de Brunner-Munzel
          # (BM devient la référence car plus fiable avec distributions différentes)
          synth$p.value <- synth_BM$p.value

          k <- .vbse(
            "Brunner-Munzel test completed. p-values updated (synth$p.value now contains BM p-values matrix).",
            "Test de Brunner-Munzel terminé. p-values mises à jour (synth$p.value contient maintenant la matrice BM).",
            verbose = verbose, k = k, cpt = "off"
          )

          # Ajouter le bootstrap de BM si disponible
          if (boot == TRUE && !is.null(synth_BM$bootstrap)) {

            # Compléter synth$bootstrap$groups avec les résultats BM
            if (is.null(synth$bootstrap$groups)) {
              synth$bootstrap$groups <- data.frame(categories = lev)
            }

            ind_temp_boot <- match(synth$bootstrap$groups[, 1], synth_BM$bootstrap$groups[, 1])
            synth$bootstrap$groups$`Brunner-Munzel_Holm_bootstrapped` <- synth_BM$bootstrap$groups[ind_temp_boot, 2]

            # Vérifier si divergence entre les tests précédents et BM bootstrap
            # Comparer avec la première colonne de résultats (généralement Wilcoxon_Holm ou Dunn_Holm)
            if (ncol(synth$bootstrap$groups) >= 2) {
              if (any(synth$bootstrap$groups[, 2] != synth_BM$bootstrap$groups[ind_temp_boot, 2], na.rm = TRUE)) {
                k <- .vbse(
                  "\tWarning! Bootstrap on Brunner-Munzel shows differences with previous bootstrap results.",
                  "\tAttention ! Bootstrap sur Brunner-Munzel montre des différences avec les résultats bootstrap précédents.",
                  verbose = verbose, k = k, cpt = "off"
                )
              } else {
                k <- .vbse(
                  "\tBootstrap confirms consistency between tests.",
                  "\tBootstrap confirme la cohérence entre les tests.",
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }
          }

          # Comparaison des conclusions entre les différents tests
          # Vérifier si Wilcoxon/Dunn et BM donnent les mêmes groupes
          n_cols <- ncol(synth$groups)
          if (n_cols >= 3) {
            # Comparer la dernière colonne BM avec les précédentes (Wilcoxon, Dunn, etc.)
            differences_detected <- FALSE

            for (col_idx in 2:(n_cols - 1)) {
              if (any(synth$groups[, col_idx] != synth$groups$`Brunner-Munzel_Holm`, na.rm = TRUE)) {
                differences_detected <- TRUE
                break
              }
            }

            if (differences_detected) {
              k <- .vbse(
                paste0("Warning! Brunner-Munzel gives different groupings than other tests.\n\t",
                       "→ Brunner-Munzel is more reliable with different distributions.\n\t",
                       "→ Prioritize Brunner-Munzel results (rightmost column)."),
                paste0("Attention ! Brunner-Munzel donne des groupements différents des autres tests.\n\t",
                       "→ Brunner-Munzel est plus fiable avec des distributions différentes.\n\t",
                       "→ Prioriser les résultats de Brunner-Munzel (colonne la plus à droite)."),
                verbose = verbose, k = k, cpt = "off"
              )
            } else {
              k <- .vbse(
                "All tests agree on group classifications.",
                "Tous les tests concordent sur les classifications de groupes.",
                verbose = verbose, k = k, cpt = "off"
              )
            }
          }

        } else {
          k <- .vbse(
            "Error! Brunner-Munzel pairwise test failed. Keeping previous test results (Wilcoxon/Dunn/Lincon).",
            "Erreur ! Test de Brunner-Munzel par paires échoué. Conservation des résultats des tests précédents (Wilcoxon/Dunn/Lincon).",
            verbose = verbose, k = k, cpt = "off"
          )
        }

        # Note académique
        k <- .vbse(
          paste0("Academic note: With different distributions detected, Brunner-Munzel test is preferable.\n\t",
                 "Wilcoxon-Mann-Whitney and Dunn tests assume identical distributions (location shift only).\n\t",
                 "Brunner-Munzel makes no such assumption and remains valid (Brunner & Munzel, 2000)."),
          paste0("Note académique : Avec des distributions différentes détectées, le test de Brunner-Munzel est préférable.\n\t",
                 "Les tests de Wilcoxon-Mann-Whitney et Dunn supposent des distributions identiques (décalage de position uniquement).\n\t",
                 "Brunner-Munzel ne fait pas cette hypothèse et reste valide (Brunner & Munzel, 2000)."),
          verbose = verbose, k = k, cpt = "off"
        )

      } else {
        # Test de Dunn avec correction de Holm
        if (code == TRUE) {
          message("# Test de Dunn\nlibrary(FSA)\ndunnTest(x ~ g, method = 'holm')\n")
        }

        k <- .vbse(
          "Posthoc - c) Dunn test [dunnTest() of {FSA}].",
          "Posthoc - c) Test de Dunn [dunnTest() de {FSA}].",
          verbose = verbose, k = k, cpt = "off"
        )

        if (trimmage > 0) {
          k <- .vbse(
            paste0("\tLower and upper trimming level of ", round(trimmage * 100, 1), "%."),
            paste0("\tAvec troncature inférieure et supérieure de ", round(trimmage * 100, 1), "%."),
            verbose = verbose, k = k, cpt = "off"
          )

          x_temp <- c()
          g_temp <- c()
          for (categor in lev) {
            x_categor <- x[g == categor]
            x_categor <- x_categor[x_categor > quantile(x_categor, p = trimmage) & x_categor < quantile(x_categor, p = (1 - trimmage))]
            g_categor <- rep(categor, length(x_categor))
            x_temp <- c(x_temp, x_categor)
            g_temp <- c(g_temp, g_categor)
          }
        } else {
          x_temp <- x
          g_temp <- g
        }

        synth2 <- pairwise(x_temp, g_temp, type = "dunn", alpha = alpha, control = control, boot = boot,
                           conf = conf, iter = iter, debug = debug)

        ind_temp <- match(synth$groups[, 1], synth2$groups[, 1])
        synth$groups <- data.frame(synth$groups, "Dunn_Holm" = synth2$groups[ind_temp, 2])

        if (boot == TRUE) {
          if (any(synth$bootstrap$groups[, 2] != synth$groups[, 2])) {
            k <- .vbse(
              "\tWarning! Bootstrap on Dunn test detects weaknesses.",
              "\tAttention ! Bootstrap sur le test de Dunn détecte des faiblesses.",
              verbose = verbose, k = k, cpt = "off"
            )
          }
          synth$bootstrap$groups <- data.frame(synth$bootstrap$groups, "Dunn_Holm_bootstrapped" = synth2$bootstrap$groups[ind_temp, 2])
        }
      }
    }
  }

  synth <- structure(
    synth,
    class = "posthoc"
  )

  .dbg(NULL, "Fin de .posthoc().",
       debug = debug)

  return(synth)
}
