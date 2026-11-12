# Socle bayesien commun : valeurs calculees a la main, commentees.
# Ces tests n'exigent ni BayesFactor ni bayestestR.

# Jeu de log BF connu, repris de dev/mbay/scripts/resultats/agent_02_bayesfactor.txt
# (y ~ A*B, BF contre le modele nul).
bf_null <- 1
bf_a    <- 18.1343746
bf_b    <- 0.3145700
bf_ab   <- 6.4448757   # A + B
bf_full <- 2.9305522   # A + B + A:B
space_2f <- list(character(0), "A", "B", c("A", "B"), c("A", "B", "A:B"))
logs_2f  <- log(c(bf_null, bf_a, bf_b, bf_ab, bf_full))

row_of <- function(res, term) res[res$term == term, , drop = FALSE]

# ---- .bay_post_prob --------------------------------------------------------

test_that("P(H1|D) : valeurs a la main du paragraphe 8 (BF = 10)", {
  # O*BF/(1+O*BF) : 2,5/3,5 ; 10/11 ; 40/41.
  expect_equal(.bay_post_prob(10, 1/4), 2.5 / 3.5)
  expect_equal(.bay_post_prob(10, 1), 10 / 11)
  expect_equal(.bay_post_prob(10, 4), 40 / 41)
  expect_equal(round(.bay_post_prob(10, c(1/4, 1, 4)), 3), c(0.714, 0.909, 0.976))
})

test_that("P(H1|D) : BF = 1 rend l'a priori (pas de mise a jour)", {
  # O/(1+O) = pi, quel que soit O.
  odds <- c(1/4, 1, 3, 9)
  expect_equal(.bay_post_prob(rep(1, 4), odds), odds / (1 + odds))
})

test_that("P(H1|D) : cas limites BF = 0, Inf, NA, NaN", {
  expect_identical(.bay_post_prob(0), 0)
  expect_identical(.bay_post_prob(Inf), 1)
  expect_true(is.na(.bay_post_prob(NA)))
  expect_true(is.na(.bay_post_prob(NaN)))
  expect_false(is.nan(.bay_post_prob(NaN)))
  expect_equal(.bay_post_prob(c(0, 1, NA, Inf)), c(0, 0.5, NA, 1))
  expect_identical(.bay_post_prob(numeric(0)), numeric(0))
})

test_that("P(H1|D) : log_bf equivaut a bf, et les noms sont conserves", {
  bf <- c(a = 0.1, b = 3, c = 250)
  expect_equal(.bay_post_prob(log_bf = log(bf), prior_odds = 2),
               .bay_post_prob(bf, prior_odds = 2))
  expect_named(.bay_post_prob(bf), c("a", "b", "c"))
})

test_that("P(H1|D) : pas de debordement pour d'enormes log BF", {
  # exp(1000) deborde : l'echelle naturelle donnerait Inf/Inf = NaN.
  expect_true(is.nan(exp(1000) / (1 + exp(1000))))
  expect_identical(.bay_post_prob(log_bf = 1000), 1)
  expect_identical(.bay_post_prob(log_bf = -1e6), 0)
  # A log BF = -700, P = exp(-700)/(1+exp(-700)) ~ exp(-700).
  expect_equal(log(.bay_post_prob(log_bf = -700)), -700)
  # La petite queue P(H0|D) reste exacte : a log BF = 40, 1 - p donnerait 0.
  expect_equal(.bay_post_prob(log_bf = -40), exp(-40) / (1 + exp(-40)))
  expect_identical(1 - .bay_post_prob(log_bf = 40), 0)
})

test_that("P(H1|D) : entrees invalides refusees", {
  expect_error(.bay_post_prob(-1))
  expect_error(.bay_post_prob(c(1, -0.1, NA)))
  expect_error(.bay_post_prob(1, 0))
  expect_error(.bay_post_prob(1, -2))
  expect_error(.bay_post_prob(1, NA))
  expect_error(.bay_post_prob(1, Inf))
  expect_error(.bay_post_prob(1, numeric(0)))
  expect_error(.bay_post_prob(1:3, c(1, 2)))
  expect_error(.bay_post_prob())
  expect_error(.bay_post_prob(1, log_bf = 0))
  expect_error(.bay_post_prob("a"))
})

test_that("P(H1|D) + P(H0|D) = 1, y compris par la voie symetrique", {
  log_bf <- c(-800, -50, -1, 0, 1, 50, 800)
  for (o in c(1/9, 1, 7)) {
    p1 <- .bay_post_prob(log_bf = log_bf, prior_odds = o)
    p0 <- .bay_post_prob(log_bf = -log_bf, prior_odds = 1 / o)
    expect_equal(p1 + p0, rep(1, length(log_bf)))
  }
})

# ---- ROPE ------------------------------------------------------------------

test_that("masse a priori du ROPE : valeurs du paragraphe 4", {
  m <- .bay_prior_mass_rope(r = 0.707, w = 0.1, prior_h1 = 1/2)
  # 2/pi*atan(0,1/0,707) ; 0,5 + 0,5*0,0894.
  expect_equal(round(m$mass_h1, 3), 0.089)
  expect_equal(round(m$mass_total, 3), 0.545)
  # La Cauchy a pour fonction de repartition pcauchy : autre ecriture de la masse.
  expect_equal(m$mass_h1, pcauchy(0.1, 0, 0.707) - pcauchy(-0.1, 0, 0.707))
})

test_that("masse a priori du ROPE : valeurs remarquables et limites", {
  # w = r : atan(1) = pi/4, donc masse = 1/2 (quartiles de la Cauchy).
  expect_equal(.bay_prior_mass_rope(0.5, 0.5)$mass_h1, 1 / 2)
  # w = r, pi = 1/2 : 1/2 + 1/2*1/2 = 3/4.
  expect_equal(.bay_prior_mass_rope(0.5, 0.5)$mass_total, 3 / 4)
  expect_equal(.bay_prior_mass_rope(1, 0)$mass_h1, 0)
  expect_equal(.bay_prior_mass_rope(1, 0)$mass_total, 1 / 2)
  expect_equal(.bay_prior_mass_rope(1, Inf)$mass_h1, 1)
  expect_equal(.bay_prior_mass_rope(1, 1, prior_h1 = 1)$mass_total, 1 / 2)
  expect_equal(.bay_prior_mass_rope(1, 1, prior_h1 = 0)$mass_total, 1)
  expect_equal(nrow(.bay_prior_mass_rope(c(0.5, 1), 0.1)), 2L)
  expect_error(.bay_prior_mass_rope(0, 0.1))
  expect_error(.bay_prior_mass_rope(-1, 0.1))
  expect_error(.bay_prior_mass_rope(1, -0.1))
  expect_error(.bay_prior_mass_rope(1, 0.1, prior_h1 = 1.2))
  expect_error(.bay_prior_mass_rope(NA, 0.1))
})

test_that("P(ROPE|D) : formule et cas limites", {
  # P(H0|D) + P(H1|D)*q : 0,09 + 0,91*0,3 = 0,363.
  expect_equal(.bay_rope_prob(0.91, 0.3), 0.363)
  expect_equal(.bay_rope_prob(1, 0.25), 0.25)
  expect_equal(.bay_rope_prob(0.5, 1), 1)
  expect_equal(.bay_rope_prob(0.5, 0), 0.5)
  # P(H1|D) = 0 : le conditionnel ne compte pas, meme manquant.
  expect_identical(.bay_rope_prob(0, NA), 1)
  expect_true(is.na(.bay_rope_prob(0.5, NA)))
  expect_true(is.na(.bay_rope_prob(NA, 0.5)))
  expect_equal(.bay_rope_prob(c(0.2, 0.8), 0.5), c(0.9, 0.6))
  expect_error(.bay_rope_prob(1.1, 0.5))
  expect_error(.bay_rope_prob(0.5, -0.1))
  expect_error(.bay_rope_prob(c(0.1, 0.2), c(0.1, 0.2, 0.3)))
})

# ---- .bay_inclusion : cas A*B ----------------------------------------------

test_that("inclusion appariee A*B : recalcul independant en echelle naturelle", {
  res <- .bay_inclusion(logs_2f, space_2f, "matched")
  # Ecrit autrement : moyennes arithmetiques sur M+ et M-, sans passer par les logs.
  # A : M+ = {A, A+B}, M- = {nul, B}.
  # B : M+ = {B, A+B}, M- = {nul, A}.
  # A:B : M+ = {A+B+A:B}, M- = {A+B}.
  hand <- c(A    = mean(c(bf_a, bf_ab)) / mean(c(bf_null, bf_b)),
            B    = mean(c(bf_b, bf_ab)) / mean(c(bf_null, bf_a)),
            `A:B` = bf_full / bf_ab)
  expect_equal(res$term, c("A", "B", "A:B"))
  expect_equal(setNames(res$bf10, res$term), hand)
  # Valeurs publiees dans agent_02 : 18,7 ; 0,3533 ; 0,4547.
  expect_equal(round(hand, 4), c(A = 18.6970, B = 0.3533, `A:B` = 0.4547), tolerance = 1e-4)
  # pi = 1/2 par terme : p = BF/(1+BF), a priori 1/2 pour TOUS les termes.
  expect_equal(res$prior_h1, rep(1 / 2, 3))
  expect_equal(res$p_h1, unname(hand / (1 + hand)))
  expect_equal(res$p_h1 + res$p_h0, rep(1, 3))
  expect_equal(res$n_plus, c(2L, 2L, 1L))
  expect_equal(res$n_minus, c(2L, 2L, 1L))
  expect_equal(res$analogue, rep("II", 3))
  expect_true(attr(res, "marginal_ok"))
})

test_that("inclusion 'all' A*B : a priori uniforme et recalcul par cotes", {
  res <- .bay_inclusion(logs_2f, space_2f, "all")
  # Table du paragraphe 3 : 3/5 et 1/5.
  expect_equal(res$prior_h1, c(3 / 5, 3 / 5, 1 / 5))
  # Probabilites de modeles, puis cotes a posteriori / cotes a priori.
  post <- c(bf_null, bf_a, bf_b, bf_ab, bf_full)
  post <- post / sum(post)
  p_a <- post[2] + post[4] + post[5]
  p_b <- post[3] + post[4] + post[5]
  p_i <- post[5]
  expect_equal(res$p_h1, c(p_a, p_b, p_i))
  expect_equal(res$bf10[1], (p_a / (1 - p_a)) / ((3 / 5) / (2 / 5)))
  expect_equal(res$bf10[3], (p_i / (1 - p_i)) / ((1 / 5) / (4 / 5)))
  expect_true(all(is.na(res$analogue)))
})

test_that("inclusion 'all' avec la loi de Chipman : 1/2, 1/2, 1/8 (paragraphe 3)", {
  # pi(A) = pi(B) = 1/2 ; heredite forte : pi(A:B | A, B) = 1/2.
  chipman <- c(1/4, 1/4, 1/4, 1/8, 1/8)
  res <- .bay_inclusion(logs_2f, space_2f, "all", prior_prob = chipman)
  expect_equal(res$prior_h1, c(1 / 2, 1 / 2, 1 / 8))
  # Avec pi(A) = 1/2, BF = P/(1-P).
  post <- chipman * c(bf_null, bf_a, bf_b, bf_ab, bf_full)
  post <- post / sum(post)
  p_a <- post[2] + post[4] + post[5]
  expect_equal(res$p_h1[1], p_a)
  expect_equal(res$bf10[1], p_a / (1 - p_a))
  # BF de A:B : modele plein contre moyenne ponderee des quatre autres.
  w <- chipman[1:4]
  expect_equal(res$bf10[3], bf_full / (sum(w * c(bf_null, bf_a, bf_b, bf_ab)) / sum(w)))
})

test_that("inclusion appariee avec loi fournie : formule ponderee du paragraphe 3", {
  w <- c(0.4, 0.1, 0.2, 0.2, 0.1)
  res <- .bay_inclusion(logs_2f, space_2f, "matched", prior_prob = w)
  # A : (0,1*A + 0,2*(A+B))/0,3 contre (0,4*nul + 0,2*B)/0,6.
  hand_a <- ((0.1 * bf_a + 0.2 * bf_ab) / 0.3) / ((0.4 * bf_null + 0.2 * bf_b) / 0.6)
  expect_equal(res$bf10[1], hand_a)
  # A priori renormalise sur M+ union M- : 0,3 / (0,3 + 0,6) = 1/3.
  expect_equal(res$prior_h1[1], 1 / 3)
  # Cotes 1/2 : p = BF/2 / (1 + BF/2).
  expect_equal(res$p_h1[1], (hand_a / 2) / (1 + hand_a / 2))
})

# ---- .bay_inclusion : invariances -------------------------------------------

test_that("inclusion : invariante au modele de reference et a l'ordre", {
  for (type in c("matched", "all")) {
    ref <- .bay_inclusion(logs_2f, space_2f, type)
    # Passer du nul a A+B comme reference retire la meme constante partout.
    shifted <- .bay_inclusion(logs_2f - log(bf_ab), space_2f, type)
    expect_equal(shifted$log_bf, ref$log_bf)
    expect_equal(shifted$p_h1, ref$p_h1)
    # Enorme decalage : aucun debordement.
    huge <- .bay_inclusion(logs_2f + 5000, space_2f, type)
    expect_equal(huge$log_bf, ref$log_bf)
    perm <- c(4, 1, 5, 3, 2)
    swapped <- .bay_inclusion(logs_2f[perm], space_2f[perm], type)
    expect_equal(swapped[match(ref$term, swapped$term), "log_bf"], ref$log_bf)
    # Les poids ne comptent qu'a une constante pres.
    w <- c(1, 2, 3, 4, 5)
    expect_equal(.bay_inclusion(logs_2f, space_2f, type, prior_prob = w)$log_bf,
                 .bay_inclusion(logs_2f, space_2f, type, prior_prob = w * 7)$log_bf)
  }
})

test_that("inclusion : l'ecriture des termes ne change rien", {
  ref <- .bay_inclusion(logs_2f, space_2f, "matched")
  strings <- .bay_inclusion(logs_2f, c("1", "A", "B", "A + B", "A + B + A:B"), "matched")
  expect_equal(strings$log_bf, ref$log_bf)
  rev_space <- list(NULL, "A", "B", c("B", "A"), c("B:A", "B", "A"))
  flipped <- .bay_inclusion(logs_2f, rev_space, "matched")
  expect_equal(flipped$log_bf, ref$log_bf)
  # Le libelle affiche reste celui de l'utilisateur (memes noms que m.test).
  expect_equal(flipped$term[3], "B:A")
})

test_that("inclusion : log BF enormes, stabilite log-sum-exp", {
  # Modeles : nul 0 ; A 5000 ; B 10 ; A+B 5002 ; A+B+A:B 5003.
  logs <- c(0, 5000, 10, 5002, 5003)
  res <- .bay_inclusion(logs, space_2f, "matched")
  # A:B : 5003 - 5002 = 1.
  expect_equal(row_of(res, "A:B")$log_bf, 1)
  # A : log[(e^5000 + e^5002)/2] - log[(1 + e^10)/2]
  #   = 5002 + log(1 + e^-2) - 10 - log(1 + e^-10) (le facteur 1/2 s'elimine).
  expect_equal(row_of(res, "A")$log_bf, 5002 + log1p(exp(-2)) - 10 - log1p(exp(-10)))
  expect_identical(row_of(res, "A")$p_h1, 1)
  expect_identical(row_of(res, "A")$bf10, Inf)
  expect_equal(row_of(res, "A")$p_h0, 0)
})

test_that("inclusion : BF nuls, infinis et manquants", {
  # A+B+A:B de BF nul : BF(A:B) = 0, donc p = 0.
  res <- .bay_inclusion(c(logs_2f[1:4], -Inf), space_2f, "matched")
  expect_identical(row_of(res, "A:B")$bf10, 0)
  expect_identical(row_of(res, "A:B")$p_h1, 0)
  # Les termes qui n'utilisent pas ce modele ne changent pas.
  ref <- .bay_inclusion(logs_2f, space_2f, "matched")
  expect_equal(row_of(res, "A")$log_bf, row_of(ref, "A")$log_bf)
  # Modele a BF infini dans M+ : BF = Inf.
  res <- .bay_inclusion(c(logs_2f[1:4], Inf), space_2f, "matched")
  expect_identical(row_of(res, "A:B")$p_h1, 1)
  # Infini des deux cotes : indetermine, NA et surtout pas NaN.
  res <- .bay_inclusion(c(logs_2f[1:3], Inf, Inf), space_2f, "matched")
  expect_true(is.na(row_of(res, "A:B")$log_bf))
  expect_false(is.nan(row_of(res, "A:B")$log_bf))
  expect_true(is.na(row_of(res, "A:B")$p_h1))
  # NA dans un modele utilise : NA pour ce terme seulement, jamais de calcul sans lui.
  res <- .bay_inclusion(c(logs_2f[1:4], NA), space_2f, "matched")
  expect_true(is.na(row_of(res, "A:B")$log_bf))
  expect_equal(row_of(res, "A")$log_bf, row_of(ref, "A")$log_bf)
  res <- .bay_inclusion(c(logs_2f[1:3], NA, logs_2f[5]), space_2f, "matched")
  expect_true(all(is.na(res$log_bf)))
})

test_that("inclusion : modele de poids nul ecarte, meme a BF infini", {
  res <- .bay_inclusion(c(logs_2f[1:4], Inf), space_2f, "all", prior_prob = c(1, 1, 1, 1, 0))
  # A:B n'a plus aucune masse a priori : P(A:B) = 0, BF non defini.
  expect_true(is.na(row_of(res, "A:B")$log_bf))
  expect_identical(row_of(res, "A:B")$prior_h1, 0)
  expect_identical(row_of(res, "A:B")$p_h1, 0)
  expect_identical(row_of(res, "A:B")$p_h0, 1)
  # A : (A, A+B) contre (nul, B), le modele plein ne compte plus.
  expect_equal(row_of(res, "A")$bf10, mean(c(bf_a, bf_ab)) / mean(c(bf_null, bf_b)))
})

# ---- .bay_inclusion : marginalite et espace incomplet ------------------------

test_that("inclusion appariee : refus d'un espace qui viole la marginalite", {
  # Espace de whichModels = "all" : 8 modeles dont A:B seul, A + A:B, B + A:B.
  terms_all <- list(character(0), "A", "B", "A:B", c("A", "B"), c("A", "A:B"),
                    c("B", "A:B"), c("A", "B", "A:B"))
  logs <- log(c(1, bf_a, bf_b, 0.4048160, bf_ab, 7.9178313, 0.1298609, bf_full))
  expect_error(.bay_inclusion(logs, terms_all, "matched"), "marginal")
  # En "all" : calcul possible, mais signale clairement.
  expect_warning(res <- .bay_inclusion(logs, terms_all, "all"), "marginal", ignore.case = TRUE)
  expect_false(attr(res, "marginal_ok"))
  # Valeurs de agent_02, 8 modeles : 19,16 ; 0,3576 ; 0,4396.
  expect_equal(round(row_of(res, "A:B")$bf10, 4), 0.4396)
  expect_equal(round(row_of(res, "A")$bf10, 2), 19.16)
})

test_that("inclusion : termes de nuisance (id, covariable, pente aleatoire)", {
  # ANCOVA (5.4) : x partout. BF(A + x vs x) = BF(A+x)/BF(x).
  res <- .bay_inclusion(log(c(2, 12)), list("x", c("x", "A")), nuisance = "x")
  expect_equal(res$term, "A")
  expect_equal(res$bf10, 6)
  expect_equal(res$prior_h1, 1 / 2)
  # Mesures repetees avec id et pente A:id (5.8) : la pente ne compte pas comme
  # terme d'ordre superieur de A, et A + id + A:id contre id + A:id est valide.
  res <- .bay_inclusion(log(c(3, 15)), list(c("id", "A:id"), c("id", "A:id", "A")),
                        "matched", nuisance = c("id", "A:id"))
  expect_equal(res$term, "A")
  expect_equal(res$bf10, 5)
  # Pentes heterogenes (5.5) : x nuisance, A:x teste, son parent x est present.
  res <- .bay_inclusion(log(c(1, 4, 8)), list("x", c("x", "A"), c("x", "A", "A:x")),
                        "matched", nuisance = "x")
  expect_equal(res$bf10[res$term == "A:x"], 2)
  expect_equal(res$bf10[res$term == "A"], 4)
  # Nuisance absente d'un modele : refus.
  expect_error(.bay_inclusion(log(c(1, 2)), list(character(0), "A"), nuisance = "id"))
  # Terme partout mais non declare : refus, sinon M- serait vide.
  expect_error(.bay_inclusion(log(c(2, 12)), list("x", c("x", "A"))))
})

test_that("inclusion appariee : modele partenaire manquant, signale et compte", {
  # Sans B : la paire {A+B, B} est ecartee, A se lit sur {A} contre {nul}.
  expect_warning(
    res <- .bay_inclusion(logs_2f[-3], space_2f[-3], "matched"),
    "A"
  )
  expect_equal(row_of(res, "A")$bf10, bf_a / bf_null)
  expect_equal(row_of(res, "A")$n_unpaired, 1L)
  expect_equal(row_of(res, "A")$n_plus, 1L)
  # A:B garde sa paire, complete.
  expect_equal(row_of(res, "A:B")$n_unpaired, 0L)
})

test_that("inclusion : entrees invalides", {
  expect_error(.bay_inclusion(logs_2f[-1], space_2f))
  expect_error(.bay_inclusion(logs_2f, space_2f[-1]))
  # Pas de modele nul.
  expect_error(.bay_inclusion(logs_2f[-1], space_2f[-1]))
  # Doublon.
  expect_error(.bay_inclusion(log(c(1, 2, 3)), list(character(0), "A", "A")))
  expect_error(.bay_inclusion(log(c(1, 2, 3)), list(character(0), "A", "A*B")))
  expect_error(.bay_inclusion(logs_2f, space_2f, type = "top"))
  expect_error(.bay_inclusion(logs_2f, space_2f, prior_prob = c(1, 1)))
  expect_error(.bay_inclusion(logs_2f, space_2f, prior_prob = c(1, 1, -1, 1, 1)))
  expect_error(.bay_inclusion(logs_2f, space_2f, prior_prob = rep(0, 5)))
  expect_error(.bay_inclusion(logs_2f, space_2f, prior_prob = c(NA, 1, 1, 1, 1)))
  expect_error(.bay_inclusion(0, list(character(0))))
})

# ---- .bay_inclusion : trois facteurs -----------------------------------------

# Enumeration independante des 19 modeles marginaux a 3 facteurs. Filtre ecrit
# sur des indices, pas sur les noms de termes.
enumerate_3f <- function() {
  all_terms <- c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C")
  parents <- list(`A:B` = 1:2, `A:C` = c(1, 3), `B:C` = 2:3, `A:B:C` = 4:6)
  out <- list()
  for (code in 0:127) {
    use <- as.logical(bitwAnd(code, 2^(0:6)))
    ok <- TRUE
    for (nm in names(parents)) {
      if (use[match(nm, all_terms)] && !all(use[parents[[nm]]])) ok <- FALSE
    }
    if (ok) out[[length(out) + 1L]] <- all_terms[use]
  }
  out
}

test_that("trois facteurs : comptes et a priori du paragraphe 3", {
  space <- enumerate_3f()
  # Table : 19 modeles ; 14/19 = 0,737 ; 6/19 = 0,316 ; 1/19 = 0,053.
  expect_length(space, 19L)
  res <- .bay_inclusion(rep(0, 19), space, "all")
  expect_equal(row_of(res, "A")$prior_h1, 14 / 19)
  expect_equal(row_of(res, "A:B")$prior_h1, 6 / 19)
  expect_equal(row_of(res, "A:B:C")$prior_h1, 1 / 19)
  # BF tous egaux a 1 : pas d'evidence, BF d'inclusion = 1, p = a priori.
  expect_equal(res$bf10, rep(1, 7))
  expect_equal(res$p_h1, res$prior_h1)
})

test_that("trois facteurs : ensembles apparies et recalcul independant", {
  space <- enumerate_3f()
  set.seed(101)
  logs <- rnorm(19, sd = 3)
  logs[lengths(space) == 0] <- 0
  res <- .bay_inclusion(logs, space, "matched")
  # A la main : A dans {A, A+B, A+C, A+B+C, A+B+C+B:C}, 5 modeles ;
  # A:B dans les 5 modeles avec A:B sans A:B:C ; A:B:C dans le seul modele plein.
  expect_equal(row_of(res, "A")$n_plus, 5L)
  expect_equal(row_of(res, "A:B")$n_plus, 5L)
  expect_equal(row_of(res, "A:B:C")$n_plus, 1L)
  expect_true(all(res$n_plus == res$n_minus))
  expect_true(all(res$n_unpaired == 0L))
  expect_equal(res$prior_h1, rep(1 / 2, 7))

  # Recalcul : on retrouve les modeles par leur composition exacte.
  value <- function(...) {
    target <- c(...)
    hit <- vapply(space, function(m) setequal(m, target), logical(1))
    exp(logs[hit])
  }
  hand_a <- mean(c(value("A"), value("A", "B"), value("A", "C"), value("A", "B", "C"),
                   value("A", "B", "C", "B:C"))) /
    mean(c(value(), value("B"), value("C"), value("B", "C"), value("B", "C", "B:C")))
  expect_equal(row_of(res, "A")$bf10, hand_a)
  hand_abc <- value("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C") /
    value("A", "B", "C", "A:B", "A:C", "B:C")
  expect_equal(row_of(res, "A:B:C")$bf10, hand_abc)
  hand_ab <- mean(c(value("A", "B", "A:B"), value("A", "B", "C", "A:B"),
                    value("A", "B", "C", "A:B", "A:C"), value("A", "B", "C", "A:B", "B:C"),
                    value("A", "B", "C", "A:B", "A:C", "B:C"))) /
    mean(c(value("A", "B"), value("A", "B", "C"), value("A", "B", "C", "A:C"),
           value("A", "B", "C", "B:C"), value("A", "B", "C", "A:C", "B:C")))
  expect_equal(row_of(res, "A:B")$bf10, hand_ab)
})

# ---- .bay_sensitivity --------------------------------------------------------

test_that("sensibilite : grille r x cotes, paragraphe 8", {
  # BF = 2 a r = 0,25 ; 10 a r = 0,5 ; 5 a r = 1.
  s <- .bay_sensitivity(c(`0.25` = log(2), `0.5` = log(10), `1` = log(5)))
  expect_equal(dim(s$grid), c(3L, 3L))
  expect_equal(colnames(s$grid), c("1:4", "1:1", "4:1"))
  expect_equal(rownames(s$grid), c("0.25", "0.5", "1"))
  expect_equal(unname(round(s$grid["0.5", ], 3)), c(0.714, 0.909, 0.976))
  # Case r = 0,25 et cotes 1:4 : 0,25*2/(1+0,25*2) = 1/3.
  expect_equal(s$grid["0.25", "1:4"], 1 / 3)
  # Case r = 1 et cotes 4:1 : 20/21.
  expect_equal(s$grid["1", "4:1"], 20 / 21)
  expect_equal(s$range, c(1 / 3, 40 / 41))
  # 1/3 est sous 1/2 et 40/41 au-dessus : instable.
  expect_false(s$stable)
  expect_equal(s$side, "mixed")
})

test_that("sensibilite : drapeau de stabilite", {
  # BF >= 40 : meme a 1:4, O*BF >= 10, tous au-dessus de 1/2.
  s <- .bay_sensitivity(log(c(40, 100, 60)))
  expect_true(s$stable)
  expect_equal(s$side, "h1")
  # BF <= 1/40 : tous dessous.
  s <- .bay_sensitivity(-log(c(40, 100, 60)))
  expect_true(s$stable)
  expect_equal(s$side, "h0")
  # Exactement 1/2 sur une case (BF = 4 a 1:4) : ni d'un cote ni de l'autre.
  s <- .bay_sensitivity(log(c(4, 10, 8)))
  expect_false(s$stable)
  # NA : stabilite non prouvee, jamais NA.
  s <- .bay_sensitivity(c(log(40), NA, log(60)))
  expect_identical(s$stable, FALSE)
  expect_equal(s$n_missing, 3L)
  expect_equal(s$range, c(10 / 11, 240 / 241))
  s <- .bay_sensitivity(c(NA_real_, NA_real_))
  expect_identical(s$stable, FALSE)
  expect_true(all(is.na(s$range)))
  # Cas extremes et autres cotes.
  s <- .bay_sensitivity(c(1000, -Inf), prior_odds = 1)
  expect_equal(s$range, c(0, 1))
  expect_false(s$stable)
  s <- .bay_sensitivity(0, prior_odds = c(1/3, 2))
  expect_equal(colnames(s$grid), c("1:3", "2:1"))
  expect_equal(unname(s$grid[1, ]), c(1 / 4, 2 / 3))
  expect_error(.bay_sensitivity(0, prior_odds = c(1, 0)))
  expect_error(.bay_sensitivity(numeric(0)))
})

# ---- .bay_prior_label ---------------------------------------------------------

test_that("annonce des a priori : contenu, bilinguisme, longueur des lignes", {
  lab <- .bay_prior_label(1/2, 0.5, "fixed")
  expect_named(lab, c("en", "fr"))
  expect_match(lab$en, "1/2", fixed = TRUE)
  expect_match(lab$en, "r = 0.5", fixed = TRUE)
  expect_match(lab$fr, "r = 0,5", fixed = TRUE)
  expect_match(lab$en, "1:1", fixed = TRUE)
  expect_false(identical(lab$en, lab$fr))
  # Autre pi : cotes 1:4 pour pi = 1/5.
  lab <- .bay_prior_label(1/5, 0.707, "delta")
  expect_match(lab$en, "0.2", fixed = TRUE)
  expect_match(lab$en, "1:4", fixed = TRUE)
  expect_match(lab$en, "delta", fixed = TRUE)
  # Plusieurs termes : aucune ligne au-dela de 80 colonnes, dans le pire des cas
  # d'affichage : .vbse cpt = "sub" met deux tabulations de 8 colonnes devant.
  lab <- .bay_prior_label(1/2, c(0.5, sqrt(2) / 4, 1, 0.5), c("fixed", "covariate", "random", "delta"))
  for (txt in lab) {
    lines <- gsub("\t", "", strsplit(txt, "\n", fixed = TRUE)[[1]], fixed = TRUE)
    expect_gt(length(lines), 1L)
    expect_true(all(nchar(lines) + 16L <= 80L))
  }
  expect_match(lab$en, "0.354", fixed = TRUE)
  expect_match(lab$fr, "0,354", fixed = TRUE)
  # Sans r : seul pi est annonce.
  expect_no_match(.bay_prior_label(1/2)$en, "r =", fixed = TRUE)
  expect_error(.bay_prior_label(0, 0.5))
  expect_error(.bay_prior_label(1, 0.5))
  expect_error(.bay_prior_label(1/2, -1))
  expect_error(.bay_prior_label(1/2, 0.5, "autre"))
})

test_that("annonce des a priori : passe par .vbse et reste ASCII a la source", {
  lab <- .bay_prior_label(1/2, 0.5)
  out <- capture.output(k <- .vbse(lab$en, lab$fr, verbose = TRUE, k = 0))
  expect_equal(k, 1)
  expect_true(any(grepl("1/2", out, fixed = TRUE)))
  # Aucune reference bibliographique dans un message (bp.log, 7.4).
  expect_no_match(paste(unlist(lab), collapse = " "), "[0-9]{4}\\)|et al")
})

# ---- .bay_effects_row ---------------------------------------------------------

test_that("ligne d'effet : colonnes reutilisables, aucune p-value", {
  row <- .bay_effects_row("A", log_bf = log(10), prior_h1 = 1/2, rscale = 0.5,
                          error_pct = 1.2, analogue = "II")
  expect_named(row, c("term", "bf10", "log_bf", "p_h1", "p_h0", "prior_h1",
                      "rscale", "error_pct", "analogue"))
  expect_false(any(grepl("p.?val", names(row), ignore.case = TRUE)))
  expect_false(any(names(row) %in% c("p", "P", "pval", "p.value", "Pr")))
  expect_equal(row$bf10, 10)
  expect_equal(row$p_h1, 10 / 11)
  expect_equal(row$p_h0, 1 / 11)
  expect_equal(row$rscale, 0.5)
  expect_equal(row$error_pct, 1.2)
  expect_equal(row$analogue, "II")
})

test_that("ligne d'effet : bf10 ou log_bf, vectorisation, a priori", {
  expect_equal(.bay_effects_row("A", bf10 = 4), .bay_effects_row("A", log_bf = log(4)))
  rows <- .bay_effects_row(c("A", "B"), log_bf = c(0, log(3)), prior_h1 = c(1/2, 1/4))
  expect_equal(nrow(rows), 2L)
  # B : cotes 1/3 et BF = 3, donc cotes a posteriori 1 et p = 1/2.
  expect_equal(rows$p_h1, c(1 / 2, 1 / 2))
  expect_equal(rows$p_h1 + rows$p_h0, c(1, 1))
  # Les petites queues survivent ; bf10 deborde, log_bf reste juste.
  row <- .bay_effects_row("A", log_bf = 1000)
  expect_identical(row$bf10, Inf)
  expect_equal(row$log_bf, 1000)
  expect_equal(log(.bay_effects_row("A", log_bf = 700)$p_h0), -700)
  # NA.
  row <- .bay_effects_row("A", log_bf = NA)
  expect_true(is.na(row$p_h1) && is.na(row$p_h0) && is.na(row$bf10))
  # pi = 0 ou 1 : les donnees ne peuvent rien y changer.
  expect_identical(.bay_effects_row("A", log_bf = 5, prior_h1 = 0)$p_h1, 0)
  expect_identical(.bay_effects_row("A", log_bf = -5, prior_h1 = 1)$p_h0, 0)
  expect_error(.bay_effects_row("A"))
  expect_error(.bay_effects_row("A", log_bf = 0, bf10 = 1))
  expect_error(.bay_effects_row("A", bf10 = -1))
  expect_error(.bay_effects_row(c("A", "B"), log_bf = 0))
  expect_error(.bay_effects_row("A", log_bf = 0, prior_h1 = 1.5))
  expect_error(.bay_effects_row("A", log_bf = 0, rscale = -1))
  expect_error(.bay_effects_row("A", log_bf = 0, error_pct = -1))
})

test_that(".bay_inclusion partage les colonnes de .bay_effects_row", {
  res <- .bay_inclusion(logs_2f, space_2f)
  expect_equal(names(res)[1:9], names(.bay_effects_row("A", log_bf = 0)))
  expect_false(any(grepl("p.?val", names(res), ignore.case = TRUE)))
})
