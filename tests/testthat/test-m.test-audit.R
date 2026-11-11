# Audit de m.test() (novembre 2026) : un test par defaut corrige.
# Chaque test echoue avant la correction et passe apres ; voir dev/mtest_audit/RAPPORT.md.

audit_quiet <- function(expr) {
  value <- NULL
  output <- utils::capture.output(value <- suppressMessages(suppressWarnings(expr)))
  list(value = value, output = output)
}

make_split_plot <- function() {
  set.seed(11)
  n <- 12
  ids <- factor(1:n)
  G <- factor(rep(c("A", "B"), each = n / 2))
  d <- expand.grid(id = ids, T = factor(c("t1", "t2", "t3")))
  d$G <- G[as.integer(d$id)]
  u <- rnorm(n, 0, 1)
  d$y <- 10 + u[as.integer(d$id)] + ifelse(d$G == "B" & d$T == "t3", 1.2, 0) +
    rnorm(nrow(d), 0, 0.8)
  d
}

# p-value d'un terme intra-sujet corrigee de la sphericite, regle de Girden (1992),
# calculee par car::Anova() sur le format large, sans passer par m.test()
girden_ref <- function(d, term) {
  w <- reshape(d[, c("id", "G", "T", "y")], idvar = c("id", "G"), timevar = "T", direction = "wide")
  ycols <- grep("^y\\.", names(w), value = TRUE)
  idata <- data.frame(T = factor(sub("^y\\.", "", ycols)))
  fit <- lm(as.matrix(w[, ycols]) ~ G, data = w, contrasts = list(G = "contr.sum"))
  sm <- summary(car::Anova(fit, idata = idata, idesign = ~T, type = 3), multivariate = FALSE)
  adj <- sm$pval.adjustments
  if (adj[term, "GG eps"] <= 0.75) adj[term, "Pr(>F[GG])"] else {
    uni <- sm$univariate.tests
    eps <- min(1, adj[term, "HF eps"])
    pf(uni[term, "F value"], eps * uni[term, "num Df"], eps * uni[term, "den Df"], lower.tail = FALSE)
  }
}

test_that("a. split-plot : l'interaction G:T est gardee, avec les p-values d'aov", {
  d <- make_split_plot()
  ref <- summary(aov(y ~ G * T + Error(id / T), data = d))
  p_G <- ref[["Error: id"]][[1]]["G", "Pr(>F)"]
  res <- audit_quiet(m.test(y ~ G * T + Error(id / T), data = d, plot = FALSE))
  eff <- res$value$effects
  expect_true("G:T" %in% eff$term)
  # G:T est intra-sujet : p-value d'aov corrigee de la sphericite (correctif l)
  expect_equal(eff$p.value[eff$term == "G:T"], girden_ref(d, "G:T"), tolerance = 1e-6)
  expect_equal(eff$p.value[eff$term == "G"], p_G, tolerance = 1e-6)
  # G est inter-sujet : il ne doit pas etre annonce comme intra-sujet
  expect_false(any(grepl("intra-sujet", res$output) & grepl("(s) : G", res$output, fixed = TRUE)))
  expect_false(any(grepl("within-subject factor(s): G", res$output, fixed = TRUE)))
})

test_that("b. split-plot avec verbose = FALSE : pas de plantage", {
  d <- make_split_plot()
  res <- audit_quiet(m.test(y ~ G * T + Error(id / T), data = d, plot = FALSE, verbose = FALSE))
  expect_true("G:T" %in% res$value$effects$term)
})

test_that("b. modele mixte force (plan desequilibre) avec verbose = FALSE : pas de plantage", {
  d <- make_split_plot()
  d <- d[-c(1, 20), ]
  expect_no_error(audit_quiet(m.test(y ~ G * T + Error(id / T), data = d,
                                     plot = FALSE, verbose = FALSE)))
  res <- audit_quiet(m.test(y ~ G * T + Error(id / T), data = d, plot = FALSE))
  expect_true("G:T" %in% res$value$effects$term)
})

test_that("f. un message anglais d'une seule chaine s'affiche en locale anglaise", {
  old_ctype <- Sys.getlocale("LC_CTYPE")
  withr::defer(Sys.setlocale("LC_CTYPE", old_ctype))
  skip_if(identical(Sys.setlocale("LC_CTYPE", "C"), ""), "locale C indisponible")
  expect_identical(.msg("en", "fr"), "en")
  out <- utils::capture.output(.vbse("Hello", "Bonjour", verbose = TRUE, k = 0))
  expect_match(out[1], "Hello")
  expect_error(.exit("Stop here", "Arret ici"), "Stop here")
  # ang absent ou vide : on retombe sur le francais
  out <- utils::capture.output(.vbse("", "Bonjour", verbose = TRUE, k = 0))
  expect_match(out[1], "Bonjour")
  # bout en bout : m.test() parle anglais
  set.seed(3)
  res <- audit_quiet(m.test(rnorm(20), rep(c("A", "B"), 10), plot = FALSE))
  expect_true(any(grepl("Fisher-Snedecor test", res$output, fixed = TRUE)))
  expect_false(any(grepl("Test de Fisher-Snedecor", res$output, fixed = TRUE)))
})

test_that("e. m.test(formula = y ~ F, data = d) est accepte, comme m.test(y ~ F, d)", {
  set.seed(5)
  d <- data.frame(F = rep(c("A", "B", "C"), each = 10), y = rnorm(30))
  ref <- audit_quiet(m.test(y ~ F, data = d, return = FALSE))$value
  p_named <- audit_quiet(m.test(formula = y ~ F, data = d, return = FALSE))$value
  expect_equal(p_named, ref)
  f <- y ~ F
  expect_equal(audit_quiet(m.test(formula = f, data = d, return = FALSE))$value, ref)
  expect_equal(audit_quiet(m.test(formula = "y ~ F", data = d, return = FALSE))$value, ref)
})

test_that("c. la confirmation de Levene porte sur les variances (Brown-Forsythe = mediane)", {
  set.seed(1)
  n <- c(5, 5, 40)
  x <- unlist(mapply(function(n, s) rnorm(n, 0, s), n, c(4, 4, 1)))
  g <- factor(rep(c("A", "B", "C"), n))
  p_mean <- car::leveneTest(x, g, center = mean)[1, 3]
  p_median <- car::leveneTest(x, g, center = median)[1, 3]
  expect_lt(p_mean, 0.05)
  expect_lt(p_median, 0.05)
  out <- utils::capture.output(v <- .variance(x, g, check_normality = FALSE, verbose = TRUE, k = 0))
  # Les deux tests de variances concluent a l'heterogeneite : .variance() doit le suivre
  expect_false(v[[1]])
  expect_true(any(grepl(.format_pval(p_median), out, fixed = TRUE)))
})

test_that("d. (1|id) donne des post-hocs apparies, comme Error(id/T)", {
  set.seed(7)
  n <- 15
  d <- expand.grid(id = factor(1:n), T = factor(c("t1", "t2", "t3")))
  u <- rnorm(n, 0, 3)
  d$y <- 20 + u[as.integer(d$id)] + c(0, 0.8, 1.2)[as.integer(d$T)] + rnorm(nrow(d), 0, 0.7)
  r_lmer <- audit_quiet(m.test(y ~ T + (1 | id), data = d, plot = FALSE, verbose = FALSE))$value
  r_err <- audit_quiet(m.test(y ~ T + Error(id / T), data = d, plot = FALSE, verbose = FALSE))$value
  expect_equal(r_lmer$p.value, r_err$p.value)
  # t1-t2 : apparie p ~ 0,02 ; non apparie (ancien comportement) p = 1
  expect_lt(r_lmer$p.value["t2", "t1"], 0.1)
})

test_that("g. la MANOVA non normale aboutit vite, par le Wilks sur les rangs", {
  set.seed(2)
  d <- data.frame(g = rep(c("A", "B", "C"), each = 30), a = rexp(90), b = rexp(90)^2)
  d$a[1:3] <- 40
  duree <- system.time(r <- audit_quiet(m.test(cbind(a, b) ~ g, data = d, plot = FALSE,
                                                 verbose = FALSE))$value)[["elapsed"]]
  expect_identical(r$method, "Robust MANOVA (Wilks, ranks)")
  ref <- rrcov::Wilks.test(cbind(d$a, d$b), grouping = factor(d$g), method = "rank")$p.value
  expect_equal(r$global_pvalue, unname(ref))
  expect_lt(duree, 120)
})

test_that("h. mesures repetees : lignes non triees par sujet, post-hocs toujours apparies", {
  set.seed(7)
  n <- 15
  d <- expand.grid(id = factor(1:n), T = factor(c("t1", "t2", "t3")))
  u <- rnorm(n, 0, 3)
  d$y <- 20 + u[as.integer(d$id)] + c(0, 0.8, 1.2)[as.integer(d$T)] + rnorm(nrow(d), 0, 0.7)
  set.seed(1)
  ds <- d[sample(nrow(d)), ]
  ref <- stats::pairwise.t.test(d$y, d$T, paired = TRUE)$p.value
  r_sorted <- audit_quiet(m.test(y ~ T + Error(id / T), data = d, plot = FALSE, verbose = FALSE))$value
  r_shuffled <- audit_quiet(m.test(y ~ T + Error(id / T), data = ds, plot = FALSE, verbose = FALSE))$value
  expect_equal(r_sorted$p.value, ref)
  expect_equal(r_shuffled$p.value, ref)
  # appel par vecteurs avec id, plus de 2 modalites : pas de plantage, memes post-hocs
  r_vec <- audit_quiet(m.test(ds$y, ds$T, id = ds$id, plot = FALSE, verbose = FALSE))$value
  expect_equal(unname(r_vec$p.value), unname(ref))
})

test_that("i. split-plot : post-hocs par emmeans, sans apparier des sujets differents", {
  d <- make_split_plot()
  res <- audit_quiet(m.test(y ~ G * T + Error(id / T), data = d, plot = FALSE, boot = FALSE))
  expect_false(any(grepl("pairwise.t.test(paired=TRUE)", res$output, fixed = TRUE)))
  expect_true(any(grepl("emmeans", res$output, fixed = TRUE)))
  mm <- lmerTest::lmer(y ~ G * T + (1 | id), data = d)
  ref <- as.data.frame(summary(pairs(emmeans::emmeans(mm, ~ G:T), adjust = "tukey")))
  p_ref <- ref$p.value[ref$contrast == "B t1 - B t3"]
  ligne <- grep("B t1 - B t3 : diff", res$output, value = TRUE)
  expect_length(ligne, 1)
  expect_match(ligne, paste0("p = ", .format_pval(p_ref)), fixed = TRUE)
})

test_that("j. lettres d'une interaction en modele mixte : la cellule entiere est nommee", {
  d <- make_split_plot()
  res <- audit_quiet(m.test(y ~ G * T + Error(id / T), data = d, plot = FALSE, boot = FALSE))
  expect_true(any(grepl("B t3 : EMM = ", res$output, fixed = TRUE)))
})

test_that("k. non normal et variances heterogenes confirmees : Welch, pas Kruskal-Wallis", {
  set.seed(1)
  n <- c(5, 5, 40)
  x <- unlist(mapply(function(n, s) rnorm(n, 0, s), n, c(4, 4, 1)))
  g <- factor(rep(c("A", "B", "C"), n))
  res <- audit_quiet(m.test(x, g, plot = FALSE, boot = FALSE))
  expect_identical(res$value$method, "Welch ANOVA")
  expect_equal(res$value$global_pvalue, oneway.test(x ~ g)$p.value)
  expect_equal(audit_quiet(m.test(x, g, return = FALSE, boot = FALSE))$value, oneway.test(x ~ g)$p.value)
})

make_ar1 <- function(seed, n = 12, k = 6, rho = 0.9) {
  set.seed(seed)
  sig <- rho^abs(outer(1:k, 1:k, "-"))
  y <- matrix(rnorm(n * k), n, k) %*% chol(sig)
  d <- data.frame(id = factor(rep(1:n, k)), T = factor(rep(paste0("t", 1:k), each = n)),
                  y = as.vector(y))
  list(d = d, wide = y)
}

test_that("l. mesures repetees : p-value toujours corrigee de la sphericite (GG ou HF)", {
  ar <- make_ar1(3)
  d <- ar$d
  k <- ncol(ar$wide); n <- nrow(ar$wide)
  # Calcul independant : epsilon de Greenhouse-Geisser sur contrastes orthonormes
  C <- stats::contr.poly(k)
  S <- t(C) %*% stats::cov(ar$wide) %*% C
  eps_gg <- sum(diag(S))^2 / ((k - 1) * sum(diag(S %*% S)))
  aov_tab <- summary(aov(y ~ T + Error(id / T), data = d))[["Error: id:T"]][[1]]
  f <- aov_tab["T", "F value"]
  # Huynh-Feldt (Lecoutre) pour la regle de Girden
  eps_hf <- min(1, ((n - 1) * (k - 1) * eps_gg - 2) / ((k - 1) * (n - 1 - (k - 1) * eps_gg)))
  eps <- if (eps_gg <= 0.75) eps_gg else eps_hf
  p_ref <- stats::pf(f, eps * (k - 1), eps * (k - 1) * (n - 1), lower.tail = FALSE)
  res <- audit_quiet(m.test(y ~ T + Error(id / T), data = d, plot = FALSE, boot = FALSE))
  expect_equal(res$value$effects$p.value[res$value$effects$term == "T"], p_ref, tolerance = 1e-6)
  expect_equal(res$value$global_pvalue, p_ref, tolerance = 1e-6)
  expect_true(any(grepl("Greenhouse-Geisser|Huynh-Feldt", res$output)))
  expect_false(any(grepl("Aucune correction n.cessaire", res$output)))
})

test_that("m. le recit du Welch heteroscedastique ne promet pas la robustesse", {
  set.seed(1)
  n <- c(5, 5, 40)
  x <- unlist(mapply(function(n, s) rnorm(n, 0, s), n, c(4, 4, 1)))
  g <- factor(rep(c("A", "B", "C"), n))
  res <- audit_quiet(m.test(x, g, plot = FALSE, boot = FALSE))
  expect_true(any(grepl("Prudence", res$output, fixed = TRUE)))
  expect_false(any(grepl("robuste \u00e0 une non-normalit\u00e9 mod\u00e9r\u00e9e", res$output, fixed = TRUE)))
})
