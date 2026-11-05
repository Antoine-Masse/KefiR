options(warn = 1, width = 120)
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE))

KEYS <- "\\[[a-zA-Z_.0-9:]+\\(|t2way|t3way|Friedman|Scheirer|aovp|lmer|Kruskal|Wilcoxon|Welch|med1way|t1way|rlm|HC3|Pillai|Wilks|Formule :|supprim|Sommes des Carr|Type I|ROBUSTE|robuste|ATTENTION|Attention|filtr|cart"

rt <- function(label, expr) {
  cat("\n\n## CAS :", label, "\n")
  err <- NULL; res <- NULL; tb <- NULL
  out <- capture.output({
    res <- tryCatch(withCallingHandlers(
      suppressMessages(eval(expr, envir = parent.frame())),
      warning = function(w) invokeRestart("muffleWarning"),
      error = function(e) { tb <<- sys.calls() }),
      error = function(e) { err <<- conditionMessage(e); NULL })
  })
  if (!is.null(err)) {
    cat("ERREUR:", substr(err, 1, 200), "\n")
    if (!is.null(tb)) {
      fn <- vapply(tb, function(cl) substr(deparse(cl[[1]])[1], 1, 40), "")
      cat("pile:", paste(tail(fn[!grepl("tryCatch|withCalling|doTry|eval|signal|stop$|^function", fn)], 8), collapse = " > "), "\n")
    }
  } else {
    cat("method:", paste(res$method, collapse = "|"), " chosen_test:", paste(res$chosen_test, collapse = "|"), " p_glob:", format(res$global_pvalue, digits = 4), "\n")
    if (!is.null(res$groups)) cat("colonnes groups:", paste(names(res$groups), collapse = ", "), "\n")
    if (!is.null(res$effects)) { cat("effects:\n"); print(res$effects, digits = 3) }
  }
  hits <- grep(KEYS, out, value = TRUE)
  hits <- unique(trimws(substr(hits, 1, 130)))
  cat("trace verbose (extraits):\n"); cat(paste0("  ", head(hits, 14)), sep = "\n")
  invisible(res)
}

A <- function(f, d, ...) m.test(f, data = d, verbose = TRUE, plot = FALSE, boot = FALSE, ...)

set.seed(21)
d2n <- data.frame(F = factor(rep(c("A","B"), each = 40))); d2n$y <- rexp(80, c(1, .3)[as.numeric(d2n$F)])
set.seed(22)
d2d <- data.frame(F = factor(rep(c("A","B"), each = 60))); d2d$y <- rpois(120, 1)
set.seed(23)
d3out <- data.frame(F = factor(rep(c("A","B","C"), each = 30))); d3out$y <- rnorm(90); d3out$y[1:4] <- 40
set.seed(24)
dabn <- expand.grid(A = c("a1","a2"), B = c("b1","b2"), r = 1:25); dabn$y <- rexp(100, 1 + (dabn$A == "a2")); dabn$y[1:3] <- 60
set.seed(25)
d3f <- expand.grid(A = c("a1","a2"), B = c("b1","b2"), C = c("c1","c2"), r = 1:12); d3f$y <- rexp(nrow(d3f))^2
set.seed(26)
dsev <- expand.grid(A = c("a1","a2"), B = c("b1","b2"), r = 1:30); dsev$y <- rnorm(120)
dsev <- dsev[-which(dsev$A == "a1" & dsev$B == "b1")[1:20], ]
set.seed(27)
dext <- expand.grid(A = c("a1","a2"), B = c("b1","b2"), r = 1:44); dext$y <- rnorm(nrow(dext))
dext <- dext[-which(dext$A == "a1" & dext$B == "b1")[1:40], ]
set.seed(28)
drmn <- data.frame(id = factor(rep(1:20, 3)), F = factor(rep(c("t1","t2","t3"), each = 20))); drmn$y <- rexp(60)^2
set.seed(29)
dp <- data.frame(id = factor(rep(1:15, 2)), F = factor(rep(c("t1","t2"), each = 15))); dp$y <- rnorm(30)
dpn <- dp; set.seed(30); dpn$y <- rexp(30)^3
set.seed(31)
dcn <- data.frame(F = factor(rep(c("A","B","C"), each = 30)), cov = runif(90)); dcn$y <- rexp(90)^2 + dcn$cov
set.seed(32)
dc2 <- expand.grid(F = c("a","b"), G = c("u","v"), r = 1:20); dc2$cov <- runif(80); dc2$y <- rexp(80)^2 + dc2$cov
set.seed(33)
dman <- data.frame(F = factor(rep(c("A","B","C"), each = 25))); dman$y1 <- rexp(75)^2; dman$y2 <- rexp(75)^2; dman$cov <- runif(75)
set.seed(34)
dsp <- expand.grid(id = 1:24, T = c("t1","t2","t3")); dsp$G <- factor(ifelse(dsp$id <= 12, "g1", "g2")); dsp$id <- factor(dsp$id)
dsp$y <- rep(rnorm(24, 0, 2), 3) + as.numeric(dsp$T) + (dsp$G == "g2") + rnorm(72, 0, 0.7)
set.seed(35)
dcorr <- data.frame(F = factor(rep(c("A","B"), each = 40))); dcorr$cov <- rnorm(80, ifelse(dcorr$F == "A", 0, 3)); dcorr$y <- dcorr$cov + rnorm(80)
set.seed(36)
dtiny <- data.frame(F = factor(rep(c("A","B","C"), each = 3))); dtiny$y <- rnorm(9)
set.seed(37)
dab2 <- expand.grid(A = c("a1","a2"), B = c("b1","b2"), r = 1:15); dab2$y <- rnorm(60)

rt("A. x,g data.frame sans data (reproduit cas 8)", quote(m.test(dab2$y, dab2[, c("A","B")], verbose = TRUE, plot = FALSE, boot = FALSE)))
rt("A2. x='y', g=c('A','B'), data", quote(m.test("y", c("A","B"), data = dab2, verbose = TRUE, plot = FALSE, boot = FALSE)))
rt("B. y~G*T+Error(id/T) split-plot (reproduit cas 18)", quote(A(y ~ G * T + Error(id/T), dsp)))
rt("B2. y~G*T+(1|id)", quote(A(y ~ G * T + (1|id), dsp)))
rt("B3. y~G*T+Error(id) ", quote(A(y ~ G * T + Error(id), dsp)))
rt("B4. y~G+T, id=id", quote(A(y ~ G + T, dsp, id = "id")))
rt("C. 2 groupes non normaux", quote(A(y ~ F, d2n)))
rt("D. 2 groupes discrets", quote(A(y ~ F, d2d)))
rt("E. 3 groupes outliers extremes", quote(A(y ~ F, d3out)))
rt("F. A*B non normal + outliers", quote(A(y ~ A * B, dabn)))
rt("G. A*B*C non normal", quote(A(y ~ A * B * C, d3f)))
rt("H. A*B desequilibre ratio 3", quote(A(y ~ A * B, dsev)))
rt("I. A*B desequilibre ratio 11", quote(A(y ~ A * B, dext)))
rt("J. RM 3 temps non normal, id", quote(A(y ~ F, drmn, id = "id")))
rt("K. paire normal, id", quote(A(y ~ F, dp, id = "id")))
rt("L. paire non normal, id", quote(A(y ~ F, dpn, id = "id")))
rt("M. ANCOVA 1 facteur non normale", quote(A(y ~ F + cov, dcn)))
rt("N. ANCOVA 2 facteurs non normale", quote(A(y ~ F + G + cov, dc2)))
rt("O. ANCOVA covariable correlee au facteur", quote(A(y ~ F + cov, dcorr)))
rt("P. ANCOVA avec id", quote(A(y ~ F + cov, cbind(dcn, id = factor(rep(1:30, 3))), id = "id")))
rt("Q. MANOVA non normale", quote(A(cbind(y1, y2) ~ F, dman)))
rt("R. MANCOVA", quote(A(cbind(y1, y2) ~ F + cov, dman)))
rt("S. n=3 par groupe", quote(A(y ~ F, dtiny)))
rt("T. x,g,id vecteur, 3 temps", quote(m.test(drmn$y, drmn$F, id = drmn$id, verbose = TRUE, plot = FALSE, boot = FALSE)))
rt("U. formula= nomme", quote(m.test(formula = y ~ F, data = dp, verbose = FALSE, plot = FALSE, boot = FALSE)))
rt("V. control", quote(m.test(y ~ F, data = d3out, control = "A", verbose = FALSE, plot = FALSE, boot = FALSE)))

cat("\n\n## Fonctions internes appelees sur (1|id)\n")
print(KefiR:::.drop_error_term(y ~ F + (1|id)))
print(KefiR:::.drop_error_term(y ~ G * T + Error(id/T)))
print(KefiR:::.formulator_safe(y ~ G * T + Error(id/T), mode = "linear"))
print(KefiR:::.formulator_safe(y ~ G * T + (1|id), mode = "linear"))
print(KefiR:::.detect_model_type(y ~ F + cov, dcn))
print(class(KefiR:::.detect_model_type(y ~ F + cov, dcn)))
print(try(KefiR:::.detect_model_type(y ~ c1, data.frame(y = 1:5, c1 = rnorm(5))), silent = TRUE))
nf <- KefiR:::.normalize_from_formula(y ~ G * T + Error(id/T), dsp, id = "id")
str(nf, max.level = 1, give.attr = FALSE)
print(identical(environment(KefiR:::.normalize_from_formula), asNamespace("KefiR")))
print(names(formals(KefiR:::.normalize_from_formula)))
print(KefiR:::.normalize_from_formula(y ~ T, dsp)$predictors)
print(try(names(KefiR:::.normalize_from_formula(y ~ T + (1|id), dsp)$mf), silent = TRUE))
print(try(KefiR:::.select_ss_type(y ~ A * B, dab2, model = lm(y ~ A * B, dab2), analysis_type = "ANOVA")$ss_type))
print(try(KefiR:::.select_ss_type(y ~ A * B, dsev, model = lm(y ~ A * B, dsev), analysis_type = "ANOVA")$ss_type))
print(try(KefiR:::.select_ss_type(y ~ A + B + C, d3f, model = lm(y ~ A + B + C, d3f), analysis_type = "ANOVA")$ss_type))
print(try(KefiR:::.select_ss_type(y ~ F + cov, dcn, model = lm(y ~ F + cov, dcn), analysis_type = "ANCOVA")$ss_type))
print(discret.test(rpois(120, 2))); print(discret.test(rnorm(100)))
