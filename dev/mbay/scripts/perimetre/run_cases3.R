options(warn = 1, width = 120)
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE))
set.seed(5)
dsp <- expand.grid(id = 1:24, T = c("t1","t2","t3")); dsp$G <- factor(ifelse(dsp$id <= 12, "g1", "g2")); dsp$id <- factor(dsp$id)
dsp$y <- rep(rnorm(24, 0, 2), 3) + as.numeric(dsp$T) + (dsp$G == "g2") + rnorm(72, 0, 0.7)
set.seed(8)
dhet <- data.frame(F = factor(rep(c("A","B","C"), each = 30))); dhet$y <- rnorm(90, 10, c(0.5, 2, 5)[as.numeric(dhet$F)])
set.seed(4); n <- 20
drm <- data.frame(id = factor(rep(seq_len(n), times = 3)), F = factor(rep(c("t1","t2","t3"), each = n)))
drm$y <- rep(rnorm(n, 0, 2), 3) + as.numeric(drm$F) + rnorm(3 * n, 0, 0.7)
trace_err <- function(label, expr) {
  cat("\n## ", label, "\n")
  tb <- NULL; err <- NULL; res <- NULL
  invisible(capture.output(res <- tryCatch(withCallingHandlers(suppressMessages(suppressWarnings(eval(expr))),
    error = function(e) tb <<- sys.calls()), error = function(e) err <<- conditionMessage(e))))
  if (!is.null(err)) {
    cat("ERREUR:", err, "\n")
    fn <- vapply(tb, function(cl) substr(deparse(cl)[1], 1, 70), "")
    cat(paste(" ", tail(fn[!grepl("tryCatch|withCalling|doTry|^eval|signal|^\\.handle|^h\\(", fn)], 6)), sep = "\n")
  } else cat("OK method=", res$method, "p=", res$global_pvalue, "\n")
}
vb <- function(label, expr, pat = NULL, n = 45) {
  cat("\n## ", label, "\n")
  o <- capture.output(r <- suppressMessages(suppressWarnings(eval(expr))))
  o <- o[nzchar(trimws(o))]
  if (!is.null(pat)) o <- grep(pat, o, value = TRUE)
  cat(substr(head(o, n), 1, 135), sep = "\n")
  cat("=> method:", r$method, " p:", r$global_pvalue, "\n")
  invisible(r)
}
trace_err("18 verbose=FALSE boot defaut", quote(m.test(y ~ G * T + Error(id/T), data = dsp, verbose = FALSE, plot = FALSE)))
trace_err("18 verbose=TRUE boot defaut", quote(m.test(y ~ G * T + Error(id/T), data = dsp, verbose = TRUE, plot = FALSE)))
trace_err("18 verbose=FALSE boot=FALSE", quote(m.test(y ~ G * T + Error(id/T), data = dsp, verbose = FALSE, plot = FALSE, boot = FALSE)))
vb("31 heteroscedastique", quote(m.test(y ~ F, data = dhet, plot = FALSE, boot = FALSE)), "variance|Variance|Welch|ANOVA|normal|Bartlett|Levene|Brown", 18)
vb("16 id+wt", quote(m.test(y ~ G, data = dsp, id = "id", wt = "T", plot = FALSE, boot = FALSE)), NULL, 22)
vb("17 within", quote(m.test(y ~ T, data = droplevels(dsp[dsp$G == "g1", ]), id = "id", within = "T", plot = FALSE, boot = FALSE)), "Friedman|Formule|within|intra|robuste|Mauchly|Wilcoxon|Nemenyi", 14)
vb("6a RM Error(id/F)", quote(m.test(y ~ F + Error(id/F), data = drm, plot = FALSE, boot = FALSE)), "Mauchly|Greenhouse|Formule|aov|ANOVA|Type|apparie|appari|Student", 14)
r <- m.test(y ~ F + (1|id), data = drm, verbose = FALSE, plot = FALSE, boot = FALSE)
cat("6b p.value (1|id) :\n"); print(r$p.value)
cat("reference pairwise.t.test non apparie Holm:\n"); print(pairwise.t.test(drm$y, drm$F, p.adjust.method = "holm", pool.sd = FALSE)$p.value)
cat("reference aov Error(id):", summary(aov(y ~ F + Error(id), data = drm))[[2]][[1]][1, "Pr(>F)"], "\n")
