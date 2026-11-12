.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_load_bay.R"); .kefir_bay()
suppressMessages(library(BayesFactor))
Sys.setenv(LANG = "en")
try1 <- function(lab, expr) {
	cat("\n=====", lab, "=====\n")
	r <- tryCatch(expr, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
	invisible(r)
}
show <- function(r) if (!is.null(r)) { print(r$effects[, c("term", "p_h1", "prior_h1", "log_bf")]); cat("p_null", r$p_null <- r$hypotheses$posterior[r$hypotheses$null], " n_models", r$n_models, "\n") }

## 1. Prior HOP by hand, 2 and 3 factors
sp <- .bay_fact_space(c("A", "B", "A:B"))
print(cbind(model = vapply(sp$models, function(m) paste(sp$keys[m], collapse = "+"), ""), p = round(exp(.bay_fact_prior(sp, "hop")), 5)))
sp3 <- .bay_fact_space(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
pm <- exp(.bay_fact_prior(sp3, "hop"))
# by hand: order 1 k=3 Beta(1,1): P(set of size s) = s!(3-s)!/4!
hand <- vapply(sp3$models, function(m) {
	k <- sp3$keys[m]; s1 <- sum(!grepl(":", k)); s2 <- sum(lengths(strsplit(k, ":")) == 2); s3 <- sum(lengths(strsplit(k, ":")) == 3)
	main <- sp3$keys[1:3][sp3$keys[1:3] %in% k]
	e2 <- choose(length(main), 2)
	e3 <- as.integer(s2 == 3)
	beta(1 + s1, 1 + 3 - s1) * beta(1 + s2, 1 + e2 - s2) * beta(1 + s3, 1 + e3 - s3)
}, 0)
cat("max |HOP - hand| 3 factors:", max(abs(pm - hand)), " sum:", sum(pm), "\n")
inc <- .bay_fact_terms(sp3, rep(0, 19), log(pm)); print(inc[, c("term", "prior_h1")])

## 2. Main dataset 2x3, compare to anovaBF whichModels all
set.seed(8)
d2 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 1:6)
d2$y <- rnorm(nrow(d2), 1.5 * (d2$A == "a2"))
r <- try1("base y~A*B", m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE)); show(r)
print(r$models)
set.seed(3)
ab <- anovaBF(y ~ A * B, data = d2, whichModels = "all", rscaleFixed = 0.5, iterations = 1e5, progress = FALSE)
print(extractBF(ab, logbf = TRUE)[, 1:2])

## 3. graine
set.seed(123); s0 <- .Random.seed
r2 <- m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE)
cat("seed intact:", identical(s0, .Random.seed), " repro:", identical(r$effects$p_h1, r2$effects$p_h1), "\n")
rm(.Random.seed, envir = globalenv())
r3 <- m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE)
cat("no seed before -> exists after:", exists(".Random.seed", envir = globalenv()), "\n")
rs <- m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = TRUE)
cat("sens TRUE vs FALSE central p_h1:\n"); print(rbind(r$effects$p_h1, rs$effects$p_h1)); print(rbind(r$models$log_bf, rs$models$log_bf))
r4 <- m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE, draws = 2000)
r5 <- m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = FALSE, draws = 9000)
cat("draws 2000 vs 4000 vs 9000 identical:", identical(r4$effects$p_h1, r$effects$p_h1), identical(r5$effects$p_h1, r$effects$p_h1), "\n")
cat("error_pct:", r$effects$error_pct[1], "\n")

## 4. Order of terms, other spellings
rb <- try1("B*A", m.bay(y ~ B * A, data = d2, verbose = FALSE, sensitivity = FALSE)); show(rb)
rc <- try1("A+B+A:B", m.bay(y ~ A + B + A:B, data = d2, verbose = FALSE, sensitivity = FALSE)); show(rc)
rd <- try1("A+B+B:A", m.bay(y ~ A + B + B:A, data = d2, verbose = FALSE, sensitivity = FALSE)); show(rd)
re <- try1("x,g df", m.bay(d2$y, d2[c("A", "B")], verbose = FALSE, sensitivity = FALSE)); show(re); if (!is.null(re)) cat(re$label, "\n")
rf <- try1("A:B only", m.bay(y ~ A:B, data = d2, verbose = FALSE, sensitivity = FALSE))
rg <- try1("A + A:B", m.bay(y ~ A + A:B, data = d2, verbose = FALSE, sensitivity = FALSE))
rh <- try1("A/B", m.bay(y ~ A / B, data = d2, verbose = FALSE, sensitivity = FALSE))

## 5. Names with spaces / special characters
d3 <- d2; names(d3)[1:2] <- c("my A", "B-2")
ri <- try1("backtick additive", m.bay(y ~ `my A` + `B-2`, data = d3, verbose = FALSE, sensitivity = FALSE)); show(ri)
rj <- try1("backtick interaction", m.bay(y ~ `my A` * `B-2`, data = d3, verbose = FALSE, sensitivity = FALSE)); show(rj)
rk <- try1("x,g df with odd names", m.bay(d3$y, d3[c("my A", "B-2")], verbose = FALSE, sensitivity = FALSE)); show(rk)
d4 <- d2; names(d4)[1:2] <- c("y", "B"); names(d4)[4] <- "score"
rl <- try1("factor named y", m.bay(score ~ y * B, data = d4, verbose = FALSE, sensitivity = FALSE)); show(rl); if (!is.null(rl)) cat("sd_within", rl$sd_within, "\n")
d5 <- d2; names(d5)[1] <- ".y"
rm2 <- try1("factor named .y", m.bay(y ~ .y * B, data = d5, verbose = FALSE, sensitivity = FALSE)); show(rm2)
d6 <- d2; names(d6)[1:2] <- c("f2", "f1")
rn <- try1("factors named f2 f1", m.bay(y ~ f2 * f1, data = d6, verbose = FALSE, sensitivity = FALSE)); show(rn)
cat("compare with A*B:\n"); show(r)
