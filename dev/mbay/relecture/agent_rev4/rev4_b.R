.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_load_bay.R"); .kefir_bay()
suppressMessages(library(BayesFactor))
try1 <- function(lab, expr) {
	cat("\n=====", lab, "=====\n")
	r <- tryCatch(expr, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
	invisible(r)
}
show <- function(r) if (!is.null(r)) { print(r$effects[, c("term", "p_h1", "prior_h1", "error_pct")]); cat("p_null", r$hypotheses$posterior[r$hypotheses$null], " n_models", r$n_models, " notes:", r$notes$en, "\n") }
q <- function(...) m.bay(..., verbose = FALSE, sensitivity = FALSE)

set.seed(8)
d2 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 1:6)
d2$y <- rnorm(nrow(d2), 1.5 * (d2$A == "a2"))

## Aliased / nested factors
da <- d2; da$C <- factor(ifelse(da$A == "a1", "c1", "c2"))   # C == A relabelled
show(try1("aliased A + C", q(y ~ A + C, data = da)))
dn <- d2; dn$N <- factor(paste(dn$A, dn$B))                    # N nested-ish: N determines A
show(try1("A + N (N determines A)", q(y ~ A + N, data = dn)))
## Empty cell, additive model
e <- d2[!(d2$A == "a1" & d2$B == "b3"), ]
show(try1("empty cell additive", q(y ~ A + B, data = e)))
## one per cell
d1 <- d2[d2$r == 1, ]
show(try1("one per cell A*B", q(y ~ A * B, data = d1)))
show(try1("one per cell A+B", q(y ~ A + B, data = d1)))
d1b <- rbind(d1, d1[1, ])
show(try1("one per cell + 1 A*B", q(y ~ A * B, data = d1b)))
## unbalanced
set.seed(2); du <- d2[sample(nrow(d2), 25), ]
show(try1("unbalanced", q(y ~ A * B, data = du)))
## numeric factor
dnum <- d2; dnum$B <- as.integer(dnum$B)
show(try1("numeric B 3 values", q(y ~ A * B, data = dnum)))
## NA
dna <- d2; dna$A[3] <- NA; dna$y[5] <- NA
show(try1("NA", q(y ~ A * B, data = dna)))
## 3 factors
set.seed(4)
d3 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), C = c("c1", "c2", "c3"), r = 1:4)
d3$y <- rnorm(nrow(d3), (d3$A == "a2") * (d3$B == "b2") * 1.5)
t3 <- system.time(r3 <- try1("A*B*C", q(y ~ A * B * C, data = d3))); show(r3); print(t3)
show(try1("A*B + C", q(y ~ A * B + C, data = d3)))
show(try1("A*B*C - A:B:C", q(y ~ A * B * C - A:B:C, data = d3)))
## interaction only effect: crossing, no main effects
set.seed(5)
dx <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), r = 1:10)
dx$y <- rnorm(nrow(dx), ifelse(dx$A == dx$A[1] & dx$B == "b1" | dx$A != dx$A[1] & dx$B == "b2", 1, -1))
rx <- try1("pure crossover", m.bay(y ~ A * B, data = dx, verbose = FALSE)); show(rx)
cat("sd_within reported:", rx$sd_within, " cell-means sd:", sigma(lm(y ~ A * B, dx)), " additive sd:", sigma(lm(y ~ A + B, dx)), "\n")
print(rx)
## big effect: printing
set.seed(6)
db <- expand.grid(A = c("a1", "a2", "a3"), B = c("b1", "b2"), r = 1:8)
db$y <- rnorm(nrow(db), 3 * as.integer(db$A))
rbig <- try1("big effect", m.bay(y ~ A * B, data = db, verbose = FALSE)); show(rbig)
print(rbig, detail = TRUE)
## null data
set.seed(7)
d0 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), r = 1:15); d0$y <- rnorm(nrow(d0))
r0 <- m.bay(y ~ A + B, data = d0, verbose = FALSE); show(r0); print(r0)
## return = FALSE
cat("return FALSE:", m.bay(y ~ A * B, data = d2, return = FALSE, sensitivity = FALSE), "\n")
