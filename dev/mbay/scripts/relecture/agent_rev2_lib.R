# Commun a la relecture adverse du 2026-11-06. Aucune modification de R/.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
root <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
for (f in c("sys_msg.R", "sys_exit.R", "sys_bay_core.R", "sys_bay_oneway.R")) source(file.path(root, "R", f))

# Reference deterministe de log m(y|H1) pour k = 2 ou 3 : quadrature produit sur les log s2_i,
# independante de la proposition, avec la meme fonction interne de g (a pas fin).
ref_h1 <- function(st, r = 0.5, h = 1/16, span_hi = 60, step = 0.2) {
  axes <- lapply(seq_len(st$k), function(i) {
    nu <- st$n[i] - 1; ss <- st$ss[i]
    lo <- log(ss) - log(qchisq(1 - 1e-14, nu))
    hi <- if (nu <= 8) log(ss) + span_hi else log(ss) - log(qchisq(1e-14, nu))
    st_ <- min(step, 0.5 * sqrt(2 / nu))
    seq(lo, hi, by = st_)
  })
  gr <- as.matrix(expand.grid(axes))
  stp <- vapply(axes, function(a) a[2] - a[1], 0)
  lq <- 0
  for (i in seq_len(st$k)) {
    nu <- st$n[i] - 1; ss <- st$ss[i]
    # densite de log s2 sous q
    lq <- lq + (nu/2) * log(ss/2) - lgamma(nu/2) - (nu/2) * gr[, i] - ss / (2 * exp(gr[, i]))
  }
  inner <- numeric(nrow(gr))
  chunk <- split(seq_len(nrow(gr)), ceiling(seq_len(nrow(gr)) / 20000))
  for (ix in chunk) inner[ix] <- .bay_hetero_log_h1_inner(exp(gr[ix, , drop = FALSE]), st, r, h)
  .bay_hetero_log_const(st) + .bay_log_sum_exp(lq + inner) + sum(log(stp))
}
ref_bf <- function(st, r = 0.5, ...) ref_h1(st, r, ...) - .bay_hetero_h0_exact(st)

# st direct a partir de resumes (n, m, ss), standardise comme .bay_hetero_prep
st_of <- function(n, m, ss) {
  m <- m - sum(n * m) / sum(n)
  unit <- sqrt(mean(ss / (n - 1)))
  list(k = length(n), n = n, m = m / unit, ss = ss / unit^2, N = sum(n), levels = as.character(seq_along(n)))
}
