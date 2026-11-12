# agent_48 : KefiR contre la reference Monte-Carlo naive independante d'agent_45
# (t_het et normal_het, k = 3, jeux 4501 a 4503). Reference sur y brut ; KefiR standardise y
# par l'ecart-type intra-groupe poole s_p : log m(y_std) = log m(y) + (N - 1) log s_p.
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"; .libPaths(c(LIB, .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R")); .kefir_bay(BASE)
ref <- read.table(file.path(BASE, "dev/mbay/scripts/resultats/agent_45.txt"), header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
donnees <- function(j) {
  if (j == 1) { set.seed(4501); y <- c(rnorm(4,0,1), rnorm(5,0.3,1), rnorm(4,0,2)); g <- rep(1:3,c(4,5,4)) }
  if (j == 2) { set.seed(4502); y <- c(rnorm(6,0,1), rnorm(6,2,1), rnorm(6,0,1)); g <- rep(1:3,each=6) }
  if (j == 3) { set.seed(4503); y <- c(rnorm(5,0,0.5), rnorm(5,0,1), rnorm(5,1,3)); g <- rep(1:3,each=5) }
  list(y = y, g = g)
}
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
out <- NULL
for (j in 1:3) {
  d <- donnees(j); st <- .bay_laws_stats(d$y, d$g); N <- st$N; k <- st$k
  sp <- sqrt(sum(tapply(d$y, d$g, function(v) sum((v - mean(v))^2))) / (N - k))
  off <- (N - 1) * log(sp)
  R <- function(loi, f) ref$log_m[ref$jeu == j & ref$loi == loi & ref$forme == f] + off
  Rse <- function(loi, f) ref$se[ref$jeu == j & ref$loi == loi & ref$forme == f]
  mixref <- function(loi) lse(c(log(.5) + R(loi, "tous"), log(.5 / k) + sapply(paste0("s", 1:k), function(f) R(loi, f))))
  for (rep in 1:3) {
    set.seed(100 * j + rep)
    nm <- .bay_law_normal(st, 0.5, 1.5, 4000L, 100 * j + rep)
    set.seed(100 * j + rep)
    t0 <- .bay_t_ml(st, FALSE, TRUE, 1.5, 0.5, 4000L)
    th <- .bay_t_h1_mix(st, TRUE, 1.5, 0.5, 4000L)
    ts <- sapply(1:k, function(s) .bay_t_ml(st, TRUE, TRUE, 1.5, 0.5, 4000L, s = s)$log_ml)
    add <- function(loi, forme, kef, re, se) out <<- rbind(out, data.frame(jeu = j, rep = rep, loi = loi, forme = forme, kefir = kef, ref = re, ecart = kef - re, se_ref = se))
    add("normal_het", "H0", nm$het_h0$log_ml, R("normal_het", "H0"), Rse("normal_het", "H0"))
    add("normal_het", "tous", nm$het_h1$log_ml_zs, R("normal_het", "tous"), Rse("normal_het", "tous"))
    add("normal_het", "melange", nm$het_h1$log_ml, mixref("normal_het"), NA)
    add("t_het", "H0", t0$log_ml, R("t_het", "H0"), Rse("t_het", "H0"))
    add("t_het", "tous", th$log_ml_zs, R("t_het", "tous"), Rse("t_het", "tous"))
    for (s in 1:k) add("t_het", paste0("s", s), ts[s], R("t_het", paste0("s", s)), Rse("t_het", paste0("s", s)))
    add("t_het", "melange", th$log_ml, mixref("t_het"), NA)
  }
}
out$ecart <- round(out$ecart, 4)
agg <- aggregate(ecart ~ jeu + loi + forme, out, function(x) c(moy = mean(x), min = min(x), max = max(x)))
sink(file.path(BASE, "dev/mbay/scripts/resultats/agent_48.txt"))
cat(R.version.string, "\nKefiR (4000 tirages, 3 graines) moins reference agent_45, en log m ; y standardise.\n")
print(do.call(data.frame, agg), digits = 4)
lbf <- aggregate(ecart ~ jeu + loi + rep, subset(out, forme %in% c("H0", "melange")), function(x) x[2] - x[1])
cat("\nEcart sur le log BF du melange contre H0 (KefiR - reference) :\n"); print(lbf, digits = 3)
sink()
