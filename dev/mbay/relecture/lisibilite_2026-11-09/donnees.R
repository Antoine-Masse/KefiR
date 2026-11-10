set.seed(42)
d2 <- data.frame(y = c(rnorm(20, 10, 2), rnorm(20, 11.5, 2)), g = factor(rep(c("A","B"), each = 20)))
d4 <- data.frame(y = c(rnorm(15, 10, 2), rnorm(15, 12, 2), rnorm(15, 10.5, 2), rnorm(15, 13, 2)),
                 g = factor(rep(c("a","b","c","d"), each = 15)))
d22 <- expand.grid(rep = 1:10, A = c("A1","A2"), B = c("B1","B2","B3"))
d22$y <- 10 + 1.5*(d22$A=="A2") + 1*(d22$B=="B3") + rnorm(nrow(d22), 0, 2)
d22$A <- factor(d22$A); d22$B <- factor(d22$B)
ns <- 15; ind <- rnorm(ns, 0, 2)
dp <- data.frame(id = factor(rep(1:ns, 2)), t = factor(rep(c("avant","apres"), each = ns), levels=c("avant","apres")))
dp$y <- 20 + ind[as.integer(dp$id)] + 1.2*(dp$t=="apres") + rnorm(2*ns, 0, 1.2)
dr <- data.frame(id = factor(rep(1:ns, 3)), t = factor(rep(c("T1","T2","T3"), each = ns)))
dr$y <- 20 + ind[as.integer(dr$id)] + c(0, 1, 1.8)[as.integer(dr$t)] + rnorm(3*ns, 0, 1.2)
da <- data.frame(g = factor(rep(c("t1","t2","t3"), each = 15)), x = runif(45, 0, 10))
da$y <- 5 + 0.8*da$x + c(0, 1.5, 3)[as.integer(da$g)] + rnorm(45, 0, 1.5)
cas <- list(
  "1_deux_groupes" = function(f) f(y ~ g, data = d2),
  "2_quatre_groupes" = function(f) f(y ~ g, data = d4),
  "3_deux_facteurs" = function(f) f(y ~ A * B, data = d22),
  "4_apparie" = function(f) f(y ~ t, data = dp, id = "id", paired = TRUE),
  "5_repete3" = function(f) f(y ~ t | id, data = dr),
  "6_ancova" = function(f) f(y ~ g + x, data = da))
