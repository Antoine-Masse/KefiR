# Taux de rejet de m.test() sous H0, par scenario.
# Usage : Rscript sim_h0.R <pkg> <N> <out.rds> <workers> <scenarios...>
args <- commandArgs(TRUE)
pkg <- args[1]; N <- as.integer(args[2]); out <- args[3]; workers <- as.integer(args[4])
scen_names <- args[-(1:4)]
lib <- "C:/Users/masse/Desktop/KefiR/mbay-lib"

gen <- list(
  # un facteur, groupes independants
  norm3     = function() { g <- rep(c("A","B","C"), each = 10); list(x = rnorm(30), g = g) },
  hetero_c  = function() { n <- c(5,5,40); g <- rep(c("A","B","C"), n)
                           list(x = unlist(mapply(function(n,s) rnorm(n,0,s), n, c(4,4,1))), g = g) },
  hetero_c2 = function() { n <- c(5,5,40); g <- rep(c("A","B","C"), n)
                           list(x = unlist(mapply(function(n,s) s*(rexp(n)-1), n, c(4,4,1))), g = g) },
  skew3     = function() { g <- rep(c("A","B","C"), each = 15); list(x = rexp(45), g = g) },
  lnorm4    = function() { g <- rep(c("A","B","C","D"), each = 12); list(x = rlnorm(48), g = g) },
  norm2     = function() { g <- rep(c("A","B"), each = 12); list(x = rnorm(24), g = g) },
  hetero2   = function() { n <- c(6,30); g <- rep(c("A","B"), n)
                           list(x = c(rnorm(6,0,3), rnorm(30,0,1)), g = g) },
  skew2     = function() { g <- rep(c("A","B"), each = 15); list(x = rexp(30), g = g) },
  # appariement
  paired2   = function() { n <- 15; u <- rnorm(n, 0, 2); d <- data.frame(id = factor(rep(1:n, 2)),
                           T = rep(c("t1","t2"), each = n)); d$y <- u[as.integer(d$id)] + rnorm(2*n); d },
  rm3       = function() { n <- 15; u <- rnorm(n, 0, 2); d <- expand.grid(id = factor(1:n), T = c("t1","t2","t3"))
                           d$y <- u[as.integer(d$id)] + rnorm(3*n); d },
  rm3_lmer  = function() { n <- 15; u <- rnorm(n, 0, 2); d <- expand.grid(id = factor(1:n), T = c("t1","t2","t3"))
                           d$y <- u[as.integer(d$id)] + rnorm(3*n); d },
  split     = function() { n <- 12; u <- rnorm(n, 0, 1); d <- expand.grid(id = factor(1:n), T = c("t1","t2","t3"))
                           d$G <- rep(c("A","B"), each = n/2)[as.integer(d$id)]
                           d$y <- u[as.integer(d$id)] + rnorm(3*n, 0, 0.8); d },
  # deux facteurs inter
  twoway    = function() { d <- expand.grid(r = 1:8, A = c("a1","a2"), B = c("b1","b2","b3")); d$y <- rnorm(nrow(d)); d },
  twoway_sk = function() { d <- expand.grid(r = 1:8, A = c("a1","a2"), B = c("b1","b2","b3")); d$y <- rexp(nrow(d)); d },
  ancova    = function() { d <- data.frame(F = rep(c("A","B","C"), each = 15), x = runif(45, 0, 10))
                           d$y <- 2 + 0.5 * d$x + rnorm(45); d }
)

call_mtest <- function(s, dat) {
  switch(s,
    paired2  = m.test(y ~ T + Error(id), data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    rm3      = m.test(y ~ T + Error(id/T), data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    rm3_lmer = m.test(y ~ T + (1|id), data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    split    = m.test(y ~ G * T + Error(id/T), data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    twoway   = , twoway_sk = m.test(y ~ A * B, data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    ancova   = m.test(y ~ F + x, data = dat, plot = FALSE, verbose = FALSE, boot = FALSE),
    m.test(dat$x, dat$g, plot = FALSE, verbose = FALSE, boot = FALSE))
}

one <- function(i, s) {
  set.seed(1000 * match(s, names(gen)) + i)
  dat <- gen[[s]]()
  r <- tryCatch(suppressWarnings(suppressMessages({
    utils::capture.output(v <- call_mtest(s, dat)); v })), error = function(e) conditionMessage(e))
  if (is.character(r)) return(data.frame(s = s, i = i, method = "ERROR", term = "global", p = NA, err = r))
  rows <- data.frame(s = s, i = i, method = paste(r$method, collapse = "|"), term = "global",
                     p = if (is.null(r$global_pvalue)) NA else as.numeric(r$global_pvalue[1]), err = "")
  if (!is.null(r$effects) && is.data.frame(r$effects))
    rows <- rbind(rows, data.frame(s = s, i = i, method = rows$method, term = r$effects$term,
                                   p = r$effects$p.value, err = ""))
  # plus petite p-value des post-hocs (erreur de famille)
  pv <- r$p.value
  if (is.matrix(pv) && any(!is.na(pv)))
    rows <- rbind(rows, data.frame(s = s, i = i, method = rows$method[1], term = "posthoc_min",
                                   p = min(pv, na.rm = TRUE), err = ""))
  rows
}

cl <- parallel::makeCluster(workers)
parallel::clusterExport(cl, c("lib", "pkg", "gen", "call_mtest", "one"))
invisible(parallel::clusterEvalQ(cl, { .libPaths(c(lib, .libPaths()))
  suppressMessages(devtools::load_all(pkg, quiet = TRUE)) }))
res <- list()
for (s in scen_names) {
  t0 <- Sys.time()
  res[[s]] <- do.call(rbind, parallel::parLapply(cl, seq_len(N), one, s = s))
  cat(s, "done in", format(Sys.time() - t0), "\n")
  saveRDS(do.call(rbind, res), out)
}
parallel::stopCluster(cl)
