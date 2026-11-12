# Correction de B2 (session m.bay) : ref.R mesurait les partitions non nulles contre l intercept seul,
# alors que la partition nulle (sujet seul) etait a 0. Ici, toutes sont lues contre lmBF(rt ~ sujet).
## agent_55 : calcul de reference independant (JZS, Rouder 2009 / 2012)
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressPackageStartupMessages(library(BayesFactor))

## ---------- B ----------
set.seed(22)
sujet <- factor(rep(1:12, 4)); session <- factor(rep(paste0("S", 1:4), each = 12))
rt <- 450 + rnorm(12, 0, 40)[sujet] + c(0, -5, -25, -30)[session] + rnorm(48, 0, 15)
dB <- data.frame(rt, sujet, session)
## ---------- B2 : partitions ----------
parts <- list(); lv <- paste0("S", 1:4)
rec <- function(i, lab) { if (i > 4) { parts[[length(parts) + 1]] <<- lab; return() }
  for (k in 1:(max(c(0, lab)) + 1)) rec(i + 1, c(lab, k)) }
rec(1, integer(0))
cat("=== B2 === nombre de partitions :", length(parts), "\n")
set.seed(2201)
lb0 <- extractBF(lmBF(rt ~ sujet, data = dB, whichRandom = "sujet", iterations = 50000, progress = FALSE), logbf = TRUE)$bf
res <- do.call(rbind, lapply(parts, function(p) {
  nom <- paste(sapply(split(lv, p), paste, collapse = "="), collapse = " | ")
  sizes <- table(p); w <- prod(factorial(sizes - 1))
  if (length(sizes) == 1) lbf <- 0 else {
    dB$bloc <- factor(p[as.integer(dB$session)])
    lbf <- extractBF(lmBF(rt ~ bloc + sujet, data = dB, whichRandom = "sujet", rscaleFixed = 0.5,
                          iterations = 50000, progress = FALSE), logbf = TRUE)$bf - lb0 }
  data.frame(partition = nom, k = length(sizes), w = w, logBF = lbf, stringsAsFactors = FALSE) }))
res$prior <- ifelse(res$k == 1, 0.5, 0.5 * res$w / sum(res$w[res$k > 1]))
lp <- log(res$prior) + res$logBF; res$post <- exp(lp - max(lp)); res$post <- res$post / sum(res$post)
res <- res[order(-res$post), ]
print(res, row.names = FALSE, digits = 4)
cat("\nTrois meilleures :\n"); print(head(res[, c("partition", "post")], 3), row.names = FALSE)
cat("\nP(les deux sessions different) :\n")
labs <- lapply(parts, identity); postv <- setNames(res$post, res$partition)
for (i in 1:3) for (j in (i + 1):4) {
  pr <- sum(sapply(seq_along(parts), function(q) {
    p <- parts[[q]]; nom <- paste(sapply(split(lv, p), paste, collapse = "="), collapse = " | ")
    if (p[i] != p[j]) postv[[nom]] else 0 }))
  cat(sprintf("  %s vs %s : %.4f\n", lv[i], lv[j], pr)) }

## ---------- C ----------
set.seed(23)
dC <- expand.grid(id = factor(1:20), seance = factor(c("s1", "s2", "s3")))
dC$groupe <- factor(ifelse(as.integer(dC$id) <= 10, "temoin", "entraine"))
dC$score <- 50 + rnorm(20, 0, 6)[dC$id] + 4 * (dC$seance != "s1") * (dC$groupe == "entraine") + rnorm(nrow(dC), 0, 3)
