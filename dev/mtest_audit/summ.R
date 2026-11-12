args <- commandArgs(TRUE)
x <- readRDS(args[1])
x$rej <- x$p < 0.05
cat("Erreurs :\n"); print(table(x$s[x$method == "ERROR"])); print(unique(x$err[x$method == "ERROR"]))
agg <- aggregate(cbind(rej = rej, n = !is.na(p)) ~ s + term, data = transform(x, rej = as.numeric(rej)),
                 FUN = sum, na.action = na.pass)
agg$taux <- round(100 * agg$rej / agg$n, 1)
print(agg[order(agg$s, agg$term), c("s", "term", "n", "taux")], row.names = FALSE)
cat("\nMethodes par scenario :\n")
g <- x[x$term == "global", ]
print(lapply(split(g$method, g$s), function(m) sort(table(m), decreasing = TRUE)))
cat("\nRejet global par methode :\n")
r <- aggregate(rej ~ s + method, data = g, FUN = function(v) round(100 * mean(v), 1))
n <- aggregate(rej ~ s + method, data = g, FUN = length)
r$n <- n$rej
print(r[order(r$s), ], row.names = FALSE)
# post-hocs : au moins une paire significative alors que le test global ne l'est pas
ph <- merge(x[x$term == "posthoc_min", c("s", "i", "p")], g[, c("s", "i", "p")], by = c("s", "i"),
            suffixes = c(".ph", ".gl"))
cat("\nPost-hoc min p < 0.05 (toutes reps) et conditionne au global non significatif :\n")
print(do.call(rbind, lapply(split(ph, ph$s), function(d) data.frame(s = d$s[1],
  ph_any = round(100 * mean(d$p.ph < 0.05, na.rm = TRUE), 1),
  ph_and_gl_ns = round(100 * mean(d$p.ph < 0.05 & d$p.gl >= 0.05, na.rm = TRUE), 1)))), row.names = FALSE)
