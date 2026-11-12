# agent_63 : temps typique de m.bay sur un jeu de TP (N ~ 40), sensibilite et post-hoc par defaut.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
set.seed(1); g2 <- factor(rep(c("A","B"), each = 20)); g4 <- factor(rep(letters[1:4], each = 10))
id <- factor(rep(1:20, 2)); t2 <- factor(rep(c("t1","t2"), each = 20))
id3 <- factor(rep(1:14, 3)); t3 <- factor(rep(c("T1","T2","T3"), each = 14))
d22 <- expand.grid(rep = 1:10, A = factor(c("a1","a2")), B = factor(c("b1","b2")))
cas <- list(
	"2 groupes" = function() m.bay(rnorm(40), g2, verbose = FALSE),
	"4 groupes" = function() m.bay(rnorm(40), g4, verbose = FALSE),
	"2 x 2 facteurs" = function() m.bay(rnorm(40) ~ A * B, data = d22, verbose = FALSE),
	"apparie" = function() m.bay(rnorm(40), t2, id = id, verbose = FALSE),
	"repete 3 temps" = function() m.bay(rnorm(42), t3, id = id3, verbose = FALSE),
	"ANCOVA" = function() m.bay(y ~ g + x, data = data.frame(y = rnorm(40), g = g4, x = runif(40)), verbose = FALSE),
	"oui/non" = function() m.bay(factor(sample(c("o","n"), 40, TRUE)), g4, verbose = FALSE),
	"comptages" = function() m.bay(rpois(40, 3), g4, verbose = FALSE),
	"Likert 1-5" = function() m.bay(sample(1:5, 40, TRUE), g4, verbose = FALSE))
out <- "agent_63 : temps typique, N ~ 40, valeurs par defaut (secondes, une mesure)"
for (nm in names(cas)) out <- c(out, sprintf("%-16s %5.1f s", nm, system.time(suppressWarnings(cas[[nm]]()))[["elapsed"]]))
writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_63.txt")); cat(out, sep = "\n")
