# Sorties verbose de m.bay sur les 6 plans de l'audit de lisibilite et un plan discret,
# en francais et en anglais, largeur 80 ; mesures de longueur (tabulation comptee 8).
# Usage : Rscript run_mbay_7plans.R <sous-dossier> (avant, apres)
args <- commandArgs(TRUE)
sous <- if (length(args)) args[1L] else "apres"
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
D <- file.path(BASE, "dev/mbay/relecture/lisibilite_2026-11-09")
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
source(file.path(D, "donnees.R"))
set.seed(41)
dj <- data.frame(traitement = factor(rep(c("temoin", "froid", "chaud"), each = 30), levels = c("temoin", "froid", "chaud")))
dj$germe <- factor(ifelse(rbinom(90, 1, c(.5, .5, .8)[dj$traitement]) == 1, "oui", "non"), levels = c("non", "oui"))
cas[["7_discret"]] <- function(f) f(germe ~ traitement, data = dj)
OUT <- file.path(D, sous)
dir.create(OUT, showWarnings = FALSE)
options(width = 80)
mes <- character(0)
for (lang in c("fr", "en")) {
	if (lang == "en") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8")) else invisible(Sys.setlocale("LC_CTYPE", "French_France.utf8"))
	for (nm in names(cas)) {
		txt <- capture.output(cas[[nm]](function(...) m.bay(..., verbose = TRUE)))
		writeLines(txt, file.path(OUT, paste0("m.bay_", nm, "_", lang, ".txt")), useBytes = FALSE)
		w <- nchar(gsub("\t", "        ", txt), type = "width")
		mes <- c(mes, sprintf("%s %-18s lignes %3d | > 80 : %3d | max %3d | moyenne %5.1f | mots %4d", lang, nm, length(txt), sum(w > 80), max(w),
			mean(w[w > 0]), sum(lengths(strsplit(txt, "[[:space:]]+")))))
	}
}
writeLines(mes, file.path(OUT, "mesures.txt"))
cat(mes, sep = "\n")
