# agent_load_bay : charge seulement les fonctions bayesiennes de KefiR, sans ses
# nombreuses dependances. Un processus pese ainsi environ 60 Mo au lieu de 600 Mo : a
# utiliser dans les simulations a plusieurs coeurs, sur une machine de 8 Go.
# Usage : source("dev/mbay/scripts/agent_load_bay.R"); .kefir_bay(src)
.kefir_bay <- function(src = Sys.getenv("KEFIR_SRC", "C:/Users/masse/Desktop/KefiR/KefiR-mbay"),
                       envir = globalenv()) {
	files <- c("sys_msg.R", "sys_exit.R", "sys_vbse.R", "discret.test.R", "sys_bay_core.R",
		"sys_bay_oneway.R", "sys_bay_laws.R", "sys_bay_input.R", "sys_bay_bilan.R", "m.bay.R",
		"sys_normalize_formula_dollar.R", list.files(file.path(src, "R"), "^sys_bay_.*[.]R$"))
	for (f in unique(files)) {
		p <- file.path(src, "R", f)
		if (file.exists(p)) sys.source(p, envir = envir, keep.source = FALSE)
	}
	invisible(TRUE)
}
