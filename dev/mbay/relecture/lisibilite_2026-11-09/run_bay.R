.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
Sys.setenv(LANGUAGE="fr")
suppressMessages(devtools::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE))
D <- "C:/Users/masse/AppData/Local/Temp/claude/dd1494da25364fcc/4e26779c-190e-41bb-89d1-f07a6f6766b4/scratchpad/audit_lisibilite"
source(file.path(D, "donnees.R"))
options(width = 120)
for (nm in names(cas)) {
  out <- file.path(D, paste0("m.bay", "_", nm, ".txt"))
  con <- file(out, open = "wt", encoding = "UTF-8")
  sink(con); sink(con, type = "message")
  t0 <- Sys.time()
  res <- tryCatch({ v <- withVisible(cas[[nm]](function(...) m.bay(..., verbose = TRUE, plot = FALSE))); if (v$visible) print(v$value) },
                  error = function(e) cat("ERREUR :", conditionMessage(e), "\n"))
  sink(type = "message"); sink(); close(con)
  cat(nm, "fait en", round(as.numeric(Sys.time()-t0, units="secs")), "s\n")
}
