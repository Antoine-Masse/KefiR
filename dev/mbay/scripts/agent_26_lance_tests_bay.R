# agent_26 : lance les tests du socle bayesien et de la voie a un facteur.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet = TRUE))
for (f in commandArgs(TRUE)) {
  d <- as.data.frame(testthat::test_file(file.path("C:/Users/masse/Desktop/KefiR/KefiR-mbay/tests/testthat", f), reporter = "silent", stop_on_failure = FALSE))
  cat(f, ": tests", nrow(d), "expect", sum(d$nb), "failed", sum(d$failed), "errors", sum(d$error), "skipped", sum(d$skipped), "warnings", sum(d$warning), "\n")
}
