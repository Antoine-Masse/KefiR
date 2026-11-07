.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet=TRUE, helpers=FALSE, attach_testthat=TRUE))
for (f in c("test-m.bay.R","test-bay-laws.R")) { r <- as.data.frame(testthat::test_file(file.path("C:/Users/masse/Desktop/KefiR/KefiR-mbay/tests/testthat",f),reporter="silent",stop_on_failure=FALSE,env=asNamespace("KefiR")))
 cat(f,": tests",nrow(r),"echecs",sum(r$failed),"erreurs",sum(r$error),"temps",round(sum(r$real)),"s\n"); print(r[r$failed>0|r$error,c("test","failed")]) }
