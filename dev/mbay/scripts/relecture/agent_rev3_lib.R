.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet=TRUE, helpers=FALSE, attach_testthat=FALSE))
tab_of <- function(y,g,r=.5,u=1.5,draws=2000L,seed=1L,laws=.bay_law_names()) {
  st <- .bay_laws_stats(y,g); .bay_laws_table(st,r,u,draws,seed,laws=laws)
}
ans_of <- function(...) { a <- list(...); r <- if(is.null(a$r)) .5 else a$r; .bay_laws_average(tab_of(...), r) }
show <- function(lab,a) cat(sprintf("%-34s logBF %8.2f P(H1) %.4f err %6.2f%% | lois: %s | logBF/loi: %s\n",lab,a$log_bf,a$p_h1,a$error_pct,
  paste(sprintf("%s=%.3f",names(a$p_law),a$p_law),collapse=" "),paste(sprintf("%.1f",a$log_bf_law),collapse=" ")))
