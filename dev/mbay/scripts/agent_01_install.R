# Installation chronometree (binaire, dans mbay-lib uniquement). Ordre : leger -> lourd.
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=900, install.packages.compile.from.source="never")
dsz <- function() { f <- list.files(LIB, recursive=TRUE, full.names=TRUE, all.files=TRUE); sum(file.size(f),na.rm=TRUE)/1024^2 }
for (p in c("BayesFactor","bayestestR","BAS","bridgesampling","rstanarm","brms")) {
  before <- rownames(installed.packages(lib.loc=LIB)); s0 <- dsz(); t0 <- Sys.time()
  err <- tryCatch({install.packages(p, lib=LIB, repos="https://cloud.r-project.org", type="binary", quiet=TRUE); "ok"}, error=function(e) conditionMessage(e), warning=function(w) paste("WARN:",conditionMessage(w)))
  el <- as.numeric(difftime(Sys.time(),t0,units="secs"))
  new <- setdiff(rownames(installed.packages(lib.loc=LIB)),before)
  cat(sprintf("\n### %s : %.0f s ; statut=%s ; %d nouveaux paquets ; +%.1f Mo ; total lib %.1f Mo\n  nouveaux: %s\n", p, el, err, length(new), dsz()-s0, dsz(), paste(new,collapse=" ")))
  flush.console()
}
ip <- installed.packages(lib.loc=LIB)
print(ip[,c("Package","Version")], quote=FALSE)
# taille par paquet
sz <- sapply(rownames(ip), function(p) sum(file.size(list.files(file.path(LIB,p),recursive=TRUE,full.names=TRUE)),na.rm=TRUE)/1024^2)
print(round(sort(sz,decreasing=TRUE),1))
cat("FIN INSTALL\n")
