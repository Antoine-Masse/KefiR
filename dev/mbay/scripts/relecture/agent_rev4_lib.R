.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
setwd("C:/Users/masse/Desktop/KefiR/KefiR-mbay")
source("dev/mbay/scripts/agent_load_bay.R"); .kefir_bay()
try_it <- function(lab, expr) {
  t0 <- Sys.time()
  out <- tryCatch({ v <- suppressWarnings(expr); list(ok=TRUE, v=v) }, error=function(e) list(ok=FALSE, m=conditionMessage(e), cl=class(e), call=deparse1(conditionCall(e))))
  el <- round(as.numeric(difftime(Sys.time(),t0,units="secs")),1)
  if (out$ok) { v <- out$v
    if (inherits(v,"mbay")) cat(sprintf("%-42s OK P=%.4f logBF=%.2f err=%.2f%% lois=%s notes=%d (%.1fs)\n",lab,v$prob_effect,v$effects$log_bf[1],v$effects$error_pct[1],paste(sprintf("%s:%.3f",v$laws$law,v$laws$posterior),collapse=" "),length(v$notes$fr),el))
    else cat(sprintf("%-42s OK (%.1fs)\n",lab,el))
  } else cat(sprintf("%-42s ERREUR [%s] %s (appel: %s)\n",lab,paste(out$cl,collapse="/"),out$m,substr(out$call,1,70)))
  invisible(out)
}
