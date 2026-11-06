# Agent 28 : tables du rapport a partir des rds agent_24 a 27.
source("agent_24_commun.R")
md <- function(df) { d <- df; d[] <- lapply(df, function(v) if (is.numeric(v)) ifelse(is.na(v), "NA", format(v, trim = TRUE, drop0trailing = TRUE, digits = 7)) else as.character(v))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(l) paste0("| ", paste(l, collapse = " | "), " |"))) }
rd <- function(f) { p <- paste0(SORTIES, f); if (file.exists(p)) readRDS(p) else NULL }
V1 <- rd("agent_24_v1.rds"); B <- rd("agent_25_v3abc.rds"); D <- rd("agent_26_v3def.rds"); T27 <- rd("agent_27_troncature.rds"); R21 <- rd("agent_21_heterosc.rds")
L <- c("# Agent 24 : verification de .bay_bf_oneway_hetero et .bay_bf_oneway_zs", "",
  paste("Genere le", format(Sys.time()), "; scripts agent_24_v1v2.R, agent_25_v3abc.R, agent_26_v3def.R, agent_27_troncature.R, agent_28_rapport.R."), "")
add <- function(...) L <<- c(L, ...)
P <- plogis
if (!is.null(V1$tab)) { t <- V1$tab
  add("## V1 Calcul independant (grille deterministe sur log s_i, log g ; hetero n_draws=20000, seed=1)", "",
    md(data.frame(jeu=t$jeu, lBF_grille=round(t$lBF_grille,5), err_quad_pct=signif(t$err_quad_pct,2), lBF_fn=round(t$lBF_fn,5), error_pct=round(t$err_pct_fn,3), ecart_pct=round(t$ecart_pct,3), z=round(t$z,2), depasse=ifelse(abs(t$z)>1,"oui","non"), depasse_2=ifelse(abs(t$z)>2,"oui","non"))), "",
    "Marginales : ecart h1 et h0 (fonction - grille), (N-1) log(tau_hat) ; la route B (k=2) passe par la Cauchy et z, sans g.", "",
    md(data.frame(jeu=t$jeu, lh1_grille=round(t$lh1_grille,5), lh1_routeB=round(t$lh1_routeB,5), lh1_fn=round(t$lh1_fn,5), lh0_grille=round(t$lh0_grille,5), lh0_fn=round(t$lh0_fn,5), d_h1=round(t$d_h1_brut,5), d_h0=round(t$d_h0_brut,5), N1_logtau=round(t$N1logtau,5))), "")
  if (!is.null(V1$reduc)) add("Reduction analytique k=3 contre integrale numerique (mu, z1, z2) :", "", md(as.data.frame(round(V1$reduc,6))), "")
  if (!is.null(V1$couv)) add("Couverture 30 graines (biais contre grille, en %) :", "", md(data.frame(lapply(V1$couv, function(x) if (is.numeric(x)) signif(x,3) else x))), "") }
if (!is.null(V1$V2)) { v <- V1$V2
  add("## V2 zs a F grand (log BF)", "", md(data.frame(n=v$n, F_vise=signif(v$F_vise,2), F_obs=signif(v$F_obs,4), zs=round(v$zs,6), integrate=round(v$integrate,6), mpfr=round(v$mpfr,6), anovaBF=round(v$anovaBF,6), err_anovaBF=signif(v$err_anovaBF,2), zs_m_int=signif(v$ecart_zs_int,2), zs_m_mpfr=signif(v$ecart_zs_mpfr,2), zs_m_BF=signif(v$ecart_zs_BF,2))), "") }
if (!is.null(B)) {
  if (!is.null(B$ndraws)) { nd <- as.data.frame(B$ndraws)
    f <- function(col, s) { d <- nd[nd$sc==s,]; c(sd_lbf=sd(d[[col]]-d$nd20000), max_dP=max(abs(P(d[[col]])-P(d$nd20000)))) }
    x <- expand.grid(col=c("nd500","nd2000","nd4000"), sc=1:2, stringsAsFactors=FALSE)
    x$sd_ecart_logBF <- NA; x$max_abs_dP <- NA; for (i in seq_len(nrow(x))) { a <- f(x$col[i], x$sc[i]); x$sd_ecart_logBF[i] <- signif(a[1],2); x$max_abs_dP[i] <- signif(a[2],2) }
    x$sc <- c("A(4,4,1)","B(1,1,4)")[x$sc]
    add("## V3 a) H0, n=(5,5,40), appariement agent_21", "", "Effet de n_draws (25 jeux par scenario, ecart a n_draws=20000) :", "", md(x), "") }
  for (sc in c("A","B")) if (!is.null(B[[paste0("a_",sc)]])) { a <- as.data.frame(B[[paste0("a_",sc)]]); a <- a[!is.na(a$gj),]; r <- R21[[sc]]; r <- r[match(a$gj, r$graine_jeu),]; n <- nrow(a)
    Ph <- P(a$log_bf); Pg <- r$P; ligne <- function(nom, p) data.frame(methode=nom, N=n, `P>0.5`=fmt_ic(sum(p>.5),n), `P>0.75`=fmt_ic(sum(p>.75),n), `P>0.95`=fmt_ic(sum(p>.95),n), P_mediane=round(median(p),3), check.names=FALSE)
    fq <- function(nom,p) data.frame(test=nom, `p<0.05`=fmt_ic(sum(p<.05),n), check.names=FALSE)
    add(paste0("Scenario ",sc," sd=(",paste(R21$plan[[paste0("sd_",sc)]],collapse=","),"), hetero n_draws=2000 :"), "",
      md(rbind(ligne("hetero",Ph), ligne("anovaBF agent_21",Pg), ligne("zs",P(a$zs)))), "", md(rbind(fq("aov",a$p_aov), fq("Welch",a$p_welch), fq("Kruskal",a$p_kw))), "",
      sprintf("Appariement : hetero>0.5 et gauss<=0.5 : %d ; hetero<=0.5 et gauss>0.5 : %d ; hetero>0.95 et gauss<=0.95 : %d ; inverse : %d. max|zs-log_bf agent_21| = %.2g ; max|p_welch-agent_21| = %.2g ; erreur hetero mediane %.2f %%, max %.2f %%, ESS min %.0f ; temps median %.3f s ; low_ess %d.",
        sum(Ph>.5&Pg<=.5), sum(Ph<=.5&Pg>.5), sum(Ph>.95&Pg<=.95), sum(Ph<=.95&Pg>.95), max(abs(a$zs-r$log_bf)), max(abs(a$p_welch-r$p_welch)), median(a$err), max(a$err), min(a$ess), median(a$t), sum(a$low_ess)), "") }
  if (!is.null(B$b)) { rows <- list(); for (ci in seq_along(B$b)) if (!is.null(B$b[[ci]])) { m <- as.data.frame(B$b[[ci]]); n<-nrow(m); Ph<-P(m$het); Pz<-P(m$zs)
      rows[[ci]] <- data.frame(k=B$b_cells$k[ci], n=B$b_cells$n[ci], N=n, h_P05=fmt_ic(sum(Ph>.5),n), h_P95=fmt_ic(sum(Ph>.95),n), h_med=round(median(Ph),3), z_P05=fmt_ic(sum(Pz>.5),n), z_P95=fmt_ic(sum(Pz>.95),n), z_med=round(median(Pz),3), aov_05=round(mean(m$aov<.05),3), welch_05=round(mean(m$welch<.05),3), corr_lBF=round(cor(m$het,m$zs),3)) }
    add("## V3 b) Homoscedastique H0, equilibre, hetero n_draws=500", "", md(do.call(rbind,rows)), "") }
  if (!is.null(B$c)) { rows <- list(); for (ci in seq_along(B$c)) if (!is.null(B$c[[ci]])) { m <- as.data.frame(B$c[[ci]]); n<-nrow(m); Ph<-P(m$het); Pz<-P(m$zs); cc<-B$c_cells[ci,]
      rows[[ci]] <- data.frame(k=cc$k, d=cc$d, n=cc$n, cond=cc$cond, N=n, h_P95=round(mean(Ph>.95),3), h_P05=round(mean(Ph<.05),3), h_med=round(median(Ph),3), z_P95=round(mean(Pz>.95),3), z_P05=round(mean(Pz<.05),3), z_med=round(median(Pz),3), welch=round(mean(m$welch<.05),3), aov=round(mean(m$aov<.05),3)) }
    add("## V3 c) Puissance (delta = d x tau, tau = racine de la moyenne des sd^2 ; effet porte par le dernier groupe), hetero n_draws=500", "", "h_ = hetero, z_ = zs ; P95 : P>0.95 ; P05 : P<0.05 ; welch, aov : p<0.05. IC binomial exact : demi-largeur de l'ordre de 0.10 a 0.13 a N=60-100 au milieu de l'echelle.", "", md(do.call(rbind,rows)), "")
    tb <- do.call(rbind,rows); h <- tb[tb$cond=="egal",]; add(sprintf("Homoscedastique (%d cases, N total par methode %d) : moyenne h_P95 %.3f, z_P95 %.3f, welch %.3f, aov %.3f ; difference appariee h_P95-z_P95 moyenne %.3f, min %.3f, max %.3f.", nrow(h), sum(h$N), mean(h$h_P95), mean(h$z_P95), mean(h$welch), mean(h$aov), mean(h$h_P95-h$z_P95), min(h$h_P95-h$z_P95), max(h$h_P95-h$z_P95)), "") } }
cl <- function(p, h) { br <- c(-1e-9,.2,.4,.6,.8,1); ci <- cut(p, br); do.call(rbind, lapply(levels(ci), function(l) { i <- ci==l & !is.na(ci); if (!any(i)) return(data.frame(classe=l,m=0,P_moy=NA,freq_H1=NA,z=NA)); e<-sum(p[i]*(1-p[i])); data.frame(classe=l, m=sum(i), P_moy=round(mean(p[i]),3), freq_H1=round(mean(h[i]),3), z=round((sum(h[i])-sum(p[i]))/sqrt(e),2)) })) }
if (!is.null(D)) {
  if (!is.null(D$d)) { add("## V3 d) Calibration sous l'a priori (r=1/2, n_draws=500)", ""); for (ci in seq_along(D$d)) if (!is.null(D$d[[ci]])) { m<-as.data.frame(D$d[[ci]]); m<-m[!is.na(m$log_bf),]; p<-P(m$log_bf)
      add(sprintf("Plage sigma [1/%g;%g], n=(%s), N=%d, echecs %d, H1 %.3f, z global %.2f, err_pct mediane %.2f :", D$d_cells$R[ci], D$d_cells$R[ci], D$d_cells$des[ci], nrow(m), nrow(D$d[[ci]])-nrow(m), mean(m$H1), (sum(m$H1)-sum(p))/sqrt(sum(p*(1-p))), median(m$err)), "", md(cl(p,m$H1)), "") } }
  if (!is.null(T27)) { add("Ecart lie a la troncature : P de hetero (a priori impropre) contre P du meme modele avec a priori tronque a la plage (grille agent_24, trapeze avec demi-poids aux bornes, memes jeux) :", ""); rows<-list(); for (ci in seq_along(T27$res)) if (!is.null(T27$res[[ci]])) { m<-as.data.frame(T27$res[[ci]]); m<-m[is.finite(m$lbf_tr)&is.finite(m$lbf_fn),]; pt<-P(m$lbf_tr); pf<-P(m$lbf_fn)
      rows[[ci]]<-data.frame(R=T27$cells$R[ci],n=T27$cells$des[ci],N=nrow(m),moy_dP=round(mean(pf-pt),4),moy_abs_dP=round(mean(abs(pf-pt)),4),max_abs_dP=round(max(abs(pf-pt)),3),z_hetero=round((sum(m$H1)-sum(pf))/sqrt(sum(pf*(1-pf))),2),z_tronque=round((sum(m$H1)-sum(pt))/sqrt(sum(pt*(1-pt))),2),change_cote=sum((pf>.5)!=(pt>.5))) }
    add(md(do.call(rbind,rows)),""); for (ci in seq_along(T27$res)) if (!is.null(T27$res[[ci]])) { m<-as.data.frame(T27$res[[ci]]); m<-m[is.finite(m$lbf_tr),]; add(sprintf("Classes, P tronque, R=%g, n=(%s) :",T27$cells$R[ci],T27$cells$des[ci]),"",md(cl(P(m$lbf_tr),m$H1)),"") } }
  if (!is.null(D$e)) { e<-as.data.frame(D$e); colnames(e)<-c("sc","i","r025","r05","r1"); nm<-c("a) A sd=(4,4,1) H0","a) B sd=(1,1,4) H0","c) k=3 d=.5 n=(90,90,30) sd=(1,1,3)","c) k=2 d=1 n=30 egal")
    rows<-lapply(1:4,function(s){d<-e[e$sc==s,];a<-P(d$r025)>.5;b<-P(d$r05)>.5;c<-P(d$r1)>.5; data.frame(scenario=nm[s],N=nrow(d),P_r025_gt_05=round(mean(a),3),P_r05_gt_05=round(mean(b),3),P_r1_gt_05=round(mean(c),3),change_025_05=round(mean(a!=b),3),change_05_1=round(mean(c!=b),3),change_au_moins_un=round(mean(a!=b|c!=b),3),med_dlBF_05_1=round(median(d$r1-d$r05),3))})
    add("## V3 e) Sensibilite a r (n_draws=500)","",md(do.call(rbind,rows)),"") }
  if (!is.null(D$f)) { f<-D$f; if (!is.data.frame(f)) f<-do.call(rbind,f); num<-vapply(f,is.numeric,TRUE); f[num]<-lapply(f[num],signif,4); add("## V3 f) Degeneres (10 graines, n_draws=4000)","",md(f),"") } }
writeLines(L, paste0(SORTIES,"agent_24_rapport.md"), useBytes=TRUE); cat(length(L), "lignes\n")
