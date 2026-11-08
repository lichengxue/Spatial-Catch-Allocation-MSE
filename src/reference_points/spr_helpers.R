make_spr_data<-function(state,years){
 r<-state$rep;d<-state$data
 annual<-length(years)==1
 list(SPR_weights=as.numeric(d$SPR_weights),log_M=r$log_M,FAA=r$FAA,
  spawn_seasons=d$spawn_seasons,spawn_regions=d$spawn_regions,fleet_regions=d$fleet_regions,
  fleet_seasons=d$fleet_seasons,fracyr_seasons=d$fracyr_seasons,
  can_move=d$can_move,must_move=d$must_move,mig_type=d$mig_type,
  trans_mu_base=r$trans_mu_base,L=r$L,
  which_F_age=as.integer(if(annual)d$which_F_age[years+1] else d$which_F_age_static),
  waa_ssb=r$waa_ssb,waa_catch=r$waa_catch,mature=r$mature_all,
  percentSPR=d$percentSPR,NAA=r$NAA,fracyr_SSB=r$fracyr_SSB_all,
  years=as.integer(years),Fbar_ages=d$Fbar_ages,
  R_XSPR=as.numeric(r$R_XSPR[if(annual)years+1 else 85,]),
  small_dim=as.integer(d$n_regions_is_small),SPR_weight_type=as.integer(d$SPR_weight_type),
  bias_correct=as.integer(d$bias_correct_brps),marg_NAA_sigma=r$marg_NAA_sigma)
}
spr_evaluator<-function(state,years){
 obj<-TMB::MakeADFun(make_spr_data(state,years),list(logF=0),DLL='spr_evaluate',
    random=NULL,map=NULL,silent=TRUE,type='Fun')
 # Explicit parameter avoids TMB's AD-derived default parameter vector when
 # the DLL is compiled for numerical reports only.
 function(logF)obj$report(as.numeric(logF))
}
solve_spr<-function(state,years){
 evaluate<-spr_evaluator(state,years)
 objective<-function(logF)as.numeric(evaluate(logF)$ratio)-.4
 low<--20;high<-20
 stopifnot(objective(low)>0,objective(high)<0)
 root<-uniroot(objective,c(low,high),tol=1e-10)$root
 result<-evaluate(root)
 stopifnot(is.finite(result$ratio),abs(result$ratio-.4)<1e-8,
     all(is.finite(c(result$Fbar,result$SSBref))),all(c(result$Fbar,result$SSBref)>0))
 list(Fbar=result$Fbar,SSBref=result$SSBref,ratio=as.numeric(result$ratio),FXSPR=exp(root))
}
