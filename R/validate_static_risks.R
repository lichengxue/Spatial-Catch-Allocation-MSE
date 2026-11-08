# Independently validate risk probabilities from compact archived OM trajectories.
# Rscript --vanilla R/validate_static_risks.R
options(digits=15)
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
project<-dirname(dirname(script));folder<-file.path(project,'data/reference_points')
scenarios<-c('NB','NBH','NU','NUH','NBI','NBO','NBC','NUC','NBHS','MB','MBH','MU','MUH','MBI','MBO','MBHS')
results<-list()
for(s in scenarios){
 cache<-readRDS(file.path(folder,paste0('reference_cache_repaired_',s,'.RDS')))
 expected<-readRDS(file.path(project,'data/performance',paste0('performance2_',s,'.RDS')))
 expected<-expected[as.character(expected$Period)=='all',]
 stopifnot(length(cache)==500,nrow(expected)==1500)
 calculated<-do.call(rbind,lapply(cache,function(z){
  stopifnot(all(is.finite(c(z$F,z$SSB,z$static_F,z$static_SSB))),all(c(z$static_F,z$static_SSB)>0),abs(z$static_spr_ratio-.4)<1e-5)
  data.frame(Realization=z$Realization,Model=z$Model,Region=c('R1','R2','G'),pof=colMeans(sweep(z$F[36:85,,drop=FALSE],2,z$static_F,'/')>1),pofd=colMeans(sweep(z$SSB[36:85,,drop=FALSE],2,z$static_SSB,'/')<.5))
 }))
 key<-function(x)paste(x$Realization,x$Model,x$Region,sep='|')
 idx<-match(key(expected),key(calculated));stopifnot(!anyNA(idx))
 err<-max(abs(expected$p_overfishing-calculated$pof[idx]),abs(expected$p_overfished-calculated$pofd[idx]))
 stopifnot(err<1e-12)
 results[[s]]<-data.frame(Scenario=s,Realizations=100,Strategies=5,Spatial_scales=3,Max_probability_error=err)
}
dest<-file.path(project,'results');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,results),file.path(dest,'Static_Risk_Validation.csv'),row.names=FALSE)
cat('All 48,000 risk probabilities (24,000 strategy/scale cases, two risks each) match the static input data exactly.\n')
