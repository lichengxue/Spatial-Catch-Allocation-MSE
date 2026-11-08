# Rebuild metrics from saved simulations and the matched repaired static cache.
# No MSE is rerun. Archives are large and are not part of the compact release.
# Rscript --vanilla R/extract_archived_performance.R <archive_folder> <destination> <raw_scenario>
extract_archived_performance<-function(mods,cache,model_names=c('Equal','Gear-Equal','Fleet-Catch','Index-Equal','Index-Gear')){
 regions<-c('R1','R2','G');rows<-list()
 for(i in seq_along(mods)){
  temporal<-vector('list',length(model_names))
  for(m in seq_along(model_names)){
   om<-mods[[i]][[m]]$om;r<-om$rep;d<-om$input$data
   stopifnot(nrow(r$SSB)==85,d$n_fleets==4,d$n_regions==2)
   C<-sapply(1:2,function(k)rowSums(r$pred_catch[,d$fleet_regions==k,drop=FALSE]));C<-cbind(C,rowSums(C))
   S<-cbind(r$SSB,rowSums(r$SSB));F<-r$Fbar[,5:7,drop=FALSE]
   hit<-which(vapply(cache,function(z)z$Realization==i&&z$Model==model_names[m],logical(1)));stopifnot(length(hit)==1)
   z<-cache[[hit]]
   stopifnot(max(abs(F-z$F))<1e-12,max(abs(S-z$SSB))<1e-9,all(is.finite(c(z$static_F,z$static_SSB))),all(c(z$static_F,z$static_SSB)>0),abs(z$static_spr_ratio-.4)<1e-5)
   aav<-function(a)colMeans(abs(a[37:85,,drop=FALSE]-a[36:84,,drop=FALSE]))/colMeans(a[36:85,,drop=FALSE])
   temporal[[m]]<-list(short_C=colMeans(C[36:45,,drop=FALSE]),short_S=colMeans(S[36:45,,drop=FALSE]),long_C=colMeans(C[76:85,,drop=FALSE]),long_S=colMeans(S[76:85,,drop=FALSE]),AAV_C=aav(C),AAV_S=aav(S),POF=colMeans(sweep(F[36:85,,drop=FALSE],2,z$static_F,'/')>1),POFD=colMeans(sweep(S[36:85,,drop=FALSE],2,z$static_SSB,'/')<.5))
  }
  for(m in seq_along(model_names))for(period in c('short','long','all')){
   z<-temporal[[m]];b<-temporal[[1]];raw<-period!='all'
   C<-if(raw)z[[paste0(period,'_C')]] else rep(NA_real_,3);S<-if(raw)z[[paste0(period,'_S')]] else rep(NA_real_,3)
   rows[[length(rows)+1]]<-data.frame(Realization=i,Model=model_names[m],Period=period,Scale=c('regional','regional','system'),Region=regions,SSB=S,Catch=C,
    rel_SSB=if(raw)(S-b[[paste0(period,'_S')]])/b[[paste0(period,'_S')]] else NA_real_,rel_Catch=if(raw)(C-b[[paste0(period,'_C')]])/b[[paste0(period,'_C')]] else NA_real_,
    AAV_SSB=if(raw)NA_real_ else z$AAV_S,AAV_Catch=if(raw)NA_real_ else z$AAV_C,p_overfishing=if(raw)NA_real_ else z$POF,p_overfished=if(raw)NA_real_ else z$POFD)
  }
 }
 do.call(rbind,rows)
}
if(sys.nframe()==0){
 args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop('Provide archive_folder, destination, raw_scenario.')
 script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
 project<-dirname(dirname(script));s<-args[3]
 cache<-readRDS(file.path(project,'data/reference_points',paste0('reference_cache_repaired_',s,'.RDS')))
 mods<-readRDS(file.path(args[1],paste0('Mods_',s,'.RDS')))
 result<-extract_archived_performance(mods,cache);result$Scenario<-s
 dir.create(args[2],recursive=TRUE,showWarnings=FALSE);saveRDS(result,file.path(args[2],paste0('performance2_',s,'.RDS')))
 cat('Extracted',nrow(result),'rows for',s,'without rerunning the MSE.\n')
}
