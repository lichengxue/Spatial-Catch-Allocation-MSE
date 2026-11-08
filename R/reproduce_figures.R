# Reproduce the final static-reference figure suite from compact performance data.
# Rscript --vanilla R/reproduce_figures.R [output_directory]
options(width=160,digits=15)
script_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
stopifnot(length(script_arg)==1)
code_dir<-dirname(normalizePath(sub('^--file=','',script_arg),winslash='/',mustWork=TRUE))
project<-dirname(code_dir)
args<-commandArgs(trailingOnly=TRUE)
output<-if(length(args))args[1] else project
for(p in c('figures/main','figures/supplementary','results'))dir.create(file.path(output,p),recursive=TRUE,showWarnings=FALSE)
# Windows user libraries are optional; no machine-specific path is required.
user_lib<-file.path(Sys.getenv('LOCALAPPDATA'),'R','win-library',paste(R.version$major,strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1],sep='.'))
if(dir.exists(user_lib)).libPaths(c(user_lib,.libPaths()))
required<-c('dplyr','tidyr','ggplot2','scales','patchwork')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Install required packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(grid);library(scales);library(patchwork)})
source(file.path(code_dir,'plot_helpers.R'))
base_theme<-theme_manuscript
SAVE_PDF<-TRUE
perf<-bind_rows(lapply(SCENARIOS_ALL_RAW,function(s){
 x<-readRDS(file.path(project,'data/performance',paste0('performance2_',s,'.RDS')))
 x$Scenario_raw<-s;x$Scenario<-recode_scenario(s);x
}))%>%mutate(Scenario=factor(Scenario,levels=SCENARIOS_ALL),Region=factor(as.character(Region),levels=REGION_LEVELS),Model=factor(as.character(Model),levels=MODEL_LEVELS),Period=factor(as.character(Period),levels=c('short','long','all')))
stopifnot(nrow(perf)==72000,!anyDuplicated(perf[c('Scenario','Realization','Region','Model','Period')]))
main<-filter(perf,Scenario%in%SCENARIOS_MAIN)
metric_real<-perf%>%transmute(Scenario,Realization,Region,Model,
 Catch_ST=if_else(Period=='short',Catch,NA_real_),SSB_ST=if_else(Period=='short',SSB,NA_real_),
 Catch_LT=if_else(Period=='long',Catch,NA_real_),SSB_LT=if_else(Period=='long',SSB,NA_real_),
 AAV_Catch_All=if_else(Period=='all',AAV_Catch,NA_real_),AAV_SSB_All=if_else(Period=='all',AAV_SSB,NA_real_),
 p_overfishing_All=if_else(Period=='all',p_overfishing,NA_real_),p_overfished_All=if_else(Period=='all',p_overfished,NA_real_))%>%
 pivot_longer(all_of(METRIC_LEVELS),names_to='Metric',values_to='Value_real')%>%filter(!is.na(Value_real))%>%mutate(Metric=factor(Metric,levels=METRIC_LEVELS))
stopifnot(nrow(metric_real)==192000,all(is.finite(metric_real$Value_real)),!anyDuplicated(metric_real[c('Scenario','Realization','Region','Model','Metric')]))
# No rel_Catch or rel_SSB enters the scores. Mean raw metrics first, then scale.
scores<-metric_real%>%group_by(Scenario,Region,Model,Metric)%>%summarise(Value=mean(Value_real),N_realizations=n(),.groups='drop')%>%group_by(Scenario,Region,Metric)%>%mutate(Score=score_metric_0_1(Value,Metric))%>%ungroup()
totals<-scores%>%group_by(Scenario,Region,Model)%>%summarise(Summed_score=sum(Score),Average_score=mean(Score),N_metrics=n(),.groups='drop')
stopifnot(nrow(scores)==1920,all(scores$N_realizations==100),all(is.finite(scores$Score)),all(scores$Score>=0&scores$Score<=1),nrow(totals)==240,all(totals$N_metrics==8))
write.csv(scores,file.path(output,'results/Metric_Scores.csv'),row.names=FALSE)
write.csv(totals,file.path(output,'results/Composite_Scores.csv'),row.names=FALSE)
write.csv(metric_real%>%group_by(Scenario,Region,Model,Metric)%>%summarise(Mean=mean(Value_real),Median=median(Value_real),Q1=quantile(Value_real,.25),Q3=quantile(Value_real,.75),N=n(),.groups='drop'),file.path(output,'results/Raw_Metric_Summaries.csv'),row.names=FALSE)
write.csv(scores%>%group_by(Scenario,Region,Metric)%>%summarise(Range=max(Value)-min(Value),Common_value=first(Value),Assigned_score=first(Score),.groups='drop')%>%filter(Range==0),file.path(output,'results/Zero_Range_Groups.csv'),row.names=FALSE)
winners<-totals%>%group_by(Scenario,Region)%>%summarise(Winner=paste(as.character(Model)[Average_score>=max(Average_score)-1e-12],collapse=' / '),Best_score=max(Average_score),.groups='drop')
write.csv(winners,file.path(output,'results/Highest_Scoring_Strategies.csv'),row.names=FALSE)

# Main and supplementary raw metrics: medians; IQR for the all-scenario figures.
for(all_om in c(FALSE,TRUE)){
 dat<-if(all_om)perf else main
 scen<-if(all_om)SCENARIOS_ALL_Y else SCENARIOS_MAIN_Y
 specifications<-if(all_om)list(c('6','Catch','Catch','Catch'),c('7','SSB','SSB','SSB'),c('8','AAV_Catch','AAV_Catch','Average annual variation in catch'),c('9','AAV_SSB','AAV_SSB','Average annual variation in SSB'),c('10','p_overfishing','POF','Probability of overfishing'),c('11','p_overfished','POFD','Probability of being overfished')) else list(c('3','Catch','Catch','Catch'),c('4','SSB','SSB','SSB'),c('S2','AAV_Catch','AAV_Catch','Average annual variation in catch'),c('S3','AAV_SSB','AAV_SSB','Average annual variation in SSB'),c('S4','p_overfishing','POF','Probability of overfishing'),c('S5','p_overfished','POFD','Probability of being overfished'))
 for(spec in specifications){
  number<-if(all_om)paste0('S',spec[1]) else spec[1]
  raw<-spec[2]%in%c('Catch','SSB');risk<-startsWith(spec[2],'p_')
  sumdat<-make_metric_summary(dat,spec[2],spec[3],scenario_levels_y=scen,period_keep=if(raw)c('short','long') else 'all',drop_equal=FALSE)
  labels<-if(raw)setNames(paste(spec[4],c('ST','LT')),paste0(spec[3],c('_ST','_LT'))) else setNames(if(risk)if(spec[2]=='p_overfishing')'P(overfishing)' else 'P(overfished)' else if(spec[2]=='AAV_Catch')'AAV catch' else 'AAV SSB',paste0(spec[3],'_All'))
  limits<-if(risk)c(0,1) else if(!raw)c(0,max(if(all_om)sumdat$Q3 else sumdat$Median)*1.08) else NULL
  filename<-file.path(output,if(grepl('^S',number))'figures/supplementary' else 'figures/main',paste0('Figure_',number,'.png'))
  theme_manuscript<-if(risk)function()base_theme()+theme(panel.spacing.x=unit(.22,'in')) else base_theme
  plot_metric_iqr(sumdat,filename,spec[4],labels,xlim=limits,auto_xlim=FALSE,facet_scales=if(raw)'free_x' else 'fixed',percent_axis=FALSE,add_zero_line=FALSE,show_iqr=all_om,median_segment_from=if(all_om)NULL else 0,width=if(all_om)15 else 14,height=if(all_om)24 else 18)
  write.csv(sumdat,file.path(output,'results',paste0('Figure_',number,'_Summary.csv')),row.names=FALSE)
  cat('Saved Figure',number,'\n');flush.console()
 }
}
for(region in c('G','R1','R2'))for(all_om in c(FALSE,TRUE)){
 number<-paste0(if(all_om)'S' else '',match(region,c('G','R1','R2'))+if(all_om)11 else 4)
 filename<-file.path(output,if(all_om)'figures/supplementary' else 'figures/main',paste0('Figure_',number,'.png'))
 plot_score_heatmap(region,scores,totals,if(all_om)SCENARIOS_ALL else SCENARIOS_MAIN,filename,paste0('Average scaled performance, ',switch(region,G='system level',R1='region 1',R2='region 2'),if(all_om)', all OM scenarios' else ''),ncol=if(all_om)4 else 3,width=if(all_om)15 else 12.5,height=if(all_om)17 else 10.5,score_col='Average_score',score_ylim=c(0,1),score_breaks=seq(0,1,by=.25))
 cat('Saved Figure',number,'\n');flush.console()
}
# Figures 1-2 are archived illustrations; numerical figures are regenerated.
for (asset in c('Figure_1.png', 'Figure_2.png')) {
  original <- file.path(project, 'figures/main', asset)
  destination <- file.path(output, 'figures/main', asset)
  if (normalizePath(original, winslash='/', mustWork=TRUE) !=
      normalizePath(destination, winslash='/', mustWork=FALSE)) {
    stopifnot(file.copy(original, destination, overwrite=TRUE))
  }
}
source(file.path(code_dir,'plot_example_trajectory.R'))
capture.output(sessionInfo(),file=file.path(output,'results/R_session_info.txt'))
cat('Completed 19 regenerated figures and 2 archived illustrations.\n')
