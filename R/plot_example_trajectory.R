# Uses project and output paths set by reproduce_figures.R.
trajectory<-read.csv(file.path(project,'data/trajectories/NBL_Example_Trajectory.csv'))%>%
 pivot_longer(c(Catch,SSB),names_to='Metric',values_to='Value')%>%
 mutate(Strategy=factor(Strategy,levels=MODEL_LEVELS),Spatial_scale=factor(Spatial_scale,levels=c('Region 1','Region 2','System')),Metric=factor(Metric,levels=c('Catch','SSB')))
historical<-filter(trajectory,Year_index<36,Strategy=='Equal')
feedback<-filter(trajectory,Year_index>=36)
p<-ggplot()+geom_line(data=historical,aes(x=Year,y=Value,group=interaction(Metric,Spatial_scale)),color='black',linewidth=.85,alpha=.95)+
 geom_vline(xintercept=2024,linetype='dashed',color='gray45',linewidth=.6)+
 geom_line(data=feedback,aes(x=Year,y=Value,color=Strategy),linewidth=.85,alpha=.72)+
 facet_grid(Metric~Spatial_scale,scales='free_y')+scale_color_manual(values=STRATEGY_COLORS,drop=FALSE)+labs(x='Year',y=NULL,color='Strategy')+
 theme_bw(base_size=14)+theme(strip.text=element_text(face='bold',size=13),axis.title.x=element_text(face='bold',size=15),axis.text=element_text(size=12),legend.position='right',legend.title=element_text(face='bold',size=13),legend.text=element_text(size=12),panel.grid.minor=element_blank(),panel.grid.major=element_line(color='gray88',linewidth=.35))
save_figure(p,file.path(output,'figures/supplementary/Figure_S1.png'),width=11,height=7,dpi=300)
