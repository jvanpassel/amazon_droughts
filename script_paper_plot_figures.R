library(ggplot2)
library(scales)
library(maptools)
library(scico)
library(khroma)
library(ggeffects)
library(sp)
library(precrec)
library(ggpattern)

#Load input data
amazon<-shapefile('amapoly_ivb.shp')
amazon_df<-fortify(amazon)
data("wrld_simpl")
worldmap<-fortify(wrld_simpl)

#Figure 1: Plot TAC linear trend ####
evi_lag1_60m_01_19_slope<-raster('evi_01_19_lag1_60m_slope.tif')
evi_lag1_36m_01_19_slope<-raster('evi_01_19_lag1_36m_slope.tif')
evi_lag1_84m_01_19_slope<-raster('evi_01_19_lag1_84m_slope.tif')

evi_lag1_slope_df<-as.data.frame(evi_lag1_60m_01_19_slope,xy=T)
colnames(evi_lag1_slope_df)[3]<-'slope'
evi_lag1_slope_df_nona<-evi_lag1_slope_df
evi_lag1_slope_df_nona<-evi_lag1_slope_df_nona[!is.na(evi_lag1_slope_df_nona$slope),]
#get percentage of positive trend values
sum(evi_lag1_slope_df_nona$slope>0)/nrow(evi_lag1_slope_df_nona)*100

#plot trend
vikO<-colour('vikO')
plot_tacslope<-ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = evi_lag1_slope_df,aes(x,y,fill=slope),
              show.legend = T) +
  scale_fill_scico(palette = 'vikO', direction = 1,
                   name='EVI TAC trend \n2001-2019',
                   na.value = 'transparent',
                   limits=c(-0.006,0.006)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

#histogram of trend values
plot_hist<-ggplot(data=evi_lag1_slope_df) + 
  geom_histogram(aes(x = slope),fill=vikO(30)) +
  xlab('EVI TAC trend 2001-2019') +
  ylab('Number of pixels') +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept=-0.0004458606,linetype='dashed',color="#3B3263")
#get mean slope to add vertical line in histogram
mean(evi_lag1_slope_df$slope,na.rm=T)

#Figure 1: Plot mean time series with random sample of 1000 pixels ####
evi_lag1_60m<-read.csv('evi_01_19_lag1_60m.csv')[,-1]
evi_lag1_84m<-read.csv('evi_01_19_lag1_84m.csv')[,-1]
evi_lag1_36m<-read.csv('evi_01_19_lag1_36m.csv')[,-1]

evi_lag1_nona<-evi_lag1_60m
evi_lag1_nona<-evi_lag1_nona[,!is.na(evi_lag1_nona[1,])]
evi_lag1_nona_df_mean<-rowMeans(evi_lag1_nona,na.rm=T)

set.seed(1)
evi_lag1_sample<-evi_lag1_nona[,sample(ncol(evi_lag1_nona),1000)]
evi_lag1_sample_t<-t(evi_lag1_sample)
rownames(evi_lag1_sample_t)<-paste0('pixel_',seq(1000))
colnames(evi_lag1_sample_t)<-paste0('time_',seq(as.Date('2005-12-01'),as.Date('2019-12-01'),by='month'))
evi_lag1_ts_sample_df<-as.data.frame(evi_lag1_sample_t)
evi_lag1_ts_sample_df$pixel<-rownames(evi_lag1_sample_t)
evi_lag1_ts_sample_melt<-melt(evi_lag1_ts_sample_df,id.vars='pixel')
evi_lag1_ts_sample_melt$time<-as.Date(gsub("time_","", evi_lag1_ts_sample_melt$variable))
#find timestep of minimal mean TAC value
which(evi_lag1_nona_df_mean==min(evi_lag1_nona_df_mean))

#create plot with density for overlapping time series
dates<-seq(as.Date('2005-12-01'),as.Date('2019-12-01'),by='month')
df_all<-evi_lag1_ts_sample_melt
x <- densCols(df_all$time,df_all$value, colramp=colorRampPalette(c("black", "white")))
df_all$'Line density'<-col2rgb(x)[1,] + 1L
bamako<-colour('bamako')

plot_meants<-ggplot()+
  geom_line(data=df_all[order(df_all$`Line density`),],aes(x = time,y = value, group = pixel,colour=`Line density`),
            size=0.5, show.legend = T) +
  #scale_colour_gradient(low=vikO(30)[14],high = vikO(30)[8]) +
  scale_colour_scico(palette = 'bamako') +
  geom_line(aes(x=dates,y=evi_lag1_nona_df_mean),col='black') +
  xlab('Time')+
  ylab('EVI TAC') +
  ylim(-0.5,0.8) +
  geom_vline(xintercept=dates[118],linetype='dashed',col='black',linewidth=0.5) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Figure 1: Combine plots figure 1 ####
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(3, 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot_tacslope, vp=define_region(1:2, 1:2))
print(plot_hist, vp = define_region(3, 1))
print(plot_meants, vp = define_region(3, 2))


#Figure 2: Plot interaction effects of model ####
tac_slope_nona<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]
seas<-raster('seasonality_index_80_19.tif')
interannvar<-raster('betweenyearvariation_80_19.tif')
seas_df<-raster::extract(seas, tac_slope_nona[,1:2])
interannvar_df<-raster::extract(interannvar, tac_slope_nona[,1:2])
droughtlegacy_df<-read.csv('droughtlegacy_01_19_nonaslope.csv')[,-1]
droughtlegacy_df[is.na(droughtlegacy_df)]<-0
droughtlegacy_clim_df<-cbind(droughtlegacy_df,seas_df,interannvar_df)
droughtlegacy_clim_df_st<-scale(droughtlegacy_clim_df[,3:8])

model_var_01_19<-cbind(tac_slope_nona,droughtlegacy_clim_df_st)
colnames(model_var_01_19)[c(3)]<-c('slope')

f <- as.formula(slope ~ d_amount_01_19 + d_avint_01_19 + d_avdur_01_19 
                + d_avint_01_19:d_amount_01_19 + d_avdur_01_19:d_amount_01_19 
                + w_totsev_01_19 + w_totsev_01_19:d_amount_01_19
                + seas_df + interannvar_df
                + d_amount_01_19:seas_df + d_amount_01_19:interannvar_df)
model_01_19_lm<-lm(f,data = model_var_01_19)

#amount * average intensity
ef2<-effect(term='d_amount_01_19*d_avint_01_19', mod=model_01_19_lm)
efdata2<-as.data.frame(ef2) #convert the effects list to a data frame
effect_avint<-ggplot(efdata2, aes(x=d_amount_01_19, y=fit, color=d_avint_01_19,group=d_avint_01_19)) + 
  geom_point() + 
  geom_line(size=1.2) +
  scale_colour_gradient(low= 'darksalmon',high='darkred',
                        breaks = c(0,2,4,6),
                        labels = c(0,2,4,6)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=d_avint_01_19),alpha=0.3) + 
  scale_fill_gradient(low= 'darksalmon',high='darkred',
                      breaks = c(0,2,4,6),
                      labels = c(0,2,4,6)) +
  labs(y="TAC trend", color="Average drought\nintensity in\n2001-2019", 
       fill="Average drought\nintensity in\n2001-2019") + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.002,0.001)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#amount * severity wet
ef3<-effect(term='d_amount_01_19*w_totsev_01_19', mod=model_01_19_lm)
efdata3<-as.data.frame(ef3) #convert the effects list to a data frame
effect_totsev<-ggplot(efdata3, aes(x=d_amount_01_19, y=fit, color=w_totsev_01_19,group=w_totsev_01_19)) + 
  geom_point() + 
  geom_line(size=1.2) +
  scale_colour_gradient(low= 'lightsteelblue',high='steelblue4',
                        breaks = c(0,25,50,75,100),
                        labels = c(0,25,50,75,100)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=w_totsev_01_19),alpha=0.3) + 
  scale_fill_gradient(low= 'lightsteelblue',high='steelblue4',
                      breaks = c(0,25,50,75,100),
                      labels = c(0,25,50,75,100)) +
  labs(y="TAC trend", color="Total wet\nseverity in\n2001-2019", 
       fill="Total wet\nseverity in\n2001-2019") + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.002,0.001)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#amount * seasonality
ef4<-effect(term='d_amount_01_19*seas_df', mod=model_01_19_lm)
efdata4<-as.data.frame(ef4) #convert the effects list to a data frame
effect_seas<-ggplot(efdata4, aes(x=d_amount_01_19, y=fit, color=seas_df,group=seas_df)) + 
  geom_point() + 
  geom_line(size=1.2) +
  scale_colour_gradient(low= 'palegoldenrod',high='darkorange2',limits=c(0,0.8),
                        breaks = c(0,0.2,0.4,0.6,0.8),
                        labels = c(0,0.2,0.4,0.6,0.8)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=seas_df),alpha=0.3) + 
  scale_fill_gradient(low= 'palegoldenrod',high='darkorange2',limits=c(0,0.8),
                      breaks = c(0,0.2,0.4,0.6,0.8),
                      labels = c(0,0.2,0.4,0.6,0.8)) +
  labs(y="TAC trend", color="Seasonality", 
       fill="Seasonality") + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.002,0.001)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#amount * interannual variability
ef1<-effect(term='d_amount_01_19*interannvar_df', mod=model_01_19_lm)
efdata1<-as.data.frame(ef1) #convert the effects list to a data frame
effect_interann<-ggplot(efdata1, aes(x=d_amount_01_19, y=fit, color=interannvar_df,group=interannvar_df)) + 
  geom_point() + 
  geom_line(size=1.2) +
  scale_colour_gradient(low= 'lemonchiffon',high='darkgreen',limits=c(0,0.4),
                        breaks = c(0,0.1,0.2,0.3,0.4),
                        labels = c(0,0.1,0.2,0.3,0.4)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=interannvar_df),alpha=0.3) + 
  scale_fill_gradient(low= 'lemonchiffon',high='darkgreen',limits=c(0,0.4),
                      breaks = c(0,0.1,0.2,0.3,0.4),
                      labels = c(0,0.1,0.2,0.3,0.4)) +
  labs(x= "Number of occurred droughts in 2001-2019", 
       y="TAC trend", color="Interannual\nvariability\nof precipitation", 
       fill="Interannual\nvariability\nof precipitation") + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.002,0.001)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#widths of plots the same, based on size of largest legend
g1 <- ggplotGrob(effect_avint)
g2 <- ggplotGrob(effect_totsev)
g3 <- ggplotGrob(effect_seas)
g4 <- ggplotGrob(effect_interann)
g2$widths <- g1$widths
g3$widths <- g1$widths
g4$widths <- g1$widths
#legends all at same horizontal position
leg1 <- convertX(sum(with(g1$grobs[[15]], grobs[[1]]$widths)), "mm")
leg2 <- convertX(sum(with(g2$grobs[[15]], grobs[[1]]$widths)), "mm")
leg3 <- convertX(sum(with(g3$grobs[[15]], grobs[[1]]$widths)), "mm")
leg4 <- convertX(sum(with(g4$grobs[[15]], grobs[[1]]$widths)), "mm")
#add empty column of width (difference with largest width) to right of legend box
g2$grobs[[15]] <- gtable_add_cols(g2$grobs[[15]], unit(abs(diff(c(leg1, leg2))), "mm"))
g3$grobs[[15]] <- gtable_add_cols(g3$grobs[[15]], unit(abs(diff(c(leg1, leg3))), "mm"))
g4$grobs[[15]] <- gtable_add_cols(g4$grobs[[15]], unit(abs(diff(c(leg1, leg4))), "mm"))

#Figure 2: Plot drivers of CSD and combine plots ####
evi_lag1_60m_01_19_slope<-raster('evi_01_19_lag1_60m_slope.tif')
ndroughts_01_19<-raster('droughtlegacy_01_19_ndroughts_all.tif')
ndroughts_01_19[is.na(ndroughts_01_19)]<-0
ndroughts_01_19<-mask(ndroughts_01_19,amazon)
ndroughts_01_19<-mask(ndroughts_01_19,evi_lag1_60m_01_19_slope)
avintdroughts_01_19<-raster('droughtlegacy_01_19_avintdroughts_all.tif')
avintdroughts_01_19[is.na(avintdroughts_01_19)]<-0
avintdroughts_01_19<-mask(avintdroughts_01_19,amazon)
avintdroughts_01_19<-mask(avintdroughts_01_19,evi_lag1_60m_01_19_slope)
avintdroughts_01_19[avintdroughts_01_19>6]<-6
avdurdroughts_01_19<-raster('droughtlegacy_01_19_avdurdroughts_all.tif')
avdurdroughts_01_19[is.na(avdurdroughts_01_19)]<-0
avdurdroughts_01_19<-mask(avdurdroughts_01_19,amazon)
avdurdroughts_01_19<-mask(avdurdroughts_01_19,evi_lag1_60m_01_19_slope)
avdurdroughts_01_19[avdurdroughts_01_19>20]<-20
totsevwets_01_19<-raster('droughtlegacy_01_19_totsevwets_all.tif')
totsevwets_01_19[is.na(totsevwets_01_19)]<-0
totsevwets_01_19<-mask(totsevwets_01_19,amazon)
totsevwets_01_19<-mask(totsevwets_01_19,evi_lag1_60m_01_19_slope)
totsevwets_01_19[totsevwets_01_19>100]<-100
seas_amazon<-raster::mask(seas,amazon)
interannvar_amazon<-raster::mask(interannvar,amazon)
seas_amazon<-raster::mask(seas_amazon,evi_lag1_60m_01_19_slope)
interannvar_amazon<-raster::mask(interannvar_amazon,evi_lag1_60m_01_19_slope)

avintdroughts_01_19_df<-as.data.frame(avintdroughts_01_19,xy=T)
plot_avint<-ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = avintdroughts_01_19_df,
              aes(x,y,fill=droughtlegacy_01_19_avintdroughts_all),show.legend=F) +
  scale_fill_gradient(low= 'darksalmon',high='darkred',
                      name='Average drought\nintensity in\n2001-2019',
                      na.value = 'transparent',
                      breaks = c(0,2,4,6),
                      labels = c(0,2,4,6)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(-80,-44),expand=c(0,0)) +
  xlab("")+
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

totsevwets_01_19_df<-as.data.frame(totsevwets_01_19,xy=T)
plot_totwetsev<-ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = totsevwets_01_19_df,
              aes(x,y,fill=droughtlegacy_01_19_totsevwets_all),show.legend=F) +
  scale_fill_gradient(low= 'lightsteelblue',high='steelblue4',name='Total wet\nseverity in\n2001-2019',
                      na.value = 'transparent',
                      breaks = c(0,25,50,75,100),
                      labels = c(0,25,50,75,100)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(-80,-44),expand=c(0,0)) +
  xlab("")+
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

seas_df<-as.data.frame(si_amazon,xy=T)
plot_seas<-ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = seas_df,aes(x,y,fill=seasonality_index_80_19),show.legend=F) +
  scale_fill_gradient(low= 'palegoldenrod',high='darkorange2',name='Seasonality',
                      na.value = 'transparent',
                      limits=c(0,0.8),
                      breaks = c(0,0.2,0.4,0.6,0.8),
                      labels = c(0,0.2,0.4,0.6,0.8)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(-80,-44),expand=c(0,0)) +
  xlab("")+
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

interann_df<-as.data.frame(interannvar_amazon,xy=T)
plot_interann<-ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = interann_df,
              aes(x,y,fill=interannvar_amazon),show.legend=F) +
  scale_fill_gradient(low= 'lemonchiffon',high='darkgreen',name='Interannual\nvariability\nof precipitation',
                      na.value = 'transparent',
                      limits=c(0,0.4),
                      breaks = c(0,0.1,0.2,0.3,0.4),
                      labels = c(0,0.1,0.2,0.3,0.4)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(-80,-44),expand=c(0,0)) +
  xlab("")+
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

#widths of plots the same, based on size of largest legend
g0 <- ggplotGrob(plot_ndroughts)
g10 <- ggplotGrob(plot_avdur)
g6 <- ggplotGrob(plot_avint)
g7 <- ggplotGrob(plot_totwetsev)
g8 <- ggplotGrob(plot_seas)
g9 <- ggplotGrob(plot_interann)
g0$widths <- g10$widths
g10$widths <- g10$widths
g6$widths <- g10$widths
g7$widths <- g10$widths
g8$widths <- g10$widths
g9$widths <- g10$widths
#legends all at same horizontal position
leg0 <- convertX(sum(with(g0$grobs[[15]], grobs[[1]]$widths)), "mm")
leg10 <- convertX(sum(with(g10$grobs[[15]], grobs[[1]]$widths)), "mm")
leg6 <- convertX(sum(with(g6$grobs[[15]], grobs[[1]]$widths)), "mm")
leg7 <- convertX(sum(with(g7$grobs[[15]], grobs[[1]]$widths)), "mm")
leg8 <- convertX(sum(with(g8$grobs[[15]], grobs[[1]]$widths)), "mm")
leg9 <- convertX(sum(with(g9$grobs[[15]], grobs[[1]]$widths)), "mm")
#add empty column of width (difference with largest width) to right of legend box
g0$grobs[[15]] <- gtable_add_cols(g0$grobs[[15]], unit(abs(diff(c(leg10, leg0))), "mm"))
g6$grobs[[15]] <- gtable_add_cols(g6$grobs[[15]], unit(abs(diff(c(leg10, leg6))), "mm"))
g7$grobs[[15]] <- gtable_add_cols(g7$grobs[[15]], unit(abs(diff(c(leg10, leg7))), "mm"))
g8$grobs[[15]] <- gtable_add_cols(g8$grobs[[15]], unit(abs(diff(c(leg10, leg8))), "mm"))
g9$grobs[[15]] <- gtable_add_cols(g9$grobs[[15]], unit(abs(diff(c(leg10, leg9))), "mm"))

#combine in figure 2
grid.newpage()
grid.arrange(g6, g1, g7, g2, g8, g3, g9, g4, nrow = 4)


#Figure 3 ####
#create polygons with regional extents
#change coordinates with 0.2° so lines don't overlap
nw_extent1<-extent(c(-79.4,-66.2,-2.8,8.66))
nw_p1<-as(nw_extent1,'SpatialPolygons')
crs(nw_p1)<-crs(evi)
nc_extent1<-extent(c(-65.8,-59.2,-2.8,8.66))
nc_p1<-as(nc_extent1,'SpatialPolygons')
crs(nc_p1)<-crs(evi)
ne_extent1<-extent(c(-58.8,-45,-2.8,8.66))
ne_p1<-as(ne_extent1,'SpatialPolygons')
crs(ne_p1)<-crs(evi)
sw_extent1<-extent(c(-79.4,-70.2,-18,-3.2))
sw_p1<-as(sw_extent1,'SpatialPolygons')
crs(sw_p1)<-crs(evi)
sc_extent1<-extent(c(-69.8,-64.2,-18,-3.2))
sc_p1<-as(sc_extent1,'SpatialPolygons')
crs(sc_p1)<-crs(evi)
se_extent1<-extent(c(-63.8,-45,-18,-3.2))
se_p1<-as(se_extent1,'SpatialPolygons')
crs(se_p1)<-crs(evi)

tac_slope_drivers_nw<-crop(amazon,nw_p1)
tac_slope_drivers_nw_df<-fortify(tac_slope_drivers_nw)
tac_slope_drivers_nc<-crop(amazon,nc_p1)
tac_slope_drivers_nc_df<-fortify(tac_slope_drivers_nc)
tac_slope_drivers_ne<-crop(amazon,ne_p1)
tac_slope_drivers_ne_df<-fortify(tac_slope_drivers_ne)
tac_slope_drivers_sw<-crop(amazon,sw_p1)
tac_slope_drivers_sw_df<-fortify(tac_slope_drivers_sw)
tac_slope_drivers_sc<-crop(amazon,sc_p1)
tac_slope_drivers_sc_df<-fortify(tac_slope_drivers_sc)
tac_slope_drivers_se<-crop(amazon,se_p1)
tac_slope_drivers_se_df<-fortify(tac_slope_drivers_se)
si_amazon<-raster('seasonality_index_80_19.tif')
si_amazon<-raster::mask(si_amazon,amazon)
si_amazon_df<-as.data.frame(si_amazon,xy=T)

ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey95',color='grey30') +
  geom_raster(data = si_amazon_df,aes(x,y,fill=seasonality_index_80_19)) +
  scale_fill_gradient(low = 'white',high = 'grey20',na.value = 'transparent',
                      name='Precipitation\nseasonality',
                      limits=c(0.2,0.9)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',size=0.3) +
  geom_polygon(data = tac_slope_drivers_nw,aes(x=long,y=lat,group=group),color='darkred',size=2,fill=NA) +
  geom_polygon_pattern(data = tac_slope_drivers_nc,aes(x=long,y=lat,group=group),color='grey20',size=2,fill=NA,
                       pattern='stripe',pattern_density=0.05) +
  geom_polygon(data = tac_slope_drivers_ne,aes(x=long,y=lat,group=group),color='steelblue4',size=2,fill=NA) +
  geom_polygon(data = tac_slope_drivers_sw,aes(x=long,y=lat,group=group),color='darkorange2',size=2,fill=NA) +
  geom_polygon_pattern(data = tac_slope_drivers_sc,aes(x=long,y=lat,group=group),color='grey20',size=2,fill=NA,
                       pattern='stripe',pattern_density=0.05) +
  geom_polygon(data = tac_slope_drivers_se,aes(x=long,y=lat,group=group),color='darkgreen',size=2,fill=NA) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 
