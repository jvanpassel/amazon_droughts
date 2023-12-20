library(ggplot2)
library(scales)
library(maptools)
library(scico)
library(khroma)
library(ggeffects)
library(sp)
library(precrec)


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
colnames(evi_lag1_sample_t)<-paste0('time_',seq(as.Date('2003-07-01'),as.Date('2017-07-01'),by='month'))
evi_lag1_ts_sample_df<-as.data.frame(evi_lag1_sample_t)
evi_lag1_ts_sample_df$pixel<-rownames(evi_lag1_sample_t)
evi_lag1_ts_sample_melt<-melt(evi_lag1_ts_sample_df,id.vars='pixel')
evi_lag1_ts_sample_melt$time<-as.Date(gsub("time_","", evi_lag1_ts_sample_melt$variable))
#find timestep of minimal mean TAC value
which(evi_lag1_nona_df_mean==min(evi_lag1_nona_df_mean))

#create plot with density for overlapping time series
dates<-seq(as.Date('2003-07-01'),as.Date('2017-07-01'),by='month')
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

#Figure 2: Plot drought-related effects of model ####
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

#average intensity
ef1<-effect(term='d_avint_01_19', mod=model_01_19_lm)
efdata1<-as.data.frame(ef1) #convert the effects list to a data frame
ggplot(efdata1,aes(x=d_avint_01_19, y=fit)) + 
  geom_point() + 
  geom_ribbon(mapping=aes(ymin=fit-se, ymax=fit+se),alpha=0.3,color='darkred') + 
  geom_line(size=1,color='darkred') +
  labs(y="TAC trend",x='Average drought intensity') + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.0008,0.0005)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#average duration
ef2<-effect(term='d_avdur_01_19', mod=model_01_19_lm)
efdata2<-as.data.frame(ef2) #convert the effects list to a data frame
ggplot(efdata2, aes(x=d_avdur_01_19, y=fit)) + 
  geom_point() + 
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=0.3, fill='darkorange2') + 
  geom_line(size=1,color='darkorange2') +
  labs(y="TAC trend",x='Average drought duration') + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.0008,0.0005)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#drought frequency
ef3<-effect(term='d_amount_01_19', mod=model_01_19_lm)
efdata3<-as.data.frame(ef3) #convert the effects list to a data frame
ggplot(efdata3, aes(x=d_amount_01_19, y=fit)) + 
  geom_point() + 
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=0.3, fill='darkgreen') + 
  geom_line(size=1,color='darkgreen',linetype='dashed') +
  labs(y="TAC trend",x='Drought frequency') + 
  scale_y_continuous(labels = scales::comma, limits = c(-0.0008,0.0005)) +
  theme_bw() + 
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")

#Figure 2: Plot drought-related drivers of CSD ####
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

ndroughts_01_19_df<-as.data.frame(ndroughts_01_19,xy=T)
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = ndroughts_01_19_df,
              aes(x,y,fill=droughtlegacy_01_19_ndroughts_all),show.legend=T) +
  scale_fill_gradient(low= 'lemonchiffon',high='darkgreen',
                      name='Droughts',
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
  scale_y_continuous(limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

avintdroughts_01_19_df<-as.data.frame(avintdroughts_01_19,xy=T)
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = avintdroughts_01_19_df,
              aes(x,y,fill=droughtlegacy_01_19_avintdroughts_all),show.legend=T) +
  scale_fill_gradient(low= 'darksalmon',high='darkred',
                      name='Drought',
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
  scale_y_continuous(limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

avdurdroughts_01_19_df<-as.data.frame(avdurdroughts_01_19,xy=T)
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey90',color='grey30') +
  geom_raster(data = avdurdroughts_01_19_df,
              aes(x,y,fill=droughtlegacy_01_19_avdurdroughts_all),show.legend=T) +
  scale_fill_gradient(low='palegoldenrod',high='darkorange2',
                      name='Drought',
                      na.value = 'transparent',
                      breaks = c(0,10,20,30),
                      labels = c(0,10,20,'30+')) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

#Figure 3 ####
#create smaller polygons by making 20 km buffer inside woody plant subregion polygons
shrinkIfPossible <- function(sf, size) {
  # compute inward buffer
  sg <- st_buffer(st_geometry(sf), -size)
  # update geometry only if polygon is not degenerate
  st_geometry(sf)[!st_is_empty(sg)] = sg[!st_is_empty(sg)]
  # return updated dataset
  return(sf)
}
pol_list<-list()
for(i in 1:13){ 
  #1. convert spfd to sf
  subregions_i <- st_as_sf(subregions[i,])
  #2. run function to create buffer inside polygon
  subregions_i_buff <- shrinkIfPossible(subregions_i, 20000) 
  #3. convert back to spdf
  subregions_i_buff<-as(subregions_i_buff,'Spatial')
  #4. transform to correct CRS
  subregions_i_buff_crs <-spTransform(subregions_i_buff,crs(amazon))
  i#5. fortify to plot
  subregions_crs_df_i<-fortify(subregions_i_buff_crs)
  pol_list[[i]]<-subregions_crs_df_i
}

subregions_crs_df_5<-fortify(subregions_crs[5,])
subregions_crs_df_7<-fortify(subregions_crs[7,])
subregions_crs_df_8<-fortify(subregions_crs[8,])
subregions_crs_df_11<-fortify(subregions_crs[11,])
subregions_crs_df_12<-fortify(subregions_crs[12,])
subregions_crs_df_13<-fortify(subregions_crs[13,])
subregions_crs_df<-fortify(subregions_crs)

si_amazon<-raster('seasonality_index_80_19.tif')
si_amazon<-raster::mask(si_amazon,amazon)
si_amazon_df<-as.data.frame(si_amazon,xy=T)
#Nothing = grey20; Frequency = darkgreen; Intensity = darkred; Duration = darkorange2
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey95',color='grey30') +
  geom_raster(data = si_amazon_df,aes(x,y,fill=seasonality_index_80_19)) +
  scale_fill_gradient(low = 'white',high = 'grey20',na.value = 'transparent',
                      name='Precipitation\nseasonality',
                      limits=c(0.2,0.9)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',size=0.3) +
  geom_polygon(data = pol_list[[5]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[7]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[8]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[11]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[12]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[13]],aes(x=long,y=lat,group=group),color='grey20',size=1.2,fill=NA) +
  geom_polygon(data = pol_list[[12]],aes(x=long,y=lat,group=group),color='darkgreen',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[6]],aes(x=long,y=lat,group=group),color='darkred',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[1]],aes(x=long,y=lat,group=group),color='darkgreen',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[1]],aes(x=long,y=lat,group=group),color='darkred',size=1.4,linetype='dashed',fill=NA) +
  geom_polygon(data = pol_list[[2]],aes(x=long,y=lat,group=group),color='darkgreen',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[2]],aes(x=long,y=lat,group=group),color='darkred',size=1.4,linetype='dashed',fill=NA) +
  geom_polygon(data = pol_list[[2]],aes(x=long,y=lat,group=group),color='darkorange2',size=1.4,linetype='dotted',fill=NA) +
  geom_polygon(data = pol_list[[3]],aes(x=long,y=lat,group=group),color='darkgreen',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[4]],aes(x=long,y=lat,group=group),color='darkgreen',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[4]],aes(x=long,y=lat,group=group),color='darkred',size=1.4,linetype='dashed',fill=NA) +
  geom_polygon(data = pol_list[[9]],aes(x=long,y=lat,group=group),color='darkorange2',size=1.4,fill=NA) +
  geom_polygon(data = pol_list[[10]],aes(x=long,y=lat,group=group),color='darkorange2',size=1.4,fill=NA) +
  geom_path(data = subregions_crs_df,aes(x=long,y=lat,group=group),color='black',size=0.3) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-20,10),expand=c(0,0)) +
  #guides(linetype=guide_legend("Region")) +
  coord_equal() 

#plot different subregions
palette_blues <- colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Blues'))(13)
ggplot() +  
  geom_polygon(data = subregions_crs_df,aes(x=long,y=lat,group=group,fill=id),color='black',size=0.3) +
  scale_fill_manual(values = palette_blues) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position = 'none') +
  coord_equal() 
