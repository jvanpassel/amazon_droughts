library(raster)
library(zoo)
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
evi<-brick('evi_resMODIS_01_19_burn_tree_mask.tif')
get_slope<-function(x){
  if(!all(is.na(x))){
    lm_slope<-lm(x~t)
    slope<-lm_slope$coefficients[[2]]
  } else {
    slope <- NA
  }
  return(slope)
}

#Create VOD dataset similar to Boulton et al, 2022####
amazon_wgs<-amazon
crs(amazon_wgs)<-'+proj=longlat +datum=WGS84 +no_defs'
#download ZIP from https://zenodo.org/record/2575599#.Yuj2hnZBzcs
#create monthly composites per year
rastlist<-list.files(path = 'VODCA/2000',pattern = '.nc',all.files = T,full.names = T)
raststack<-stack(rastlist)
raststack_am<-crop(raststack,amazon_wgs)
mi <- c(rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31)))
mi <- c(rep(1:12, c(31,29,31,30,31,30,31,31,30,31,30,31))) #leap year
length(mi)==nlayers(raststack_am)
raststack_am_monthlymean <- monthlyComposite(raststack_am, mi, fun = mean)
writeRaster(raststack_am_monthlymean,'VODCA/vodca_ku_2000.tif')

#Create broadleaf evergreen forest mask 
modis_lct1<-raster('MCD12Q1.006_LC_Type1_doy2001001_aid0001.tif')
#lct = 2 is broadleaf evergreen forest
modis_lct1_bef<-modis_lct1
modis_lct1_bef[modis_lct1_bef!=2]<-0
modis_lct1_bef[modis_lct1_bef!=0]<-1
modis_lct1_bef_fraction<-projectRaster(modis_lct1_bef,vodcastack,method='bilinear')
modis_lct1_bef_fraction_crop<-mask(modis_lct1_bef_fraction,amazon_wgs)
bef_mask<-modis_lct1_bef_fraction_crop
bef_mask[bef_mask<0.8]<-NA
bef_mask_df<-getValues(bef_mask)
#remove human land-use grid cells (lct=12,13,14)
modis_lct1_hlu<-modis_lct1
modis_lct1_hlu[modis_lct1_hlu==12]<-100
modis_lct1_hlu[modis_lct1_hlu==13]<-100
modis_lct1_hlu[modis_lct1_hlu==14]<-100
modis_lct1_hlu[modis_lct1_hlu!=100]<-0
modis_lct1_hlu[modis_lct1_hlu==100]<-1
modis_lct1_hlu_fraction<-projectRaster(modis_lct1_hlu,vodcastack,method = 'bilinear')
modis_lct1_hlu_fraction_crop<-mask(modis_lct1_hlu_fraction,amazon_wgs)
hlu_mask<-modis_lct1_hlu_fraction_crop
hlu_mask[hlu_mask>0]<-NA
bef_hlu_mask<-mask(bef_mask,hlu_mask)
bef_hlu_mask_df<-getValues(bef_hlu_mask)
writeRaster(bef_hlu_mask,'LCT1_BEF_80_mask.tif',overwrite=T)

#Combine yearly stacks and mask with broadleaf evergreen forest mask
vodcalist<-list.files(path = 'VODCA',pattern = '.tif$',all.files = T,full.names = T)
vodcastack<-stack(vodcalist)
vodcastack_crop<-mask(vodcastack,amazon_wgs)
bef_mask<-raster('LCT1_BEF_80_mask.tif')
vodcastack_mask<-mask(vodcastack_crop,bef_mask)
writeRaster(vodcastack_mask,'vodca_ku_88_16_befmask.tif',overwrite=T)

#Calculate VOD TAC and trend in TAC ####
vodcastack_mask<-stack('vodca_ku_88_16_befmask.tif')
vodca_df<-as.data.frame(vodcastack_mask,xy=T)
write.csv(vodca_df,'vodca_ku_88_16_befmask_df.tif')

#Detrend VOD time series
vodca_co<-vodca_df[,1:2]
vodca_val<-vodca_df[,3:350] #for 1988-2016
vodca_val<-vodca_df[,159:350] #for 2001-2016

vodca_nona_index<-which(!is.na(vodca_val[,1]))
write.csv(vodca_nona_index,'./lag1_ac/vodca_88_16_nona_index.csv')
vodca_na_index<-which(is.na(vodca_val[,1]))
vodca_val_nona<-vodca_val[vodca_nona_index,]
vodca_val_t<-as.data.frame(t(vodca_val_nona))
vodca_co_nona<-vodca_co[vodca_nona_index,]
#convert df into ts
dates<-seq(as.Date('2001-01-01'),as.Date('2016-12-01'),by='month')
mo <- as.numeric(format(dates[1], "%m"))
yr <- as.numeric(format(dates[1], "%Y"))
vodca_val_ts<-ts(vodca_val_t, start = c(yr, mo), freq = 12)
#apply STL decomposition
vodca_seas<-as.data.frame(matrix(nrow=nrow(vodca_val_ts),ncol=ncol(vodca_val_ts)))
vodca_trend<-as.data.frame(matrix(nrow=nrow(vodca_val_ts),ncol=ncol(vodca_val_ts)))
vodca_remain<-as.data.frame(matrix(nrow=nrow(vodca_val_ts),ncol=ncol(vodca_val_ts)))
for(i in 1:ncol(vodca_val_t)){
  vodca_val_ts_i<-vodca_val_ts[,i]
  if(all(!is.na(vodca_val_ts_i))){
    vodca_val_stl_i<-stl(vodca_val_ts_i,s.window = 'periodic',t.window = 19,l.window = 13)
    vodca_seas[,i]<-vodca_val_stl_i$time.series[,1]
    vodca_trend[,i]<-vodca_val_stl_i$time.series[,2]
    vodca_remain[,i]<-vodca_val_stl_i$time.series[,3]
  } else {
    vodca_seas[,i]<-rep(NA,192)
    vodca_trend[,i]<-rep(NA,192)
    vodca_remain[,i]<-rep(NA,192)
  }
  print(i)
}
write.csv(vodca_co_nona,'vodca_01_16_stl_co.csv')
write.csv(vodca_seas,'vodca_01_16_stl_seas.csv')
write.csv(vodca_trend,'vodca_01_16_stl_trend.csv')
write.csv(vodca_remain,'vodca_01_16_stl_remain.csv')

#Calculate lag-1 autocorrelation
vodca_co_nona<-read.csv('vodca_01_16_stl_co.csv')[,-1]
vodca_remain<-read.csv('vodca_01_16_stl_remain.csv')[,-1]

vodca_lag1_60m<-as.data.frame(matrix(nrow=133,ncol=ncol(vodca_remain))) 
for(i in 1:ncol(vodca_remain)){
  vodca_remain_i<-vodca_remain[,i]
  if(all(!is.na(vodca_remain_i))){
    TS <- rollapply(vodca_remain_i, width = 60, #width of rolling window in months (3,5,7 years)
                    FUN = function(z) acf(z,na.action=na.pass,lag.max= 1,plot=FALSE)$acf[2],
                    by.column = FALSE, align = "right")
    vodca_lag1_60m[,i]<-TS
  } else {
    vodca_lag1_60m[,i]<-rep(NA,133)
  }
  print(i)
}
write.csv(vodca_lag1_60m,'vodca_01_16_lag1_60m.csv')

#Create rasterstack of TAC
vodca_nona_index<-read.csv('vodca_01_16_nona_index.csv')[,-1]
vodca_lag1_60m<-read.csv('vodca_01_16_lag1_60m.csv')[,-1]
vodca_lag1_60m_df<-as.data.frame(matrix(nrow=133,ncol=14980))
vodca_lag1_60m_df[,vodca_nona_index]<-vodca_lag1_60m
vodca_lag1_60m_df<-t(vodca_lag1_60m_df)
vodca_lag1_60m_rasterstack <- brick(extent(vodcastack_mask), nl=133, nrows=nrow(vodcastack_mask),
                                    ncols=ncol(vodcastack_mask), crs=proj4string(vodcastack_mask))
for(k in 1:133){
  vodca_lag1_60m_rasterstack[[k]] <- vodca_lag1_60m_df[,k]
}
writeRaster(vodca_lag1_60m_rasterstack,'vodca_01_16_lag1_60m.tif',overwrite=T)
vodca_lag1_60m_amazon<-mask(vodca_lag1_60m_rasterstack,amazon)
writeRaster(vodca_lag1_60m_amazon,'vodca_01_16_lag1_60m.tif',overwrite=T)
vodca_lag1_60m_amazon_df<-as.data.frame(vodca_lag1_60m_rasterstack,xy=T)

#calculate TAC slope
t<-c(1:133)
vodca_lag1_60m_slope<-calc(vodca_lag1_60m_rasterstack,fun=get_slope)
writeRaster(vodca_lag1_60m_slope,'vodca_01_16_lag1_60m_slope.tif',overwrite=T)

#Compare EVI and VOD TAC 01-16 ####
evi_lag1_01_16_slope<-raster('evi_01_16_lag1_60m_slope.tif')
vod_lag1_01_16_slope<-raster('vodca_01_16_lag1_60m_slope.tif')
evi_lag1_01_16_slope_df<-as.data.frame(evi_lag1_01_16_slope,xy=T)

#Change resolution VOD to resolution EVI
vod_lag1_01_16_slope_resevi<-disaggregate(vod_lag1_01_16_slope,fact=5)
vod_lag1_01_16_slope_resevi<-projectRaster(vod_lag1_01_16_slope_resevi,evi_lag1_01_16_slope)
vod_lag1_01_16_slope_resevi_df<-as.data.frame(vod_lag1_01_16_slope_resevi,xy=T)

#Plot EVI vs VOD
evi_lag1_01_16_slope_df_nona<-which(!is.na(evi_lag1_01_16_slope_df$evi_01_16_lag1_60m_slope))
evi_lag1_01_16_slope_df_nona_val<-evi_lag1_01_16_slope_df[evi_lag1_01_16_slope_df_nona,]
vod_lag1_01_16_slope_resevi_df_nona_val<-vod_lag1_01_16_slope_resevi_df[evi_lag1_01_16_slope_df_nona,]
vod_lag1_01_16_slope_resevi_df_nona_val2<-vod_lag1_01_16_slope_resevi_df_nona_val[!is.na(vod_lag1_01_16_slope_resevi_df_nona_val$vodca_01_16_lag1_60m_slope),]
evi_lag1_01_16_slope_df_nona_val2<-evi_lag1_01_16_slope_df_nona_val[!is.na(vod_lag1_01_16_slope_resevi_df_nona_val$vodca_01_16_lag1_60m_slope),]
plot(evi_lag1_01_16_slope_df_nona_val2$evi_01_16_lag1_60m_slope,vod_lag1_01_16_slope_resevi_df_nona_val2$vodca_01_16_lag1_60m_slope)
df_all<-as.data.frame(cbind(evi_lag1_01_16_slope_df_nona_val2$evi_01_16_lag1_60m_slope,vod_lag1_01_16_slope_resevi_df_nona_val2$vodca_01_16_lag1_60m_slope))
colnames(df_all)<-c('EVI_slope','VOD_slope')
x <- densCols(df_all$EVI_slope,df_all$VOD_slope, colramp=colorRampPalette(c("black", "white")))
df_all$xdens<-col2rgb(x)[1,] + 1L
batlow<-colour('batlow')
cols <-  batlow(256)
df_all$xcol <- cols[df_all$xdens]
plot(df_all$EVI_slope~df_all$VOD_slope,data=df_all[order(df_all$xdens),],
     pch=20,col=df_all$xcol)
summary(lm(df_all$EVI_slope~df_all$VOD_slope))

bamako<-colour('bamako')
ggplot()+
  geom_point(data=df_all[order(df_all$xdens),],aes(x=VOD_slope,y=EVI_slope,colour=xdens),
             size=0.5, show.legend = FALSE) +
  scale_colour_scico(palette = 'bamako') +
  xlab('VOD TAC slope') +
  ylab('EVI TAC slope') +
  geom_abline(intercept=-5.313e-04,slope=3.552e-02,col="#001260",linetype='dashed') +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

evi_lag1_06_19_slope_df<-as.data.frame(vod_lag1_01_16_slope_resevi,xy=T)
colnames(evi_lag1_06_19_slope_df)[3]<-'slope'
evi_lag1_06_19_slope_df_nona<-evi_lag1_06_19_slope_df
evi_lag1_06_19_slope_df_nona<-evi_lag1_06_19_slope_df_nona[!is.na(evi_lag1_06_19_slope_df_nona$slope),]
ggplot() +  
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill='grey50',color='grey30') +
  geom_raster(data = evi_lag1_06_19_slope_df,aes(x,y,fill=slope)) +
  scale_fill_scico(palette = 'vik',name='VOD TAC trend \n2001-2016',
                   na.value = 'transparent',
                   limits=c(-0.006,0.006)) +
  geom_map(data = worldmap,map = worldmap,aes(x=long,y=lat,map_id=id),fill=NA,color='grey30') +
  geom_path(data = amazon_df,aes(x=long,y=lat,group=group),color='black',fill='white',size=0.3) +
  theme_bw() +
  scale_x_continuous(name=expression(paste("Longitude [deg]")),
                     limits=c(-80,-44),expand=c(0,0)) +
  scale_y_continuous(name=expression(paste("Latitude [deg]")),
                     limits=c(-19,10),expand=c(0,0)) +
  coord_equal() 

#Compare time series
vodca_lag1_60m_df<-read.csv('vodca_01_16_lag1_60m.csv')[,-1]
evi_lag1_60m_df<-read.csv('evi_01_16_lag1_60m.csv')[,-1]

evi_lag1_ts<-t(evi_lag1_60m_df)
evi_lag1_60m_df_mean<-colMeans(evi_lag1_ts,na.rm=T)
evi_lag1_60m_df_sd<-apply(evi_lag1_ts,2,function(x){sd(x,na.rm=T)})
evi_lag1_60m_df_meanplussd<-evi_lag1_60m_df_mean+evi_lag1_60m_df_sd
evi_lag1_60m_df_meanminussd<-evi_lag1_60m_df_mean-evi_lag1_60m_df_sd
vod_lag1_ts<-vodca_lag1_60m_df[,-c(1:2)]
vod_lag1_60m_df_mean<-as.numeric(colMeans(vod_lag1_ts,na.rm=T))
vod_lag1_60m_df_sd<-as.numeric(apply(vod_lag1_ts,2,function(x){sd(x,na.rm=T)}))
vod_lag1_60m_df_meanplussd<-vod_lag1_60m_df_mean+vod_lag1_60m_df_sd
vod_lag1_60m_df_meanminussd<-vod_lag1_60m_df_mean-vod_lag1_60m_df_sd

vikO<-colour('vikO')
vikO(7)(1)
vik[7](7)
dates<-seq(as.Date('2005-12-01'),as.Date('2016-12-01'),by='month')
ggplot()+
  geom_line(aes(x=dates,y=evi_lag1_60m_df_mean,col='EVI'))+
  geom_ribbon(aes(x=dates,y=evi_lag1_60m_df_mean, ymin=evi_lag1_60m_df_meanminussd,ymax=evi_lag1_60m_df_meanplussd),alpha=0.2,fill="#334B7F")+
  geom_line(aes(x=dates,y=vod_lag1_60m_df_mean,col='VOD'))+
  geom_ribbon(aes(x=dates,y=vod_lag1_60m_df_mean, ymin=vod_lag1_60m_df_meanminussd,ymax=vod_lag1_60m_df_meanplussd),alpha=0.2,fill="#7B261E")+
  xlab('Time')+
  ylab('mean TAC') +
  scale_colour_manual('',breaks = c('EVI','VOD'),
                      values = c("#334B7F","#7B261E")) +
  geom_vline(xintercept=dates[118],linetype='dashed',color="#334B7F") +
  geom_vline(xintercept=dates[106],linetype='dashed',color="#7B261E") +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dates<-seq(as.Date('2005-12-01'),as.Date('2016-12-01'),by='month')
which(evi_lag1_60m_df_mean==min(evi_lag1_60m_df_mean))
which(vod_lag1_60m_df_mean==min(vod_lag1_60m_df_mean))
dates[113]
#September 2014 for VOD TAC, September 2015 for EVI TAC

#Create plot of mean time series and random sample of 1000 pixels
vodca_lag1_60m_df<-read.csv('vodca_01_16_lag1_60m.csv')[,-1]
evi_lag1_60m_df<-read.csv('evi_01_16_lag1_60m.csv')[,-1]
evi_lag1_60m_nona<-evi_lag1_60m_df
evi_lag1_60m_nona<-evi_lag1_60m_nona[,!is.na(evi_lag1_60m_nona[1,])]
evi_lag1_60m_df_mean<-rowMeans(evi_lag1_60m_nona,na.rm=T)
set.seed(1)
evi_lag1_sample<-evi_lag1_60m_nona[,sample(ncol(evi_lag1_60m_nona),1000)]
evi_lag1_sample_t<-t(evi_lag1_sample)
rownames(evi_lag1_sample_t)<-paste0('pixel_',seq(1000))
colnames(evi_lag1_sample_t)<-paste0('time_',seq(as.Date('2005-12-01'),as.Date('2016-12-01'),by='month'))
evi_lag1_ts_sample_df<-as.data.frame(evi_lag1_sample_t)
evi_lag1_ts_sample_df$pixel<-rownames(evi_lag1_sample_t)
evi_lag1_ts_sample_melt<-melt(evi_lag1_ts_sample_df,id.vars='pixel')
evi_lag1_ts_sample_melt$time<-as.Date(gsub("time_","", evi_lag1_ts_sample_melt$variable))
which(evi_lag1_60m_df_mean==min(evi_lag1_60m_df_mean))

vod_lag1_60m_nona<-vodca_lag1_60m_df[,3:135]
vod_lag1_60m_nona<-vod_lag1_60m_nona[!is.na(vod_lag1_60m_nona[,1]),]
vod_lag1_60m_df_mean<-colMeans(vod_lag1_60m_nona,na.rm=T)
set.seed(1)
vod_lag1_sample_t<-vod_lag1_60m_nona[sample(nrow(vod_lag1_60m_nona),1000),]
rownames(vod_lag1_sample_t)<-paste0('pixel_',seq(1000))
colnames(vod_lag1_sample_t)<-paste0('time_',seq(as.Date('2005-12-01'),as.Date('2016-12-01'),by='month'))
vod_lag1_ts_sample_df<-as.data.frame(vod_lag1_sample_t)
vod_lag1_ts_sample_df$pixel<-rownames(vod_lag1_sample_t)
vod_lag1_ts_sample_melt<-melt(vod_lag1_ts_sample_df,id.vars='pixel')
vod_lag1_ts_sample_melt$time<-as.Date(gsub("time_","", vod_lag1_ts_sample_melt$variable))
which(vod_lag1_60m_df_mean==min(vod_lag1_60m_df_mean))

dates<-seq(as.Date('2005-12-01'),as.Date('2016-12-01'),by='month')
ggplot()+
  geom_line(data = evi_lag1_ts_sample_melt,
            aes(x = time,y = value, group = pixel),linewidth=0.1, alpha=0.01,
            col="#2D7CA5") +
  geom_line(data = vod_lag1_ts_sample_melt,
            aes(x = time,y = value, group = pixel),linewidth=0.1, alpha=0.01,
            col="#AF893D") +
  geom_line(aes(x=dates,y=evi_lag1_60m_df_mean,col="EVI")) +
  geom_line(aes(x=dates,y=vod_lag1_60m_df_mean,col="VOD")) +
  scale_colour_manual('',breaks = c("EVI","VOD"),
                      values = c("#001260","#601200")) +
  xlab('Time')+
  ylab('mean TAC') +
  ylim(-0.4,0.6) +
  geom_vline(xintercept=dates[118],linetype='dashed',col="#001260",linewidth=0.5) +
  geom_vline(xintercept=dates[106],linetype='dashed',col="#601200",linewidth=0.5) +
  theme_classic()


