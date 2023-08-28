library(raster)
library(zoo)
library(ggplot2)

amazon<-shapefile('amapoly_ivb.shp')
evi<-brick('evi_resMODIS_01_19_burn_tree_mask.tif')

#nir is the rasterstack of the monthly mean of the near-infrared MODIS MCD43C4 bands 
#red is the rasterstack of the monthly mean of the red MODIS MCD43C4 bands 

# Calculate kNDVI from NIR and red bands ####
# following https://github.com/IPL-UV/kNDVI
sigma1 <- (nir + red)/2
sigma2 <- calc(sigma1,fun=function(x) median(x,na.rm=T))
knr <- exp(-(nir-red)^2/(2*sigma2^2))
kndvi <- (1-knr) / (1+knr)

# Mask pixels with more than 4 consecutive missing values:
getna<-is.na(kndvi) #FALSE means non-NA, TRUE means NA
fn<-function(a) {
  seq <- rle(a)
  n=length(seq$lengths[seq$lengths>4 & seq$values==TRUE])
  return(n)
}
kndvi_consecutive_na<-calc(getna,fn)
kndvi_consecutive_na[kndvi_consecutive_na > 0] <-NA
kndvi_mask<-mask(kndvi,kndvi_consecutive_na)

# Fill gaps and remove outliers
remove_outliers <- function(x){
  if(!all(is.na(x))){
    x_initial <- ts(as.numeric(as.vector(x)), start = c(2001, 1), freq = 12)
    x_approx <- na.approx(x_initial, rule = 2,maxgap=4)
    ts_approx1<-x_approx
    ts_approx2<-x_approx
    
    window_width<-round(2*12/7,digits=0) #half window width is number of values per year divided by 7
    cutoff<-sd(x_approx)*3
    
    median_ts<-rollmedian(x_approx,k=window_width, fill = c(first(x_approx),NA,last(x_approx)))
    y1<-median_ts+cutoff
    y2<-median_ts-cutoff
    
    # Data value is defined as an outlier following two criteria: (1) it deviates more than the cutoff value
    # from the median in a moving window, and (2) it is lower than the mean value of its immediate neighbors 
    # minus the cutoff or it is larger than the highest value of its immediate neighbor plus the cutoff.
    
    for(i in 2:length(x_approx)-1){
      is.na(ts_approx1[i])<- x_approx[i] < y2[i] | x_approx[i] > y1[i] #criteria 1
      is.na(ts_approx2[i])<- x_approx[i] < (mean(x_approx[i-1],x_approx[i+1]) - cutoff) | x_approx[i] > (mean(x_approx[i-1],x_approx[i+1]) + cutoff) #criteria 2
    }
    
    crit1<-which(is.na(ts_approx1))
    crit2<-which(is.na(ts_approx2))
    outliers_index<-Reduce(intersect, list(crit1,crit2))
    x_approx[outliers_index]<-NA
    x_approx_final <- na.approx(x_approx, rule = 2) #peeks replaced by NA instead of interpolation
    return(x_approx_final)
  } else {
    return(x)
  }
}
kndvi_clean <- calc(kndvi_mask, fun=function(x){
  y<-as.data.frame(x)
  res <- apply(y, 2, remove_outliers)
  return(res)
})

# Use same mask as EVI
kndvi_clean<-mask(kndvi_clean, evi[[1]])

writeRaster(kndvi_clean,'kndvi_resMODIS_01_19_burn_tree_mask.tif')

# Detrend kNDVI and calculate TAC trend script_paper_calc_tac ####
# Compare EVI and kNDVI ####
# Compare mean time series
evi_lag1_60m<-read.csv('evi_01_19_lag1_60m.csv')[,-1]
kndvi_lag1_60m<-read.csv('kndvi_01_19_lag1_60m.csv')[,-1]

evi_lag1_ts<-t(evi_lag1_60m)
evi_lag1_60m_df_mean<-colMeans(evi_lag1_ts,na.rm=T)
evi_lag1_60m_df_sd<-apply(evi_lag1_ts,2,function(x){sd(x,na.rm=T)})
evi_lag1_60m_df_meanplussd<-evi_lag1_60m_df_mean+evi_lag1_60m_df_sd
evi_lag1_60m_df_meanminussd<-evi_lag1_60m_df_mean-evi_lag1_60m_df_sd

kndvi_lag1_ts<-t(kndvi_lag1_60m)
kndvi_lag1_60m_df_mean<-colMeans(kndvi_lag1_ts,na.rm=T)
kndvi_lag1_60m_df_sd<-apply(kndvi_lag1_ts,2,function(x){sd(x,na.rm=T)})
kndvi_lag1_60m_df_meanplussd<-kndvi_lag1_60m_df_mean+kndvi_lag1_60m_df_sd
kndvi_lag1_60m_df_meanminussd<-kndvi_lag1_60m_df_mean-kndvi_lag1_60m_df_sd

dates<-seq(as.Date('2005-12-01'),as.Date('2019-12-01'),by='month')
ggplot()+
  geom_line(aes(x=dates,y=evi_lag1_60m_df_mean,col='EVI'))+
  geom_ribbon(aes(x=dates,y=evi_lag1_60m_df_mean, ymin=evi_lag1_60m_df_meanminussd,ymax=evi_lag1_60m_df_meanplussd),alpha=0.2,fill="#334B7F")+
  geom_line(aes(x=dates,y=kndvi_lag1_mean,col='kNDVI'))+
  geom_ribbon(aes(x=dates,y=kndvi_lag1_mean, ymin=kndvi_lag1_60m_df_meanminussd,ymax=kndvi_lag1_60m_df_meanplussd),alpha=0.2,fill="#7B261E")+
  xlab('Time')+
  ylab('mean TAC') +
  scale_colour_manual('',breaks = c('EVI','kNDVI'),
                      values = c("#334B7F","#7B261E")) +
  geom_vline(xintercept=dates[118],linetype='dashed',color="black") +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none')

# Compare trends
evi_lag1_60m_slope<-raster('evi_01_19_lag1_60m_slope.tif')
kndvi_lag1_60m_slope<-raster('kndvi_01_19_lag1_60m_slope.tif')

evi_slope_df<-as.data.frame(evi_lag1_60m_slope,xy=T)
kndvi_slope_df<-as.data.frame(kndvi_lag1_60m_slope,xy=T)
colnames(evi_slope_df)[3]<-'slope'
colnames(kndvi_slope_df)[3]<-'slope'

df_all<-as.data.frame(cbind(evi_slope_df$slope,kndvi_slope_df$slope))
colnames(df_all)<-c('EVI_slope','kNDVI_slope')
x <- densCols(df_all$EVI_slope,df_all$kNDVI_slope, colramp=colorRampPalette(c("black", "white")))
df_all$xdens<-col2rgb(x)[1,] + 1L
df_all$xcol <- cols[df_all$xdens]
summary(lm(df_all$EVI_slope~df_all$kNDVI_slope))
bamako<-colour('bamako')
ggplot()+
  geom_point(data=df_all[order(df_all$xdens),],aes(x=kNDVI_slope,y=EVI_slope,colour=xdens),
             size=0.5, show.legend = F) +
  scale_colour_scico(palette = 'bamako') +
  xlab('kNDVI TAC slope') +
  ylab('EVI TAC slope') +
  geom_abline(intercept=-6.861e-05,slope=9.069e-01,col="#001260",linetype='dashed') +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

