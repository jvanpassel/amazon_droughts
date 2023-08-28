library(raster)
library(zoo)

#load input data
evi<-brick('evi_resMODIS_01_19_burn_tree_mask.tif')
kndvi<-brick('kndvi_resMODIS_01_19_burn_tree_mask.tif')
get_slope<-function(x){
  if(!all(is.na(x))){
    lm_slope<-lm(x~t)
    slope<-lm_slope$coefficients[[2]]
  } else {
    slope <- NA
  }
  return(slope)
}

#Detrend EVI/kNDVI time series ####
vi<-evi #change to evi/kndvi
vi_df<-as.data.frame(vi,xy=T)
vi_co<-vi_df[,1:2]
vi_val<-vi_df[,3:230]
nona_index<-which(!is.na(vi_val[,1]))
na_index<-which(is.na(vi_val[,1]))
vi_val_nona<-vi_val[nona_index,]
vi_val_t<-as.data.frame(t(vi_val_nona))
vi_co_nona<-vi_co[nona_index,]
#convert dataframe into time series
dates<-seq(as.Date('2001-01-01'),as.Date('2019-12-01'),by='month')
mo <- as.numeric(format(dates[1], "%m"))
yr <- as.numeric(format(dates[1], "%Y"))
vi_val_ts<-ts(vi_val_t, start = c(yr, mo), freq = 12)
#apply STL decomposition
vi_seas<-as.data.frame(matrix(nrow=nrow(vi_val_ts),ncol=ncol(vi_val_ts)))
vi_trend<-as.data.frame(matrix(nrow=nrow(vi_val_ts),ncol=ncol(vi_val_ts)))
vi_remain<-as.data.frame(matrix(nrow=nrow(vi_val_ts),ncol=ncol(vi_val_ts)))

for(i in 1:ncol(vi_val_t)){
  vi_val_ts_i<-vi_val_ts[,i]
  if(all(!is.na(vi_val_ts_i))){
    vi_val_stl_i<-stl(vi_val_ts_i,s.window = 13,t.window = 19,l.window = 13)
    vi_seas[,i]<-vi_val_stl_i$time.series[,1]
    vi_trend[,i]<-vi_val_stl_i$time.series[,2]
    vi_remain[,i]<-vi_val_stl_i$time.series[,3]
  } else {
    vi_seas[,i]<-rep(NA,228)
    vi_trend[,i]<-rep(NA,228)
    vi_remain[,i]<-rep(NA,228)
  }
  print(i)
}
write.csv(vi_co_nona,'evi_01_19_stl_co.csv')
write.csv(vi_seas,'evi_01_10_stl_seas.csv')
write.csv(vi_trend,'evi_01_19_stl_trend.csv')
write.csv(vi_remain,'evi_01_19_stl_remain.csv')

#Calculate lag-1 autocorrelation from detrended EVI/kNDVI remainder ####
#EVI
vi_co_nona<-read.csv('evi_01_19_stl_co.csv')[,-1]
vi_remain<-read.csv('evi_01_19_stl_remain.csv')[,-1]
#kNDVI
vi_co_nona<-read.csv('kndvi_01_19_stl_co.csv')[,-1]
vi_remain<-read.csv('kndvi_01_19_stl_remain.csv')[,-1]

#with moving window of 3 years/36 months
vi_lag1_36m<-as.data.frame(matrix(nrow=193,ncol=ncol(vi_remain)))
for(i in 1:ncol(vi_remain)){
  vi_remain_i<-vi_remain[,i]
  if(all(!is.na(vi_remain_i))){
    TS <- rollapply(vi_remain_i, width = 36, #width of rolling window in months (3,5,7 years)
                    FUN = function(z) acf(z,na.action=na.pass,lag.max= 1,plot=FALSE)$acf[2],
                    by.column = FALSE, align = "right")
    vi_lag1_36m[,i]<-TS
  } else {
    vi_lag1_36m[,i]<-rep(NA,193)
  }
  print(i)
}
write.csv(vi_lag1_36m,'evi_01_19_lag1_36m.csv')

#with moving window of 5 years/60 months
vi_lag1_60m<-as.data.frame(matrix(nrow=169,ncol=ncol(vi_remain))) 
for(i in 1:ncol(vi_remain)){
  vi_remain_i<-vi_remain[,i]
  if(all(!is.na(evi_remain_i))){
    TS <- rollapply(vi_remain_i, width = 60, #width of rolling window in months (3,5,7 years)
                    FUN = function(z) acf(z,na.action=na.pass,lag.max= 1,plot=FALSE)$acf[2],
                    by.column = FALSE, align = "right")
    vi_lag1_60m[,i]<-TS
  } else {
    vi_lag1_60m[,i]<-rep(NA,169)
  }
  print(i)
}
write.csv(vi_lag1_60m,'evi_01_19_lag1_60m.csv')

#with moving window of 7 years/84 months
evi_lag1_84m<-as.data.frame(matrix(nrow=145,ncol=ncol(evi_remain)))
for(i in 1:ncol(vi_remain)){
  vi_remain_i<-vi_remain[,i]
  if(all(!is.na(vi_remain_i))){
    TS <- rollapply(vi_remain_i, width = 84, #width of rolling window in months (3,5,7 years)
                    FUN = function(z) acf(z,na.action=na.pass,lag.max= 1,plot=FALSE)$acf[2],
                    by.column = FALSE, align = "right")
    vi_lag1_84m[,i]<-TS
  } else {
    vi_lag1_84m[,i]<-rep(NA,145)
  }
  print(i)
}
write.csv(vi_lag1_84m,'evi_01_19_lag1_84m.csv')


#Add NA rows back into dataframes and convert into EVI/kNDVI TAC rasterstack and linear trend ####
#EVI
nona_index<-read.csv('evi_01_19_nona_index.csv')[,-1]
vi_lag1_36m<-read.csv('evi_01_19_lag1_36m.csv')[,-1]
vi_lag1_60m<-read.csv('evi_01_19_lag1_60m.csv')[,-1]
vi_lag1_84m<-read.csv('evi_01_19_lag1_84m.csv')[,-1]

#kNDVI
nona_index<-read.csv('kndvi_01_19_nona_index.csv')[,-1]
vi_lag1_60m<-read.csv('kndvi_01_19_lag1_60m.csv')[,-1]

#for moving window of 36 months
vi_lag1_36m_df<-as.data.frame(matrix(nrow=193,ncol=372732))
vi_lag1_36m_df[,nona_index]<-vi_lag1_36m
vi_lag1_36m_df<-t(vi_lag1_36m_df)
vi_lag1_36m_rasterstack <- brick(extent(evi), nl=193, nrows=nrow(evi),
                                  ncols=ncol(evi), crs=proj4string(evi))
for(k in 1:193){
  vi_lag1_36m_rasterstack[[k]] <- vi_lag1_36m_df[,k]
}
writeRaster(vi_lag1_36m_rasterstack,'evi_01_19_lag1_36m.tif')

t<-c(1:193)
vi_lag1_36m_01_19_slope<-calc(vi_lag1_36m_rasterstack,fun=get_slope)
writeRaster(vi_lag1_36m_01_19_slope,'evi_01_19_lag1_36m_slope.tif')

#for moving window of 60 months
vi_lag1_60m_df<-as.data.frame(matrix(nrow=169,ncol=372732))
vi_lag1_60m_df[,nona_index]<-vi_lag1_60m
vi_lag1_60m_df<-t(vi_lag1_60m_df)
vi_lag1_60m_rasterstack <- brick(extent(evi), nl=169, nrows=nrow(evi),
                                  ncols=ncol(evi), crs=proj4string(evi))
for(k in 1:169){
  vi_lag1_60m_rasterstack[[k]] <- vi_lag1_60m_df[,k]
}
writeRaster(vi_lag1_60m_rasterstack,'evi_01_19_lag1_60m.tif',overwrite=T)

t<-c(1:169)
vi_lag1_60m_01_19_slope<-calc(vi_lag1_60m_rasterstack,fun=get_slope)
writeRaster(vi_lag1_60m_01_19_slope,'evi_01_19_lag1_60m_slope.tif',overwrite=T)

#for moving window of 84 months
vi_lag1_84m_df<-as.data.frame(matrix(nrow=145,ncol=372732))
vi_lag1_84m_df[,nona_index]<-vi_lag1_84m
vi_lag1_84m_df<-t(vi_lag1_84m_df)
vi_lag1_84m_rasterstack <- brick(extent(evi), nl=145, nrows=nrow(evi),
                                  ncols=ncol(evi), crs=proj4string(evi))
for(k in 1:145){
  vi_lag1_84m_rasterstack[[k]] <- vi_lag1_84m_df[,k]
}
writeRaster(vi_lag1_84m_rasterstack,'evi_01_19_lag1_84m.tif')

t<-c(1:145)
vi_lag1_84m_01_19_slope<-calc(vi_lag1_84m_rasterstack,fun=get_slope)
writeRaster(vi_lag1_84m_01_19_slope,'evi_01_19_lag1_84m_slope.tif')

