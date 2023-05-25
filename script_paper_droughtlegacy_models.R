library(raster)
library(ggeffects)

#Load input data
amazon<-shapefile('amapoly_ivb.shp')

#Create drought legacy variables for all pixels ####
droughts <- brick('cwd_80_19_anomalies_droughts.tif')
wets <- brick('cwe_80_19_anomalies_wets.tif')
droughts_01_19<-droughts[[253:480]]
wets_01_19<-wets[[253:480]]
droughts_01_19_am<-mask(droughts_01_19,amazon)
wets_01_19_am<-mask(wets_01_19,amazon)

getndroughts<-function(x){
  if(!all(is.na(x))){
    y<-which(!is.na(x))
    droughts<-split(y, cumsum(c(1, diff(y) != 1)))
    n_droughts<-length(droughts)
    return(n_droughts)
  } else {
    return(NA)
  }
}
ndroughts_01_19<-calc(droughts_01_19,getndroughts)
ndroughts_01_19_am<-mask(ndroughts_01_19,amazon)
writeRaster(ndroughts_01_19_am,'droughtlegacy_01_19_ndroughts_all.tif',overwrite=T)

getavdurdroughts<-function(x){
  if(!all(is.na(x))){
    y<-which(!is.na(x))
    droughts<-split(y, cumsum(c(1, diff(y) != 1)))
    duration<-lapply(droughts,length)
    avdur_droughts<-mean(unlist(duration))
    return(avdur_droughts)
  } else {
    return(NA)
  }
}
avdurdroughts_01_19<-calc(droughts_01_19_am,getavdurdroughts)
writeRaster(avdurdroughts_01_19,'droughtlegacy_01_19_avdurdroughts_all.tif',overwrite=T)

getavintdroughts<-function(x){
  if(!all(is.na(x))){
    y<-which(!is.na(x))
    yval<-na.omit(x)
    droughts<-split(y, cumsum(c(1, diff(y) != 1)))
    duration<-lapply(droughts,length)
    droughtval<-split(yval, rep(1:length(droughts),as.vector(unlist(duration))))
    intensity<-lapply(droughtval,function(x){abs(min(x))})
    avint_droughts<-mean(unlist(intensity))
    return(avint_droughts)
  } else {
    return(NA)
  }
}
avintdroughts_01_19<-calc(droughts_01_19_am,getavintdroughts)
writeRaster(avintdroughts_01_19,'droughtlegacy_01_19_avintdroughts_all.tif',overwrite=T)

gettotsevwets<-function(x){
  if(!all(is.na(x))){
    y<-which(!is.na(x))
    yval<-na.omit(x)
    wets<-split(y, cumsum(c(1, diff(y) != 1)))
    duration<-lapply(wets,length)
    wetval<-split(yval, rep(1:length(wets),as.vector(unlist(duration))))
    severity<-lapply(wetval,function(x){abs(sum(x))})
    totsev_wets<-sum(unlist(severity))
    return(totsev_wets)
  } else {
    return(NA)
  }
}
totsevwets_01_19<-calc(wets_01_19_am,gettotsevwets)
writeRaster(totsevwets_01_19,'droughtlegacy_01_19_totsevwets_all.tif',overwrite=T)

#Create drought legacy variables for all non-NA pixels in TAC slope ####
#load TAC files
evi_lag1_60m_01_19<-brick('evi_01_19_lag1_60m.tif')
evi_lag1_60m_01_19_slope<-brick('evi_01_19_lag1_60m_slope.tif')
evi_lag1_60m_01_19_slope_df<-as.data.frame(evi_lag1_60m_01_19_slope,xy=T)
slope_nona<-evi_lag1_60m_01_19_slope_df[!is.na(evi_lag1_60m_01_19_slope_df$evi_lag1_60m_01_19_slope),]
co_nona<-slope_nona[,1:2]
nona_index<-which(!is.na(evi_lag1_60m_01_19_slope_df$evi_lag1_60m_01_19_slope))

#load droughts and wets files (1980-2019)
droughts <- brick('cwd_80_19_anomalies_droughts.tif')
wets <- brick('cwe_80_19_anomalies_wets.tif')

droughts_nona<-as.data.frame(raster::extract(droughts,co_nona))
wets_nona<-as.data.frame(raster::extract(wets,co_nona))

#For period 2001-2019
droughts_01_19<-droughts_nona[,c(253:480)]
wets_01_19<-wets_nona[,c(253:480)]
#For extra 5 years of drought compared to 2001-2019
droughts_96_19<-droughts_nona[,c(193:480)]
wets_96_19<-wets_nona[,c(193:480)]
#For extra 10 years of drought compared to 2001-2019
droughts_91_19<-droughts_nona[,c(133:480)]
wets_91_19<-wets_nona[,c(133:480)]
#For extra 15 years of drought compared to 2001-2019
droughts_86_19<-droughts_nona[,c(73:480)]
wets_86_19<-wets_nona[,c(73:480)]

#Run for loop
#decide on period 
drought_period<-droughts_01_19
wet_period<-wets_01_19

droughtlegacy_char<-as.data.frame(matrix(nrow=nrow(co_nona),ncol = 4))
names(droughtlegacy_char)<-c('d_amount','d_avdur','d_avint','w_totsev')
for(i in 1:nrow(co_nona)){
  #Define drought periods and derive characteristics
  drought_i<-as.numeric(drought_period[i,])
  if(!all(is.na(drought_i))){
    #drought periods
    NonNAindex <- which(!is.na(drought_i)) 
    droughtval<-na.omit(drought_i)
    diffdroughts<-split(NonNAindex, cumsum(c(1, diff(NonNAindex) != 1)))
    droughtstart<-lapply(diffdroughts, min)
    droughtend<-lapply(diffdroughts, max)
    droughtend<-c(1,droughtend)
    duration<-lapply(diffdroughts,length)
    diffdroughtval<-split(droughtval,rep(1:length(diffdroughts),as.vector(unlist(duration))))
    maxintensity<-lapply(diffdroughtval,function(x){abs(min(x))})
    
    #Amount of droughts
    droughtlegacy_char[i,1]<-as.numeric(length(diffdroughts))
    #Average drought characterization
    droughtlegacy_char[i,2]<-mean(unlist(duration))
    droughtlegacy_char[i,3]<-mean(unlist(maxintensity))
  } 
  #Define wet periods and derive characteristics
  wet_i<-as.numeric(wet_period[i,])
  if(!all(is.na(wet_i))){
    #drought periods
    NonNAwetindex <- which(!is.na(wet_i)) 
    wetval<-na.omit(wet_i)
    diffwets<-split(NonNAwetindex, cumsum(c(1, diff(NonNAwetindex) != 1)))
    wetstart<-lapply(diffwets, min)
    wetend<-lapply(diffwets, max)
    wetend<-c(1,wetend)
    wetduration<-lapply(diffwets,length)
    diffwetval<-split(wetval,rep(1:length(diffwets),as.vector(unlist(wetduration))))
    maxwetintensity<-lapply(diffwetval,function(x){abs(min(x))})
    wetseverity<-lapply(diffwetval,function(x){abs(sum(x))})
    
    #Total wet severity
    droughtlegacy_char[i,4]<-sum(unlist(wetseverity))
  }
}
colnames(droughtlegacy_char)<-paste(colnames(droughtlegacy_char),'01_19',sep = '_')
droughtlegacy_char<-cbind(co_nona,droughtlegacy_char)
write.csv(droughtlegacy_char,'droughtlegacy_01_19_nonaslope.csv')

#convert to rasterstack
droughtlegacy_char_df<-as.data.frame(matrix(nrow=372732,ncol=6))
droughtlegacy_char_df[nona_index,]<-droughtlegacy_char[,3:6]
colnames(droughtlegacy_char_df)<-colnames(droughtlegacy_char)[3:6]
droughtlegacy_char_rasterstack <- brick(extent(evi), nl=4, nrows=nrow(evi),
                                        ncols=ncol(evi), crs=proj4string(evi))
for(k in 1:4){
  droughtlegacy_char_rasterstack[[k]] <- droughtlegacy_char_df[,k]
}
names(droughtlegacy_char_rasterstack)<-colnames(droughtlegacy_char_df)
writeRaster(droughtlegacy_char_rasterstack,'droughtlegacy_01_19_nonaslope.tif',overwrite=T)

#Load input data model ####
seas<-raster('seasonality_index_80_19.tif')
interannvar<-raster('betweenyearvariation_80_19.tif')
tac_slope_nona<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]
seas_df<-raster::extract(seas, tac_slope_nona[,1:2])
interannvar_df<-raster::extract(interannvar, tac_slope_nona[,1:2])

#Create model trend 01-19 (60m) ~ drought 01-19 + climate####
#change tac_slope to _36m or _84m for impact different moving window in TAC calculation
tac_slope_nona<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]
#change droughtlegacy to _96_19, _91_19, _86_19 for impact longer drought period included
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
summary(model_01_19_lm)
VIF(model_01_19_lm)
## Variogram of the residuals
coords<- data.frame(lon=model_var_01_19$x, lat=model_var_01_19$y)
coordinates(coords)<- c("lon", "lat")
proj4string(coords)<- "+proj=longlat +ellps=clrk66 +no_defs"
coords_projected<- spTransform(coords, CRS("+proj=utm +zone=21 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
rs<-rstandard(model_01_19_lm)
spdata2 <- data.frame(resid = rs, x = coordinates(coords_projected)[,1], y = coordinates(coords_projected)[,2])
coordinates(spdata2)<-c("x","y")
proj4string(spdata2)<-"+proj=utm +zone=21 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
vario <- gstat::variogram(spdata2$resid ~ 1, locations=coordinates(spdata2), data=spdata2, width=5000, cutoff=80000)
exp.vario<-fit.variogram(vario, vgm("Exp"))
rang<-exp.vario$range[2] #range of 8 km?
nug<-exp.vario$psill[1] #nugget
## Replacing coordinates by projected coords
model_var_01_19$x <- coordinates(spdata2)[,1]
model_var_01_19$y <- coordinates(spdata2)[,2]
## Weights matrix
nb_dosi<-dnearneigh(as.matrix(model_var_01_19[,1:2]), 0, rang, longlat = FALSE)
listw_dosi<-nb2listw(nb_dosi, style="W", zero.policy = TRUE)
summary(unlist(listw_dosi$weights))
#Moran's I test
lm.morantest(model_01_19_lm,listw_dosi, zero.policy = T)
#do Lagrange multiplier diagnostics test
model_01_19_lm_test<-lm.LMtests(model_01_19_lm,listw_dosi, zero.policy = T,test="all")
summary(model_01_19_lm_test)
#Both spatial error and spatial lag are significant
W_dosi <- as(listw_dosi, "CsparseMatrix")
tr_dosi<-trW(W_dosi)
## Running  spatial lag model
model_01_19_lm_lag<-lagsarlm(f,data = model_var_01_19,
                             listw=listw_dosi, trs=tr_dosi,
                             method="Matrix", zero.policy = TRUE)
#running spatial error model
model_01_19_lm_err<-errorsarlm(f,data = model_var_01_19,
                               listw=listw_dosi, trs=tr_dosi,
                               method="Matrix", zero.policy = TRUE)
AIC(model_01_19_lm,model_01_19_lm_lag,model_01_19_lm_err)
summary(model_01_19_lm_lag)

#get info on average conditions in Amazon
droughtlegacy_clim_df_mean<-colMeans(droughtlegacy_clim_df[,3:8],na.rm=T)
droughtlegacy_clim_df_sd<-apply(droughtlegacy_clim_df[,3:8],2,function(x){sd(x,na.rm=T)})

#Create model for parts of Amazon: slope 01-19 (60m) ~ drought 01-19 + climate####
tac_slope_nona<-read.csv('evi_01_19_lag1_60m_nonaslope.csv')[,-1]
droughtlegacy_df<-read.csv('droughtlegacy_01_19_nonaslope.csv')[,-1]
droughtlegacy_df[is.na(droughtlegacy_df)]<-0
droughtlegacy_clim_df<-cbind(droughtlegacy_df,seas_df,interannvar_df)

#select pixels in NW Amazon: x < -66 and y > -3 (12309 pixels)
tac_slope_nona_nw<-tac_slope_nona[tac_slope_nona$x < (-66) & tac_slope_nona$y > (-3),]
droughtlegacy_clim_df_nw<-droughtlegacy_clim_df[tac_slope_nona$x < (-66) & tac_slope_nona$y > (-3),]
#select pixels in NE Amazon: x > -59 and y > -3 (11485 pixels)
tac_slope_nona_ne<-tac_slope_nona[tac_slope_nona$x > (-59) & tac_slope_nona$y > (-3),]
droughtlegacy_clim_df_ne<-droughtlegacy_clim_df[tac_slope_nona$x > (-59) & tac_slope_nona$y > (-3),]
#select pixels in SW Amazon: x < -70 and y < -5 (11031 pixels)
tac_slope_nona_sw<-tac_slope_nona[tac_slope_nona$x < (-70) & tac_slope_nona$y < (-3),]
droughtlegacy_clim_df_sw<-droughtlegacy_clim_df[tac_slope_nona$x < (-70) & tac_slope_nona$y < (-3),]
#select pixels in SE Amazon: x > -64 and y < -3 (12181 pixels)
tac_slope_nona_se<-tac_slope_nona[tac_slope_nona$x > (-64) & tac_slope_nona$y < (-3),]
droughtlegacy_clim_df_se<-droughtlegacy_clim_df[tac_slope_nona$x > (-64) & tac_slope_nona$y < (-3),]
#select pixels in S Amazon: 70 < x < -64 and y < -3 (12347 pixels)
tac_slope_nona_s<-tac_slope_nona[tac_slope_nona$x < (-64) & tac_slope_nona$x > (-70) & tac_slope_nona$y < (-3),]
droughtlegacy_clim_df_s<-droughtlegacy_clim_df[tac_slope_nona$x < (-64) & tac_slope_nona$x > (-70) & tac_slope_nona$y < (-3),]
#select pixels in N Amazon: -60 > x > -66 and y > -3 (12444 pixels)
tac_slope_nona_n<-tac_slope_nona[tac_slope_nona$x < (-59) & tac_slope_nona$x > (-66) & tac_slope_nona$y > (-3),]
droughtlegacy_clim_df_n<-droughtlegacy_clim_df[tac_slope_nona$x < (-59) & tac_slope_nona$x > (-66) & tac_slope_nona$y > (-3),]

#choose region
droughtlegacy_clim_df_st<-scale(droughtlegacy_clim_df_nw[,3:8]) 
model_var_01_19<-cbind(tac_slope_nona_nw,droughtlegacy_clim_df_st)
colnames(model_var_01_19)[c(3)]<-c('slope')

f <- as.formula(slope ~ d_amount_01_19 + d_avint_01_19 + d_avdur_01_19 
                + d_avint_01_19:d_amount_01_19 + d_avdur_01_19:d_amount_01_19 
                + w_totsev_01_19 + w_totsev_01_19:d_amount_01_19 
                + seas_df + interannvar_df
                + d_amount_01_19:seas_df + d_amount_01_19:interannvar_df)
model_01_19_lm<-lm(f,data = model_var_01_19)
summary(model_01_19_lm)
VIF(model_01_19_lm)
## Variogram of the residuals
coords<- data.frame(lon=model_var_01_19$x, lat=model_var_01_19$y)
coordinates(coords)<- c("lon", "lat")
proj4string(coords)<- "+proj=longlat +ellps=clrk66 +no_defs"
coords_projected<- spTransform(coords, CRS("+proj=utm +zone=21 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
rs<-rstandard(model_01_19_lm)
spdata2 <- data.frame(resid = rs, x = coordinates(coords_projected)[,1], y = coordinates(coords_projected)[,2])
coordinates(spdata2)<-c("x","y")
proj4string(spdata2)<-"+proj=utm +zone=21 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
vario <- gstat::variogram(spdata2$resid ~ 1, locations=coordinates(spdata2), data=spdata2, width=5000, cutoff=80000)
exp.vario<-fit.variogram(vario, vgm("Exp"))
rang<-exp.vario$range[2] #range of 8 km?
nug<-exp.vario$psill[1] #nugget
## Replacing coordinates by projected coords
model_var_01_19$x <- coordinates(spdata2)[,1]
model_var_01_19$y <- coordinates(spdata2)[,2]
## Weights matrix
nb_dosi<-dnearneigh(as.matrix(model_var_01_19[,1:2]), 0, rang, longlat = FALSE)
listw_dosi<-nb2listw(nb_dosi, style="W", zero.policy = TRUE)
summary(unlist(listw_dosi$weights))
#Moran's I test
lm.morantest(model_01_19_lm,listw_dosi, zero.policy = T)
#do Lagrange multiplier diagnostics test
model_01_19_lm_test<-lm.LMtests(model_01_19_lm,listw_dosi, zero.policy = T,test="all")
summary(model_01_19_lm_test)
#Both spatial error and spatial lag are significant
W_dosi <- as(listw_dosi, "CsparseMatrix")
tr_dosi<-trW(W_dosi)
## Running  spatial lag model
model_01_19_lm_lag<-lagsarlm(f,data = model_var_01_19,
                             listw=listw_dosi, trs=tr_dosi,
                             method="Matrix", zero.policy = TRUE)
summary(model_01_19_lm_lag)

#get info on average conditions in region of Amazon
colMeans(droughtlegacy_clim_df_nw[,3:8],na.rm=T)
apply(droughtlegacy_clim_df_nw[,3:8],2,function(x){sd(x,na.rm=T)})

