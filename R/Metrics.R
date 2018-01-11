#'@title Class Area
#'
#'@description Get the area in hectares (ha) of raster classes.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value in hectares(ha) of the class.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@import raster
#'@import rgdal
#'@import rgeos
#' @import igraph
#'@export
CA <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment_table=data.frame(cell_fragment)
  cell_fragment_table=cell_fragment_table[!is.na(cell_fragment_table[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(cell_fragment_table$count)
  total_area_class_ha=sum_area * area_pixel/10000
  return(total_area_class_ha)
}


#'@title Number of Patches
#'
#'@description Get the number of patches of the corresponding patch type (class).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the number of patches of the respective class.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#' @export
NP <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment_table=data.frame(cell_fragment)
  cell_fragment_table=cell_fragment_table[!is.na(cell_fragment_table[,1]),]
  NP=max(cell_fragment_table$value)
  return(NP)
}


#'@title Mean Patch Size
#'
#'@description Get the mean patch size in hectares (ha).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean patch size in hectares (ha).
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#' @export
MPS <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment=cell_fragment[!is.na(cell_fragment[,1]),]
  area_pixel=res(x)[1] * res(x)[2]
  table_fragm=data.frame(cell_fragment)
  table_fragm$size_ha=(table_fragm$count*area_pixel/10000)
  MPS=sum(table_fragm$size_ha)/length(cell_fragment[,"count"])
  return(MPS)
}


#'@title Median Patch Size
#'
#'@description Get the fragment size found in the median in hectares (ha).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the patch size in the median in hectares (ha).
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#' @export
MedPS <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment=cell_fragment[!is.na(cell_fragment[,1]),]
  area_pixel=res(x)[1]*res(x)[2]
  table_fragm=data.frame(cell_fragment)
  table_fragm$size_ha=(table_fragm$count*area_pixel/10000)
  MedPS=median(table_fragm$size_ha)
  return(MedPS)
}


#'@title Patch Size Standard Deviation
#'
#'@description Calculates the standard deviation of the patch size in hectares (ha).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the patch size standard deviation in hectares (ha).
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#' @export
PSSD <- function(x){
  class_forest=x== 1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment=cell_fragment[!is.na(cell_fragment[,1]),]
  area_pixel=res(x)[1] * res(x)[2]
  a=cell_fragment[ ,"count"]
  b=sqrt(sum((a - mean(a))^2)/(length(a)))
  PSSD=b*area_pixel/10000
  return(PSSD)
}


#'@title Patch Size Coefficient of Variation
#'
#'@description Calculates of the patch size coefficient of variation in percentage.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the patch size coefficient of variation in percentage.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
PSCoV <- function(x) {
  class_forest=x== 1
  fragment_forest=clump(class_forest)
  cell_fragment=freq(fragment_forest)
  cell_fragment=cell_fragment[!is.na(cell_fragment[,1]),]
  area_pixel=res(x)[1]*res(x)[2]
  a=cell_fragment[ ,"count"]
  b=sqrt(sum((a - mean(a))^2)/(length(a)))
  PSSD=b * area_pixel/10000
  table_fragm=data.frame(cell_fragment)
  table_fragm$size_ha=(table_fragm$count*area_pixel/10000)
  MPS=sum(table_fragm$size_ha)/length(cell_fragment[,"count"])
  PSCoV=(PSSD/MPS)*100
  return(PSCoV)
}


#'@title Total Edge
#'
#'@description Calculates the lengths, in meters (m), of all edge segm ents involving the corresponding patch type.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The sum of the lengths, in meters (m), of all edge segments involving the corresponding patch type.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TE <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  TE=gLength(polygons)
  return(TE)
}

#'@title Edge Density
#'
#'@description Calculates the total weighted class edges by the landscape area in meters per hectare (m/ha).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@param a Total area of the landscape.
#'@return The edge density in meteres per hectare (m/ha).
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export

ED <- function(x,a){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  TE=gLength(polygons)
  TE_res=TE/a
  return(TE_res)
}


#'@title Landscape Shape Index
#'
#'@description Calculates the landscape shape index.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the landscape shape index.
#'@details The value must be greater than or equal to 1.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
LSI <- function(x) {
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  area_pixel=res(x)[1] * res(x)[2]
  a= gArea(polygons)/area_pixel
  TE=gLength(polygons)/res(x)[1]
  b=trunc(sqrt(a))
  c=a - b^2
  minTE=rep(0,length(c))
  for (ii in 1:length(c)){
    if (c[ii]==0) minTE[ii]=4*b[ii]
    if (b[ii]^2<a[ii] & a[ii]<=b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 2
    if (a[ii] > b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 4
  }
  return(TE/minTE)
}


#'@title Mean Shape Index
#'
#'@description Calculates the mean shape index.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean shape index.
#'@details The value must be greater than or equal to 1.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MSI <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  areas_m2=gArea(polygons, byid = T)
  area_pixel=res(x)[1] * res(x)[2]
  areas_pixels=areas_m2/area_pixel
  perims_m=gLength(polygons, byid = T)
  perims_edges=perims_m/res(x)[1]
  d=data.frame(id = names(areas_m2), perims_edges, perims_m, areas_pixels, areas_m2)
  n=length(polygons)
  {
    LSI = function(a,TE) {
      b=trunc(sqrt(a))
      c=a - b^2
      minTE=rep(0,length(c))
      for (ii in 1:length(c)){
        if (c[ii]==0) minTE[ii]=4*b[ii]
        if (b[ii]^2<a[ii] & a[ii]<=b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 2
        if (a[ii] > b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 4
      }
      return(TE/minTE)
    }
  }
  d$SHAPE=LSI(d$areas_pixels, d$perims_edges)
  MSI=(sum(d$SHAPE))/(n)
  return(MSI)
}


#'@title Area-Weighted Mean Shape Index
#'
#'@description Calculates the area-weighted mean shape index.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the area-weighted mean shape index.
#'@details The value must be greater than or equal to 1.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
AWMSI <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  areas_m2=gArea(polygons, byid = T)
  area_pixel=res(x)[1] * res(x)[2]
  areas_pixels=areas_m2/area_pixel
  perims_m=gLength(polygons, byid = T)
  perims_edges=perims_m/res(x)[1]
  d=data.frame(id = names(areas_m2), perims_edges, perims_m, areas_pixels, areas_m2)
  n=length(polygons)
  {
    LSI = function(a,TE) {
      b=trunc(sqrt(a))
      c=a - b^2
      minTE=rep(0,length(c))
      for (ii in 1:length(c)){
        if (c[ii]==0) minTE[ii]=4*b[ii]
        if (b[ii]^2<a[ii] & a[ii]<=b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 2
        if (a[ii] > b[ii]*(1+b[ii])) minTE[ii]=4 * b[ii] + 4
      }
      return(TE/minTE)
    }
  }
  d$SHAPE=LSI(d$areas_pixels, d$perims_edges)
  d$area_weighted=d$areas_m2/sum(areas_m2)
  d$AWMSI=(d$SHAPE)*(d$area_weighted)
  AWMSI=sum(d$AWMSI)
  return(AWMSI)
}


#'@title Mean Patch Fractal Dimension
#'
#'@description Calculates the mean patch fractal dimension.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean patch fractal dimension.
#'@details The MPFD value varies between 1 and 2.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MPFD <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  areas_m2=gArea(polygons, byid = T)
  perims_m=gLength(polygons, byid = T)
  d=data.frame(id = names(areas_m2), perims_m, areas_m2)
  n=length(polygons)
  MPFD=(2*log(0.25 * d$perims_m))/log(d$areas_m2)
  MPFD=sum(MPFD)/n
  return(MPFD)
}


#'@title Area-Weighted Mean Patch Fractal Dimension
#'
#'@description Calculates the area-weighted mean patch fractal dimension.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the area-weighted mean patch fractal dimension.
#'@details The AWMPFD value varies between 1 and 2.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
AWMPFD <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  areas_m2=gArea(polygons, byid = T)
  perims_m=gLength(polygons, byid = T)
  d=data.frame(id = names(areas_m2), perims_m, areas_m2)
  n=length(polygons)
  d$area_weighted=d$areas_m2/sum(areas_m2)
  AWMPFD=(2*log(0.25 * d$perims_m))/log(d$areas_m2)
  AWMPFD=AWMPFD * d$area_weighted
  AWMPFD=sum(AWMPFD)
  return(AWMPFD)
}


#'@title Mean Perimeter-Area Ratio
#'
#'@description Calculates de mean perimeter area ratio in hectares (ha).
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean perimeter area ratio in hectares (ha).
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MPAR <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  areas_m2=gArea(polygons, byid = T)
  perims_m=gLength(polygons, byid = T)
  d=data.frame(id = names(areas_m2), perims_m, areas_m2)
  n=length(polygons)
  a=data.frame(id = names(areas_m2), areas_m2/10000, perims_m)
  MPAR <- sum(a$perims_m/a$areas_m2.10000)
  MPAR <- MPAR/n
  return(MPAR)
}


#'@title Total Core Area with edge of 20 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export

TCA_20 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_20=(sum_area*area_pixel/10000)
  return(TCA_20)
}


#'@title Number of Core Area with edge of 20 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_20 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_20=length(polygons_2)
  return(NCA_20)
}


#'@title Mean Core Area with edge of 20 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_20 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_20=(sum_area*area_pixel/10000)
  NCA_20=length(polygons_2)
  MCA_20=TCA_20/NCA_20
  return(MCA_20)
}


#'@title Total Core Area Index with edge of 20 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_20 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_20=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_20=TCA_20/area_class*100
  return(TCAI_20)
}


#'@title Core Area Standard Deviation with edge of 20 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_20 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_20=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_20=CASD_20 * area_pixel/10000
  return(CASD_20)
}


#'@title Core Area Coefficient of Variation with edge of 20 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 20 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 20 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_20 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_20=gBuffer(polygons, byid=TRUE, width = -20.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_20)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_20=(sum(a))*area_pixel/10000
  NCA_20=length(frequency_fragm[,"count"])
  MCA_20=TCA_20/NCA_20
  CASD_20=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_20=CASD_20 * area_pixel/10000
  CACV_20=(CASD_20/MCA_20)*100
  return(CACV_20)
}


#'@title Total Core Area with edge of 40 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_40 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_40=(sum_area*area_pixel/10000)
  return(TCA_40)
}

#'@title Number of Core Area with edge of 40 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_40 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_40=length(polygons_2)
  return(NCA_40)
}


#'@title Mean Core Area with edge of 40 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#' @export
MCA_40 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_40=(sum_area*area_pixel/10000)
  NCA_40=length(polygons_2)
  MCA_40=TCA_40/NCA_40
  return(MCA_40)
}


#'@title Total Core Area Index with edge of 40 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_40 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_40=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_40=TCA_40/area_class*100
  return(TCAI_40)
}


#'@title Core Area Standard Deviation with edge of 40 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_40 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_40=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_40=CASD_40 * area_pixel/10000
  return(CASD_40)
}


#'@title Core Area Coefficient of Variation with edge of 40 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 40 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 40 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_40 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_40=gBuffer(polygons, byid=TRUE, width = -40.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_40)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_40=(sum(a))*area_pixel/10000
  NCA_40=length(frequency_fragm[,"count"])
  MCA_40=TCA_40/NCA_40
  CASD_40=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_40=CASD_40 * area_pixel/10000
  CACV_40=(CASD_40/MCA_40)*100
  return(CACV_40)
}


#'@title Total Core Area with edge of 60 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_60 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_60=(sum_area*area_pixel/10000)
  return(TCA_60)
}

#'@title Number of Core Area with edge of 60 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_60 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_60=length(polygons_2)
  return(NCA_60)
}


#'@title Mean Core Area with edge of 60 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_60 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_60=(sum_area*area_pixel/10000)
  NCA_60=length(polygons_2)
  MCA_60=TCA_60/NCA_60
  return(MCA_60)
}


#'@title Total Core Area Index with edge of 60 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_60 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_60=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_60=TCA_60/area_class*100
  return(TCAI_60)
}


#'@title Core Area Standard Deviation with edge of 60 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_60 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_60=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_60=CASD_60 * area_pixel/10000
  return(CASD_60)
}


#'@title Core Area Coefficient of Variation with edge of 60 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 60 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 60 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_60 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_60=gBuffer(polygons, byid=TRUE, width = -60.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_60)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_60=(sum(a))*area_pixel/10000
  NCA_60=length(frequency_fragm[,"count"])
  MCA_60=TCA_60/NCA_60
  CASD_60=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_60=CASD_60 * area_pixel/10000
  CACV_60=(CASD_60/MCA_60)*100
  return(CACV_60)
}


#'@title Total Core Area with edge of 80 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_80 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_80=(sum_area*area_pixel/10000)
  return(TCA_80)
}


#'@title Number of Core Area with edge of 80 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_80 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_80=length(polygons_2)
  return(NCA_80)
}


#'@title Mean Core Area with edge of 80 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_80 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_80=(sum_area*area_pixel/10000)
  NCA_80=length(polygons_2)
  MCA_80=TCA_80/NCA_80
  return(MCA_80)
}


#'@title Total Core Area Index with edge of 80 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_80 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_80=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_80=TCA_80/area_class*100
  return(TCAI_80)
}


#'@title Core Area Standard Deviation with edge of 80 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_80 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_80=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_80=CASD_80 * area_pixel/10000
  return(CASD_80)
}


#'@title Core Area Coefficient of Variation with edge of 80 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 80 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 80 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_80 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_80=gBuffer(polygons, byid=TRUE, width = -80.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_80)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_80=(sum(a))*area_pixel/10000
  NCA_80=length(frequency_fragm[,"count"])
  MCA_80=TCA_80/NCA_80
  CASD_80=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_80=CASD_80 * area_pixel/10000
  CACV_80=(CASD_80/MCA_80)*100
  return(CACV_80)
}


#'@title Total Core Area with edge of 100 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_100 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_100=(sum_area*area_pixel/10000)
  return(TCA_100)
}


#'@title Number of Core Area with edge of 100 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_100 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_100=length(polygons_2)
  return(NCA_100)
}


#'@title Mean Core Area with edge of 100 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_100 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_100=(sum_area*area_pixel/10000)
  NCA_100=length(polygons_2)
  MCA_100=TCA_100/NCA_100
  return(MCA_100)
}


#'@title Total Core Area Index with edge of 100 meters
#'
#'@description Measurement of relative quantity of core areain the landscape, in percentage, when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_100 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_100=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_100=TCA_100/area_class*100
  return(TCAI_100)
}


#'@title Core Area Standard Deviation with edge of 100 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_100 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_100=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_100=CASD_100 * area_pixel/10000
  return(CASD_100)
}


#'@title Core Area Coefficient of Variation with edge of 100 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 100 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 100 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_100 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_100=gBuffer(polygons, byid=TRUE, width = -100.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_100)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_100=(sum(a))*area_pixel/10000
  NCA_100=length(frequency_fragm[,"count"])
  MCA_100=TCA_100/NCA_100
  CASD_100=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_100=CASD_100 * area_pixel/10000
  CACV_100=(CASD_100/MCA_100)*100
  return(CACV_100)
}


#'@title Total Core Area with edge of 140 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_140 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_140=(sum_area*area_pixel/10000)
  return(TCA_140)
}


#'@title Number of Core Area with edge of 140 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_140 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_140=length(polygons_2)
  return(NCA_140)
}


#'@title Mean Core Area with edge of 140 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_140 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_140=(sum_area*area_pixel/10000)
  NCA_140=length(polygons_2)
  MCA_140=TCA_140/NCA_140
  return(MCA_140)
}


#'@title Total Core Area Index with edge of 140 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_140 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_140=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_140=TCA_140/area_class*100
  return(TCAI_140)
}


#'@title Core Area Standard Deviation with edge of 140 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_140 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_140=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_140=CASD_140 * area_pixel/10000
  return(CASD_140)
}


#'@title Core Area Coefficient of Variation with edge of 140 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 140 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 140 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_140 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_140=gBuffer(polygons, byid=TRUE, width = -140.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_140)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_140=(sum(a))*area_pixel/10000
  NCA_140=length(frequency_fragm[,"count"])
  MCA_140=TCA_140/NCA_140
  CASD_140=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_140=CASD_140 * area_pixel/10000
  CACV_140=(CASD_140/MCA_140)*100
  return(CACV_140)
}


#'@title Total Core Area with edge of 200 meters
#'
#'@description Calculates the total core area, in hectare (ha), when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core are in hectare (ha), with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCA_200 <- function (x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_200=(sum_area*area_pixel/10000)
  return(TCA_200)
}

#'@title Number of Core Area with edge of 200 meters
#'
#'@description Calculates the total number of core area when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total number of core area, with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
NCA_200 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  NCA_200=length(polygons_2)
  return(NCA_200)
}


#'@title Mean Core Area with edge of 200 meters
#'
#'@description Calculates the mean core area, in hectare (ha), when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean core are in hectare (ha), with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MCA_200 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  polygons_2=rasterToPolygons(Fragments, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(x)[1] * res(x)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_200=(sum_area*area_pixel/10000)
  NCA_200=length(polygons_2)
  MCA_200=TCA_200/NCA_200
  return(MCA_200)
}


#'@title Total Core Area Index with edge of 200 meters
#'
#'@description Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the total core area index in percentage, with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
TCAI_200 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm_table=data.frame(frequency_fragm)
  frequency_fragm_table[is.na(frequency_fragm_table)]=0
  frequency_fragm_table[frequency_fragm_table$value==0,]=0
  area_pixel=res(raster)[1] * res(raster)[2]
  sum_area=sum(frequency_fragm_table$count)
  TCA_200=(sum_area*area_pixel/10000)
  area_class=gArea(polygons)
  area_class=area_class/10000
  TCAI_200=TCA_200/area_class*100
  return(TCAI_200)
}


#'@title Core Area Standard Deviation with edge of 200 meters
#'
#'@description Calculates the patch core area standard deviation in hectare (ha), when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area standard deviation in hectare (ha), with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CASD_200 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  CASD_200=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_200=CASD_200 * area_pixel/10000
  return(CASD_200)
}


#'@title Core Area Coefficient of Variation with edge of 200 meters
#'
#'@description Calculates the patch core area coefficient of variation in pecentage, when excluded the edge from 200 meters.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the core area coefficient of variation in percentage, with subtracted edge of 200 meters.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
CACV_200 <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  Buffer_200=gBuffer(polygons, byid=TRUE, width = -200.0, capStyle = "SQUARE", joinStyle = "MITRE", mitreLimit = 5.0)
  {
    clip<-function(raster,shape) {
      a1_crop<-crop(raster,shape)
      step1<-rasterize(shape,a1_crop)
      a1_crop*step1}
  }
  clip=clip(fragment_forest,Buffer_200)
  Fragments=clump(clip)
  frequency_fragm=freq(Fragments)
  frequency_fragm=frequency_fragm[!is.na(frequency_fragm[,1]),]
  area_pixel=res(raster)[1] * res(raster)[2]
  a=frequency_fragm[ ,"count"]
  TCA_200=(sum(a))*area_pixel/10000
  NCA_200=length(frequency_fragm[,"count"])
  MCA_200=TCA_200/NCA_200
  CASD_200=sqrt(sum((frequency_fragm[,"count"] - mean(frequency_fragm[,"count"]))^2) / (length(frequency_fragm[,"count"])))
  CASD_200=CASD_200 * area_pixel/10000
  CACV_200=(CASD_200/MCA_200)*100
  return(CACV_200)
}


#'@title Mean Nearest-Neighbor Distance
#'
#'@description Calculates the mean nearest neighbor distance in meters (m) between patchs.
#'@param x Raster of data with patches identified as classes in RasterLayer format.
#'@return The value of the mean nearest neighbor distance in meters (m) between patchs.
#'@references MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
#'@export
MNN <- function(x){
  class_forest=x==1
  fragment_forest=clump(class_forest)
  polygons=rasterToPolygons(fragment_forest, fun=NULL, na.rm=TRUE, dissolve=TRUE)
  dist_pol=gDistance(polygons, byid = TRUE)
  dist_min_value=apply(dist_pol, 2, FUN = function(x) {min(x[x > 0])})
  table_data=cbind(polygons@data, dist_min_value)
  table_data$NEAR=(table_data$dist_min_value)+(min(dist_pol[dist_pol > 0]))
  MNN=(sum(table_data$NEAR))/(length(polygons))
  return(MNN)
}


