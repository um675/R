library(sp)
library(raster)
library(Rgbp)
library(rgdal)
library(sn)
library(stats4)
library(mnormt)
library(GEOmap)
library(PythonInR)
library(maptools)
library(mapproj)
library(geomapdata)
library(shapefiles)
library(tiff)
library(Rquake)
library(RSEIS)
library(cluster)
library(randomForest)


setClass("Classificazione",
         representation(immagine="RasterBrick",nclassi="numeric")
)
setGeneric("k_means",
           function(x)
             standardGeneric("k_means")
)
setMethod("k_means","Classificazione",
  function(x){
    lsat=x@immagine
    stk_lsat=stack(lsat)
    v_lsat=getValues(lsat)
    ind=which(!is.na(v_lsat))
    v_lsat=na.omit(v_lsat)
    n=x@nclassi
    ####################################
    ###### K means classification ######
    ####################################
    E=kmeans(v_lsat,n,iter.max = 100,nstart=n)
    kmeans_raster=raster(stk_lsat)
    kmeans_raster[ind]=E$cluster
    #plot(kmeans_raster)
    #sort(unique(values(kmeans_raster)))
    return(kmeans_raster)
  }
)


setGeneric("kmeansRandomForest",
           function(x)
             standardGeneric("kmeansRandomForest")
)

setMethod("kmeansRandomForest","Classificazione",
 function(x){
    lsat=x@immagine
    stk_lsat=stack(lsat)
    v_lsat=getValues(lsat)
    ind=which(!is.na(v_lsat))
    v_lsat=na.omit(v_lsat)
    n=x@nclassi
    ###########################
    ##### Random Forest #######
    ###########################
    
    ## unsupervised randomForest classification using kmeans
    vx<-v_lsat[sample(nrow(v_lsat), 500),]
    rf = randomForest(vx)
    rf_prox <- randomForest(vx,ntree = 1000, proximity = TRUE)$proximity
    
    E_rf <- kmeans(rf_prox, n, iter.max = 100, nstart = n)
    rf <- randomForest(vx,as.factor(E_rf$cluster),ntree = 500)
    rf_raster<- predict(stk_lsat,rf)
    return(rf_raster)
  }
)






