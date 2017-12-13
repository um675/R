library(sp)
library(raster)
library(Rgbp)
library(rgdal)
library(sn)
library(stats4)
library(mnormt)
library(GEOmap)
library(maptools)
library(mapproj)
library(geomapdata)
library(shapefiles)
library(tiff)


setClass("RasterPlus",
         representation(r1="RasterLayer",r2="RasterLayer",t1="numeric",t2="numeric")
)

setGeneric("maskApply",
           function(x)
            standardGeneric("maskApply")
)

setMethod("maskApply","RasterPlus",
          function(x){
            mask=x@r1
            victim=x@r2
            result=victim
            values(result)=values(victim)*values(mask)
            return(result)
          }
            
              
)

setGeneric("maskOUT",
           function(x)
             standardGeneric("maskOUT")
)

setMethod("maskOUT","RasterPlus",
          function(x){
            ORraster=x@r2
            tr1=x@t1
            tr2=x@t2
            indici=which((values(ORraster)<tr1) | (values(ORraster)>tr2))
            values(ORraster)[indici]=NaN
            result=ORraster
            return(result)
          }
          
          
)