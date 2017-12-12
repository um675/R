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
rm(list=ls())

##################
# Crea la Classe #
##################

setClass("RLmath",
         representation(r="RasterLayer")
)

setGeneric("SMC2dielconst",
           function(x)
             standardGeneric("SMC2dielconst")
)

setMethod("SMC2dielconst","RLmath",
          function(object){
            mappaSMC=object@r
            dielconst=mappaSMC
            y=values(mappaSMC)
            values(dielconst)=3.03+(9.3*y)+(146*(y^2))+(76.7*(y^3))
            return(dielconst)
          }
)

setGeneric("SMC2diel_orgsoil",
           function(x)
             standardGeneric("SMC2diel_orgsoil")
)

setMethod("SMC2diel_orgsoil","RLmath",
          function(object){
            mappaSMC=object@r
            dielconst=mappaSMC
            y=values(mappaSMC)
            values(dielconst)=1.74-(0.34*y)+(135*(y^2))-(55*(y^3))
            return(dielconst)
          }
)
#Il seguente è il migliore
setGeneric("SMC2diel_ulaby",
           function(x)
             standardGeneric("SMC2diel_ulaby")
)

setMethod("SMC2diel_ulaby","RLmath",
          function(object){
            mappaSMC=object@r
            Clay=0.74
            Silt=0.26
            
            #### Per Banda L valgono questi coefficienti ####
            a0=2.862
            a1=-0.012
            a2=0.001
            
            b0=3.803
            b1=0.462
            b2=-0.341
            
            c0=119.006
            c1=-0.500
            c2=0.633
            ##################################################
            A=a0+(a1*Silt)+(a2*Clay)
            B=b0+(b1*Silt)+(b2*Clay)
            C=c0+(c1*Silt)+(c2*Clay)
            
            dielconst=mappaSMC
            y=values(mappaSMC)
            values(dielconst)=A+(B*y)+(C*(y^2))
            return(dielconst)
          }
)






################################
#Costante dielettrica nel vuoto#
#    eps=8.85418781762E-12     #
################################
