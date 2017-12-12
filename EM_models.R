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

#Crea la Classe
setClass("EM_models",
         representation(teta="RasterLayer",epsrel="RasterLayer",smc="RasterLayer",k="numeric",rms="numeric",lmbd="numeric",lcorr="numeric")
)

setGeneric("Dubois",
           function(object)
             standardGeneric("Dubois")
)

setMethod("Dubois","EM_models",
          function(object){
            map_epsrel=object@epsrel
            map_teta=object@teta
            k_nonda=object@k
            Hrms=object@rms
            lambda=object@lmbd
            

            
            #definisco gli elementi delle equazioni di dubois
            costeta=map_teta
            values(costeta)=cos(values(map_teta))
            sinteta=map_teta
            values(sinteta)=sin(values(map_teta))
            tanteta=map_teta
            values(tanteta)=values(sinteta)/values(costeta)
            krms=k_nonda*Hrms
            
            #definisco sottoparti delle equazioni
            
            #parte 1
            
            ratio_HH=map_teta
            values(ratio_HH)=((values(costeta))^(1.5))/((values(sinteta))^(5))
            parte1_HH=map_teta
            values(parte1_HH)=(10^(-2.75))*values(ratio_HH)
            
            ratio_VV=map_teta
            values(ratio_VV)=((values(costeta))^(3))/((values(sinteta))^(3))
            parte1_VV=map_teta
            values(parte1_VV)=(10^(-2.35))*values(ratio_VV)
            
            
            #parte 2
            
            esponente_HH=map_teta
            values(esponente_HH)=0.028*(values(map_epsrel))*(values(tanteta))
            parte2_HH=map_teta
            values(parte2_HH)=10^(values(esponente_HH))
            
            esponente_VV=map_teta
            values(esponente_VV)=0.046*(values(map_epsrel))*(values(tanteta))            
            parte2_VV=map_teta
            values(parte2_VV)=10^(values(esponente_VV))           
            
            #parte3
            
            parte3_HH=map_teta
            values(parte3_HH)=((krms*(values(sinteta)))^(1.4))*(lambda^(0.7))
            
            parte3_VV=map_teta
            values(parte3_VV)=((krms*(values(sinteta)))^(1.1))*(lambda^(0.7)) 
            
            #risultato
            sigmaHH=map_teta
            values(sigmaHH)=values(parte1_HH)*values(parte2_HH)*values(parte3_HH)
            sigmaVV=map_teta
            values(sigmaVV)=values(parte1_VV)*values(parte2_VV)*values(parte3_VV)
            
            sigma=brick(sigmaHH,sigmaVV)
            
            return(sigma)
          }
)


setGeneric("Oh",
           function(object)
             standardGeneric("Oh")
)


setMethod("Oh","EM_models",
          function(object){
            map_epsrel=object@epsrel
            map_teta=object@teta
            k_nonda=object@k
            Hrms=object@rms
            lambda=object@lmbd
            sm=object@smc
            
            
            
            #definisco gli elementi delle equazioni di dubois
            costeta=map_teta
            values(costeta)=cos(values(map_teta))
            sinteta=map_teta
            values(sinteta)=sin(values(map_teta))
            tanteta=map_teta
            values(tanteta)=values(sinteta)/values(costeta)
            krms=k_nonda*Hrms
            
            #definisco le 3 equazioni
            
            #equazione p (HH/VV)
            
            pimezzi=pi/2
            parteP1=values(map_teta)/pimezzi
            esponenteP1=(values(sm)^(-0.65))*0.35
            parteP2=exp((krms^(1.4))*(-0.4))
            num_p=1-((parteP1^esponenteP1)*parteP2)
            p=sm
            values(p)=num_p
            
            #equazione q (HV/VV)
            
            esponente_Q=exp((-1.3)*(krms^0.9))
            parteQ1=(0.13+sin(1.5*(values(map_teta))))^1.4
            num_Q=0.095*parteQ1*(1-esponente_Q)
            q=sm
            values(q)=num_Q

            #equazione HV

            parte1_HV=0.11*((values(sm))^(0.7))
            esponente_HV=exp((-0.32)*(krms^1.8))
            parte2_HV=((values(costeta))^(2.2))*(1-esponente_HV)
            num_HV=parte1_HV*parte2_HV
            HV=sm
            values(HV)=num_HV
            
            #risultato

            VV=sm
            values(VV)=values(HV)/values(q)
            HH=sm
            values(HH)=(values(VV))*(values(p))
            
            sigma_HH_VV_HV=brick(HH,VV,HV)
            
            return(sigma_HH_VV_HV)
          }
)

setGeneric("IEM",
           function(object)
             standardGeneric("IEM")
)


setMethod("IEM","EM_models",
          function(object){
            
            epsrel=object@epsrel
            map_teta=object@teta
            #k_nonda=object@k
            Hrms=object@rms
            lambda=object@lmbd
            sm=object@smc
            L_corr=object@lcorr
            
            y=Hrms
            l=L_corr #lunghezza correlazione
            x=epsrel #cost dielettrica
            #t1=30 # teta 
            #t2=t1*pi/180 # teta in radianti
            t2=map_teta
            k_s=(2*pi/lambda)
            
            #IEM backscattering coefficient
            
            
            
            
            
            
            ######## COMPONENTE HH #########
            
            A=((k_s)^2/2)*exp(-2*(k_s*cos(t2)*y)^2)
            Rpr=((cos(t2)-(x-sin(t2)^2)^0.5)/(cos(t2)+(x-sin(t2)^2)^0.5))
            fhh=(-2/cos(t2))*Rpr
            Shh=-1*(((2*sin(t2)^2)*(1+Rpr)^2)/cos(t2))*((x-sin(t2)^2-cos(t2)^2)/((cos(t2))^2))
            
            Ih=brick(Shh,Shh,Shh,Shh,Shh,Shh,Shh,Shh,Shh,Shh)
            values(Ih)=0
            W=Ih
            iterazioni=1:10
            
            for (i in iterazioni){
              
              Ih_b=(2*k_s*cos(t2)*y)^i*fhh*exp(-1*(k_s*cos(t2)*y)^2)+0.5*(k_s*cos(t2)*y)^i*Shh
              W_b=((l/i)^2*(1+(-2*k_s*sin(t2)*l/i)^2)^(-1.5))/(factorial(i))
              
              Ih[[i]]=Ih_b
              W[[i]]=W_b
              
              rm(Ih_b,W_b)  
            }
            
            
            ######## COMPONENTE VV #########
            
            Rpl=((x*cos(t2)-(x-sin(t2)^2)^0.5)/(x*cos(t2)+(x-sin(t2)^2)^0.5))
            fvv=(2/cos(t2))*Rpl
            Svv=((2*sin(t2)^2*(1+Rpl)^2)/cos(t2))*((1-(1/x))+(x-sin(t2)^2-x*cos(t2)^2)/((x*cos(t2))^2))
            
            Iv=Ih
            
            for (i in iterazioni){
              
              Iv_b=(2*k_s*cos(t2)*y)^i*fvv*exp(-1*(k_s*cos(t2)*y)^2)+0.5*(k_s*cos(t2)*y)^i*Svv
              
              
              Iv[[i]]=Iv_b
              
              
              rm(Iv_b)  
            }
            
            
            
            sum_iemh=Ih[[1]]^2*W[[1]]
            sum_iemv=Iv[[1]]^2*W[[1]]
            
            iter=2:10
            for (i in iter){
              
              sum_iemh=sum_iemh+Ih[[i]]^2*W[[i]]
              sum_iemv=sum_iemv+Ih[[i]]^2*W[[i]]
              
            }
            
            iemh=(A*sum_iemh)
            iemv=(A*sum_iemv)
            
            iem_HH_VV=brick(iemh,iemv)
            
            
            return(iem_HH_VV)
          }
)





