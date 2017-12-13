############ FUNZIONE CHE RESTITUISCE IL RASTER in Db (sia raster che brick) ###############
rasterDecibel=function(datoraster){           
  
  decibellato=10*(log(values(datoraster),base=10))
  rdecibellato=datoraster
  values(rdecibellato)=decibellato
  
  return(rdecibellato)
  
}
############################################################################################