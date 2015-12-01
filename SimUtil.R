PlotMSD <- function(positions,lags,doplot){
  # Calculates the MSD curve
  #
  # Args:
  #   positions: x/y positions of the particle [µm]
  #   lags: vector of timelages [unitless, integer]

  
  dx= mapply(diff,lag=lags,MoreArgs=list(x=positions[,1]))
  dy= mapply(diff,lag=lags,MoreArgs=list(x=positions[,2]))
  
  msd <- vector(length=length(dx));
  weights <- vector(length=length(dx));
  for(i in 1:length(dx)) {
    msd[i] <- mean((dx[[i]]^2 + dy[[i]]^2));
    weights[i] <- length(dx[[i]]);
  }
  if(doplot){
    plot(lags,msd,type = "b");
  }
  
  return(list(MSD=msd,WEIGHTS=weights,LAG=lags))
}



############## NOT TESTED YET #####################
HydrDiameterToDiffucionCoefficient = function(dia, temp, visk) {
  # Converts to hydrodynamic diameter to the diffusion coeffcient
  #
  # Args:
  #   dia: hydrodynamic diameter [µm]
  #   temp: temperature in kelvin [k]
  #   visk: viscosity 
  #
  # Returns:
  #   Diameter in [µm] oder DC 
  kBoltz = 1.3806488*10^-23 #Boltzmann-Konstante 
  2*kBoltz*temp/(dia*6*pi*visk)*10^24
}

