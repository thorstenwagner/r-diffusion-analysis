FreeDiffusion <- function(n,dc,dt,dr, x=0,y=0){
  # Simulates a free diffusion of particle with diffusion coefficient 
  # dc and a timelag between two sucessive steps of dt.
  #
  # Args:
  #   n: Number of steps.
  #   dc: Diffusion coefficient [µm²/s]
  #   dt: time lag between two steps [s]
  #   dr: Drift in [µm/s]. The drift direction will be selected randomly.
  #   x: start position x axis
  #   y: start position y axis
  #
  # Returns:
  #   Vector of positions
  
  #Random steplengths
  u <- runif(n,min=0,max=1)
  steplengths <- sqrt(-4*dc*dt*log(1-u))
  
  #Random directions
  alpha <- runif(n,min=0,max=1)*2*pi;
  
  #Drift
  beta <- runif(1,min=0,max=1)*2*pi;
  driftStepLength <- dt*dr;
  
  driftX <- cos(beta)*driftStepLength
  driftY <- sin(beta)*driftStepLength
  
  #Random Positions
  stepsX = x+cumsum(cos(alpha)*steplengths + driftX)
  stepsY = y+cumsum(sin(alpha)*steplengths + driftY)
  stepsX = c(x,stepsX)
  stepsY = c(y,stepsY)
  positions <- cbind(stepsX,stepsY)
  
  return(positions)
}