ConfinedDiffusion <- function(n,dc,dt,dr,r,x=0,y=0){
  # Simulates a confined diffusion of particle inside a circle
  #
  # Args:
  #   n: Number of stebs.
  #   dc: Diffusion coefficient [µm²/s]
  #   dt: time lag between two steps [s]
  #   dr: Drift of the circle in [µm/s]. The direction will be choosen randomly.
  #   r: radius of the circular region.
  #   x: x position of the circle
  #   y: y position of the circle
  #
  # Returns:
  #   Vector of positions
  
  posX <- vector(length=n);
  posY <- vector(length=n);
  posX[1] <- 0
  posY[1] <- 0
  
  for(i in 2:n){
    p <- ConfindedDiffusionSingleStep(x=posX[i-1],y=posY[i-1],dc=dc,dt=dt,r = r)
    posX[i] <- p[1]
    posY[i] <- p[2]
  }
  #Add and shift to  start position
  posX = posX + x
  posY = posY + y
  
  #Add drift to circle
  driftDirection <- runif(1,min=0,max=1)*2*pi;
  driftStepLength <- dt*dr;
  driftX <- cos(driftDirection)*driftStepLength
  driftY <- sin(driftDirection)*driftStepLength
  posX <- posX + cumsum(rep(driftX,n))
  posY <- posY + cumsum(rep(driftY,n))
  posX = c(x,posX)
  posY = c(y,posY)
  pos <- cbind(posX,posY);
  return(pos);
}

ConfindedDiffusionSingleStep <- function(dc,dt,r,x,y){
  # Simulates a single step (for dt) of a confined diffusion inside of a circle. 
  # Therefore each step is split up in 100 substeps. A substep which collidates 
  # with an object is set to the previous position.
  #
  # Args:
  #   n: Number of stebs.
  #   dc: Diffusion coefficient [µm²/s]
  #   dt: time lag between two steps [s]
  #   r: radius of the circular region.
  #   x: x position of the circle
  #   y: y position of the circle
  #
  # Returns:
  #   Vector of position
  
  numberOfSubsteps = 100;
  dt = dt / numberOfSubsteps
  
  stepCounter <- 0;
  repeat{
    #Random steplength
    u <- runif(1,min=0,max=1)
    steplength <- sqrt(-4*dc*dt*log(1-u))
    
    #Random direction
    alpha <- runif(1,min=0,max=1)*2*pi;
    
    #Random Positions
    stepX <- cos(alpha)*steplength
    stepY <- sin(alpha)*steplength
    newX <- x + stepX;
    newY <- y + stepY;
    
    
    #Bezugspunkt ist immer (0,0)
    stepCounter <- stepCounter+1;
    distanceFromZero = sqrt((newX)^2 + (newY)^2)
    if(distanceFromZero<=r){
      x<-newX;
      y<-newY;
      
    }
    if(stepCounter==numberOfSubsteps){
      break;
    }
  }
  
  return(c(x,y))
}
