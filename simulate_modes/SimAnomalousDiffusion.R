AnomalousDiffusion <- function(n,dc,dt,dr,scene,sceneWidth,sceneHeight,
                               x=sceneWidth/2,y=sceneHeight/2){
  # Simulates a anomalous diffusion of particle with diffusion coefficient 
  # dc and a timelag between two sucessive steps of dt.
  #
  # Args:
  #   n: Number of stebs.
  #   dc: Diffusion coefficient [µm²/s]
  #   dt: time lag between two steps [s]
  #   dr: Drift in [µm/s]. The direction will be choosen randomly.
  #   scene: binary matrix with where 1=object and 0 = no object
  #   scenewidth: Width of the scene in µm
  #   sceneheight: Height of the scene in µm
  #   x: start position x axis
  #   y: start position y axis
  #
  # Returns:
  #   Vector of positions
  
  # Split a single step for the time lag dt in 100 steps a dt/100
  posX <- vector(length=n);
  posY <- vector(length=n);
  posX[1] <- x
  posY[1] <- y
  driftDirection <- pi#runif(1,min=0,max=1)*2*pi;
  for(i in 2:(n+1)){
    p <- AnomalousDiffusionSingleStep(x=posX[i-1],y=posY[i-1],dc=dc,dt =dt,
                                      dr=dr,beta=driftDirection,scene = scene,
                                      sceneWidth = sceneWidth, sceneHeight = sceneHeight)
    posX[i] <- p[1]
    posY[i] <- p[2]
  }
  pos <- cbind(posX,posY);
  return(pos);
}

AnomalousDiffusionSingleStep <- function(x,y,dc,dt,dr,beta,scene,sceneWidth,
                                         sceneHeight){
  # Simulates a anomalous diffusion of a particle with diffusion coefficient 
  # dc and a timelag between two sucessive steps of dt. Therefore each step
  # is split up in 100 substeps. A substep which collidates with an object is
  # set to the previous position.
  #
  # Args:
  #   n: Number of stebs.
  #   dc: Diffusion coefficient [µm²/s]
  #   dt: time lag between two steps [s]
  #   dr: Drift in µm/s.
  #   beta: Direction of the drift [0,2pi]
  #   scene: binary matrix with where 1=object and 0 = no object
  #   scenewidth: Width of the scene in µm
  #   sceneheight: Height of the scene in µm
  #   x: start position x axis
  #   y: start position y axis
  #
  # Returns:
  #   Vector with a single (x,y) Position in µm
  
  #Check if the startposition is valid
  l = length(scene);
  wlRatio <- sceneWidth/l;
  hlRatio <- sceneHeight/l; 
  sceneIndexX = as.integer(x*wlRatio)+1 #as.integer((x/sceneWidth)*length(scene)) +1;
  sceneIndexY = as.integer(y*hlRatio)+1  #as.integer((y/sceneHeight)*length(scene)) +1 ;
  
  if(scene[sceneIndexX,sceneIndexY]>0){
    stop("Start position is not a free position")
  }
  
  numberOfSubsteps = 100;
  dt = dt / numberOfSubsteps
  
  #Drift
  
  driftStepLength <- dt*dr;
  driftX <- cos(beta)*driftStepLength
  driftY <- sin(beta)*driftStepLength
  
  stepCounter <- 0;
  repeat{
    #Random steplength
    u <- runif(1,min=0,max=1)
    steplength <- sqrt(-4*dc*dt*log(1-u))
    
    #Random direction
    alpha <- runif(1,min=0,max=1)*2*pi;
    
    #Random Positions
    stepX <- cos(alpha)*steplength+driftX
    stepY <- sin(alpha)*steplength+driftY
    newX <- x + stepX;
    newY <- y + stepY;
    sceneIndexX = as.integer(newX*wlRatio)+1 #as.integer(newX/sceneWidth * length(scene)) +1;
    sceneIndexY = as.integer(newY*hlRatio)+1 #as.integer(newY/sceneHeight * length(scene)) + 1;
    stepCounter=stepCounter+1;
    if(scene[sceneIndexX,sceneIndexY]==0){
      x=newX;
      y=newY;
    }
    
    if(stepCounter==numberOfSubsteps){
      break;
    }
  }
  
  return(c(x,y))
}