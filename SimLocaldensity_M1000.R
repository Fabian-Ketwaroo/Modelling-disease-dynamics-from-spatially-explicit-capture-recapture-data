#=====================================================================================================================================================================================

# 1. High disease transmission case 
library(nimble)
#library(ggmcmc)
library(coda)

load("~/SCR/R Models y1/Rmodsy1/simTrapGrid.Rdata")
load("~/SCR/R Models y1/Rmodsy1/simHabitatGrid.Rdata")

load("~/OSCR_server/simTrapGrid.Rdata")
load("~/OSCR_server/simHabitatGrid.Rdata")



# skip unnecessary calculations + vectorised detection probabilities 
# When zi = 0, the likelihood is one if the individual was never observed—always the case for augmented individuals—which can be calculated without
# any distances or detection probabilities
dBernoulli_skip_unnecessary = nimbleFunction(
  run = function(x = double(1), z = double(0),  p0 = double(0), sigma = double(0), s = double(1), traplocs = double(2), log = integer(0, default =0) ){
    returnType(double(0))
    
    ## Shortcut if the current individual is not available for detection
    if(z == 0){
      
      prob <- numeric(length(x))
      
    } else {
      
      ## Calculate the log-probability of the vector of detections
      alpha <- -1.0 / (2.0 * sigma * sigma)
      
      D2 <- pow(s[1]- traplocs[ ,1], 2) + pow(s[2]-traplocs[ ,2],2) 
      prob <-  p0*exp( alpha*D2) 
      
    }
    
    logProb <- sum(dbinom(x, prob = prob, size = 1, log = TRUE))
    
    if(log)return(logProb)
    return(exp(logProb))
    
  }
)

#disease data likelihood only when individual is captured
dcat_alive = nimbleFunction(
  run = function(x = double(0), a = double(1),  disease = double(0), p = double(2), log = integer(0, default =0) ){
    returnType(double(0))
    
    ## Shortcut if the current individual is not available for detection
    if( any(a==1) ){
      
      logProb <-  dcat(x , p[disease,1:8], log = TRUE) # disease observation process  
      
    } else {
      
      logProb <- 0 #dcat(x , p[1,1:8], log = TRUE) # d
      
    }
    
    if(log)return(logProb)
    return(exp(logProb))
    
  }
)


rcat_alive = nimbleFunction(
  run = function(n = double(0), a = double(1),  disease = double(0), p = double(2) ){
    returnType(double(0))
    if(n != 1) print("rcat_alive only allows n = 1; using n = 1.")
    y <-  rcat(n , p[disease,1:8]) # disease observation process
    
    return(y)
    
  }
)


Psilden = nimbleFunction(
  run= function(S= double(2), habitatGrid = double(2), Den = double(1), beta0 = double(0), beta1= double(0)){
    
    M = dim(S)[1]
    idpsi = numeric(M)
    
    for (i in 1:M) {
      
      ## EXTRACT LOCATION OF THE ID
      sID <- habitatGrid[trunc(S[i,2])+1, trunc(S[i,1])+1]
      idpsi[i] = ilogit(beta0 + beta1*Den[sID])  
      
    }
    return(idpsi)
    returnType(double(1))
  }
)

#Psilden = compileNimble(Psilden)


dbern_disease <- nimbleFunction(
  run = function(x = double(0),S = double(1), habitatGrid = double(2), beta0 = double(0), beta1= double(0), Den= double(1), z= double(0), d = double(0), log = integer(0, default =0)  ){
    returnType(double(0))
    
    if(z==1 & d ==0){
      sID <- habitatGrid[trunc(S[2])+1, trunc(S[1])+1]
      psi.disease = ilogit(beta0 + beta1*Den[sID])  
      logProb = dbinom( x, prob = psi.disease, size = 1, log = TRUE) 
      
    } else  logProb = dbinom( x, prob = d, size = 1, log = TRUE) 
    
    
    if(log)return(logProb)
    return(exp(logProb))
    
  }
)

rbern_disease <- nimbleFunction(
  run = function(n = double(0),S = double(1), habitatGrid = double(2), beta0 = double(0), beta1= double(0), Den= double(1), z= double(0), d = double(0) ){
    returnType(double(0))
    if(n != 1) print("rbern disease only allows n = 1; using n = 1.")
    y <-  rbinom(1,1,0.5)
    
    return(y)
    
  }
)

binary_state_sampler <- nimbleFunction(
  name = 'binary_state_sampler',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## checks
    if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
    # if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    model[[target]] <<- 1 - model[[target]]
    otherLogProbPrior <- model$calculate(target)
    if(otherLogProbPrior == -Inf) {
      otherLogProb <- otherLogProbPrior
    } else {
      otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
    }
    acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
    jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
    if(jump) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    }
  },
  methods = list(
    reset = function() { }
  )
)


calculateDensity <- nimbleFunction(
  run = function( s = double(2), habitatGrid = double(2) , indicator = double(1), numWindows = double(0), nIndividuals = double(0)){
    # Return type declaration
    returnType(double(1))
    
    dens <- numeric(length = numWindows, value = 0)
    for(i in 1:nIndividuals){
      if(indicator[i]==1){
        windowInd <- habitatGrid[trunc(s[i,2]) + 1, trunc(s[i,1]) + 1]
        dens[windowInd] <- 1 + dens[windowInd]
      }
    }
    
    return(dens)
    
  })


e2dist <- function (x, y)
{ i <- sort(rep(1:nrow(y), nrow(x)))
dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


IdLden = function(S, habitatGrid, lden){
  
  M = dim(S)[1]
  idlen = numeric(M)
  
  for (i in 1:M) {
    
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(S[i,2])+1, trunc(S[i,1])+1]
    idlen[i] =  lden[sID]
    
  }
  
  return(idlen)
  
}


simOSCR.fn <- function(N,phi, gamma ,p0,M, T, delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, grid, habitatGrid, xl, xu, yl, yu, sigma, K){ # M is total ever alive
  
  ntraps <- dim(grid)[1]
  nreps<- K
  pmat <- array(NA, dim = c(M,ntraps, K, T))
  dmat <- array(NA, dim = c(M,ntraps, T))
  Nuninfected <- Ninfected <- matrix(NA, K,T)
  z <- al <- r <- d <- od <- array(0, dim = c(M, K, T ))
  lden <- array(0, dim = c(prod(dim(habitatGrid)), K, T ))
  DenCov <-  array(0, dim = c(M, K, T ))
  surv = numeric(M)
  
  sx <-runif(M,xl,xu)
  sy <-runif(M,yl,yu)
  S <- cbind(sx, sy)
  
  # first year 
  # Initial states 
  r[,1,1] <- rbinom(M,1,gamma[1,1])
  z[ ,1, 1] <- r[,1,1]
  d[,1,1] <- rbinom(M, z[ ,1, 1] ,delta) #  delta - probability an individual is infected at the start of the study.
  od[,1,1] <- d[,1,1]+1
  Ninfected[1,1]  <- sum(z[1:M,1,1]*d[1:M,1,1]) # number of individual that are infected 
  Nuninfected[1,1] <- sum(z[1:M,1,1]*(1-d[1:M,1,1] )) # number of individuals that are uninfected
  
  # number of individual AC in each habitat cell 
  lden[, 1,1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,1, 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
  
  # transition between states 
  for (k in 2:K) {
    
    
    # alive/death/recruit transistion 
    for (i in 1:M) {
      surv[i] <- rbinom(1,1,z[i,k-1,1]*phi[od[i,k-1,1]]) # need to choose phi
      
    }
    
    al[, k,1] <- apply(matrix(z[ , 1:(k-1), 1],nrow=M,byrow=FALSE),1,sum)>0 # was the individual ever recruited before k
    idx <- 1- as.numeric(al[,k,1])
    r[,k,1] <- rbinom(M,idx,gamma[2,1]) #  recruitment  
    z[, k,1] <- surv + r[,k,1]
    
    # disease transition
    lden[, k,1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,k, 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
    DenCov[, k, 1] <-  IdLden(S, habitatGrid, lden[,k,1])
    psi.disease <- Psilden(S, habitatGrid, lden[,k,1], beta0, beta1)
    
    
    
    new_inf <- rbinom(M,1, z[,k,1]*(1-d[ ,k-1,1])*psi.disease) # transition from uninfected to infected 
    d[,k,1] <-   (d[,k-1,1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
    od[,k,1] <- d[,k,1]+1 
    
    #Test whether DenCov have power to identify effect
    #df = data.frame(y=new_inf,den = DenCov[, k, 1] )
    #glm( y~den,data=df,family="binomial")
    
    Ninfected[k,1]  <- sum(z[1:M,k,1]*d[1:M,k,1]) # number of individual that are infected 
    Nuninfected[k,1] <- sum(z[1:M,k,1]*(1-d[1:M,k,1] )) # number of individuals that are uninfected
  }
  
  
  for (t in 2:T){
    
    
    
    # alive/death/recruit transition 
    for (i in 1:M) {
      surv[i] <- rbinom(1,1,z[i,K,t-1]*phi[od[i,K,t-1]]) # need to choose phi
      
    }
    
    al[,1,t] <- apply(matrix(z[ , 1:K, 1:(t-1)],nrow=M,byrow=FALSE),1,sum)>0  # was the individual ever recruited
    idx <- 1- as.numeric(al[,1,t])
    r[,1,t] <- rbinom(M,idx,gamma[1,t])  #  recruitment  
    z[, 1,t] <- surv + r[,1,t]
    
    # disease transition
    lden[, 1,t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,1, t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
    DenCov[, 1, t] <-  IdLden(S, habitatGrid, lden[,1,t])
    psi.disease <- Psilden(S, habitatGrid, lden[,1,t], beta0, beta1)
    
    new_inf <- rbinom(M,1,  z[, 1,t]*(1-d[,K,t-1])*psi.disease) # transition from uninfected to infected 
    d[,1,t] <- (d[,K,t-1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
    od[,1,t] <- d[,1,t]+1 
    Ninfected[1,t]  <- sum(z[1:M,1,t]*d[1:M,1,t]) # number of individual that are infected 
    Nuninfected[1,t] <- sum(z[1:M,1,t]*(1-d[1:M,1,t])) # number of individuals that are uninfected
    
    for (k in 2:K) {
      
      # alive/death/recruit transistion 
      for (i in 1:M) {
        surv[i] <- rbinom(1,1,z[i,k-1,t]*phi[od[i,k-1,t]]) # need to choose phi
        
      }
      
      al[, k,t] <- (al[,1,t] + apply(matrix(z[ , 1:(k-1), t],nrow=M,byrow=FALSE),1,sum)) >0 # was the individual ever recruited
      idx <- 1- as.numeric(al[,k,t])
      r[,k,t] <- rbinom(M,idx,gamma[2,t]) #  recruitment  
      z[, k,t] <- surv + r[,k,t]
      
      # disease transition
      lden[, k,t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,k, t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
      DenCov[, k, t] <-  IdLden(S, habitatGrid, lden[,k,t])
      psi.disease <- Psilden(S, habitatGrid, lden[,k,t], beta0, beta1)
      
      new_inf <- rbinom(M,1,z[, k,t]*(1-d[,k-1,t])*psi.disease) # transition from uninfected to infected 
      d[,k,t] <- (d[,k-1,t]) + new_inf # infected remain infected and newly infected individuals are accounted for 
      od[,k,t] <- d[,k,t]+1 
      
      Ninfected[k,t]  <- sum(z[1:M,k,t]*d[1:M,k,t]) # number of individual that are infected 
      Nuninfected[k,t] <- sum(z[1:M,k,t]*(1-d[1:M,k,t] )) # number of individuals that are uninfected
      
    }
    
  }
  
  for (t in 1:T){
    dmat[,,t] <- e2dist(S,grid)
    
    for (i in 1:M) {
      for (k in 1:K) {
        pmat[i, ,k,t] <- p0[od[i,k,t]]*exp(-(1/(2*sigma[od[i,k, t]]^2))*dmat[i,,t]*dmat[i,,t]) 
      }
    }
    
  }
  
  
  Nmat= matrix(0, K,T)
  for (t in 1:T) {
    for (k in 1:K) {
      Nmat[k,t] = sum(z[1:M, k, t])
    }
    
  }
  #Nmat   
  
  Rmat= matrix(0, K,T)
  for (t in 1:T) {
    for (k in 1:K) {
      Rmat[k,t] = sum(r[1:M, k, t])
    }
    
  }
  
  {
    p <- matrix(NA, nrow = 2, ncol = 8) 
    p[1,1] <-  sp1*sp2*sp3   #probability of being uninfected and giving result 1 (all negative) 
    p[1,2] <- (1-sp1)*sp2*sp3
    p[1,3] <-  sp1*(1-sp2)*sp3
    p[1,4] <-  sp1*sp2*(1-sp3)
    p[1,5] <- (1-sp1)*(1-sp2)*sp3
    p[1,6] <- (1-sp1)*sp2*(1-sp3)
    p[1,7] <- sp1*(1-sp2)*(1-sp3)
    p[1,8] <- (1-sp1)*(1-sp2)*(1-sp3)  
    
    p[2,1] <-  (1-se1)*(1-se2)*(1-se3)   #probability of being infected and giving result 1 (all negative) 
    p[2,2] <- se1*(1-se2)*(1-se3)
    p[2,3] <- (1-se1)*se2*(1-se3)
    p[2,4] <- (1-se1)*(1-se2)*se3
    p[2,5] <- se1*se2*(1-se3)
    p[2,6] <- se1*(1-se2)*se3
    p[2,7] <- (1-se1)*se2*se3
    p[2,8] <- se1*se2*se3  
  } 
  
  
  
  y<-array(0,dim=c(M,ntraps,K, T))
  obs_disease <- array(NA, dim = c(M,K,T))
  for(i in 1:M){
    for(t in 1:T){
      for (k in 1:K) {
        obs_disease[i,k, t] <- rcat(1, p[od[i,k,t], ] ) # disease testing observation process
        y[i,1:ntraps ,k,t]<-rbinom(ntraps,1, pmat[i, 1:ntraps,k, t]*z[i,k,t] ) # per occasion 
      }
    }
    
  }
  ycapt=y[which(rowSums(y[, , , ])>0), , ,]
  disease_observe = obs_disease[which(rowSums(y[, , , ])>0), ,]
  list(y=ycapt, obs_disease = disease_observe, ninfect= Ninfected, nunifect = Nuninfected, z=z,r=r,gamma=gamma,N=Nmat,R=Rmat, S=S, lden = lden, DenCov= DenCov)
}


# Trap locations
head(grid)
#save(grid, file = "simTrapGrid.Rdata")
J=dim(grid)[1]

# Define state space - 11x11 grid
xl<- 0
yl<- 0
xu<- 11
yu<- 11

#plot(grid, xlim=c(xl,xu), ylim=c(yl,yu), pch=3, cex=0.75, xlab = "sx", ylab = "sy") # trapping grid 

T <- 2 #number of years/seasons
K <- 4 # number of sampling occasions
sigma=c(0.5,1)
p0= c(0.5, 0.2)
N= 200 # N=80; M= 200 works
M= 1000#1000#500
phi= c(0.9,0.8)#c(0.9,0.8)
phi[1]; phi[2]
#gamma = matrix(c(0.15,0.005,0.10,0.005,0.1,0.005,0.1,0.005, 0.08,0.005,0.06,0.005,0.005,0.005,0.005) ,2,T);gamma
gamma = matrix(c(0.40,0.1,0.20,0.15) ,2,T);gamma
psi.disease = 0.15 # prob of  
beta0 = -2.5
beta1 = 0.25
delta= 0.15
se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999

nsim = 5
numbers = 1:nsim

# matrices to store results. 
psi_ldenc_time <- matrix(0, nsim, 3)
N_sum = array(0, dim = c(K*T,5, nsim))
Ninf_sum = array(0, dim = c(K*T,5, nsim))
Nunif_sum = array(0, dim = c(K*T,5, nsim))
gamma11 <- matrix(0,nsim,5)
gamma12 <- matrix(0,nsim,5)
gamma21 <- matrix(0,nsim,5)
gamma22 <- matrix(0,nsim,5)
phi_uninfected <- matrix(0,nsim,5)
phi_infected <- matrix(0,nsim,5)
beta0_sum <- matrix(0,nsim,5)
beta1_sum <- matrix(0,nsim,5)
delta_sum <- matrix(0,nsim,5)
p0_uninfected <- matrix(0,nsim,5)
p0_infected <- matrix(0,nsim,5)
sigma_uninfected <- matrix(0,nsim,5)
sigma_infected <- matrix(0,nsim,5)
se1_sum <- matrix(0,nsim,5)
se2_sum <- matrix(0,nsim,5)
se3_sum <- matrix(0,nsim,5)
sp1_sum <- matrix(0,nsim,5)
sp2_sum <- matrix(0,nsim,5)
sp3_sum <- matrix(0,nsim,5)
ess <- matrix(0, nsim , 43 )
converge = matrix(0, nsim, 19)
colnames(ess) = c("N[1, 1]", "N[2, 1]",  "N[3, 1]",  "N[4, 1]",  "N[1, 2]",  "N[2, 2]", "N[3, 2]" , "N[4, 2]" , "Nnegative[1, 1]", "Nnegative[2, 1]", "Nnegative[3, 1]", "Nnegative[4, 1]", "Nnegative[1, 2]", "Nnegative[2, 2]",  "Nnegative[3, 2]", "Nnegative[4, 2]",  "Npostive[1, 1]",  "Npostive[2, 1]", 
                  "Npostive[3, 1]",  "Npostive[4, 1]",  "Npostive[1, 2]",  "Npostive[2, 2]",  "Npostive[3, 2]",  "Npostive[4, 2]", "beta0", "beta1", "delta_I", "gamma[1, 1]", "gamma[2, 1]",   "gamma[1, 2]" , "gamma[2, 2]" , "p0[1]" , "p0[2]", "phi[1]",  "phi[2]", "se1", "se2",  "se3","sigma[1]", "sigma[2]",  "sp1", "sp2", "sp3")


colnames(converge) = c( "beta0", "beta1","delta_I", "gamma[1, 1]", "gamma[2, 1]",   "gamma[1, 2]" , "gamma[2, 2]" , "p0[1]" , "p0[2]", "phi[1]",  "phi[2]", "se1", "se2",  "se3","sigma[1]", "sigma[2]",  "sp1", "sp2", "sp3")


trials = rep(1, J )


for (i in 1:nsim) {
  
  print(i)
  set.seed(numbers[i])
  simdat<- simOSCR.fn(N,phi,gamma ,p0,M, T, delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, grid,habitatGrid, xl, xu, yl, yu, sigma, K)
  simdat$N
  simdat$ninfect
  simdat$nunifect
  # simdat$Den
  # simdat$psi_transit
  # colSums(simdat$psi_transit, na.rm = T)
  # 
  
  (off_den = mean(simdat$lden))
  
  ntot=dim(simdat$y)[1] #total ever observed in this simulated dataset
  #add M-ntot zero encounter histories (data augmentation)
  dataug <- array(0, dim=c(M, J, K, T))
  dataug[1:ntot, , ,] <- simdat$y
  dim(dataug)
  
  zinit<- array(0, dim = c(M,K,T))
  zinit[1:ntot, ,]<-1
  
  sst <- cbind(runif(M, xl, xu), runif(M, yl, yu))
  dim(sst)
  
  dest= array(0, dim= c(M, K, T))   #matrix(0, M,T)
  dest[1:ntot, ,]<-1
  dim(dest)
  
  disease_data = array(-1, dim = c(M,K,T))
  
  for (t in 1:T) {
    for (k in 1:K) {
      for (j in 1:ntot) {
        if(sum(dataug[j,,k,t])>0 ){
          disease_data[j,k,t] <- simdat$obs_disease[j,k,t] 
        } else  disease_data[j,k,t] <- -1
      }
    }
  }
  
  # head(disease_data)
  # tail(disease_data)
  
  area <- (xu - xl) * (yu - yl)
  
  ntraps= J
  
  win.data = list( M= M, n.cells = prod(dim(habitatGrid)), habitatGrid= habitatGrid, y.maxDet = dim(habitatGrid)[1], x.maxDet = dim(habitatGrid)[2],  T =T, J = ntraps, traplocs = as.matrix(grid), K=K, xlim= xl, xupp= xu, ylim = yl, yupp = yu, area = area, trials = trials)
  str(win.data)
  
  if(i==1){
    
    
    jsmmod_markov = nimbleCode({
      
      # priors
      # psi.disease ~ dunif(0,1) # probability of transitioning from unifected to infected
      delta_I ~ dunif(0,1) # probability of being infected at the start of the study.
      
      se1 ~ dbeta(127.02,131.12) # sensitivity for statpak testting
      sp1 ~ dbeta(10.22,1.68) # specificity for stat pak
      
      se2 ~ dbeta(26.41,7) # ifn 
      sp2 ~ dbeta(9.95,1.61)
      
      se3 ~ dbeta(2.25,12.26) # culture
      sp3 ~ dbeta(60.61,1.06)
      
      
      for (t in 1:T) {
        for (k in 1:2) {
          gamma[k,t] ~ dunif(0,1)     
        }
      }
      
      
      for (s in 1:2) {
        
        alpha0[s] ~ dnorm(0, 0.1)
        logit(phi[s]) <- alpha0[s]
        
        alpha1[s] ~ dnorm(0, 0.1)
        logit(p0[s]) <- alpha1[s]
        sigma[s] ~ dunif(0,15)
        #  alpha2[s] <- 1/(2*sigma[s]^2)
      }
      
      
      # need to run this bit only for real analysis only 
      for (t in 1:T) {
        
        for (k in 1:K) {
          # Model summaries
          N[k,t] <- sum(z[1:M,k,t])    # N per year  
          # D[k,t] <- N[k,t]/ area        #  density per year 
          Rc[k,t] <- sum(R[1:M,k,t])    #number recruited 
          Npostive[k,t]  <- sum(z[1:M,k,t]*disease[1:M,k,t]) # number of individual that are infected 
          Nnegative[k,t] <- sum(z[1:M,k,t]*(1-disease[1:M,k,t] )) # number of individuals that are uninfected 
          # Dpositive[k,t] <- Npostive[t]/area # density of individuals that are infected 
          #  Dnegative[k,t] <- Nnegative[t]/area # density of individuals that are uninfected 
          
        }
      }
      
      # disease2; 1- uninfected , 2= infected
      p[1,1] <-  sp1*sp2*sp3   #probability of being uninfected and giving result 1 (all negative) 
      p[1,2] <- (1-sp1)*sp2*sp3
      p[1,3] <-  sp1*(1-sp2)*sp3
      p[1,4] <-  sp1*sp2*(1-sp3)
      p[1,5] <- (1-sp1)*(1-sp2)*sp3
      p[1,6] <- (1-sp1)*sp2*(1-sp3)
      p[1,7] <- sp1*(1-sp2)*(1-sp3)
      p[1,8] <- (1-sp1)*(1-sp2)*(1-sp3)  
      
      p[2,1] <-  (1-se1)*(1-se2)*(1-se3)   #probability of being infected and giving result 1 (all negative) 
      p[2,2] <- se1*(1-se2)*(1-se3)
      p[2,3] <- (1-se1)*se2*(1-se3)
      p[2,4] <- (1-se1)*(1-se2)*se3
      p[2,5] <- se1*se2*(1-se3)
      p[2,6] <- se1*(1-se2)*se3
      p[2,7] <- (1-se1)*se2*se3
      p[2,8] <- se1*se2*se3  
      
      
      beta0 ~ dnorm(0, 0.1)
      beta1 ~ dnorm(0, 0.1)
      
      #  Compute local density at each grid at each time point
      for (k in 2:K) {
        dens[1:n.cells,k,  1] <- calculateDensity( s = S[1:M,1:2], habitatGrid = habitatGrid[1:y.maxDet, 1:x.maxDet], indicator = z[1:M,k,1]  , numWindows = n.cells , nIndividuals = M)
      }
      
      # 2nd year overlap
      for (t in 2:T) {
        
        dens[1:n.cells,1,  t] <- calculateDensity( s = S[1:M,1:2], habitatGrid = habitatGrid[1:y.maxDet, 1:x.maxDet], indicator = z[1:M,1,t]  , numWindows = n.cells , nIndividuals = M)
        
        for (k in 2:K) {
          
          dens[1:n.cells,k,  t] <- calculateDensity( s = S[1:M,1:2], habitatGrid = habitatGrid[1:y.maxDet, 1:x.maxDet], indicator = z[1:M,k,t]  , numWindows = n.cells , nIndividuals = M)
          
        }
      }
      
      
      # first occasion
      for (i in 1:M){    #loop over M individuals (includes the augmented data) 
        
        # first year
        # 1st entry into the population
        z[i,1,1] ~ dbin(gamma[1,1], 1) # initial state
        a[i,1,1] <-(1-z[i,1,1]) 
        R[i,1,1] <- z[i,1,1]      # Calculate recruits
        
        # independent activity centers upon first entry into the population
        S[i,1] ~ dunif(xlim, xupp)  # set priors for the X and Y coordinates of each individual  
        S[i,2] ~ dunif(ylim, yupp)
        
        disease[i,1,1] ~ dbern(z[i,1,1]*delta_I) # disease: 0= unifected, 1 =  infected
        disease2[i,1,1] <- disease[i,1,1] + 1 # disease2; 1- uninfected , 2= infected
        
        
        for (k in 2:K) {
          
          disease[i,k, 1] ~ dbern_disease(S[i, 1:2],habitatGrid[1:y.maxDet, 1:x.maxDet], beta0, beta1, dens[1:n.cells,k,  1] - off_den, z[i,k,1], disease[i,k-1,1])#dbern(disease[i,k-1,1]+ (z[i,k,1]*(1-disease[i,k-1,1])*psi.disease[i, k,1]))  # an individual infected at t-1 will remain infected and an individual not infected at t-1 will become infected with probability psi.disease at t. 
          disease2[i,k,1] <- disease[i,k,1] + 1 # disease2; 1- uninfected , 2= infected
          
          a1[i,k,1] <- sum(z[i, 1:k , 1])  # have you ever been alive (0 = no, >1 = yes), basically a check to see if individual was ever recruited in the population    
          a[i,k,1] <- 1-nimStep(a1[i,k,1] - 1)   #use the step function to make a1 binary. a indicates if individual is available to be recruited at occasion k,time t.  
          mu[i,k,1] <- (phi[disease2[i,k-1,1]]*z[i,k-1,1]) + (gamma[2,1]*a[i,k-1,1]) 
          z[i,k,1] ~ dbern(mu[i,k,1]) 
          R[i,k,1] <- z[i,k,1]*a[i,k-1, 1] #
          
        }
        
        
        # other years
        for (t in 2:T) {
          
          disease[i,1, t] ~ dbern_disease(S[i, 1:2],habitatGrid[1:y.maxDet, 1:x.maxDet], beta0, beta1, dens[1:n.cells,1,  t] - off_den, z[i,1,t], disease[i,K,t-1])#dbern(disease[i,K,t-1]+ (z[i,1,t]*(1-disease[i,K,t-1])*psi.disease[i, 1,t]))  # an individual infected at t-1 will remain infected and an individual not infected at t-1 will become infected with probability psi.disease at t. 
          disease2[i,1,t] <- disease[i,1,t] + 1 # disease2; 1- uninfected , 2= infected
          
          # Model for transition between years
          a1[i,1,t] <- sum(z[i, 1:K , 1:(t-1)]) + z[i,1,t]  # have you ever been alive (0 = no, >1 = yes), basically a check to see if individual was ever recruited in the population    
          a[i,1,t] <- 1- nimStep(a1[i,1,t] - 1)   #use the step function to make a1 binary. a indicates if individual is available to be recruited at occasion k,time t.  
          mu[i,1,t] <- (phi[disease2[i,K,t-1]]*z[i,K,t-1]) + (gamma[1,t]*a[i,K,t-1]) 
          z[i,1,t] ~ dbern(mu[i,1,t]) 
          R[i,1,t] <- z[i,1,t]*a[i,K, t-1] # whether individual i is recruited at time t; a[i,K,t] =0 and z[i,1,t] =1 - individual is recruited,
          
          for (k in 2:K) {
            
            disease[i,k, t] ~ dbern_disease(S[i, 1:2],habitatGrid[1:y.maxDet, 1:x.maxDet], beta0, beta1, dens[1:n.cells,k,  t] - off_den, z[i,k,t], disease[i,k-1,t])#dbern(disease[i,k-1,t]+ (z[i,k,t]*(1-disease[i,k-1,t])*psi.disease[i, k,t]))  # an individual infected at t-1 will remain infected and an individual not infected at t-1 will become infected with probability psi.disease at t. 
            disease2[i,k,t] <- disease[i,k,t] + 1 # disease2; 1- uninfected , 2= infected
            
            a1[i,k,t] <- sum(z[i, 1:K , 1:(t-1)]) + sum(z[i, 1:k , t])  # have you ever been alive (0 = no, >1 = yes), basically a check to see if individual was ever recruited in the population    
            a[i,k,t] <- 1-nimStep(a1[i,k,t] - 1)   #use the step function to make a1 binary. a indicates if individual is available to be recruited at occasion k,time t.  
            mu[i,k,t] <- (phi[disease2[i,k-1,t]]*z[i,k-1,t]) + (gamma[2,t]*a[i,k-1,t]) 
            z[i,k,t] ~ dbern(mu[i,k,t]) 
            R[i,k,t] <- z[i,k,t]*a[i,k-1, t]
            
          }
        }
        
        
        for(t in 1:T){ 
          
          
          # D2[i,1:J,t] <- pow(S[i,1]- traplocs[1:J,1], 2) + pow(S[i,2]-traplocs[1:J,2],2) 
          
          for(k in 1:K) {    #loop over all traps 
            
            obs_disease[i,k,t] ~ dcat_alive(y[i,1:J,k,t], disease2[i,k,t] , p[1:2,1:8])#dcat(p[disease2[i,k,t],1:8]) # disease observation process
            
            #detect[i,1:J,k,t] <- detection(z[i,k,t], p0[disease2[i,k,t]], sigma[disease2[i,k,t]], s=S[i,1:2], traplocs = traplocs[1:J,1:2]) 
            #y[i,1:J,k,t] ~ dBernoulliVector(detect[i,1:J,k,t], trials[1:J])
            y[i,1:J,k,t] ~ dBernoulli_skip_unnecessary( z[i,k,t], p0[disease2[i,k,t]], sigma[disease2[i,k,t]], s=S[i,1:2], traplocs = traplocs[1:J,1:2] )
            
            
          }
          
        }
        
        
      } 
      
      
      
    })
    
    set.seed(i)
    oscr = nimbleModel(jsmmod_markov, constants =  win.data, data = list(y = dataug, obs_disease = disease_data, off_den = off_den), inits = list( z=zinit, S = sst, disease = dest, beta0 = beta0-0.25, beta1 = beta1+0.1 , alpha0 = c(11.7,1.1), alpha1=c(0.4, -0.8), sigma = sigma+0.5, gamma = gamma+0.1 ,delta_I = delta+0.1 , se1 = rbeta(1,127,131) ,sp1=rbeta(1,10,2), se2=rbeta(1,26,7),sp2 = rbeta(1,9,2), se3= rbeta(1,2,12), sp3= rbeta(1,60,1) ), calculate =T) 
    oscr$calculate()
    # oscr$calculate("y")
    # oscr$calculate("obs_disease")
    # oscr$calculate("disease")
    # 
    #length(oscr$getDependencies("D[2,1]"))
    
    # oscr$disease[1,2,1]
    # oscr$disease[1,1,1]
    # oscr$disp[1,2,1]
    # oscr$psi.disease[1,2,1]
    # 
    coscr <- compileNimble(oscr) 
    oscroconf = configureMCMC(oscr)  
    oscroconf$resetMonitors()
    oscroconf$addMonitors(c("phi", "p0", "gamma", "sigma", "beta0", "beta1", "delta_I","se1", "sp1", "se2", "sp2", "se3", "sp3", "N", "Npostive", "Nnegative"))
    
    oscroconf$removeSamplers(c('alpha1[1]', 'alpha1[2]', 'sigma[1]', 'sigma[2]'))
    oscroconf$addSampler(target = c('alpha1[1]','sigma[1]'), type = 'RW_block')
    oscroconf$addSampler(target = c('alpha1[2]','sigma[2]'), type = 'RW_block')
    
    
    oscroconf$removeSamplers('disease', print = FALSE)
    ## get the vector of nodes for the new sampler
    dnodes <- oscr$expandNodeNames("disease")
    
    ## add a sampler for each znode
    for(dnode in dnodes) oscroconf$addSampler(target = dnode, type = binary_state_sampler, print=FALSE)
    #oscroconf$printSamplers("disease")
    
    
    ## remove the default samplers for s
    oscroconf$removeSamplers('S', print = FALSE)
    
    ## add block samplers for pairs of coordinates
    ## first get the pairs:
    sNodePairs <- split( matrix(oscr$expandNodeNames('S'), ncol = 2), 1:M )
    
    for(s in seq_along(sNodePairs)) oscroconf$addSampler(target = sNodePairs[[s]], type = 'RW_block', control = list(adaptScaleOnly = TRUE), print=FALSE, silent = TRUE)
    
    oscromcmc = buildMCMC(oscroconf)
    oscromod = compileNimble(oscromcmc, project = oscr, resetFunctions = TRUE)
    #oscromcmc.out = runMCMC(oscromod, niter = 25000, nburnin = 15000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE)
    lden_time <- system.time( oscromcmc.out <- runMCMC(oscromod, niter = 25000, nburnin = 15000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    #lden_time <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 400, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    
    psi_ldenc_time[i,1:3] = lden_time[1:3]
    
    N_sum[,,i] <- oscromcmc.out$summary$all.chains[1:(K*T),] 
    Nunif_sum[,,i] <- oscromcmc.out$summary$all.chains[9:16,] 
    Ninf_sum[,,i] <- oscromcmc.out$summary$all.chains[17:24,] 
    delta_sum[i, ] <- oscromcmc.out$summary$all.chains["delta_I",]
    beta0_sum[i,] <- oscromcmc.out$summary$all.chains["beta0",]
    beta1_sum[i,] <- oscromcmc.out$summary$all.chains["beta1",]
    gamma11[i, ] <- oscromcmc.out$summary$all.chains["gamma[1, 1]",] # gamma[1:2,t]
    gamma21[i, ] <- oscromcmc.out$summary$all.chains["gamma[2, 1]",]
    gamma12[i, ] <- oscromcmc.out$summary$all.chains["gamma[1, 2]",]
    gamma22[i, ] <- oscromcmc.out$summary$all.chains["gamma[2, 2]",]
    phi_infected[i, ] <- oscromcmc.out$summary$all.chains["phi[2]",]
    phi_uninfected[i, ] <- oscromcmc.out$summary$all.chains["phi[1]",]
    p0_infected[i,] <- oscromcmc.out$summary$all.chains["p0[2]",]
    p0_uninfected[i,] <- oscromcmc.out$summary$all.chains["p0[1]",]
    sigma_infected[i, ] <- oscromcmc.out$summary$all.chains["sigma[2]",]
    sigma_uninfected[i, ] <- oscromcmc.out$summary$all.chains["sigma[1]",]
    se1_sum[i, ] <- oscromcmc.out$summary$all.chains["se1",] 
    se2_sum[i, ] <- oscromcmc.out$summary$all.chains["se2",] 
    se3_sum[i, ] <- oscromcmc.out$summary$all.chains["se3",] 
    sp1_sum[i, ] <- oscromcmc.out$summary$all.chains["sp1",] 
    sp2_sum[i, ] <- oscromcmc.out$summary$all.chains["sp2",] 
    sp3_sum[i, ] <- oscromcmc.out$summary$all.chains["sp3",] 
    ess[i,] = effectiveSize(oscromcmc.out$samples)
    a = gelman.diag(oscromcmc.out$samples[, 25:43])
    converge[i,] <- a$psrf[,1] # converges
    
    save(psi_ldenc_time, file= "psilden_time_ldenm.Rdata")
    save(N_sum, file= "psilden_Nsum_ldenm.Rdata")
    save(Ninf_sum, file= "psilden_Ninfsum_psildenm.Rdata")
    save(Nunif_sum, file= "psilden_Nunifsum_psildenm.Rdata")
    save(delta_sum, file= "psilden_delta_ldenm.Rdata")
    save(beta0_sum, file = "psilden_beta0_ldenm.Rdata")
    save(beta1_sum, file = "psilden_beta1_ldenm.Rdata")
    save(gamma11, file = "psilden_gamma11_ldenm.Rdata")
    save(gamma21, file = "psilden_gamma21_ldenm.Rdata")
    save(gamma12, file = "psilden_gamma12_ldenm.Rdata")
    save(gamma22, file = "psilden_gamma22_ldenm.Rdata")
    save(phi_infected, file = "psilden_phi_infected_ldenm.Rdata")
    save(phi_uninfected, file = "psilden_phi_uninfected_ldenm.Rdata")
    save(p0_infected, file = "psilden_p0_infected_ldenm.Rdata")
    save(p0_uninfected, file = "psilden_p0_uninfected_ldenm.Rdata")
    save(sigma_infected, file = "psilden_sigma_infected_ldenm.Rdata")
    save(sigma_uninfected, file = "psilden_sigma_uninfected_ldenm.Rdata")
    save(se1_sum, file= "psilden_se1_ldenm.Rdata")
    save(se2_sum, file= "psilden_se2_ldenm.Rdata")
    save(se3_sum, file= "psilden_se3_ldenm.Rdata")
    save(sp1_sum, file= "psilden_sp1_ldenm.Rdata")
    save(sp2_sum, file= "psilden_sp2_ldenm.Rdata")
    save(sp3_sum, file= "psilden_sp3_ldenm.Rdata")
    save(ess, file = "psilden_dess_ldenm.Rdata")
    save(converge, file= "psilden_converge_ldenm.Rdata")
    
    
    
  } else {
    
    coscr$y <- dataug
    coscr$obs_disease  <- disease_data
    coscr$off_den = off_den
    coscr$z = zinit
    coscr$S = sst
    coscr$disease = dest
    coscr$beta0 = beta0-0.25; coscr$beta1= beta1+0.1
    coscr$alpha0 = c(11.7,1.1)
    coscr$alpha1= c(0.4, -0.8)
    coscr$sigma = sigma+0.5
    coscr$gamma = gamma+0.1
    coscr$delta_I = delta+0.1
    coscr$se1 = rbeta(1,127,131) ; coscr$sp1=rbeta(1,10,2); coscr$se2=rbeta(1,26,7); coscr$sp2 = rbeta(1,9,2); coscr$se3= rbeta(1,2,12); coscr$sp3= rbeta(1,60,1)
    
    #lden_time <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 20000, nburnin = 10000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    lden_time <- system.time( oscromcmc.out <- runMCMC(oscromod, niter = 25000, nburnin = 15000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    
    psi_ldenc_time[i,1:3] = lden_time[1:3]
    
    N_sum[,,i] <- oscromcmc.out$summary$all.chains[1:(K*T),] 
    Nunif_sum[,,i] <- oscromcmc.out$summary$all.chains[9:16,] 
    Ninf_sum[,,i] <- oscromcmc.out$summary$all.chains[17:24,] 
    delta_sum[i, ] <- oscromcmc.out$summary$all.chains["delta_I",]
    beta0_sum[i,] <- oscromcmc.out$summary$all.chains["beta0",]
    beta1_sum[i,] <- oscromcmc.out$summary$all.chains["beta1",]
    gamma11[i, ] <- oscromcmc.out$summary$all.chains["gamma[1, 1]",] # gamma[1:2,t]
    gamma21[i, ] <- oscromcmc.out$summary$all.chains["gamma[2, 1]",]
    gamma12[i, ] <- oscromcmc.out$summary$all.chains["gamma[1, 2]",]
    gamma22[i, ] <- oscromcmc.out$summary$all.chains["gamma[2, 2]",]
    phi_infected[i, ] <- oscromcmc.out$summary$all.chains["phi[2]",]
    phi_uninfected[i, ] <- oscromcmc.out$summary$all.chains["phi[1]",]
    p0_infected[i,] <- oscromcmc.out$summary$all.chains["p0[2]",]
    p0_uninfected[i,] <- oscromcmc.out$summary$all.chains["p0[1]",]
    sigma_infected[i, ] <- oscromcmc.out$summary$all.chains["sigma[2]",]
    sigma_uninfected[i, ] <- oscromcmc.out$summary$all.chains["sigma[1]",]
    se1_sum[i, ] <- oscromcmc.out$summary$all.chains["se1",] 
    se2_sum[i, ] <- oscromcmc.out$summary$all.chains["se2",] 
    se3_sum[i, ] <- oscromcmc.out$summary$all.chains["se3",] 
    sp1_sum[i, ] <- oscromcmc.out$summary$all.chains["sp1",] 
    sp2_sum[i, ] <- oscromcmc.out$summary$all.chains["sp2",] 
    sp3_sum[i, ] <- oscromcmc.out$summary$all.chains["sp3",] 
    ess[i,] = effectiveSize(oscromcmc.out$samples)
    a = gelman.diag(oscromcmc.out$samples[, 25:43])
    converge[i,] <- a$psrf[,1] # converges
    
    save(psi_ldenc_time, file= "psilden_time_ldenm.Rdata")
    save(N_sum, file= "psilden_Nsum_ldenm.Rdata")
    save(Ninf_sum, file= "psilden_Ninfsum_psildenm.Rdata")
    save(Nunif_sum, file= "psilden_Nunifsum_psildenm.Rdata")
    save(delta_sum, file= "psilden_delta_ldenm.Rdata")
    save(beta0_sum, file = "psilden_beta0_ldenm.Rdata")
    save(beta1_sum, file = "psilden_beta1_ldenm.Rdata")
    save(gamma11, file = "psilden_gamma11_ldenm.Rdata")
    save(gamma21, file = "psilden_gamma21_ldenm.Rdata")
    save(gamma12, file = "psilden_gamma12_ldenm.Rdata")
    save(gamma22, file = "psilden_gamma22_ldenm.Rdata")
    save(phi_infected, file = "psilden_phi_infected_ldenm.Rdata")
    save(phi_uninfected, file = "psilden_phi_uninfected_ldenm.Rdata")
    save(p0_infected, file = "psilden_p0_infected_ldenm.Rdata")
    save(p0_uninfected, file = "psilden_p0_uninfected_ldenm.Rdata")
    save(sigma_infected, file = "psilden_sigma_infected_ldenm.Rdata")
    save(sigma_uninfected, file = "psilden_sigma_uninfected_ldenm.Rdata")
    save(se1_sum, file= "psilden_se1_ldenm.Rdata")
    save(se2_sum, file= "psilden_se2_ldenm.Rdata")
    save(se3_sum, file= "psilden_se3_ldenm.Rdata")
    save(sp1_sum, file= "psilden_sp1_ldenm.Rdata")
    save(sp2_sum, file= "psilden_sp2_ldenm.Rdata")
    save(sp3_sum, file= "psilden_sp3_ldenm.Rdata")
    save(ess, file = "psilden_dess_ldenm.Rdata")
    save(converge, file= "psilden_converge_ldenm.Rdata") 
    
  }
  
}

#=====================================================================================================================================================================================


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Posterior inference results 
library(nimble)
#-----------------------------------------------------------------------------------------------------------------------------------------------------

{
  T <- 2 #number of years/seasons
  K <- 4 # number of sampling occasions
  sigma=c(0.5,1)
  p0= c(0.5, 0.2)
  N= 200 # N=80; M= 200 works
  M= 1000#1000#500
  phi= c(0.9,0.8)#c(0.9,0.8)
  phi[1]; phi[2]
  #gamma = matrix(c(0.15,0.005,0.10,0.005,0.1,0.005,0.1,0.005, 0.08,0.005,0.06,0.005,0.005,0.005,0.005) ,2,T);gamma
  gamma = matrix(c(0.40,0.1,0.20,0.15) ,2,T);gamma
  psi.disease = 0.15 # prob of  
  beta0 = -2.5
  beta1 = 0.25
  delta= 0.15
  se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
  sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999
  
  nsim = 10
  numbers = 1:nsim
  
  sumSim <- function(samples, true_value, nsim, gamma = F ){
    lam_sum = samples
    lambda = true_value
    if (gamma == TRUE){
      mean1 = mean(lam_sum[,1]) # mean
      return(mean1)
    } else{
      med = median(lam_sum[,2]) # mean
      cov = mean(lam_sum[,4] < lambda & lam_sum[,5]>lambda) # coverage
      rmse= sqrt(mean((lam_sum[,2]-lambda)^2))/lambda
      
      lam = lam_sum[,2]
      N= nsim
      lam_bias = CV_bias= CI_width= numeric(N)
      
      for (i in 1:N) {
        lam_bias[i] = (lam[i] - lambda)/lambda
        CV_bias[i] = lam_sum[i,3]/ abs(lam_sum[i,2])
        CI_width[i] = lam_sum[i,5] - lam_sum[i,4]
      }
      
      bias = matrix(0, 1,3)
      colnames(bias) = c("Median", "2.5% Quantile", "97.5% Quantile")
      bias[1,1] = median(lam_bias)
      bias[, 2:3] = quantile(lam_bias, c(0.025, 0.975))
      
      
      CV = matrix(0, 1,3)
      colnames(CV) = c("Median", "2.5% Quantile", "97.5% Quantile")
      CV[1,1] = median(CV_bias)
      CV[, 2:3] = quantile(CV_bias, c(0.025, 0.975))
      
      CIwidth= matrix(0, 1,3)
      colnames(CIwidth) = c("Median", "2.5% Quantile", "97.5% Quantile")
      CIwidth[1,1] = median(CI_width)
      CIwidth[, 2:3] = quantile(CI_width, c(0.025, 0.975))
      
      
      out = list(med = med, cov = cov, rmse = rmse, bias = bias, CV = CV, CIwidth= CIwidth)
      return(out)
      
    }
    
  }
  
}


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_gamma21_ldenM1K.Rdata")
eg = gamma21
tg = gamma[2,1]; tg

sumSim(eg, tg, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_phi_infected_ldenM1K.Rdata")
p <- phi_infected
ph <- phi[2]

sumSim(p, ph, nsim)



load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_phi_uninfected_ldenM1K.Rdata")
p <- phi_uninfected
ph <- phi[1]

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sigma_infected_ldenM1K.Rdata")
ph = sigma[2]
p = sigma_infected

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sigma_uninfected_ldenM1K.Rdata")
ph = sigma[1]
p = sigma_uninfected

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_p0_infected_ldenM1K.Rdata")

p <- p0_infected
ph <- p0[2]

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_p0_uninfected_ldenM1K.Rdata")
p <- p0_uninfected
ph <- p0[1]

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_delta_ldenM1K.Rdata")
ph = delta
p = delta_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_beta0_ldenM1K.Rdata")
ph = beta0
p = beta0_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_beta1_ldenM1K.Rdata")
ph = beta1
p = beta1_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se1_ldenM1K.Rdata")

ph = se1
p = se1_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se2_ldenM1K.Rdata")

ph = se2
p = se2_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se3_ldenM1K.Rdata")

ph = se3
p = se3_sum

sumSim(p, ph, nsim)


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp1_ldenM1K.Rdata")

ph = sp1
p = sp1_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp2_ldenM1K.Rdata")

ph = sp2
p = sp2_sum

sumSim(p, ph, nsim)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp3_ldenM1K.Rdata")

ph = sp3
p = sp3_sum

sumSim(p, ph, nsim)


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_converge_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_ess_ldenM1K.Rdata")

converge
#ess
colMeans(ess)


#---------------------------------------------------------------------------------------------------------------------------------------

# Posterior catepillar plots 


{
  T <- 2 #number of years/seasons
  K <- 4 # number of sampling occasions
  sigma=c(0.5,1)
  p0= c(0.5, 0.2)
  N= 200 # N=80; M= 200 works
  M= 1000#1000#500
  phi= c(0.9,0.8)#c(0.9,0.8)
  phi[1]; phi[2]
  #gamma = matrix(c(0.15,0.005,0.10,0.005,0.1,0.005,0.1,0.005, 0.08,0.005,0.06,0.005,0.005,0.005,0.005) ,2,T);gamma
  gamma = matrix(c(0.40,0.1,0.20,0.15) ,2,T);gamma
  psi.disease = 0.15 # prob of  
  beta0 = -2.5
  beta1 = 0.25
  delta= 0.15
  se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
  sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999
  
  nsim = 10
  numbers = 1:nsim
}


library(ggmcmc)
library(coda)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_phi_infected_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_phi_uninfected_ldenM1K.Rdata")

sam_phi = as.mcmc(cbind(phi_uninfected[,2], phi_infected[,2]))
colnames(sam_phi) = c("uninfected", "infected")
head(sam_phi)

surv = data.frame(value = c(phi[2],phi[1]), Parameter = c("infected", "uninfected"))
ggs_caterpillar(ggs(sam_phi), horizontal = F, sort = F, thick_ci = c(0.05, 0.95),thin_ci = c(0.5, 0.5)) + geom_point(data = surv, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Survival Probability"))

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_gamma11_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_gamma12_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_gamma21_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_gamma22_ldenM1K.Rdata")

sam_gamma = as.mcmc(cbind(gamma11[,2], gamma21[,2],gamma12[,2],gamma22[,2] ))
colnames(sam_gamma) = c("gamma11", "gamma21", "gamma12", "gamma22")
head(sam_gamma)

ga = data.frame(value = c(gamma), Parameter = c("gamma11", "gamma21", "gamma12", "gamma22") )
ggs_caterpillar(ggs(sam_gamma, keep_original_order = TRUE), horizontal = F, sort = F, thick_ci = c(0.05, 0.95),thin_ci = c(0.5, 0.5)) + geom_point(data = ga, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Recruitment Probability"))

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_p0_infected_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_p0_uninfected_ldenM1K.Rdata")

sam_p0 = as.mcmc(cbind(p0_uninfected[,2], p0_infected[,2]))
colnames(sam_p0) = c("uninfected", "infected")
head(sam_p0)

surv = data.frame(value = c(p0[2],p0[1]), Parameter = c("infected", "uninfected"))
ggs_caterpillar(ggs(sam_p0), horizontal = F, sort = F) + geom_point(data = surv, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Baseline detection probabiliy"))

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sigma_infected_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sigma_uninfected_ldenM1K.Rdata")


sam_p0 = as.mcmc(cbind(sigma_uninfected[,2], sigma_infected[,2]))
colnames(sam_p0) = c("uninfected", "infected")
head(sam_p0)

surv = data.frame(value = c(sigma[2],sigma[1]), Parameter = c("infected", "uninfected"))
ggs_caterpillar(ggs(sam_p0), horizontal = F, sort = F) + geom_point(data = surv, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Sigma"))

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_delta_ldenM1K.Rdata")


#load("~/SCR/R Models y1/Rmodsy1/efficient_KOSCR_simulation_results/eff_psidisease.Rdata")
delta_sum
sam_dis = as.mcmc(delta_sum[,2])
#colnames(sam_dis) = c("delta", "psi.disease")
head(sam_dis)

boxplot(delta_sum[,2], xlab = "Delta")
points(delta, col = "red", pch= 15, cex = 2)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_beta0_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_beta1_ldenM1K.Rdata")


sam_beta = as.mcmc(cbind(beta0_sum[,2], beta1_sum[,2]))
colnames(sam_beta) = c("beta0", "beta1")
#head(sam_p0)

surv = data.frame(value = c(beta0,beta1), Parameter = c("beta0", "beta1"))
ggs_caterpillar(ggs(sam_beta), horizontal = F, sort = F) + geom_point(data = surv, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Beta's"))


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se1_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se2_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_se3_ldenM1K.Rdata")

sam_gamma = as.mcmc(cbind(se1_sum[,2], se2_sum[,2], se3_sum[,2] ))
colnames(sam_gamma) = c("se1", "se2", "se3")
head(sam_gamma)

ga = data.frame(value = c(se1,se2, se3), Parameter = c("se1", "se2", "se3") )
ggs_caterpillar(ggs(sam_gamma), horizontal = F, sort = F) + geom_point(data = ga, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Sensitivity"))

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp1_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp2_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_sp3_ldenM1K.Rdata")

sam_gamma = as.mcmc(cbind(sp1_sum[,2], sp2_sum[,2], sp3_sum[,2] ))
colnames(sam_gamma) = c("sp1", "sp2", "sp3")
head(sam_gamma)

ga = data.frame(value = c(sp1,sp2, sp3), Parameter = c("sp1", "sp2", "sp3") )
ggs_caterpillar(ggs(sam_gamma), horizontal = F, sort = F) + geom_point(data = ga, mapping = aes(y= Parameter, x = value), shape = 23, fill="red", size =4 )+ ylab(c("Specificity"))

#=====================================================================================================================================================

# Comparing N's

{
  library(nimble)
  library(ggmcmc)
  library(coda)
  
  
  load("~/SCR/R Models y1/Rmodsy1/simTrapGrid.Rdata")
  load("~/SCR/R Models y1/Rmodsy1/simHabitatGrid.Rdata")
  
  
  Psilden = nimbleFunction(
    run= function(S= double(2), habitatGrid = double(2), Den = double(1), beta0 = double(0), beta1= double(0)){
      
      M = dim(S)[1]
      idpsi = numeric(M)
      
      for (i in 1:M) {
        
        ## EXTRACT LOCATION OF THE ID
        sID <- habitatGrid[trunc(S[i,2])+1, trunc(S[i,1])+1]
        idpsi[i] = ilogit(beta0 + beta1*Den[sID])  
        
      }
      return(idpsi)
      returnType(double(1))
    }
  )
  
  #Psilden = compileNimble(Psilden)
  
  
  
  calculateDensity <- nimbleFunction(
    run = function( s = double(2), habitatGrid = double(2) , indicator = double(1), numWindows = double(0), nIndividuals = double(0)){
      # Return type declaration
      returnType(double(1))
      
      dens <- numeric(length = numWindows, value = 0)
      for(i in 1:nIndividuals){
        if(indicator[i]==1){
          windowInd <- habitatGrid[trunc(s[i,2]) + 1, trunc(s[i,1]) + 1]
          dens[windowInd] <- 1 + dens[windowInd]
        }
      }
      
      return(dens)
      
    })
  
  
  e2dist <- function (x, y)
  { i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  
  
  IdLden = function(S, habitatGrid, lden){
    
    M = dim(S)[1]
    idlen = numeric(M)
    
    for (i in 1:M) {
      
      ## EXTRACT LOCATION OF THE ID
      sID <- habitatGrid[trunc(S[i,2])+1, trunc(S[i,1])+1]
      idlen[i] =  lden[sID]
      
    }
    
    return(idlen)
    
  }
  
  
  simOSCR.fn <- function(N,phi, gamma ,p0,M, T, delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, grid, habitatGrid, xl, xu, yl, yu, sigma, K){ # M is total ever alive
    
    ntraps <- dim(grid)[1]
    nreps<- K
    pmat <- array(NA, dim = c(M,ntraps, K, T))
    dmat <- array(NA, dim = c(M,ntraps, T))
    Nuninfected <- Ninfected <- matrix(NA, K,T)
    z <- al <- r <- d <- od <- array(0, dim = c(M, K, T ))
    lden <- array(0, dim = c(prod(dim(habitatGrid)), K, T ))
    DenCov <-  array(0, dim = c(M, K, T ))
    surv = numeric(M)
    
    sx <-runif(M,xl,xu)
    sy <-runif(M,yl,yu)
    S <- cbind(sx, sy)
    
    # first year 
    # Initial states 
    r[,1,1] <- rbinom(M,1,gamma[1,1])
    z[ ,1, 1] <- r[,1,1]
    d[,1,1] <- rbinom(M, z[ ,1, 1] ,delta) #  delta - probability an individual is infected at the start of the study.
    od[,1,1] <- d[,1,1]+1
    Ninfected[1,1]  <- sum(z[1:M,1,1]*d[1:M,1,1]) # number of individual that are infected 
    Nuninfected[1,1] <- sum(z[1:M,1,1]*(1-d[1:M,1,1] )) # number of individuals that are uninfected
    
    # number of individual AC in each habitat cell 
    lden[, 1,1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,1, 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
    
    # transition between states 
    for (k in 2:K) {
      
      
      # alive/death/recruit transistion 
      for (i in 1:M) {
        surv[i] <- rbinom(1,1,z[i,k-1,1]*phi[od[i,k-1,1]]) # need to choose phi
        
      }
      
      al[, k,1] <- apply(matrix(z[ , 1:(k-1), 1],nrow=M,byrow=FALSE),1,sum)>0 # was the individual ever recruited before k
      idx <- 1- as.numeric(al[,k,1])
      r[,k,1] <- rbinom(M,idx,gamma[2,1]) #  recruitment  
      z[, k,1] <- surv + r[,k,1]
      
      # disease transition
      lden[, k,1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,k, 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
      DenCov[, k, 1] <-  IdLden(S, habitatGrid, lden[,k,1])
      psi.disease <- Psilden(S, habitatGrid, lden[,k,1], beta0, beta1)
      
      
      
      new_inf <- rbinom(M,1, z[,k,1]*(1-d[ ,k-1,1])*psi.disease) # transition from uninfected to infected 
      d[,k,1] <-   (d[,k-1,1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
      od[,k,1] <- d[,k,1]+1 
      
      #Test whether DenCov have power to identify effect
      #df = data.frame(y=new_inf,den = DenCov[, k, 1] )
      #glm( y~den,data=df,family="binomial")
      
      Ninfected[k,1]  <- sum(z[1:M,k,1]*d[1:M,k,1]) # number of individual that are infected 
      Nuninfected[k,1] <- sum(z[1:M,k,1]*(1-d[1:M,k,1] )) # number of individuals that are uninfected
    }
    
    
    for (t in 2:T){
      
      
      
      # alive/death/recruit transition 
      for (i in 1:M) {
        surv[i] <- rbinom(1,1,z[i,K,t-1]*phi[od[i,K,t-1]]) # need to choose phi
        
      }
      
      al[,1,t] <- apply(matrix(z[ , 1:K, 1:(t-1)],nrow=M,byrow=FALSE),1,sum)>0  # was the individual ever recruited
      idx <- 1- as.numeric(al[,1,t])
      r[,1,t] <- rbinom(M,idx,gamma[1,t])  #  recruitment  
      z[, 1,t] <- surv + r[,1,t]
      
      # disease transition
      lden[, 1,t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,1, t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
      DenCov[, 1, t] <-  IdLden(S, habitatGrid, lden[,1,t])
      psi.disease <- Psilden(S, habitatGrid, lden[,1,t], beta0, beta1)
      
      new_inf <- rbinom(M,1,  z[, 1,t]*(1-d[,K,t-1])*psi.disease) # transition from uninfected to infected 
      d[,1,t] <- (d[,K,t-1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
      od[,1,t] <- d[,1,t]+1 
      Ninfected[1,t]  <- sum(z[1:M,1,t]*d[1:M,1,t]) # number of individual that are infected 
      Nuninfected[1,t] <- sum(z[1:M,1,t]*(1-d[1:M,1,t])) # number of individuals that are uninfected
      
      for (k in 2:K) {
        
        # alive/death/recruit transistion 
        for (i in 1:M) {
          surv[i] <- rbinom(1,1,z[i,k-1,t]*phi[od[i,k-1,t]]) # need to choose phi
          
        }
        
        al[, k,t] <- (al[,1,t] + apply(matrix(z[ , 1:(k-1), t],nrow=M,byrow=FALSE),1,sum)) >0 # was the individual ever recruited
        idx <- 1- as.numeric(al[,k,t])
        r[,k,t] <- rbinom(M,idx,gamma[2,t]) #  recruitment  
        z[, k,t] <- surv + r[,k,t]
        
        # disease transition
        lden[, k,t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ ,k, t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
        DenCov[, k, t] <-  IdLden(S, habitatGrid, lden[,k,t])
        psi.disease <- Psilden(S, habitatGrid, lden[,k,t], beta0, beta1)
        
        new_inf <- rbinom(M,1,z[, k,t]*(1-d[,k-1,t])*psi.disease) # transition from uninfected to infected 
        d[,k,t] <- (d[,k-1,t]) + new_inf # infected remain infected and newly infected individuals are accounted for 
        od[,k,t] <- d[,k,t]+1 
        
        Ninfected[k,t]  <- sum(z[1:M,k,t]*d[1:M,k,t]) # number of individual that are infected 
        Nuninfected[k,t] <- sum(z[1:M,k,t]*(1-d[1:M,k,t] )) # number of individuals that are uninfected
        
      }
      
    }
    
    for (t in 1:T){
      dmat[,,t] <- e2dist(S,grid)
      
      for (i in 1:M) {
        for (k in 1:K) {
          pmat[i, ,k,t] <- p0[od[i,k,t]]*exp(-(1/(2*sigma[od[i,k, t]]^2))*dmat[i,,t]*dmat[i,,t]) 
        }
      }
      
    }
    
    
    Nmat= matrix(0, K,T)
    for (t in 1:T) {
      for (k in 1:K) {
        Nmat[k,t] = sum(z[1:M, k, t])
      }
      
    }
    #Nmat   
    
    Rmat= matrix(0, K,T)
    for (t in 1:T) {
      for (k in 1:K) {
        Rmat[k,t] = sum(r[1:M, k, t])
      }
      
    }
    
    {
      p <- matrix(NA, nrow = 2, ncol = 8) 
      p[1,1] <-  sp1*sp2*sp3   #probability of being uninfected and giving result 1 (all negative) 
      p[1,2] <- (1-sp1)*sp2*sp3
      p[1,3] <-  sp1*(1-sp2)*sp3
      p[1,4] <-  sp1*sp2*(1-sp3)
      p[1,5] <- (1-sp1)*(1-sp2)*sp3
      p[1,6] <- (1-sp1)*sp2*(1-sp3)
      p[1,7] <- sp1*(1-sp2)*(1-sp3)
      p[1,8] <- (1-sp1)*(1-sp2)*(1-sp3)  
      
      p[2,1] <-  (1-se1)*(1-se2)*(1-se3)   #probability of being infected and giving result 1 (all negative) 
      p[2,2] <- se1*(1-se2)*(1-se3)
      p[2,3] <- (1-se1)*se2*(1-se3)
      p[2,4] <- (1-se1)*(1-se2)*se3
      p[2,5] <- se1*se2*(1-se3)
      p[2,6] <- se1*(1-se2)*se3
      p[2,7] <- (1-se1)*se2*se3
      p[2,8] <- se1*se2*se3  
    } 
    
    
    
    y<-array(0,dim=c(M,ntraps,K, T))
    obs_disease <- array(NA, dim = c(M,K,T))
    for(i in 1:M){
      for(t in 1:T){
        for (k in 1:K) {
          obs_disease[i,k, t] <- rcat(1, p[od[i,k,t], ] ) # disease testing observation process
          y[i,1:ntraps ,k,t]<-rbinom(ntraps,1, pmat[i, 1:ntraps,k, t]*z[i,k,t] ) # per occasion 
        }
      }
      
    }
    ycapt=y[which(rowSums(y[, , , ])>0), , ,]
    disease_observe = obs_disease[which(rowSums(y[, , , ])>0), ,]
    list(y=ycapt, obs_disease = disease_observe, ninfect= Ninfected, nunifect = Nuninfected, z=z,r=r,gamma=gamma,N=Nmat,R=Rmat, S=S, lden = lden, DenCov= DenCov)
  }
  
  
  # Trap locations
  head(grid)
  #save(grid, file = "simTrapGrid.Rdata")
  J=dim(grid)[1]
  
  # Define state space - 11x11 grid
  xl<- 0
  yl<- 0
  xu<- 11
  yu<- 11
  
  #plot(grid, xlim=c(xl,xu), ylim=c(yl,yu), pch=3, cex=0.75, xlab = "sx", ylab = "sy") # trapping grid 
  
  T <- 2 #number of years/seasons
  K <- 4 # number of sampling occasions
  sigma=c(0.5,1)
  p0= c(0.5, 0.2)
  N= 200 # N=80; M= 200 works
  M= 1000#1000#500
  phi= c(0.9,0.8)#c(0.9,0.8)
  phi[1]; phi[2]
  #gamma = matrix(c(0.15,0.005,0.10,0.005,0.1,0.005,0.1,0.005, 0.08,0.005,0.06,0.005,0.005,0.005,0.005) ,2,T);gamma
  gamma = matrix(c(0.40,0.1,0.20,0.15) ,2,T);gamma
  psi.disease = 0.15 # prob of  
  beta0 = -2.5
  beta1 = 0.25
  delta= 0.15
  se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
  sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999
  
  nsim = 10
  numbers = 1:nsim
  
  
}


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_Nsum_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_Ninf_sum_ldenM1K.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/SCR/LocalDensity_SimResults/M1000_high/psilden_Nunifsum_ldenM1K.Rdata")

N_bias = N_cov = NCV = NCI= matrix(0, nsim, T*K)
Ninf_bias = inf_cov = NinfCV = NinfCI= matrix(0, nsim, T*K)
Nuninf_bias = uninf_cov = NunifCV = NunifCI= matrix(0, nsim, T*K)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  # simulate data from overlap model
  simdat<- simOSCR.fn(N,phi,gamma ,p0,M, T, delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, grid,habitatGrid, xl, xu, yl, yu, sigma, K)
  Ninfected = c(simdat$ninfect)
  Nuninfected = c(simdat$nunifect)
  N =  c(simdat$N)
  for (j in 1:(T*K)) {
    
    #bias
    N_bias[i,j] =  (N_sum[j, 2 ,i] - N[j]) / N[j] 
    Ninf_bias[i,j] =  (Ninf_sum[j, 2 ,i] - Ninfected[j]) / Ninfected[j] 
    Nuninf_bias[i,j] =  (Nunif_sum[j, 2 ,i] - Nuninfected[j]) / Nuninfected[j] 
    
    # CV
    NCV[i,j] = N_sum[j,3,i]/ abs(N_sum[j,2,i])
    NinfCV[i,j] = Ninf_sum[j,3,i]/ abs(Ninf_sum[j,2,i])
    NunifCV[i,j] = Nunif_sum[j,3,i]/ abs(Nunif_sum[j,2,i])
    
    # CI width
    NCI[i,j] =   N_sum[j,5,i] -  N_sum[j,4,i]
    NinfCI[i,j] = Ninf_sum[j, 5 ,i] - Ninf_sum[j, 4 ,i]
    NunifCI[i,j] = Nunif_sum[j, 5 ,i] - Nunif_sum[j, 4 ,i]
    
    # coverage
    if( N_sum[j, 4 ,i] <= N[j] & N_sum[j, 5 ,i] >=N[j]) N_cov[i,j] = 1
    else N_cov[i,j] = 0
    
    if( Ninf_sum[j, 4 ,i] <= Ninfected[j] & Ninf_sum[j, 5 ,i] >=Ninfected[j]) inf_cov[i,j] = 1
    else inf_cov[i,j] = 0
    
    if( Nunif_sum[j, 4 ,i] <= Nuninfected[j] & Nunif_sum[j, 5 ,i] >=Nuninfected[j]) uninf_cov[i,j] = 1
    else uninf_cov[i,j] = 0
    
  }
  
}


Nbias = colMeans(N_bias)
plot(Nbias, pch = 15, xlab = "Occasions")
ggs_caterpillar(ggs(as.mcmc(N_bias)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))
ggs_caterpillar(ggs(as.mcmc(NCV)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))
ggs_caterpillar(ggs(as.mcmc(NCI)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))

# Relative bias summary
N_RBsummay = Ninf_RBsummay= Nunif_RBsummay= array(0, dim=c(1,3,T*K)) # median, 95% qunatile
N_CVsummay = Ninf_CVsummay= Nunif_CVsummay= array(0, dim=c(1,3,T*K)) # median, 95% qunatile
N_CIsummay = Ninf_CIsummay= Nunif_CIsummay= array(0, dim=c(1,3,T*K)) # median, 95% qunatile

for (i in 1:(T*K)) {
  N_RBsummay[,1,i] = median(N_bias[,i])
  N_RBsummay[,2:3,i] = quantile(N_bias[,i], c(0.025,0.975))
  
  N_CVsummay[,1,i] = median(NCV[,i])
  N_CVsummay[,2:3,i] = quantile(NCV[,i], c(0.025,0.975))
  
  N_CIsummay[,1,i] = median(NCI[,i])
  N_CIsummay[,2:3,i] = quantile(NCI[,i], c(0.025,0.975))
  
  N_RBsummay[,1,i] = median(N_bias[,i])
  N_RBsummay[,2:3,i] = quantile(N_bias[,i], c(0.025,0.975))
  
  Ninf_RBsummay[,1,i] = median(Ninf_bias[,i])
  Ninf_RBsummay[,2:3,i] = quantile(Ninf_bias[,i], c(0.025,0.975))
  
  Ninf_CVsummay[,1,i] = median(NinfCV[,i])
  Ninf_CVsummay[,2:3,i] = quantile(NinfCV[,i], c(0.025,0.975))
  
  Ninf_CIsummay[,1,i] = median(NinfCI[,i])
  Ninf_CIsummay[,2:3,i] = quantile(NinfCI[,i], c(0.025,0.975))
  
  Nunif_RBsummay[,1,i] = median(Nuninf_bias[,i])
  Nunif_RBsummay[,2:3,i] = quantile(Nuninf_bias[,i], c(0.025,0.975))
  
  Nunif_CVsummay[,1,i] = median(NunifCV[,i])
  Nunif_CVsummay[,2:3,i] = quantile(NunifCV[,i], c(0.025,0.975))
  
  Nunif_CIsummay[,1,i] = median(NunifCI[,i])
  Nunif_CIsummay[,2:3,i] = quantile(NunifCI[,i], c(0.025,0.975))
  
}


Ncov = colMeans(N_cov)
plot(Ncov, pch = 15, xlab = "Occasions")

NInf_bias = colMeans(Ninf_bias)
plot(NInf_bias, pch = 15, xlab = "Occasions")
ggs_caterpillar(ggs(as.mcmc(Ninf_bias)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))

NInf_cov = colMeans(inf_cov)
plot(NInf_cov, pch = 15, xlab = "Occasions")


NUnInf_bias = colMeans(Nuninf_bias)
plot(NUnInf_bias, pch = 15, xlab = "Occasions")
ggs_caterpillar(ggs(as.mcmc(Nuninf_bias)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))

NUnInf_cov = colMeans(uninf_cov)
plot(NUnInf_cov, pch = 15, xlab = "Occasions")


#====================================================================================================================================================


