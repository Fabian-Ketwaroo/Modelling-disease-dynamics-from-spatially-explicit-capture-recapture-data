
rm(list = ls())
library(nimble)
#library(nimbleSCR)
library(MCMCvis)
library(coda)

dcat_hazard = nimbleFunction(
  run = function(x = double(0), z = double(0),  lambda = double(0), sigma = double(0), s = double(1), traplocs = double(2), id = double(1), log = integer(0, default =0) ){
    returnType(double(0))
    
    h_trap <- nimNumeric(dim(traplocs)[1]+1 ) 
    
    if (z ==0) { 
      
      h_trap[1] <- 1
      logProb <- dcat(x, h_trap, log=TRUE)  #dbinom(x, prob = prob[y], size = 1, log = TRUE) 
    }
    else {
      
      alpha <- -1.0 / (2.0 * sigma * sigma)
      D2 <- pow(s[1]- traplocs[id ,1], 2) + pow(s[2]-traplocs[id ,2],2) 
      pmat <- lambda*exp( alpha*D2) #lambda*exp(-(1/(2*sigma^2))*dmat[i, id]*dmat[i, id])  
      h_star <- sum(pmat)
      h_trap[nimC(1, id+1)] <- nimC(exp(-h_star*z), (1-exp(-h_star*z) )*(pmat/h_star)  ) 
      logProb <- dcat(x, h_trap, log=TRUE)  #dbinom(x, prob = prob[y], size = 1, log = TRUE) 
    }
    
    
    if(log)return(logProb)
    return(exp(logProb))
    
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


dcat_alive = nimbleFunction(
  run = function(x = double(0), a = double(0),  disease = double(0), p = double(2), log = integer(0, default =0) ){
    returnType(double(0))
    
    ## Shortcut if the current individual is not available for detection
    if( a > 1 ){
      
      logProb <-  dcat(x , p[disease,1:8], log = TRUE) # disease observation process  
      
    } else {
      
      logProb <- 0 #dcat(x , p[1,1:8], log = TRUE) # d
      
    }
    
    if(log)return(logProb)
    return(exp(logProb))
    
  }
)

#compileNimble(dcat_alive)

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


e2dist <- function (x, y)
{ i <- sort(rep(1:nrow(y), nrow(x)))
dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

simKOSCR.fn <- function(phi, gamma,lambda,delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3,M, T, R,  grid,habitatGrid, xl, xu, yl, yu, sigma, id, nid){ # M is total ever alive
  
  ntraps <- dim(grid)[1]
  pmat <- array(0, dim = c(M,R, ntraps, T))
  dmat <- array(NA, dim = c(M, ntraps))
  h_trap <- array(0, dim = c(M,R,ntraps+1, T))
  h_star <- array(0, dim = c(M,R, T))
  Nuninfected <- Ninfected <- numeric(T)
  lden <- array(0, dim = c(prod(dim(habitatGrid)), T ))
  DenCov <-  array(0, dim = c(M, T ))
  
  sx <-runif(M,xl,xu)
  sy <-runif(M,yl,yu)
  S <- cbind(sx, sy)
  
  z <- al <- r <-  d <- od <- array(0, dim = c(M, T ))
  r[,1] <- rbinom(M,1,gamma[1]) 
  z[ , 1] <- r[,1]
  d[,1] <- rbinom(M, z[ , 1] ,delta) #  delta - probability an individual is infected at the start of the study.
  od[,1] <- d[,1]+1
  Ninfected[1]  <- sum(z[1:M,1]*d[1:M,1]) # number of individual that are infected 
  Nuninfected[1] <- sum(z[1:M,1]*(1-d[1:M,1] )) # number of individuals that are uninfected
  
  # number of individual AC in each habitat cell 
  lden[, 1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ , 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
  
  
  for (t in 2:T) {
    
    surv <- rbinom(M ,1 , z[ ,t-1]*phi[od[,t-1]])
    al[,t] <- apply(matrix(z[ ,1:(t-1)],nrow=M,byrow=FALSE),1,sum)>0  
    idx <- 1- as.numeric(al[,t])
    r[,t] <- rbinom(M,idx,gamma[t]) #  recruitment  
    z[, t] <- surv + r[,t]
    
    # disease transition
    # local density at each grid
    lden[, t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ , t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
    
    # Local density of each individual.   
    DenCov[, t] <-  IdLden(S, habitatGrid, lden[,t])
    
    # Disease transition probabiliy as a function of local density  
    psi.disease <- Psilden(S, habitatGrid, lden[,t], beta0, beta1)
    
    # disease transition
    new_inf <- rbinom(M, 1, z[,t]*(1-d[,t-1])*psi.disease) # transition from uninfected to infected 
    d[,t] <- (d[,t-1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
    od[,t] <- d[,t]+1 
    Ninfected[t]  <- sum(z[1:M,t]*d[1:M,t]) # number of individual that are infected 
    Nuninfected[t] <- sum(z[1:M,t]*(1-d[1:M,t])) # number of individuals that are uninfected
    
  }
  
  
  
  Nmat= numeric(T)#matrix(0, K,T)
  for (t in 1:T) {
    Nmat[t] = sum(z[1:M, t])
  }
  
  # For the 3 area, different number of traps are active at each time point 
  dmat <- e2dist(S,grid)
  for (i in 1:M) {
    for (t in 1:T) {
      for (r in 1:R) {
        pmat[i, r, id[1:nid[r,t],r,t], t] = lambda[od[i ,t]]*exp(-(1/(2*sigma[od[i ,t]]^2))*dmat[i, id[1:nid[r,t],r,t] ]*dmat[i, id[1:nid[r,t],r,t]])  
        h_star[i,r, t] = sum(pmat[i,r, id[1:nid[r,t],r,t] ,t])
        h_trap[i,r,1 ,t] = exp(-h_star[i,r,t]*z[i,t])  
        h_trap[i,r,id[1:nid[r,t],r,t]+1 ,t] = (1-exp(-h_star[i,r,t]*z[i,t]) )*(pmat[i,r, id[1:nid[r,t],r,t] ,t]/h_star[i,r, t]) 
        
        
      }
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
  
  
  y<-array(0,dim=c(M,R,T))
  obs_disease <- array(NA, dim = c(M,R, T))
  for (t in 1:T) {
    for (i in 1:M) {
      for (r in 1:R) {
        y[i,r,t] <- rcat(1, h_trap[i,r, ,t])
        obs_disease[i,r, t] <- rcat(1, p[od[i,t], ] )
      }
    }
  }
  
  ycapt=y[which(rowSums(y[,,])>T*R),,]
  disease_observe = obs_disease[which(rowSums(y[, , ])>T*R), ,  ]
  list(y=ycapt, obs_disease = disease_observe, h_trap = h_trap, ninfect= Ninfected, nunifect = Nuninfected, z=z,N=Nmat, lden = lden, DenCov= DenCov, S =S, SX=sx, SY=sy)
}


# Load data
load("~/data/id_oscrdata_T2014_2018.Rdata")
load("~/data/simBadgerHabitatGrid.Rdata")
load("~/data/simrescale_traplocs_2014_2018.Rdata")

M = 500
#nind = ntot = 259
T <- 8 # 2 years with 4 sampling occasions in each year
#K <- 4 #number of years/seasons
R <- 3 # 3 areas
sigma= c(0.25,0.6)
lambda = c(1.5,0.25)
#M = 650
phi = c(0.9,0.8) #matrix(c(0.6,0.9,0.8,0.7,0.6,0.9,0.7,0.8  ) ,K,T);phi
#gamma = c(0.4, 0.25,0.1, 0.08)
gamma = c(0.4, rep(0.1,3), 0.25, rep(0.15,3))# c(matrix(c(0.40,0.1,0.20,0.15) ,2,2))
#psi.disease = 0.15#0.15 # prob of  
beta0 = -3
beta1 = 0.5# 0.3
delta= 0.15
se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999



nid = matrix(0, nrow = R, ncol = T)

for (t in 1:T) {
  for (r in 1:R) {
    nid[r,t]= sum(id[,r,t] != 0)
  }
}

nid

buffer = 0 #  
ntraps <- nrow(traplocs)
xl <- min(traplocs[, 1] - buffer)
xu <- max(traplocs[, 1] + buffer )
yl <- min(traplocs[, 2] - buffer)
yu <- max(traplocs[, 2] + buffer ) # these values here are already scaled.
#plot(traplocs, xlim = c(Xl,Xu), ylim = c(Yl, Yu) )
area <- (xu - xl) * (yu - yl)/4; area 
# sst <- cbind(runif(M, xl, xu), runif(M, yl, yu))
# dim(sst)

grid = traplocs
J=dim(grid)[1]


nsim = 40
numbers = 1:nsim

# matrices to store results. 
psi_ldenc_time <- matrix(0, nsim, 3)
N_sum = array(0, dim = c(T,5, nsim))
Ninf_sum = array(0, dim = c(T,5, nsim))
Nunif_sum = array(0, dim = c(T,5, nsim))
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
ess <- matrix(0, nsim , 49 )
converge = matrix(0, nsim, 49)
#colnames(ess) = c("N[1]", "N[2]",  "N[3]",  "N[4]",  "N[5]",  "N[6]", "N[7]" , "N[8]" , "Nnegative[1]", "Nnegative[2]", "Nnegative[3]", "Nnegative[4]", "Nnegative[5]", "Nnegative[6]",  "Nnegative[7]", "Nnegative[8]",  "Npostive[1]",  "Npostive[2]", 
#                  "Npostive[3]",  "Npostive[4]",  "Npostive[5]",  "Npostive[6]",  "Npostive[7]",  "Npostive[8]", "beta0", "beta1", "delta_I", "gamma[1, 1]", "gamma[2, 1]",   "gamma[1, 2]" , "gamma[2, 2]" , "p0[1]" , "p0[2]", "phi[1]",  "phi[2]", "se1", "se2",  "se3","sigma[1]", "sigma[2]",  "sp1", "sp2", "sp3")


#colnames(converge) = c( "beta0", "beta1","delta_I", "gamma11", "gamma[2, 1]",   "gamma[1, 2]" , "gamma[2, 2]" , "lambda[1]" , "lambda[2]", "phi[1]",  "phi[2]", "se1", "se2",  "se3","sigma[1]", "sigma[2]",  "sp1", "sp2", "sp3")





for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  mysimdat<- simKOSCR.fn(phi=phi, gamma= gamma, lambda ,delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, M=M, T=T,R = R, grid=grid,habitatGrid, xl=xl, xu=xu, yl=yl, yu=yu, sigma=sigma, id, nid)
  mysimdat$N
  mysimdat$ninfect
  # as time goes by more individuals can become infected. Hence, there there wiil be more infected than unifnected. I guess to counter this can lower the survival of infected individuals.
  mysimdat$nunifect
  
  ntot=dim(mysimdat$y)[1] #total ever observed in this simulated dataset
  ntot
  
  (off_den = mean(mysimdat$lden))
  #add Mc-ntot zero encounter histories (data augmentation)
  dataug <- array(1, dim=c(M, R, T))
  dataug[1:ntot,  ,] <- mysimdat$y
  dim(dataug)
  
  zinit<- array(0, dim = c(M,T))
  zinit[1:ntot ,]<-1
  
  
  SX<-runif(M, xl, xu)
  SY<-runif(M, yl, yu)  
  S= cbind(SX, SY)
  
  # Observed disease data
  disease_data = array(-1, dim = c(M, R, T))
  
  for (t in 1:T) {
    for (r in 1:R) {
      for (j in 1:M) {
        if(dataug[j,r,t]>1) disease_data[j,r,t] <- mysimdat$obs_disease[j,r,t] # individual capture 
        else disease_data[j,r,t] <- -1 # individual not capture so no data available
      }
    }
  }
  
  
  # starting values for disease status
  dest= array(0, dim= c(M, T))   # says individuals to be recruited are uninfected
  #head(disease_data[,,1])
  dim(disease_data)
  #dest[1:ntot ,]<-1 # this assumption that individuals caught are all infected may not be the best place for the MCMC to start
  #dim(dest)
  id_uninfect = 1:4
  id_infect = 5:8
  
  for (t in 1:T) {
    if(t==1){
      for (i in 1:ntot) {
        
        if(all(is.na(disease_data[i, , t]))){  # no disease data 
          dest[i,t] = 0 # assumed to be uninfect so that they can transition 
        } else {
          
          # extract disease data value
          dobs <- disease_data[i, which(!is.na(disease_data[i, ,t])), t] 
          
          # if only one observation
          if(length(dobs) == 1) {
            if(any(dobs == id_uninfect)) dest[i,t] = 0
            else dest[i, t] = 1 
          } else{ # more than one observation
            
            unif <- sum(dobs %in% id_uninfect) # number of suggested  uninfected 
            inf <-  sum(dobs %in% id_infect) # number of suggested infected 
            
            if(unif > inf) dest[i,t] = 0
            if(unif < inf) dest[i,t] = 1
            if(unif == inf) dest[i,t] = sample( c(0,1),1,  0.5)
            
          }
          
        }
        
        
        
      }
    } else {
      
      for (i in 1:ntot) {
        
        if(all(is.na(disease_data[i, , t]))){  # no disease data 
          
          dest[i,t] = dest[i, t-1]# take the previous time point status
          l = dbinom(dest[i,t], 1, dest[i, t-1], log = 0)
          if(l == 0) dest[i,t] = 1  #print("-Inf")
          
        } else {
          
          # extract disease data value
          dobs <- disease_data[i, which(!is.na(disease_data[i, ,t])), t] 
          
          # if only one observation
          if(length(dobs) == 1) {
            if(any(dobs == id_uninfect)) dest[i,t] = 0
            else dest[i, t] = 1 
            l = dbinom(dest[i,t], 1, dest[i, t-1], log = 0)
            if(l == 0) dest[i,t] = 1  #print("-Inf")
          } else{ # more than one observation
            
            unif <- sum(dobs %in% id_uninfect) # number of suggested  uninfected 
            inf <-  sum(dobs %in% id_infect) # number of suggested infected 
            
            if(unif > inf) dest[i,t] = 0
            if(unif < inf) dest[i,t] = 1
            if(unif == inf) dest[i,t] = dest[i,t-1]
            l = dbinom(dest[i,t], 1, dest[i, t-1], log = 0)
            if(l == 0) dest[i,t] = 1  #print("-Inf")
          }
          
        }
        
        
        
      }
      
      
    }
    
    
  }
  


  win.data = list(id = id, nid = nid,  n.cells = prod(dim(habitatGrid)), habitatGrid= habitatGrid, y.maxDet = dim(habitatGrid)[1], x.maxDet = dim(habitatGrid)[2], M= M, T =T,  J = J , traplocs = as.matrix(grid), R1= R, xlim= xl, xupp= xu, ylim = yl, yupp = yu, area = area)
  str(win.data)
  
  if(n==1){
    
    jsmmod_markov = nimbleCode( {
      
      # priors
      beta0 ~ dnorm(0, sd= 1)
      beta1 ~ dnorm(0, sd= 1)
      delta_I ~ dunif(0,1) # probability of being infected at the start of the study.
      
      se1 ~ dbeta(127.02,131.12) # sensitivity for statpak testting
      sp1 ~ dbeta(10.22,1.68) # specificity for stat pak
      
      se2 ~ dbeta(26.41,7) # ifn 
      sp2 ~ dbeta(9.95,1.61)
      
      se3 ~ dbeta(2.25,12.26) # culture
      sp3 ~ dbeta(60.61,1.06)
      
      for (s in 1:2) {
        phi[s] ~ dunif(0,1)
        lambda[s] ~ dgamma(0.1,0.1) # baseline encounter rate
        sigma[s] ~ dunif(0,5)
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
      
      # Recruitment at the beginning of each year
      gamma[1] ~ dunif(0,1)
      gamma[5] ~ dunif(0,1)
      
      # Assume recruitment is constant for other occasions after the 1st
      gamma11 ~ dunif(0,1)
      gamma22 ~ dunif(0,1)
      
      for (t in 2:4) {
        gamma[t] <- gamma11
      }
      
      for (t in 6:T) {
        gamma[t] <- gamma22
      }
      
      for(t in 1:T){
        #gamma[t] ~ dunif(0,1)
        N[t] <- sum(z[1:M,t])    # N per year   
        Npostive[t]  <- sum(z[1:M,t]*disease[1:M,t]) # number of individual that are infected 
        Nnegative[t] <- sum(z[1:M,t]*(1-disease[1:M,t] )) # number of individuals that are uninfected 
        
      }
      
      # 2nd year overlap
      for (t in 1:T) {
        dens[1:n.cells,  t] <- calculateDensity( s = S[1:M,1:2], habitatGrid = habitatGrid[1:y.maxDet, 1:x.maxDet], indicator = z[1:M,t]  , numWindows = n.cells , nIndividuals = M)
        
      }
      
      # state model  
      for (i in 1:M){           #loop over M individuals (includes the augmented data) 
        z[i,1] ~ dbin(gamma[1], 1)
        a[i,1] <-(1-z[i,1])  
        R[i,1] <- z[i,1]      #Calculate recruits
        
        
        S[i,1] ~ dunif(xlim, xupp)  # set priors for the X and Y coordinates of each individual  
        S[i,2] ~ dunif(ylim, yupp) 
        
        disease[i,1] ~ dbern(z[i,1]*delta_I) # disease: 0= unifected, 1 =  infected
        disease2[i,1] <- disease[i,1] + 1 # disease2; 1- uninfected , 2= infected
        
        for(t in 2:T){
          
          a1[i,t] <- sum(z[i, 1:t]) #have you ever been alive (0 = no, >1 = yes)    
          a[i,t] <- 1-step(a1[i,t] - 1)     #use the step function to make a1 binary. a indicates if individual is available to be recruited at time t.  
          mu[i,t] <- (phi[disease2[i,t-1]]*z[i,t-1]) + (gamma[t]*a[i,t-1])   # alive state at t is dependent on phi and gamma
          z[i,t] ~ dbern(mu[i,t])    # transistion  of individual from state t to t+1 
          
          disease[i, t] ~ dbern_disease(S[i, 1:2],habitatGrid[1:y.maxDet, 1:x.maxDet], beta0, beta1, dens[1:n.cells, t] - off_den, z[i,t], disease[i,t-1])#dbern(disease[i,k-1,t]+ (z[i,k,t]*(1-disease[i,k-1,t])*psi.disease[i, k,t]))  # an individual infected at t-1 will remain infected and an individual not infected at t-1 will become infected with probability psi.disease at t. 
          disease2[i,t] <- disease[i,t] + 1 # disease2; 1- uninfected , 2= infected 
          
          
        }
      }
      
      # Observation model
      for (t in 1:T) {
        for (i in 1:M) { # loop over M individuals (includes the augmented data
          for (r in 1:R1) {
            y[i,r,t] ~ dcat_hazard(z[i,t], lambda[disease2[i,t]], sigma[disease2[i,t]], S[i, 1:2], traplocs[1:J,1:2], id[1:nid[r,t], r, t])
            obs_disease[i,r,t] ~ dcat_alive(y[i,r,t], disease2[i,t] , p[1:2,1:8]) # disease observation process
          }
          
        }
      }
      
    })
    
    
    set.seed(n)
    gamma_init = gamma + runif(T, 0,0.15)
    oscr = nimbleModel(jsmmod_markov, constants =  win.data, data = list(y = dataug,  obs_disease = disease_data, off_den = off_den), inits = list(z = zinit, S = S, disease = dest, beta0 = beta0-0.25, beta1 = beta1+0.1  , lambda = runif(2,0,3) ,sigma = runif(2,1,3), phi = runif(2), gamma = gamma_init, gamma11 = gamma[2], gamma22 = gamma[6], delta_I = runif(1,0.1,0.3), se1 = rbeta(1,127,131) ,sp1=rbeta(1,10,2), se2=rbeta(1,26,7),sp2 = rbeta(1,9,2), se3= rbeta(1,2,12), sp3= rbeta(1,60,1)), calculate = FALSE)
    #oscr$calculate()
    # logProb sensitive to sigma
    
    
    coscr <- compileNimble(oscr,showCompilerOutput = FALSE) 
    oscroconf = configureMCMC(oscr,monitors = c( "phi", "lambda", "gamma", "sigma", "beta0", "beta1", "delta_I","se1", "sp1", "se2", "sp2", "se3", "sp3", "N", "Npostive", "Nnegative"))  
    oscroconf$monitors
    
    oscroconf$removeSamplers('disease', print = FALSE)
    ## get the vector of nodes for the new sampler
    dnodes <- oscr$expandNodeNames("disease")
    
    ## add a sampler for each znode
    for(dnode in dnodes) oscroconf$addSampler(target = dnode, type = binary_state_sampler, print=FALSE)
    
    ## remove the default samplers for s
    oscroconf$removeSamplers('S', print = FALSE)
    sNodePairs <- split( matrix(oscr$expandNodeNames('S'), ncol = 2), 1:M )
    ## then add the block samplers
    for(s in seq_along(sNodePairs)) oscroconf$addSampler(target = sNodePairs[[s]], type = 'RW_block', control = list(adaptScaleOnly = TRUE), print=FALSE,silent = TRUE)
    
    
    oscromcmc = buildMCMC(oscroconf)
    oscromod = compileNimble(oscromcmc, project = oscr, resetFunctions = TRUE)
    #lden_time  <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 1000, nburnin = 400, thin = 1, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    #lden_time  <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 10000, nburnin = 4000, thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    lden_time  <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 25000, nburnin = 15000, thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    psi_ldenc_time[n,1:3] = lden_time[1:3]
    
    N_sum[,,n] <- oscromcmc.out$summary$all.chains[1:(T),] 
    Nunif_sum[,,n] <- oscromcmc.out$summary$all.chains[9:16,] 
    Ninf_sum[,,n] <- oscromcmc.out$summary$all.chains[17:24,] 
    delta_sum[n, ] <- oscromcmc.out$summary$all.chains["delta_I",]
    beta0_sum[n,] <- oscromcmc.out$summary$all.chains["beta0",]
    beta1_sum[n,] <- oscromcmc.out$summary$all.chains["beta1",]
    gamma11[n, ] <- oscromcmc.out$summary$all.chains["gamma[1]",] # gamma[1:2,t]
    gamma21[n, ] <- oscromcmc.out$summary$all.chains["gamma[2]",]
    gamma12[n, ] <- oscromcmc.out$summary$all.chains["gamma[5]",]
    gamma22[n, ] <- oscromcmc.out$summary$all.chains["gamma[6]",]
    phi_infected[n, ] <- oscromcmc.out$summary$all.chains["phi[2]",]
    phi_uninfected[n, ] <- oscromcmc.out$summary$all.chains["phi[1]",]
    p0_infected[n,] <- oscromcmc.out$summary$all.chains["lambda[2]",]
    p0_uninfected[n,] <- oscromcmc.out$summary$all.chains["lambda[1]",]
    sigma_infected[n, ] <- oscromcmc.out$summary$all.chains["sigma[2]",]
    sigma_uninfected[n, ] <- oscromcmc.out$summary$all.chains["sigma[1]",]
    se1_sum[n, ] <- oscromcmc.out$summary$all.chains["se1",] 
    se2_sum[n, ] <- oscromcmc.out$summary$all.chains["se2",] 
    se3_sum[n, ] <- oscromcmc.out$summary$all.chains["se3",] 
    sp1_sum[n, ] <- oscromcmc.out$summary$all.chains["sp1",] 
    sp2_sum[n, ] <- oscromcmc.out$summary$all.chains["sp2",] 
    sp3_sum[n, ] <- oscromcmc.out$summary$all.chains["sp3",] 
    ess[n,] = effectiveSize(oscromcmc.out$samples)
    converge[n,] <- MCMCsummary(oscromcmc.out$samples)[,6]
    
    save(psi_ldenc_time, file= "psilden_time_ldenh.Rdata")
    save(N_sum, file= "psilden_Nsum_ldenh.Rdata")
    save(Ninf_sum, file= "psilden_Ninfsum_psildenh.Rdata")
    save(Nunif_sum, file= "psilden_Nunifsum_psildenh.Rdata")
    save(delta_sum, file= "psilden_delta_ldenh.Rdata")
    save(beta0_sum, file = "psilden_beta0_ldenh.Rdata")
    save(beta1_sum, file = "psilden_beta1_ldenh.Rdata")
    save(gamma11, file = "psilden_gamma11_ldenh.Rdata")
    save(gamma21, file = "psilden_gamma21_ldenh.Rdata")
    save(gamma12, file = "psilden_gamma12_ldenh.Rdata")
    save(gamma22, file = "psilden_gamma22_ldenh.Rdata")
    save(phi_infected, file = "psilden_phi_infected_ldenh.Rdata")
    save(phi_uninfected, file = "psilden_phi_uninfected_ldenh.Rdata")
    save(p0_infected, file = "psilden_p0_infected_ldenh.Rdata")
    save(p0_uninfected, file = "psilden_p0_uninfected_ldenh.Rdata")
    save(sigma_infected, file = "psilden_sigma_infected_ldenh.Rdata")
    save(sigma_uninfected, file = "psilden_sigma_uninfected_ldenh.Rdata")
    save(se1_sum, file= "psilden_se1_ldenh.Rdata")
    save(se2_sum, file= "psilden_se2_ldenh.Rdata")
    save(se3_sum, file= "psilden_se3_ldenh.Rdata")
    save(sp1_sum, file= "psilden_sp1_ldenh.Rdata")
    save(sp2_sum, file= "psilden_sp2_ldenh.Rdata")
    save(sp3_sum, file= "psilden_sp3_ldenh.Rdata")
    save(ess, file = "psilden_dess_ldenh.Rdata")
    save(converge, file= "psilden_converge_ldenh.Rdata")
    
  }  else {
    
    coscr$y <- dataug
    coscr$obs_disease  <- disease_data
    coscr$off_den = off_den
    coscr$z = zinit
    coscr$S = S
    coscr$disease = dest
    coscr$beta0 = beta0-0.25; coscr$beta1= beta1+0.1
    coscr$lambda = runif(2,0,3) 
    coscr$sigma = runif(2,1,3)
    coscr$phi = runif(2)
    coscr$gamma = gamma_init
    coscr$gamma11 = gamma_init[2]
    coscr$gamma22 = gamma_init[6]
    coscr$delta_I = delta+0.1
    coscr$se1 = rbeta(1,127,131) ; coscr$sp1=rbeta(1,10,2); coscr$se2=rbeta(1,26,7); coscr$sp2 = rbeta(1,9,2); coscr$se3= rbeta(1,2,12); coscr$sp3= rbeta(1,60,1)
    
    #lden_time <- system.time(oscromcmc.out <- runMCMC(oscromod, niter = 20000, nburnin = 10000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    lden_time <- system.time( oscromcmc.out <- runMCMC(oscromod, niter = 25000, nburnin = 15000,thin = 5, samplesAsCodaMCMC = TRUE, nchains = 2,  summary = TRUE))
    
    psi_ldenc_time[n,1:3] = lden_time[1:3]
    
    N_sum[,,n] <- oscromcmc.out$summary$all.chains[1:(T),] 
    Nunif_sum[,,n] <- oscromcmc.out$summary$all.chains[9:16,] 
    Ninf_sum[,,n] <- oscromcmc.out$summary$all.chains[17:24,] 
    delta_sum[n, ] <- oscromcmc.out$summary$all.chains["delta_I",]
    beta0_sum[n,] <- oscromcmc.out$summary$all.chains["beta0",]
    beta1_sum[n,] <- oscromcmc.out$summary$all.chains["beta1",]
    gamma11[n, ] <- oscromcmc.out$summary$all.chains["gamma[1]",] # gamma[1:2,t]
    gamma21[n, ] <- oscromcmc.out$summary$all.chains["gamma[2]",]
    gamma12[n, ] <- oscromcmc.out$summary$all.chains["gamma[5]",]
    gamma22[n, ] <- oscromcmc.out$summary$all.chains["gamma[6]",]
    phi_infected[n, ] <- oscromcmc.out$summary$all.chains["phi[2]",]
    phi_uninfected[n, ] <- oscromcmc.out$summary$all.chains["phi[1]",]
    p0_infected[n,] <- oscromcmc.out$summary$all.chains["lambda[2]",]
    p0_uninfected[n,] <- oscromcmc.out$summary$all.chains["lambda[1]",]
    sigma_infected[n, ] <- oscromcmc.out$summary$all.chains["sigma[2]",]
    sigma_uninfected[n, ] <- oscromcmc.out$summary$all.chains["sigma[1]",]
    se1_sum[n, ] <- oscromcmc.out$summary$all.chains["se1",] 
    se2_sum[n, ] <- oscromcmc.out$summary$all.chains["se2",] 
    se3_sum[n, ] <- oscromcmc.out$summary$all.chains["se3",] 
    sp1_sum[n, ] <- oscromcmc.out$summary$all.chains["sp1",] 
    sp2_sum[n, ] <- oscromcmc.out$summary$all.chains["sp2",] 
    sp3_sum[n, ] <- oscromcmc.out$summary$all.chains["sp3",] 
    ess[n,] = effectiveSize(oscromcmc.out$samples)
    converge[n,] <- MCMCsummary(oscromcmc.out$samples)[,6]
    
    save(psi_ldenc_time, file= "psilden_time_ldenh.Rdata")
    save(N_sum, file= "psilden_Nsum_ldenh.Rdata")
    save(Ninf_sum, file= "psilden_Ninfsum_psildenh.Rdata")
    save(Nunif_sum, file= "psilden_Nunifsum_psildenh.Rdata")
    save(delta_sum, file= "psilden_delta_ldenh.Rdata")
    save(beta0_sum, file = "psilden_beta0_ldenh.Rdata")
    save(beta1_sum, file = "psilden_beta1_ldenh.Rdata")
    save(gamma11, file = "psilden_gamma11_ldenh.Rdata")
    save(gamma21, file = "psilden_gamma21_ldenh.Rdata")
    save(gamma12, file = "psilden_gamma12_ldenh.Rdata")
    save(gamma22, file = "psilden_gamma22_ldenh.Rdata")
    save(phi_infected, file = "psilden_phi_infected_ldenh.Rdata")
    save(phi_uninfected, file = "psilden_phi_uninfected_ldenh.Rdata")
    save(p0_infected, file = "psilden_p0_infected_ldenh.Rdata")
    save(p0_uninfected, file = "psilden_p0_uninfected_ldenh.Rdata")
    save(sigma_infected, file = "psilden_sigma_infected_ldenh.Rdata")
    save(sigma_uninfected, file = "psilden_sigma_uninfected_ldenh.Rdata")
    save(se1_sum, file= "psilden_se1_ldenh.Rdata")
    save(se2_sum, file= "psilden_se2_ldenh.Rdata")
    save(se3_sum, file= "psilden_se3_ldenh.Rdata")
    save(sp1_sum, file= "psilden_sp1_ldenh.Rdata")
    save(sp2_sum, file= "psilden_sp2_ldenh.Rdata")
    save(sp3_sum, file= "psilden_sp3_ldenh.Rdata")
    save(ess, file = "psilden_dess_ldenh.Rdata")
    save(converge, file= "psilden_converge_ldenh.Rdata")
    
    
  }
    
  
}


#=====================================================================================================================================================================================


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Posterior inference results 
library(nimble)
#-----------------------------------------------------------------------------------------------------------------------------------------------------

{
  M = 500
  #nind = ntot = 259
  T <- 8 # 2 years with 4 sampling occasions in each year
  #K <- 4 #number of years/seasons
  R <- 3 # 3 areas
  sigma= c(0.25,0.6)
  lambda = c(1.5,0.25)
  #M = 650
  phi = c(0.9,0.8) #matrix(c(0.6,0.9,0.8,0.7,0.6,0.9,0.7,0.8  ) ,K,T);phi
  #gamma = c(0.4, 0.25,0.1, 0.08)
  gamma = c(0.4, rep(0.1,3), 0.25, rep(0.15,3))# c(matrix(c(0.40,0.1,0.20,0.15) ,2,2))
  #psi.disease = 0.15#0.15 # prob of  
  beta0 = -3
  beta1 = 0.5# 0.3
  delta= 0.15
  se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
  sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999
  
  nsim = 30
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


load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_gamma11_ldenh.Rdata")
eg = gamma11
tg = gamma[1]; tg

sumSim(eg, tg, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_gamma21_ldenh.Rdata")
eg = gamma21
tg = gamma[2]; tg

sumSim(eg, tg, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_gamma12_ldenh.Rdata")
eg = gamma12
tg = gamma[5]; tg

sumSim(eg, tg, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_gamma22_ldenh.Rdata")
eg = gamma22
tg = gamma[6]; tg

sumSim(eg, tg, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_phi_infected_ldenh.Rdata")
p <- phi_infected
ph <- phi[2]

sumSim(p, ph, nsim)


load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_phi_uninfected_ldenh.Rdata")
p <- phi_uninfected
ph <- phi[1]

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_sigma_infected_ldenh.Rdata")
ph = sigma[2]
p = sigma_infected

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_sigma_uninfected_ldenh.Rdata")
ph = sigma[1]
p = sigma_uninfected

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_p0_infected_ldenh.Rdata")
p <- p0_infected
ph <- lambda[2]

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_p0_uninfected_ldenh.Rdata")
p <- p0_uninfected
ph <- lambda[1]

sumSim(p, ph, nsim)



load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_delta_ldenh.Rdata")
ph = delta
p = delta_sum

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_beta0_ldenh.Rdata")
ph = beta0
p = beta0_sum

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_beta1_ldenh.Rdata")
ph = beta1
p = beta1_sum
p[19,3] = 0.004
sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_se1_ldenh.Rdata")
ph = se1
p = se1_sum

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_se2_ldenh.Rdata")
ph = se2
p = se2_sum
sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_se3_ldenh.Rdata")
ph = se3
p = se3_sum

sumSim(p, ph, nsim)


load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_sp1_ldenh.Rdata")
ph = sp1
p = sp1_sum
sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_sp2_ldenh.Rdata")
ph = sp2
p = sp2_sum

sumSim(p, ph, nsim)

load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_sp3_ldenh.Rdata")
ph = sp3
p = sp3_sum

sumSim(p, ph, nsim)


load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_dess_ldenh.Rdata")
load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_converge_ldenh.Rdata")

converge
#ess
colMeans(ess)


#---------------------------------------------------------------------------------------------------------------------------------------
# Dont need to update these plots but will leave it just in case.
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
  
  
  dcat_alive = nimbleFunction(
    run = function(x = double(0), a = double(0),  disease = double(0), p = double(2), log = integer(0, default =0) ){
      returnType(double(0))
      
      ## Shortcut if the current individual is not available for detection
      if( a > 1 ){
        
        logProb <-  dcat(x , p[disease,1:8], log = TRUE) # disease observation process  
        
      } else {
        
        logProb <- 0 #dcat(x , p[1,1:8], log = TRUE) # d
        
      }
      
      if(log)return(logProb)
      return(exp(logProb))
      
    }
  )
  
  #compileNimble(dcat_alive)
  
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
  
  
  e2dist <- function (x, y)
  { i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  
  simKOSCR.fn <- function(phi, gamma,lambda,delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3,M, T, R,  grid,habitatGrid, xl, xu, yl, yu, sigma, id, nid){ # M is total ever alive
    
    ntraps <- dim(grid)[1]
    pmat <- array(0, dim = c(M,R, ntraps, T))
    dmat <- array(NA, dim = c(M, ntraps))
    h_trap <- array(0, dim = c(M,R,ntraps+1, T))
    h_star <- array(0, dim = c(M,R, T))
    Nuninfected <- Ninfected <- numeric(T)
    lden <- array(0, dim = c(prod(dim(habitatGrid)), T ))
    DenCov <-  array(0, dim = c(M, T ))
    
    sx <-runif(M,xl,xu)
    sy <-runif(M,yl,yu)
    S <- cbind(sx, sy)
    
    z <- al <- r <-  d <- od <- array(0, dim = c(M, T ))
    r[,1] <- rbinom(M,1,gamma[1]) 
    z[ , 1] <- r[,1]
    d[,1] <- rbinom(M, z[ , 1] ,delta) #  delta - probability an individual is infected at the start of the study.
    od[,1] <- d[,1]+1
    Ninfected[1]  <- sum(z[1:M,1]*d[1:M,1]) # number of individual that are infected 
    Nuninfected[1] <- sum(z[1:M,1]*(1-d[1:M,1] )) # number of individuals that are uninfected
    
    # number of individual AC in each habitat cell 
    lden[, 1] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ , 1] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
    
    
    for (t in 2:T) {
      
      surv <- rbinom(M ,1 , z[ ,t-1]*phi[od[,t-1]])
      al[,t] <- apply(matrix(z[ ,1:(t-1)],nrow=M,byrow=FALSE),1,sum)>0  
      idx <- 1- as.numeric(al[,t])
      r[,t] <- rbinom(M,idx,gamma[t]) #  recruitment  
      z[, t] <- surv + r[,t]
      
      # disease transition
      # local density at each grid
      lden[, t] <-  calculateDensity(s = S,habitatGrid = habitatGrid,  indicator = z[ , t] , numWindows = prod(dim(habitatGrid)),  nIndividuals = M)
      
      # Local density of each individual.   
      DenCov[, t] <-  IdLden(S, habitatGrid, lden[,t])
      
      # Disease transition probabiliy as a function of local density  
      psi.disease <- Psilden(S, habitatGrid, lden[,t], beta0, beta1)
      
      # disease transition
      new_inf <- rbinom(M, 1, z[,t]*(1-d[,t-1])*psi.disease) # transition from uninfected to infected 
      d[,t] <- (d[,t-1]) + new_inf # infected remain infected and newly infected individuals are accounted for 
      od[,t] <- d[,t]+1 
      Ninfected[t]  <- sum(z[1:M,t]*d[1:M,t]) # number of individual that are infected 
      Nuninfected[t] <- sum(z[1:M,t]*(1-d[1:M,t])) # number of individuals that are uninfected
      
    }
    
    
    
    Nmat= numeric(T)#matrix(0, K,T)
    for (t in 1:T) {
      Nmat[t] = sum(z[1:M, t])
    }
    
    # For the 3 area, different number of traps are active at each time point 
    dmat <- e2dist(S,grid)
    for (i in 1:M) {
      for (t in 1:T) {
        for (r in 1:R) {
          pmat[i, r, id[1:nid[r,t],r,t], t] = lambda[od[i ,t]]*exp(-(1/(2*sigma[od[i ,t]]^2))*dmat[i, id[1:nid[r,t],r,t] ]*dmat[i, id[1:nid[r,t],r,t]])  
          h_star[i,r, t] = sum(pmat[i,r, id[1:nid[r,t],r,t] ,t])
          h_trap[i,r,1 ,t] = exp(-h_star[i,r,t]*z[i,t])  
          h_trap[i,r,id[1:nid[r,t],r,t]+1 ,t] = (1-exp(-h_star[i,r,t]*z[i,t]) )*(pmat[i,r, id[1:nid[r,t],r,t] ,t]/h_star[i,r, t]) 
          
          
        }
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
    
    
    y<-array(0,dim=c(M,R,T))
    obs_disease <- array(NA, dim = c(M,R, T))
    for (t in 1:T) {
      for (i in 1:M) {
        for (r in 1:R) {
          y[i,r,t] <- rcat(1, h_trap[i,r, ,t])
          obs_disease[i,r, t] <- rcat(1, p[od[i,t], ] )
        }
      }
    }
    
    ycapt=y[which(rowSums(y[,,])>T*R),,]
    disease_observe = obs_disease[which(rowSums(y[, , ])>T*R), ,  ]
    list(y=ycapt, obs_disease = disease_observe, h_trap = h_trap, ninfect= Ninfected, nunifect = Nuninfected, z=z,N=Nmat, lden = lden, DenCov= DenCov, S =S, SX=sx, SY=sy)
  }
  
  
  # Load data
  load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Data/id_oscrdata_T2014_2018.Rdata")
  #load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Data/dcat_oscrdata_T2014_2018.Rdata")
  #load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Data/disease_dpp_data_2014_2018.RData")
  load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/simBadgerHabitatGrid.Rdata")
  load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/simrescale_traplocs_2014_2018.Rdata")

  
  
  #  load("~/SCR_revisions/id_oscrdata_T2014_2018.Rdata")
  # # load("~/SCR_revisions/dcat_oscrdata_T2014_2018.Rdata")
  # # load("~/SCR_revisions/disease_dpp_data_2014_2018.RData")
  #  load("~/SCR_revisions/BadgerHabitatGrid.Rdata")
  #  load("~/SCR_revisions/rescale_traplocs_2014_2018.Rdata")
  
  # noether
  # load("~/SCR/id_oscrdata_T2014_2018.Rdata")
  # load("~/SCR/simBadgerHabitatGrid.Rdata")
  # load("~/SCR/simrescale_traplocs_2014_2018.Rdata")
  
  M = 500
  #nind = ntot = 259
  T <- 8 # 2 years with 4 sampling occasions in each year
  #K <- 4 #number of years/seasons
  R <- 3 # 3 areas
  sigma= c(0.25,0.6)
  lambda = c(1.5,0.25)
  #M = 650
  phi = c(0.9,0.8) #matrix(c(0.6,0.9,0.8,0.7,0.6,0.9,0.7,0.8  ) ,K,T);phi
  #gamma = c(0.4, 0.25,0.1, 0.08)
  gamma = c(0.4, rep(0.1,3), 0.25, rep(0.15,3))# c(matrix(c(0.40,0.1,0.20,0.15) ,2,2))
  #psi.disease = 0.15#0.15 # prob of  
  beta0 = -3
  beta1 = 0.5# 0.3
  delta= 0.15
  se1 <- 0.492; se2 <- 0.809; se3 <- 0.1
  sp1 <- 0.931; sp2 <- 0.936; sp3 <- 0.999
  
  
  
  nid = matrix(0, nrow = R, ncol = T)
  
  for (t in 1:T) {
    for (r in 1:R) {
      nid[r,t]= sum(id[,r,t] != 0)
    }
  }
  
  nid
  
  buffer = 0 #  
  ntraps <- nrow(traplocs)
  xl <- min(traplocs[, 1] - buffer)
  xu <- max(traplocs[, 1] + buffer )
  yl <- min(traplocs[, 2] - buffer)
  yu <- max(traplocs[, 2] + buffer ) # these values here are already scaled.
  #plot(traplocs, xlim = c(Xl,Xu), ylim = c(Yl, Yu) )
  area <- (xu - xl) * (yu - yl)/4; area 
  # sst <- cbind(runif(M, xl, xu), runif(M, yl, yu))
  # dim(sst)
  
  grid = traplocs
  J=dim(grid)[1]
  
  
  nsim = 30
  numbers = 1:nsim
  
  
  
}


load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_Nsum_ldenh.Rdata")
load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_Nunifsum_psildenh.Rdata")
load("C:/Users/fke/OneDrive - Vogelwarte/SCR/Revisions/Simulation_results/psilden_Ninfsum_psildenh.Rdata")


N_bias = N_cov = NCV = NCI= matrix(0, nsim, T)
Ninf_bias = inf_cov = NinfCV = NinfCI= matrix(0, nsim, T)
Nuninf_bias = uninf_cov = NunifCV = NunifCI= matrix(0, nsim, T)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  # simulate data from overlap model
  mysimdat<- simKOSCR.fn(phi=phi, gamma= gamma, lambda ,delta,beta0, beta1, sp1,sp2,sp3, se1,se2,se3, M=M, T=T,R = R, grid=grid,habitatGrid, xl=xl, xu=xu, yl=yl, yu=yu, sigma=sigma, id, nid)
  N = c(mysimdat$N)
  Ninfected  = c(mysimdat$ninfect)
  Nuninfected = c(mysimdat$nunifect)
  
  for (j in 1:(T)) {
    
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
mean(N_bias)
quantile(N_bias, c(0.025, 0.975))
plot(Nbias, pch = 15, xlab = "Occasions")

mean(NCV)
quantile(NCV, c(0.025, 0.975))

mean(Ncov)



mean(Nuninf_bias)
quantile(Nuninf_bias, c(0.025, 0.975))

mean(NunifCV)
quantile(NunifCV, c(0.025, 0.975))
mean(uninf_cov)

mean(Ninf_bias)
quantile(Ninf_bias, c(0.025, 0.975))

mean(NinfCV)
quantile(NinfCV, c(0.025, 0.975))
mean(inf_cov)


ggs_caterpillar(ggs(as.mcmc(N_bias)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))
ggs_caterpillar(ggs(as.mcmc(NCV)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))
ggs_caterpillar(ggs(as.mcmc(NCI)), horizontal = F, sort = F, thick_ci = c(0.025, 0.975),thin_ci = c(0.5, 0.5))

# Relative bias summary
N_RBsummay = Ninf_RBsummay= Nunif_RBsummay= array(0, dim=c(1,3,T)) # median, 95% qunatile
N_CVsummay = Ninf_CVsummay= Nunif_CVsummay= array(0, dim=c(1,3,T)) # median, 95% qunatile
N_CIsummay = Ninf_CIsummay= Nunif_CIsummay= array(0, dim=c(1,3,T)) # median, 95% qunatile

for (i in 1:(T)) {
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


