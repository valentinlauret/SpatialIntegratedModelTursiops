# Spatial capture recapture

# This code run SCR analyze for Bottlenose dolphins in the NW Mediterranean Sea
# It requires "SCR5.rdata" file
# Analyzes are run with NIMBLE, please visit https://r-nimble.org/ for information



#--- load required packages  -----
library(nimble)
library(dplyr)


# ---- function to load data  -----

dataSCR <- function(M = double(0)){
  
  # load data 
  load("SCR/SCR5.rdata")
  
  
  y.scr <- y.scr[,samp.sites,] # subset the detection history matrix to sampled sites only
  
  
  # get UTM coordinates of sampled sites and of all sites of the study area 
  traploc <- cbind(g2$X[samp.sites],g2$Y[samp.sites])
  sites <- cbind(g2$X,g2$Y)
  
  
  # 
  ybis <- y.scr
  
  
  nsites <- nrow(sites)  # number of sites
  # full analysis will be run with nsites = 4356
  nind <- dim(ybis)[1]  # nb of ind = 927
  nocc <- dim(ybis)[3]  # nb of occ = 8
  ntrap <- dim(ybis)[2]  # nb of trap = 327 (i.e. nb of sampled sites)
  
  
  # Data Augmentation 
  M <- M
  if(M < nind) stop()
  
  
  ## DT: here's the main change we're going to make
  ## change the representation of ybis, so instead of a 3-dimensional
  ## array (individuals x traps x occasions) of 0's and 1's,
  ## instead it will be a 3-dimensional array of
  ## (individuals x occasions x trap IDs with captures)
  ## specifically: for individual i, on occasion t, ybis[i, t, 1:MAX] will contain
  ## the trap IDs of these captures (with 0's following the final capture trap ID,
  ## or all 0's if individual i was never captured on occasion t.
  ## we'll have to set the size MAX correctly, based on the maximum
  ## number of sightings of any individual
  
  MAX <- max(apply(ybis, c(1,3), sum))
  ybis2 <- array(0, c(nind, nocc, MAX))
  for(i in 1:nind) {
    for(t in 1:nocc) {
      trapIDs <- which(ybis[i,,t] > 0)
      numSightings <- length(trapIDs)
      if(numSightings > 0) ybis2[i, t, 1:numSightings] <- trapIDs
    }
  }
  
  
  ## DT: new data augmentation using ybis2
  ## augmented individuals just have all 0's in ybis2
  ## we'll call this new variable y.da2, to distinguish it
  
  y.da2 <- array(0, c(M, nocc, MAX))
  y.da2[1:nind,,] <- ybis2
  dim(y.da2)
  max(y.da2) ## must be < ntrap
  
  
  ## DT: create z initial values using y.da2
  z <- as.numeric(apply(y.da2, 1, sum) > 0)
  
  # draw starting locations from the list of sites
  # for known individual, take the mean location where they have been seen
  set.seed(32)   ## DT: moved set.seed() here
  head(sites)
  
  id <- sample(1:nsites, M, replace= TRUE)
  
  for(i in 1:nind) {
    
    if(z[i] == 1) {
      siteVector <- as.numeric(y.da2[i,,1:MAX])
      if(max(siteVector) == 0) stop()
      sitesSeen <- siteVector[which(siteVector>0)]
      id[i] <- samp.sites[median(sitesSeen)]   
      ## if an ind is seen at multiple locations, we kept oen location as the starting activity center
    }
  }
  
  # environmental covariate
  habitat <- as.numeric(g2$bathy.sc)
  
  # effort covariate log-transformed
  eff.cov <- seff[samp.sites,] 
  eff.cov <- eff.cov[,c("occ1", "occ2","occ3","occ4", "occ5", "occ6", "occ7", "occ8")]
  eff.cov$geometry <-  NULL
  eff.cov <- log(eff.cov)
  eff.cov <- as.matrix(eff.cov)
  
  # as_tibble() %>% 
  # select(occ1, occ2,occ3,occ4, occ5, occ6, occ7, occ8) %>%
  # log() %>% 
  # as.matrix()
  
  eff.cov[eff.cov == - Inf] <-  0
  
  # make a binary cov with 1 if sites i is sampled during occ j, and 0 otherwise
  eff.ind <- eff.cov
  eff.ind[eff.ind!=0] <- 1
  
  
  dataSCR <- list(nsites = nsites, nind = nind, nocc = nocc, ntrap = ntrap, traploc = traploc,
                  habitat = habitat, seff = eff.cov, eff.ind = eff.ind, MAX = MAX, M = M, sites = sites, z=z,
                  y.scr = y.da2, id = id)
  
  return(dataSCR)
}

# --- load data ---

M = 5000 # data augmentation parameter
datSCR <- dataSCR(M = M)

# --- NIMBLE model ----
SCRcode <-  nimbleCode({
  
  # priors for the density regression
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  # priors for the sigma 
  sigma ~ dunif(10,10^5)
  sig2 <- 2*sigma*sigma
  # priors for the pbar regression
  p0 ~ dnorm(0,1) 
  p1 ~ dnorm(0,1) 
  # prior psi
  #psi ~ dunif(0,1)
  
  log(mu[1:nsites]) <- mu0 + mu1 * habitat[1:nsites]  # abundance model via IPP 
  
  EN <- sum(mu[1:nsites]) # expected number of individual for the study area
  
  psi <-  EN/M
  
  
  # logit regression of pbar
  logit(pbar[1:ntrap, 1:nocc]) <- p0 + p1 * seff[1:ntrap,1:nocc]
  pbar_mask[1:ntrap, 1:nocc] <- pbar[1:ntrap, 1:nocc]*eff.ind[1:ntrap,1:nocc] 
  
  # cell prob
  probs[1:nsites] <- mu[1:nsites]/EN
  
  # data augmentation
  for(i in 1:M) {
    
    z[i] ~ dbern(psi)  # is i alive or not ?
    
    # induce the uniform distribution of activity centers over the discrete set 1:nsites
    id[i] ~ dunif(0.5, nsites+0.5) 
    
    # location of AC
    SX[i] <- sites[round(id[i]),1]
    SY[i] <- sites[round(id[i]),2]
    
    #zeros's trick
    nll[i] <- - log(probs[round(id[i])])
    zeros[i] ~ dpois(nll[i])
    
    ## DT: likelihood using new dBernoulliVector4 distribution:
    y.scr[i, 1:nocc, 1:MAX] ~ dBernoulliVector5(
      MAX = MAX, ntrap = ntrap, nocc = nocc,
      sigma = sig2, pbar = pbar_mask[1:ntrap,1:nocc], 
      eff.ind = eff.ind[1:ntrap, 1:nocc], SX = SX[i], SY = SY[i],
      traploc = traploc[1:ntrap,1:2], indicator = z[i])
    
  }
  N <- sum(z[1:M]) 
})

# ---- v5 -----
# Val mpodif to account for multiple sightings in same occasion 
dBernoulliVector5 <- nimbleFunction(
  run = function(
    x = double(2), MAX = double(0), ntrap = double(0), nocc = double(0),
    sigma = double(0), pbar = double(2), eff.ind = double(2),  SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0), log = integer(0, default = 0)) {
    if(indicator == 0) { if(sum(x) == 0) return(0) else return(-Inf) }
    lp <- 0
    for(s in 1:ntrap) {
      dist <- (SX-traploc[s,1])^2 + (SY-traploc[s,2])^2
      for(j in 1:nocc) {
        p <- pbar[s,j]*exp(-dist/(sigma))
        for(i in 1:MAX){
          if(x[j,i] == s) {
            lp <- lp + log(p)
          } else { lp <- lp + log(1-p) }
        }
      }
    }
    returnType(double())
    if(log) return(lp)
    return(exp(lp))
  })

rBernoulliVector5 <- nimbleFunction(
  run = function(
    n = integer(), MAX = double(0), ntrap = double(0), nocc = double(0),
    sigma = double(0), pbar = double(2), eff.ind = double(2), SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0)) {
    stop('should never call rBernoulliVector5')
    returnType(double(2))
    return(array(0, c(2, 2)))
  })

registerDistributions(list(
  dBernoulliVector5 = list(
    BUGSdist = 'dBernoulliVector5(MAX,ntrap,nocc,sigma,pbar,eff.ind,SX,SY,traploc,indicator)',
    types = c('value = double(2)', 'sigma = double(0)', 'pbar = double(2)','eff.ind = double(2)', 'traploc = double(2)'),
    mixedSizes = TRUE   ## DT: prevents a warning about size/dimension mismatch
  )))
# ---- Build data ----

# constants
constants <-  list(nsites = datSCR$nsites, nocc = datSCR$nocc, M = datSCR$M, ntrap = datSCR$ntrap,
                   traploc = datSCR$traploc/100000, habitat = datSCR$habitat,
                   seff = datSCR$seff, eff.ind = datSCR$eff.ind,
                   MAX = datSCR$MAX)  

# data
data <- list(sites= datSCR$sites/100000, y.scr = datSCR$y.scr, zeros = rep(0,M))

# Inits
inits <-  list(z=datSCR$z, id = c(ACi, sample(1:4356,1)),#datSCR$id,
               mu0=0,mu1=0,sigma = 1000, p0 = 0, p1= 0)

# NIMBLE RUN ----

# Build model
Rmodel <- nimbleModel(SCRcode, constants, data, inits)
Rmodel$calculate() # - 127631
# configure model
conf <- configureMCMC(Rmodel)

conf$printMonitors()
conf$resetMonitors()
conf$addMonitors(c('EN','mu0','mu1',"sigma", "p0","p1", "id","z"))

# custom samplers OPTIONNAL
conf$printSamplers(byType= TRUE)

conf$removeSamplers(c("mu0","mu1","p0","p1"))

conf$addSampler(c("mu0","mu1"), type ="RW_block")
conf$addSampler(c("p0","p1"), type ="RW_block")
conf$printSamplers(byType= TRUE)

conf$removeSamplers('id')
conf$addSampler('id', type = 'RW', scalarComponents = TRUE,
                control = list(adaptive = FALSE, reflective = TRUE, scale = datSCR$nsites/3)) ## DT


# Build and compile MCMC
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

# Run
# ART: 20-30 h 
t <- system.time(samples <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, samplesAsCodaMCMC = TRUE) ) ## DT: use runMCMC


save(samples,t, file ="SCR_res.Rdata")
