# Spatial capture recapture

# This code run SIM analyze for Bottlenose dolphins in the NW Mediterranean Sea
# It requires "SCR5.rdata", "DS2.rdata", "sea.grid.rdata" files
# Analyzes are run with NIMBLE, please visit https://r-nimble.org/ for information


#--- load packages  -----
library(nimble)
library(dplyr)
library(sf)

# ----function to load data  -----

# SCR
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
  eff.cov <- seff[samp.sites,] %>% 
    as_tibble() %>% 
    select(occ1, occ2,occ3,occ4, occ5, occ6, occ7, occ8) %>%
    log() %>% 
    as.matrix()
  
  eff.cov[eff.cov == - Inf] <-  0
  
  # make a binary cov with 1 if sites i is sampled during occ j, and 0 otherwise
  eff.ind <- eff.cov
  eff.ind[eff.ind!=0] <- 1
  
  
  dataSCR <- list(nsites = nsites, nind = nind, nocc = nocc, ntrap = ntrap, traploc = traploc,
                  habitat = habitat, seff = eff.cov, eff.ind = eff.ind, MAX = MAX, M = M, sites = sites, z=z,
                  y.scr = y.da2, id = id)
  
  return(dataSCR)
}

# DS 
dataDS2 <- function(){
  
  #load DS datasets
  load("DS/DS2.Rdata")
  load("sea.grid.rdata")
  
  g2 <- sea
  
  #---- DS datasets ----
  
  dim(sites) # "sites" contient les donnees sur les 1335 sites echantillonnés par SAMM avec les cov
  dim(indiv) # "indiv" contient les 129 detections de groupe de dauphins avec les cov
  
  # create variables for DS
  sitesDS <- sites[]
  nsitesDS <- nrow(sites) # nb de sites prospectes : 1335
  nobs <- nrow(indiv) # nb of detections : 129 # 
  nocc <-  2
  
  site <- indiv$DSid #   site id of every detection
  occ <- indiv$occ # sampling occasion of every detection
  
  # sampled sites by aerial surveys 
  gDS <- g2 %>% mutate(DSid = 1:nrow(g2))
  
  sampDS <- rep(NA, nrow(sites))
  for(i in 1:nrow(sites)){
    index <-  sites$objectid[i]
    sampDS[i] <-  which(g2$objectid == index)
  }
  
  
  # distance class format
  indiv$dec_angle
  indiv$d_perp[indiv$d_perp == max(indiv$d_perp)] <-  130.2
  d <- indiv$d_perp /max(indiv$d_perp) # distance perpendiculaire au transect. 0< d < 1
  B <- max(d)
  delta <- 0.1 # distance bin width, we make 10 bins
  midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
  dclass <- d %/% delta + 1 # convert distances to cat. distances
  nD <- length(midpt) # Number of distance intervals
  
  # group size format
  gpW <- sites$groupsizeW ## sitesDS$groupsize # taille du groupe detecte
  gpS <- sites$groupsizeS
  groupsize <- cbind(gpW,gpS) # Input groupsize as data
  groupsize[is.na(groupsize)] <- 0
  
  
  # covariates 
  
  # -- new sampling effort as % of sampling area
  sampAW <- (sitesDS$seffW * 1200) / st_area(sitesDS) # calculate % sampled area within each grid-cell
  sampAS <- (sitesDS$seffS * 1200) / st_area(sitesDS) # according to a 600m width bandstrip
  DSarea <- cbind(sampAW, sampAS)
  DSarea[DSarea> 1] <- 1
  #DSarea[DSarea == 0] <- 0.000001
  
  #
  seffDSW <- log(as.numeric(sitesDS$seffW)) # sapmpling effort covariate
  seffDSS <- log(as.numeric(sitesDS$seffS))
  seffDS <- cbind(seffDSW, seffDSS) # sapmpling effort covariate
  # 
  effind <- seffDS
  effind[effind>0] <- 1
  effind[effind<0] <- 0
  #habitat <- as.numeric(sites$bathy.sc) # depth as habitat covariate
  sea.stateDSW <- as.numeric(scale(sitesDS$seg_subj_numW)) # sea state as detection covariate
  sea.stateDSS <- as.numeric(scale(sitesDS$seg_subj_numS)) # sea state as detection covariate
  sea.stateDS <- cbind(sea.stateDSW, sea.stateDSS)
  #
  seffDS[effind ==0] <- 0
  sea.stateDS[effind ==0] <- 0
  # covariates
  habitat <- as.numeric(g2$bathy.sc)[sampDS] # cov
  
  
  # inits
  Nst <- groupsize + 1 
  Nst[is.na(Nst)] <- 0
  
  
  # save list
  dataDS <- list(nobs = nobs, B = B, nsitesDS = nsitesDS, midpt = midpt, ## DS
                 delta = delta, nD = nD, site = site,sampDS= sampDS,
                 seffDS = seffDS, 
                 sea.stateDS = sea.stateDS,
                 habitat = habitat,dclass= dclass,  groupsize = groupsize, Nst = Nst,
                 occ = occ, nocc = nocc, effind = effind,
                 DSarea = DSarea)
  
  return(dataDS)
}


# --- NIMBLE model ---- 

SIMcode <-  nimbleCode({
  # SCR
  # priors for the mu-density regression  
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  # priors for the sigma SCR regression
  sigma.scr ~ dunif(1,1000) 
  sig2 <- 2*sigma.scr*sigma.scr
  # priors for the sigma DS regression
  sigma.ds ~ dunif(0,10) 
  # priors for the pbar SCR regression
  p0 ~ dnorm(0,1) # Intercept
  p1 ~ dnorm(0,1) # slope
  # priors for the p0 DS regression
  alpha0 ~ dnorm(0,1)
  alpha1 ~ dnorm(0,1)
  
  # Lantent abundance via the IPP
  mu[1:nsites] <- exp(mu0 + mu1 * habitat[1:nsites])
  
  EN <-  sum(mu[1:nsites]) # expected number of individual for the study area
  
  psi <-  EN / M 
  
  probs[1:nsites] <- mu[1:nsites]/EN
  
  # SCR Observation model
  
  # logit regression of pbar
  logit(pbar[1:ntrap, 1:nocc]) <- p0 + p1 * seffSCR[1:ntrap,1:nocc]
  pbar_mask[1:ntrap, 1:nocc] <- pbar[1:ntrap,1:nocc]*eff.ind[1:ntrap,1:nocc]
  
  # data augmentation
  for( i in 1:M){ #
    
    z[i] ~ dbern(psi)  # is i alive or not ?
    
    id[i] ~ dunif(0.5, nsites+0.5)
    
    # location of AC 
    SX[i] <- sites[round(id[i]),1]
    SY[i] <- sites[round(id[i]),2]
    
    # zero's trick
    nll[i] <- - log(probs[round(id[i])])
    zeros[i] ~ dpois(nll[i])
    
    
    ## DT: likelihood using new dBernoulliVector4 distribution:
    y.scr[i, 1:nocc, 1:MAX] ~ dBernoulliVector4(
      MAX = MAX, ntrap = ntrap, nocc = nocc,
      sigma = sig2, pbar = pbar_mask[1:ntrap,1:nocc],
      eff.ind = eff.ind[1:ntrap, 1:nocc], SX = SX[i], SY = SY[i],
      traploc = traploc[1:ntrap,1:2], indicator = z[i])
  }# M
  
  # derived parameter for SCR 
  NtotSCR <- sum(z[1:M]) 
  
  # DS observation model
  
  # logit regression on p0 for DS
  logit(p0.ds[1:nsitesDS,1:noccDS]) <- alpha0 + alpha1 * sea.stateDS[1:nsitesDS,1:noccDS] 
  
  
  for(i in 1:nobs){
    dclass[i] ~ dcat(fc[site[i],1:nD,occ[i]]) # Part 1 of Hierachical Model : distance of observation of obs i
  }
  
  pi <- delta / B # Probability per interval
  
  for(s in 1:nsitesDS){
    
    lambda[s] <- mu[sampDS[s]] # Abundance estimated via the IPP at s
    
    for(t in 1:2){
      
      for( g in 1:nD){
        
        p[s,g,t] <- (p0.ds[s,t] * exp(- midpt[g] * midpt[g] / (2*sigma.ds*sigma.ds)))
        
        f[s,g,t] <- p[s,g,t] * pi
        
        fc[s,g,t] <- f[s,g,t] / pcap[s,t]
      }
      
      pcap[s,t] <- sum(f[s,1:nD,t]) # Pr(capture): sum of rectangular areas
      
      gpsize[s,t] ~ dbin(pcap[s,t], N_ds[s,t]) # Part 2 of HM
      
      N_ds[s,t] ~ dpois(lambda[s]*DSarea[s,t]*effindDS[s,t]) # account for sampling area
      
    }
    
  }
  
  # Derived parameters
  NtotDS <- sum(N_ds[1:nsitesDS,1:2])
  
})

##
# ---- Nimble custom distribution ----
dBernoulliVector4 <- nimbleFunction(
  run = function(
    x = double(2), MAX = double(0), ntrap = double(0), nocc = double(0),
    sigma = double(0), pbar = double(2), eff.ind = double(2),  SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0), log = integer(0, default = 0)) {
    if(indicator == 0) { if(sum(x) == 0) return(0) else return(-Inf) }
    lp <- 0
    xInd <- rep(1, nocc)
    for(s in 1:ntrap) {
      dist <- (SX-traploc[s,1])^2 + (SY-traploc[s,2])^2
      for(j in 1:nocc) {
        p <-  pbar[s,j]*exp(-dist/(sigma))
        if(x[j,xInd[j]] == s) {
          lp <- lp + log(p)
          if(xInd[j] < MAX) xInd[j] <- xInd[j] + 1
        } else {
          lp <- lp + log(1-p)
        }
      }
    }
    returnType(double())
    if(log) return(lp)
    return(exp(lp))
  })

rBernoulliVector4 <- nimbleFunction(
  run = function(
    n = integer(), MAX = double(0), ntrap = double(0), nocc = double(0),
    sigma = double(0), pbar = double(2), eff.ind = double(2), SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0)) {
    stop('should never call rBernoulliVector4')
    returnType(double(2))
    return(array(0, c(2, 2)))
  })

registerDistributions(list(
  dBernoulliVector4 = list(
    BUGSdist = 'dBernoulliVector4(MAX,ntrap,nocc,sigma,pbar,eff.ind,SX,SY,traploc,indicator)',
    types = c('value = double(2)', 'sigma = double(0)', 'pbar = double(2)','eff.ind = double(2)', 'traploc = double(2)'),
    mixedSizes = TRUE   ## DT: prevents a warning about size/dimension mismatch
  )))

# ---- load DATA ----
M = 5000

datSCR <- dataSCR(M = M)

datDS <- dataDS2()

# RUN ---- 

constants <- list(nobs = datDS$nobs, B = datDS$B, nsitesDS = datDS$nsitesDS, midpt = datDS$midpt, ## DS
                  delta = datDS$delta, nD = datDS$nD, site = datDS$site,sampDS= datDS$sampDS,
                  seffDS = datDS$seffDS, sea.stateDS = datDS$sea.stateDS, noccDS = datDS$nocc,
                  occ = datDS$occ, DSarea = datDS$DSarea, effindDS = datDS$effind,
                  habitat = datSCR$habitat, ## SCR
                  nsites = datSCR$nsites, nocc = datSCR$nocc, M = datSCR$M, ntrap = datSCR$ntrap,
                  traploc = datSCR$traploc/100000, seffSCR = datSCR$seff, eff.ind = datSCR$eff.ind, 
                  MAX = datSCR$MAX)

data <- list(sites= datSCR$sites/100000,dclass= datDS$dclass,  gpsize = datDS$groupsize, y.scr = datSCR$y.scr,
             zeros = rep(0,M))


inits <- list(alpha0=0, alpha1=0,mu0=0,mu1=0,p0=0,p1=0,
              sigma.ds = runif(1,0,10), sigma.scr = runif(1,1,1000),
              id = datSCR$id, z =datSCR$z ,  N_ds = datDS$Nst*datDS$effind)

# Build model 

Imodel <- nimbleModel(SIMcode, constants, data, inits)
Imodel$calculate() #-10267035
# Configure monitors
Iconf <- configureMCMC(Imodel)

Iconf$printMonitors() # see wht parameters are monitored
Iconf$resetMonitors()
Iconf$addMonitors('EN',"mu0","mu1",
                  "p0", "p1","sigma.scr", "sigma.ds",
                  "alpha0","alpha1","id") # add monitors 


# Custom samplers OPTIONNAL
# add Random Walk Block Samples for correlated parameters on sigma log-linear regression, and the IPP
Iconf$printSamplers(byType= TRUE)

Iconf$removeSamplers(target = c("mu0","mu1","alpha0","alpha1",
                                "p0","p1"))

Iconf$addSampler(target = c("mu0","mu1"), type ="RW_block")
Iconf$addSampler(target = c("p0","p1"), type ="RW_block")
Iconf$addSampler(target = c("alpha0","alpha1"), type ="RW_block")
Iconf$printSamplers(byType= TRUE)

Iconf$removeSamplers('id')
Iconf$addSampler('id', type = 'RW', scalarComponents = TRUE,
                 control = list(adaptive = FALSE, reflective = TRUE, scale = datSCR$nsites/3))


# Build and compile MCMC
IRmcmc <- buildMCMC(Iconf)
ICmodel <- compileNimble(Imodel)
ICmcmc <- compileNimble(IRmcmc, project = ICmodel)

# Run 
# ART : 20-30 h
t <- system.time(samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 10000, nchains = 3,thin = 10 , samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC

save(samples,t, file ="SIM_res.rdata")


