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
dataDS <- function(){
  
  #load DS datasets
  load("DS/DS2.Rdata")
  load("sea.grid.rdata")
  
  g2 <- sea
  
  #---- DS datasets ----
  
  dim(sites) # "sites" includes the coordinates of 1335 sampled sites with associated covariates
  dim(indiv) # "indiv" includes the locations of 129 dolphins detection with associated covariates
  
  # create variables for DS
  sitesDS <- sites[]
  nsitesDS <- nrow(sites) # nb of sampled sites : 1335
  nobs <- nrow(indiv) # nb of detections : 129  
  
  site <- indiv$DSid #   site id of every detection
  
  gDS <- g2 %>% mutate(DSid = 1:nrow(g2))  # sampled sites by aerial surveys 
  
  sampDS <- rep(NA, nrow(sites))
  for(i in 1:nrow(sites)){
    index <-  sites$objectid[i]
    sampDS[i] <-  which(g2$objectid == index)
  }
  
  # format detection data with 10 Bins
  d <- indiv$d_perp / max(indiv$d_perp) # perpendicular distance to the transect
  B <- max(d) # max distance
  delta <- 0.1 # distance bin width: 0.1, we created 10 bins
  midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
  dclass <- d %/% delta + 1 # convert distances to cat. distances
  nD <- length(midpt) # Number of distance intervals
  
  # group size format
  gp <- sites$taille_grp  # group size detected
  groupsize <- gp # Input groupsize as data
  groupsize[is.na(groupsize)] <- 0
  
  
  # sampling effort covariates 
  seffDS <- log(as.numeric(sitesDS$samm.effort)) # sampling effort covariate
  sea.stateDS <- as.numeric(scale(sitesDS$seg_subj_num)) # sea state as detection covariate
  
  # environmental covariate
  habitat <- as.numeric(g2$bathy.sc)[sampDS]
  
  # inits
  Nst <- groupsize + 1 
  Nst[is.na(Nst)] <- 0
  
  # save data
  dataDS <- list(nobs = nobs, B = B, nsitesDS = nsitesDS, midpt = midpt,
                 delta = delta, nD = nD, site = site,sampDS= sampDS,
                 seffDS = seffDS, sea.stateDS = sea.stateDS,
                 habitat = habitat,dclass= dclass,  groupsize = groupsize, Nst = Nst)
  
  return(dataDS)
}


# --- NIMBLE model ---- 

SIMcode <-  nimbleCode({
  # SCR
  # priors for the mu-density regression  
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  # priors for the sigma SCR regression
  beta0 ~ dnorm(0,1) # Intercept of mu-density regression
  beta1 ~ dnorm(0,1)
  # priors for the pbar SCR regression
  p0 ~ dnorm(0,1) # Intercept
  p1 ~ dnorm(0,1) # slope
  # priors for the sigma DS regression
  alpha0 ~ dnorm(0,1)
  alpha1 ~ dnorm(0,1)
  alpha2 ~ dnorm(0,1)
  
  # Lantent abundance via the IPP
  mu[1:nsites] <- exp(mu0 + mu1 * habitat[1:nsites])
 
  EN <-  sum(mu[1:nsites]) # expected number of individual for the study area
  
  psi <-  EN / M 
  
  # SCR Observation model
  
  # log-linear regression on sigma 
  log(sigmaSCR[1:ntrap,1:nocc]) <- beta0 + beta1 * seffSCR[1:ntrap,1:nocc]
  
  # logit regression of pbar
  logit(pbar[1:ntrap, 1:nocc]) <- p0 + p1 * seffSCR[1:ntrap,1:nocc] 
  
  # cell prob
  # cellprob[1:nsites] <- mu[1:nsites]/(sum(mu[1:nsites]))
  
  # data augmentation
  for( i in 1:M){ #
    
    z[i] ~ dbern(psi)  # is i alive or not ?
    
    # induce the uniform distribution of activity centers over the discrete set 1:nsites
    id[i] ~  dcat(mu[1:nsites]) # dunif(0.5, nsites+0.5) #
    
    # location of AC
    SX[i] <- sites[round(id[i]),1]
    SY[i] <- sites[round(id[i]),2]
    
   
    ## DT: likelihood using new dBernoulliVector4 distribution:
    y.scr[i, 1:nocc, 1:MAX] ~ dBernoulliVector4(
      MAX = MAX, ntrap = ntrap, nocc = nocc,
      sigma = sigmaSCR[1:ntrap,1:nocc], pbar = pbar[1:ntrap,1:nocc],
      eff.ind = eff.ind[1:ntrap, 1:nocc], SX = SX[i], SY = SY[i],
      traploc = traploc[1:ntrap,1:2], indicator = z[i])
  }# M
  
  # derived parameter for SCR 
  NtotSCR <- sum(z[1:M]) 

  # DS observation model
   
   for(i in 1:nobs){
     dclass[i] ~ dcat(fc[site[i],1:nD]) # Part 1 of Hierachical Model : distance of observation of obs i
   }
   
   for(s in 1:nsitesDS){
     
     # Construct bin probabilities for nD multinomial bins
     for(g in 1:nD){ # midpt = mid-point of each bin
       
       log(p[s,g]) <- - midpt[g] * midpt[g] / (2*sigmaDS[s]*sigmaDS[s])
       
       pi[s,g] <- delta / B # Probability per interval
       
       f[s,g] <- p[s,g] * pi[s,g]
       
       fc[s,g] <- f[s,g] / pcap[s]
     }
     
     pcap[s] <- sum(f[s,1:nD]) # Pr(capture): sum of rectangular areas
     gpsize[s] ~ dbin(pcap[s], N[s]) # Part 2 of HM : observed abundance at s
     N[s] ~ dpois(lambda[s]) # Part 3 of HM : latent abundance at s 
     lambda[s] <- mu[sampDS[s]] # Abundance estimated via the IPP at s
     
     log(sigmaDS[s])<- alpha0 + alpha1*seffDS[s]  + alpha2*sea.stateDS[s] # log-linear model detection
   }
  
   # Derived parameters
   NtotDS <- sum(N[1:nsitesDS])
  
})

# ---- Nimble custom distribution ----
dBernoulliVector4 <- nimbleFunction(
  run = function(
    x = double(2), MAX = double(0), ntrap = double(0), nocc = double(0),
    sigma = double(2), pbar = double(2), eff.ind = double(2),  SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0), log = integer(0, default = 0)) {
    if(indicator == 0) { if(sum(x) == 0) return(0) else return(-Inf) }
    lp <- 0
    xInd <- rep(1, nocc)
    for(s in 1:ntrap) {
      dist <- (SX-traploc[s,1])^2 + (SY-traploc[s,2])^2
      for(j in 1:nocc) {
        p <-  pbar[s,j]*max(min(exp(-dist/(2*sigma[s,j]^2)),0.999),0.001)*eff.ind[s,j]
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
    sigma = double(2), pbar = double(2), eff.ind = double(2), SX = double(0), SY = double(0),
    traploc = double(2), indicator = double(0)) {
    stop('should never call rBernoulliVector4')
    returnType(double(2))
    return(array(0, c(2, 2)))
  })

registerDistributions(list(
  dBernoulliVector4 = list(
    BUGSdist = 'dBernoulliVector4(MAX,ntrap,nocc,sigma,pbar,eff.ind,SX,SY,traploc,indicator)',
    types = c('value = double(2)', 'sigma = double(2)', 'pbar = double(2)','eff.ind = double(2)', 'traploc = double(2)'),
    mixedSizes = TRUE   ## DT: prevents a warning about size/dimension mismatch
  )))

# ---- load DATA ----

datSCR <- dataSCR(M = 5000)

datDS <- dataDS()

# RUN ---- 

constants <- list(nobs = datDS$nobs, B = datDS$B, nsitesDS = datDS$nsitesDS, midpt = datDS$midpt, ## DS
                  delta = datDS$delta, nD = datDS$nD, site = datDS$site,sampDS= datDS$sampDS,
                  seffDS = datDS$seffDS, sea.stateDS = datDS$sea.stateDS,
                  habitat = datSCR$habitat, ## SCR
                  nsites = datSCR$nsites, nocc = datSCR$nocc, M = datSCR$M, ntrap = datSCR$ntrap,
                  traploc = datSCR$traploc, seffSCR = datSCR$seff, eff.ind = datSCR$eff.ind, MAX = datSCR$MAX)

data <- list(sites= datSCR$sites,dclass= datDS$dclass,  gpsize = datDS$groupsize, y.scr = datSCR$y.scr)


inits <- list(alpha0=0, alpha1=0,alpha2= 0,beta0=0, beta1 = 0,mu0=0,mu1=0,p0=0,p1=0,
              id = datSCR$id, z =datSCR$z , N =datDS$Nst)
  
# Build model 

  Rmodel <- nimbleModel(SIMcode, constants, data, inits)
 Rmodel$calculate() #_121041
  # Configure monitors
  conf <- configureMCMC(Rmodel)
  
  conf$printMonitors() # see wht parameters are monitored
  conf$resetMonitors()
  conf$addMonitors('EN','NtotSCR',"NtotDS","psi","mu0","mu1",
                   "beta0","beta1","p0", "p1",
                   "alpha0","alpha1","alpha2") # add monitors 
   
  
  # Custom samplers OPTIONNAL
  # add Random Walk Block Samples for correlated parameters on sigma log-linear regression, and the IPP
  conf$printSamplers(byType= TRUE)
  
  conf$removeSamplers(target = c("mu0","mu1","beta0","beta1","alpha0","alpha1","alpha2",
                                 "p0","p1"))

  conf$addSampler(target = c("mu0","mu1"), type ="RW_block")
  conf$addSampler(target = c("beta0","beta1"), type ="RW_block")
  conf$addSampler(target = c("alpha0","alpha1","alpha2"), type ="RW_block")
  conf$printSamplers(byType= TRUE)
  
  conf$removeSamplers('id')
  conf$addSampler('id', type = 'RW', scalarComponents = TRUE,
                 control = list(adaptive = FALSE, reflective = TRUE, scale = datSCR$nsites/3))
  
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
  
  # Run 
  # ART : 2.5 days
  t <- system.time(samples <- runMCMC(Cmcmc, niter = 100000, nburnin = 35000, nchains = 3, samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC
  # ART: 500s
  t <- system.time(samples <- runMCMC(Cmcmc, niter = 1000, nburnin = 350, nchains = 1, samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC
  
save(samples,t, file ="SIMres.rdata")

