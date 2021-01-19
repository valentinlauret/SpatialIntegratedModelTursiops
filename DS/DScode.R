# Distance Sampling 

# This code run DS analyze for Bottlenose dolphins in the NW Mediterranean Sea
# It requires "DS2.rdata", and "sea.grid.rdata" files
# Analyzes are run with NIMBLE, please visit https://r-nimble.org/ for information

#--- load required packages  -----
library(nimble)
library(tidyverse)
library(sf)

# ---- function to load data  -----
load("sea.grid.rdata")
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

# --- load data --- 
datDS <- dataDS()

# ----- NIMBLE model -----
DScode <-  nimbleCode({
  
  # priors 
  # priors for the mu-density regression  
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  
  # priors for the sigma DS regression
  alpha0 ~ dnorm(0,1)
  alpha1 ~ dnorm(0,1)
  alpha2 ~ dnorm(0,1)
  
  # Abundance model via an inhomogeneous point process for all sites of the study area
  mu[1:nlam] <- exp(mu0 + mu1 * habitat[1:nlam])   
  
  # DS obervation model 
   for(i in 1:nobs){
     dclass[i] ~ dcat(fc[site[i],1:nD]) # Part 1 of Hierachical Model : distance of observation of obs i
   }
   
   for(s in 1:nsites){
     
     # Construct bin probabilities for nD multinomial bins
     for(g in 1:nD){ # midpt = mid-point of each bin
       
       log(p[s,g]) <- - midpt[g] * midpt[g] / (2*sigma[s]*sigma[s])
       
       pi[s,g] <- delta / B # Probability per interval
       
       f[s,g] <- p[s,g] * pi[s,g]
       
       fc[s,g] <- f[s,g] / pcap[s]
     }
     
     pcap[s] <- sum(f[s,1:nD]) # Pr(capture): sum of rectangular areas
     
     gpsize[s] ~ dbin(pcap[s], N[s]) # Part 2 of HM : observed abundance at s
     N[s] ~ dpois(lambda[s]) # Part 3 of HM : Latent abundance at s
     #log(lambda[s]) <- mu0 + mu1 * habitat[s] # Eq 2 Linear model abundance
     
     lambda[s] <- mu[sampDS[s]] # Abundance estimated via the IPP at s
     
     log(sigma[s])<- alpha0 + alpha1*seff[s]  + alpha2*sea.state[s] # log-linear model detection
   }
  
  # Derived parameters
   Ntot <- sum(N[1:nsites])
   EN <- sum(mu[1:nlam])
  
})

# ---- NIMBLE RUN ---- 

constants <- list(nobs = datDS$nobs, B = datDS$B, nsites = datDS$nsites,midpt = datDS$midpt,
                 seff = datDS$seff, sea.state = datDS$sea.state, 
                  delta = datDS$delta, nD = datDS$nD, site = datDS$site,
                  sampDS = datDS$sampDS, nlam = nrow(sea),  habitat = as.numeric(sea$bathy.sc))

data <- list(dclass =datDS$dclass, gpsize = datDS$groupsize)


inits <- list(alpha0=0, alpha1=0,alpha2= 0,mu0=0, mu1 = 0, N = datDS$Nst)

# Build model
  Rmodel <- nimbleModel(DScode, constants, data, inits)
  
  
# Configure monitors
  conf <- configureMCMC(Rmodel)

  conf$printMonitors() # see what parameters are monitored
  conf$addMonitors("EN", "Ntot") # add derived parameters to be monitored 

# Custom samplers OPTIONNAL
  # add Random Walk Block Samples for correlated parameters on sigma log-linear regression, and the IPP
  conf$printSamplers(byType= TRUE)

  conf$removeSamplers(target = c("mu0","mu1","alpha0","alpha1","alpha2"))
  
  conf$addSampler(target = c("mu0","mu1"), type ="RW_block")
  conf$addSampler(target = c("alpha0","alpha1","alpha2"), type ="RW_block")
  conf$printSamplers(byType= TRUE)

# Build and compile MCMC
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

# --- Run --- 
t <- system.time(samples <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC

save(samples,t, file ="DSres.Rdata")


