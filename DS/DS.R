# SIPM 1

# Combining Distance sampling & SCR data
# Distance sampling --> SAMM data
# SCR --> GDEGeM data

#--- load packages  -----
library(nimble)
library(tidyverse)
library(sf)

# ----function to load data  -----

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


# NIMBLE code 
# ---- BUGS model -----

DScode <-  nimbleCode({
  
  # priors 
  # priors for the mu-density regression  
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  # prior on sigma
  sigma ~ dunif(0,10)
  # p0 ~ dunif(0,1)
  alpha0 ~ dnorm(0,1) # Intercept of mu-density regression
  alpha1 ~ dnorm(0,1) # slope of unique covariate for now

  logit(p0[1:nsites,1:nocc]) <- alpha0 + alpha1 * sea.state[1:nsites,1:nocc] 
  # DS
  
  for(i in 1:nobs){
    dclass[i] ~ dcat(fc[site[i],1:nD,occ[i]]) # Part 1 of HM
  }
  
  pi <- delta / B # Probability per interval
  
  for(s in 1:nsites){
    
    for(t in 1:2){
      
      for( g in 1:nD){
        
        p[s,g,t] <- (p0[s,t] * exp(- midpt[g] * midpt[g] / (2*sigma*sigma)))
        
        f[s,g,t] <- p[s,g,t] * pi
        
        fc[s,g,t] <- f[s,g,t] / pcap[s,t]
      }
      
      pcap[s,t] <- sum(f[s,1:nD,t]) # Pr(capture): sum of rectangular areas
      
      gpsize[s,t] ~ dbin(pcap[s,t], N_ds[s,t]) # Part 2 of HM
      
      N_ds[s,t] ~ dpois(lambda[s]*DSarea[s,t]*effind[s,t]) # account for sampling area
      
    }
    
    log(lambda[s]) <- mu0 + mu1 * habitat[s] # Eq 2 Linear model abundance
    
  }
  
  ##  # Derived parameters
  EN <-  sum(lambda[1:nsites])     # expected number of individual for the study area
  
})


# ---- load DATA ----


datDS <- dataDS2()

# RUN ---- 

constants <- list(nobs = datDS$nobs, B = datDS$B, nsites = datDS$nsites,midpt = datDS$midpt,
                  habitat = datDS$habitat, seff = datDS$seff, sea.state = datDS$sea.state, 
                  delta = datDS$delta, nD = datDS$nD, site = datDS$site,
                  occ = datDS$occ, nocc = datDS$nocc, DSarea = datDS$DSarea, effind = datDS$effind)

data <- list(dclass =datDS$dclass, gpsize = datDS$groupsize)


inits <- list(mu0=0, mu1 = 0, alpha0 = 0, alpha1 = 0,   sigma = runif(1,0,1), N_ds = datDS$Nst*datDS$effind)


library(nimble)

# Build model


Rmodel <- nimbleModel(DScode, constants, data, inits)
Rmodel$initializeInfo()
Rmodel$calculate() #  -4792

# Configure monitors
conf <- configureMCMC(Rmodel)

conf$printMonitors() # see wht parameters are monitored
#conf$resetMonitors()
conf$addMonitors("EN") #,"mu0", "mu1","sigma") # add monitors 

# custom samplers OPTIONNAL
conf$printSamplers(byType= TRUE)

conf$removeSamplers(target = c("mu0","mu1"))
conf$removeSamplers(target = c("alpha0","alpha1"))
conf$addSampler(target = c("mu0","mu1"), type ="RW_block")
conf$addSampler(target = c("alpha0","alpha1"), type ="RW_block")
conf$printSamplers(byType= TRUE)

# Build and compile MCMC
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

# Run 
t <- system.time(samples <- runMCMC(Cmcmc, niter = 100000, nburnin = 10000, nchains = 3, samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC

library(mcmcplots)
denplot(samples)
save(samples, t, file = "DS_res.rdata")

