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
dataDS <- function(){
  
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
  
  site <- indiv$DSid #   site id of every detection
  
  # sampled sites by aerial surveys 
  gDS <- g2 %>% mutate(DSid = 1:nrow(g2))
  
  sampDS <- rep(NA, nrow(sites))
  for(i in 1:nrow(sites)){
    index <-  sites$objectid[i]
    sampDS[i] <-  which(g2$objectid == index)
  }
  
  
  # distance class format
  d <- indiv$d_perp /max(indiv$d_perp) # distance perpendiculaire au transect. 0< d < 1
  B <- max(d)
  delta <- 0.1 # distance bin width, we make 10 bins
  midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
  dclass <- d %/% delta + 1 # convert distances to cat. distances
  nD <- length(midpt) # Number of distance intervals
  
  # group size format
  gp <- sites$taille_grp ## sitesDS$groupsize # taille du groupe detecte
  groupsize <- gp # Input groupsize as data
  groupsize[is.na(groupsize)] <- 0
  
  
  # covariates 
  seffDS <- log(as.numeric(sitesDS$samm.effort)) # sapmpling effort covariate
  #habitat <- as.numeric(sites$bathy.sc) # depth as habitat covariate
  sea.stateDS <- as.numeric(scale(sitesDS$seg_subj_num)) # sea state as detection covariate
  
  # covariates
  habitat <- as.numeric(g2$bathy.sc)[sampDS] # cov
  
  
  # inits
  Nst <- groupsize + 1 
  Nst[is.na(Nst)] <- 0
  
  
  # save list
  dataDS <- list(nobs = nobs, B = B, nsitesDS = nsitesDS, midpt = midpt, ## DS
                 delta = delta, nD = nD, site = site,sampDS= sampDS,
                 seffDS = seffDS, sea.stateDS = sea.stateDS,
                 habitat = habitat,dclass= dclass,  groupsize = groupsize, Nst = Nst)
  
  return(dataDS)
}


# NIMBLE code 
# ---- BUGS model -----

DScode <-  nimbleCode({
  
  # priors 
  # priors for the mu-density regression  
  mu0 ~ dnorm(0,1) # Intercept of mu-density regression
  mu1 ~ dnorm(0,1) # slope of unique covariate for now
  # priors for the sigma DS regression
  alpha0 ~ dnorm(0,1)
  alpha1 ~ dnorm(0,1)
  alpha2 ~ dnorm(0,1)
  # prior on sigma
  sigma ~ dunif(0,5)
  
  # DS
   
   for(i in 1:nobs){
     dclass[i] ~ dcat(fc[site[i],1:nD]) # Part 1 of HM
   }
   
   for(s in 1:nsites){
     
     # Construct cell probabilities for nD multinomial cells
     for(g in 1:nD){ # midpt = mid-point of each cell
       
       p[s,g] <- p0[s] * exp(- midpt[g] * midpt[g] / (2*sigma*sigma))
       
       pi[s,g] <- delta / B # Probability per interval
       
       f[s,g] <- p[s,g] * pi[s,g]
       
       fc[s,g] <- f[s,g] / pcap[s]
     }
     
     pcap[s] <- sum(f[s,1:nD]) # Pr(capture): sum of rectangular areas
     gpsize[s] ~ dpois(pcap[s]*lambda[s])
     #gpsize[s] ~ dbin(pcap[s], N[s]) # Part 2 of HM
     N[s] ~ dpois(lambda[s]) # Part 3 of HM
     log(lambda[s]) <- mu0 + mu1 * habitat[s] # Eq 2 Linear model abundance
     logit(p0[s])<- alpha0 + alpha1*seff[s]  + alpha2*sea.state[s] # Linear model detection
   }
 ##  # Derived parameters
   #Ntot <- sum(N[1:nsites])
   EN <-  sum(lambda[1:nsites])     # expected number of individual for the study area
  # area <- nsites*1*2*B # Unit length == 1, half-width = B
  # D <- Ntotal/area
  
})


# ---- load DATA ----


datDS <- dataDS()

# RUN ---- 

constants <- list(nobs = datDS$nobs, B = datDS$B, nsites = datDS$nsites,midpt = datDS$midpt,
                  habitat = datDS$habitat, seff = datDS$seff, sea.state = datDS$sea.state, 
                  delta = datDS$delta, nD = datDS$nD, site = datDS$site)

data <- list(dclass =datDS$dclass, gpsize = datDS$groupsize)


inits <- list(alpha0=0, alpha1=0,alpha2= 0,mu0=0, mu1 = 0, sigma = runif(1,0,1), N = datDS$Nst)

# parralellize 
# 
# library(parallel)
# 
# this_cluster <- makeCluster(3)
# 
# set.seed(32)
# 
# DS_para <- function(seed, code, constants, data, inits){
  
  library(nimble)


# nimbleMCMC
out <- nimbleMCMC(code = DScode,
           data = data,
           constants = constants,
           inits = inits,
           monitors = c("EN","mu0", "mu1","sigma","alpha0","alpha1","alpha2"),
           niter = 1000,
           nburnin = 200,
           nchains = 2)
# Build model
#tSIPM1 <- proc.time()

Rmodel <- nimbleModel(DScode, constants, data, inits)
Rmodel$initializeInfo()
Rmodel$calculate() #  -4792

# Configure monitors
conf <- configureMCMC(Rmodel)

conf$printMonitors() # see wht parameters are monitored
conf$resetMonitors()
conf$addMonitors("EN","mu0", "mu1","sigma","alpha0","alpha1","alpha2") # add monitors 

# custom samplers OPTIONNAL
conf$printSamplers(byType= TRUE)

conf$removeSamplers(target = c("mu0","mu1","alpha0","alpha1","alpha2"))
#conf$removeSamplers(target = c("mu0","mu1","beta0","beta1"))
conf$addSampler(target = c("mu0","mu1"), type ="RW_block")
conf$addSampler(target = c("alpha0","alpha1","alpha2"), type ="RW_block")
conf$printSamplers(byType= TRUE)

# Build and compile MCMC
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

# Run 
t <- system.time(samples <- runMCMC(Cmcmc, niter = 50000, nburnin = 5000, nchains = 3, samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC

save(samples, t, file = "DS_sig0res2.rdata")

