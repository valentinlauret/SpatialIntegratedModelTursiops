# test prior sigma et p 
denplot(samples[,c("mu0","EN","mu1","sigma","p0","p1")])

# try to see what scale for sigma.ds and sigma.scr

# DS
sigma.ds <- 1.03
p.ds <- exp(- datDS$midpt / (2*sigma.ds^2))
p.ds

ggplot() + 
  geom_line(aes(x = datDS$midpt, y = p.ds))

# SCR
traploc <-  datSCR$traploc
sites <-  datSCR$sites
head(sites)

sigma.scr <- 100000
p.scr <- dist2.scr <- rep(NA,1000)

for(i in 1:1000){
trap <- traploc[sample(nrow(traploc),1),]

id <-  runif(1, 0.5, nrow(sites)+0.5)

# location of AC 
SX <- sites[round(id),1]
SY <- sites[round(id),2]

dist2.scr[i] <- (SX-trap[1])^2 + (SY-trap[2])^2

p.scr[i] <- exp(-dist2.scr[i]/(2*sigma.scr^2))
}

ggplot() + 
  geom_line(aes(x = sqrt(dist2.scr)/1000, y = p.scr))

denplot(sqrt(dist.scr)/1000)
