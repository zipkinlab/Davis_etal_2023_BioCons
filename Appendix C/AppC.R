########################################################################################
# Integrated population model (IPM) for St. Clair Flats Black Terns, 2013 - 2022
# Kayla Davis, Sarah Saunders, Stephanie Beilke, Erin Ford, Jenni Fuller, Ava Landgraf, and Elise Zipkin

# Adapted from original scripts by Michael Schaub & Marc Kéry (2021)
# Modified by K. Davis, 2022

# Code to recreate the results from Appendix D: BPVA with small, positive trend

########################################################################################



##################################################################
# BPVA with small positive trend
# Appendix C

# requires data (SCF-IPM_github.xlsx, SCF_AHY13to22_update.txt, SCF_HY13to22_update.txt, YearCovs.csv) from main text analysis to fit the IPM
# requires Hist_NAO_Means.csv to complete the bootstrap method
##################################################################


# setup ------------------------------------------------------------------------
#load libraries
library(openxlsx)
library(tidyverse)
library(lubridate)
library(coda)
library(ggpubr)
library(MCMCvis)
library(jagsUI)
library(scales)
library(plotrix)

#load function to create a m-array based on capture-recapture data (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  
  #Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  
  #Calculate the number of individuals that are never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}



# bootstrap method ------------------------------------------------------------------------

###################################################
## Bootstrap method to obtain temporal variation
## in NAO values over time
##################################################

## Read in data
nao_all <- read.csv("Hist_NAO_Means.csv", stringsAsFactors = FALSE)
st.mean <- mean(nao_all$NAO_Jan_Jun_Mean[114:123])
lt.mean <- mean(nao_all$NAO_Jan_Jun_Mean[24:123])

## Bootstrap method for pulling 10 random values of historical NAO
## to calculate SD [temporal variability of NAO]
library(boot)

sdstats <- function(data, indices){
  dat <- data[indices,] #allows boot to select sample
  s1 <- sample(dat[,2], size=10, replace=FALSE)
  sd.s <- sd(s1)
  return(sd.s)
}

# bootstrapping wtih 1000 reps
set.seed(1234)
results <- boot(data=nao_all, statistic=sdstats,
                R=100000)

# view results
results
mnao.sd <- mean(results$t) #mean of bootstrapped samples: 0.58
sdnao.sd <- sd(results$t) #sd of bootstrapped samples: 0.13



# the IPM data ---------------------------------------------------------------------

########################################################################
# Capture-recapture data: m-array of juveniles (HY) and adults (AHY)
########################################################################

#First read in capture histories for birds marked as chicks during 2013-2022
CH.J <- read.table("SCF_HY13to22_update.txt")

#convert to matrix
CH.J <- data.matrix(CH.J)

#read in capture histories for birds marked as adults during 2013-2022
CH.A <- read.table("SCF_AHY13to22_update.txt")

#convert to matrix
CH.A <- data.matrix(CH.A)

#create two m-arrays, one for juveniles and one for adults
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,] # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,] # Juvenile CH never recaptured

#Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

#Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)

#Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

#Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

#The last column should show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

#Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray

#outputs: CH.A.marray and CH.J.marray
#convert outputs to names of m-arrays used in models

#delete last 2 rows of juvenile m-array (can't recap birds released in last 2 occassions)
marray.j <- CH.J.marray[-c(7:8),]
marray.a <- CH.A.marray

# #Create NA marrays
# marray.j <- matrix(0, nrow = 7, ncol = 10)
# marray.a <- matrix(0, nrow = 9, ncol = 10)


########################################################################
# Population count data
########################################################################

dat1 <- read.xlsx("SCF-IPM.xlsx")

# breeding pairs and productivity
year <- dat1$Year
nyears <- length(year)
y <- dat1$MinPairs # the number of breeding pairs per year
j1 <- dat1$MinFledges # the number of fledglings recorded flying 
j2 <- dat1$Nanotag # the number of fledglings recorded via nanotag
R <- y # number of pairs/broods monitored
R_tag <- dat1$NumberTagged # number of fledglings tagged

########################################################################
# Covariates
########################################################################

# read in covariate data (top supported model only)
yeareffects = read.csv("YearCovs.csv", header=T, sep=',', na.strings=T)

# adult survival covs
nao_hur = yeareffects$nao_jan.jun

# make into vector for jags
nao_hur = as.numeric(as.vector(nao_hur[1:9]))

#calculate mean nao
mnao <- mean(nao_hur)
sdnao <- sd(nao_hur)

# standardize effects
znao_hur <- as.vector(scale(nao_hur))


# the Population Viability Analysis-------------------------------------------------------------

###################################################
# Combined mgmt scenario: 
# increase adult survival 5%
# increase juv survival 10%
# double fecundity
##################################################

sink("blte_nao_bpva_all_trend")
cat("
model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(5, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois(0.5 * f1[1] * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)          # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)
    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation

    
    # Priors for beta coefficients
    beta.phia1 ~ dnorm(0, 0.01)


    # Precision of standard deviations of temporal variability
    sig.phij ~ dexp(1)
    tau.phij <- pow(sig.phij, -2) 
    sig.phia ~ dexp(1)
    tau.phia <- pow(sig.phia, -2) 
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.13 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 0.58 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    
    mu.nao[9] <- st.mean                    #use long-term (lt.mean) or short-term (st.mean) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] + 0.004        #use small, positive trend from 100-yr regression 
    nao.pre[t] ~ dnorm(mu.nao[t],tau.nao)      
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao 
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    }
    
    # Future: use priors for NAO on phia and increase adult survival ~5% (0.85 --> 0.90)
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- log(0.9/(1-0.9)) + beta.phia1*nao.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }
    
    
    #Past fecundity
    for (t in 1:(nyears)){
    log(f1[t]) <- l.mfec1
    }
    
    #Future fecundity
    for (t in (nyears+1):(nyears+K)){
    log(f1[t]) <- l.mfec1 + log(2)
    }
    
    #Past juv survival
    for (t in 1:(nyears-3)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]
    }

    #Future juv survival (increase mean juv survival by ~10% --> .08 to .18)
    for (t in (nyears-2):(nyears-3+K)){
    logit(phij[t]) <- log(0.18/(1-0.18)) + epsilon.phij[t]
    }

    for (t in 1:(nyears-1+K)){    # Same for past and future
    logit(p[t]) <- l.p + epsilon.res[t]
    }
    
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts

    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:(nyears+K)){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:(nyears+K)){
    mean1[t] <- 0.5 * f1[t] * Ntot[t]    # all adults can breed
    N0[t] ~ dpois(mean1[t])                 # Fledglings 
    Nad[t] ~ dbin(phia[t-1], Ntot[t-1])
    }
     
    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(Ntot[t], tau.obs)          # all adults are counted
    }
    
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    

    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears) {
    j1[t] ~ dpois(rho1[t])    # number young observed as fledged
    rho1[t] <- R[t] * f1[t]    # number of pairs and fecundity
    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
K<-10
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, 
                  marray.a = marray.a, r.j = rowSums(marray.j), 
                  r.a = rowSums(marray.a), y = y, j1 = j1, 
                  R = R, mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao = mnao, sdnao=sdnao, K=K, st.mean = st.mean, lt.mean=lt.mean)

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1", "mlam.tot",
                "mlam.hist", "l.mphij", "l.mphia","l.mfec1", "p",
                "sig.phij", "sig.phia", "sig.fec","sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1",
                "nao.fut", "nao.pre") 

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

# Call JAGS from R
nao.s.trend.all <- jags(jags.data, inits, parameters, "blte_nao_bpva_all_trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)


###################################################
# Adult survival mgmt scenario: 
# increase adult survival 5%
##################################################

sink("blte_nao_bpva_phia_trend")
cat("
model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(5, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois(0.5 * f1[1] * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)          # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)
    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation

    
    # Priors for beta coefficients
    beta.phia1 ~ dnorm(0, 0.01)


    # Precision of standard deviations of temporal variability
    sig.phij ~ dexp(1)
    tau.phij <- pow(sig.phij, -2) 
    sig.phia ~ dexp(1)
    tau.phia <- pow(sig.phia, -2) 
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.13 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 0.58 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[9] <- st.mean                    #use long-term (lt.mean) or short-term (st.mean) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] + 0.004        #use small, positive trend from 100-yr regression 
    nao.pre[t] ~ dnorm(mu.nao[t],tau.nao)      
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao 
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    }
    
    # Future: use priors for NAO on phia and increase adult survival ~5% (0.85 --> 0.90)
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- log(0.9/(1-0.9)) + beta.phia1*nao.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }
    
    # Fecundity same for past and future
    for (t in 1:(nyears+K)){
    log(f1[t]) <- l.mfec1
    }
    
    #juv survival same for past and future
    for (t in 1:(nyears-3+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]
    }

    
    for (t in 1:(nyears-1+K)){    # Same for past and future
    logit(p[t]) <- l.p + epsilon.res[t]
    }
    
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts


    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:(nyears+K)){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:(nyears+K)){
    mean1[t] <- 0.5 * f1[t] * Ntot[t]    # all adults can breed
    N0[t] ~ dpois(mean1[t])                 # Fledglings 
    Nad[t] ~ dbin(phia[t-1], Ntot[t-1])
    }
     

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(Ntot[t], tau.obs)          # all adults are counted
    }
    
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    

    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears) {
    j1[t] ~ dpois(rho1[t])    # number young observed as fledged
    rho1[t] <- R[t] * f1[t]    # number of pairs and fecundity
    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
###################################################################
# Bundle data
K<-10
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, 
                  marray.a = marray.a, r.j = rowSums(marray.j), 
                  r.a = rowSums(marray.a), y = y, j1 = j1, 
                  R = R, mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao = mnao, sdnao=sdnao, K=K, st.mean = st.mean, lt.mean=lt.mean)

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1", "mlam.tot",
                "mlam.hist", "l.mphij", "l.mphia","l.mfec1", "p",
                "sig.phij", "sig.phia", "sig.fec","sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1",
                "nao.fut", "nao.pre")

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

nao.s.trend.phia <- jags(jags.data, inits, parameters, "blte_nao_bpva_phia_trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# Juvenile survival mgmt scenario: 
# increase juv survival 10%
##################################################

sink("blte_nao_bpva_phij_trend")
cat("
model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(5, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois(0.5 * f1[1] * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)          # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)
    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation

    
    # Priors for beta coefficients
    beta.phia1 ~ dnorm(0, 0.01)


    # Precision of standard deviations of temporal variability
    sig.phij ~ dexp(1)
    tau.phij <- pow(sig.phij, -2) 
    sig.phia ~ dexp(1)
    tau.phia <- pow(sig.phia, -2) 
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.13 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 0.58 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[9] <- st.mean                    #use long-term (lt.mean) or short-term (st.mean) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] + 0.004        #use small, positive trend from 100-yr regression 
    nao.pre[t] ~ dnorm(mu.nao[t],tau.nao)      
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao 
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past  adult survival
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    }
    
    # Future: use priors for NAO on phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*nao.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }

    
    # Fecundity same for past and future
    for (t in 1:(nyears+K)){
    log(f1[t]) <- l.mfec1
    }
    
    
    #Past juv survival
    for (t in 1:(nyears-3)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]
    }

    #Future juv survival (increase mean juv survival by ~10% --> .08 to .18)
    for (t in (nyears-2):(nyears-3+K)){
    logit(phij[t]) <- log(0.18/(1-0.18)) + epsilon.phij[t]
    }

    
    for (t in 1:(nyears-1+K)){    # Same for past and future
    logit(p[t]) <- l.p + epsilon.res[t]
    }
    
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts

  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:(nyears+K)){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:(nyears+K)){
    mean1[t] <- 0.5 * f1[t] * Ntot[t]    # all adults can breed
    N0[t] ~ dpois(mean1[t])                 # Fledglings 
    Nad[t] ~ dbin(phia[t-1], Ntot[t-1])
    }
     

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(Ntot[t], tau.obs)          # all adults are counted
    }
    
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    

    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears) {
    j1[t] ~ dpois(rho1[t])    # number young observed as fledged
    rho1[t] <- R[t] * f1[t]    # number of pairs and fecundity
    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
###################################################################
# Bundle data
K<-10
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, 
                  marray.a = marray.a, r.j = rowSums(marray.j), 
                  r.a = rowSums(marray.a), y = y, j1 = j1, 
                  R = R, mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao = mnao, sdnao=sdnao, K=K, st.mean = st.mean, lt.mean=lt.mean)

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1", "mlam.tot",
                "mlam.hist", "l.mphij", "l.mphia","l.mfec1", "p",
                "sig.phij", "sig.phia", "sig.fec","sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1",
                "nao.fut", "nao.pre")

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

nao.s.trend.phij <- jags(jags.data, inits, parameters, "blte_nao_bpva_phij_trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# Fecundity mgmt scenario: 
# double annual fecundity
##################################################

sink("blte_nao_bpva_fec_trend")
cat("
model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(5, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois(0.5 * f1[1] * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)          # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)
    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation

    
    # Priors for beta coefficients
    beta.phia1 ~ dnorm(0, 0.01)


    # Precision of standard deviations of temporal variability
    sig.phij ~ dexp(1)
    tau.phij <- pow(sig.phij, -2) 
    sig.phia ~ dexp(1)
    tau.phia <- pow(sig.phia, -2) 
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.13 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 0.58 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[9] <- st.mean                    #use long-term (lt.mean) or short-term (st.mean) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] + 0.004        #use small, positive trend from 100-yr regression 
    nao.pre[t] ~ dnorm(mu.nao[t],tau.nao)      
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao 
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past  adult survival
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    }
    
    # Future: use priors for NAO on phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*nao.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }

    
    #Past fecundity
    for (t in 1:(nyears)){
    log(f1[t]) <- l.mfec1
    }
    
    #Future fecundity
    for (t in (nyears+1):(nyears+K)){
    log(f1[t]) <- l.mfec1 + log(2)
    }
    
    #juv survival same for past and future
    for (t in 1:(nyears-3+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]
    }

    
    for (t in 1:(nyears-1+K)){    # Same for past and future
    logit(p[t]) <- l.p + epsilon.res[t]
    }
    
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts

  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:(nyears+K)){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:(nyears+K)){
    mean1[t] <- 0.5 * f1[t] * Ntot[t]    # all adults can breed
    N0[t] ~ dpois(mean1[t])                 # Fledglings 
    Nad[t] ~ dbin(phia[t-1], Ntot[t-1])
    }
     

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(Ntot[t], tau.obs)          # all adults are counted
    }
    
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    

    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears) {
    j1[t] ~ dpois(rho1[t])    # number young observed as fledged
    rho1[t] <- R[t] * f1[t]    # number of pairs and fecundity
    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
###################################################################
# Bundle data
K<-10
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, 
                  marray.a = marray.a, r.j = rowSums(marray.j), 
                  r.a = rowSums(marray.a), y = y, j1 = j1, 
                  R = R, mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao = mnao, sdnao=sdnao, K=K, st.mean = st.mean, lt.mean=lt.mean)

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1", "mlam.tot",
                "mlam.hist", "l.mphij", "l.mphia","l.mfec1", "p",
                "sig.phij", "sig.phia", "sig.fec","sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1",
                "nao.fut", "nao.pre")

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

nao.s.trend.fec <- jags(jags.data, inits, parameters, "blte_nao_bpva_fec_trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# No additional mgmt scenario: 
# Same for past and future
##################################################

sink("blte_nao_bpva_null_trend")
cat("
model {

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(5, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois(0.5 * f1[1] * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)          # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)
    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation

    
    # Priors for beta coefficients
    beta.phia1 ~ dnorm(0, 0.01)


    # Precision of standard deviations of temporal variability
    sig.phij ~ dexp(1)
    tau.phij <- pow(sig.phij, -2) 
    sig.phia ~ dexp(1)
    tau.phia <- pow(sig.phia, -2) 
    sig.res ~ dunif(0, 10)
    tau.res <- pow(sig.res, -2)
    
    sig.obs ~ dunif(0.5, 50)   # residual variance
    tau.obs <- pow(sig.obs, -2)
    
    sig.sig <- sdnao.sd  #sdnao.sd is 0.13 from bootstrap results 
    tau.sig <- pow(sig.sig, -2)
    sig.nao ~ dnorm(mnao.sd, tau.sig) #mnao.sd is 0.58 from bootstrap results
    tau.nao <- pow(sig.nao, -2)
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    
    for (t in 1:(nyears-3+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    mu.nao[9] <- st.mean                    #use long-term (lt.mean) or short-term (st.mean) mean here
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    mu.nao[t] <- mu.nao[t-1] + 0.004        #use small, positive trend from 100-yr regression 
    nao.pre[t] ~ dnorm(mu.nao[t],tau.nao)      
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao 
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    # Past  adult survival
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    }
    
    # Future: use priors for NAO on phia
    for (t in nyears:(nyears-1+K)){
    logit(phia[t]) <- l.mphia + beta.phia1*nao.fut[t] + epsilon.phia[t]              # Adult apparent survival (3+ y old survival)
    }

    
    # Fecundity same for past and future
    for (t in 1:(nyears+K)){
    log(f1[t]) <- l.mfec1
    }

    
    #juv survival same for past and future
    for (t in 1:(nyears-3+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]
    }

    
    for (t in 1:(nyears-1+K)){    # Same for past and future
    logit(p[t]) <- l.p + epsilon.res[t]
    }
    
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts

  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:(nyears+K)){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:(nyears+K)){
    mean1[t] <- 0.5 * f1[t] * Ntot[t]    # all adults can breed
    N0[t] ~ dpois(mean1[t])                 # Fledglings 
    Nad[t] ~ dbin(phia[t-1], Ntot[t-1])
    }
     

    # 4.1.2 Observation process
    for (t in 1:nyears){      
    y[t] ~ dnorm(Ntot[t], tau.obs)          # all adults are counted
    }
    
    
    # 4.2 Likelihood for capture-recapture data: CJS model
    # Multinomial likelihood
    for (t in 1:(nyears-3)){                  
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t]) 
    }
    
    for (t in 1:(nyears-1)){
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])  
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-3)){    
    # Main diagonal
    pr.j[t,(t+2)] <- phij[t]*p[t+2]
    # Above main diagonal
    for (j in (t+3):(nyears-1)){   
    pr.j[t,j] <- phij[t]*prod(phia[(t+3):j])*prod(q[(t+2):(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t+1)){  
    pr.j[t,j] <- 0   
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults  
    for (t in 1:(nyears-1)){      
    q[t] <- 1-p[t]     
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    

    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:nyears) {
    j1[t] ~ dpois(rho1[t])    # number young observed as fledged
    rho1[t] <- R[t] * f1[t]    # number of pairs and fecundity
    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
###################################################################
# Bundle data
K<-10
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, 
                  marray.a = marray.a, r.j = rowSums(marray.j), 
                  r.a = rowSums(marray.a), y = y, j1 = j1, 
                  R = R, mnao.sd = mnao.sd, sdnao.sd = sdnao.sd,
                  mnao = mnao, sdnao=sdnao, K=K, st.mean = st.mean, lt.mean=lt.mean)

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1", "mlam.tot",
                "mlam.hist", "l.mphij", "l.mphia","l.mfec1", "p",
                "sig.phij", "sig.phia", "sig.fec","sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1",
                "nao.fut", "nao.pre")

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

nao.s.trend <- jags(jags.data, inits, parameters, "blte_nao_bpva_null_trend", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------
# ~~~~ save output for use later ~~~~
 # save(nao.trend, file="nao.trend.Rdata", compress = "xz", compression_level = 9)
 # save(nao.trend.fec, file="nao.trend.fec.Rdata", compress = "xz", compression_level = 9)
 # save(nao.trend.phij, file="nao.trend.phij.Rdata", compress = "xz", compression_level = 9)
 # save(nao.trend.phia, file="nao.trend.phia.Rdata", compress = "xz", compression_level = 9)
 # save(nao.trend.all, file="nao.trend.all.Rdata", compress = "xz", compression_level = 9)
 # save(nao.s.trend, file="nao.s.trend.Rdata", compress = "xz", compression_level = 9)
 # save(nao.s.trend.fec, file="nao.s.trend.fec.Rdata", compress = "xz", compression_level = 9)
 # save(nao.s.trend.phij, file="nao.s.trend.phij.Rdata", compress = "xz", compression_level = 9)
 # save(nao.s.trend.phia, file="nao.s.trend.phia.Rdata", compress = "xz", compression_level = 9)
 # save(nao.s.trend.all, file="nao.s.trend.all.Rdata", compress = "xz", compression_level = 9)



# BPVA results data for Table C1  ----------------------------------------------------------------------------------------------------------------------------------------

# load saved data
load(file="nao.trend.Rdata")
load(file="nao.trend.fec.Rdata")
load(file="nao.trend.phij.Rdata")
load(file="nao.trend.phia.Rdata")
load(file="nao.trend.all.Rdata")
load(file="nao.s.trend.Rdata")
load(file="nao.s.trend.fec.Rdata")
load(file="nao.s.trend.phij.Rdata")
load(file="nao.s.trend.phia.Rdata")
load(file="nao.s.trend.all.Rdata")


# Create table of pop size after 10 years, prob that pop size after
# 10 years is smaller than in last year of data collection, and quasi-extinction probability
# long-term mean
# Population size after 10 years
nao.trend$mean$Ntot[20]
nao.trend.phia$mean$Ntot[20]
nao.trend.phij$mean$Ntot[20]
nao.trend.fec$mean$Ntot[20]
nao.trend.all$mean$Ntot[20]


# credible intervals of projected pop size
c(nao.trend$q2.5$Ntot[20], nao.trend$q97.5$Ntot[20])
c(nao.trend.phia$q2.5$Ntot[20], nao.trend.phia$q97.5$Ntot[20])
c(nao.trend.phij$q2.5$Ntot[20], nao.trend.phij$q97.5$Ntot[20])
c(nao.trend.fec$q2.5$Ntot[20], nao.trend.fec$q97.5$Ntot[20])
c(nao.trend.all$q2.5$Ntot[20], nao.trend.all$q97.5$Ntot[20])


# probability that pop size in 2032 < pop size in 2022
sum(nao.trend$sims.list$Ntot[,20] < nao.trend$sims.list$Ntot[,10])/30000
sum(nao.trend.phia$sims.list$Ntot[,20] < nao.trend.phia$sims.list$Ntot[,10])/30000
sum(nao.trend.phij$sims.list$Ntot[,20] < nao.trend.phij$sims.list$Ntot[,10])/30000
sum(nao.trend.fec$sims.list$Ntot[,20] < nao.trend.fec$sims.list$Ntot[,10])/30000
sum(nao.trend.all$sims.list$Ntot[,20] < nao.trend.all$sims.list$Ntot[,10])/30000



# quasi extinction probability (prob that <=6 pairs in 2032)
sum(nao.trend$sims.list$Ntot[,20] < 7)/30000
sum(nao.trend.phia$sims.list$Ntot[,20] < 7)/30000
sum(nao.trend.phij$sims.list$Ntot[,20] < 7)/30000
sum(nao.trend.fec$sims.list$Ntot[,20] < 7)/30000
sum(nao.trend.all$sims.list$Ntot[,20] < 7)/30000

#short-term mean
# Population size after 10 years
nao.s.trend$mean$Ntot[20]
nao.s.trend.phia$mean$Ntot[20]
nao.s.trend.phij$mean$Ntot[20]
nao.s.trend.fec$mean$Ntot[20]
nao.s.trend.all$mean$Ntot[20]


# credible intervals of projected pop size
c(nao.s.trend$q2.5$Ntot[20], nao.s.trend$q97.5$Ntot[20])
c(nao.s.trend.phia$q2.5$Ntot[20], nao.s.trend.phia$q97.5$Ntot[20])
c(nao.s.trend.phij$q2.5$Ntot[20], nao.s.trend.phij$q97.5$Ntot[20])
c(nao.s.trend.fec$q2.5$Ntot[20], nao.s.trend.fec$q97.5$Ntot[20])
c(nao.s.trend.all$q2.5$Ntot[20], nao.s.trend.all$q97.5$Ntot[20])


# probability that pop size in 2032 < pop size in 2022
sum(nao.s.trend$sims.list$Ntot[,20] < nao.s.trend$sims.list$Ntot[,10])/30000
sum(nao.s.trend.phia$sims.list$Ntot[,20] < nao.s.trend.phia$sims.list$Ntot[,10])/30000
sum(nao.s.trend.phij$sims.list$Ntot[,20] < nao.s.trend.phij$sims.list$Ntot[,10])/30000
sum(nao.s.trend.fec$sims.list$Ntot[,20] < nao.s.trend.fec$sims.list$Ntot[,10])/30000
sum(nao.s.trend.all$sims.list$Ntot[,20] < nao.s.trend.all$sims.list$Ntot[,10])/30000



# quasi extinction probability (prob that <=6 pairs in 2032)
sum(nao.s.trend$sims.list$Ntot[,20] < 7)/30000
sum(nao.s.trend.phia$sims.list$Ntot[,20] < 7)/30000
sum(nao.s.trend.phij$sims.list$Ntot[,20] < 7)/30000
sum(nao.s.trend.fec$sims.list$Ntot[,20] < 7)/30000
sum(nao.s.trend.all$sims.list$Ntot[,20] < 7)/30000


# Quasi-extinction plot using BPVA results (Figure C1) ----------------------------------------------------------------------------------------------------------------------------------------

par(mfrow=c(1,1))
#quasi-extinction probability for 10 year projections for all scenarios
## long-term
### null
trend.null <- as.vector(c(sum(nao.trend$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.trend$sims.list$Ntot[,20] < 7)/30000))
plot(trend.null)
### adutrend survival
trend.phia <- as.vector(c(sum(nao.trend.phia$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.trend.phia$sims.list$Ntot[,20] < 7)/30000))
plot(trend.phia)
### juvenile survival
trend.phij <- as.vector(c(sum(nao.trend.phij$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.trend.phij$sims.list$Ntot[,20] < 7)/30000))
plot(trend.phij)
### fecundity
trend.fec <- as.vector(c(sum(nao.trend.fec$sims.list$Ntot[,11] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,12] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,13] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,14] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,15] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,16] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,17] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,18] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,19] < 7)/30000,
                         sum(nao.trend.fec$sims.list$Ntot[,20] < 7)/30000))
plot(trend.fec)
### combined
trend.all <- as.vector(c(sum(nao.trend.all$sims.list$Ntot[,11] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,12] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,13] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,14] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,15] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,16] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,17] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,18] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,19] < 7)/30000,
                         sum(nao.trend.all$sims.list$Ntot[,20] < 7)/30000))
plot(trend.all)

## short-term
### null
trend.s.null <- as.vector(c(sum(nao.s.trend$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.s.trend$sims.list$Ntot[,20] < 7)/30000))

### adutrend survival
trend.s.phia <- as.vector(c(sum(nao.s.trend.phia$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.s.trend.phia$sims.list$Ntot[,20] < 7)/30000))

### juvenile survival
trend.s.phij <- as.vector(c(sum(nao.s.trend.phij$sims.list$Ntot[,11] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,12] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,13] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,14] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,15] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,16] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,17] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,18] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,19] < 7)/30000,
                          sum(nao.s.trend.phij$sims.list$Ntot[,20] < 7)/30000))

### fecundity
trend.s.fec <- as.vector(c(sum(nao.s.trend.fec$sims.list$Ntot[,11] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,12] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,13] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,14] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,15] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,16] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,17] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,18] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,19] < 7)/30000,
                         sum(nao.s.trend.fec$sims.list$Ntot[,20] < 7)/30000))

### combined
trend.s.all <- as.vector(c(sum(nao.s.trend.all$sims.list$Ntot[,11] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,12] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,13] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,14] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,15] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,16] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,17] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,18] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,19] < 7)/30000,
                         sum(nao.s.trend.all$sims.list$Ntot[,20] < 7)/30000))

#----------------------------------------------------------------------------------------------------------------------------------------
# Code for fig quasi-extinction probabilities under all 10 climate/mgmt scenarios
#save plot
# units1 <- "BLTE_"
# fname <- paste0(units1[1],"BPVA_quasi-extinct_trend.pdf")
# pdf(fname, height=5, width=5)

year <- 1:10
nyears <- 10
plot(0, 0, ylim = c(0, 0.8), xlim = c(1,10), xlab = "Year", ylab = "Quasi-extinction probability", col = "black", type = "l",lwd = 2,  axes = F, frame = F)
axis(2)
axis(1, at = 1:(nyears), labels = year)
lines(x=c(1:(nyears)), trend.null, col = "#fde725", lwd = 2)
lines(x=c(1:(nyears)), trend.fec, col = "#5ec962", lwd = 2)
lines(x=c(1:(nyears)), trend.phij, col = "#21918c", lwd = 2)
lines(x=c(1:(nyears)), trend.phia, col = "#3b528b", lwd = 2)
lines(x=c(1:(nyears)), trend.all, col = "#440154", lwd = 2)
lines(x=c(1:(nyears)), trend.s.null, col = "#fde725", lty =2, lwd = 2)
lines(x=c(1:(nyears)), trend.s.fec, col = "#5ec962", lty =2,lwd = 2)
lines(x=c(1:(nyears)), trend.s.phij, col = "#21918c",lty =2, lwd = 2)
lines(x=c(1:(nyears)), trend.s.phia, col = "#3b528b",lty =2, lwd = 2)
lines(x=c(1:(nyears)), trend.s.all, col = "#440154", lty =2,lwd = 2)
legend(x = 0.5, y = 0.8, legend = c("No additional management", "Increase fecundity", "Increase juvenile survival", "Increase adult survival", "Combined management"), lty = c(1, 1),lwd = c(2, 2), col = c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"), bty = "n", cex = 0.75)
legend(x = 0.5, y = 0.5, legend = c("100-yr NAOI", "10-yr NAOI"), lty = c(1, 2),lwd = c(2, 2), bty = "n", cex = 0.75)

# dev.off()
