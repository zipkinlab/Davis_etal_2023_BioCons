########################################################################################
# Integrated population model (IPM) for St. Clair Flats Black Terns, 2013 - 2022
# Kayla Davis, Sarah Saunders, Stephanie Beilke, Erin Ford, Jenni Fuller, Ava Landgraf, and Elise Zipkin

# Adapted from original scripts by Michael Schaub & Marc Kéry (2021)
# Modified by K. Davis, 2022

########################################################################################

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


# the data ---------------------------------------------------------------------

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
yeareffects = read.csv("CovData/YearCovs.csv", header=T, sep=',', na.strings=T)

# adult survival covs
nao_hur = yeareffects$nao_jan.jun

# make into vector for jags
nao_hur = as.numeric(as.vector(nao_hur[1:9]))

#calculate mean nao
mnao <- mean(nao_hur)
sdnao <- sd(nao_hur)

# standardize effects
znao_hur <- as.vector(scale(nao_hur))


# the model ---------------------------------------------------------------------


#############################################################################
# Integrated population model (IPM) for St. Clair Flats Black Terns
# Code by Kayla Davis, Michigan State University, 2022
# Data provided by Detroit Audubon
# Adapted from original scripts by Michael Schaub & Marc Kéry (2021)
# Modified by K. Davis, 2022
# See main text for full description of modeling framework
#
# Notations:
# nyears = 10
# marray.j is an m-array of capture histories for individuals first banded as
# chicks during 2013 - 2022
# marray.a is an m-array of capture histories for individuals first banded as
# adults during 2013 - 2022
# y = number of breeding pairs annually
# j1 = number of fledglings observed annually from on-the-ground counts
# j2 = number of fledglings observed annually (2019 and 2021) from nanotag observations
# R = number of pairs/broods monitored annually
# R_tag = number of pre-fledged chicks tagged with nanotags in 2019 and 2021
# znao_hur = North Atlantic Oscillation index preceeding the hurricane season of year t-1
# (i.e. NAOI in January-June before breeding season of year t-1)

#############################################################################

sink("blte_ipm")
cat("
model {
    #------------------------------------------------------------------------
    #  Integrated population model
    #  - Stage-structured model with 2 stages: 
    #       juvenile (age 0-3), adults (age 4+)
    #  - Age at first breeding is third year
    #  - 1 and 2 y olds are 'invisible' (don't return to breeding grounds
    #       until age 3)
    #  - Prebreeding census, female-based
    #  - All vital rates are assumed to be time-dependent
    #  - Includes env. stochasticity thru random time effects for all params except fecundity
    #------------------------------------------------------------------------

    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes

     for (t in 1:3){
     n1[t] ~ dnorm(10, 0.1)I(0,)               # New 3 year olds; needs prior for first 3 years
     N1[t] <- round(n1[t])               
     Ntot[t] <- Nad[t] + N1[t]                 # Initial population size of adults in first 3 years is sum of new 3yos and returning adults
     }#t  
    
     N0[1] ~ dpois((f1[1]*0.5) * Ntot[1])      # Initial pop size of fledglings

     nadSurv ~ dnorm(300, 0.01)I(0,)           # Adults in year 1
     Nad[1] <- round(nadSurv)
  
    
    # Mean demographic parameters (on appropriate scale)

    mphi.juv ~ dbeta(3, 12)
    mphi.ad ~ dbeta(13.5, 4)
    l.mphij <- log(mphi.juv / (1-mphi.juv)) # Logit transformation
    l.mphia <- log(mphi.ad / (1-mphi.ad))   # Logit transformation
    
    fec1 ~ dunif(0,3)                    # productivity
    l.mfec1 <- log(fec1)                 # Log transform
    fec2 ~ dunif(0,3)                    # productivity
    l.mfec2 <- log(fec2)                 # Log transformation
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
    

    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1)){
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.res[t] ~ dnorm(0, tau.res)T(-5,5)
    }
    
    for (t in 1:(nyears-3)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)
    }
    
    
    #---------------------------------------------
    # 2. Constrain parameters (temp variability)
    #---------------------------------------------
    for (t in 1:(nyears-1)){
    logit(phia[t]) <- l.mphia + beta.phia1 * znao_hur[t] + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    logit(p[t]) <- l.p + epsilon.res[t] 
    }
    
    for (t in 1:nyears){
    log(f1[t]) <- l.mfec1     # f1 = fecundity from fledgling counts (this is used in the IPM)
    log(f2[t]) <- l.mfec2     # f2 = fecundity from nanotag counts (this is estimated outside of the IPM for comparison to fledgling count)
    }

    
    for (t in 1:(nyears-3)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]         
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))        # Mean juvenile survival 
    mphia <- exp(l.mphia)/(1+exp(l.mphia))        # Mean adult survival 
    mfec1 <- exp(l.mfec1)                         # Mean productivity fledgling counts
    mfec2 <- exp(l.mfec2)                         # Mean productivity nanotag counts

    # Population growth rate (total adult breeders [3+ y olds])
    for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / (Ntot[t] + 0.001)   
    logla[t] <- log(lambda[t])
    }

    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geo mean all yrs
  
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (k in 4:nyears){
    N1[k] ~ dbin(phij[k-3], N0[k-3])
    Ntot[k] <- N1[k] + Nad[k]
    }#k
    
    for (t in 2:nyears){
    mean1[t] <- (f1[t]*0.5) * Ntot[t]       # all adults can breed
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
    j1[t] ~ dpois(rho1[t])         # number young observed as fledged
    rho1[t] <- f1[t] * R[t]        # number of pairs and fecundity
    
    j2[t] ~ dpois(rho2[t])         # number young fledged with nanotag
    rho2[t] <- R_tag[t] * f2[t]    # number tagged and fecundity

    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
jags.data <- list(znao_hur = znao_hur, nyears = nyears, marray.j = marray.j, marray.a = marray.a, r.j = rowSums(marray.j), r.a = rowSums(marray.a), y = y, j1 = j1, j2 = j2, R = R, R_tag = R_tag)  

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), fec2 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sig.fec1 = runif(1, 0.1, 10), sig.fec2 = runif(1, 0.1, 10), sigma.obs = runif(1, 0, 1), beta.phia1 = rnorm(1,0,1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "f2", "lambda", 
                "mphij", "mphia","mfec1","mfec2", "mlam",
                "l.mphij", "l.mphia","l.mfec1","l.mfec2", "p",
                "sig.phij", "sig.phia", "sig.obs", 
                "N1", "Nad", "Ntot", "beta.phia1") 

# MCMC settings
ni <- 500000    
nt <- 10
nb <- 400000
nc <- 3

# Call JAGS from R
scf <- jags(jags.data, inits, parameters, "blte_ipm", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)

# ~~~~ save output for use later ~~~~
#save(scf, file="scf.Rdata")

#use MCMC vis to look at different credible intervals (80% and 95%)
scf

MCMCsummary(scf,
            params = c("beta.phia1"),
            probs = c(0.025, 0.05, 0.075, 0.1,0.25, 0.5, 0.75, 0.9, 0.925, 0.95, 0.975),
            round = 2)


# BPVA ---------------------------------------------------------------------


###############################################################################################################################
## Bayesian population viability analysis with top-supported model covariates
#############################################################################################################################


###################################################
## Bootstrap method to obtain temporal variation
## in NAO values over time
##################################################

## Read in data
nao_all <- read.csv("CovData/Hist_NAO_Means.csv", stringsAsFactors = FALSE)
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



###################################################
# Combined mgmt scenario: 
# increase adult survival 5%
# increase juv survival 10%
# double fecundity
##################################################

sink("blte_nao_bpva_all")
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
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    nao.pre[t] ~ dnorm(lt.mean,tau.nao)     #use short term mean(st.mean) or long term mean (lt.mean) here   
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
nao.lt.all <- jags(jags.data, inits, parameters, "blte_nao_bpva_all", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)


###################################################
# Adult survival mgmt scenario: 
# increase adult survival 5%
##################################################

sink("blte_nao_bpva_phia")
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
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    nao.pre[t] ~ dnorm(lt.mean,tau.nao)     #use short term mean(st.mean) or long term mean (lt.mean) here   
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

nao.lt.phia <- jags(jags.data, inits, parameters, "blte_nao_bpva_phia", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# Juvenile survival mgmt scenario: 
# increase juv survival 10%
##################################################

sink("blte_nao_bpva_phij")
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
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    nao.pre[t] ~ dnorm(lt.mean,tau.nao)     #use short term mean(st.mean) or long term mean (lt.mean) here   
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

nao.lt.phij <- jags(jags.data, inits, parameters, "blte_nao_bpva_phij", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# Fecundity mgmt scenario: 
# double annual fecundity
##################################################

sink("blte_nao_bpva_fec")
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
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    nao.pre[t] ~ dnorm(lt.mean,tau.nao)     #use short term mean(st.mean) or long term mean (lt.mean) here   
    nao.fut[t] <- (nao.pre[t] - mnao)/sdnao #standardize?
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

nao.lt.fec <- jags(jags.data, inits, parameters, "blte_nao_bpva_fec", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



###################################################
# No additional mgmt scenario: 
# Same for past and future
##################################################

sink("blte_nao_bpva_null")
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
    
    for (t in nyears:(nyears-1+K)){         #pull future NAO values from prior
    nao.pre[t] ~ dnorm(st.mean,tau.nao)     #use short term mean(st.mean) or long term mean (lt.mean) here   
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

nao.st <- jags(jags.data, inits, parameters, "blte_nao_bpva_null", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)



#----------------------------------------------------------------------------------------------------------------------------------------
# ~~~~ save output for use later ~~~~
# save(nao.st, file="nao.st.Rdata")
# save(nao.lt, file="nao.lt.Rdata")
# save(nao.st.fec, file="nao.st.fec.Rdata")
# save(nao.lt.fec, file="nao.lt.fec.Rdata")
# save(nao.lt.phij, file="nao.lt.phij.Rdata")
# save(nao.st.phij, file="nao.st.phij.Rdata")
# save(nao.st.phia, file="nao.st.phia.Rdata")
# save(nao.lt.phia, file="nao.lt.phia.Rdata")
# save(nao.lt.all, file="nao.lt.all.Rdata")
# save(nao.st.all, file="nao.st.all.Rdata")


