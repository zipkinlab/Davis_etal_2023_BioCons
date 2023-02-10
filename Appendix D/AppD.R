########################################################################################
# Integrated population model (IPM) for St. Clair Flats Black Terns, 2013 - 2022
# Kayla Davis, Sarah Saunders, Stephanie Beilke, Erin Ford, Jenni Fuller, Ava Landgraf, and Elise Zipkin

# Adapted from original scripts by Michael Schaub & Marc Kéry (2021)
# Modified by K. Davis, 2022

# Appendix D code 

########################################################################################
getwd()
setwd("G:/My Drive/BLTE/GenericIPM_Practice")
options(max.print=99999)


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



# the data ---------------------------------------------------------------------

########################################################################
# Capture-recapture data: m-array of juveniles (HY) and adults (AHY)
########################################################################

#Create NA marrays
marray.j <- matrix(0, nrow = 7, ncol = 10)
marray.a <- matrix(0, nrow = 9, ncol = 10)

########################################################################
# Population count data
########################################################################

dat1 <- read.xlsx("NA-IPM.xlsx")

# breeding pairs and productivity
year <- dat1$Year
nyears <- length(year)
y <- dat1$MinPairs # the number of breeding pairs per year
j <- dat1$MinFledges # the number of fledglings recorded flying 
R <- rep(0,nyears) # number of pairs/broods monitored

###############################################################################################################################
## Check what priors look like
#############################################################################################################################

#define range
p = seq(0,1, length=100)

#plot beta priors for adult and juvenile survival

#save plot
# units1 <- "priors_"
# fname <- paste0(units1,"BetaSurvival.pdf")
# pdf(fname, height=3.5, width=5)

plot(p, dbeta(p, 13.5, 4), ylab='density', type ='l', col='#5ec962', lwd = 2)
lines(p, dbeta(p, 3, 12), col='#481f70', lwd = 2) 
legend(x = 0.3, y = 4, 
       legend = c("Beta prior phi[ad]", "Beta prior phi[juv]"), 
       lty = c(1, 1),lwd = c(2, 2), col = c("#5ec962", '#481f70'),
       bty = "n", cex = 0.75)

# dev.off()


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
# marray.j is an empty m-array of capture histories for individuals first banded as
# chicks during 2013 - 2022
# marray.a is an empty m-array of capture histories for individuals first banded as
# adults during 2013 - 2022
# y = number of breeding pairs annually set to Na
# j1 = number of fledglings observed annually from on-the-ground counts set to NA
# R = number of pairs/broods monitored annually set to NA


#############################################################################

sink("NA_ipm")
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
    # fec2 ~ dunif(0,3)                    # productivity
    # l.mfec2 <- log(fec2)                 # Log transformation
    res ~ dunif(0,1)                     # mean detection probability
    l.p <- log(res / (1-res))            # Logit transformation
   

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
    logit(phia[t]) <- l.mphia + epsilon.phia[t]  # epsilon.phia is random temporal effect for env. stoch.                           
    logit(p[t]) <- l.p + epsilon.res[t] 
    }
    
    for (t in 1:nyears){
    log(f1[t]) <- l.mfec1       # f1 = fecundity from fledgling counts (this is used in the IPM)
    # log(f2[t]) <- l.mfec2     # f2 = fecundity from nanotag counts (this is estimated outside of the IPM for comparison to fledgling count)
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
    # mfec2 <- exp(l.mfec2)                       # Mean productivity nanotag counts

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
    
    # j2[t] ~ dpois(rho2[t])         # number young fledged with nanotag
    # rho2[t] <- R_tag[t] * f2[t]    # number tagged and fecundity

    }
    
    }
    ",fill = TRUE)
sink()

###################################################################
# Bundle data
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, r.j = rowSums(marray.j), r.a = rowSums(marray.a), y = y, j1 = j1,R = R)  

# Initial values
inits <- function(){list(mphi.juv = runif(0.1, 0.2, 0.4), mphi.ad = runif(1, 0.7, 0.9), fec1 = runif(1, 0, 2), res = runif(1, 0, 0.5), sig.phij = runif(1, 0.1, 5), sig.phia = runif(1, 0.1, 5), 
                         sigma.obs = runif(1, 0, 1))} 


# Parameters monitored
parameters <- c("phij", "phia","f1", "lambda", 
                "mphij", "mphia","mfec1","mlam",
                "l.mphij", "l.mphia","l.mfec1","p",
                "sig.phij", "sig.phia", "sig.obs", 
                "N1", "Nad", "Ntot") 

# MCMC settings
ni <- 10000   
nt <- 10
nb <- 4000
nc <- 3

# Call JAGS from R
null <- jags(jags.data, inits, parameters, "NA_ipm", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)


