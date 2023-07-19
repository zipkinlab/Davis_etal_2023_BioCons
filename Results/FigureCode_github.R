########################################################################################
# Integrated population model (IPM) for St. Clair Flats Black Terns, 2013 - 2022
# Kayla Davis, Sarah Saunders, Stephanie Beilke, Erin Ford, Jenni Fuller, Ava Landgraf, and Elise Zipkin

# Adapted from original scripts by Michael Schaub & Marc Kéry (2021)
# Modified by K. Davis, 2022

# Post-processing to create figures (2 & 3) and table 2 data from the main text

########################################################################################

# setup ------------------------------------------------------------------------
#load libraries
library(openxlsx)
library(tidyverse)
library(coda)
library(ggpubr)
library(jagsUI)
library(scales)
library(plotrix)

# Population count data ----------------------------------------------------------------------------------------------------------------------------------------

dat1 <- read.xlsx("SCF-IPM_github.xlsx") # use data file from analysis

# breeding pairs and productivity
year <- dat1$Year
nyears <- length(year)
y <- dat1$MinPairs # the number of breeding pairs per year
j1 <- dat1$MinFledges # the number of fledglings recorded flying 
j2 <- dat1$Nanotag # the number of fledglings recorded via nanotag
R <- y # number of pairs/broods monitored
R_tag <- dat1$NumberTagged # number of fledglings tagged

# IPM results figure (Figure 2) ----------------------------------------------------------------------------------------------------------------------------------------

#load data
load(file = "scf.Rdata")

# Code for fig of counts vs. ests, annual adult and juv surv. prob, and fecundity

#save plot
# units1 <- "BLTE_"
# fname <- paste0(units1[1],"IPM_Plot.pdf")
# pdf(fname, height=5.2, width=8.5)

par(mfrow = c(2, 2), cex.axis = 1, cex.lab = 1, las = 1, mar = c(5, 5, 1, 5), mgp=c(3, 1, 0))

nyears <- 10
lower <- upper <- numeric()
year <- 2013:2022
for (i in 1:nyears){
  lower[i] <- quantile(scf$sims.list$Ntot[,i], 0.025)
  upper[i] <- quantile(scf$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(scf$mean$Ntot, y, lower), na.rm = T)
m2 <- max(c(scf$mean$Ntot, y, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, nyears), ylab = "Population size (pairs)", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:nyears, labels = year)
polygon(x = c(1:nyears, nyears:1), y = c(lower, upper[nyears:1]), col = alpha("#5ec962", 0.3), border = alpha("#5ec962", 0.3))
points(scf$mean$Ntot, type = "l", col = "#5ec962", lwd = 1.4)
points(y, type = "l", col = "#481f70", lwd = 1.4)
legend(x = 1, y = 100, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("#481f70", "#5ec962"), bty = "n", cex = 0.75)
corner.label("(a)", font=2, cex=1.5, xoff=0.2)


lower <- upper <- numeric()
T <- nyears
for (t in 1:T){
  lower[t] <- quantile(scf$sims.list$f1[,t], 0.025)
  upper[t] <- quantile(scf$sims.list$f1[,t], 0.975)}
lowerf2 <- upperf2 <- numeric()
for (t in 1:T){
  lowerf2[t] <- quantile(scf$sims.list$f2[,t], 0.025)
  upperf2[t] <- quantile(scf$sims.list$f2[,t], 0.975)}
plot(0,0, ylim = c(0, 0.6), xlim=c(1,10), ylab = "Fecundity (fledgling / female)", xlab = "", axes = F, cex = 1.1, frame = F, lwd = 1.3, col = '#b4de2c')
axis(2, ylim = c(0, 0.6))
axis(1, at = 1:(T), labels = 2013:2022)
par(new = T)
polygon(x = c(7:9, 9:7), y = c(lowerf2[7:9], upperf2[7:9]), col = alpha("#481f70", 0.3), border = alpha("#481f70", 0.3))
segments(1, mean(scf$sims.list$mfec1), T, mean(scf$sims.list$mfec1), lty = 1, lwd = 1.4, col = "#5ec962")
polygon(x = c(1:nyears, nyears:1), y = c(lower[1:nyears], upper[1:nyears]), col = alpha("#5ec962", 0.3), border = alpha("#5ec962", 0.3))
points(7, (dat1$Nanotag[7]/dat1$NumberTagged[7]), col = "#481f70", pch = 0)
points(9, (dat1$Nanotag[9]/dat1$NumberTagged[9]), col = "#481f70", pch = 0)
points(dat1$Productivity, col = '#5ec962', pch = 0)
legend(x = 0.5, y = 0.55, legend = c("Fledgling Mean", "Nanotag Mean"), lty = c(1, 1),lwd = c(2, 2), col = c("#5ec962", "#481f70"), bty = "n", cex = 0.75)
legend(x = 0.5, y = 0.44, legend = c("Observed Fledgling", "Observed Nanotag"), pch = c(0, 0),col = c("#5ec962", "#481f70"), bty = "n", cex = 0.75)
corner.label("(b)", font=2, cex=1.5, xoff=0.2)


lower <- upper <- numeric()
T <- nyears-3
for (t in 1:T){
  lower[t] <- quantile(scf$sims.list$phij[,t], 0.025)
  upper[t] <- quantile(scf$sims.list$phij[,t], 0.975)}
plot(0,0, ylim = c(0, 0.6), xlim = c(1,8), xaxt = "n", ylab = "Juvenile survival probability", xlab = "", axes = F, cex = 1.1, frame = F, lwd = 1.3)
axis(2)
axis(1, at = 1:(T+1), labels = 2015:2022)
polygon(x = c(1:(T+1), (T+1):1), y = c(rep(quantile(scf$sims.list$mphij, 0.025), (T+1)), rep(quantile(scf$sims.list$mphij, 0.975), (T+1))), col = alpha("#5ec962", 0.3), border = alpha("#5ec962", 0.3))
segments(1, scf$mean$mphij, T+1, scf$mean$mphij, lty = 1, lwd = 1.4, col = "#5ec962")
segments((1:T)+0.5, lower, (1:T)+0.5, upper, col = "#481f70")
points(y = scf$mean$phij, x = (1:T)+0.5, type = "b", pch = 16, col = "#481f70")
corner.label("(c)", font=2, cex=1.5, xoff=0.2)


lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
  lower[t] <- quantile(scf$sims.list$phia[,t], 0.025)
  upper[t] <- quantile(scf$sims.list$phia[,t], 0.975)}
plot(0,0, cex = 1.1, lwd = 1.3, ylim = c(0, 1.0), xlim = c(1,10), ylab = "Adult survival probability", xaxt = "n", xlab = "", frame = F)
axis(2)
axis(1, at = 1:(nyears), labels = 2013:2022)
polygon(x = c(1:(T+1), (T+1):1), y = c(rep(quantile(scf$sims.list$mphia, 0.025), (T+1)), rep(quantile(scf$sims.list$mphia, 0.975), (T+1))), col = alpha("#5ec962", 0.3), border = alpha("#5ec962", 0.3))
segments(1, scf$mean$mphia, T+1, scf$mean$mphia, lty = 1, lwd = 1.4, col = "#5ec962")
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
points(y = scf$mean$phia, x = (1:T)+0.5, type = "b", pch = 16, col = "#481f70")
corner.label("(d)", font=2, cex=1.5, xoff=0.2)

# dev.off()

# BPVA results data for Table 1  ----------------------------------------------------------------------------------------------------------------------------------------

# load saved data
load(file="nao.st.Rdata")
load(file="nao.lt.Rdata")
load(file="nao.st.fec.Rdata")
load(file="nao.lt.fec.Rdata")
load(file="nao.lt.phij.Rdata")
load(file="nao.st.phij.Rdata")
load(file="nao.st.phia.Rdata")
load(file="nao.lt.phia.Rdata")
load(file="nao.lt.all.Rdata")
load(file="nao.st.all.Rdata")


# Create table of pop size after 10 years, prob that pop size after
# 10 years is smaller than in last year of data collection, and quasi-extinction probability

# Population size after 10 years
nao.lt$mean$Ntot[20]
nao.st$mean$Ntot[20]
nao.lt.phia$mean$Ntot[20]
nao.st.phia$mean$Ntot[20]
nao.lt.phij$mean$Ntot[20]
nao.st.phij$mean$Ntot[20]
nao.lt.fec$mean$Ntot[20]
nao.st.fec$mean$Ntot[20]
nao.lt.all$mean$Ntot[20]
nao.st.all$mean$Ntot[20]

# credible intervals of projected pop size
c(nao.lt$q2.5$Ntot[20], nao.lt$q97.5$Ntot[20])
c(nao.st$q2.5$Ntot[20], nao.st$q97.5$Ntot[20])
c(nao.lt.phia$q2.5$Ntot[20], nao.lt.phia$q97.5$Ntot[20])
c(nao.st.phia$q2.5$Ntot[20], nao.st.phia$q97.5$Ntot[20])
c(nao.lt.phij$q2.5$Ntot[20], nao.lt.phij$q97.5$Ntot[20])
c(nao.st.phij$q2.5$Ntot[20], nao.st.phij$q97.5$Ntot[20])
c(nao.lt.fec$q2.5$Ntot[20], nao.lt.fec$q97.5$Ntot[20])
c(nao.st.fec$q2.5$Ntot[20], nao.st.fec$q97.5$Ntot[20])
c(nao.lt.all$q2.5$Ntot[20], nao.lt.all$q97.5$Ntot[20])
c(nao.st.all$q2.5$Ntot[20], nao.st.all$q97.5$Ntot[20])


# probability that pop size in 2032 < pop size in 2022
sum(nao.lt$sims.list$Ntot[,20] < nao.lt$sims.list$Ntot[,10])/30000
sum(nao.lt.phia$sims.list$Ntot[,20] < nao.lt.phia$sims.list$Ntot[,10])/30000
sum(nao.lt.phij$sims.list$Ntot[,20] < nao.lt.phij$sims.list$Ntot[,10])/30000
sum(nao.lt.fec$sims.list$Ntot[,20] < nao.lt.fec$sims.list$Ntot[,10])/30000
sum(nao.lt.all$sims.list$Ntot[,20] < nao.lt.all$sims.list$Ntot[,10])/30000
sum(nao.st$sims.list$Ntot[,20] < nao.st$sims.list$Ntot[,10])/30000
sum(nao.st.phia$sims.list$Ntot[,20] < nao.st.phia$sims.list$Ntot[,10])/30000
sum(nao.st.phij$sims.list$Ntot[,20] < nao.st.phij$sims.list$Ntot[,10])/30000
sum(nao.st.fec$sims.list$Ntot[,20] < nao.st.fec$sims.list$Ntot[,10])/30000
sum(nao.st.all$sims.list$Ntot[,20] < nao.st.all$sims.list$Ntot[,10])/30000

# population growth rate from 2022 to 2032
mean(nao.lt$sims.list$Ntot[,20] / nao.lt$sims.list$Ntot[,10])
mean(na.omit(nao.lt.phia$sims.list$Ntot[,20] / nao.lt.phia$sims.list$Ntot[,10]))
mean(nao.lt.phij$sims.list$Ntot[,20] / nao.lt.phij$sims.list$Ntot[,10])
mean(nao.lt.fec$sims.list$Ntot[,20] / nao.lt.fec$sims.list$Ntot[,10])
mean(nao.lt.all$sims.list$Ntot[,20] / nao.lt.all$sims.list$Ntot[,10])
mean(nao.st$sims.list$Ntot[,20] / nao.st$sims.list$Ntot[,10])
mean(nao.st.phia$sims.list$Ntot[,20] / nao.st.phia$sims.list$Ntot[,10])
mean(nao.st.phij$sims.list$Ntot[,20] / nao.st.phij$sims.list$Ntot[,10])
mean(nao.st.fec$sims.list$Ntot[,20] / nao.st.fec$sims.list$Ntot[,10])
mean(nao.st.all$sims.list$Ntot[,20] / nao.st.all$sims.list$Ntot[,10])
# lower intervals
quantile(nao.lt$sims.list$Ntot[,20] / nao.lt$sims.list$Ntot[,10], 0.025)
quantile(na.omit(nao.lt.phia$sims.list$Ntot[,20] / nao.lt.phia$sims.list$Ntot[,10]),0.025)
quantile(nao.lt.phij$sims.list$Ntot[,20] / nao.lt.phij$sims.list$Ntot[,10],0.025)
quantile(nao.lt.fec$sims.list$Ntot[,20] / nao.lt.fec$sims.list$Ntot[,10],0.025)
quantile(nao.lt.all$sims.list$Ntot[,20] / nao.lt.all$sims.list$Ntot[,10],0.025)
quantile(nao.st$sims.list$Ntot[,20] / nao.st$sims.list$Ntot[,10],0.025)
quantile(nao.st.phia$sims.list$Ntot[,20] / nao.st.phia$sims.list$Ntot[,10],0.025)
quantile(nao.st.phij$sims.list$Ntot[,20] / nao.st.phij$sims.list$Ntot[,10],0.025)
quantile(nao.st.fec$sims.list$Ntot[,20] / nao.st.fec$sims.list$Ntot[,10],0.025)
quantile(nao.st.all$sims.list$Ntot[,20] / nao.st.all$sims.list$Ntot[,10],0.025)
# upper intervals
quantile(nao.lt$sims.list$Ntot[,20] / nao.lt$sims.list$Ntot[,10], 0.975)
quantile(na.omit(nao.lt.phia$sims.list$Ntot[,20] / nao.lt.phia$sims.list$Ntot[,10]),0.975)
quantile(nao.lt.phij$sims.list$Ntot[,20] / nao.lt.phij$sims.list$Ntot[,10],0.975)
quantile(nao.lt.fec$sims.list$Ntot[,20] / nao.lt.fec$sims.list$Ntot[,10],0.975)
quantile(nao.lt.all$sims.list$Ntot[,20] / nao.lt.all$sims.list$Ntot[,10],0.975)
quantile(nao.st$sims.list$Ntot[,20] / nao.st$sims.list$Ntot[,10],0.975)
quantile(nao.st.phia$sims.list$Ntot[,20] / nao.st.phia$sims.list$Ntot[,10],0.975)
quantile(nao.st.phij$sims.list$Ntot[,20] / nao.st.phij$sims.list$Ntot[,10],0.975)
quantile(nao.st.fec$sims.list$Ntot[,20] / nao.st.fec$sims.list$Ntot[,10],0.975)
quantile(nao.st.all$sims.list$Ntot[,20] / nao.st.all$sims.list$Ntot[,10],0.975)

# quasi extinction probability (prob that <=6 pairs in 2032)
sum(nao.lt$sims.list$Ntot[,20] < 7)/30000
sum(nao.st$sims.list$Ntot[,20] < 7)/30000
sum(nao.lt.phia$sims.list$Ntot[,20] < 7)/30000
sum(nao.lt.phij$sims.list$Ntot[,20] < 7)/30000
sum(nao.lt.fec$sims.list$Ntot[,20] < 7)/30000
sum(nao.lt.all$sims.list$Ntot[,20] < 7)/30000
sum(nao.st.phia$sims.list$Ntot[,20] < 7)/30000
sum(nao.st.phij$sims.list$Ntot[,20] < 7)/30000
sum(nao.st.fec$sims.list$Ntot[,20] < 7)/30000
sum(nao.st.all$sims.list$Ntot[,20] < 7)/30000



# Quasi-extinction plot using BPVA results (Figure 3) ----------------------------------------------------------------------------------------------------------------------------------------

par(mfrow=c(1,1))
#quasi-extinction probability for 10 year projections for all scenarios
## long-term
### null
lt.null <- as.vector(c(sum(nao.lt$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.lt$sims.list$Ntot[,20] < 7)/30000))
plot(lt.null)
### adult survival
lt.phia <- as.vector(c(sum(nao.lt.phia$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.lt.phia$sims.list$Ntot[,20] < 7)/30000))
plot(lt.phia)
### juvenile survival
lt.phij <- as.vector(c(sum(nao.lt.phij$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.lt.phij$sims.list$Ntot[,20] < 7)/30000))
plot(lt.phij)
### fecundity
lt.fec <- as.vector(c(sum(nao.lt.fec$sims.list$Ntot[,11] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,12] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,13] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,14] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,15] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,16] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,17] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,18] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,19] < 7)/30000,
                      sum(nao.lt.fec$sims.list$Ntot[,20] < 7)/30000))
plot(lt.fec)
### combined
lt.all <- as.vector(c(sum(nao.lt.all$sims.list$Ntot[,11] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,12] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,13] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,14] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,15] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,16] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,17] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,18] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,19] < 7)/30000,
                      sum(nao.lt.all$sims.list$Ntot[,20] < 7)/30000))
plot(lt.all)

## short-term
### null
st.null <- as.vector(c(sum(nao.st$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.st$sims.list$Ntot[,20] < 7)/30000))
plot(st.null)
### adult survival
st.phia <- as.vector(c(sum(nao.st.phia$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.st.phia$sims.list$Ntot[,20] < 7)/30000))
plot(st.phia)
### juvenile survival
st.phij <- as.vector(c(sum(nao.st.phij$sims.list$Ntot[,11] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,12] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,13] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,14] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,15] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,16] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,17] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,18] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,19] < 7)/30000,
                       sum(nao.st.phij$sims.list$Ntot[,20] < 7)/30000))
plot(st.phij)
### fecundity
st.fec <- as.vector(c(sum(nao.st.fec$sims.list$Ntot[,11] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,12] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,13] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,14] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,15] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,16] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,17] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,18] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,19] < 7)/30000,
                      sum(nao.st.fec$sims.list$Ntot[,20] < 7)/30000))
plot(st.fec)
### combined
st.all <- as.vector(c(sum(nao.st.all$sims.list$Ntot[,11] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,12] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,13] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,14] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,15] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,16] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,17] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,18] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,19] < 7)/30000,
                      sum(nao.st.all$sims.list$Ntot[,20] < 7)/30000))
plot(st.all)
#----------------------------------------------------------------------------------------------------------------------------------------
# Code for fig quasi-extinction probabilities under all 10 climate/mgmt scenarios
#save plot
# units1 <- "BLTE_"
# fname <- paste0(units1[1],"BPVA_quasi-extinct_final.pdf")
# pdf(fname, height=5, width=5)

year <- 1:10
nyears <- 10
plot(0, 0, ylim = c(0, 0.8), xlim = c(1,10), xlab = "Year", ylab = "Quasi-extinction probability", col = "black", type = "l",lwd = 2,  axes = F, frame = F)
axis(2)
axis(1, at = 1:(nyears), labels = year)
lines(x=c(1:(nyears)), lt.null, col = "#fde725", lwd = 2)
lines(x=c(1:(nyears)), lt.fec, col = "#5ec962", lwd = 2)
lines(x=c(1:(nyears)), lt.phij, col = "#21918c", lwd = 2)
lines(x=c(1:(nyears)), lt.phia, col = "#3b528b", lwd = 2)
lines(x=c(1:(nyears)), lt.all, col = "#440154", lwd = 2)
lines(x=c(1:(nyears)), st.null, col = "#fde725", lty = 2, lwd = 2)
lines(x=c(1:(nyears)), st.fec, col = "#5ec962", lty = 2, lwd = 2)
lines(x=c(1:(nyears)), st.phij, col = "#21918c", lty = 2, lwd = 2)
lines(x=c(1:(nyears)), st.phia, col = "#3b528b", lty = 2, lwd = 2)
lines(x=c(1:(nyears)), st.all, col = "#440154", lty =2, lwd = 2)
legend(x = 0.5, y = 0.8, legend = c("No additional management", "Increase fecundity", "Increase juvenile survival", "Increase adult survival", "Combined management"), lty = c(1, 1),lwd = c(2, 2), col = c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154"), bty = "n", cex = 0.75)
legend(x = 0.5, y = 0.5, legend = c("100-yr NAOI", "10-yr NAOI"), lty = c(1, 2),lwd = c(2, 2), bty = "n", cex = 0.75)
# dev.off()
