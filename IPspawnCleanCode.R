# SK 05/08/2015

###############################################################################
###############################################################################


###############################################################################
###############################################################################

library(CoralSpawn)
library(reshape)
library(plyr)


###############################################################################
###############################################################################

## 1. LOAD DATA

###############################################################################
###############################################################################

# raw data has been scaled to mean = 0, sd = 1
load("DataForAnalysisCoralSpawnJune2015.RData")
head(d)

d$current.speed <- sqrt(d$currentu.abs^2 + d$currentv.abs^2)
# standardise current speed data
d$current.speed <- sqrt(d$current.speed)  # transform first
d$current.speed <- scale(d$current.speed) # scale/centre

d <- d[,c(1:18,23,19:22)]


###############################################################################
###############################################################################

## 2. CHECK MULTICOLINEARITY

###############################################################################
###############################################################################

cor(d[,7:18])
pairs(d[,7:18])
# FYI - also check using VIF in next step

# only the variables included in the final models
cor(d[,c(8,10,11,13,14,15,16,19)])
pairs(d[,c(8,10,11,13,14,15,16,17,18,23)])
# no correlations >0.6

###############################################################################
###############################################################################

## 3. SPAWNING SEASONALITY MODEL 

###############################################################################
###############################################################################

# reduce to just sites with data for year-round spawning
ds <- d[!is.na(d$Spawn),]
colnames(ds)
nrow(ds)

## LINEAR OR QUADRATIC TERMS?
lin.quad(7:19,ds,binomial)
# quadratic wind.avg

# GLM with full model
# Choose variables that are most expected to make sense
# Want the simplest possible model
form.spawn <- Spawn ~ sst.avg.rel+sst.rel.change+par10+rain.avg+par.rel.change+
                         current.speed+sunset+poly(wind.avg,2)
glm.mod <- glm.diagnostic(form.spawn,binomial,ds)
# PAR relative change VIF >2.5, drop this variable

## Check interactions that biological/climatic make sense
spawn.int <- glm(Spawn ~ sst.avg.rel+sst.rel.change+par10+rain.avg+sunset+
                    current.speed+poly(wind.avg,2)+
                    current.speed*sst.rel.change+current.speed*poly(wind.avg,2)+
                    par10*sst.rel.change+sst.avg.rel*rain.avg+
                    sst.rel.change*rain.avg,
                    family=binomial(link="logit"),data=ds,na.action=na.fail)
summary(spawn.int)
confint(spawn.int)
# based on confidence intervals not overlapping zero, significant interactions:
# currentv.abs * wind
# currentu.abs * wind^2

## MIXED EFFECTS MODEL WITH SIGNIFICANT INTERACTIONS & QUADRATICS
form.spawn.mix <- Spawn ~ sst.avg.rel+sst.rel.change+rain.avg+par10+current.speed+
                           poly(wind.avg,2)+sunset+
                           current.speed*poly(wind.avg,2)+
                           (1|month)+(1|Site.number)

mixed.mod.output(form.spawn.mix,binomial,ds,c("month","Site.number"))

# model produces convergence errors.
# check convergence with alternative method from B. Bolker
# If < 0.001 should be ok
fit <- glmer(form.spawn.mix,binomial,data=ds,na.action=na.fail)
relgrad <- with(fit@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# result is 0.01 so definitely not converging.
# N.B. result is actually the as in error message after "max!grad!"
# Model needs to be simplified.

# remove interactions &/or quadratic
form.spawn.mix <- Spawn ~ sst.avg.rel+sst.rel.change+par10+
                           current.speed+poly(wind.avg,2)+
                           (1|month)+(1|Site.number)
mixed.mod.output(form.spawn.mix,binomial,ds,c("month","Site.number"))
# able to add currentv.abs*wind back in (interaction was with linear term of wind)

###############################################################################
## MODEL SELECTION AND AVERAGING

# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects.
# Interactions not significant in full mixed model were removed before this step.
# Ranked by BIC so strong penalisation for complexity.

form.spawn.avg <-  Spawn ~ sst.rel.change+par10+current.speed+poly(wind.avg,2)+
                            sst.avg.rel+(1|month)+(1|Site.number)

### USE NEW RELATIVE VAR THINGY!
sa.fit.spawn <- select.average("mixed",form.spawn.avg,binomial,ds)

coefTable(sa.fit.spawn$averaged.model,full=T,adjust.se=T)

## CONFIDENCE INTERVALS AND ODDS RATIOS
# Get 95% confidence intervals and odds ratios for coefficients.
oddsspawn <- cbind(exp(summary(sa.fit.spawn$averaged.mod)$coefmat.full[,1]),
                   exp(summary(sa.fit.spawn$averaged.mod)$coefmat.full[,1]-1.96*summary(sa.fit.spawn$averaged.mod)$coefmat.full[,3]), 
                   exp(summary(sa.fit.spawn$averaged.mod)$coefmat.full[,1]+1.96*summary(sa.fit.spawn$averaged.mod)$coefmat.full[,3]))
colnames(oddsspawn) <- c("or","l95","u95")
oddsspawn

###############################################################################


# PARTIAL COEFFICIENT PLOTS for significant variables 
# doesn't work with model-averaged result so generate best model
# first based on model selection step above
best.spawn.mod <- glmer(Spawn ~ sst.rel.change+sst.avg.rel+poly(wind.avg,2)+
                           (1|month)+(1|Site.number),
                           family=binomial(link="logit"),data=ds,na.action=na.fail)

summary(best.spawn.mod)
# R SQAURED & AICc OF THE BEST PAR MODEL
# R2m = variation explained by fixed effects; 
# R2c = variance explained by fixed + random
r.squaredGLMM(best.spawn.mod)
BIC(best.spawn.mod)
dotplot(ranef(best.spawn.mod,condVar=TRUE,whichel="month"))
dotplot(ranef(best.spawn.mod,condVar=TRUE,whichel="Site.number"))

# plot all partial coefficients on the same screen - cool :)
zspawn <- allEffects(mod=best.spawn.mod,default.levels=50,confidence.level=0.95)
# pdf("AllEffectsSpawnMo1delAllVarsFromAveragedModel.pdf")
plot(zspawn,rug=F,rescale.axis=F,ylab="p(spawning)",ylim=c(0,1),multiline=T,
     ci.style="bands",band.transparency=0.1,colors=c("black","white",
                     "blue","white","purple","white","green","white"),lines=1)
plot(zspawn,rug=F,rescale.axis=F,ylab="p(spawning)",ylim=c(0,1))
# dev.off()


save(ds,sa.fit.spawn,oddsspawn,best.spawn.mod,zspawn,file="ModelOutSpawn.RData")




#######################################################
#######################################################

## 3. PEAK SPAWNING MODEL

#######################################################
#######################################################

# cloglog link is good for dependent variable with rare events
# because it is asymmetric so will use here

# reduce to just sites with data for year-round spawning
dp <- d[!is.na(d$Peak),]

## LINEAR OR QUADRATIC TERMS?
lin.quad.peak(7:19,dp,binomial)
cor(dp[,7:19])
# linear for all models

# GLM with full model
form.peak <- Peak ~ sst.avg.rel+sst.rel.change+par10+rain.avg+
                         wind.avg+sunset+current.speed
glm.modp <- glm.diagnostic(form.peak,binomial(link="cloglog"),dp)

## check for interactions that are make sense
peak.int <- glm(Peak ~ sst.avg.rel+sst.rel.change+par10+rain.avg+wind.avg+sunset+
                    current.speed+current.speed*sst.rel.change+current.speed*wind.avg+
                    par10*sst.rel.change+sst.avg.rel*rain.avg+sst.rel.change*rain.avg,
                    family=binomial(link="cloglog"),data=dp,na.action=na.fail)
summary(peak.int)
confint(peak.int)
# based on confidence intervals not overlapping zero, no significant interactions

## MIXED EFFECTS MODEL WITH SIGNIFICANT INTERACTIONS 
form.peak.mix <- Peak ~ sst.avg.rel+sst.rel.change+par10+rain.avg+wind.avg+current.speed+
                         (1|month)

mixed.mod.output(form.peak.mix,binomial(link="cloglog"),dp,"month")


###############################################################################
## MODEL SELECTION AND AVERAGING

# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects.
# Interactions not significant in full mixed model were removed before this step.
# Ranked by BIC so strong penalisation for complexity.

form.peak.avg <- Peak ~ sst.avg.rel+sst.rel.change+par10+rain.avg+wind.avg+
                         current.speed+(1|month)

# model averaging uses "revised.var=TRUE" so that CIs work with shirnkage
sa.fitp <- select.average("mixed",form.peak.avg,binomial(link="cloglog"),dp)



###############################################################################

# PARTIAL COEFFICIENT PLOTS for significant variables 
# doesn't work with model-averaged result so generate best model
# first based on model selection step above
best.peak.mod <- glmer(Peak ~ sst.rel.change+current.speed+(1|month),
                        family=binomial(link="cloglog"),data=dp,na.action=na.fail)

summary(best.peak.mod)
# R SQAURED & AICc OF THE BEST PAR MODEL
# R2m = variation explained by fixed effects; 
# R2c = variance explained by fixed + random
r.squaredGLMM(best.peak.mod)
BIC(best.peak.mod)
dotplot(ranef(best.peak.mod,condVar=TRUE,whichel="month"))

# plot all partial coefficients on the same screen - cool :)
# make sure it uses cloglog and its inverse
zpeak <- allEffects(mod=best.peak.mod,default.levels=50,confidence.level=0.95,
                    transformation=list(link=function(x) log(-log(1-x)),inverse=function(x) 1-exp(-exp(x))))
# pdf("AllEffectsSpawnModelAllVarsFromAveragedModel.pdf")
plot(zpeak,rug=F,rescale.axis=F,ylab="p(peak)",ylim=c(0,1))
# dev.off()

save(dp,sa.fitp,best.peak.mod,zpeak,file="ModelOutPeak.RData")


#########################################################
#########################################################

## 4. FIGURES

#########################################################
#########################################################

library(plotrix)

load("ModelOutPeak.RData")
load("ModelOutSpawn.RData")
load("ReducedSitesEnvRaw.RData")

d2 <- d[order(d$Lat),]
d3 <- d[order(d$Long),]


#########################################################
## histogram of spawning duration

# How long does spawning last in each site?
spawn.m <- cast(d2,Site.number~month,value="Spawn")
peak.m <- cast(d2,Site.number~month,value="Peak")
rownames(spawn.m) <- spawn.m[,1]
rownames(peak.m) <- peak.m[,1]
spawn.m <- spawn.m[,-1]; peak.m <- peak.m[,-1]
month.name <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
colnames(spawn.m) <- month.name
colnames(peak.m) <- month.name
rowSums(spawn.m)
spawn.m
hist(rowSums(spawn.m),xlim=c(0,12),main="spawning season length",
     xlab="number of months",ylab="number of sites")


#########################################################
## wheel of spawning times

# scaling factor to plot circles concentrically
site.ring <- seq(0.25,1,0.9/34)
# colours for line segments
months <- factor(rep(1:12, each = 100))
# matrix to hold segments for each site
col.sites <- matrix(NA,1200,34)
# generate a bunch of angles to draw the circle over
angles <- rev(seq(0, 2*pi, length=1200))
angles2 <- angles[c(901:1200,1:900)]

# create empty plot
plot(cos(angles),sin(angles),col=months,type="n",axes="n",frame.plot=F)

for(z in 1:nrow(spawn.m)){
   # fill in colours of each site
   col.sites[,z] <- as.numeric(rep(spawn.m[z,], each = 100))   
   for(i in 1:1199){
      # multiply the angles by a scaling factor for different sites
      segments(cos(angles2)[i]*site.ring[z], sin(angles2)[i]*site.ring[z], cos(angles2)[i+1]*site.ring[z], 
               sin(angles2)[i+1]*site.ring[z], col=col.sites[i,z], lwd=1)
   }
}

month.colours <- c("blue","powderblue","aquamarine","seagreen3","seagreen","darkolivegreen4","greenyellow","thistle1","plum2","maroon3","purple","purple4")
# add a circle around the outside with different colours for each month
for(i in 1:1199){
      segments(cos(angles2)[i]*1.05, sin(angles2)[i]*1.05, cos(angles2)[i+1]*1.05, 
               sin(angles2)[i+1]*1.05, col=month.colours[months[i]], lwd=2)
}

# prepare positions for text labels for each month
month.mid <- seq(50,1150,100)
text.month <- angles2[month.mid] # midpoints of the text
cwise <- c(1,2,3,10,11,12) # which labels will be plotted clockwise?

# add points for peak spawning
for(i in 1:34){
      points(cos(text.month[peak.m[i,]>0])*site.ring[i],
             sin(text.month[peak.m[i,]>0])*site.ring[i],pch=19,col=2)
}

# add month names
for(i in 1:length(month.name)){
   if(i %in% cwise) arctext(month.name[i], center=c(0,0), radius=0.9, middle=text.month[i],"clockwise"=T)
   else arctext(month.name[i], center=c(0,0), radius=0.9, middle=text.month[i],"clockwise"=F)
}


##################################################################
## plot odds ratios from both averaged regression models


cilor <- log(oddsspawn[c(4,2,3,5),])

# plot empty plot to add points and lines to
# make sure margins will be large enough to fit in labels
par(mar=c(5,15,2,2))
plot(c(-3,3),c(1,5),cex.lab=1.5,col="white",xlab="Log Odds Ratio",ylab="",yaxt="n")
# create labels for environmental variables
yax.lab<-c("SST relative change","wind speed","wind speed^2",
           "longitudinal current speed")
# create y axis
axis(2,at=c(1,2,3,4),labels=rev(yax.lab),pos=-3.5,las=2,col="white")

# plot points for log odds ratios
# These are odds ratio for x and place of category for y
y<-c(1,2,3,4)

for(i in 1:(nrow(cilor))){
   points(cilor[i,1],rev(y)[i])
}

# add lower and upper confidence limits as lines
for(i in 1:(nrow(cilor))){
	lines(x=c(cilor[i,2],cilor[i,3]),y=c(rev(y)[i],rev(y)[i]))
}

# add a dashed line at zero for reference
abline(v=0,lty="dashed")

# add a legend for colour and barrier
legend("topleft",as.character(c("spawning season","peak month")),col=c(1,2),cex=1,lty=1,bg="white",title.adj=0.5,title="Model")



#########################################################
## plot odds ratios from peak averaged model

cilorp <- log(oddspeak[c(2,3),])

# plot points for log odds ratios
points(cilorp[1,1],3.25,col=2)
points(cilorp[2,1],1.25,col=2)

# add lower and upper confidence limits as lines
lines(x=c(cilorp[1,2],cilorp[1,3]),y=c(3.25,3.25),col=2)
lines(x=c(cilorp[2,2],cilorp[2,3]),y=c(1.25,1.25),col=2)



#########################################################
## temperature profiles per site per month,
## peak spawning indicated by filled bar

## PLOTS TO EXPLORE DATA, NOT FOR MANUSCRIPT

load("ReducedSitesEnvRaw.RData")

par(mfcol=c(3,4))
for(i in 1:length(unique(dp$Site.number))){
   x <- subset(dp,dp$Site.number==unique(dp$Site.number)[i])
   cols.spawn <- x$Peak
   barplot(x$sst.avg.rel,col=cols.spawn,main=x$Location[1])   
}

all.raw <- cbind(spawn.data,env.raw.reduced[,7:18])
all.raw$currentv.abs <- abs(all.raw$currentv)
all.raw$currentu.abs <- abs(all.raw$currentu)
all.raw <- all.raw[,-(18:20)]

par(mfcol=c(3,4))
# Find mean annual SST for each site
site.index2 <- seq(1,34*12,12)
mean.ann.sst <- c()
for(i in 1:length(site.index2)){
   sst12 <- all.raw$sst.avg[site.index2[i]:(site.index2[i]+11)]
   mean.ann.sst[site.index2[i]:(site.index2[i]+11)] <- rep(mean(sst12),12)    
}
# add to dataframe
all.raw <- cbind(all.raw,mean.ann.sst)

# reduce to sites that have peak spawning available
peak.raw <- all.raw[!is.na(all.raw$Peak),]

#pdf("SST avg profiles by site.pdf")
par(mfcol=c(3,4))
for(i in 1:length(unique(peak.raw$Site.number))){
   x <- subset(peak.raw,peak.raw$Site.number==unique(peak.raw$Site.number)[i])
   cols.spawn <- x$Peak
   barplot(x$sst.avg,col=cols.spawn,main=x$country[1],ylim=c(0,32),
           ylab="degrees C",xlab="month",names=month.name,las=2)   
}
#dev.off()


#pdf("SST rel change profiles by site SN.pdf")
par(mfcol=c(3,4))
for(i in 1:length(unique(peak.raw$Site.number))){
   x <- subset(peak.raw,peak.raw$Site.number==unique(peak.raw$Site.number)[i])
   cols.spawn <- x$Peak
   barplot(x$sst.rel.change,col=cols.spawn,main=x$Site.number[1],ylim=c(0,0.5),
           ylab="rel change",xlab="month",names=month.name,las=2)   
}
#dev.off()


# PLOT SST RELATIVE CHANGE IN MONTH OF PEAK SPAWNING BY LATITUDE
par(mfrow=c(1,2))
peak1 <- peak.raw[peak.raw$Peak==1,]
plot(peak1$Lat,peak1$sst.rel.change,col=1,pch=19,cex=0.5)
peaklat <- peak1[order(peak1$Lat),]
lines(loess.smooth(peaklat$Lat,peaklat$sst.rel.change),lty=2)

# PLOT SST RELATIVE CHANGE IN MONTH OF PEAK SPAWNING BY LONGITUDE
peak1 <- peak.raw[peak.raw$Peak==1,]
plot(peak1$Long,peak1$sst.rel.change,col=1,pch=19,cex=0.5)
peaklong <- peak1[order(peak1$Long),]
lines(loess.smooth(peaklong$Long,peaklong$sst.rel.change),lty=2)

# PLOT SST RELATIVE CHANGE IN MONTH OF PEAK SPAWNING BY MEAN ANNUAL SST
peak1 <- peak.raw[peak.raw$Peak==1,]
plot(peak1$mean.ann.sst,peak1$sst.rel.change,col=1,pch=19,cex=0.5)
peak.ann.sst <- peak1[order(peak1$mean.ann.sst),]
lines(loess.smooth(peak.ann.sst$mean.ann.sst,peak.ann.sst$sst.rel.change),lty=2)

# SST relative change per site
month.rel.change <- matrix(NA,length(unique(peak.raw$Site.number)),5) 
colnames(month.rel.change) <- c("-2","1","PEAK","+1","+2")
for(i in 2:length(unique(peak.raw$Site.number))){
   x <- subset(peak.raw,peak.raw$Site.number==unique(peak.raw$Site.number)[i])
   peak.row <- which(x$Peak==1)
   month.rel.change <- x$sst.rel.change[(peak.row-2):(peak.row+2)] 
}
# THIS DOESNT WORK BECAUSE OF SPAWNING AT END OF YEAR!
# cheated and did it manually

month.rel.change <- read.csv("peakANDsstrelchange-2+2.csv")
colnames(month.rel.change) <- c("Site.number","-2","1","PEAK","+1","+2")
matplot(t(month.rel.change[,2:6]),type="l",col=1,lty=1)
abline(v=3,lty=2)
boxplot(month.rel.change[,2:6],col=c("grey","grey","green","grey","grey"))

# repeat for SST relative average
month.rel.avg <- read.csv("peakANDsstrelavg-2+2.csv")
colnames(month.rel.avg) <- c("Site.number","-2","1","PEAK","+1","+2")
matplot(t(month.rel.avg[,2:6]),type="l",col=1,lty=1)
abline(v=3,lty=2)
boxplot(month.rel.avg[,2:6],col=c("grey","grey","green","grey","grey"))

# repeat for SST average
month.avg <- read.csv("peakANDsstavg-2+2.csv")
colnames(month.avg) <- c("Site.number","-2","1","PEAK","+1","+2")
matplot(t(month.avg[,2:6]),type="l",col=1,lty=1)
abline(v=3,lty=2)
boxplot(month.avg[,2:6],col=c("grey","grey","green","grey","grey"))




#################################################
## MANUAL PLOTS OF PARTIAL COEFFICIENTS
# 01/07/2015

pdf("PartialCoefMAvg.pdf")

# PEAK SPAWNING
######

par(mfrow=c(2,2))
coef95p <- coefTable(sa.fitp$averaged.mod,full=T)
cil <- coef95p[,1]-1.96*coef95p[,2]
ciu <- coef95p[,1]+1.96*coef95p[,2]
coef95p <- cbind(coef95p,cil,ciu)

xsrc <- seq(min(dp$sst.rel.change),max(dp$sst.rel.change),(max(dp$sst.rel.change)-min(dp$sst.rel.change))/335)
xc <- seq(min(dp$current.speed),max(dp$current.speed),(max(dp$current.speed)-min(dp$current.speed))/335)
# SST relative change response plot       
plot(xsrc,plogis(coef95p[1,1] + coef95p[2,1]*xsrc + coef95p[3,1] * mean(dp$current.speed)),lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(peak spawn)",xlab="SST relative change")       
lines(xsrc,plogis(coef95p[1,3] + coef95p[2,3]*xsrc + coef95p[3,3] * mean(dp$current.speed)),lwd=1,lty=3)          
lines(xsrc,plogis(coef95p[1,4] + coef95p[2,4]*xsrc + coef95p[3,4] * mean(dp$current.speed)),lwd=1,lty=3) 
# Current speed response plot
plot(xc,plogis(coef95p[1,1] + coef95p[2,1]*mean(dp$sst.rel.change) + coef95p[3,1] * xc),lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(peak spawn)",xlab="current speed")       
lines(xc,plogis(coef95p[1,3] + coef95p[2,3]*mean(dp$sst.rel.change) + coef95p[3,3] * xc),lwd=1,lty=3)          
lines(xc,plogis(coef95p[1,4] + coef95p[2,4]*mean(dp$sst.rel.change) + coef95p[3,4] * xc),lwd=1,lty=3) 



######
# SPAWNING
#####
par(mfrow=c(2,2))

coef95s <- coefTable(sa.fit.spawn$averaged.mod,full=T)
cils <- coef95s[,1]-1.96*coef95s[,2]
cius <- coef95s[,1]+1.96*coef95s[,2]
coef95s <- cbind(coef95s,cils,cius)

xsrcs <- seq(min(ds$sst.rel.change),max(ds$sst.rel.change),(max(ds$sst.rel.change)-min(ds$sst.rel.change))/335)
xsar <- seq(min(ds$sst.avg.rel),max(ds$sst.avg.rel),(max(ds$sst.avg.rel)-min(ds$sst.avg.rel))/335)
xw <- seq(min(ds$wind.avg),max(ds$wind.avg),(max(ds$wind.avg)-min(ds$wind.avg))/335)

# wind speed response plot       
plot(xw,plogis(coef95s[1,1] + coef95s[2,1]*xw + coef95s[3,1]*(xw^2)),
                 lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(spawn)",xlab="current speed")  
lines(xw,plogis(coef95s[1,3] + coef95s[2,3]*xw + coef95s[3,3]*(xw^2)),lty=3,type="l",ylim=c(0,1)) 
lines(xw,plogis(coef95s[1,4] + coef95s[2,4]*xw + coef95s[3,4]*(xw^2)),lty=3,type="l",ylim=c(0,1)) 

# SST average (relative) response plot       
plot(xsar,plogis(coef95s[1,1] + coef95s[5,1]*xsar),lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(spawn)",xlab="SST average (relative)")  
lines(xsar,plogis(coef95s[1,3] + coef95s[5,3]*xsar),lty=3,type="l",ylim=c(0,1)) 
lines(xsar,plogis(coef95s[1,4] + coef95s[5,4]*xsar),lty=3,type="l",ylim=c(0,1)) 

# SST relative change response plot       
plot(xsrcs,plogis(coef95s[1,1] +  coef95s[4,1]*xsrcs),lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(spawn)",xlab="SST average (relative)")  
lines(xsrcs,plogis(coef95s[1,3] +  coef95s[4,3]*xsrcs),lty=3,type="l",ylim=c(0,1)) 
lines(xsrcs,plogis(coef95s[1,4] + coef95s[4,4]*xsrcs),lty=3,type="l",ylim=c(0,1)) 



# #####
# # Wind speed x current V response plot
# windquant <- quantile(dp$wind.avg,probs=seq(0,1,1/5))
# par(mfrow=c(2,2))
# for(i in 2:(length(windquant)-1)){
#    plot(xcu,plogis(coef95s[1,1] + coef95s[2,1]*mean(ds$sst.rel.change) + coef95s[3,1]*windquant[i] + coef95s[4,1]*windquant[i]^2 + coef95s[5,1]*mean(ds$currentv.abs)+ 
#                       coef95s[6,1]*xcu*windquant[i] + coef95s[7,1]*xcu),lwd=1,lty=1,type="l",ylim=c(0,1),ylab="P(spawn)",xlab="current V",main=paste("wind = ",round(windquant[i],2),sep=""))  
#    lines(xcu,plogis(coef95s[1,2] + coef95s[2,2]*mean(ds$sst.rel.change) + coef95s[3,2]*windquant[i] + coef95s[4,2]*windquant[i]^2 + coef95s[5,2]*mean(ds$currentv.abs)+ 
#                        coef95s[6,2]*xcu*windquant[i] + coef95s[7,2]*xcu),lwd=1,lty=3,type="l",ylim=c(0,1)) 
#    lines(xcu,plogis(coef95s[1,3] + coef95s[2,3]*mean(ds$sst.rel.change) + coef95s[3,3]*windquant[i] + coef95s[4,3]*windquant[i]^2 + coef95s[5,3]*mean(ds$currentv.abs)+ 
#                        coef95s[6,3]*xcu*windquant[i] + coef95s[7,3]*xcu),lwd=1,lty=3,type="l",ylim=c(0,1)) 
# }
# 
# dev.off()

##############



