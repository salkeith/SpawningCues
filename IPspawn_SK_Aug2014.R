# SK 02/08/2014

###########################################
###########################################

# A Baird Coral Spawning Cues Indo-Pacific

###########################################
###########################################

# Selected variables based on hypotheses

library(MuMIn)
library(lme4)
library(visreg)
library(effects)
load("Normalized data used for models IP.RData")

################################################################################

### 1. SST MODEL

#################################################

# no DL because correlated with PAR_avg
sst <- glm(Spawning ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
                       SST_avg*Precip_avg+SST_avg*SST_deriv+SST_avg*diff_ss+
                       SST_deriv*Precip_avg+SST_deriv*diff_ss+Precip_avg*diff_ss+Month+Site.number,
                       family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sst)
# remove non-significant interaction terms                       
sst <- glm(Spawning ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
                       SST_deriv*Precip_avg+Precip_avg*diff_ss+Month+Site.number,
                       family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sst)
# Month and site not significant - mixed model not required 
sst <- glm(Spawning ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
                       SST_deriv*Precip_avg+Precip_avg*diff_ss,
                       family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sst)                       

# check GLM assumptions
hist(resid(sst))
par(mfcol=c(2,2))
plot(sst)
# Test for overdispersion 
X2 <- sum(residuals(sst,type="pearson")^2)
phi <- X2/sst$df.residual; phi    # if >1 suggests over-dispersion. Want it to be just below 1
# if output below significant, indicates over-dispersion  
p.val <- pchisq(phi * sst$df.residual,sst$df.residual, lower = F); p.val        
# if output below significant, the model doesn't fit well
pp <- sum(resid(sst, type="pearson")^2);  print(1-pchisq(pp,sst$df.resid))

# plot response (shows change in probability of spawning for each variable 
# within the model, with 95% confidence intervals)
par(mfcol=c(3,2))
visreg(sst,partial=F,scale="response",ylim=c(0,1))

# plot significant interactions
visreg(sst,"diff_ss",by="Precip_avg",scale="response",breaks=3)
visreg(sst,"SST_deriv",by="Precip_avg",scale="response",breaks=3)


# model selection and avergaing
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects
msubset <- expression(c(Wind_avg|!`I(Wind_avg^2)`),(diff_ss|!Precip_avg*diff_ss),
                      (Precip_avg|!Precip_avg*diff_ss),(Precip_avg|!Precip_avg*SST_deriv),
                      (SST_deriv|!Precip_avg*SST_deriv))
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(sst,subset=msubset)
sst.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
sst.sel


#######################################
## REMOVE ALL INTERACTIONS
sst.no.int <- glm(Spawning ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss,
         family = binomial(link="logit"), data = d, na.action = na.fail)      
# model selection and avergaing
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects
msubset <- expression(c(Wind_avg|!`I(Wind_avg^2)`))
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(sst.no.int,subset=msubset)
sst.no.int.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
sst.no.int.sel

sst.no.int.best <- glm(Spawning ~ SST_avg+SST_deriv,family = binomial(link="logit"), data = d, na.action = na.fail)      

#######################################

sst.best <- glm(Spawning ~ Precip_avg+SST_avg+SST_deriv+diff_ss+Precip_avg*diff_ss+diff_ss*Precip_avg,
                family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sst.best)
sst.best.no.int <- glm(Spawning ~ SST_avg+SST_deriv,family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sst.best.no.int)

# model-averaged confidence intervals
cisstb <- confint(sst.best)
# use FULL model-averaged coefficients because not all models in the set contain
# both interaction terms and main effects, or quadratic wind term
cosstb <- coef(sst.best,full=T)
CIsstb <- cbind(cosstb,cisstb) 
CIsstb
z <- allEffects(sst.best,transformation=list(link=logit,
                inverse=family(sst.best)$linkinv)
plot(exp(z),rug=F)
par(mfcol=c(2,2))
visreg(sst.best,scale="response",ylim=c(0,1))

# model-averaged confidence intervals
cisst <- confint(sst.sel)
# use FULL model-averaged coefficients because not all models in the set contain
# both interaction terms and main effects, or quadratic wind term
cosst <- coef(sst.sel,full=T)
CIsst <- cbind(cosst,cisst) 
CIsst
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. not sure how to treat the interactions here...
CIsst[!(CIsst[,2]<0 & CIsst[,3]>0),]

##################################################
# RESPONSE PLOTS for averaged model 
# manual plotting required, visreg doesn't work on model-averged object
pdf("SpawningSSTModelPartialCoefs.pdf")
par(mfrow=c(2,2),oma=c(0,1,0,0))

# SST AVERAGE
x <- seq(min(d$SST_avg),max(d$SST_avg),length=100)
sst.pred <- with(d,(predict(sst.sel,type="response",se=T,newdata=data.frame(SST_avg=x,
            SST_deriv=mean(SST_deriv),diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~SST_avg,data=d,col="lightgray",pch=16,ylab="P(Spawning)",
       xlab="SST average",cex.lab=1.5)
points(sst.pred$fit~x,type="l",col=1,lwd=2)
lines(sst.pred$fit+sst.pred$se.fit~x,lty=2,type="l")
lines(sst.pred$fit-sst.pred$se.fit~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
sst.pred <- with(d,(predict(sst.sel,type="response",se=T,newdata=data.frame(SST_deriv=x,
            SST_avg=mean(SST_avg),diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~SST_deriv,data=d,col="lightgray",pch=16,ylab="P(Spawning)",
       xlab="SST derivative",cex.lab=1.5)
points(sst.pred$fit~x,type="l",col=1,lwd=2)
lines(sst.pred$fit+sst.pred$se.fit~x,lty=2,type="l")
lines(sst.pred$fit-sst.pred$se.fit~x,lty=2,type="l")

# plot interaction 
# plot diff_ss as binned variable, Precip_avg as continuous
# DIFF_SS * PRECIP_AVG
par(mfcol=c(2,2))
for(i in -2:1){
   diff_ssint <- i
   x <- seq(min(d$Precip_avg),max(d$Precip_avg),length=100)
   precipint.pred <- with(d,(predict(sst.sel,type="response",se=T,newdata=data.frame(Precip_avg=x,
                 diff_ss=diff_ssint,SST_avg=mean(SST_avg),Wind_avg=mean(Wind_avg),
                 SST_deriv=mean(SST_deriv)))))
   plot(Spawning~Precip_avg,data=d,col="lightgray",pch=16,ylab="P(Spawning)",
                  xlab="Mean Precipitation",cex.lab=1.5,
                  main=paste("Difference in Sunset Mean ",diff_ssint,sep=""))
   points(precipint.pred$fit~x,type="l",col=1,lwd=2)
   lines(precipint.pred$fit+precipint.pred$se.fit~x,lty=2,type="l")
   lines(precipint.pred$fit-precipint.pred$se.fit~x,lty=2,type="l")
}


#################################
# RESPONSE PLOTS - BEST MODEL
# manual plotting required, visreg doesn't work on model-averged object
pdf("SpawningCoefsBestModelNoInteractions.pdf")
par(mfcol=c(2,2)) 

# SST Average
x <- seq(min(d$SST_avg),max(d$SST_avg),length=100)
ssta.pred <- with(d,(predict(sst.best.no.int,type="response",se.fit=T,newdata=data.frame(SST_avg=x,
                     SST_deriv=mean(SST_deriv)))))
ssta.pred2 <- with(par.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~SST_avg,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="SST mean",cex.lab=1.5)
points(ssta.pred2[,1]~x,type="l",col=1,lwd=2)
lines(ssta.pred2[,2]~x,lty=2,type="l")
lines(ssta.pred2[,3]~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
sstd.pred <- with(d,(predict(sst.best.no.int,type="response",se.fit=T,newdata=data.frame(SST_deriv=x,
                    SST_avg=mean(SST_avg)))))
sstd.pred2 <- with(sstd.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~SST_deriv,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="SST derivative",cex.lab=1.5)
points(sstd.pred2[,1]~x,type="l",col=1,lwd=2)
lines(sstd.pred2[,2]~x,lty=2,type="l")
lines(sstd.pred2[,3]~x,lty=2,type="l")

dev.off()

################################################################################

### 2. PAR MODEL

############################################

light <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon+
                         PAR_avg*PAR_10.mon+PAR_avg*Precip_avg+PAR_10.mon*Precip_avg+
                         PAR_avg*diff_ss+PAR_10.mon*diff_ss+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(light)
# remove non-significant interaction terms 
light <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(light)                         
# Month and site not significant - mixed model not required  
light <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(light)
hist(resid(light))                         
windows()
par(mfcol=c(2,2))
plot(light) 
X2 <- sum(residuals(light,type="pearson")^2)
phi <- X2/light$df.residual; phi    
p.val <- pchisq(phi * light$df.residual,light$df.residual, lower = F); p.val          
pp <- sum(resid(light, type="pearson")^2); print(1-pchisq(pp,light$df.resid))  

# plot response
windows()
par(mfcol=c(3,2))
visreg(light,partial=F,scale="response",ylim=c(0,1)) 

# model selection and avergaing
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects
msubset <- expression(Wind_avg|!`I(Wind_avg^2)`)
                      
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(light,subset=msubset)
light.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
light.sel

# model-averaged confidence intervals 
cilight <- confint(light.sel)
colight <- coef(light.sel)
CIlight <- cbind(colight,cilight) 
CIlight
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. not sure how to treat the interactions here...
CIlight[!(CIlight[,2]<0 & CIlight[,3]>0),]

# RESPONSE PLOTS  
# manual plotting required, visreg doesn't work on model-averged object
par(mfcol=c(2,1))
# PAR AVERAGE
x <- seq(min(d$PAR_avg),max(d$PAR_avg),length=100)
par.pred <- with(d,(predict(light.sel,type="response",se=T,newdata=data.frame(PAR_avg=x,
            PAR_10.mon=mean(PAR_10.mon),diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~PAR_avg,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="PAR average",cex.lab=1.5)
points(par.pred$fit~x,type="l",col=1,lwd=2)
lines(par.pred$fit+par.pred$se.fit~x,lty=2,type="l")
lines(par.pred$fit-par.pred$se.fit~x,lty=2,type="l")
# WIND AVERAGE
x <- seq(min(d$Wind_avg),max(d$Wind_avg),length=100)
par.pred <- with(d,(predict(light.sel,type="response",se=T,newdata=data.frame(Wind_avg=x,
            PAR_avg=mean(PAR_avg),diff_ss=mean(diff_ss),PAR_10.mon=mean(PAR_10.mon),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~Wind_avg,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="Wind average",cex.lab=1.5)
points(par.pred$fit~x,type="l",col=1,lwd=2)
lines(par.pred$fit+par.pred$se.fit~x,lty=2,type="l")
lines(par.pred$fit-par.pred$se.fit~x,lty=2,type="l")



#############
# PAR & SST
#############

sstlight <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                         SST_deriv*PAR_10.mon+SST_deriv*Precip_avg+PAR_10.mon*Precip_avg+
                         SST_deriv*diff_ss+PAR_10.mon*diff_ss+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sstlight)
# remove non-significant interaction terms
sstlight <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                         Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sstlight)
# Month and site not significant, no mixed model required
sstlight <- glm(Spawning ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv,
                     family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sstlight)
hist(resid(sstlight))                         
par(mfcol=c(2,2))
plot(sstlight) 
X2 <- sum(residuals(sstlight,type="pearson")^2)
phi <- X2/sstlight$df.residual; phi    
p.val <- pchisq(phi * sstlight$df.residual,sstlight$df.residual, lower = F); p.val          
pp <- sum(resid(sstlight, type="pearson")^2); print(1-pchisq(pp,sstlight$df.resid)) 

# plot response
par(mfcol=c(3,2))
visreg(sstlight,partial=F,scale="response",ylim=c(0,1)) 

# model selection and avergaing
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects
msubset <- expression(Wind_avg|!`I(Wind_avg^2)`)
                      
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(sstlight,subset=msubset)
sstlight.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
sstlight.sel

# model-averaged confidence intervals 
cisstlight <- confint(sstlight.sel)
cosstlight <- coef(sstlight.sel,full=T)
CIsstlight <- cbind(cosstlight,cisstlight) 
CIsstlight
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. not sure how to treat the interactions here...
CIsstlight[!(CIsstlight[,2]<0 & CIsstlight[,3]>0),]

# best model
sstlight.best <- glm(Spawning ~ Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv,
                     family = binomial(link="logit"), data = d, na.action = na.fail)
cisstlight <- confint(sstlight.best)
cosstlight <- coef(sstlight.best,full=T)
CIsstlight.best <- cbind(cosstlight,cisstlight) 
CIsstlight.best

####################################
# RESPONSE PLOTS - MODEL AVERAGED
# manual plotting required, visreg doesn't work on model-averged object
par(mfcol=c(2,2)) 
# PAR AVERAGE
x <- seq(min(d$PAR_10.mon),max(d$PAR_10.mon),length=100)
parsst.pred <- with(d,(predict(sstlight.sel,type="response",interval="confidence",se=T,newdata=data.frame(PAR_10.mon=x,
            diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),SST_deriv=mean(SST_deriv),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~PAR_10.mon,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="PAR 10 month integral",cex.lab=1.5)
points(parsst.pred$fit~x,type="l",col=1,lwd=2)
lines(parsst.pred$fit+parsst.pred$se.fit~x,lty=2,type="l")
lines(parsst.pred$fit-parsst.pred$se.fit~x,lty=2,type="l")
# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
parsst.pred <- with(d,(predict(sstlight.sel,type="response",se=T,newdata=data.frame(PAR_10.mon=mean(PAR_10.mon),
            diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),SST_deriv=x,
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~SST_deriv,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="SST derivative",cex.lab=1.5)
points(parsst.pred$fit~x,type="l",col=1,lwd=2)
lines(parsst.pred$fit+parsst.pred$se.fit~x,lty=2,type="l")
lines(parsst.pred$fit-parsst.pred$se.fit~x,lty=2,type="l")
# DIFFERENCE SUNSET
x <- seq(min(d$diff_ss),max(d$diff_ss),length=100)
parsst.pred <- with(d,(predict(sstlight.sel,type="response",se=T,newdata=data.frame(PAR_10.mon=mean(PAR_10.mon),
            diff_ss=x,Wind_avg=mean(Wind_avg),SST_deriv=mean(SST_deriv),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~diff_ss,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="Difference in Sunset",cex.lab=1.5)
points(parsst.pred$fit~x,type="l",col=1,lwd=2)
lines(parsst.pred$fit+parsst.pred$se.fit~x,lty=2,type="l")
lines(parsst.pred$fit-parsst.pred$se.fit~x,lty=2,type="l")
# WIND AVERAGE
x <- seq(min(d$Wind_avg),max(d$Wind_avg),length=100)
parsst.pred <- with(d,(predict(sstlight.sel,type="response",se=T,newdata=data.frame(PAR_10.mon=mean(PAR_10.mon),
            diff_ss=mean(diff_ss),Wind_avg=x,SST_deriv=mean(SST_deriv),
            Precip_avg=mean(Precip_avg)))))
plot(Spawning~Wind_avg,data=d,col="lightgray",pch=16,ylab="Spawning probability",
       xlab="Wind average",cex.lab=1.5)
points(parsst.pred$fit~x,type="l",col=1,lwd=2)
lines(parsst.pred$fit+parsst.pred$se.fit~x,lty=2,type="l")
lines(parsst.pred$fit-parsst.pred$se.fit~x,lty=2,type="l")

#################################
# RESPONSE PLOTS - BEST MODEL
# manual plotting required, visreg doesn't work on model-averged object
par(mfcol=c(2,2)) 

# PAR AVERAGE
x <- seq(min(d$PAR_10.mon),max(d$PAR_10.mon),length=100)
par.pred <- with(d,(predict(sstlight.best,type="response",se.fit=T,newdata=data.frame(PAR_10.mon=x,
                        diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),SST_deriv=mean(SST_deriv)))))
par.pred2 <- with(par.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~PAR_10.mon,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="PAR 10 month integral",cex.lab=1.5)
points(par.pred2[,1]~x,type="l",col=1,lwd=2)
lines(par.pred2[,2]~x,lty=2,type="l")
lines(par.pred2[,3]~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
sstd.pred <- with(d,(predict(sstlight.best,type="response",se.fit=T,newdata=data.frame(SST_deriv=x,
                     diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),PAR_10.mon=mean(PAR_10.mon)))))
sstd.pred2 <- with(sstd.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~SST_deriv,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="SST derivative",cex.lab=1.5)
points(sstd.pred2[,1]~x,type="l",col=1,lwd=2)
lines(sstd.pred2[,2]~x,lty=2,type="l")
lines(sstd.pred2[,3]~x,lty=2,type="l")

# DIFFERENCE SUNSET
x <- seq(min(d$diff_ss),max(d$diff_ss),length=100)
ss.pred <- with(d,(predict(sstlight.best,type="response",se.fit=T,newdata=data.frame(diff_ss=x,
                  PAR_10.mon=mean(PAR_10.mon),Wind_avg=mean(Wind_avg),SST_deriv=mean(SST_deriv)))))
ss.pred2 <- with(ss.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~diff_ss,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="Difference in sunset",cex.lab=1.5)
points(ss.pred2[,1]~x,type="l",col=1,lwd=2)
lines(ss.pred2[,2]~x,lty=2,type="l")
lines(ss.pred2[,3]~x,lty=2,type="l")

# WIND AVERAGE
x <- seq(min(d$Wind_avg),max(d$Wind_avg),length=100)
w.pred <- with(d,(predict(sstlight.best,type="response",se.fit=T,newdata=data.frame(Wind_avg=x,
                      diff_ss=mean(diff_ss),PAR_10.mon=mean(PAR_10.mon),SST_deriv=mean(SST_deriv)))))
w.pred2 <- with(w.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Spawning~Wind_avg,data=d,col="lightgray",pch=16,ylab="Spawning probability",
     xlab="Averaege wind speed",cex.lab=1.5)
points(w.pred2[,1]~x,type="l",col=1,lwd=2)
lines(w.pred2[,2]~x,lty=2,type="l")
lines(w.pred2[,3]~x,lty=2,type="l")

### Which is the best model out of the three peak spawning mixed models?
AICc(light,sst,sstlight,sstlight.best)


################################################################################

### 3. ALL THREE PREVIOUS MODELS WITH "PEAK" AS DEPENDENT VARIABLE

################################################################

########
# SST
########

psst <- glm(Peak ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
                       SST_avg*Precip_avg+SST_avg*SST_deriv+SST_avg*diff_ss+
                       SST_deriv*Precip_avg+SST_deriv*diff_ss+Precip_avg*diff_ss+Month+Site.number,
                       family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psst)
# remove non-significant interaction terms
psst <- glm(Peak ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
                       SST_deriv*diff_ss+Month+Site.number,
                       family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psst)
# Month significant - mixed model required                      
# check assumptions as GLM first                  
par(mfcol=c(2,2))
plot(psst)
# Test for overdispersion 
X2 <- sum(residuals(psst,type="pearson")^2)
phi <- X2/psst$df.residual; phi    
p.val <- pchisq(phi * psst$df.residual,psst$df.residual, lower = F);  p.val         
pp <- sum(resid(psst, type="pearson")^2); print(1-pchisq(pp,psst$df.resid))  

# mixed model
psstm <- glmer(Peak ~ Wind_avg+I(Wind_avg^2)+Precip_avg+SST_avg+SST_deriv+diff_ss+
            SST_deriv*diff_ss+(1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstm)

# check whether intercepts vary a lot from the mean for each site/month
dotplot(ranef(psstm,condVar=TRUE,whichel="Month"))

# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
month <- 0.63/(0.63+3.29)
month
# VPC is 0.16 so 16% of the variance is accounted for by month

# is there a big improvment with the addition of random effects?
# (i.e., does the AICc drop > 3 AIC?)
AICc(psst,psstm)

# model selection and avergaing
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects
msubset <- expression(c(Wind_avg|!`I(Wind_avg^2)`),(diff_ss|!SST_deriv*diff_ss),
                      (SST_deriv|!SST_deriv*diff_ss))
# model selection by delta < 3 (Bolker 2009) with polynomials
dmod <- dredge(psstm,subset=msubset)
psstm.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
psstm.sel

#################################
# RESPONSE PLOTS - BEST MODEL
# manual plotting required, visreg doesn't work on model-averged object
pdf("PeakSSTonlyPartialCoefs.pdf")
par(mfrow=c(2,2),oma=c(0,2,0,0)) 

# SST AVERAGE
x <- seq(min(d$SST_avg),max(d$SST_avg),length=100)
ssta.pred <- with(d,(predict(psstm.sel,type="response",se.fit=T,newdata=data.frame(SST_avg=x,
                     diff_ss=mean(diff_ss),SST_deriv=mean(SST_deriv),Precip_avg=mean(Precip_avg),
                     Wind_avg=mean(Wind_avg)))))
ssta.pred2 <- with(ssta.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Peak~SST_avg,data=d,col="lightgray",pch=16,ylab="P(Peak spawning)",
     xlab="Mean SST",cex.lab=1.5)
points(ssta.pred2[,1]~x,type="l",col=1,lwd=2)
lines(ssta.pred2[,2]~x,lty=2,type="l")
lines(ssta.pred2[,3]~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
sstd.pred <- with(d,(predict(psstm.sel,type="response",se.fit=T,newdata=data.frame(SST_deriv=x,
                     diff_ss=mean(diff_ss),SST_avg=mean(SST_avg),Precip_avg=mean(Precip_avg),
                     Wind_avg=mean(Wind_avg)))))
sstd.pred2 <- with(sstd.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Peak~SST_deriv,data=d,col="lightgray",pch=16,ylab="P(Peak spawning)",
     xlab="SST derivative",cex.lab=1.5)
points(sstd.pred2[,1]~x,type="l",col=1,lwd=2)
lines(sstd.pred2[,2]~x,lty=2,type="l")
lines(sstd.pred2[,3]~x,lty=2,type="l")

dev.off()


########
# PAR
########

plight <- glm(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon+
                         PAR_avg*PAR_10.mon+PAR_avg*Precip_avg+PAR_10.mon*Precip_avg
                         +PAR_avg*diff_ss+PAR_10.mon*diff_ss+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(plight)
# remove non-significant interaction terms   
# month significant and model under-dispersed so require mixed model
plight <- glm(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon,
                         family = binomial(link="logit"), data = d, na.action = na.fail)
# check assumptions. Will need mixed model because Month is significant
# but assumptions must be checked as GLM first                       
hist(resid(plight))                         
par(mfcol=c(2,2))
plot(plight)
X2 <- sum(residuals(plight,type="pearson")^2)
phi <- X2/plight$df.residual
phi    # if >1 suggests over-dispersion. Want it to be just below 1
p.val <- pchisq(phi * plight$df.residual,plight$df.residual, lower = F)
p.val       # if significant, indicates over-dispersion   
pp <- sum(resid(plight, type="pearson")^2)
print(1-pchisq(pp,plight$df.resid))  # if significant, the model doesn't fit well

# month significant and model under-dispersed so require mixed model
plightm <- glmer(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_avg+PAR_10.mon+
                         (1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)

summary(plightm)

# check whether intercepts vary a lot from the mean for each site/month
dotplot(ranef(plightm,condVar=TRUE,whichel="Month"))

# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
month <- 1.24/(1.24+3.29)
month
# VPC is 0.27 so 27% of the variance is accounted for by month

# is there a big improvment with the addition of random effects?
# (i.e., does the AICc drop > 3 AIC?)
AICc(plight,plightm)
                 
# model selection by delta < 3 (Bolker 2009) with polynomials
msubset <- expression(Wind_avg|!`I(Wind_avg^2)`)   
dmod <- dredge(plightm,subset=msubset)
plightm.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
plightm.sel



##############
# PAR & SST
##############

psstlight <- glm(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                         SST_deriv*PAR_10.mon+SST_deriv*Precip_avg+PAR_10.mon*Precip_avg+
                         SST_deriv*diff_ss+PAR_10.mon*diff_ss+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstlight)
# remove non-significant interaction terms
psstlight <- glm(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                         PAR_10.mon*Precip_avg+SST_deriv*diff_ss+Month+Site.number, 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstlight)
# Month marginally significant, check mixed model
psstlightm <- glmer(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                         PAR_10.mon*Precip_avg+SST_deriv*diff_ss+(1|Month), 
                         family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstlightm)
# Month as a random effect accounts for a fair amount of variation so mixed model appropriate

# check whether intercepts vary a lot from the mean for each site/month
dotplot(ranef(psstlightm,condVar=TRUE,whichel="Month"))

# Variance partition coefficient 
# VPC = variance for random effect/(variance for random effect + 3.29)
# The 3.29 is because it uses a binomial distribution
month <- 1.08/(1.08+3.29); month
# VPC is 0.25 so 25% of the variance is accounted for by month

# is there a big improvment with the addition of random effects?
# (i.e., does the AICc drop > 3 AIC?)
AICc(psstlight,psstlightm)
                 
# model selection by delta < 3 (Bolker 2009) with polynomials
msubset <- expression(c(Wind_avg|!`I(Wind_avg^2)`),(diff_ss|!diff_ss*SST_deriv),
                      (SST_deriv|!diff_ss*SST_deriv),(Precip_avg|!Precip_avg*PAR_10.mon),
                      (PAR_10.mon|!Precip_avg*PAR_10.mon))   
dmod <- dredge(psstlightm,subset=msubset)
psstlightm.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
psstlightm.sel


# model-averaged confidence intervals 
cipsstlightm <- confint(psstlightm.sel)
copsstlightm <- coef(psstlightm.sel,full=T)
CIpsstlightm <- cbind(copsstlightm,cipsstlightm) 
CIpsstlightm
# return variables that do not have coefficients that overlap zero
# (i.e., those that are significant)
# N.B. not sure how to treat the interactions here...
CIpsstlightm[!(CIpsstlightm[,2]<0 & CIpsstlightm[,3]>0),]
# Precip_avg is significant but also part of a significant interaction
# so view as interaction only. Main effect meaningless.

#######
# RESPONSE PLOTS for significant variables 
# manual plotting required, visreg doesn't work on model-averged object
pdf("PeakSpawningSSTPARModelPartialCoefs.pdf")
par(mfrow=c(2,2),oma=c(0,1,0,0))

# DIFF_SS
x <- seq(min(d$diff_ss),max(d$diff_ss),length=100)
diffss.pred <- with(d,(predict(psstlightm.sel,type="response",se=T,newdata=data.frame(diff_ss=x,
               SST_deriv=mean(SST_deriv),PAR_10.mon=mean(PAR_10.mon),Wind_avg=mean(Wind_avg),
               Precip_avg=mean(Precip_avg)))))
plot(Spawning~diff_ss,data=d,col="lightgray",pch=16,ylab="P(Peak Spawning)",
     xlab="Difference in Sunset",cex.lab=1.5)
points(diffss.pred$fit~x,type="l",col=1,lwd=2)
lines(diffss.pred$fit+diffss.pred$se.fit~x,lty=2,type="l")
lines(diffss.pred$fit-diffss.pred$se.fit~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
SST_deriv.pred <- with(d,(predict(psstlightm.sel,type="response",se=T,newdata=data.frame(SST_deriv=x,
               diff_ss=mean(diff_ss),PAR_10.mon=mean(PAR_10.mon),Wind_avg=mean(Wind_avg),
               Precip_avg=mean(Precip_avg)))))
plot(Spawning~SST_deriv,data=d,col="lightgray",pch=16,ylab="P(Peak Spawning)",
     xlab="SST Derivative",cex.lab=1.5)
points(SST_deriv.pred$fit~x,type="l",col=1,lwd=2)
lines(SST_deriv.pred$fit+SST_deriv.pred$se.fit~x,lty=2,type="l")
lines(SST_deriv.pred$fit-SST_deriv.pred$se.fit~x,lty=2,type="l")

# plot interaction 
# plot diff_ss as binned variable, Precip_avg as continuous
# PAR_10.MON * PRECIP_AVG
par(mfrow=c(2,2),oma=c(0,1,0,0))
for(i in -2:1){
   PAR10monint <- i+0.5
   x <- seq(min(d$Precip_avg),max(d$Precip_avg),length=100)
   precipint.pred <- with(d,(predict(psstlightm.sel,type="response",se=T,newdata=data.frame(Precip_avg=x,
                           PAR_10.mon=PAR10monint,diff_ss=mean(diff_ss),Wind_avg=mean(Wind_avg),
                           SST_deriv=mean(SST_deriv)))))
   plot(Spawning~Precip_avg,data=d,col="lightgray",pch=16,ylab="P(Spawning)",
        xlab="Mean Precipitation",cex.lab=1.5,
        main=paste("PAR 10 Months Mean ", PAR10monint,sep=""))
   points(precipint.pred$fit~x,type="l",col=1,lwd=2)
   lines(precipint.pred$fit+precipint.pred$se.fit~x,lty=2,type="l")
   lines(precipint.pred$fit-precipint.pred$se.fit~x,lty=2,type="l")
}

dev.off()

### Which is the best model out of the three peak spawning mixed models?
AICc(plightm,psst,psstlightm)


####################################
## REMOVE INTERACTIONS
psstlightm.no.int <- glmer(Peak ~ Precip_avg+Wind_avg+I(Wind_avg^2)+diff_ss+PAR_10.mon+SST_deriv+
                              (1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstlightm.no.int)
# model selection by delta < 3 (Bolker 2009) with polynomials
msubset <- expression(c(Wind_avg|!`I(Wind_avg^2)`))   
dmod <- dredge(psstlightm.no.int,subset=msubset)
psstlightm.no.int.sel <- summary(model.avg(get.models(subset(dmod,delta < 3))))
psstlightm.no.int.sel

psstlightm.no.int.best <- glmer(Peak ~ Precip_avg+diff_ss+SST_deriv+
                              (1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(psstlightm.no.int.best)

#################################
# RESPONSE PLOTS - BEST MODEL
# manual plotting required, visreg doesn't work on model-averged object
par(mfcol=c(2,2)) 

# PRECIPITATION AVERAGE
x <- seq(min(d$Precip_avg),max(d$Precip_avg),length=100)
precip.pred <- with(d,(predict(psstlightm.no.int.best,type="response",se.fit=T,newdata=data.frame(Precip_avg=x,
                 diff_ss=mean(diff_ss),SST_deriv=mean(SST_deriv)))))
precip.pred2 <- with(precip.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Peak~Precip_avg,data=d,col="lightgray",pch=16,ylab="Peak spawning probability",
     xlab="Mean precipitation",cex.lab=1.5)
points(precip.pred2[,1]~x,type="l",col=1,lwd=2)
lines(precip.pred2[,2]~x,lty=2,type="l")
lines(precip.pred2[,3]~x,lty=2,type="l")

# SST DERIVATIVE
x <- seq(min(d$SST_deriv),max(d$SST_deriv),length=100)
sstd.pred <- with(d,(predict(psstlightm.no.int.best,type="response",se.fit=T,newdata=data.frame(SST_deriv=x,
                  diff_ss=mean(diff_ss)))))
sstd.pred2 <- with(sstd.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Peak~SST_deriv,data=d,col="lightgray",pch=16,ylab="Peak spawning probability",
     xlab="SST derivative",cex.lab=1.5)
points(sstd.pred2[,1]~x,type="l",col=1,lwd=2)
lines(sstd.pred2[,2]~x,lty=2,type="l")
lines(sstd.pred2[,3]~x,lty=2,type="l")

# DIFFERENCE SUNSET
x <- seq(min(d$diff_ss),max(d$diff_ss),length=100)
ss.pred <- with(d,(predict(psstlightm.no.int.best,type="response",se.fit=T,newdata=data.frame(diff_ss=x,
                  SST_deriv=mean(SST_deriv)))))
ss.pred2 <- with(ss.pred,cbind(fit,low=fit-1.96*se.fit,up=fit+1.96*se.fit))
plot(Peak~diff_ss,data=d,col="lightgray",pch=16,ylab="Peak spawning probability",
     xlab="Difference in sunset",cex.lab=1.5)
points(ss.pred2[,1]~x,type="l",col=1,lwd=2)
lines(ss.pred2[,2]~x,lty=2,type="l")
lines(ss.pred2[,3]~x,lty=2,type="l")





#########################################

## PEAK SPAWNING WITH PAR 10 MONTH MEAN ONLY

par10 <- glmer(Peak ~ PAR_10.mon+(1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(par10)
sstd <- glmer(Peak ~ SST_deriv+(1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(sstd)
precip <- glmer(Peak ~ Precip_avg+(1|Month),family = binomial(link="logit"), data = d, na.action = na.fail)
summary(precip)



