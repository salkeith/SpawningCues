# SK 14/02/2014

## DIAGNOSTICS FOR GLM

# OUTPUTS
# Regression diagnostic plots, VIF, histogram of residuals,
# over-dispersion, goodness-of-fit

glm.diagnostic <- function(formula,error.distribution,dataset){
   
   require("car")
   fit <- glm(formula,family=error.distribution,data=dataset,na.action=na.fail)
   summary(fit)
   par(mfcol=c(2,2))
   plot(fit)
   hist(resid(fit))
   # test for collinearity amongst predictors
   # Logan book suggests VIF >5 as a cut-off for colinearity problems
   # and tolerance (inverse of VIF) of <0.2
   print("VIF")
   print(vif(fit))
   
   # Test for overdispersion 
   X2 <- sum(residuals(fit,type="pearson")^2)
   phi <- X2/fit$df.residual
   if(phi>1) print(paste("Model is overdispersed. Phi =",round(phi,3))) else
      print(paste("Model is NOT overdispersed. Phi =",round(phi,3)))
   # if output below significant, indicates over-dispersion  
   p.val <- pchisq(phi * fit$df.residual,fit$df.residual, lower = F)
   if(p.val<0.05) print(paste("Model is overdispersed. p.val =",round(p.val,3))) else
      print(paste("Model is NOT overdispersed. p.val =",round(p.val,3)))
   # if output below significant, the model doesn't fit well
   pp <- sum(resid(fit, type="pearson")^2)
   gof <- 1-pchisq(pp,fit$df.resid)
   if(gof<0.05) print(paste("Model fit is poor. pp =",round(gof,3))) else
      print(paste("Model fit is OK. pp =",round(gof,3)))
}

