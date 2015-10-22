# SK 14/02/2014

## MIXED EFFECTS MODEL WITH SIGNIFICANT INTERACTIONS & QUADRATICS

# OUTPUTS
# Regression summary printed, dotplots of random effects, 
# variance partitioning coefficients

mixed.mod.output <- function(formula,error.distribution,dataset,random){
   
   require(lme4)
   fit <- glmer(formula,family=error.distribution,data=dataset,na.action=na.fail)
   print(summary(fit))
   dp <- dotplot(ranef(fit,condVar=TRUE,whichel=random[1]))
   print(dp)
   if(length(random)>1){      
       dp2 <- dotplot(ranef(fit,condVar=TRUE,whichel=random[2]))
       print(dp2)
   }
   
   # VARIANCE PARTITION COEFFICIENT 
   # VPC = variance for random effect/(variance for random effect + 3.29)
   # The 3.29 is because it uses a binomial distribution
   random.vars <- as.data.frame(VarCorr(fit)) 
   month.vpc <- random.vars[which(random.vars[,1]==random[1]),4]
   month.vpc/(month.vpc+3.29)
   site.vpc <- random.vars[which(random.vars[,1]==random[2]),4]
   site.vpc/(site.vpc+3.29)
   print(paste(random[1]," VPC = ",round(month.vpc,3),
               ". ",random[1]," accounts for ",round(month.vpc,5)*100," of the variance in the data",sep=""))
   print(paste(random[2]," VPC = ",round(site.vpc,3),
               ". ",random[2]," accounts for ",round(site.vpc,5)*100," of the variance in the data",sep=""))
}
