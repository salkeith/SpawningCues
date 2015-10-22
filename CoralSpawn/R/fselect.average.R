# SK 14/02/2014

## MODEL SELECTION AND AVERAGING
# expression so that quadratic term cannot be included unless linear term present
# and interaction terms cannot be included without the main effects.

# OUTPUTS
# Dredged model list, averaged model 

select.average <- function(mod.type,formula,error.distribution,dataset,msubset=NULL){
   
   require(lme4)
   if(mod.type=="mixed") 
      fit <- glmer(formula,family=error.distribution,data=dataset,na.action=na.fail)
   else
      fit <- glm(formula,family=error.distribution,data=dataset,na.action=na.fail)
   print("Patience... the next step can be slow...")
   
   # model selection by delta < 3 (Bolker 2009) with polynomials
   dmod <- dredge(fit,rank="BIC")
   fit.avg <- summary(model.avg(dmod,revised.var=TRUE,subset=delta<3))
   print(fit.avg) 
   
   return(list("dredge.res"=dmod,"averaged.model"=fit.avg))
   
}