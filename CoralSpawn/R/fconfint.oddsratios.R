# SK 14/02/2014

## CONFIDENCE INTERVALS AND ODDS RATIOS
# Get 95% confidence intervals and odds ratios for coefficients.

# OUTPUTS
# 95% confidence intervals, odds ratios

CI.OR <- function(fit.avg){
   
   ## MODEL-AVERAGED CONFIDENCE INTERVALS
   cimod <- confint(fit.avg)
   # full model-averaged coefficients not needed because all models in the set contain
   # both interaction terms and main effects, or cubic term
   comod <- fit.avg$coef.shrinkage
   CImod <- cbind(comod,cimod) 
   CImod
   # return variables that do not have coefficients that overlap zero
   # (i.e., those that are significant)
   # N.B. doesn't apply to interactions or cubic 
   CIORmod <- exp(CImod)
   
   return(list("CI95"=data.frame(CImod),"odds.ratios"=data.frame(CIORmod)))
   
}