# SK 14/02/2014

## ARE QUADRATIC TERMS REQUIRED FOR ANY PREDICTOR VARIABLES?

# OUTPUTS
# pasted message of whether quadratic is needed

lin.quad <- function(col.numbers,dataset,error.distribution){

   for(i in min(col.numbers):max(col.numbers)){
      diff.BIC <- BIC(glm(Spawn~poly(dataset[,i],1),data=dataset,
                              family=error.distribution))-
                     BIC(glm(Spawn~poly(dataset[,i],2),data=dataset,
                              family=error.distribution))
      print(paste(colnames(dataset)[i]," difference in BIC =",round(diff.BIC,2)))
      if (diff.BIC > 3) print("Use quadratic") else print("Use linear")
   }

}