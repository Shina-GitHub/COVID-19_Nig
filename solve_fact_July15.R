#library(rootSolve)
model<- function(x,p){
  F1 = rbinom(1,x,p)
  return (F1)
}


#ss <- multiroot(f = model, start = 100,parms=c(p=0.5, report= 30, rand_var= runif(1)))
local_adjusted<-function(report, prob){
  upper_value = 5000
  solu_ideas = report:upper_value  #solution guess
  solu = rep(0,1) # 1000 realization
  solu_ideas = matrix(solu_ideas,ncol =1)
  for (i in seq(1,1)){
      myguess = apply(solu_ideas,1,model, p= prob)
      
      solu[i] = median(solu_ideas[myguess < (report+2) & myguess>= (report-2)])

  }

  myresult = round(median(solu))

  return(myresult)
}




