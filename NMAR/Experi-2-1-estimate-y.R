rm(list=ls())
n <- 50000
ite <- 100

y_list <- seq(0,0,length.out = ite)

x_cate <- 2
for(iteite in 1:ite){
  print(iteite)
  x <- rmultinom(n,1,seq(1.0/x_cate,1.0/x_cate,length.out = x_cate))
  x1 <- rbinom(n,1,0.5)
  for(i in 1:n){
    x1[i] <- x[,i] %*% seq(0,x_cate-1.0,length.out = x_cate)
  }
  x2 <- runif(n,-1,1)
  y <- seq(0,0,length.out = n)
  delta <- seq(0,0,length.out = n)
  marge_pipi<- seq(0,0,length.out = n)
  mu <- -1.0-0.4*x1+0.5*x2*x2
  sigma_ <- 1.0
  alpha <- -0.3
  beta <-  -0.4
  gamma <- 0.5
  true_gamma <- gamma
  
  
  ### Marge PS
  phi <- sigma_*sigma_
  aaa <- (0.5*(gamma*phi+mu)*(gamma*phi+mu)-0.5*mu*mu)/phi
  for(i in 1:n){
    marge_pipi[i] <- 1.0/(1.0+exp(alpha+beta*x1[i]+aaa[i]))
    delta[i] <- rbinom(1,1,marge_pipi[i])
  }
  #### PS model
  sigmoid <- function(xxx){
    return(1.0/(1.0+exp(alpha+beta*x[i]+gamma*y[i])))
  }
  
  ### 
  for(i in 1:n){
    if(delta[i]==1){
      y[i] <- mu[i]+ rnorm(1,0,sigma_)
    }
    else{
      y[i] <- mu[i]+phi*gamma+rnorm(1,0,sigma_)
    }
  }
  y_list[iteite] <- mean(y)
}
print(mean(y_list))