n <-  100000
ite <- 200
x_cate <- 4
gamma_estimate <- seq(1,1,length.out = ite)
gamma_estimate2 <- seq(1,1,length.out = ite)
gamma_estimate3 <- seq(1,1,length.out = ite)
gamma_estimate4 <- seq(1,1,length.out = ite)
gamma_estimate5 <- seq(1,1,length.out = ite)
gamma_estimate6 <- seq(1,1,length.out = ite)
gamma_estimate7 <- seq(1,1,length.out = ite)
gamma_estimate8 <- seq(1,1,length.out = ite)
gamma_estimate9<- seq(1,1,length.out = ite)
y_true_mean <- seq(1,1,length.out = ite)

for(iteite in 1:ite){
  print(iteite)
  x <- rmultinom(n,1,seq(1.0/x_cate,1.0/x_cate,length.out = x_cate))
  y <- matrix(1, nrow = 2, ncol = n)
  z <- matrix(1, nrow = 2, ncol = n)
  x_transform <- seq(1,1,length.out = n)
  z_transform <- seq(1,1,length.out = n)
  y_transform <- seq(1,1,length.out = n)
  comb_x_z <- matrix(1, nrow = x_cate*2, ncol = n)
  comb_x_z_transform  <- seq(1,1,length.out = n)
  ccc1 <- seq(1,1,length.out = n)
  ccc2 <- seq(1,1,length.out = n)
  
  pipi <- seq(1,1,length.out = n)
  delta <- seq(1,1,length.out = n)
  
  for(i in 1:n){
    x_transform[i] <- x[,i] %*% seq(0,x_cate-1.0,length.out = x_cate)
    z_transform[i] <- rbinom(1,1,0.5)
  }
  for(i in 1:n){
    aaa <- 1/(1+exp(1.3-(x_transform[i]-1.6)*1.0-z_transform[i]*1.5))
    ###aaa <- 0.9-0.12*x_transform[i]-z_transform[i]*0.4
    y[,i] <- rmultinom(1,1,c(aaa,1.0-aaa))
    y_transform[i] <- y[,i]%*% c(0,1)
  }
  
  for(i in 1:n){
    comb_x_z[,i] <- seq(0,0,length.out = x_cate)
    kkkk <- x_transform[i]+(z_transform[i])*x_cate+1
    comb_x_z[kkkk,i] <- 1 
  }
  for(i in 1:n){
    ttt <- (2*x_cate)
    comb_x_z_transform[i] <- comb_x_z[,i] %*% seq(0,2*x_cate-1,length.out = ttt)
  }
  
  alpha <- 0.3
  beta <- 0.2
  beta1 <- 0.4
  beta2 <- 0.0
  beta3 <- 0.8
  gamma <- -0.6
  true_gamma <- -gamma
  
  #### PS model
  sigmoid <- function(xxx,yyy){
    return(1.0/(1.0+exp(-beta-beta2*(xxx-1.3)*(xxx-1.3)-beta3*(xxx)-gamma*yyy)))
  }
  
  for(i in 1:n){
    pipi[i] <- sigmoid(x_transform[i],y_transform[i])
    delta[i] <- rbinom(1,1,pipi[i])
  }
  gamma <- 0.0
  y_true_mean[iteite] <- mean(y_transform)
}
mean(y_true_mean)