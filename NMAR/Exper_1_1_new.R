### Setting (b) Compare all estiamtors + Estimate y

rm(list=ls())
n <- 4000
ite <- 300
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
true_y <- 0.6159
count_list2 <- 0.0
count_list4 <- 0.0
var_list2_out <- seq(1,1,length.out = ite)
var_list4_out <- seq(1,1,length.out = ite)

y_estimate <- seq(1,1,length.out = ite)
y_estimate2 <- seq(1,1,length.out = ite)
y_estimate2_ipw <- seq(1,1,length.out = ite)
y_estimate2_out <- seq(1,1,length.out = ite)
y_estimate2_db <- seq(1,1,length.out = ite)
var_list2_out <- seq(1,1,length.out = ite)
y_estimate2_var <- seq(1,1,length.out = ite)
####y_estimate3 <- seq(1,1,length.out = ite)
y_estimate4 <- seq(1,1,length.out = ite)
y_estimate4_ipw <- seq(1,1,length.out = ite)
y_estimate4_out <- seq(1,1,length.out = ite)
y_estimate5 <- seq(1,1,length.out = ite)
y_estimate5_ipw <- seq(1,1,length.out = ite)
y_estimate5_out <- seq(1,1,length.out = ite)
y_estimate7_ipw <- seq(1,1,length.out = ite)
y_estimate9_out <- seq(1,1,length.out = ite)
y_estimate9_ipw <- seq(1,1,length.out = ite)
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
  y_transform[delta==0] <- 0.0
  ## Profile ML method
  
  likelihood1 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    
    
    ccc5 <- delta*exp(gamma.*y_transform)*y_transform
    ccc6 <-  delta*exp(gamma.*y_transform)
    g_list5 <- seq(1,1, length.out = x_cate)
    for(kkkk in 0:x_cate-1){
      g_list5[kkkk+1] <- sum(ccc5[x_transform==kkkk])/sum(ccc6[x_transform==kkkk])
    }
    g_list6 <- apply(apply(comb_x_z, 2, function(x) x*g_list5),2,sum)
    
    ccc3 <- delta*exp(gamma.*y_transform)*(pseudo_pi*(y_transform-g_list6))
    ccc4 <-  delta*exp(gamma.*y_transform)
    
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
    aaa <- delta*(1-pseudo_pi)*(y_transform-g_list6)-(1-delta)*g_list4
    return(sum(aaa)*sum(aaa))
  }
  
  boku9 <- optimize(likelihood1,c(-2.0,2.0))
  gamma_estimate9[iteite] <- boku9$minimum
  #### IPW 
  g_list <- seq(1,1, length.out = x_cate)
  ccc1 <- delta*exp(boku9$minimum*y_transform)
  ccc2 <- (1-delta)
  for(kkkk in 0:x_cate-1){
    g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  }
  g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku9$minimum*y_transform)))
  y_estimate9_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
  ### Outcome
  ccc3 <- delta*exp(boku9$minimum*y_transform)*y_transform
  ccc4 <-  delta*exp(boku9$minimum*y_transform)
  g_list3 <- c(1,1,2*x_cate)
  for(i in 0:x_cate-1){
    for(j in 0:1){
      g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
    }
  }
  g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  y_estimate9_out[iteite] <- mean(g_list4[delta==0])*sum(delta==0)/n+mean(y_transform[delta==1])*sum(delta==1)/n
  
  
  
  #### Wang&Shao method (Optimal IV)
  
  likelihood2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    ccc3 <- delta*exp(gamma.*y_transform)*(-pseudo_pi*y_transform)
    ccc4 <-  delta*exp(gamma.*y_transform)
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
    aaa <- (delta/pseudo_pi-1.0)*g_list4
    return(sum(aaa)*sum(aaa))
  }
  likelihood2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
  
    
    ccc5 <- delta*exp(gamma.*y_transform)*y_transform
    ccc6 <-  delta*exp(gamma.*y_transform)
    g_list5 <- seq(1,1, length.out = x_cate)
    for(kkkk in 0:x_cate-1){
      g_list5[kkkk+1] <- sum(ccc5[x_transform==kkkk])/sum(ccc6[x_transform==kkkk])
    }
    g_list6 <- apply(apply(x, 2, function(x) x*g_list5),2,sum)
    
    ccc3 <- delta*exp(gamma.*y_transform)*(-pseudo_pi*(y_transform-g_list6))
    ccc4 <-  delta*exp(gamma.*y_transform)
    
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
    aaa <- (delta/pseudo_pi-1.0)*g_list4
    return(sum(aaa)*sum(aaa))
  }
  boku2<- optimize(likelihood2,c(-2.0,2.0))
  gamma_estimate2[iteite] <- boku2$minimum
  print(boku2$minimum)
  
  ##### Y estimate
  
  #### IPW 
  g_list <- seq(1,1, length.out = x_cate)
  ccc1 <- delta*exp(boku2$minimum*y_transform)
  ccc2 <- (1-delta)
  for(kkkk in 0:x_cate-1){
    g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  }
  g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku2$minimum*y_transform)))
  y_estimate2_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
  ### Outcome
  ccc3 <- delta*exp(boku2$minimum*y_transform)*y_transform
  ccc4 <-  delta*exp(boku2$minimum*y_transform)
  g_list3 <- c(1,1,2*x_cate)
  for(i in 0:x_cate-1){
    for(j in 0:1){
      g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
    }
  }
  g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  y_estimate2_out[iteite] <- mean(g_list4[delta==0])*sum(delta==0)/n+mean(y_transform[delta==1])*sum(delta==1)/n
  
  ### Estimave_variance
  cf2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    ###pseudo_pi <- 1.0/(1.0+exp(-beta-beta2*(x_transform-1.4)*(x_transform-1.4)-beta3*(x_transform)-gamma.*y_transform))
    
    ccc3 <- delta*exp(gamma.*y_transform)*(y_transform)*pseudo_pi
    ccc4 <-  delta*exp(gamma.*y_transform)
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
    
    ccc5 <- delta*exp(gamma.*y_transform)*(y_transform)*pseudo_pi
    ccc6 <- delta*exp(gamma.*y_transform)
    g_list5 <- c(1,1,x_cate)
    for(kkkk in 0:x_cate-1){
      g_list5[kkkk+1] <- sum(ccc5[x_transform==kkkk])/sum(ccc6[x_transform==kkkk])
    }
    g_list6 <- apply(apply(x, 2, function(x) x*g_list5),2,sum)
    
    ##aaa <- (delta/pseudo_pi-1.0)*(g_list6-apply(apply(x, 2, function(x) x*g_list3),2,sum))
    aaa <- (delta/pseudo_pi-1)*(g_list6-g_list4)
    bbb <- (1-pseudo_pi)*(y_transform*delta/pseudo_pi-g_list6)*(g_list6-g_list4)
    return(aaa/(mean(bbb)))
    ###return(var(aaa))
  }
  
  cf2_2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    ###pseudo_pi <- 1.0/(1.0+exp(-beta-beta2*(x_transform-1.4)*(x_transform-1.4)-beta3*(x_transform)-gamma.*y_transform))
    
    ccc3 <- delta*exp(gamma.*y_transform)*(y_transform)
    ccc4 <-  delta*exp(gamma.*y_transform)
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  
    ##aaa <- (delta/pseudo_pi-1.0)*(g_list6-apply(apply(x, 2, function(x) x*g_list3),2,sum))
    ###aaa <- (dela/pseudo_pi-1)*(g_list6-g_list4)
    ccc1 <- (1-pseudo_pi)/pseudo_pi*delta*(y_transform- g_list4)**2
    return(ccc1)
    ###return(var(aaa))
  }
  print("coverage")
  ###print(cf2_3(boku2$minimum)*1000/n)
  auxi <- cf2(true_gamma)
  auxi2 <- cf2_2(true_gamma)
  ###print(var(auxi))
   ### print(auxi2)
  ###print(auxi3)
  print(var(g_list4[delta==0]+(delta/pseudo_pi)*(y_transform-g_list4[delta==0])+auxi2*auxi)*1000/n)
  print(var(g_list4[delta==0]+(delta/pseudo_pi)*(y_transform-g_list4[delta==0]))*1000/n)
  var_list2_out[iteite] <-  var(g_list4[delta==0]+(delta/pseudo_pi)*(y_transform-g_list4[delta==0])+auxi2*auxi)/n
  
  ### Calculate CI
  if(-1.96*sqrt(var_list2_out[iteite])+y_estimate2_out[iteite]<true_y){
    if(1.96*sqrt(var_list2_out[iteite])+y_estimate2_out[iteite]>true_y){
      count_list2 <-count_list2 + 1.0 
    }
  }
  ### DB
  g_list <- seq(1,1, length.out = x_cate)
  ccc1 <- delta*exp(boku2$minimum*y_transform)
  ccc2 <- (1-delta)
  for(kkkk in 0:x_cate-1){
    g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  }
  g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku2$minimum*y_transform)))
  
  ccc3 <- delta*exp(boku2$minimum*y_transform)*y_transform
  ccc4 <-  delta*exp(boku2$minimum*y_transform)
  g_list3 <- c(1,1,2*x_cate)
  for(i in 0:x_cate-1){
    for(j in 0:1){
      g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
    }
  }
  g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  y_estimate2_db[iteite] <- mean(delta*y_transform/pseudo_pi)+mean((1-delta/pseudo_pi)*g_list4)
  
  
  ### Not optimal IV
  likelihood4 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    aaa <- (delta/pseudo_pi-1.0)*z_transform
    return(sum(aaa)*sum(aaa))
  }
  boku4 <- optimize(likelihood4,c(-2.0,2.0))
  ###print(boku4$minimum)
  gamma_estimate4[iteite] <- boku4$minimum
  
  #### IPW 
  g_list <- seq(1,1, length.out = x_cate)
  ccc1 <- delta*exp(boku4$minimum*y_transform)
  ccc2 <- (1-delta)
  for(kkkk in 0:x_cate-1){
    g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  }
  g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku4$minimum*y_transform)))
  y_estimate4_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
  ### Outcome
  ccc3 <- delta*exp(boku4$minimum*y_transform)*y_transform
  ccc4 <-  delta*exp(boku4$minimum*y_transform)
  g_list3 <- c(1,1,2*x_cate)
  for(i in 0:x_cate-1){
    for(j in 0:1){
      g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
    }
  }
  g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  y_estimate4_out[iteite] <- mean(g_list4[delta==0])*sum(delta==0)/n+mean(y_transform[delta==1])*sum(delta==1)/n
  
  #### Wang&Shao method (Optimal IV)
  
  likelihood5 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    ccc3 <- delta*(1-pseudo_pi)/(pseudo_pi)*(y_transform)
    ccc4 <-  delta*(1-pseudo_pi)/(pseudo_pi*pseudo_pi)
    g_list3 <- c(1,1,2*x_cate)
    for(i in 0:x_cate-1){
      for(j in 0:1){
        g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
      }
    }
    g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
    aaa <- (delta/pseudo_pi-1.0)*g_list4
    return(sum(aaa)*sum(aaa))
  }
  boku5 <- optimize(likelihood5,c(-2.0,2.0))
  gamma_estimate5[iteite] <- boku5$minimum
  
  #### IPW 
  g_list <- seq(1,1, length.out = x_cate)
  ccc1 <- delta*exp(boku5$minimum*y_transform)
  ccc2 <- (1-delta)
  for(kkkk in 0:x_cate-1){
    g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  }
  g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku5$minimum*y_transform)))
  y_estimate5_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
  ### Outcome
  ccc3 <- delta*exp(boku5$minimum*y_transform)*y_transform
  ccc4 <-  delta*exp(boku5$minimum*y_transform)
  g_list3 <- c(1,1,2*x_cate)
  for(i in 0:x_cate-1){
    for(j in 0:1){
      g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
    }
  }
  g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  y_estimate5_out[iteite] <- mean(g_list4[delta==0])*sum(delta==0)/n+mean(y_transform[delta==1])*sum(delta==1)/n
  
  ### Estimave_variance
  cf4 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    
    ccc5 <- delta*exp(gamma.*y_transform)*z_transform
    ccc6 <- delta*exp(gamma.*y_transform)
    g_list5 <- c(1,1,x_cate)
    for(kkkk in 0:x_cate-1){
      g_list5[kkkk+1] <- sum(ccc5[x_transform==kkkk])/sum(ccc6[x_transform==kkkk])
    }
    g_list6 <- apply(apply(x, 2, function(x) x*g_list5),2,sum)
    
    ccc7 <- delta*exp(gamma.*y_transform)*y_transform
    ccc8 <- delta*exp(gamma.*y_transform)
    g_list7 <- c(1,1,x_cate)
    for(kkkk in 0:x_cate-1){
      g_list7[kkkk+1] <- sum(ccc7[x_transform==kkkk])/sum(ccc8[x_transform==kkkk])
    }
    g_list8 <- apply(apply(x, 2, function(x) x*g_list7),2,sum)
    
    aaa <- (delta/pseudo_pi-1)*(z_transform-g_list6)
    bbb <- (1-pseudo_pi)*(y_transform*delta/pseudo_pi-g_list8)*(z_transform-g_list6)
    return(aaa/(mean(bbb)))
  }

  print("coverage")
  ###print(cf2_3(boku2$minimum)*1000/n)
  auxi <- cf4(true_gamma)
  auxi2 <- cf2_2(true_gamma)
  
  ###print(var(auxi))
  print(var(g_list4[delta==0]+(delta/pseudo_pi)*(y_transform-g_list4[delta==0])+auxi2*auxi)*1000/n)
  var_list4_out[iteite] <-  var(g_list4[delta==0]+(delta/pseudo_pi)*(y_transform-g_list4[delta==0])-auxi2*auxi)/n

  ### Calculate CI
  if(-1.96*sqrt(var_list4_out[iteite])+y_estimate4_out[iteite]<true_y){
    if(1.96*sqrt(var_list4_out[iteite])+y_estimate4_out[iteite]>true_y){
      count_list4 <-count_list4 + 1.0
    }
  }
    
  likelihood7 <- function(gamma.){
    g_list2 <- exp(gamma.[1]+gamma.[2]*x_transform)
    ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.[3]*y_transform)))
    aaa <- (delta/pseudo_pi-1.0)*1.0
    aaa1 <- (delta/pseudo_pi-1.0)*x_transform
    aaa2 <- (delta/pseudo_pi-1.0)*z_transform
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  ### boku7 <- optim(c(0.4,0.4,0.3),likelihood7)
  ###gamma_estimate7[iteite] <- boku7$par[3]
  #### IPW 
  ##g_list2 <- exp(boku7$par[1]+boku7$par[2]*x_transform)
  ##pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku7$par[3]*y_transform)))
  ###y_estimate7_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
  # likelihood8 <- function(gamma.){
  #   g_list2 <- exp(0.2+0.4*x_transform)
  #   pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
  #   ccc3 <- delta*exp(gamma.*y_transform)*(-pseudo_pi*y_transform)
  #   ccc4 <-  delta*exp(gamma.*y_transform)
  #   g_list3 <- c(1,1,2*x_cate)
  #   for(i in 0:x_cate-1){
  #     for(j in 0:1){
  #       g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #     }
  #   }
  #   g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #   aaa <- delta*(delta-pseudo_pi)*y_transform+(1-delta)*g_list4
  #   return(sum(aaa)*sum(aaa))
  # }
  # boku8 <- optimise(likelihood8,c(-2,2))
  # gamma_estimate8[iteite] <- boku8$minimum
}

# print(mean((gamma_estimate-true_gamma)*(gamma_estimate-true_gamma)))
print(mean((gamma_estimate2-true_gamma)*(gamma_estimate2-true_gamma)))
# print(mean((gamma_estimate3-true_gamma)*(gamma_estimate3-true_gamma)))
print(mean((gamma_estimate4-true_gamma)*(gamma_estimate4-true_gamma)))
print(mean((gamma_estimate5-true_gamma)*(gamma_estimate5-true_gamma)))
# print(mean((gamma_estimate6-true_gamma)*(gamma_estimate6-true_gamma)))
# print(mean((gamma_estimate7-true_gamma)*(gamma_estimate7-true_gamma)))
# print(mean((gamma_estimate8-true_gamma)*(gamma_estimate8-true_gamma)))
 print(mean((gamma_estimate9-true_gamma)*(gamma_estimate9-true_gamma)))
# print(mean(gamma_estimate))
# print(mean(gamma_estimate2))
# print(mean(gamma_estimate3))
# print(mean(gamma_estimate4))
# print(mean(gamma_estimate5))
# print(mean(gamma_estimate6))
# print(mean(gamma_estimate7))
# print(mean(gamma_estimate8))
# print(mean(gamma_estimate9))

print(mean(y_estimate2_ipw))
print(mean(y_estimate2_out))
print(mean(y_estimate2_db))
print(mean(y_estimate4_ipw))
print(mean(y_estimate4_out))
print(mean(y_estimate5_ipw))
print(mean(y_estimate5_out))
true_y <- 0.6159
print(mean((y_estimate2_out-true_y)*(y_estimate2_out-true_y))*1000)
print(mean((y_estimate2_ipw-true_y)*(y_estimate2_ipw-true_y))*1000)
print(mean((y_estimate2_db-true_y)*(y_estimate2_db-true_y))*1000)
print(mean((y_estimate4_out-true_y)*(y_estimate4_out-true_y))*1000)
print(mean((y_estimate4_ipw-true_y)*(y_estimate4_ipw-true_y))*1000)
print(mean((y_estimate5_out-true_y)*(y_estimate5_out-true_y))*1000)
print(mean((y_estimate9_ipw-true_y)*(y_estimate9_ipw-true_y))*1000)
print(mean((y_estimate9_out-true_y)*(y_estimate9_out-true_y))*1000)
print("aaa")
print(count_list2/ite)
print(count_list4/ite)
boxplot(gamma_estimate4,gamma_estimate2,gamma_estimate9,gamma_estimate5)

abline(h=0.6)
boxplot(y_estimate4_ipw,y_estimate4_out,y_estimate9_ipw, y_estimate9_out,y_estimate2_ipw, y_estimate2_out,y_estimate5_ipw, y_estimate5_out)
abline(h=0.6159)
###print(mean((y_estimate7_ipw-true_y)*(y_estimate7_ipw-true_y))*1000)

png("1_1_MU.png")
par(mar=c(7,3,2,1))
ccc <- data.frame(x= c(c(y_estimate4_ipw, y_estimate4_out,y_estimate4_out),c(y_estimate9_ipw,y_estimate9_out, y_estimate9_out),c(y_estimate2_ipw,y_estimate2_out, y_estimate2_out),c(y_estimate5_ipw,y_estimate5_out,y_estimate5_out)),
                  y=rep(c("GMM","SCORE","CA1","CA2_A"),each=ite*3),
                  z=rep(rep(c('IPW',"MP","DR"), each=ite),4),
                  stringsAsFactors = FALSE
)

boxplot(ccc$x ~ factor(ccc$z,levels=c('IPW',"MP","DR")) + factor(ccc$y,levels=c("GMM","SCORE","CA1","CA2_A")) ,las = 2,at=c(1:3,5:7,9:11,13:15),main="M1: MU",las = 2)
###colnames(ccc) <- c("IPW_GMM","MP_GMM","DR_GMM","IPW_SCORE","MP_SCORE","DR_SCORE","IPW_CA1","MP_CA1","DR_CA1","IPW_CA2_A","MP_CA2_A","DR_CA2_A","IPW_CA2_S","MP_CA2_S","DR_CA2_S")
###ccc <- data.frame(IPW_SCORE=y_estimate4_ipw,DR_SCORE=y_estimate4_out,IPW_GMM=y_estimate9_ipw, DR_GMM=y_estimate9_out,y_estimate4_ipw, y_estimate2_out,y_estimate5_ipw, y_estimate5_out)
abline(h=0.6159)
###boxplot(ccc,las = 2)
dev.off()


png("1_1_GAMMA.png")
ccc <- data.frame(gamma_estimate4 , gamma_estimate9 ,gamma_estimate2 ,gamma_estimate5)
colnames(ccc) <- c("GMM","SCORE","CA1","CA2")
###ccc <- data.frame(IPW_SCORE=y_estimate4_ipw,DR_SCORE=y_estimate4_out,IPW_GMM=y_estimate9_ipw, DR_GMM=y_estimate9_out,y_estimate4_ipw, y_estimate2_out,y_estimate5_ipw, y_estimate5_out)
###boxplot(ccc,las = -0.8862813)
boxplot(ccc,las = 2,main="M1: GAMMA")
abline(h=0.6)
dev.off()
