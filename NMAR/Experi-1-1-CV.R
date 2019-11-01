#### 5.1 Situation 1 (Construction of CV)

rm(list=ls())
n <- 1000
ite <- 500
x_cate <- 4
gamma_estimate <- seq(1,1,length.out = ite)
gamma_estimate2 <- seq(1,1,length.out = ite)
var_list2 <- seq(1,1,length.out = ite)
count_list2 <- 0.0

gamma_estimate3 <- seq(1,1,length.out = ite)
gamma_estimate4 <- seq(1,1,length.out = ite)
var_list4 <- seq(1,1,length.out = ite)
count_list4 <- 0.0

gamma_estimate5 <- seq(1,1,length.out = ite)
y_estimate <- seq(1,1,length.out = ite)
y_estimate_2 <- seq(1,1,length.out = ite)
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
  print("confirm_coverage_probability")
  print(sum(delta==1)/n)
  ### Profile ML method
  # likelihood <- function(gamma.){
  #   g_list <- seq(1,1,length.out =  x_cate)
  #   ccc1 <- delta*exp(gamma.*y_transform)
  #   ccc2 <- (1-delta)
  #   for(kkkk in 0:x_cate-1){
  #     g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  #   }
  #   g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  #   pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
  #   ccc3 <- delta*(1.0-pseudo_pi)/pseudo_pi
  #   ccc4 <-  delta/pseudo_pi
  #   g_list3 <- c(1,1,2*x_cate)
  #   for(i in 0:x_cate-1){
  #     for(j in 0:1){
  #       g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #     }
  #   }
  #   g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #   aaa <- delta*log(pseudo_pi)+(1-delta)*log(g_list4)
  #   return(-sum(aaa))
  # }
  # 
  # boku <- optimize(likelihood,c(-2.0,2.0))
  # gamma_estimate[iteite] <- boku$minimum
  # 
  # all <- data.frame(x_transform,z_transform,y_transform,pipi,comb_x_z_transform)
  known_gamma <- 0.8
  
  
  #### Wang&Shao method (Optimal IV)
  
  likelihood2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
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
  boku2 <- optimize(likelihood2,c(-2.0,3.0))
  gamma_estimate2[iteite] <- boku2$minimum
  
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
    aaa <- (1-pseudo_pi)/pseudo_pi*(g_list6-g_list4)**2
    bbb <- (1-pseudo_pi)/pseudo_pi*delta*y_transform*(g_list6-g_list4)
    return(mean(aaa)/(n*mean(bbb)*mean(bbb)))
    ###return(var(aaa))
  }
  print(cf2(true_gamma))
  var_list2[iteite] <- cf2(true_gamma)
  
  ### Calculate CI
  if(-1.96*sqrt(cf2(boku2$minimum))+boku2$minimum<true_gamma){
    if(1.96*sqrt(cf2(boku2$minimum))+boku2$minimum>true_gamma){
      count_list2 <-count_list2 + 1.0 
    }
  }
  
  #### (Score method)
  likelihood3 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
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
    aaa <- delta*(delta-pseudo_pi)*y_transform+(1-delta)*g_list4
    return(sum(aaa)*sum(aaa))
  }
  boku3 <- optimize(likelihood3,c(-2.0,2.0))
  gamma_estimate3[iteite] <- boku3$minimum
  
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
  boku4 <- optimize(likelihood4,c(-2.0,3.0))
  gamma_estimate4[iteite] <- boku4$minimum
  
  ### CF
  
  cf4_2 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y_transform)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
    ###pseudo_pi <- 1.0/(1.0+exp(-beta-beta2*(x_transform-1.4)*(x_transform-1.4)-beta3*(x_transform)-gamma*y_transform))
    g_list3 <- seq(1,1, length.out = x_cate)
    ccc3 <- delta*exp(gamma.*y_transform)*z_transform
    ccc4 <- delta*exp(gamma.*y_transform)
    for(kkkk in 0:x_cate-1){
      g_list3[kkkk+1] <- sum(ccc3[x_transform==kkkk])/sum(ccc4[x_transform==kkkk])
    }
    
    ###aaa <- (delta/pseudo_pi-1.0)*(z_transform-apply(apply(x, 2, function(x) x*g_list3),2,sum))
    aaa <- (1-pseudo_pi)/pseudo_pi*(z_transform-apply(apply(x, 2, function(x) x*g_list3),2,sum))**2
    bbb <- (1-pseudo_pi)/pseudo_pi*delta*y_transform*(z_transform-apply(apply(x, 2, function(x) x*g_list3),2,sum))
    ##aaa <- (1-pseudo_pi)/pseudo_pi*(z_transform)**2
    ###bbb <- (1-pseudo_pi)*y_transform*(z_transform)
    return(mean(aaa)/(n*mean(bbb)*mean(bbb)))
  }
  print(cf4_2(true_gamma))
  var_list4[iteite] <- cf4_2(true_gamma)
  
  ### Calculate CI
  if(-1.96*sqrt(cf4_2(boku4$minimum))+boku4$minimum<true_gamma){
    if(1.96*sqrt(cf4_2(boku4$minimum))+boku4$minimum>true_gamma){
      count_list4 <-count_list4 + 1.0 
    }
  }
  
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
}

print(mean((gamma_estimate-true_gamma)*(gamma_estimate-true_gamma)))
print(mean((gamma_estimate2-true_gamma)*(gamma_estimate2-true_gamma)))
print(mean(var_list2))
print(count_list2/ite)
print(mean((gamma_estimate3-true_gamma)*(gamma_estimate3-true_gamma)))
print(mean((gamma_estimate4-true_gamma)*(gamma_estimate4-true_gamma)))
print(mean((gamma_estimate5-true_gamma)*(gamma_estimate5-true_gamma)))
print(mean(gamma_estimate))
print(mean(gamma_estimate2))
print(mean(gamma_estimate3))
print(mean(gamma_estimate4))
print(mean(var_list4))
print(count_list4/ite)
print(mean(gamma_estimate5))