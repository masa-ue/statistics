### Setting (b) Compare all estiamtors + Estimate y

rm(list=ls())
n <- 1000
ite <-200
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

y_estimate <- seq(1,1,length.out = ite)
y_estimate2 <- seq(1,1,length.out = ite)
y_estimate2_ipw <- seq(1,1,length.out = ite)
y_estimate2_out <- seq(1,1,length.out = ite)
y_estimate2_db <- seq(1,1,length.out = ite)
####y_estimate3 <- seq(1,1,length.out = ite)
y_estimate4 <- seq(1,1,length.out = ite)
y_estimate4_ipw <- seq(1,1,length.out = ite)
y_estimate4_out <- seq(1,1,length.out = ite)
y_estimate5 <- seq(1,1,length.out = ite)
y_estimate5_ipw <- seq(1,1,length.out = ite)
y_estimate5_out <- seq(1,1,length.out = ite)
y_estimate7_ipw <- seq(1,1,length.out = ite)
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
  beta <- 1.6
  beta2 <- -0.8
  beta3 <- 0.7
  gamma <- -0.6
  true_gamma <- -gamma
  
  #### PS model
  sigmoid <- function(xxx,yyy){
    return(1.0/(1.0+exp(-beta-beta2*sin(xxx)-gamma*yyy)))
  }
  
  for(i in 1:n){
    pipi[i] <- sigmoid(x_transform[i],y_transform[i])
    delta[i] <- rbinom(1,1,pipi[i])
  }
  ### Profile ML method
  # likelihood1 <- function(gamma.){
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
  # boku1 <- optimize(likelihood1,c(-2.0,2.0))
  # gamma_estimate[iteite] <- boku1$minimum
  # 
  # gamma_estimate2_list =seq(0,0,length.out = 100)
  # gamma_known <- 0.8
  # for(tttt in 1:5){
  #   likelihood2 <- function(gamma.,gamma_known. = gamma_known){
  #     g_list <- seq(1,1, length.out = x_cate)
  #     ccc1 <- delta*exp(gamma.*y_transform)
  #     ccc2 <- (1-delta)
  #     for(kkkk in 0:x_cate-1){
  #       g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  #     }
  #     g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
  #     ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
  #     pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y_transform)))
  #     ccc3 <- delta*exp(gamma_known.*y_transform)*(-pseudo_pi*y_transform)
  #     ccc4 <-  delta*exp(gamma_known.*y_transform)
  #     g_list3 <- c(1,1,2*x_cate)
  #     for(i in 0:x_cate-1){
  #       for(j in 0:1){
  #         g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #       }
  #     }
  #     g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #     aaa <- (delta/pseudo_pi-1.0)*g_list4
  #     return(sum(aaa)*sum(aaa))
  #   }
  #   
  #   boku2 <- optimize(likelihood2,c(-2.0,2.0))
  #   gamma_known <- boku2$minimum
  #   print(gamma_known)
  # }
  # gamma_estimate9[iteite] <- gamma_known
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
  
  # #### (Score method)
  # likelihood3 <- function(gamma.){
  #   g_list <- seq(1,1, length.out = x_cate)
  #   ccc1 <- delta*exp(gamma.*y_transform)
  #   ccc2 <- (1-delta)
  #   for(kkkk in 0:x_cate-1){
  #     g_list[kkkk+1] <- sum(ccc1[x_transform==kkkk])/sum(ccc2[x_transform==kkkk])
  #   }
  #   g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
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
  # boku3 <- optimize(likelihood3,c(-2.0,2.0))
  # print(boku3$minimum)
  # gamma_estimate3[iteite] <- boku3$minimum
  # 
  
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
  print(boku4$minimum)
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
  
  # ### Paramaetirc model
  # likelihood6 <- function(gamma.){
  #   g_list2 <- exp(gamma.[1]+gamma.[2]*x_transform)
  #   ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
  #   pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.[3]*y_transform)))
  #   ccc3 <- delta*exp(gamma.[3]*y_transform)*(-pseudo_pi*y_transform)
  #   ccc4 <-  delta*exp(gamma.[3]*y_transform)
  #   g_list3 <- c(1,1,2*x_cate)
  #   for(i in 0:x_cate-1){
  #     for(j in 0:1){
  #       g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #     }
  #   }
  #   g_list4 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #   
  #   ccc3 <- delta*exp(gamma.[3]*y_transform)*(-pseudo_pi*x_transform)
  #   ccc4 <-  delta*exp(gamma.[3]*y_transform)
  #   g_list3 <- c(1,1,2*x_cate)
  #   for(i in 0:x_cate-1){
  #     for(j in 0:1){
  #       g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #     }
  #   }
  #   g_list5 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #   
  #   ccc3 <- delta*exp(gamma.[3]*y_transform)*(-pseudo_pi*1)
  #   ccc4 <-  delta*exp(gamma.[3]*y_transform)
  #   g_list3 <- c(1,1,2*x_cate)
  #   for(i in 0:x_cate-1){
  #     for(j in 0:1){
  #       g_list3[1+i+x_cate*j]<- sum(ccc3[comb_x_z_transform==(i+x_cate*j)])/sum(ccc4[comb_x_z_transform==(i+x_cate*j)])
  #     }
  #   }
  #   g_list6 <- apply(apply(comb_x_z, 2, function(x) x*g_list3),2,sum)
  #   aaa <- (delta/pseudo_pi-1.0)*g_list4
  #   aaa1 <- (delta/pseudo_pi-1.0)*g_list5
  #   aaa2 <- (delta/pseudo_pi-1.0)*g_list6
  #   return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  # }
  # 
  # boku6 <- optim(c(0.4,0.4,0.3),likelihood6)
  # gamma_estimate6[iteite] <- boku6$par[3]
  
  likelihood7 <- function(gamma.){
    g_list2 <- exp(gamma.[1]+gamma.[2]*x_transform)
    ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.[3]*y_transform)))
    aaa <- (delta/pseudo_pi-1.0)*1.0
    aaa1 <- (delta/pseudo_pi-1.0)*x_transform
    aaa2 <- (delta/pseudo_pi-1.0)*z_transform
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  boku7 <- optim(c(0.4,0.4,0.3),likelihood7)
  gamma_estimate7[iteite] <- boku7$par[3]
  #### IPW 
  g_list2 <- exp(boku7$par[1]+boku7$par[2]*x_transform)
  pseudo_pi <- as.vector(g_list2/(g_list2+exp(boku7$par[3]*y_transform)))
  y_estimate7_ipw[iteite] <- mean(y_transform/pseudo_pi)
  
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
# print(mean((gamma_estimate2-true_gamma)*(gamma_estimate2-true_gamma)))
# print(mean((gamma_estimate3-true_gamma)*(gamma_estimate3-true_gamma)))
# print(mean((gamma_estimate4-true_gamma)*(gamma_estimate4-true_gamma)))
# print(mean((gamma_estimate5-true_gamma)*(gamma_estimate5-true_gamma)))
# print(mean((gamma_estimate6-true_gamma)*(gamma_estimate6-true_gamma)))
# print(mean((gamma_estimate7-true_gamma)*(gamma_estimate7-true_gamma)))
# print(mean((gamma_estimate8-true_gamma)*(gamma_estimate8-true_gamma)))
# print(mean((gamma_estimate9-true_gamma)*(gamma_estimate9-true_gamma)))
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
true_y <- 0.6158693
print(mean((y_estimate2_out-true_y)*(y_estimate2_out-true_y))*1000)
print(mean((y_estimate2_ipw-true_y)*(y_estimate2_ipw-true_y))*1000)
print(mean((y_estimate2_db-true_y)*(y_estimate2_db-true_y))*1000)
print(mean((y_estimate4_out-true_y)*(y_estimate4_out-true_y))*1000)
print(mean((y_estimate4_ipw-true_y)*(y_estimate4_ipw-true_y))*1000)
print(mean((y_estimate5_out-true_y)*(y_estimate5_out-true_y))*1000)
print(mean((y_estimate5_ipw-true_y)*(y_estimate5_ipw-true_y))*1000)
####print(mean((y_estimate7_ipw-true_y)*(y_estimate7_ipw-true_y))*1000)



