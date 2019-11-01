rm(list=ls())
n <- 4000
ite <- 20

gamma_estimate <- seq(1,1,length.out = ite)
gamma_estimate2 <- seq(1,1,length.out = ite)
gamma_estimate3 <- seq(1,1,length.out = ite)
gamma_estimate4 <- seq(1,1,length.out = ite)
gamma_estimate5 <- seq(1,1,length.out = ite)
gamma_estimate6 <- seq(1,1,length.out = ite)
gamma_estimate7 <- seq(1,1,length.out = ite)
gamma_estimate8 <- seq(1,1,length.out = ite)

x_cate <- 5
for(iteite in 1:ite){
  print(iteite)
  x <- rmultinom(n,1,seq(1.0/x_cate,1.0/x_cate,length.out = x_cate))
  x1 <- rbinom(n,1,0.5)
  for(i in 1:n){
    x1[i] <- x[,i] %*% seq(0,x_cate-1.0,length.out = x_cate)
  }
  x2 <- runif(n,-1,1)
  ####x2 <- rnorm(n,0,1.0)
  y <- seq(0,0,length.out = n)
  delta <- seq(0,0,length.out = n)
  marge_pipi<- seq(0,0,length.out = n)
  mu <- exp(-1.7-0.4*x1+0.5*x2)
  sigma_ <- 1.0
  alpha <- -1.7
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
  
  likelihood_outcome <- function(gamma.){
    aaa <- 1.0/sqrt(gamma.[1]*gamma.[1])*exp(-0.5*(y[delta==1]-exp(gamma.[2]+gamma.[3]*x1[delta==1]+gamma.[4]*x2[delta==1]))**2/(gamma.[1]*gamma.[1]))
    return(-sum(log(aaa)))
  }
  boku2 <- optim(c(1.0,-1.7,0.4,0.5),likelihood_outcome)
  mu <- exp(boku2$par[2]+boku2$par[3]*x1+boku2$par[4]*x2)
  sigma_ <- boku2$par[1]
  
  likelihood <- function(gamma.,y.=y){
    pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y.))
    aaa <- delta*log(pseudo_pi)+(1-delta)*log(1-pseudo_pi)
    return(-sum(aaa))
  }
  ###boku2 <- optim(c(-1.7,-0.4,0.5),likelihood)
  ###print(boku2$par)
  ##gamma_estimate[iteite] <- boku2$par[3]
  
  likelihood2 <- function(gamma.,y.=y){
    ###gamma.[1] <- alpha
    ####gamma.[2] <- beta
    pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y.))
    ####pseudo_pi <- 1.0/(1.0+exp(-1.7-0.4*x1+gamma.*y.))
    aaa <- (delta/pseudo_pi-1.0)*x2
    aaa1 <- (delta/pseudo_pi-1.0)*x1
    aaa2 <- (delta/pseudo_pi-1.0)
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  ####boku2 <- optim(c(-1.7,-0.4,0.5),likelihood2)
  ####print(boku2$par)
  ####gamma_estimate2[iteite] <- boku2$par[3]
  
  likelihood3 <- function(gamma.,y.=y){
    ###gamma.[1] <- alpha
    ####gamma.[2] <- beta
    pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y.))
    pseudo_pi2 <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*(1.5*gamma.[3]*sigma_*sigma_+mu)))
    aaa <- (delta/pseudo_pi-1.0)*pseudo_pi2*(gamma.[3]*sigma_*sigma_+mu)
    aaa1 <- (delta/pseudo_pi-1.0)*pseudo_pi2
    aaa2 <- (delta/pseudo_pi-1.0)*pseudo_pi2*(x1)
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  ###boku3 <- optim(c(-1.7,-0.4,0.5),likelihood3)
  ###gamma_estimate3[iteite] <- boku3$par[3]
  
  likelihood4 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x1==kkkk])/sum(ccc2[x1==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(1.0/(1.0+exp(-log(g_list2)+gamma.*y)))
    pseudo_pi2 <- as.vector(1.0/(1.0+exp(-log(g_list2)+gamma.*(1.5*gamma.*sigma_*sigma_+mu))))
    # pseudo_pi <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*y))
    # pseudo_pi2 <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*(1.5*gamma.*sigma_*sigma_+mu)))
    ###pseudo_pi <- 1.0/(1.0+exp(-1.7-0.4*x1+gamma.*y))
    ###pseudo_pi2 <- 1.0/(1.0+exp(-1.7-0.4*x1+gamma.*(1.5*gamma.*sigma_*sigma_+mu)))
    aaa <- (delta/pseudo_pi-1.0)*pseudo_pi2*(gamma.*sigma_*sigma_+mu)
    aaa1 <- (delta/pseudo_pi-1.0)*pseudo_pi2
    aaa2 <- (delta/pseudo_pi-1.0)*pseudo_pi2*(x1)
    return(sum(aaa)*sum(aaa)+sum(aaa1)*sum(aaa1)+sum(aaa2)*sum(aaa2))
  }
  
  boku4 <- optimise(likelihood4,c(-1.0,2.0))
  gamma_estimate4[iteite] <- boku4$minimum
  
  likelihood5 <- function(gamma.){
    g_list <- seq(1,1, length.out = x_cate)
    ccc1 <- delta*exp(gamma.*y)
    ccc2 <- (1-delta)
    for(kkkk in 0:x_cate-1){
      g_list[kkkk+1] <- sum(ccc1[x1==kkkk])/sum(ccc2[x1==kkkk])
    }
    g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
    pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y)))
    ######pseudo_pi2 <- as.vector(g_list2/(g_list2+exp(gamma.*(1.5*gamma.*sigma_*sigma_+mu))))
    ###pseudo_pi <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*y))
    ###pseudo_pi2 <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*(1.5*gamma.*sigma_*sigma_+mu)))
    aaa <- (delta/pseudo_pi-1.0)*x2
    return(sum(aaa)*sum(aaa))
  }
  
  boku5 <- optimise(likelihood5,c(-2.0,2.0))
  gamma_estimate5[iteite] <- boku5$minimum
  
  gamma. <- c(-1.7,-0.4,0.5)
  pseudo_pi <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*y))
  pseudo_pi2 <- 1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*(1.5*gamma.[3]*sigma_+mu)))
  pseudo_pi3 <- 1.0-1.0/(1.0+exp(gamma.[1]+gamma.[2]*x1+gamma.[3]*(1.5*gamma.[3]*sigma_+mu)))
  y_exp <- (gamma.[3]*sigma_*sigma_+mu)*pseudo_pi2+(2*gamma.[3]*sigma_*sigma_+mu)*pseudo_pi3
  ####print(c(mean(y[delta==0]),mean(y_exp[delta==0])))
  
  print("EM-P-EE-2")
  gamma_estimate2_list =seq(0,0,length.out = 100)
  gamma_known <- 0.5
  mmmm <- 300
  large_matrix <- matrix(1,n,mmmm)
  for(iii in 1:n){
    large_matrix[iii,] <- rnorm(mmmm, mu[iii], sigma_)
  }
  for(tttt in 1:1){
    likelihood2 <- function(gamma.,gamma_known. = gamma_known){
      g_list <- seq(1,1, length.out = x_cate)
      ccc1 <- delta*exp(gamma.*y)
      ccc2 <- (1-delta)
      for(kkkk in 0:x_cate-1){
        g_list[kkkk+1] <- sum(ccc1[x1==kkkk])/sum(ccc2[x1==kkkk])
      }
      g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
      ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
      pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y)))
      ###pseudo_pi <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*y))
      g_list4 <- seq(1,1, length.out = x_cate)
      for(iii in 1:n){
        sum3 <- seq(1,1,length.out = mmmm)
        sum4 <- seq(1,1,length.out = mmmm)
        for(jjj in 1:mmmm){
          sum3[jjj] <- exp(gamma.*large_matrix[iii,jjj])*large_matrix[iii,jjj]*g_list2[iii]/(g_list2[iii]+exp(gamma.*large_matrix[iii,jjj]))
          sum4[jjj] <- exp(gamma.*large_matrix[iii,jjj])
        }
        g_list4[iii] <- mean(sum3)/mean(sum4)
      }
      aaa <- (delta/pseudo_pi-1)*g_list4
      return(sum(aaa)*sum(aaa))
    }
    boku2 <- optimize(likelihood2,c(-2.0,2.0))
    gamma_known2 <- boku2$minimum
    print(gamma_known2)
  }
  gamma_estimate6[iteite] <- gamma_known2
  
  print("EM-P-Score")
  
  gamma_estimate2_list =seq(0,0,length.out = 100)
  gamma_known <- 0.5
  mmmm <- 300
  large_matrix <- matrix(1,n,mmmm)
  for(iii in 1:n){
    large_matrix[iii,] <- rnorm(mmmm, mu[iii], sigma_)
  }
  for(tttt in 1:1){
    likelihood2 <- function(gamma.,gamma_known. = gamma_known){
      g_list <- seq(1,1, length.out = x_cate)
      ccc1 <- delta*exp(gamma.*y)
      ccc2 <- (1-delta)
      for(kkkk in 0:x_cate-1){
        g_list[kkkk+1] <- sum(ccc1[x1==kkkk])/sum(ccc2[x1==kkkk])
      }
      g_list2 <- apply(apply(x, 2, function(x) x*g_list),2,sum)
      ####g_list2 <- exp(beta+beta2*(x_transform-1.3)*(x_transform-1.3)+beta3*(x_transform))
      pseudo_pi <- as.vector(g_list2/(g_list2+exp(gamma.*y)))
      ###pseudo_pi <- exp(1.7+0.4*x1)/(exp(1.7+0.4*x1)+exp(gamma.*y))
      g_list4 <- seq(1,1, length.out = x_cate)
      for(iii in 1:n){
        sum3 <- seq(1,1,length.out = mmmm)
        sum4 <- seq(1,1,length.out = mmmm)
        for(jjj in 1:mmmm){
          sum3[jjj] <- exp(gamma.*large_matrix[iii,jjj])*large_matrix[iii,jjj]*g_list2[iii]/(g_list2[iii]+exp(gamma.*large_matrix[iii,jjj]))
          sum4[jjj] <- exp(gamma.*large_matrix[iii,jjj])
        }
        g_list4[iii] <- mean(sum3)/mean(sum4)
      }
      aaa <- delta*(1-pseudo_pi)*y-(1-delta)*g_list4
      return(sum(aaa)*sum(aaa))
    }
    boku2 <- optimize(likelihood2,c(-2.0,2.0))
    gamma_known2 <- boku2$minimum
    print(gamma_known2)
  }
  gamma_estimate7[iteite] <- gamma_known2
  
  
}
###print(mean((gamma_estimate-true_gamma)*(gamma_estimate-true_gamma)))
print(mean((gamma_estimate2-true_gamma)*(gamma_estimate2-true_gamma)))
print(mean((gamma_estimate3-true_gamma)*(gamma_estimate3-true_gamma)))
print(mean((gamma_estimate4-true_gamma)*(gamma_estimate4-true_gamma)))
print(mean((gamma_estimate5-true_gamma)*(gamma_estimate5-true_gamma)))
print(mean(gamma_estimate2-true_gamma))
print(mean(gamma_estimate3-true_gamma))
print(mean(gamma_estimate4-true_gamma))
print(mean(gamma_estimate5-true_gamma))
print(mean(gamma_estimate6-true_gamma))
print(mean(gamma_estimate7-true_gamma))
