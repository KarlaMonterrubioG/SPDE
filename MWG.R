#mwg ar ex2


# Libraries ---------------------------------------------------------------
library("MASS")
library(limSolve)
library(Matrix)
library(fields)

# Load data ---------------------------------------------------------------

load("D:/Hierarchical models + MCMC/Codes/forpaper/sin/data_sin.RData")
set.seed(100)

# Discretisation ----------------------------------------------------------

N <- 350
N.extended.domain <- 40 # This is for one-side
N.total <- N + 2*N.extended.domain

x<- seq(0,8,length.out = N)
h=x[2]-x[1]

try<-seq((0-h)*100,(8+h)*10,by=h)
x2<-try[-(1:( which(try==0)-N.extended.domain-1))]
x2<-x2[1:(N.extended.domain+350+N.extended.domain )]

theory <- matrix(0, nrow=M, ncol=N.total)

NM_tmp<-(N-1)/(M-1) 

for(i in 1:M){
  theory[i, NM_tmp*i + N.extended.domain] <- 1
}

index.theory<-seq(N.extended.domain+1,N+N.extended.domain,by=1)


# Hyperprior --------------------------------------------------------------

EXPO_L.extension <- function(xi=x,l=exp(l_hyperprior[1]),h.=h, N.extension = N.extended.domain){
  n1    <- length(xi)
  L <- matrix(0,n1+2*N.extension,n1+2*N.extension)
  alpha <- 1
  d <- h./alpha/l
  e <- l/alpha/h.
  a0 <- (sqrt(d)+sqrt(d+4*e))/sqrt(2)
  a1 <- (sqrt(d)-sqrt(d+4*e))/sqrt(2)
  for (j in 1:(n1+2*N.extension-1)){
    L[j,j:(j+1)] <- 0.5*c(a0,a1)
  } 
  L[n1+2*N.extension,n1+2*N.extension] <- 1/sqrt(alpha)
  return(list(L,a0))
}



EXPO_L2<-function(xi=x,l=exp(l_hyperprior[1]),h.=h){
  n1    <- length(xi)
  alpha <- 1
  d <- h./alpha/l
  e <- l/alpha/h.
  s1<-sqrt(d)
  s2<-sqrt(d+4*e)
  s3<-sqrt(2)
  a0 <- (s1+s2)/s3
  a1 <- (s1-s2)/s3
  diagonal<-c(rep(0.5*a0,(n1-1)),1/sqrt(alpha))
  return(list(diagonal,a0,0.5*a1))
}

# MCMC initialisation -----------------------------------------------------
nIterations <- 100000



sigma_proposal <- rep(1, N.total)
P              <- sigma_proposal*diag(1,N.total)
counter1 <- rep(0,N.total)
counter2 <- rep(0,N.total)
matrixaccpetance <- matrix(NA,nIterations, N.total)
accel <- c()
cont3 <- 0
v_save <- matrix(NA, nrow=N.total, ncol=nIterations+1)
u_save <- matrix(NA, nrow=N.total, ncol=nIterations+1)
u_new  <- vector(length=N.total)
u_old  <- vector(length=N.total)
A_old <- matrix(0,N.total, N.total)    
delta <- 0.01


l_hyperprior    <- c()
l_hyperprior[1] <- log(.5)


A_old_hp<-EXPO_L.extension(xi=x,l=exp(l_hyperprior[1]), h.=h, N.extended.domain)[[1]]
u_save[,1] <- solve(A_old_hp,rnorm(N.total))
u_old      <- u_save[,1]
tmp_llu_new<- exp(u_save[,1])
a0<-EXPO_L.extension(xi=x,l=exp(l_hyperprior[1]), h.=h, N.extended.domain)[[2]]
SEQ<-seq(N.total+1, N.total*N.total+1,N.total+1)
l.abd=matrix(0,2,N.total)
l.abd.temp=matrix(0,2,N.total)
l.abd[1,2:N.total]<-A_old_hp[SEQ]
l.abd[2,]<-diag(A_old_hp)

sdhyper <- 2
mag    <- c()
mag[1:(nIterations+1)] <- log(1)


# Prior initialisation -----------------------------------------------------

for (j in 1:N.total){
  A_old[j,(j-2):(j)%%N.total+1]<-(c(0,1,0)-c(1,-2,1)*tmp_llu_new[j]^2/h^2)*sqrt(h/ (exp(mag[1])*tmp_llu_new[j]))
}

A_old[1,N.total] <- 0
A_old[N.total,1] <- 0

#  -----------------------------------------------------

sigma    <- c()
sigma[1] <- 0.02
sd1      <- 1
accepted <- c()
accepted_hyper <- c()

m<-length(meas)
zero_n<-rep(0, N.total)

cont1 <- 0


length_bb<-N.total+m
v_save[,1] <- qr.solve(rbind(theory/sqrt(sigma[1]), A_old), (c(meas/sqrt(sigma[1]),zero_n)-rnorm(length_bb)))


acvecnoise <- c()
acvec_hyper<- c()
cont_hyper <- 0


######################MCMC
start.time<-Sys.time()
for (i in 1:nIterations){
  
  
  ################################ noise variance
  
  
  tz<-v_save[index.theory,i] 
  m_tz<-meas-tz
  quad<-sum((m_tz)^2)
  llik.old<- -m/2*log(sigma[i])-.5/sigma[i]*quad
  
  lsigma_new <- log(sigma[i])+rnorm(1, mean=0, sd=sd1)  
  llik.new<- -m/2*log(exp(lsigma_new))-.5/exp(lsigma_new)*quad
  alpha1<-llik.new-lsigma_new^2/(2*3)-llik.old+log(sigma[i])^2/(2*3)
  
  if(log(runif(1))<alpha1){
    accepted<-accepted+1
    sigma[i+1]<-exp(lsigma_new)
    acvecnoise[i]<-1
    llik.old<-llik.new
  } else {
    sigma[i+1]<-sigma[i]
    acvecnoise[i]<-0
  }
  
  
  if (i%%50==0) {
    cont1 <- cont1+1
    accenoise <- sum(acvecnoise[(seq(0, nIterations, by=50)[cont1]):i])
    if((accenoise/50)>.44){
      sd1 <- exp(log(sd1)+delta)
    } 
    else {
      sd1 <- exp(log(sd1)-delta)
    }
  }
  
  
  
  
  v_save[,i+1] <- qr.solve(rbind(theory/(sqrt(sigma[i+1])), A_old),
                           (c(meas/(sqrt(sigma[i+1])),zero_n)-rnorm(length_bb)))
  
  #STEP 2 M-H
  u_new <- u_save[,i]+sigma_proposal*rnorm(N.total) # mvrnorm(1,rep(0,N), P)
  
  
  for(j in 1:N.total){
    tmp_u_new     <- u_old
    tmp_u_new[j] <- u_new[j]
    
    temp_alpha<-A_old_hp[(j-2):j%%N.total+1,(j-2):j%%N.total+1]
    
    log_hyperprior_old <- -0.5*sum(( temp_alpha%*%u_old[(j-2):j%%N.total+1])^2)
    log_hyperprior_new <- -0.5*sum(( temp_alpha%*%tmp_u_new[(j-2):j%%N.total+1])^2)
    
    tmp_llu_new <- exp(tmp_u_new)
    if (j == 1){
      kernelw <- (c(0,1,0)-c(0,-2,1)*tmp_llu_new[j]^2/h^2)*sqrt(h/(exp(mag[i+1])*tmp_llu_new[j]))
    } else if (j == N.total) {
      kernelw <- (c(0,1,0)-c(1,-2,0)*tmp_llu_new[j]^2/h^2)*sqrt(h/(exp(mag[i+1])*tmp_llu_new[j]))
    } else {
      kernelw <- (c(0,1,0)-c(1,-2,1)*tmp_llu_new[j]^2/h^2)*sqrt(h/(exp(mag[i+1])*tmp_llu_new[j]))
    }
    
    log_prior_new <- -0.5*sum((kernelw%*%v_save[(j-2):j%%N.total+1,i+1])^2)
    log_prior_old <- -0.5*sum((A_old[j,(j-2):j%%N.total+1]%*%v_save[(j-2):j%%N.total+1,i+1])^2)
    
    A_new<- matrix(c(A_old[(j-2)%%N.total+1 ,(((j-2):j))%%N.total+1], kernelw, 
                     A_old[(j)%%N.total+1 ,(((j-2):j))%%N.total+1]),3,3,byrow=T)
    
    
    norm_cost      <- A_old[j ,(((j-2):j))%%N.total+1]%*%solve(A_new)
    log_norm_const <- log(norm_cost[2])   
    
    
    log_old     <- log_prior_old+log_hyperprior_old
    log_new     <- log_prior_new+log_hyperprior_new
    exponential <- (log_new-log_old)-log_norm_const  
    
    if(exponential>=0){
      u_old[j] <- u_new[j]
      A_old[j, ((j-2):j)%%N.total+1]<- kernelw
      counter1[j] <- counter1[j]+1 
      matrixaccpetance[i,j]<-1
    } 
    else {
      r=runif(1)
      if (r<exp(exponential)) {
        u_old[j] <- u_new[j]
        A_old[j, ((j-2):j)%%N.total+1]<- kernelw
        counter2[j] <- counter2[j]+1
        matrixaccpetance[i,j] <- 1
      }   
      else { 
        matrixaccpetance[i,j]<-0
      }
    }
  }
  
  
  if (i%%50==0) {
    accepted_l <- counter1+counter2
    cont3 <- cont3+1
    accel <- colSums(matrixaccpetance[(seq(0, nIterations, by=50)[cont3]):i,])
    
    sigma_proposal<- ifelse((accel/50)>.44,exp(log(sigma_proposal)+rep(delta,N.total)),exp(log(sigma_proposal)-rep(delta,N.total)) )
    
  }
  
  
  
  u_save[, i+1] <- u_old
  
  ###l_hyperprior with adaptive RW
  
  
  l_hyperprior_new<-(l_hyperprior[i])+rnorm(1, mean=0, sd=sdhyper) 
  
  
  
  
  L_lambda_temp<-EXPO_L2(x2,l=exp(l_hyperprior_new),h.=h)
  l.abd.temp[1,2:N.total]<-L_lambda_temp[[3]]
  l.abd.temp[2,]<-L_lambda_temp[[1]]
  a0_temp<-L_lambda_temp[[2]]
  
  llik.old <-(N.total-1)*log(.5*a0) - 
    0.5*sum((l.abd[2,]*u_save[,i+1]+ c(l.abd[1,2:N.total]*u_save[2:N.total,i+1],0) )^2)
  
  llik.new  <- (N.total-1)*log(.5*a0_temp) - 
    0.5*sum((l.abd.temp[2,]*u_save[,i+1]+c(l.abd.temp[1,2:N.total]*u_save[2:N.total,i+1],0))^2)
  
  alpha2 <- llik.new-l_hyperprior_new^2/(2*3)-llik.old+(l_hyperprior[i])^2/(2*3)
  
  if (log(runif(1)) < alpha2) {
    accepted_hyper    <- accepted_hyper+1
    l_hyperprior[i+1] <- l_hyperprior_new
    l.abd<-l.abd.temp
    acvec_hyper[i]    <- 1
    a0<-a0_temp
    diag(A_old_hp)<- l.abd.temp[2,]
    A_old_hp[SEQ]<-l.abd.temp[1,2:N.total]
  } 
  else {
    l_hyperprior[i+1]<-l_hyperprior[i]
    acvec_hyper[i]<-0
  }
  
  
  if (i%%50==0) {
    cont_hyper<-cont_hyper+1
    accepted_hyper<-sum(acvec_hyper[(seq(0, nIterations, by=50)[cont_hyper]):i])
    
    if((accepted_hyper/50)>.44){
      sdhyper <-exp(log(sdhyper)+delta)
    } else {
      sdhyper<-exp(log(sdhyper)-delta)
    }
    
  }
  
  
  print(i)
  
  
  
  
}


end.time<-Sys.time()
time.taken<-end.time-start.time
time.taken





