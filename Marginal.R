#Ex2 Marginal with domain extension

# Libraries ---------------------------------------------------------------
library("MASS")
library("fields")
library("scales")
library(limSolve)

# Load data ---------------------------------------------------------------

load("D:/Hierarchical models + MCMC/Codes/forpaper/sin/data_sin.RData")
set.seed(100)



# Discretisation ----------------------------------------------------------
N=350
M=350
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

ttheory<-t(theory)
index.theory<-seq(N.extended.domain+1,N+N.extended.domain,by=1)

# Hyperprior --------------------------------------------------------------
EXPO_L<-function(xi=x,l=exp(l_hyperprior[1]),h.=h){
  n1    <- length(xi)
  L <- matrix(0,n1,n1)
  alpha <- 1
  d <- h./alpha/l
  e <- l/alpha/h.
  a0 <- (sqrt(d)+sqrt(d+4*e))/sqrt(2)
  a1 <- (sqrt(d)-sqrt(d+4*e))/sqrt(2)
  for (j in 1:(n1-1)){
    L[j,j:(j+1)] <- 0.5*c(a0,a1)
  }
  L[n1,n1] <- 1/sqrt(alpha)
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


# Likelihood --------------------------------------------------------------
LTLfast=function(l.abd, N=nrow(A_old)){
  l2.abd=matrix(0,5,N)
  l2.abd[3,]=colSums(l.abd^2)
  dudl=l.abd[2,-N]*l.abd[1,-1]+l.abd[2,-1]*l.abd[3,-N]
  l2.abd[2,]=c(0,dudl)
  l2.abd[4,]=c(dudl,0)
  ul=l.abd[1,c(-1,-2)]*l.abd[3,c(-(N-1),-N)]
  l2.abd[1,]=c(0,0,ul)
  l2.abd[5,]=c(ul,0,0)
  return(l2.abd)
}

uppdiag<-seq(N+1, N.total*N.total+1,N.total+1)
lowdiag<-seq(2,N.total*(N.total-1),N.total+1)
l.abd.prior<-matrix(0,3,N.total)
l.abd.prior.temp<-matrix(0,3,N.total)
Lik<- function(l2.abd, s, index.theory, ATy, yTy){
  aux=Solve.banded(l2.abd, nup=2, nlow=2, B = ttheory)
  aux=aux[index.theory,]
  diag(aux)=diag(aux)+s
  l2.abd[3,index.theory]=l2.abd[3,index.theory]+s^-1
  aux2=Solve.banded(l2.abd, nup=2, nlow=2, B = ATy)
  ss=-.5/s*(yTy-sum(s^-1*ATy*aux2))
  llik.new<- -m/2*log(2*pi)-.5*determinant(aux,logarithm=T)$modulus[1] +ss
  llik.new}



# MCMC initialisation -----------------------------------------------------
nIterations<-100000
u_save<-matrix(NA, nrow=length(x2), ncol=nIterations+1)
sigma<-c()
sigma[1]<-.02

l_hyperprior<-c()
l_hyperprior[1]<-log(.5)

A_old_hp<-EXPO_L(xi=x2,l=exp(l_hyperprior[1]))[[1]]
inv_l<-solve(A_old_hp)
CC<-inv_l%*%t(inv_l)
u_save[,1]<-mvrnorm(1, rep(0, N.total), CC)
tmp_llu_new<-exp(u_save[,1])
a0<-EXPO_L(xi=x,l=exp(l_hyperprior[1]))[[2]]
l.abd=matrix(0,2,N.total)
l.abd.temp=matrix(0,2,N.total)
uppdiag<-seq(N+1, N.total*N.total+1,N.total+1)
l.abd[1,2:N.total]<-A_old_hp[uppdiag]
l.abd[1,1]<-0
l.abd[2,]<-diag(A_old_hp)




A_old<-matrix(0,N.total, N.total)  

mag<-c()
mag[1:(nIterations+1)]<-log(1)




for(j in 1:N.total){
  A_old[j, ((j-2):j%%N.total)+1]<-(c(0,1,0)-c(1,-2,1)*tmp_llu_new[j]^2/h^2)*sqrt(h/ (exp(mag[1])*tmp_llu_new[j]))
}


A_old[1,N.total]<-0
A_old[N.total,1]<-0
A_new<-matrix(0,N.total, N.total)
A_TEMP<-matrix(0,N.total, N.total)
delta<-.01
I<-diag(1, N.total)
cont1<-0
sd1<-1
sd2<-2
sdhyper<-2.5
accepted<-c()
accepted_hyper<-c()

acvecnoise<-c()
acvec_hyper<-c()
cont_hyper=0
cont2<-0

eta<-matrix(NA, nrow=length(x2), ncol=nIterations+1)
eta[,1]<-solve(t(chol(CC+.00001*diag(1, length(x2)))))%*%(u_save[,1])


m<-length(meas)
zeromean<-rep(0, m)
Im<-diag(1,m)
zeromeanN<-rep(0, N.total)
In<-diag(1,N.total)


l.abd.prior[1,]<-c(0,A_old[uppdiag])
l.abd.prior[2,]<-diag(A_old)
l.abd.prior[3,]<-c(A_old[lowdiag],0)
l2.abd<-LTLfast(l.abd.prior)
yTy=sum(meas^2)
ATy=matrix(0,N.total,1)
ATy[index.theory]=meas
llik.old<-Lik(l2.abd, s=sigma[1], index.theory, ATy, yTy)
counter<-rep(0,nIterations)

# For_loop ----------------------------------------------------------------

start.time<-Sys.time()
for (i in 1:nIterations){
  lsigma_new<-log(sigma[i])+rnorm(1, mean=0, sd=sd1)  
  
  
  llik.new<-Lik(l2.abd, s=exp(lsigma_new), index.theory, ATy, yTy)
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
    cont1<-cont1+1
    accenoise<-sum(acvecnoise[(seq(0, nIterations, by=50)[cont1]):i])
    if((accenoise/50)>.44){
      sd1 <-exp(log(sd1)+delta)
    } else {
      sd1<-exp(log(sd1)-delta)
    }
  }
  
  
  
  v_ess<-mvrnorm(1,zeromeanN,In) 
  unif<-runif(1)
  
  loglikelihood<-llik.old
  ly<-loglikelihood+log(unif)
  
  theta<-runif(1,0,2*pi)
  theta_min<-theta-(2*pi)
  theta_max<-theta
  
  
  while(TRUE){
    
    etanew<-eta[,i]*cos(theta)+v_ess*sin(theta)
    
    
    tmp_u_new<-Solve.banded(l.abd, nup=1, nlow=0,B=etanew)
    tmp_llu_new<-exp(tmp_u_new)
    
    
    
    alpha<-sapply(tmp_llu_new, function(x) x^2/h^2)
    beta<-sapply(tmp_llu_new, function(x) sqrt(h/((exp(mag[i+1]))*x)))
    ab<-alpha*beta
    l.abd.prior[3,]<-c((-ab)[2:N.total], 0)
    l.abd.prior[2,]<-beta+2*ab
    l.abd.prior[1,]<-c(0, (-ab)[1:(N.total-1)])
    
    
    l2.abd<-LTLfast(l.abd.prior)
    
    
    loglikelihood_prime<-Lik(l2.abd, s=sigma[1+i], index.theory, ATy, yTy)
    counter[i]<-counter[i]+1
    if(loglikelihood_prime>ly){break}
    
    if(theta<0){
      theta_min<-theta
    } else {
      theta_max<-theta
    }
    theta<-runif(1,theta_min,theta_max)
  }
  
  
  eta[,i+1]<-etanew
  u_save[, i+1]<-tmp_u_new
  
  
  
  llik.old<-loglikelihood_prime
  l_hyperprior_new<-l_hyperprior[i]+rnorm(1, mean=0, sd=sdhyper) 
  
  L_lambda_temp<-EXPO_L2(x2,l=exp(l_hyperprior_new))
  l.abd.temp[1,2:N.total]<-L_lambda_temp[[3]]
  l.abd.temp[1,1]<-0
  l.abd.temp[2,]<-L_lambda_temp[[1]]
  a0_temp<-L_lambda_temp[[2]]
  tmp_u_new<-c(Solve.banded(l.abd.temp, nup=1, nlow=0,B=eta[,i+1]))
  tmp_llu_new<-exp(tmp_u_new)
  
  
  
  alpha<-sapply(tmp_llu_new, function(x) x^2/h^2)
  beta<-sapply(tmp_llu_new, function(x) sqrt(h/((exp(mag[i+1]))*x)))
  ab<-alpha*beta
  l.abd.prior.temp[3,]<-c((-ab)[2:N.total], 0)
  l.abd.prior.temp[2,]<-beta+2*ab
  l.abd.prior.temp[1,]<-c(0, (-ab)[1:(N.total-1)])
  
  
  
  l2.abd.temp<-LTLfast(l.abd.prior.temp)
  llik.new<- Lik(l2.abd.temp, s=sigma[1+i], index.theory, ATy, yTy)
  
  
  alpha2<-llik.new-(l_hyperprior_new^2/(2*3))-llik.old+((l_hyperprior[i])^2/(2*3))
  
  if(log(runif(1))<alpha2){
    accepted_hyper<-accepted_hyper+1
    l_hyperprior[i+1]<-l_hyperprior_new
    acvec_hyper[i]<-1
    l.abd<-l.abd.temp
    l.abd.prior<-l.abd.prior.temp
    l2.abd<-l2.abd.temp
    u_save[, i+1]<-tmp_u_new
    llik.old<-llik.new
  } else {
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
  #if(i%%100==0) {print(i)}
}

end.time<-Sys.time()
time.taken<-end.time-start.time
time.taken



