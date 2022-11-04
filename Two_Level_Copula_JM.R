### define functions ###

lambda01=function(t,x){exp(Beta1*x)*(1/(2+t))}
lambda02=function(t,x){exp(Beta2*x)*(1/(4+t))}
#lambda01=1
#lambda02=1

Lambda01=function(t, x){integrate(lambda01, x=x,lower = 0, upper = t, subdivisions=2000)$value}
Lambda02=function(t, x){integrate(lambda02, x=x,lower = 0, upper = t, subdivisions=2000)$value}

S01=function(t, x){exp(-Lambda01(t,x))}
S02=function(t, x){exp(-Lambda02(t,x))}

F01=function(t, x){1-S01(t,x)}
F02=function(t, x){1-S02(t,x)}

f01=function(t, x){lambda01(t,x)*S01(t,x)}
f02=function(t, x){lambda02(t,x)*S02(t,x)}

# Clayton's copula:

S12=function(t,x){(S01(t,x)^(-Theta)+S02(t,x)^(-Theta)-1)^(-1/Theta)}

ft1=function(t,x){(S01(t,x)^(-Theta)+S02(t,x)^(-Theta)-1)^(-1/Theta-1)*S01(t,x)^(-Theta-1)*f01(t,x)}
ft2=function(t,x){(S01(t,x)^(-Theta)+S02(t,x)^(-Theta)-1)^(-1/Theta-1)*S02(t,x)^(-Theta-1)*f02(t,x)}

Ft1=function(t,x){integrate(Vectorize(ft1), x=x,lower = 0, upper = t, subdivisions=2000)$value}
Ft2=function(t,x){integrate(Vectorize(ft2), x=x,lower = 0, upper = t, subdivisions=2000)$value}

Ft=function(t,x){Ft1(t,x)+Ft2(t,x)}

Fy_ij=function(yij,Qt,rrho,BetaT,wij,zij,bi,fy_ij){
  return(pALD(yij,mu=(BetaT%*%wij+zij%*%bi)[1],sigma=rrho,p=Qt,lower.tail = TRUE)-fy_ij)
}

YIJ=function(Qt,rrho,BetaT,wij,zij,bi,fy_ij){
  return(uniroot(Fy_ij, c(-1000, 1000), tol = 0.0001, Qt=Qt,rrho=rrho,BetaT=BetaT,wij=wij,zij=zij,bi=bi,fy_ij=fy_ij)$root)
}




library(rstan)
library(survival)
library(copula)
library(splines2)
library(mvtnorm)
library(ald)
library(quantreg)
library(statmod)
#library(dclone)

###  competing risks with copula simulations

i=1
set.seed(i*10002)


#### data generations ####

n=500 # number of subjects
nc=1

TAU=0.5
Theta=2
Beta1=1.5 # parameter betaU true value
Beta2=-0.5

K=4 # 4 continuous, will add disparate in future works.

RhoT1=0.30
RhoT2=0.25
RhoT3=0.20
RhoT4=0.35

Rho12=0.20
Rho13=0.15
Rho14=0.40

Rho23=0.35
Rho24=0.50

Rho34=0.60

Omega=matrix(c(1, RhoT1, RhoT2, RhoT3, RhoT4,
               RhoT1, 1, Rho12, Rho13, Rho14,
               RhoT2, Rho12, 1, Rho23, Rho24,
               RhoT3, Rho13, Rho23, 1, Rho34,
               RhoT4, Rho14, Rho24, Rho34, 1),ncol = 5)

Rho=Omega[lower.tri(Omega, diag = FALSE)]


CovT=matrix(c(1), ncol = 1)
CovY=matrix(c(1, Rho12, Rho13, Rho14,
              Rho12, 1, Rho23, Rho24,
              Rho13, Rho23, 1, Rho34,
              Rho14, Rho24, Rho34, 1),ncol = 4)
CovYT=c(RhoT1, RhoT2, RhoT3, RhoT4)

CCovY=CovY-CovYT%*%CovT%*%CovYT


Qt1=0.25
BetaT1=c(1,0.5)
rrho1=1
SiG1=matrix(c(0.4,0.04,0.04,0.3), ncol=2)  # cov of Bi

Qt2=0.5
BetaT2=c(1.5,1)
rrho2=1
SiG2=matrix(c(0.4,0.07,0.07,0.3), ncol=2)  # cov of Bi


Qt3=0.75
BetaT3=c(1,1.5)
rrho3=1
SiG3=matrix(c(0.5,0.15,0.15,0.4), ncol=2)  # cov of Bi

Qt4=0.9
BetaT4=c(1.5,2)
rrho4=1
SiG4=matrix(c(0.5,0.2,0.2,0.4), ncol=2)  # cov of Bi

BetaT=rbind(BetaT1,BetaT2,BetaT3,BetaT4)
rrho=c(rrho1,rrho2,rrho3,rrho4)
SiG=c(SiG1[lower.tri(SiG1, diag = TRUE)], SiG2[lower.tri(SiG2, diag = TRUE)], SiG3[lower.tri(SiG3, diag = TRUE)], SiG4[lower.tri(SiG4, diag = TRUE)])



ID=seq(1,n)
X=matrix(rnorm(n,0,1), nrow = n, ncol = 1) # risk factor values

# generate copula
clayton.cop=claytonCopula(Theta,dim = 2)
FF=rCopula(n,clayton.cop)
S1=FF[,1]
S2=FF[,2]


# generate times

#U=-log(S1)/exp(Beta1*X)
#V=-log(S2)/exp(Beta2*X)

U=c()
for(i in 1:n){
  x=X[i,1]
  s=S1[i]
  u0=0
  u=0.01
  while(u-u0>0.001){
    u0=u
    u=u0-(1/(lambda01(u0,x)))*(integrate(lambda01, x=x,lower = 0, upper = u0, subdivisions=2000)$value+log(s))
  }
  U=as.vector(c(U,u))
}



V=c()
for(i in 1:n){
  x=X[i,1]
  s=S2[i]
  v0=0
  v=0.01
  while(v-v0>0.001){
    v0=v
    v=v0-(1/(lambda02(v0,x)))*(integrate(lambda02, x=x,lower = 0, upper = v0, subdivisions=2000)$value+log(s))
  }
  V=as.vector(c(V,v))
}

C=15*rbeta(n,3,5)

T=c()
D=c()
for(i in 1:n){
  Ti=min(U[i],V[i],C[i])
  T=c(T,Ti)
  if(Ti==U[i]){
    D=c(D,1)
  }else{
    if(Ti==V[i]){
      D=c(D,2)
    }else{
      D=c(D,0)
    }
  }
  
}


FTT=c()
for(i in 1:n){
  ti=T[i]
  # ti=min(U[i],V[i])
  xi=X[i,1]
  FTT=c(FTT,Ft(ti,xi))
}

QFTT=qnorm(FTT)


## use different columns for same variables of each outcome is just to make it expanable for future work with different observation times.
## or missing values in some longitudinal outcomes.
ID_Y1=c()
ID_Y2=c()
ID_Y3=c()
ID_Y4=c()

T1ij=c()
T2ij=c()
T3ij=c()
T4ij=c()

W1ij=c()
W2ij=c()
W3ij=c()
W4ij=c()

Z1ij=c()
Z2ij=c()
Z3ij=c()
Z4ij=c()

Y1ij=c()
Y2ij=c()
Y3ij=c()
Y4ij=c()

T1i=c()
T2i=c()
T3i=c()
T4i=c()

X1i=c()
X2i=c()
X3i=c()
X4i=c()

N1i=c()
N2i=c()
N3i=c()
N4i=c()



for (i in 1:length(FTT)) {
  #nob1=sample(c(5,6,7,8,9,10,11,12,13,14,15),1,prob = c(0.02,0.03,0.05,0.1,0.2,0.2,0.1,0.1,0.1,0.05,0.05))
  #nob1=sample(c(1,2,3,4),1,prob = c(0.1,0.3,0.4,0.2))  
  #nob2=sample(c(5,6,7,8,9,10,11,12,13,14,15),1,prob = c(0.02,0.03,0.05,0.1,0.2,0.2,0.1,0.1,0.1,0.05,0.05))
  #nob2=nob1
  nob=sample(c(5,6,7,8,9,10,11,12,13,14,15),1,prob = c(0.02,0.03,0.05,0.1,0.2,0.2,0.1,0.1,0.1,0.05,0.05))
  
  N1i=c(N1i,nob)
  N2i=c(N2i,nob)
  N3i=c(N3i,nob)
  N4i=c(N4i,nob)
  
  ID_Y1=c(ID_Y1,rep(ID[i],nob))
  t1i=rep(T[i],nob)
  T1i=c(T1i,t1i)
  X1i=rbind(X1i,do.call(rbind, replicate(nob, X[i,], simplify=FALSE)))
  t1ij=sort(runif(nob,0,T[i]))
  T1ij=c(T1ij,t1ij)
  
  ID_Y2=c(ID_Y2,rep(ID[i],nob))
  t2i=rep(T[i],nob)
  T2i=c(T2i,t2i)
  X2i=rbind(X2i,do.call(rbind, replicate(nob, X[i,], simplify=FALSE)))
  #t2ij=sort(runif(nob,0,T[i]))
  t2ij=t1ij
  T2ij=c(T2ij,t2ij)
  
  ID_Y3=c(ID_Y3,rep(ID[i],nob))
  t3i=rep(T[i],nob)
  T3i=c(T3i,t3i)
  X3i=rbind(X3i,do.call(rbind, replicate(nob, X[i,], simplify=FALSE)))
  #t3ij=sort(runif(nob,0,T[i]))
  t3ij=t1ij
  T3ij=c(T3ij,t3ij)
  
  ID_Y4=c(ID_Y4,rep(ID[i],nob))
  t4i=rep(T[i],nob)
  T4i=c(T4i,t4i)
  X4i=rbind(X4i,do.call(rbind, replicate(nob, X[i,], simplify=FALSE)))
  #t4ij=sort(runif(nob,0,T[i]))
  t4ij=t1ij
  T4ij=c(T4ij,t4ij)
  
  
  w1ij=cbind(rep(X[i,1],nob),t1ij)
  W1ij=rbind(W1ij,w1ij)
  z1ij=cbind(rep(1,nob),t1ij)
  Z1ij=rbind(Z1ij,z1ij)
  b1i=rmvnorm(1,mean = rep(0,length(z1ij[1,])),sigma = SiG1)[1,]
  
  w2ij=cbind(rep(X[i,1],nob),t2ij)
  W2ij=rbind(W2ij,w2ij)
  z2ij=cbind(rep(1,nob),t2ij)
  Z2ij=rbind(Z2ij,z2ij)
  b2i=rmvnorm(1,mean = rep(0,length(z2ij[1,])),sigma = SiG2)[1,]
  
  w3ij=cbind(rep(X[i,1],nob),t3ij)
  W3ij=rbind(W3ij,w3ij)
  z3ij=cbind(rep(1,nob),t3ij)
  Z3ij=rbind(Z3ij,z3ij)
  b3i=rmvnorm(1,mean = rep(0,length(z3ij[1,])),sigma = SiG3)[1,]
  
  w4ij=cbind(rep(X[i,1],nob),t4ij)
  W4ij=rbind(W4ij,w4ij)
  z4ij=cbind(rep(1,nob),t4ij)
  Z4ij=rbind(Z4ij,z4ij)
  b4i=rmvnorm(1,mean = rep(0,length(z4ij[1,])),sigma = SiG4)[1,]
  
  
  #QFyy=rmvnorm(max(nob1,nob2), mean = CovYT%*%CovT%*%QFTT[i], sigma = CCovY)
  QFyy=rmvnorm(nob, mean = CovYT%*%CovT%*%QFTT[i], sigma = CCovY)
  
  Y1_ij=c()
  for (j in 1:nob) {
    QFy=QFyy[j,1]
    fy_ij=pnorm(QFy)
    yij=YIJ(Qt1,rrho1,BetaT1,w1ij[j,],z1ij[j,],b1i,fy_ij)
    Y1_ij=c(Y1_ij,yij)
  }
  Y1ij=c(Y1ij,Y1_ij)
  
  
  Y2_ij=c()
  for (j in 1:nob) {
    QFy=QFyy[j,2]
    fy_ij=pnorm(QFy)
    yij=YIJ(Qt2,rrho2,BetaT2,w2ij[j,],z2ij[j,],b2i,fy_ij)
    Y2_ij=c(Y2_ij,yij)
  }
  Y2ij=c(Y2ij,Y2_ij)
  
  
  Y3_ij=c()
  for (j in 1:nob) {
    QFy=QFyy[j,3]
    fy_ij=pnorm(QFy)
    yij=YIJ(Qt3,rrho3,BetaT3,w3ij[j,],z3ij[j,],b3i,fy_ij)
    Y3_ij=c(Y3_ij,yij)
  }
  Y3ij=c(Y3ij,Y3_ij)
  
  
  Y4_ij=c()
  for (j in 1:nob) {
    QFy=QFyy[j,4]
    fy_ij=pnorm(QFy)
    yij=YIJ(Qt4,rrho4,BetaT4,w4ij[j,],z4ij[j,],b4i,fy_ij)
    Y4_ij=c(Y4_ij,yij)
  }
  Y4ij=c(Y4ij,Y4_ij)
  
}



ID=rep(ID,nc)+rep(n*seq(0,nc-1),each=length(ID))
T=rep(T,nc)
D=rep(D,nc)
X=do.call(rbind, replicate(nc, X, simplify=FALSE))

T1=subset(cbind(T,D),D==1)[,1]
T2=subset(cbind(T,D),D==2)[,1]


ID_Y1=rep(ID_Y1,nc)+rep(n*seq(0,nc-1),each=length(ID_Y1))
T1y=rep(T1ij,nc)
T1i=rep(T1i,nc)
X1i=do.call(rbind, replicate(nc, X1i, simplify=FALSE))
W1=do.call(rbind, replicate(nc, W1ij, simplify=FALSE))
Z1=do.call(rbind, replicate(nc, Z1ij, simplify=FALSE))

Y1=rep(Y1ij,nc)

ID_Y2=rep(ID_Y2,nc)+rep(n*seq(0,nc-1),each=length(ID_Y2))
T2y=rep(T2ij,nc)
T2i=rep(T2i,nc)
X2i=do.call(rbind, replicate(nc, X2i, simplify=FALSE))
W2=do.call(rbind, replicate(nc, W2ij, simplify=FALSE))
Z2=do.call(rbind, replicate(nc, Z2ij, simplify=FALSE))

Y2=rep(Y2ij,nc)


ID_Y3=rep(ID_Y3,nc)+rep(n*seq(0,nc-1),each=length(ID_Y3))
T3y=rep(T3ij,nc)
T3i=rep(T3i,nc)
X3i=do.call(rbind, replicate(nc, X3i, simplify=FALSE))
W3=do.call(rbind, replicate(nc, W3ij, simplify=FALSE))
Z3=do.call(rbind, replicate(nc, Z3ij, simplify=FALSE))

Y3=rep(Y3ij,nc)

ID_Y4=rep(ID_Y4,nc)+rep(n*seq(0,nc-1),each=length(ID_Y4))
T4y=rep(T4ij,nc)
T4i=rep(T4i,nc)
X4i=do.call(rbind, replicate(nc, X4i, simplify=FALSE))
W4=do.call(rbind, replicate(nc, W4ij, simplify=FALSE))
Z4=do.call(rbind, replicate(nc, Z4ij, simplify=FALSE))

Y4=rep(Y4ij,nc)




ID_Y=c(ID_Y1, ID_Y2, ID_Y3, ID_Y4)
Ty=c(T1y, T2y, T3y, T4y)
TI=c(T1i, T2i, T3i, T4i)
Xi=rbind(X1i, X2i, X3i, X4i)
W=rbind(W1, W2, W3, W4)
Z=rbind(Z1, Z2, Z3, Z4)
Y=c(Y1, Y2, Y3, Y4)
Qt=c(Qt1, Qt2, Qt3, Qt4)

Ny = c(sum(N1i), sum(N2i), sum(N3i), sum(N4i))




p1=length(X[1,])
p2=length(W1[1,])
p3=length(Z1[1,])



rate0=sum(as.numeric(D==0))/length(D)

rate1=sum(as.numeric(D==1))/length(D)

rate2=sum(as.numeric(D==2))/length(D)



N_L=10 # number of baseline parameters
N_D=3  # degree for Spline methods

# N_Q=seq(1,N_L-N_D)/(N_L-N_D+1) # quantile location of knots for Spline methods

# B=bSpline(T, knots = quantile(T,probs = N_Q), Boundary.knots = c(0,max(T)), degree = N_D)
# IB=ibs(T, knots = quantile(T,probs = N_Q), Boundary.knots = c(0,max(T)), degree = N_D)
# M=mSpline(T, knots = quantile(T,probs = N_Q), Boundary.knots = c(0,max(T)), degree = N_D)
# IM=iSpline(T, knots = quantile(T,probs = N_Q), Boundary.knots = c(0,max(T)), degree = N_D)


CPU=c()
CPV=c()
N1=length(T1)
N2=length(T2)
T11=sort(T1)
T22=sort(T2);
CPU[1]=0
CPV[1]=0
for (i in 2:N_L){
  CPU[i]=T11[(N1/N_L)*(i-1)]
  CPV[i]=T22[(N2/N_L)*(i-1)]
  #CPU[i]=max(T1)/N_L*(i-1)
  #CPV[i]=max(T2)/N_L*(i-1)
}

LUU=c()
LVV=c()

for (i in 1:(N_L-1)){
  LUU[i]=(log(2+CPU[i+1])-log(2+CPU[i]))/(CPU[i+1]-CPU[i])
  LVV[i]=(log(4+CPV[i+1])-log(4+CPV[i]))/(CPV[i+1]-CPV[i])
}
LUU[N_L]=(log(2+max(T1))-log(2+CPU[N_L]))/(max(T1)-CPU[N_L])
LVV[N_L]=(log(4+max(T2))-log(4+CPV[N_L]))/(max(T1)-CPU[N_L])


N_GQ=10

GQ_ti=gauss.quad(N_GQ,kind="legendre")$nodes
GQ_wi=gauss.quad(N_GQ,kind="legendre")$weights




## prepare data for stan fit##

Data=list(
  N=length(T),
  K=K,
  N1=length(T1),
  N2=length(T2),
  N_L=N_L,
  N_Y=length(Y),
  Ny=Ny[1],
  Ni=N1i,
  N_GQ=N_GQ,
  p1=p1,
  p2=p2,
  p3=p3,
  ID=ID,
  T=T,
  D=D,
  X=X,
  T1=T1,
  T2=T2,
  #B=B,
  #IB=IB,
  #M=M,
  #IM=IM,
  ID_Y=ID_Y,
  Ty=Ty,
  TI=TI,
  Xi=Xi,
  W=W,
  Z=Z,
  Y=Y,
  Qt=Qt,
  CPU=CPU,
  CPV=CPV,
  GQ_ti=GQ_ti,
  GQ_wi=GQ_wi
)




## initial values:

dat=data.frame(cbind(T,D,X))
dat1=subset(dat,dat[,2]!=2)
dat2=subset(dat,dat[,2]!=1)
dat2$D=as.numeric(dat2$D!=0)
beta1=c()
for (j in 1:p1) {
  beta1=c(beta1,unname(coxph(Surv(dat1$T,dat1$D) ~ subset(dat1,select = -c(T,D))[,j])$coef))
}

beta2=c()
for (j in 1:p1) {
  beta2=c(beta2,unname(coxph(Surv(dat2$T,dat2$D) ~ subset(dat2,select = -c(T,D))[,j])$coef))
}


lu=rep(1,Data$N_L)
lv=rep(1,Data$N_L)


k0=1
k1=0
betaT=c()
for (i in 1:K) {
  k1=k1+Ny[i]
  Qt_k=Qt[i]
  Y_k=Y[k0:k1]
  W_k=W[k0:k1,]
  
  betaT=rbind(betaT,unname(rq(Y_k~W_k-1,tau = Qt_k)$coef))
  
  k0=k1+1
}



b=matrix(rep(0,length(T)*K*p3),nrow = length(T)*K,ncol = p3)



L_omega=diag(1,K+1,K+1)

sig=array(rep(diag(1,p3,p3),K), dim = c(K,p3,p3))

sig[1,,]=diag(1,p3,p3)
sig[2,,]=diag(1,p3,p3)
sig[3,,]=diag(1,p3,p3)
sig[4,,]=diag(1,p3,p3)



init00=list(
  beta1=as.array(beta1+rnorm(length(beta1),0,0.05)),
  beta2=as.array(beta2+rnorm(length(beta2),0,0.05)),
  theta=1+0.05,
  LU=lu+rnorm(length(lu),0,0.05),
  LV=lv+rnorm(length(lu),0,0.05),
  L_omega=L_omega,
  betaT=betaT+matrix(rnorm(K*p2,0,0.05), ncol = p2),
  rrho=c(rrho1, rrho2,rrho3, rrho4)+rnorm(K,0,0.05),
  sig=sig,
  b=b
)  

init01=list(
  beta1=as.array(beta1+rnorm(length(beta1),0,0.05)),
  beta2=as.array(beta2+rnorm(length(beta2),0,0.05)),
  theta=1+0.05,
  LU=lu+rnorm(length(lu),0,0.05),
  LV=lv+rnorm(length(lu),0,0.05),
  L_omega=L_omega,
  betaT=betaT+matrix(rnorm(K*p2,0,0.05), ncol = p2),
  rrho=c(rrho1, rrho2,rrho3, rrho4)+rnorm(K,0,0.05),
  sig=sig,
  b=b
)  



init0=list(init00,init01)




#######################################
######### Fit the stan model ##########

fit=stan(file = "Two_Level_Copula_JM.stan",data = Data, warmup = 10000, iter = 20000, init = init0, chains = 2, cores =2, refresh=10)

round(summary(fit, pars = c("beta1", "beta2", "theta", "Ktau", "LU", "LV", "Rho", "betaT", "rrho", "sig", "lp__", "DEV", "DEV_J"), probs = c(0.025, 0.975))$summary,digits = 2)








# 
# 
# # latent event times:
# E_B1=mean(extract(fit)$beta1)
# SD_B1=sd(extract(fit)$beta1)
# CP_B1=as.numeric(quantile(extract(fit)$beta1, 0.025)<Beta1 & quantile(extract(fit)$beta1, 0.975)>Beta1)
# 
# E_B2=mean(extract(fit)$beta2)
# SD_B2=sd(extract(fit)$beta2)
# CP_B2=as.numeric(quantile(extract(fit)$beta2, 0.025)<Beta2 & quantile(extract(fit)$beta2, 0.975)>Beta2)
# 
# E_T=mean(extract(fit)$theta)
# SD_T=sd(extract(fit)$theta)
# CP_T=as.numeric(quantile(extract(fit)$theta, 0.025)<Theta & quantile(extract(fit)$theta, 0.975)>Theta)
# 
# E_Ktau=mean(extract(fit)$Ktau)
# SD_Ktau=sd(extract(fit)$Ktau)
# CP_Ktau=as.numeric(quantile(extract(fit)$Ktau, 0.025)<TAU & quantile(extract(fit)$Ktau, 0.975)>TAU)
# 
# 
# # Longitudinal outcomes:
# E_betaT=c()
# SD_betaT=c()
# CP_betaT=matrix(0, nrow = K, ncol = p2)
# E_rrho=colMeans(extract(fit)$rrho)
# SD_rrho=unname(sapply(data.frame(extract(fit)$rrho), sd))
# CP_rrho=c()
# for (k in 1:K){
#   E_betaT=rbind(E_betaT,colMeans(extract(fit)$betaT[,k,]))
#   
#   SD_betaT=rbind(SD_betaT, unname(sapply(data.frame(extract(fit)$betaT[,k,]), sd)))
#   
#   for (j in 1:p2){
#     CP_betaT[k,j]=as.numeric(quantile(extract(fit)$betaT[,k,j], 0.025)<BetaT[k,j] & quantile(extract(fit)$betaT[,k,j], 0.975)>BetaT[k,j])
#   }
#   
#   CP_rrho=c(CP_rrho, as.numeric(quantile(extract(fit)$rrho[,k], 0.025)<rrho[k] & quantile(extract(fit)$rrho[,k], 0.975)>rrho[k]))
#   
# }
# 
# 
# # Conditional normal correlations:
# 
# Rho_fit=c()
# for (i in 1:length(extract(fit)$Rho[,1,1])) {
#   Rho_fit=rbind(Rho_fit, extract(fit)$Rho[i,,][lower.tri(extract(fit)$Rho[i,,], diag = FALSE)])
# }
# 
# E_rho=colMeans(Rho_fit)
# SD_rho=unname(sapply(data.frame(Rho_fit), sd))
# CP_rho=c()
# for (i in 1:length(E_rho)) {
#   CP_rho=c(CP_rho, as.numeric(quantile(Rho_fit[,i], 0.025)<Rho[i] & quantile(Rho_fit[,i], 0.975)>Rho[i]))
# }
# 
# 
# 
# # Random effects:
# 
# sig_fit=c()
# for (k in 1:K) {
#   sig_k=c()
#   for (i in 1:length(extract(fit)$sig[,1,1,1])) {
#     sig_ik=extract(fit)$sig[i,k,,]
#     sig_k=rbind(sig_k, sig_ik[lower.tri(sig_ik, diag = TRUE)])
#   }
#   sig_fit=cbind(sig_fit, sig_k)
# }
# 
# E_sig=colMeans(sig_fit)
# SD_sig=unname(sapply(data.frame(sig_fit), sd))
# CP_sig=c()
# for (i in 1:length(E_sig)) {
#   CP_sig=c(CP_sig, as.numeric(quantile(sig_fit[,i], 0.025)<SiG[i] & quantile(sig_fit[,i], 0.975)>SiG[i]))
# }
# 
# True_para=c(Beta1, Beta2, Theta, TAU, c(t(BetaT)), rrho, Rho, SiG)
# True_para
# 
# #write.table(t(c(E_B1, E_B2, E_T,E_Ktau,c(t(E_betaT)), E_rrho, E_rho, E_sig)),"Para_C_J.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# #write.table(t(c(SD_B1, SD_B2, SD_T,SD_Ktau,c(t(SD_betaT)), SD_rrho, SD_rho, SD_sig)),"SD_C_J.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# #write.table(t(c(CP_B1, CP_B2, CP_T,CP_Ktau,c(t(CP_betaT)), CP_rrho, CP_rho, CP_sig)),"CP_C_J.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# 
# 
# 
# ##########  LU LV  ####
# 
# E_LU=c()
# SD_LU=c()
# CP_LU=c()
# E_LV=c()
# SD_LV=c()
# CP_LV=c()
# 
# for(i in 1:N_L){
#   E_LU=c(E_LU,mean(extract(fit)$LU[,i]))
#   SD_LU=c(SD_LU,sd(extract(fit)$LU[,i]))
#   CP_LU=c(CP_LU,as.numeric(quantile(extract(fit)$LU[,i], 0.025)<LUU[i] & quantile(extract(fit)$LU[,i], 0.975)>LUU[i]))
#   E_LV=c(E_LV,mean(extract(fit)$LV[,i]))
#   SD_LV=c(SD_LV,sd(extract(fit)$LV[,i]))
#   CP_LV=c(CP_LV,as.numeric(quantile(extract(fit)$LV[,i], 0.025)<LVV[i] & quantile(extract(fit)$LV[,i], 0.975)>LVV[i]))
# }
# 
# #write.table(t(c(E_LU, E_LV)),"Para_C_LULV.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# #write.table(t(c(SD_LU, SD_LV)),"SD_C_LULV.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# #write.table(t(c(CP_LU, CP_LV)),"CP_C_LULV.txt",append = TRUE, row.names = FALSE, col.names = FALSE)
# 
# 
# ###############
# 
# 
