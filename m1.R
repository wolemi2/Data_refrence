#paper review##04/11/2018
#############################November version
#setwd("C:\\Users\\user\\Desktop\\CALIBRATION_REVIEW\\November")
library(doParallel)
library(LaplacesDemon)
#library(dlm)
#library(MCMCpack)
source("rwish.R",echo=FALSE)
library(MASS)
library(emulator)
library(invgamma)
library(compiler)
library(abind)
library(Matrix)
library(parallel)
#load("IDYNO_train")
#load("NUFEB_train")
#load("design")
nam=c("sub","o2","nh4","muhet","muaob","Yhet","Yaob")
N=97#100
nout=4
ninp=7
T=96
dess=c(30e-3,10e-3,6e-3,0.00006944444,0.00001158333,0.63,0.24)##
############################################read in data
gap=gap1=list()
for(j in 1:N){
gap[[j]]=as.matrix(read.csv(sub("j",j,paste("data","input_j.csv",sep="/")),header=TRUE))
gap1[[j]]=as.matrix(read.csv(sub("j",j,paste("data","output_j.csv",sep="/")),header=TRUE))
}
input=abind(gap,along=3)
y=abind(gap1,along=3)
#paddle training data
temp1=input
temp2=y
for(i in 1:T){
temp1[i,,]=input[96,,,drop=FALSE]
temp2[i,,]=y[96,,,drop=FALSE]
}
input0=abind(input,temp1,along=1)
y0=abind(y,temp2,along=1)
input=input0
y=y0
#####################scale
min1=apply(y,2,min)
max1=apply(y,2,max)
Y=y
for(j in 1:nout){
Y[,j,]=(y[,j,]-min1[j])/(max1[j]-min1[j])
}
arx1=Y
#Y=Y+rnorm(96*97*4,.00001,.0001)
###############initialize
arx=Y
X=F=input[1,,]
#X=F=(rbind((input[1,,]),arx[(1),,]))
TT=96
ndes=97
ntime=150#96
r=dim(y)[3]
s=97
p=2##order
q=4##nuber of output
ninp=ninp#+q
lambda=.01
m0 = matrix(1,nrow=ninp,ncol=q)
C0 = 15*diag(1,ninp)
n0 = 20
d0 =diag(20,q)
dis=c(.5,.5)

################################
tvarm1=function(input,arx,p,dis,m0,C0,d0,n0,psi,betai,Re){
T=ntime
r=ndes
bv=dis[1]; bw=dis[2]
m=list();d=list();n=list();C=list();h=list();s=list()
length(s)=length(n)=length(m)=length(C)=T
A=R=Q=e=d=list()
m[[1]]=m0
C[[1]]=C0
d[[1]]=d0                                                                                
n[[1]]=n0
h[[1]]=n0+q-1
s[[1]]=d0/n0
for(t in p:(T)){
F=input[t,,]
#F=(rbind((input[t,,]),arx[(t-1),,]))
Re=corr.matrix(t(F),scales=betai,method=1)+diag(lambda,r)
G=diag(psi)
R[[t]]=G%*%C[[t-1]]%*%t(G)/bw
Q[[t]]=(t(F)%*%R[[t]]%*%F)+ Re#[1,1]
A[[t]]=R[[t]]%*%F%*%solve(Q[[t]]) 
at=G%*%m[[t-1]]
e[[t]]=t(arx[t,,])-t(F)%*%at#G%*%m[[t-1]]
m[[t]]=at+A[[t]]%*%e[[t]]
d[[t]]=bv*d[[t-1]]+(t(e[[t]])%*%solve(Q[[t]])%*%e[[t]])
n[[t]]=(n[[t-1]])+ndes
#h[[t]]=n[[t]]+q+1
h[[t]]=(bv*h[[t-1]])+ndes
C[[t]]=R[[t]]-A[[t]]%*%Q[[t]]%*%t(A[[t]])
C[[t]]=(C[[t]]+t(C[[t]]))/2##make it a symmetric matrix
}
d[[1]]=d[[p]]=d[[(p+1)]]
n[[1]]=n[[p]]=n[[(p+1)]]
h[[1]]=h[[p]]=h[[(p+1)]]
nn=abind(n,along=1)
hh=abind(h,along=1)
dd=d#abind(d,along=0)
list(dd=dd,nn=nn,hh=hh)
}

myV=function(hh,dd,p,dis){
VV=list()
bv=dis[1]
T=length(hh)
ht=hh#unlist(nn)+q-1
V=list();length(V)=T
r0=solve(dd[[T]])
r0=(r0+t(r0))/2
VT=rwish(ht[T],(r0))#inverse variance
VT=(VT+t(VT))/2
VV[[T]]=(VT)##variance
for(t in (T-1):(p)){
#ht=nn[[t]]+q-1
VT=(bv*VV[[t+1]])+rwish((1-bv)*ht[t],solve(dd[[t]]))	
VT=(VT+t(VT))/2
VV[[t]]=(VT)	
}
VV[[1]]=VV[[p]]=VV[[(p+1)]]
V=lapply(VV,solve)
v1=abind(V,along=0)
return(v1)
}

tvarm2=function(input,arx,p,dis,m0,C0,Vi,psi,betai,Re){
T=ntime
r=ndes
bv=dis[1]
bw=dis[2]
Theta=list()#array(m0,dim=c(ninp,q,T))
Err=array(0,dim=c(T,q,ndes))
G=diag(psi)
###forward filtering
ee=list()
aa=list()
C=list()
m=list()
W=list()
length(m)=length(C)=T
mt=m0
Ct=C0
for(t in (p):(T)){
F=(input[t,,])
#F=(input[t,,])
Re=corr.matrix(t(F),scales=betai,method=1)+diag(lambda,r)
R=G%*%Ct%*%t(G)/bw
W[[t]]=R*(1-bw)#1/((1-bw)*R)
Q=(t(F)%*%R%*%F)+ Re#[1,1]
A=R%*%F%*%solve(Q) 
at=G%*%mt
e=t(arx[t-1,,])-t(F)%*%at
mt=(at)+A%*%e
m[[t]]=mt
Ct=(R-A%*%Q%*%t(A))
Ct=(Ct+t(Ct))/2
C[[t]]=Ct#as.matrixnorm(nearPD(Ct)$mat)
ee[[t]]=e
aa[[t]]=at
}
aa[[1]]=(aa[[p]])
ee[[1]]=(ee[[p]])
C[[1]]=(C[[p]])#
m[[1]]=(m[[p]])#
W[[1]]=W[[p]]
CC=C#abind(C,along=3)
mm=abind(m,along=3)
WW=W#abind(W,along=3)
###backward sampling
vr=(Vi[T,,]+t(Vi[T,,]))/2
cc=(CC[[T]]+t(CC[[T]]))/2
Theta[[T]]=rmatrixnorm(mm[,,T],as.positive.definite(cc),as.positive.definite(vr))
for(t in (T-1):(p)){
mt=((1-bw)*mm[,,t])+(bw*Theta[[(t+1)]])##posterior
if(t<TT){
Err[t+1,,]=(arx[t+1,,])-t(Theta[[t+1]])%*%F
}else{
Err[t+1,,]=(arx[96,,])-t(Theta[[t+1]])%*%F
}	
Ct=(1-bw)*(CC[[t]]+t(CC[[t]]))/2
vr=(Vi[t,,]+t(Vi[t,,]))/2
Theta[[t]]=rmatrixnorm(mt,as.positive.definite(Ct),as.positive.definite(vr))
}
W0=abind(W,along=3)
C00=abind(CC,along=3)#[[T]]
m00=mm#[,,T]
Err[1,,]=Err[p,,]=Err[(p+1),,]
Theta[[1]]=Theta[[p]]=Theta[[(p+1)]]
Thetai=abind(Theta,along=0)
Erri=Err
list(Thetai=Thetai,Erri=Erri,WW=WW,Ct=C00,mt=m00,W0=W0)
}
#res3=tvarm2(input,arx,p,dis,m0,C0,Vi,psi)
#Thetai=res3$Thetai
#W0=res3$WW
####################psi AR
psifun=function(Thetai,W0,psi0,tau0){
theta=Thetai#abind(Thetai,along=0)
ww0=abind(W0,along=3)
#Psi=rep(NA,ninp)
#for (t in 1:T)
kk=1
T=ntime
#w=1/(ww0[,,T])## G%*%theta_precision
w=1/apply(ww0,1,rowSums)
tau= 1/((1/tau0) + diag(diag(w)) %*% crossprod(theta[-(T),,1]))
ps=diag(tau * (diag(diag(w)) %*% t(theta[-(T),,kk]) %*% theta[-1,,kk]) + psi0/tau0)
Psi=mvrnorm(1,ps,diag(diag(tau)))
return(Psi)
}
#psi=psifun(Thetai,W0,psi0,tau0)
######## Metropolis algorithm ################
MH=function(betai,input,Erri,Vi,p){
T=ntime
r=ndes
LogPrior=0 
dd=length(betai)
for(i in 1:dd){
LogPrior=LogPrior-(betai[i]/4)-0.9*(1-exp(-betai[i]/4))
}
#############
LogLikeB=0
llkk=rep(NA,T)
#for(t in p:T){
LogLike=foreach(t=p:T,.combine=sum,.packages="emulator",.export=c("lambda","ndes","q"))%do%{
F=(input[t,,])
#F=(rbind((input[t,,]),arx[(1),,]))
ee=corr.matrix(t(F),scales=betai,method=1)+diag(lambda,r)
SigmaInv=solve(ee)
Sigmainv=(SigmaInv + t(SigmaInv))/2
LogLikeA = -1.0*q*log(det(ee))/2
#LogLikeB=0
#LogLikeC=foreach(t=p:T,.combine=sum,.packages="emulator",.export=c("lambda","ndes","q"))%do%{
#for(t in p:T){
vt=Vi[t,,]
LogLikeB =LogLikeA -1.0*(ndes)*log(det(vt))/2-tr(((Erri[t,,])%*%Sigmainv%*%t(Erri[t,,])%*%solve(vt))/2)
llkk[t]=LogLikeB
}
#LogLike=LogLikeB+LogLikeA
#LogLike=sum(llkk[-1])
#########################
LogJ=sum(log(betai))
score=LogLike+LogPrior+LogJ
return(score)
}

MH2=function(betai,input,Erri,Vi,p,scorei,w,ww,sig){
#acc=list()
T=nrow(Erri)
#proposal function for correlation betai
betabeta = log(betai) + sig*rnorm(length(betai),0,1)
beta=exp(betabeta)
#alpha=alphai
score=MH(beta,input,Erri,Vi,p)
acceptance=score-scorei
if(log(runif(1,0,1))< acceptance){
beta=beta;score=score
#print("accept")
rate=0
}else{
beta=betai;score=scorei
#print("reject")
rate=1
}
list(beta=beta,score=score,rate=rate)
}

MH=cmpfun(MH)
MH2=cmpfun(MH2)
tvarm1= cmpfun(tvarm1)
tvarm2= cmpfun(tvarm2)
psifun=cmpfun(psifun)
cl <- makeCluster(4)
registerDoParallel(cl)
############initialize the fitting
Gibbs=function(arx,ninp,p,dis,m0,C0,d0,n0,input,NIter){ 
NIter =10000#00#00#1100
burn =500
thin=5
T=ntime
r=ndes
sig=.1
# initialization of parameters
WWt=array(0,c(ninp,ninp,T,NIter/thin))#list()
CCt=array(0,c(ninp,ninp,T,NIter/thin))#list()
MMt=array(0,c(ninp,q,T,NIter/thin))#list()
THETA=array(0,c(T,ninp,q,NIter/thin))
V=array(0,c(T,q,q,NIter/thin))
ERR=array(0,c(T,q,ndes,NIter/thin))
Beta=array(1,c(ninp,NIter/thin))
betai=rep(1,ninp)#c(6.0779326,4.1288388,0.2816155,28.2552092,0.6274502,1.5451606,1.1944421,1)
psi0 = 0.1; tau0 = 0.5#1
#psi <- rnorm(ninp,psi0,tau0)
psi <- mvrnorm(1,rep(psi0,ninp),diag(tau0,ninp))
PSI= array(0,c(ninp,NIter/thin))
# initialization of calibration parameters
ww=100
RATE=array(NA,c(ww,NIter/thin))
rates=rep(NA,ww)##Compute acceptance rate
for(k in 1:NIter){
h=k
ind1=k%%thin
n=h[!ind1]
print(k)
#F=(rbind((input[1,,]),arx[(1),,]))
#F=(input[1,,])
#Re=corr.matrix(t(F),scales=betai,method=1)+diag(lambda,r)
#forward filtering with unknown variance
res=tvarm1(input,arx,p,dis,m0,C0,d0,n0,psi,betai,Re)###
##sample variance
Vi=myV(res$hh,res$dd,p,dis)###produce [V]
###forward filtering/ backward sampling with estimated variance
ress=tvarm2(input,arx,p,dis,m0,C0,Vi,psi,betai,Re)#
V[,,,n/thin]=Vi
THETA[,,,n/thin]=Thetai=ress$Thetai
ERR[,,,n/thin]=Erri=ress$Erri
W0=ress$WW
WWt[,,,n/thin]=ress$W0
CCt[,,,n/thin]=ress$Ct
MMt[,,,n/thin]=ress$mt
##compute psi. AR coefficients
psi=psifun(Thetai,W0,psi0,tau0)
PSI[,n/thin]=psi
####Metropolis algorithm
scorei=MH(betai,input,Erri,Vi,p)##produces [scorei]
for(w in 1:ww){
hh=MH2(betai,input,Erri,Vi,p,scorei,w,ww,sig)
betai=hh$beta;scorei=hh$score
rates[w]=hh$rate
}
Beta[,n/thin]=hh$beta
RATE[,n/thin]=unlist(rates)
#CCt[[n/thin]]=ress$Ct
#MMt[[n/thin]]=ress$mt
}  
list(V,THETA,ERR,Beta,PSI,RATE,CCt,MMt,WWt) 
}  

Gibbs=cmpfun(Gibbs)   	
multGP1=Gibbs(arx,ninp,p,dis,m0,C0,d0,n0,input,NIter)
path=getwd()
save(multGP1,file=paste(path,"multGP1",sep="/"))
GP=multGP1
length(which(GP[[6]][20,]==0))/ncol(GP[[6]])
length(which(GP[[6]][,]==0))/(ncol(GP[[6]])*nrow(GP[[6]]))
stopCluster(cl)
#############
N2=11
ytest=inptest1=list()
for(j in 1:N2){
ytest[[j]]=as.matrix(read.csv(sub("j",j,paste("data","ytest_j.csv",sep="/")),header=TRUE))
inptest1[[j]]=as.matrix(read.csv(sub("j",j,paste("data","inptest_j.csv",sep="/")),header=TRUE))
}
#F0=(rbind((input[1,,]),arx[(1),,]))
F0=(input[1,,])
#b1=min1
#b2=(max1)###not maximum deciation
burn=500
nsim=1000
design=t(X)
input3=inptest1
y2=ytest
V0=apply(GP[[1]][,,,-(1:burn)],c(2,3),rowMeans)#96x2x2	
THETA0=apply(GP[[2]][,,,-(1:burn)],c(2,3),rowMeans)#96x7x2
ERR0=apply(GP[[3]][,,,-(1:burn)],c(2,3),rowMeans)#96x2x75
BETA0=apply(GP[[4]][,-(1:burn)],1,mean)
PSI0=apply(GP[[5]][,-(1:burn)],1,mean)
C0=apply(GP[[7]][,,,-(1:burn)],c(1,2),rowMeans)
M0=apply(GP[[8]][,,,-(1:burn)],c(1,2),rowMeans)
W0=apply(GP[[9]][,,,-(1:burn)],c(1,2),rowMeans)
#6.0779326  4.1288388  0.2816155 28.2552092  0.6274502  1.5451606  1.1944421
#6.0999550  4.1138221  0.2783141 28.2141863  0.6046648  1.5631418  1.2033580

########################
TT=ntime#96
pred=function(newinput,design){
T=min(nrow(newinput),ntime)
newoutput1=array(NA,c(T,q))
newoutput2=array(NA,c(T,q,q))
temp=rep(NA,nout)
newoutput2[1,,]=matrix(0,nout,nout)
newoutput1[1,]=(ini[1,]*(max1-min1))+min1
#r=1
for(t0 in (p):T){
F=(input[t,,])
Re=corr.matrix(t(F),scales=BETA0,method=1)+diag(lambda,r)
newdata0=newinput[t0,,drop=FALSE]
SigmaInv=solve(Re)
rho=corr.matrix((newdata0),t(F),scales=BETA0)
err1=t(ERR0[t0,,])
coef1=(THETA0[t0,,])
VV=V0[t0,,]
mut = (newdata0[,]%*%(coef1)) + (t(rho)%*%SigmaInv%*%(err1))
#mut = (newdata0[,]%*%(coef1)) + (t(rho)%*%SigmaInv%*%(err1))
sigma2t = VV*diag(diag(1,1) - t(rho)%*% SigmaInv %*% rho)
newoutput2[t0,,]= sigma2t
for(j in 1:nout){temp[j]=(mut[,j]*(max1[j]-min1[j]))+min1[j]}
newoutput1[t0,]=temp
#ini=mut
}
return(list(newoutput1,newoutput2))
}

ytest2=ytest
for(k in 1:N2){
for(j in 1:nout){
ytest2[[k]][,j]=(ytest[[k]][,j]-min1[j])/(max1[j]-min1[j])
}}
####
vv0=vv=lam=list()
dev=list()##variance
P=RMSE=matrix(NA,11,q)
for(k in 1:length(y2)){
newinput=inptest1[[k]]##validation
newoutput=ytest[[k]]
#ini=newoutput[8,,drop=FALSE]
ini=arx1[,,1]
vv[[k]]=pred(newinput,(design))
vv0[[k]]=emu=vv[[k]][[1]]
lammp=y2[[k]]
T=min(nrow(newinput),ntime)
lammp=lammp[1:T,]
lam[[k]]=lammp
dev[[k]]=vv[[k]][[2]]
##proportion of variance explained
P[k,]=1-(colSums((lammp-emu)^2)/colSums((lammp- colMeans(lammp))^2))
RMSE[k,]=sqrt(colMeans((lammp-emu)^2))##
}
P
RMSE
r1=1;r2=2;r3=3;r4=4
q1=(abind(vv0,along=1))
q2=(abind(lam,along=1))##lammp
#q2=(q2-b1)/(b2-b1)
#q3=na.omit(cbind(q1,q2))
1-(colSums((q2-q1)^2)/colSums((q2- colMeans(q2))^2))
sqrt(colMeans((q2-q1)^2))##
mean(ergMean(GP[[1]][96,r1,r1,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r2,r2,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r3,r3,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r4,r4,-c(1:burn)]))
