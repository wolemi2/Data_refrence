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
dis=c(.9,.9)

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
multGP3=Gibbs(arx,ninp,p,dis,m0,C0,d0,n0,input,NIter)
path=getwd()
save(multGP3,file=paste(path,"multGP3",sep="/"))
GP=multGP3
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
#####
TT=ntime#96
pred=function(newinput,design){
T=min(nrow(newinput),ntime)
newoutput1=array(NA,c(T,q))
newoutput2=array(NA,c(T,q))
temp=rep(NA,nout)
temp2=temp
newoutput2[1,]=matrix(0,nout)
newoutput1[1,]=(ini[1,]*(max1-min1))+min1
#r=1
for(t0 in (p):T){
F=(input[t0,,])
Re=corr.matrix(t(F),scales=BETA0,method=1)+diag(lambda,r)
newdata0=newinput[t0,,drop=FALSE]
SigmaInv=solve(Re)
rho=corr.matrix((newdata0),t(F),scales=BETA0)
err1=t(ERR0[t0,,])
coef1=(THETA0[t0,,])
VV=V0[t0,,]
mut = (newdata0[,]%*%(coef1)) + (t(rho)%*%SigmaInv%*%(err1))
sigma2t = diag(VV*diag(diag(1,1) - t(rho)%*% SigmaInv %*% rho))

for(j in 1:nout){
temp[j]=(mut[,j]*(max1[j]-min1[j]))+min1[j]
temp2[j]=(sqrt(sigma2t[j])*(max1[j]-min1[j])^2)
}
newoutput1[t0,]=temp
newoutput2[t0,]=temp2# abind(temp2,along=0)
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
1-(colSums((q2-q1)^2)/colSums((q2- colMeans(q2))^2))
sqrt(colMeans((q2-q1)^2))##
mean(ergMean(GP[[1]][96,r1,r1,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r2,r2,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r3,r3,-c(1:burn)]))
mean(ergMean(GP[[1]][96,r4,r4,-c(1:burn)]))
#######mvrnorm generation
nsim=1000
SIM=list()
for(k in 1:length(vv)){
ss=list()
resb=vv0[[k]]
res2a=dev[[k]]
for(i in 1:nrow(resb)){
mm=mvrnorm(nsim,resb[i,],Sigma=(diag(res2a[i,]))^2)
#ss[[i]]=c(colMeans(mm),apply(mm,2,sd))
ss[[i]]=c(colMeans(mm),((res2a[i,])))
}
SIM[[k]]=abind(ss,along=0)
}
#########Figure 5
###################################PLOTS
ind1=1:150
ind2=1:150
lap=c(5,10)
par(mar = c(4, 5, 3, 1) + 0.1, cex =1.2)
par(mfrow=c(2,2))
#for(i in lap){
i=5#10
ind=ind2
j=1
dev1=(SIM[[i]][,j])-2*((SIM[[i]][,j+4]))
dev2=(SIM[[i]][,j])+2*((SIM[[i]][,j+4]))
dev1=dev1[ind];dev2=dev2[ind]
y0=(y2[[i]])[ind,j]
mmu=SIM[[i]][ind,j]
ttime=0.5*(1:length(mmu))
plot(ttime,y0,"l",xlab="Time (hr)",lwd=2,ylab="Biomass (kg/m^3)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,80),ylim=c(-10,30))
lines(ttime,mmu,col="green",lwd=2)
lines(ttime,dev1,col="red",ylab="",xlab="",lwd=2,lty=1)
lines(ttime,dev2,col="red",ylab="",xlab="",lwd=2,lty=1)
legend("topleft",c("Prediction","Simulation","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

j=2
dev1=(SIM[[i]][,j])-2*((SIM[[i]][,j+4]))
dev2=(SIM[[i]][,j])+2*((SIM[[i]][,j+4]))
dev1=dev1[ind];dev2=dev2[ind]
y0=(y2[[i]])[ind,j]
mmu=SIM[[i]][ind,j]#(SIM[[i]][,j]*(b2[j]-b1[j]))+b1[j]
ttime=0.5*(1:length(mmu))
plot(ttime,y0,"l",xlab="Time (hr)",lwd=2,ylab="Total no of particles",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(0,35),xlim=c(0,80))
lines(ttime,mmu,col="green",lwd=2)
lines(ttime,dev1,col="red",ylab="",xlab="",lwd=2,lty=1)
lines(ttime,dev2,col="red",ylab="",xlab="",lwd=2,lty=1)
legend("topleft",c("Prediction","Simulation","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

##height
j=3
dev1=(SIM[[i]][,j])-2*((SIM[[i]][,j+4]))
dev2=(SIM[[i]][,j])+2*((SIM[[i]][,j+4]))
dev1=dev1[ind];dev2=dev2[ind]
y0=(y2[[i]])[ind,j]
mmu=SIM[[i]][ind,j]#*(b2[1]-b1[1]))+b1[1]
plot(ttime,y0,"l",xlab="Time (hr)",lwd=2,ylab="Biofilm height (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(-15,9),xlim=c(0,80))
lines(ttime,mmu,col="green",lwd=2)
lines(ttime,dev1,col="red",ylab="",xlab="",lwd=2,lty=1)
lines(ttime,dev2,col="red",ylab="",xlab="",lwd=2,lty=1)
legend("topleft",c("Prediction","Simulation","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

j=4
#roughness
dev1=(SIM[[i]][,j])-2*((SIM[[i]][,j+4]))
dev2=(SIM[[i]][,j])+2*((SIM[[i]][,j+4]))
dev1=dev1[ind];dev2=dev2[ind]
y0=(y2[[i]])[ind,j]
mmu=SIM[[i]][ind,j]
plot(ttime,y0,"l",xlab="Time (hr)",lwd=2,ylab="Roughness (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(-25,25),xlim=c(0,80))
lines(ttime,mmu,col="green",lwd=2)
lines(ttime,dev1,col="red",ylab="",xlab="",lwd=2,lty=1)
lines(ttime,dev2,col="red",ylab="",xlab="",lwd=2,lty=1)
legend("topleft",c("Prediction","Simulation","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

###############################################################OTHER PLOTS

##Fig3a
##observation variance for some selected time points
r1=1;r2=2
tvec=c(1:96)
lap=c(1,50,96)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(3,3))
for(j in lap){
hh=GP[[1]][j,r1,r1,-(1:burn)]
foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ', sigma[11[a]]^2),list(a=tvec[j])),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
#foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ',sigma[a]^2),list(a=tvec[r])))
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
lines(density(GP[[1]][j,r1,r1,-(1:burn)]),col="red",xlab=sub("j",j,paste("V","j",sep="_")))
plot(ergMean(GP[[1]][j,r1,r1,-(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
plot(GP[[1]][j,r1,r1,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
}
####Figure 3b
tvec=c(1:96)
lap=c(1,50,96)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(3,3))
for(j in lap){
hh=GP[[1]][j,r2,r2,-(1:burn)]
foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ', sigma[22[a]]^2),list(a=tvec[j])),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
lines(density(GP[[1]][j,r2,r2,-(1:burn)]),col="red",xlab=sub("j",j,paste("V","j",sep="_")))
plot(ergMean(GP[[1]][j,r2,r2,-(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
plot(GP[[1]][j,r2,r2,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
}
#Fig S2
#covariance
tvec=c(1:96)
r1=1;r2=2
lap=c(1,50,96)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(3,3))
for(j in lap){
hh=GP[[1]][j,r1,r2,-(1:burn)]
foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ', sigma[12[a]]),list(a=tvec[j])))
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
lines(density(GP[[1]][j,r1,r2,-(1:burn)]),col="red",xlab=sub("j",j,paste("V","j",sep="_")))
ts.plot(ergMean(GP[[1]][j,r1,r2,-(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean")
ts.plot(GP[[1]][j,r1,r2,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots")
}
###########Fig S4
##beta
tvec=1:7
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(3,3))
for(j in 1:3){
hh=GP[[4]][j,-(1:burn)]
foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="", prob=T,main=substitute(paste('Density of ',beta[a]),list(a=tvec[j])),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
lines(density(GP[[4]][j,-(1:burn)]),col="red",xlab=sub("j",j,paste("Beta","j",sep="_")))
plot(ergMean(GP[[4]][j,-c(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
plot(GP[[4]][j,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
}

##########Fig 2
##beta
tvec=1:7
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(4,2))
for(j in 1:7){
hh=(GP[[4]][j,-c(1:burn)])
foo=hist(hh,breaks=10,xaxt="n",prob=TRUE,xlab=substitute(paste(beta[a]),list(a=tvec[j])),main="",col="green",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
}
#
##psi
tvec=1:7
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(4,2))
for(j in 1:7){
m1=(GP[[5]][j,-c(1:burn)])
hist(m1,prob=TRUE,xlab=substitute(paste(Psi[a]),list(a=tvec[j])),main="",col="green",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
}

###Fig S1
##Psi=1
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(3,3))
for(j in 1:3){
hh=GP[[5]][j,-(1:burn)]
foo=hist(hh,breaks=10,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ',psi[a]),list(a=tvec[j])),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=2),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
lines(density(GP[[5]][j,-(1:burn)]),col="red",xlab=sub("j",j,paste("Psi","j",sep="_")))
plot(ergMean(GP[[5]][j,-c(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
plot(GP[[5]][j,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
}

###Figure S3a
tt=1:96
tvec=1:7
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.5)
par(mfrow=c(4,2))
for(j in 1:7){
m1=apply(GP[[2]][tt,j,r1,-c(1:burn)],1,mean)
m2=(apply(GP[[2]][tt,j,r1,-c(1:burn)],1,sd))
#m1=apply(GP[[2]][j,-1,-c(1:burn)],1,quantile,prob=c(.50))
#m2=apply(GP[[2]][j,-1,-c(1:burn)],1,quantile,prob=c(.05))
m3=m1-m2;m4=m1+m2
yylim=range(m1,m2,m3,m4)
plot(0.5*(1:96),m1,type="l",col="green",xlab="Time (hr)",ylab=substitute(paste(theta[a]),list(a=tvec[j])),ylim=yylim,cex.axis=1.5,cex.lab=1.4,cex.main=1.5)
abline(0,0)
lines(0.5*(1:96),m3,col="red",lwd=2,lty=2,ylim=yylim,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(0.5*(1:96),m4,col="red",lwd=2,lty=2,ylim=yylim,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
}

###Figure S3b
tvec=1:7
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.5)
par(mfrow=c(4,2))
for(j in 1:7){
m1=apply(GP[[2]][tt,j,r2,-c(1:burn)],1,mean)
m2=(apply(GP[[2]][tt,j,r2,-c(1:burn)],1,sd))
#m1=apply(GP[[2]][j,-1,-c(1:burn)],1,quantile,prob=c(.50))
#m2=apply(GP[[2]][j,-1,-c(1:burn)],1,quantile,prob=c(.05))
m3=m1-m2;m4=m1+m2
yylim=range(m1,m2,m3,m4)
plot(0.5*(1:96),m1,type="l",col="green",xlab="Time (hr)",ylab=substitute(paste(theta[a]),list(a=tvec[j])),ylim=yylim,cex.axis=1.5,cex.lab=1.4,cex.main=1.5)
abline(0,0)
lines(0.5*(1:96),m3,col="red",lwd=2,lty=2,ylim=yylim,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(0.5*(1:96),m4,col="red",lwd=2,lty=2,ylim=yylim,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
}
#############################################Figure 1
gap=gap1=gap2=list()
for(j in 1:97){
#gap[[j]]=as.matrix(read.csv(sub("j",j,paste("data2","idynotest_j.csv",sep="/")),header=TRUE))
gap2[[j]]=as.matrix(read.csv(sub("j",j,paste("data2","idynotrain_j.csv",sep="/")),header=TRUE))
}
arx2=abind(gap2,along=3)
idynotest=list()
for(j in 1:11){
idynotest[[j]]=as.matrix(read.csv(sub("j",j,paste("data2","idynotest_j.csv",sep="/")),header=TRUE))
}
##############################
r1=2
library(RColorBrewer)
arx=y
arr1=arx[tt,r1,]
arr2=arx2[,r1,]
par(mfrow=c(2,2))
par(mfrow=c(2,2))
foo=hist(arr1,prob=T,breaks=seq(min(arr1),max(arr1),l=12),panel.first=grid(),xaxt="n",main="Simulation",cex=1.2,xlab="Total no of particles",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(4,15),ylim=c(0,.30))
axis(side=1,at=foo$mids,labels=round(seq(min(arr1),max(arr1),l=11)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
foo=hist(arr2,prob=T,breaks=seq(min(arr2),max(arr2),l=12),panel.first=grid(),xaxt="n",main="Observation",cex=1.2,xlab="Total no of particles",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(4,15),ylim=c(0,.30))
axis(side=1,at=foo$mids,labels=round(seq(min(arr2),max(arr2),l=11)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(4,15))
n=97
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:96)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr1[,1],col=sample(color,1),main="Simulation",ylab="Total no of particles",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(4,15),xlim=c(0,48))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr1[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	par(new=TRUE)
}
#
par(new=FALSE)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:72)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr2[,1],col=sample(color,1),main="Observation",ylab="Total no of particles",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(4,15),xlim=c(0,48))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr2[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	par(new=TRUE)
}

########################################biomass
r1=1
brk=seq(-6.5,4.5,l=10)
eps=.01
arr1=arx[tt,r1,]
arr2=arx2[,r1,]
par(mar = c(4, 5, 3, 1) + 0.1, cex = 1.5)
par(mfrow=c(2,2))
foo=hist(arr1,prob=T,breaks=max(arr1)-min(arr1),panel.first=grid(),xaxt="n",main="Simulation",cex=1.2,xlab=expression(paste("Biomass concentration"," ", (kgm^-3))),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(-6,4))
axis(side=1,at=foo$mids,labels=round(seq(min(arr1),max(arr1))),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(-6,4))

foo=hist(arr2,prob=T,breaks=max(arr2)-min(arr2),panel.first=grid(),xaxt="n",main="Observation",cex=1.2,xlab=expression(paste("Biomass concentration"," ", (kgm^-3))),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(-6,4))
axis(side=1,at=foo$mids,labels=round(seq(min(arr2),max(arr2))),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(-6,4))
n <- 97
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:96)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr1[,1],col=sample(color,1),main="Simulation",ylab=expression(paste("Biomass concentration"," ", (kgm^-3))),xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48),ylim=c(-6,4))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr1[,i],col=sample(color,1),xlab="",ylab=expression(paste("Biomass concentration"," ", (kgm^-3))),type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48))
	par(new=TRUE)
}
#
par(new=FALSE)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:72)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr2[,1],col=sample(color,1),main="Observation",ylab=expression(paste("Biomass concentration"," ", (kgm^-3))),xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48),ylim=c(-6,4))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr2[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48))
	par(new=TRUE)
}
############FigS6
#######################biofilm height
r1=3
arr1=arx[tt,r1,]
arr2=arx2[,r1,]
par(mfrow=c(2,2))
par(mfrow=c(2,2))
foo=hist(arr1,prob=T,breaks=seq(min(arr1),max(arr1),l=12),panel.first=grid(),xaxt="n",main="Simulation",cex=1.2,xlab="Biofilm height (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=round(seq(min(arr1),max(arr1),l=11)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
foo=hist(arr2,prob=T,breaks=seq(min(arr2),max(arr2),l=12),panel.first=grid(),xaxt="n",main="Observation",cex=1.2,xlab="Biofilm height (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=round(seq(min(arr2),max(arr2),l=11)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(4,15))
n=97
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:96)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr1[,1],col=sample(color,1),main="Simulation",ylab="Biofilm height (m)",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr1[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48))
	par(new=TRUE)
}
#
par(new=FALSE)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:72)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr2[,1],col=sample(color,1),main="Observation",ylab="Biofilm height (m)",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr2[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	par(new=TRUE)
}

#################roughness
r1=4
arr1=arx[tt,r1,]
arr2=arx2[,r1,]
par(mfrow=c(2,2))
par(mfrow=c(2,2))
foo=hist(arr1,prob=T,breaks=seq(min(arr1),max(arr1),l=9),panel.first=grid(),xaxt="n",main="Simulation",cex=1.2,xlab="Surface roughness (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=round(seq(min(arr1),max(arr1),l=8)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
foo=hist(arr2,prob=T,breaks=seq(min(arr2),max(arr2),l=9),panel.first=grid(),xaxt="n",main="Observation",cex=1.2,xlab="Surface roughness (m)",cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
axis(side=1,at=foo$mids,labels=round(seq(min(arr2),max(arr2),l=8)),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
n=97
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:96)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr1[,1],col=sample(color,1),main="Simulation",ylab="Surface roughness (m)",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr1[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48))
	par(new=TRUE)
}
#
par(new=FALSE)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ttime=0.5*(1:72)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(ttime,arr2[,1],col=sample(color,1),main="Observation",ylab="Surface roughness (m)",xlab="Time (hr)",type="l",lwd=1,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,xlim=c(0,48),ylim=c(-21,-11))
par(new=TRUE)
for(i in 2:97){
	plot(ttime,arr2[,i],col=sample(color,1),xlab="",ylab="",type="l",lwd=1,axes=FALSE,cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	par(new=TRUE)
}
##########################################################################

#Figure 4
##E(sigma^2)
load("multGP1")#.50
load("multGP2")#.80
load("multGP3")#.90
load("multGP4")#.95
#par(mar = c(4, 5, 3, 1) + 0.1, cex =1.2)
f1=apply(multGP1[[1]][tt,,,-c(1:burn)],c(2,3),rowMeans)
f2=apply(multGP2[[1]][tt,,,-c(1:burn)],c(2,3),rowMeans)
f3=apply(multGP3[[1]][tt,,,-c(1:burn)],c(2,3),rowMeans)
f4=apply(multGP4[[1]][tt,,,-c(1:burn)],c(2,3),rowMeans)
r1=1
ttime=0.5*(1:96)
plot(ttime,f1[,r1,r1],lwd=2,"l",xlab="Time (hr)",ylab="Posterior variance",main=expression(E(sigma[11]^2)),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(ttime,f2[,r1,r1],col="green",lwd=2,ylim=c(0,1.7))
lines(ttime,f3[,r1,r1],col="blue",lwd=2,ylim=c(0,1.7))
lines(ttime,f4[,r1,r1],col="red",lwd=2,ylim=c(0,1.7))
legend("topleft",c("0.50","0.80","0.90","0.95"),title="Discount factor",
fill=c("blue","green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

##
r1=2
ttime=0.5*(1:96)
plot(ttime,f1[,r1,r1],lwd=2,"l",xlab="Time (hr)",ylab="Posterior variance",main=expression(E(sigma[22]^2)),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(ttime,f2[,r1,r1],col="green",lwd=2,ylim=c(0,2.7))
lines(ttime,f3[,r1,r1],col="blue",lwd=2,ylim=c(0,2.7))
lines(ttime,f4[,r1,r1],col="red",lwd=2,ylim=c(0,2.7))
legend("topleft",c("0.50","0.80","0.90","0.95"),title="Discount factor",
fill=c("blue","green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)
##
r1=3
ttime=0.5*(1:96)
plot(ttime,f1[,r1,r1],lwd=2,"l",xlab="Time (hr)",ylab="Posterior variance",main=expression(E(sigma[33]^2)),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(ttime,f2[,r1,r1],col="green",lwd=2,ylim=c(0,11))
lines(ttime,f3[,r1,r1],col="blue",lwd=2,ylim=c(0,11))
lines(ttime,f4[,r1,r1],col="red",lwd=2,ylim=c(0,11))
legend("topleft",c("0.50","0.80","0.90","0.95"),title="Discount factor",
fill=c("blue","green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

#
r1=4
ttime=0.5*(1:96)
plot(ttime,f1[,r1,r1],lwd=2,"l",xlab="Time (hr)",ylab="Posterior variance",main=expression(E(sigma[44]^2)),cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(ttime,f2[,r1,r1],col="green",lwd=2,ylim=c(0,21))
lines(ttime,f3[,r1,r1],col="blue",lwd=2,ylim=c(0,21))
lines(ttime,f4[,r1,r1],col="red",lwd=2,ylim=c(0,21))
legend("topright",c("0.50","0.80","0.90","0.95"),title="Discount factor",
fill=c("blue","green","black","red"),bty="0",border="black",cex=1.5,horiz=FALSE)

