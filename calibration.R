##########11-11-2018/ PAPER RESULTS

library(MASS)
library(emulator)
library(invgamma)
library(compiler)
library(abind)
library(LaplacesDemon)
#library(MCMCpack)
source("rwish.R",echo=FALSE)
nam=c("sub","o2","nh4","muhet","muaob","Yhet","Yaob")
N=97
ntime=72
nout=4
############################################read in data
gap=gap1=gap2=list()
for(j in 1:N){
gap[[j]]=as.matrix(read.csv(sub("j",j,paste("data","input_j.csv",sep="/")),header=TRUE))
gap1[[j]]=as.matrix(read.csv(sub("j",j,paste("data","output_j.csv",sep="/")),header=TRUE))
gap2[[j]]=as.matrix(read.csv(sub("j",j,paste("data2","idynotrain_j.csv",sep="/")),header=TRUE))
}
input=abind(gap,along=3)
input1=input[1:ntime,,]
design0=t(input[1,,])
y=abind(gap1,along=3)
#####################scale
min1=apply(y,2,min)
max1=apply(y,2,max)
Y=y
for(j in 1:nout){
Y[,j,]=(y[,j,]-min1[j])/(max1[j]-min1[j])
}
y=Y
y=y[1:ntime,,]
arx=y[1:ntime,,]
arx2=abind(gap2,along=3)
#####################scale idynomics
min11=apply(arx2,2,min)
max11=apply(arx2,2,max)
Y=arx2
for(j in 1:nout){
Y[,j,]=(arx2[,j,]-min11[j])/(max11[j]-min11[j])
}
arx2=Y
load("results/multGP3")
GP=multGP3
#load("data/min2")
#load("data/max2")
max2min2=as.matrix(read.csv("max2min2.txt",sep=" ",header=FALSE))
max2=c(max2min2[1,])
min2=c(max2min2[2,])
r=dim(y)[3]
#################

ndes=1
ninp=7
p=2
rr=97#75#1
q=4
sig=.05#1.38#.8#.001#.01
sig2=1.38/q#.001#.1
lambda=.001
dess=c(30e-3,10e-3,6e-3,0.00006944444,0.00001158333,0.63,0.24)##
etai=(dess-min2)/(max2-min2)##etai_0
eta0=etai[1:3]
etai=etai[4:7]
a=0.0*etai;b=2.0*etai##uniform distr parameter for etai
#inverse gamma distri parameter for sigmaz
sigmazi=sig2*diag(1,q*r)#rep(.01,ndes)


MH=function(betai,design0,Erri,Vi,p,Thetai,obsinput,arx2,etai,sigmazi){
T=ntime
St=Mt=list()
Sigma=corr.matrix(design0[1:rr,,drop=FALSE],design0[1:rr,,drop=FALSE],scales=betai,method=1)+diag(lambda,rr)
SigmaInv=solve(Sigma)
gg=design0[1:rr,,drop=FALSE]
gg[,4:7]=matrix(etai,nrow=rr,ncol=length(etai),byrow=TRUE)
rho=corr.matrix(gg,design0[1:rr,,drop=FALSE],scales=betai,method=1)
A1=0;A2=0
llk=0
for(t in p:T){
inp=input1[1,,]
inp[4:7,]=t(matrix(etai,nrow=rr,ncol=length(etai),byrow=TRUE))
Mt[[t]] = t(Thetai[1,,])%*%inp+ t(t(rho)%*%SigmaInv%*%t(Erri[t,,1:rr]))
#pp=diag(rep(diag(sigmazi),each=rr))+Vi[t,,]%x%(1-t(rho)%*% SigmaInv %*% rho)
pp=diag(diag(sigmazi+Vi[t,,]%x%(1-t(rho)%*% SigmaInv %*% rho)))
St[[t]]=(pp+t(pp))/2
A1=A1+(determinant(((St[[t]])))$modulus[1])
mut=c(Mt[[t]])#matrix(Mt[[t]],nrow=q,ncol=rr)
A2=A2+c(arx2[t,,1:rr]-mut)%*%solve(St[[t]])%*%c(arx2[t,,1:rr]-mut)##sum((arx2[2:96,r]-unlist(Mt))^2/unlist(St[2:96]))
}
LogPrior2=0
LogPrior3=0
for(i in 1:length(etai)){
LogPrior2= LogPrior2+log(dunif(etai[i],min=0,max=1))#
}
for(i in 1:(q*r)){
LogPrior3=LogPrior3+log(dinvgamma(sigmazi[i,i],shape=3,scale=.5))
}
score=c(-0.5*A1-0.5*A2+LogPrior2+LogPrior3)#+sum(log(etai))+sum(log(tr(sigmazi))))
return(score)
}

MH2=function(betai,design0,Erri,Vi,p,Thetai,obsinput,arx2,etai,sigmazi,scorei){
T=nrow(Erri)
##Proposal function for calibration eta and sigmaz
ll=length(etai)
ll2=length(sigmazi)
eta=exp(log(etai)+ sig*rnorm(ll,0,1))
sigmaz=exp(log(sigmazi)+ .05*rnorm(q*r,0,1))
score=MH(betai,design0,Erri,Vi,p,Thetai,obsinput,arx2,etai=eta,sigmazi=sigmaz)
acceptance=score-scorei
if(log(runif(1,0,1))< acceptance){
eta=eta;sigmaz=sigmaz;score=score
}else{
eta=etai;sigmaz=sigmazi;score=scorei
}
list(eta=eta,sigmaz=sigmaz,score=score)
}
MH=cmpfun(MH)
MH2=cmpfun(MH2) 

ww=5000
thin=5 
burn=2000
ETA=array(0,c(ninp-3,ww/thin))
SIGMAZ=array(0,c(q*r,q*r,ww/thin))
Vi=apply(GP[[1]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x2x2	
Thetai=apply(GP[[2]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x7x2
Erri=apply(GP[[3]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x2x75
betai=apply(GP[[4]][,-(1:burn)],1,mean) 
 
for(w in 1:ww){
s=sample(1:2000,1)
Vi=GP[[1]][1:ntime,,,s]	
Thetai=GP[[2]][1:ntime,,,s]
Erri=GP[[3]][1:ntime,,,s]
betai=GP[[4]][,s]
h=w
ind1=w%%thin
n=h[!ind1]
scorei=MH(betai,design0,Erri,Vi,p,Thetai,obsinput,arx2,etai,sigmazi)
hh=MH2(betai,design0,Erri,Vi,p,Thetai,obsinput,arx2,etai,sigmazi,scorei)
etai=hh$eta;scorei=hh$score;sigmazi=hh$sigmaz#
print(w)
ETA[,n/thin]=hh$eta
SIGMAZ[,,n/thin]=hh$sigmaz 
}
###############
length(unique(ETA[1,]))/ncol(ETA)

path=getwd()
save(ETA,file=paste(path,"ETA5",sep="/"))
save(SIGMAZ,file=paste(path,"SIGMAZ5",sep="/"))

#############################################################################PAPER PLOTS
load("ETA5")
load("SIGMAZ5")
length(unique(ETA[1,]))/ncol(ETA)
length(unique(SIGMAZ[1,1,]))/length(SIGMAZ[1,1,])

burn=500
V0=apply(GP[[1]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x2x2	
THETA0=apply(GP[[2]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x7x2
ERR0=apply(GP[[3]][1:ntime,,,-(1:burn)],c(2,3),rowMeans)#96x2x75
BETA0=apply(GP[[4]][,-(1:burn)],1,mean)

T=72
burn=4000#2000
est=apply(ETA[,-(1:burn)],1,mean)# 
est2=matrix(est,nrow=97,ncol=length(est),byrow=TRUE)
vv=apply(SIGMAZ[,,-(1:burn)],2,rowMeans)
v0=array(diag(vv),c(97,4))
emu0=list()
dev=devv=list()
for(t in 1:T){
eta00=t(input1[t,1:3,])
aa=cbind(eta00,est2)
Sigma=corr.matrix(design0[,,drop=FALSE],design0[,,drop=FALSE],scales=BETA0)+diag(lambda,rr)
SigmaInv=solve(Sigma)
rho=corr.matrix((aa),design0[,,drop=FALSE],scales=BETA0)
emu0[[t]]=t(THETA0[t,,])%*%t(aa) + t(t(rho)%*%SigmaInv%*%t(ERR0[t,,]))
#dev[[t]]=cbind(v0[1,1]+V0[t,1,1]*diag(diag(1,r) - t(rho)%*% SigmaInv %*% rho),vv[2,2]+V0[t,2,2]*diag(diag(1,r) - t(rho)%*% SigmaInv %*% rho))

#pp=diag(rep(diag(vv),each=rr))+V0[t,,]%x%(1-t(rho)%*% SigmaInv %*% rho)
pp=vv+V0[t,,]%x%(1-t(rho)%*% SigmaInv %*% rho)
devv[[t]]=matrix(diag(pp),nrow=97,ncol=4)
}
emu=abind(emu0,along=0)#t(emu0[[t]])
cbind(emu[,,56],arx2[,,56])
diag(cor(emu[,1,],arx2[,1,])^2)
##################predictions
tit=c("Biomass concentration (kg/m^3)","Total no of particles","Biofilm height (m)","Surface roughness (m)")
emuu=emu
emuvv=abind(devv,along=0)
ttime=(1:72)
TT=72
library(plotrix)
lap=c(15,10,72,80)
#k=15#72
par(mar = c(4, 5, 3, 1) + 0.1, cex = 1.5)
par(mfrow=c(2,2))
lim=list(c(-11,45),c(-5,45),c(-17,45),c(-23,45))
for(j in 1:nout){
sd=2*sqrt(emuvv[1:TT,k,j])*(max1[j]-min1[j])
y1=(arx2[1:TT,j,k])*(max11[j]-min11[j])+min11[j]#idyno
y2=(arx[1:TT,j,k])*(max1[j]-min1[j])+min1[j]##NUFEB
y3=(emuu[1:TT,j,k])*(max1[j]-min1[j])+min1[j]##emulator
dev1=(emuu[1:TT,j,k]-sd)
dev2=(emuu[1:TT,j,k]+sd)
coll=c("red","green","blue","black")
F=y3
U=dev2
L=dev1
#expression(paste(sub("j",j,tit[j])," ", (kgm^-3)))
plotCI(ttime[1:TT],F,ui=U,li=L,ylab=(paste(sub("j",j,tit[j])))
,xlab="Time (hr)",main="",cex=1.2,cex.lab=1.2,cex.axis=1.2,
cex.main=1.8,col="green",scol="red",slty=1,lwd=1,xlim=c(0,80),ylim=lim[[j]],cex.lab=1.4)
points(ttime,y2,col="black",ylab="",xlab="","o")
points(ttime,y1,col="blue",ylab="",xlab="","o")
legend("topleft",c("iDynoMiCS","NUFEB","Emulator","95% C.I"),
fill=c("blue","black","green","red"),bty="0",border="black",cex=1.5,horiz=FALSE)
}


###################Figure 6
tvec=c(1:96)
burn=4000
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(2,2))
for(j in 1:4){
gap=(ETA[j,-(1:burn)])
x=seq(0,1,length=1000)
#hist(gap,breaks=15,col="green",xlab="",prob=T,main=substitute(paste('Density of ',eta[a]),list(a=tvec[j])),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,.4))
hist(gap,breaks=15,col="green",xlab="",prob=T,main=substitute(paste('Density of ',eta[a]),list(a=tvec[j])),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
par(new=TRUE)
curve(dunif(x,0,1),add=TRUE,from=0,to=1,axes=FALSE, ylab="",xlab="",col="red",lwd=2)
}
#################or Figure 6
burn=4000
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(2,2))
for(j in 1:4){
gap=(ETA[j,-(1:burn)])
x=seq(0,1,length=1000)
foo=hist(gap,breaks=8,xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ',eta[a]),list(a=tvec[j])),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,1))
aa=c(0,foo$mids,1)
lab=c(0,signif(seq(min(gap),max(gap),l=length(foo$mids)),digits=1),1)
axis(side=1,xpd=NA,at=aa,labels=lab,cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,1))
par(new=TRUE)
curve(dunif(x,0,1),add=TRUE,from=0,to=1, ylab="",xlab="",col="red",lwd=2)
}
#####################Figure 7
#
k=0
tvec=c(1,2,3,4)
lap=c(k+1,k+98,k+195,k+292)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(2,2))
for(r in 1:4){
j=(97*r)-94
hh=SIGMAZ[j,j,-(1:burn)]
foo=hist(hh,breaks=8,panel.first=grid(),xaxt="n",col="green",xlab="",prob=T,main=substitute(paste('Density of ',sigma[a]^2),list(a=tvec[r])))
x=seq(0,1,length=1000)
axis(side=1,at=foo$mids,labels=signif(seq(min(hh),max(hh),l=length(foo$mids)),digits=3),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
par(new=TRUE)
curve(dinvgamma(x,shape=3,scale=.5),ylab="",xlab="",col="red",lwd=2,add=TRUE)
}
#############
tvec=c(1:96)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 0.6)
par(mfrow=c(4,3))
for(j in 1:4){
hist(ETA[j,-(1:burn)],col="green",xlab="",prob=T,main=substitute(paste('Density of ',eta[a]),list(a=tvec[j])))
lines(density(ETA[j,-(1:burn)],kernel="rectangular"),col="red",xlab=sub("j",j,paste("Eta","j",sep="_")))
ts.plot(ergMean(ETA[j,-c(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean")
ts.plot(ETA[j,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots")
}

par(mfrow=c(4,3))
for(j in 1:4){
hist(SIGMAZ[j,j,-(1:burn)],col="green",xlab="",prob=T,main=substitute(paste('Density of ',sigma[a]^2),list(a="z")))
lines(density(SIGMAZ[j,j,-(1:burn)]),col="red",xlab=expression(sigma^2),cex.main=1.2,cex.axis=1.2)
plot(ergMean(SIGMAZ[j,j,-c(1:burn)]),col="red",xlab="Iterations",ylab="ErgMean",cex.axis=1.3)
plot(SIGMAZ[j,j,-(1:burn)],col="blue",xlab="Iterations",ylab="Traceplots",cex.axis=1.3)
}
#########################Table 4
olu=apply(ETA[,-(1:burn)],1,mean)
u1= quantile(ETA[1,-(1:burn)],probs=c(.025,.975))
u2= quantile(ETA[2,-(1:burn)],probs=c(.025,.975))
u3= quantile(ETA[3,-(1:burn)],probs=c(.025,.975))
u4= quantile(ETA[4,-(1:burn)],probs=c(.025,.975))
uu=rbind(u1,u2,u3,u4)
(uu*(max2[4:7]-min2[4:7]))+min2[4:7]
(olu*(max2[4:7]-min2[4:7]))+min2[4:7]

