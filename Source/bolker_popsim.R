########################
# CRL - lmer model fit #
########################

### Steady growth

g <- 1.2
nt <- 20

year <- 1:nt
N <- numeric(nt)
Nobs <- numeric(nt)
N[1] <- 10000
temp <- rnorm(nt,1,0.1)
for (t in 1:(nt-1)) {
	N[t+1] <- rpois(1,lambda=g*N[t]*temp[t])
	}

par(mfrow=c(2,1),las=1)
plot(N,type="b")
plot(N[-1]/N[-length(N)],type="b",ylab=expression(lambda))

### Fit model


require(lme4)
mod <- lmer(N[-1]~(1|year[-length(N)])+offset(log(N[-length(N)])),family="poisson")
summary(mod)
exp(coef(mod))

### Hyperbolic

#b <- 0
#a <- 1
#N[1] <- 20
#for (t in 1:(nt-1)) {
#	N[t+1] <- rpois(1,lambda=a*N[t]/(b+N[t]))
#	}



###################################################
### chunk number 1: 
###################################################
library(emdbook) ## for trcoef
require(bbmle)
require(Hmisc)
require(MASS) ## for eqscplot
source("chapskel.R")

###################################################
### chunk number 2: 
###################################################
## linear simulation with process/observation error
linsim = function(nt=20,N0=2,dN=1,sd_process=sqrt(2),
  sd_obs=sqrt(2)) {
  Nobs = numeric(nt)
  N_cur=N0
  Nobs[1]=N_cur+rnorm(1,sd=sd_obs)
  for (i in 2:nt) {
    N_cur=N_cur+dN+rnorm(1,sd=sd_process)
    Nobs[i]=N_cur+rnorm(1,sd=sd_obs)
  }
  Nobs
}
N = linsim()
nsim=1000
nt=20
Nmat_obsonly = matrix(ncol=nsim,nrow=nt)
for (j in 1:nsim) {
  Nmat_obsonly[,j] = linsim(sd_process=0,sd_obs=2)
}
env_obs = t(apply(Nmat_obsonly,1,quantile,c(0.025,0.975)))
Nmat_proconly = matrix(ncol=nsim,nrow=nt)
for (j in 1:nsim) {
  Nmat_proconly[,j] = linsim(sd_process=2,sd_obs=0)
}
env_proc = t(apply(Nmat_proconly,1,quantile,c(0.025,0.975)))

###################################################
### chunk number 3: 
###################################################
## discrete pop with process/observation error
  ## equilibrium of N(t+1)=N*a/(b+N): a/(b+N)=1 or N=a-b
  ## with detection probability p we get (a-b)*p
  ## pure process: Poisson variation; p=1, mean=(a-b), var=(a-b)
  ## pure observation: Binomial variation around (a-b)*p;
  ## mean = (a-b)*p, var= (a-b)*p*(1-p)
  ## solve: M,V equal with (a1-b1) = (a2-b2)*p
  ##                       (a1-b1) = (a2-b2)*p*(1-p)
  ##                       poismean = (a2-b2)*p
  ##                       can't do it -- not possible
  ## have to settle for equal means

dsim = function(nt=20,N0=(a-b),a=6,b=1,p=1,
  proc_err=TRUE) {
  Nobs=numeric(nt)
  N_cur=N0
  Nobs[1]=rbinom(1,size=N_cur,prob=p)
  for (i in 2:nt) {
    if (proc_err) {
      N_cur = rpois(1,N_cur*a/(b+N_cur))
      Nobs[i] = rbinom(1,size=N_cur,prob=p)
    } else {
      N_cur = N_cur*a/(b+N_cur)
      Nobs[i] = rbinom(1,size=floor(N_cur)+rbinom(1,size=1,prob=N_cur-floor(N_cur)),
            prob=p)
    }
  }
  Nobs
}
nt=20
dN=dsim()
nsim=1000
nt=20
dNmat_obsonly = matrix(ncol=nsim,nrow=nt)
for (j in 1:nsim) {
  dNmat_obsonly[,j] = dsim(proc_err=FALSE,p=5/8,a=9)
}
denv_obs = t(apply(dNmat_obsonly,1,quantile,c(0.025,0.975)))
dNmat_proconly = matrix(ncol=nsim,nrow=nt)
for (j in 1:nsim) {
  dNmat_proconly[,j] = dsim(p=1)
}
denv_proc = t(apply(dNmat_proconly,1,quantile,c(0.025,0.975)))


###################################################
### chunk number 4:  eval=FALSE
###################################################
##   long_dNmat_proconly = matrix(ncol=5000,nrow=100)
##   for (j in 1:5000) {
##     long_dNmat_proconly[,j] = dsim(nt=100,p=1)
##   }
## long_denv_proc = t(apply(long_dNmat_proconly,1,quantile,c(0.025,0.975)))
## matplot(1:100,long_denv_proc,type="l",col=1,lty=1)


###################################################
### chunk number 5: 
###################################################
op = par(mfrow=c(1,2))
par(cex=1.5,mar=c(4,4,0,1)+0.1,
  mgp=c(2.5,0.75,0),las=1,bty="l")
plot(1:nt,N,ylim=c(-8,35),type="b",xlab="Time",lwd=2,
     ylab="Population density")
matlines(1:nt,env_obs,type="l",lty=2,col=1)
matlines(1:nt,env_proc,type="l",lty=3,col=1)
text(13,30,"Process",adj=0)
text(15,22.5,"Observation",adj=0)
##
plot(1:nt,dN,type="b",ylim=c(0,max(denv_proc)+1),lwd=2,xlab="Time",
     ylab="Population numbers")
matlines(1:nt,denv_obs,type="l",lty=2,col=1)
matlines(1:nt,denv_proc,type="l",lty=3,col=1)
text(3,10.5,"Process",adj=0)
text(5.5,7.5,"Observation",adj=0)
par(op)


###################################################
### chunk number 6:  eval=FALSE
###################################################
## N = numeric(nt)
## Nobs = numeric(nt)
## N[1] = a
## for (t in 1:(nt-1)) {
##   N[t+1] = b+N[t]
##   Nobs[t]=rnorm(1,mean=N[t],sd=sd.obs)
## }
## Nobs[nt] = rnorm(1,mean=N[nt])


###################################################
### chunk number 7:  eval=FALSE
###################################################
## N = numeric(nt)
## Nobs = numeric(nt)
## N[1] = a
## for (t in 1:(nt-1)) {
##   N[t+1] = rnorm(1,mean=b+N[t],sd=sd.proc)
##   Nobs[t]=N[t]
## }
## Nobs[nt]=N[nt]


###################################################
### chunk number 8:  eval=FALSE
###################################################
## N = numeric(nt)
## Nobs = numeric(nt)
## N[1] = N0
## for (t in 1:(nt-1)) {
##   N[t+1] = a*N[t]/(b+N[t])
##   Nobs[t]=rbinom(1,size=round(N[t+1]),prop=p)
## }
## Nobs[nt]=rbinom(1,size=round(N[nt]),prop=p)


###################################################
### chunk number 9:  eval=FALSE
###################################################
## N = numeric(nt)
## Nobs = numeric(nt)
## N[1] = N0
## for (t in 1:(nt-1)) {
##   N[t+1] = rpois(1,lambda=a*N[t]/(b+N[t]))
##   Nobs[t]=N[t]
## }
## Nobs[nt]=N[nt]


###################################################
### chunk number 10: 
###################################################
tvec = 1:200
a.true=5
b.true=0.05 ## 0.01
x0 = rnorm(200,mean=a.true+b.true*tvec,sd=2)
x = x0[-200]
y = x0[-1]
lm.ols = lm(x0~tvec)
lm1 = lm(y~x)
lm2 = lm(x~y)
tmpf=function(p) {
  a=p[1]
  b=p[2]
  sum((y-b*x-a)^2/(b^2+1))
}
O1 = optim(fn=tmpf,par=c(1,1))
a1 = arima(x0,c(1,0,0))


###################################################
### chunk number 11: 
###################################################
op=par(pty="s",mfrow=c(1,2),cex=1.5,mgp=c(2.5,1,0),
  mar=c(4,4,2,1)+0.1,las=1,lwd=2,bty="l")
plot(x0,xlab="Time",ylab="N(t)",type="l")
abline(lm.ols,lwd=2)
plot(x,y,
     xlab="N(t)",ylab="N(t+1)",col="gray")
xvec = seq(floor(min(x)),ceiling(max(x)),length=100)
matlines(xvec,predict(lm1,interval="prediction",newdata=data.frame(x=xvec)),
         col=1,lty=c(1,2,2))
invisible(require(ellipse,quietly=TRUE))
cov1 = cov(cbind(x,y))
lines(ellipse(cov1,centre=c(mean(x),mean(y))),lty=2)
## calculate principal axis
e1 = eigen(cov1)$values[1]
rmaslope = sqrt(coef(lm1)[2]*coef(lm2)[2])
## y = a+e1*x
##abline(a=mean(y)-e1*mean(x),b=e1)
##abline(a=mean(y)-rmaslope*mean(x),b=rmaslope,col=2)
abline(a=O1$par[1],b=O1$par[2],lty=2)
par(xpd=NA)
legend(2,25,c("process error","observation error"),
       lty=1:2,bty="n")
par(xpd=FALSE)
par(op)


###################################################
### chunk number 12: 
###################################################
## simulate logistic y with process and observation error
set.seed(1001)
r = 1
K = 10
t0 = 5
n0 = 1
tot.t = 10
dt=0.5
sd.proc= 1
sd.obs = 1
set.seed(1001)
tvec = seq(1,tot.t,by=dt)
n=length(tvec)
y = numeric(n)
ydet=numeric(n)
y[1] = n0
ydet[1] = n0
e.proc = rnorm(n,mean=0,sd=sd.proc)
e.obs = rnorm(n,mean=0,sd=sd.obs)
for (i in 2:n) {
  ydet[i] = ydet[i-1]+r*ydet[i-1]*(1-ydet[i-1]/K)*dt
  y[i] = y[i-1]+(r*y[i-1]*(1-y[i-1]/K)+e.proc[i-1])*dt ## process only
}
## sd is variance in GROWTH RATE: should translate to
## sd.proc/4 with delta-t=0.5
y.procobs = y+e.obs
y.obs = ydet+e.obs
X = cbind(ydet,y,y.procobs,y.obs)


###################################################
### chunk number 13: 
###################################################
t0 = 1
## fit to logistic by shooting
shootfun = function(n0,r,K,sigma) {
  y.pred = K/(1+(K/n0-1)*exp(-r*(tvec-t0)))
  -sum(dnorm(y.procobs,y.pred,sd=sigma,log=TRUE))
}
m.shoot = mle2(shootfun,start=list(n0=1,r=1,K=10,sigma=1),
  method="Nelder-Mead")


###################################################
### chunk number 14: 
###################################################
## calculate diagonal points???
## find the intersection of (xn,yn)-a*(x+y) and (x,K/(1+(K/n0-1)*exp(r*(x-t0))))
##  xn + a*D = x
##  yn - a*D = K/(1+(K/n0-1)*exp(r*(x-t0)))
## solve for D:
##  yn - a*D = K/(1+(K/n0-1)*exp(r*(xn+a*D-t0)))
## transcendental, I'm afraid
intdist = function(x,y,pars) {
  tmpf = function(D) { with(as.list(pars),y-D-K/(1+(K/(x+D)-1)*exp(-r*dt)))}
  D = uniroot(tmpf,c(-10,10))$root
}
D = numeric(n-1)
for (i in 1:(n-1)) {
  D[i] = intdist(y.procobs[i],y.procobs[i+1],coef(m.shoot))
}   


###################################################
### chunk number 15: 
###################################################
op = par(mfrow=c(1,2))
par(cex=1.5,mar=c(4,4,0,1)+0.1,
  mgp=c(2.5,0.75,0),las=1,bty="l",
  mfrow=c(1,2))
## N vs t
plot(tvec,y.procobs,xlab="Time",ylab="N(t)",ylim=c(0,15))
with(as.list(coef(m.shoot)),curve(K/(1+(K/n0-1)*exp(-r*(x-t0))),add=TRUE))
y.pred = with(as.list(coef(m.shoot)),K/(1+(K/n0-1)*exp(-r*(tvec-t0))))
segments(tvec,y.procobs,tvec,y.pred,lty=2)
points(tvec,y.procobs,pch=16,col="gray")
##text(6,4,paste(names(coef(m.shoot)),round(coef(m.shoot),2),sep="=",collapse="\n"))
## N(t+1) vs N(t)
## FIXME: not sure this is really an improvement
eqscplot(y.procobs[-n],y.procobs[-1],xlab="N(t)",ylab="N(t+1)",
         xlim=c(0,15),ylim=c(0,15),tol=0)
with(as.list(coef(m.shoot)),curve(K/(1+(K/x-1)*exp(-r*dt)),add=TRUE))
segments(y.procobs[-n],y.procobs[-1],y.procobs[-n]+D,y.procobs[-1]-D,lty=2)
points(y.procobs[-n],y.procobs[-1],pch=16,col="gray")
par(op)


###################################################
### chunk number 16: 
###################################################
t0 = 1
## fit to logistic by one-step-ahead
stepfun = function(r,K,sigma) {
  y2 = y.procobs[-n]
  y.pred = K/(1+(K/y2-1)*exp(-r*dt))
  -sum(dnorm(y.procobs[-1],y.pred,sd=sigma,log=TRUE))
}
m.step = mle2(stepfun,start=list(r=1,K=10,sigma=1),method="Nelder-Mead",
  control=list(parscale=c(1,10,1)))


###################################################
### chunk number 17: 
###################################################
op = par(cex=1.5,mar=c(4,4,0,1)+0.1,
  mgp=c(2.5,0.75,0),las=1,bty="l",
  mfrow=c(1,2))
plot(tvec,y.procobs,pch=16,ylim=c(0,15),
xlab="time",ylab="N(t)")
logist = function(x0,t,r=1,K=10) {
  K/(1+(K/x0-1)*exp(-r*dt))
}
y.pred = with(as.list(coef(m.step)),K/(1+(K/y.procobs[-n]-1)*exp(-r*dt)))
arrows(tvec[-n],y.procobs[-n],tvec[-1],y.pred,length=0.1,angle=20)
points(tvec[-1],y.pred)
segments(tvec[-1],y.pred,tvec[-1],y.procobs[-1],lty=2)
##text(6,4,paste(names(coef(m.step)),round(coef(m.step),2),sep="=",collapse="\n"))
legend("topleft",c("observed","predicted"),
       pch=c(16,1))
##
plot(y.procobs[-n],y.procobs[-1],xlab="N(t)",ylab="N(t+1)",xlim=c(0,15),ylim=c(0,15))
with(as.list(coef(m.step)),curve(K/(1+(K/x-1)*exp(-r*dt)),add=TRUE))
segments(y.procobs[-n],y.pred,y.procobs[-n],y.procobs[-1],lty=2)
points(y.procobs[-n],y.procobs[-1],pch=16,col="gray")
par(op)


###################################################
### chunk number 18: 
###################################################
simlogistdata = function(seed=1001,
  r=1,K=10,n0=1,t0=1,tot.t=10,dt=0.5,
  sd.proc=1,sd.obs=1) {
  if (!is.null(seed)) set.seed(seed)
  tvec = seq(1,tot.t,by=dt)
  n=length(tvec)
  y = numeric(n)
  ydet=numeric(n)
  y[1] = n0
  ydet[1] = n0
  e.proc = rnorm(n,mean=0,sd=sd.proc)
  e.obs = rnorm(n,mean=0,sd=sd.obs)
  for (i in 2:n) {
    ydet[i] = ydet[i-1]+r*ydet[i-1]*(1-ydet[i-1]/K)*dt
    y[i] = y[i-1]+(r*y[i-1]*(1-y[i-1]/K)+e.proc[i-1])*dt ## process only
  }
  y.procobs = y+e.obs
  y.obs = ydet+e.obs
  cbind(tvec,y.procobs)
}


###################################################
### chunk number 19:  eval=FALSE
###################################################
## ylast = y[-1]
## yfirst = y[-n]
## sum((ylast-yfirst)*yfirst*(1-yfirst/10))/sum(yfirst*(1-yfirst/10))
## -sum(yfirst^2)/sum(ylast-yfirst*(1+1))


###################################################
### chunk number 20: 
###################################################
## fit to logistic by one-step-ahead (assume proc. error only)
stepfun2 = function(r,K,sigma,y) {
  ystart = y[-n]
  yend = y[-1]
  y.pred = K/(1+(K/ystart-1)*exp(-r*dt))
  -sum(dnorm(yend,y.pred,sd=sigma,log=TRUE))
}
jagged = FALSE
newdata = simlogistdata(tot.t=100,sd.proc=2)
n = nrow(newdata)
y.procobs = newdata[,2]
tvec = newdata[,1]
randvals = rnorm(n,mean=0,sd=1)
simex0 = function(s,...) {
  y.simex = y.procobs+if (jagged) rnorm(n,mean=0,sd=s) else randvals*s
  coef(mle2(stepfun2,start=list(r=1,K=10,sigma=1),
                      data=list(y=y.simex),...))
}     
sdvec = seq(0,2,length=10)
simextab = t(sapply(sdvec,simex0,method="Nelder-Mead"))
predmat = apply(simextab,1,
  function(X) {
    r=X[1]; K=X[2]; n0=y.procobs[1]
    K/(1+(K/n0-1)*exp(-r*(tvec-t0)))
  })
##matplot(predmat,type="b")
rvec = simextab[,"r"]
Kvec = simextab[,"K"]
tsdvec = sqrt(sd.obs^2+sdvec^2)
tsdvec2a = seq(0,1,length=40)
tsdvec2b = seq(1,max(tsdvec),length=40)
q.r = lm(rvec~tsdvec+I(tsdvec^2))
q.K = lm(Kvec~tsdvec+I(tsdvec^2))
l.r = lm(rvec~tsdvec)
l.K = lm(Kvec~tsdvec)


###################################################
### chunk number 21: 
###################################################
op=par(mar=c(5,5,2,5)+0.1,las=1,cex=1.5,lwd=2,bty="l")
plot(tsdvec,rvec,xlim=c(0,max(tsdvec)),
     xlab="Total observation error",ylab="",ylim=c(0.9,max(rvec)))
mtext(side=2,"r",line=3,cex=1.5)
lines(tsdvec2a,predict(q.r,newdata=data.frame(tsdvec=tsdvec2a)),lty=2)
lines(tsdvec2b,predict(q.r,newdata=data.frame(tsdvec=tsdvec2b)))
abline(h=1,lty=3)
points(par("usr")[1],coef(q.r)[1],pch=16,xpd=NA)
## K plot
par(new=TRUE)
plot(tsdvec,Kvec,pch=2,xlim=c(0,max(tsdvec)),ylim=c(9,10.5),axes=FALSE,
     xlab="",ylab="",col="darkgray")
lines(tsdvec2a,predict(q.K,newdata=data.frame(tsdvec=tsdvec2a)),lty=2,col="darkgray")
lines(tsdvec2b,predict(q.K,newdata=data.frame(tsdvec=tsdvec2b)),col="darkgray")
axis(side=4,col="darkgray",col.axis="darkgray")
mtext(side=4,at=9.75,line=par("mgp")[1],"K",col="darkgray")
abline(h=10,lty=3,col="darkgray")
points(par("usr")[2],coef(q.K)[1],pch=16,col="darkgray",xpd=NA)
par(op)


###################################################
### chunk number 22:  eval=FALSE
###################################################
## #library(MASS)
## dx = 3
## dy = 3
## sc = 0.5
## xpos = rep(0:1,each=3)*dx
## ypos = rep(0:2,2)*dy
## symbols(xpos,ypos,circles=rep(1,6),inches=sc,
##         xlab="",ylab="",xlim=c(-1,5),ylim=c(-1,7))
## for (i in 1:3) {
##   text(0,ypos[i],
##        substitute(N[i],list(i=i)))
## }
## for (i in 1:3) {
##   text(dx,ypos[i],
##        substitute(N[list(obs,i)],list(i=i)))
## }


###################################################
### chunk number 23: 
###################################################
kfpred.gen <- function(A,B,C,D,U,V,N0,MX.start,WX.start,Y) {
  ## follows Schnute 1994 notation
  ## would need to allocate space for answers
  MX <- MX.start
  WX <- WX.start
  MY[,1] <- C + D %*% MX.start
  WY[,,1] <- D %*% WX.start %*% t(D) + V
  for (t in 2:nt) {
    MXi <- A+B %*% MX
    WXi <- B %*% WX %*% t(B) + U
    MY[,t] <- C+D %*% MXi
    WY[,,t] <- D %*% WXi %*% t(D) + V
    MX <- MX + WX %*% t(D) %*% solve(WY[,,t]) %*% (Y[,t]-MY[,t])
    WX <- WX - WX %*% t(D) %*% solve(WY[,,t]) %*% D %*% WX
  } 
}


###################################################
### chunk number 24:  eval=FALSE
###################################################
## kfpred <- function(a,b,procvar,obsvar,M.n.start,Var.n.start,Nobs,c=1,d=1) {
##   nt <- length(Nobs)
##   M.nobs <- numeric(nt)
##   Var.nobs <- numeric(nt)
##   M.n <- M.n.start
##   Var.n <- Var.n.start
##   M.nobs[1] <- M.n.start
##   Var.nobs[1] <- Var.n.start+obsvar
##   for (t in 2:nt) {
##     M.ni <- a+b*M.n
##     Var.ni <- b^2*Var.n + procvar
##     M.nobs[t] <- M.ni
##     Var.nobs[t] <- Var.ni + obsvar
##     M.n <- M.ni +  Var.ni/Var.nobs[t]*(Nobs[t]-M.nobs[t])
##     Var.n <- Var.ni*(1 - Var.ni/Var.nobs[t])
##   } 
##   list(mean=M.nobs,var=Var.nobs)
## }
## kflik <- function(a,b,logprocvar,logobsvar,M.n.start,logVar.n.start) {
##   pred <- kfpred(a,b,exp(logprocvar),exp(logobsvar),M.n.start,exp(logVar.n.start),
##                  Nobs=Nobs)
##   -sum(dnorm(Nobs,pred$mean,sqrt(pred$var),log=TRUE))
## }
## ## simulate data:
## set.seed(1001)
## nt <- 200
## obsvar <- 10
## procvar <- 20
## a <- 0.5
## b <- 1
## N <- numeric(nt)
## Nobs <- numeric(nt)
## N[1] <- 150
## for (t in 2:nt) {
##   N[t] <- a +b*N[t-1]+rnorm(1,sd=sqrt(procvar))
## }
## tvec <- 1:nt
## Nobs <- N+rnorm(nt,sd=sqrt(obsvar))
## 
## ## kfpred(2,0.99,procvar=procvar,obsvar=obsvar,M.n.start=Nobs[1],Var.n.start=1,Nobs=Nobs)
## kflik(2,1,logprocvar=log(procvar),logobsvar=log(obsvar),
##       M.n.start=N[1],logVar.n.start=log(procvar+obsvar))
## m1 <- mle2(minuslogl=kflik,
##           start=list(a=1,b=0.9,logprocvar=0,logobsvar=0,
##             M.n.start=150,logVar.n.start=0),
##           fixed=list(b=1),
##           method="Nelder-Mead",
##           control=list(maxit=1000))
## pred1 <- do.call(kfpred,c(as.list(trcoef(coef(m1))),list(Nobs=Nobs)))
## plot(pred1$mean,type="l")
## points(Nobs)
## matlines(1:length(Nobs),
##          pred1$mean+2*outer(sqrt(pred1$var),c(-1,1),"*"),lty=2,col=1)
## plot(pred1$mean,Nobs)
## ## get silly answers, at first ...
## ## need to limit abs(b)<1, var params > 0
## 
## -logLik(m1)
## confint(m1,method="hessian")


###################################################
### chunk number 25:  eval=FALSE
###################################################
## m1prof <- profile(m1,which="logprocvar",trace=TRUE)
## newvals <- list(a=33.067243,b=1.00000,logprocvar=3.810141,logobsvar=-31.542905,
##                 M.n.start=151.222582,logVar.n.start=-85.620510)
## do.call("kflik",newvals)
## ## something wacky going on inside -- reports new better minimum in
## ##   profiling, yet reported parameters don't do better in outside
## ##   scope ???
## ##
## a.prof <- profile(m1,which="a")
## pv <- a.prof@profile$a$par.vals
## z  <- a.prof@profile$a$z
## a.ci <- confint(m1,parm="a")
## plot(pv[,"a"],z)
## abline(v=a.ci)
## abline(v=confint(m1,method="hessian",parm="a"),col=2)
## plot(pv[,"a"],z,type="b")
## abline(v=a.ci)
## abline(v=confint(m1,method="hessian",parm="a"),col=2)
## ## agrees well enough with the Hessian-based CI
## ## over this range what happens
## rrange <- range(which(abs(z)<2.5))
## rvals <- rrange[1]:rrange[2]
## matplot(pv[rvals,"a"],exp(pv[rvals,c("logobsvar","logprocvar")]),type="l")
## plot(pv[rvals,"a"],exp(pv[rvals,"logprocvar"]),type="b")
## matplot(pv[rvals,"a"],pv[rvals,-1],ylim=c(-10,20),type="b")
## 
## m2 <- mle2(minuslogl=kflik,
##           start=list(a=1,b=0.9,logprocvar=0,logobsvar=0,
##             M.n.start=150,logVar.n.start=0),
##           fixed=list(b=1,M.n.start=150,logVar.n.start=-2),
##           method="Nelder-Mead",
##           control=list(maxit=1000))
## 
## m3 <- mle2(minuslogl=kflik,
##           start=list(a=1,b=0.9,logprocvar=0,logobsvar=0,
##             M.n.start=150,logVar.n.start=0),
##           fixed=list(b=1,logprocvar=-30,M.n.start=150,logVar.n.start=-2),
##           method="Nelder-Mead",
##           control=list(maxit=1000))


###################################################
### chunk number 26: 
###################################################
nlkfpred = function(r,K,procvar,obsvar,M.n.start,Var.n.start,Nobs) {
  nt = length(Nobs)
  M.nobs = numeric(nt)
  Var.nobs = numeric(nt)
  M.n = M.n.start
  Var.n = Var.n.start
  M.nobs[1] = M.n.start
  Var.nobs[1] = Var.n.start+obsvar
  for (t in 2:nt) {
    M.ni = M.n+r*M.n*(1-M.n/K)
    b = 1+r-2*r*M.n/K
    Var.ni = b^2*Var.n + procvar
    M.nobs[t] = M.ni
    Var.nobs[t] = Var.ni + obsvar
    M.n = M.ni +  Var.ni/Var.nobs[t]*(Nobs[t]-M.nobs[t])
    Var.n = Var.ni*(1 - Var.ni/Var.nobs[t])
  } 
  list(mean=M.nobs,var=Var.nobs)
}
nlkflik = function(logr,logK,logprocvar,logobsvar,logM.n.start,logVar.n.start,obs.data) {
  pred = nlkfpred(r=exp(logr),K=exp(logK),
                   procvar=exp(logprocvar),obsvar=exp(logobsvar),
                   M.n.start=exp(logM.n.start),Var.n.start=exp(logVar.n.start),
                   Nobs=y.procobs2)
  -sum(dnorm(obs.data,mean=pred$mean,sd=sqrt(pred$var),log=TRUE))
}


###################################################
### chunk number 27: 
###################################################
rval = 0.25
n0val = 3
procvar = 0.5
obsvar = 0.5
nt = 100
newdata2 = simlogistdata(tot.t=nt,sd.proc=sqrt(0.5),sd.obs=sqrt(0.5),
  n0=3,r=0.25,seed=1002,dt=1)
y.procobs2 = newdata2[,2]


###################################################
### chunk number 28: 
###################################################
##plot(y.procobs2,ylim=c(-1,12))




###################################################
### chunk number 29:  eval=FALSE
###################################################
## L1 = do.call("nlkfpred",c(trcoef(startvec),list(Nobs=y.procobs2)))
## head(cbind(NL1$mean,NL1$var,y.procobs2),20)
## do.call("nlkflik",startvec)


###################################################
### chunk number 30: 
###################################################
startvec = list(logr=log(0.25),logK=log(10),logprocvar=log(0.5),logobsvar=log(0.5),
                 logM.n.start=log(3),logVar.n.start=-2)


###################################################
### chunk number 31: 
###################################################
m4 = mle2(minuslogl=nlkflik,
  start=startvec,data=list(obs.data=y.procobs2),
  method="Nelder-Mead",
  control=list(maxit=2000))


###################################################
### chunk number 32: 
###################################################
ucsim = replicate(1000,with(as.list(trcoef(coef(m4))),
  simlogistdata(tot.t=100,sd.proc=sqrt(procvar),sd.obs=sqrt(obsvar),
                n0=M.n.start,r=r,K=K,seed=NULL,dt=1)[,2]))
ucsim[ucsim<0] = NA  ## clear out negative values
ucsim.mean = rowMeans(ucsim,na.rm=TRUE)
ucsim.env = t(apply(ucsim,1,quantile,c(0.025,0.975),na.rm=TRUE))
##  save("ucsim.mean","ucsim.env","m4",file="nlkf.RData")
##} else load("nlkf.RData")


###################################################
### chunk number 33: 
###################################################
restab = cbind(trcoef(unlist(startvec)),trcoef(coef(m4)),exp(confint(m4,method="quad")))
## ignoring starting values since they're kind of funky anyway
restab <- restab[1:4,]
dimnames(restab) = list(
          c("$r$","$K$","$\\ssqproc$","$\\ssqobs$"),
          c("true","fitted","2.5\\%","97.5\\%"))
latex(round(restab[1:4,],2),file="",title="",table.env=FALSE)


###################################################
### chunk number 34: 
###################################################
op = par(cex=2,lwd=2,mar=c(4,4,2,1)+0.1,mfrow=c(1,2),las=1,bty="l")
pred2 = do.call(nlkfpred,c(as.list(trcoef(coef(m4))),list(Nobs=y.procobs2)))
plot(pred2$mean,type="l",ylim=c(0,15),xlab="Time",ylab="Population density",axes=FALSE)
axis(side=2)
axis(side=1,at=c(0,50,100))
box()
points(y.procobs2,cex=0.5)
matlines(1:length(y.procobs2),
         pred2$mean+2*outer(sqrt(pred2$var),c(-1,1),"*"),lty=2,col=1)
matlines(1:length(y.procobs2),
         cbind(ucsim.mean,ucsim.env),col="darkgray",lty=c(1,2,2))
### sweeeeeeet.
legend(25,7,
       c("obs","pred.","unconditional"),
       pch=c(16,NA,NA),
       lty=c(NA,1,1),
       col=c("black","black","gray"),
       merge=TRUE)
transf = function(x) exp(x)
par(mar=c(4,4.5,2,1)+0.1)
plot(transf(ellipse(vcov(m4)[3:4,3:4],centre=coef(m4)[3:4])),type="l",
     xlab=expression(sigma[var]^2),
     ylab=expression(sigma[obs]^2),
     ylim=c(0.2,1.1))
points(transf(coef(m4)[3]),transf(coef(m4)[4]))
text(transf(coef(m4)[3]),transf(coef(m4)[4])+0.04,"estimated")
points(transf(startvec[[3]]),transf(startvec[[4]]),pch=16)
text(transf(startvec[[3]]),transf(startvec[[4]])+0.04,"true")
text(0.28,1.09,"95% conf. int")
## abline(v=exp(confint(m4,method="quad")[3,]),lty=2)
##
## curve(1/14*1/x,add=TRUE)
## curve(1/2.5*1/x,add=TRUE)
## plot(ellipse(vcov(m4),t=1.96,which=3:4,centre=coef(m4)[3:4]),type="l")
## abline(v=confint(m4,method="quad")[3,],lty=2)
par(op)


###################################################
### chunk number 35: 
###################################################
library(R2WinBUGS)
maxr <- 2
n0rate <- 1/y.procobs2[1]


###################################################
### chunk number 36: 
###################################################
o <- y.procobs2
N <- length(y.procobs2)

statespace.data <- list ("N", "o","maxr","n0rate")
## inits <- list(list(n0=y.procobs2[1],r=0.2,K=10,tau.obs=1,tau.proc=1),
##               list(n0=y.procobs2[1],r=0.1,K=10,tau.obs=1,tau.proc=1),
##               list(n0=y.procobs2[1],r=0.4,K=10,tau.obs=1,tau.proc=1),
##               list(n0=y.procobs2[1],r=0.2,K=10,tau.obs=3,tau.proc=1),
##              list(n0=y.procobs2[1],r=0.2,K=10,tau.obs=1,tau.proc=3))
inits = perturb.params(
  list(n0=y.procobs2[1],r=0.2,K=10,tau.obs=1,tau.proc=1),
  alt=list(r=c(0.1,0.4),tau.obs=3,tau.proc=3))


###################################################
### chunk number 37: 
###################################################
parameters <- c("r","K","tau.obs","tau.proc","n0")


###################################################
### chunk number 38:  eval=FALSE
###################################################
## statespace.sim <- bugs(data=statespace.data, inits, param=parameters,
##                        model="statespace.bug", n.chains=length(inits), 
##                        n.iter=15000)
## s1 = as.mcmc.bugs(statespace.sim)


###################################################
### chunk number 39: 
###################################################
library(coda)


###################################################
### chunk number 40: 
###################################################
statespace.sim <- bugs(data=statespace.data, inits, param=parameters,
                       model="statespace.bug", n.chains=length(inits), n.iter=15000,
                       n.burnin=5000,n.thin=37)
s1 = as.mcmc.bugs(statespace.sim)
sum1 = summary(s1)
## save("statespace.sim","s1","sum1",file="statespace.bugsim")


###################################################
### chunk number 41: 
###################################################
gelman.diag(s1)


###################################################
### chunk number 42: 
###################################################
tab1 <- sum1$quantiles[1:4,c(1,3,5)]
tab1[3:4,] <- 1/tab1[3:4,]
tab1[3:4,] <- tab1[4:3,3:1]
rownames(tab1)[3:4] <- c("$\\vproc$","$\\vobs$")
colnames(tab1) <- c("2.5\\%","median","97.5\\%")
latex(round(tab1,2),file="",title="",table.env=FALSE)


###################################################
### chunk number 43: 
###################################################
library(MASS)
vproc <- 1/statespace.sim$sims.list$tau.proc
vobs <- 1/statespace.sim$sims.list$tau.obs
##plot(vproc,vobs)
k1 = kde2d(vproc,vobs,n=60,lims=c(0,1,0.1,1.1))
k1z = cumsum(rev(sort(k1$z)))/sum(k1$z)
w <- which.min(abs(k1z-0.95))
critlevel <-rev(sort(k1$z))[w]
w2 <- which.min(abs(k1z-0.5))
critlevel2 <-rev(sort(k1$z))[w2]
## double-check:
## sum(rev(sort(k1$z))[(w+1):(length(k1$z))])/sum(k1$z)
transf <- function(x) exp(x)


###################################################
### chunk number 44: 
###################################################
op <- par(cex=1.5,lwd=2,mar=c(4,5,2,1)+0.1,las=1,mgp=c(2.5,1,0))
plot(transf(ellipse(vcov(m4)[3:4,3:4],centre=coef(m4)[3:4])),type="l",
     xlab=expression(sigma[proc]^2),
     ylab=expression(sigma[obs]^2),ylim=c(0,1.1),xlim=c(0,1),col="gray")
points(transf(coef(m4)[3]),transf(coef(m4)[4]),col="gray",pch=16)
#text(transf(coef(m4)[3]),transf(coef(m4)[4])+0.04,col="gray","KF")
points(transf(startvec[[3]]),transf(startvec[[4]]),pch=16)
text(transf(startvec[[3]]),transf(startvec[[4]])+0.04,"true")
contour(k1,level=critlevel,drawlabels=FALSE,add=TRUE)
contour(k1,level=critlevel2,drawlabels=FALSE,add=TRUE,lty=2)
points(mean(vproc),mean(vobs))
text(mean(vproc),mean(vobs)-0.05,"mean")
##segments(0.44,0.44,mean(vproc),mean(vobs)-0.02)
maxind <- which(k1$z==max(k1$z),TRUE)
##points(median(vproc),median(vobs),pch=2)
points(k1$x[maxind[1]],k1$y[maxind[2]],pch=3)
text(0.3,0.58,"mode")
text(0.67,0.116,"95% cred. int.")
text(0.22,0.74,"50%")
par(op)


###################################################
### chunk number 45:  eval=FALSE
###################################################
## par(mgp=c(2.5,1,0),mfrow=c(1,2))
## plot(1/as.numeric(s1[,"tau.obs"]),1/as.numeric(s1[,"tau.proc"]),log="xy",pch=".",
##      xlab=expression(sigma[obs]^2),ylab=expression(sigma[proc]^2))
## plot(1/as.numeric(s1[,"tau.obs"]),as.numeric(s1[,"r"]),log="x",pch=".",
##      ylab="r",xlab=expression(sigma[obs]^2))


###################################################
### chunk number 46: 
###################################################
trellis.par.set(canonical.theme(color=FALSE))
print(densityplot(s1,trace=FALSE))


###################################################
### chunk number 47: 
###################################################
nlkfpred = function(r,K,procvar,obsvar,M.n.start,Var.n.start,Nobs) {
  nt = length(Nobs)
  M.nobs = numeric(nt)
  Var.nobs = numeric(nt)
  M.n = M.n.start
  Var.n = Var.n.start
  M.nobs[1] = M.n.start
  Var.nobs[1] = Var.n.start+obsvar
  for (t in 2:nt) {
    M.ni = M.n+r*M.n*(1-M.n/K)
    b = 1+r-2*r*M.n/K
    Var.ni = b^2*Var.n + procvar
    M.nobs[t] = M.ni
    Var.nobs[t] = Var.ni + obsvar
    M.n = M.ni +  Var.ni/Var.nobs[t]*(Nobs[t]-M.nobs[t])
    Var.n = Var.ni*(1 - Var.ni/Var.nobs[t])
  } 
  list(mean=M.nobs,var=Var.nobs)
}


###################################################
### chunk number 48: 
###################################################
nlkflik = function(logr,logK,logprocvar,logobsvar,logM.n.start,logVar.n.start,obs.data) {
  pred = nlkfpred(r=exp(logr),K=exp(logK),
                   procvar=exp(logprocvar),obsvar=exp(logobsvar),
                   M.n.start=exp(logM.n.start),Var.n.start=exp(logVar.n.start),
                   Nobs=y.procobs2)
  -sum(dnorm(obs.data,mean=pred$mean,sd=sqrt(pred$var),log=TRUE))
}


