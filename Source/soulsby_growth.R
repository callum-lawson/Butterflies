#############
# FUNCTIONS #
#############

soulpred <- function(t,N,t0,TE,T){
	
	#derived params
	tdash <- t-t0
	X <- pi*tdash/TE
	a <- TE/(pi*T)

	n <- vector()

	# when t' is less than 0
	n[tdash<0] <- 0

	# when t' lies between 0 and T[E]
	c2 <- tdash>=0 & tdash<=TE
	A1 <- 3*N/(4*(a^2+9))
	A2 <- a*sin(X[c2])^3 - 3*sin(X[c2])^2 * cos(X[c2])
	A3 <- 6/(a^2+1)
	A4 <- a*sin(X[c2]) - cos(X[c2]) + exp(-a*X[c2])
	n[c2] <- A1 * (A2 + A3 * A4)

	# when t' is greater than T[E]
	c3 <- tdash>TE
	B1 <- 9*N*exp(-a*X[c3])*(1+exp(a*pi))
	B2 <- 2*(a^2+9)*(a^2+1)
 	n[c3] <- B1/B2
	return(n)
	}

soulLik <- function(par, y, t){
	
	# read in params
	N <- par[1]
	t0 <- par[2]
	TE <- par[3]
	T <- par[4]
	
	#calculate n
	n <- soulpred(t,N,t0,TE,T)

	#likelihood
	ndash <- n[!(n<0.1 & t<=t0)] # adjusted to prevent predictions of zero
	ydash <- y[!(n<0.1 & t<=t0)] # (likelihood of zero predictions is 1)
	Lik <- sum(ydash*log(ndash)-ndash)

	return(Lik)
	
	}

###############
# SINGLE SITE #
###############

curdat <- bspl$"Silver-spotted Skipper"

# load curdat as SSSk
( small <- subset(curdat,site==991 & year==2003) )
	
y <- small$count
t <- small$week # may need to change for other species

result1 <- optim(par=c(N=100,t0=16,TE=6,T=1), 
	fn=soulLik, 
	y=small$count, t=small$week, 
	control=list(fnscale=-1))

plot(y~t,type="b",pch=16,las=1)
with(result1,
	lines(soulpred(t,par[1],par[2],par[3],par[4])~t,
		col="red")
	)

##############################
# POOLED SITES - SINGLE YEAR #
##############################

oneyear <- subset(na.omit(curdat),year==2003) # NAs ommitted

result2 <- optim(par=c(N=100,t0=16,TE=6,T=1), 
	fn=soulLik, 
	y=oneyear$count, t=oneyear$week, 
	control=list(fnscale=-1))

predweek <- with(oneyear, min(week):max(week))
preds <- with(result2, soulpred(predweek,par[1],par[2],par[3],par[4]))
wmeans <- with(na.omit(oneyear),tapply(count,week,mean))
se.func <- function(x) sd(x)/sqrt(length(x))
wse <- with(na.omit(oneyear),tapply(count,week,se.func))
require(gplots)
barpos <- barplot2(wmeans, plot.ci=T, 
	ci.u=wmeans+1.96*wse,
	ci.l=wmeans-1.96*wse,
	las=1,
	)
with(oneyear, lines(preds ~ barpos, type="l",col="red"))

##########################
# TRYING OUT CURVE SHAPE #
##########################

predweek <- with(oneyear, min(week):max(week))
predsoul <- function(x,T){
	soulpred(x,100,16,6,T)
	}

Tseq <- seq(0.25,2,0.25)

pdf("emergence_survival_curves.pdf",width=6,height=6)
plot(1,1,type="n",
	ylab="Abundance",
	xlab="Week",
	las=1,
	xlim=c(16,26),
	ylim=c(0,60),
	bty="n"
	)
for(i in 1:length(Tseq)){
	curve(predsoul(x,Tseq[i]),
		xlim=c(16,26),add=T,col=i)
	}
legend("topright",legend=Tseq,
	title="mean lifespan (weeks)",
	lty=rep(1,length(Tseq)),
	col=1:length(Tseq),
	bty="n")
mtext("N=100",side=3,adj=0)
dev.off()

##################
# MULTIPLE SITES #
##################

soulLik.fixed <- function(par, y, t){
	
	# read in params
	N <- par[1]
	t0 <- par[2]
	TE <- par[3]
	T <- par[4]
	
	
	#calculate n
	n <- soulpred(t,N,t0,TE,T)

	#likelihood
	ndash <- n[!(n<0.1 & t<=t0)] # adjusted to prevent predictions of zero
	ydash <- y[!(n<0.1 & t<=t0)] # (likelihood of zero predictions is 1)
	Lik <- sum(ydash*log(ndash)-ndash)
	return(Lik)
	
	}




