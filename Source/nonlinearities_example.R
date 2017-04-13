#############
# Functions #
#############

varplot <- function(nlfunc){

	blankplot <- function(x,y,...){
		plot(x,y,
			type="l",bty="n",
			xaxt="n",yaxt="n",
			ann=F,
			col=distcol,
			...)
		}

	### Graphical params

	meancol <- "red"
	distcol <- "blue"
	dashw <- 1.5

	envmean <- 10
	envsd <- 5

	xmar <- c(4,1)
	ymar <- c(3,2)

	xdist <- seq(0,20,length.out=1000)
	xdens <- dnorm(xdist,envmean,envsd)
	ydist <- nlfunc(xdist)
	ydens <- nlfunc(xdens)/sum(nlfunc(xdens))
	mean_y <- sum(ydens*ydist)

	xlim <- range(xdist)
	ylim <- c(0,20)

	### Begin figure

	# X hist
	par(mar=c(1,xmar[1],0,xmar[2])+0.1)
	blankplot(xdist,xdens,xlim=xlim)
	abline(v=envmean,col=meancol,lty=2,lwd=dashw)

	# Y hist
	par(mar=c(ymar[1],1,ymar[2],0)+0.1)
	blankplot(ydens,ydist,
		xlim=rev(range(ydens)),
		ylim=ylim
		)
	abline(h=nlfunc(envmean),col=meancol,lty=2,lwd=dashw)
	abline(h=mean_y,col=distcol,lty=2,lwd=dashw)

	# Central plot
	par(mar=c(ymar[1],xmar[1],ymar[2],xmar[2])+0.1)
	curve(nlfunc(x),
		xlim=xlim,ylim=ylim,
		xaxt="n",yaxt="n",
		ann=F,
		bty="n"
		)
	lines(c(envmean,envmean),c(nlfunc(envmean),-10000),xpd=T,
		col=meancol,lty=2,lwd=dashw)
	lines(c(-10000,envmean),c(nlfunc(envmean),nlfunc(envmean)),
		xpd=T,col=meancol,lty=2,lwd=dashw)
	axis(1,labels=F,ann=F,xaxp=c(xlim[1],xlim[2],1))
	mtext("climate",side=1,line=1)
	axis(2,labels=F,ann=F,yaxp=c(ylim[1],ylim[2],1))
	mtext(expression(log(lambda)),side=2,line=1)
	
	par(mar=c(5,4,4,2)+0.1)

	}

addlet <- function(mytext){
	mtext(mytext,side=3,adj=0.1,line=1)
	}

##########
# Figure #
##########

expfunc <- function(x) exp(x*0.15)
MMfunc <- function(x) 20*x/(2.5+x)
linear <- function(x) x

mymat <- matrix(
	c(rep(c(2,rep(3,4)),4),0,rep(1,4)),
	nr=5,nc=5,byrow=T)
addmat <- function(mymat){
	mymatadd <- mymat+3
	mymatadd[mymatadd==3] <- 0
	mymatadd
	}

mymat2 <- addmat(mymat)
mymat3 <- addmat(mymat2)
mymat4 <- matrix(rep(max(mymat3)+1,5^2),nc=5,nr=5)
fullmat <- rbind(cbind(mymat,mymat2),cbind(mymat3,mymat4))

setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS")
pdf("nl_fig1.pdf",width=6,height=6)

layout(fullmat)
varplot(linear)
addlet("(a) linear")
varplot(expfunc)
addlet("(b) convex")
varplot(MMfunc)
addlet("(c) concave")

plot(1,1,type="n",ann=F,xaxt="n",yaxt="n",bty="n")
legend("center",lty=rep(2,2),lwd=1.5,
	col=c("red","blue"),
	legend=c("Mean climate","Variable climate"), 
	title="Predicted mean",
	cex=1.5,
	bty="n"
	)

dev.off()



