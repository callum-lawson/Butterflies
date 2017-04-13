require(R2WinBUGS)
curdir <- getwd()

setwd("C:/Users/Callum/Documents/WinBUGS")
sink("weeklymod.txt")
cat("model{

# Priors
	N ~ dnorm(0,0.001)
	t0 ~ dnorm(0,0.001)
	TE ~ dnorm(0,0.001)
	T ~ dnorm(0,0.001)



# Likelihood
	for (i in 1:R){
		y[i] ~ dpois(n[i])
	
		n[tdash[i]<0] <- 0
		
		# Derived params
		tdash[i] <- t[i]-t0
		X[i] <- pi*tdash[i]/TE
		a <- TE/(pi*T)

		}



	}",fill=TRUE)
sink()

# Bundle data
win.data <- with(small, 
	list(y=count, t=week, R=length(week), pi=pi)
	)

# Inits function
inits <- function(){
	list(mu.group=rnorm(1,-1,1), sigma.group=gamsim(0.8, 10), 
		beta1=rnorm(1,1,1),beta2=rnorm(1,1,1),
		beta3=rnorm(1,-1,1),beta4=rnorm(1,-1,1),
		beta5=rnorm(1,0.5,0.5),
		beta6=rnorm(1,-0.5,0.5),beta7=rnorm(1,0.5,0.5),
		beta8=rnorm(1,-1,0.5),beta9=rnorm(1,-0.5,0.5),
		beta10=rnorm(1,-0.5,0.5)
		)
	} 

# Parameters to estimate
params <- c("N","t0","TE","T")

# MCMC settings
nc <- 3
nb <- 100
ni <- 1100
nt <- 1

# Start Gibbs sampler
bcomb <- bugs(win.data, inits, params, "weeklymod.txt", 
	n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt, 
	debug = F, bugs.directory="C:/WinBUGS14", working.directory=getwd()
	)

setwd(curdir)
