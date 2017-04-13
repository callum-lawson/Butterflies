#############################################################################
# Source code belonging to the paper: van de Pol, M. & Cockburn, A. (2011). #
# "Identifying the critical climatic window that affects trait expression"  #
# The American Naturalist MS 52579                                          #
# For more details see Appendix B of the paper, which also gives a          #
# link to the accompanying datafiles                                        #
#############################################################################

setwd("C:/Users/Callum/Documents/Dropbox/PhD/Analyses/BMS")

# load the libraries 'survival' and 'evd' that are assumed to be installed
library(survival)
library(evd)

# read input data
trait <- read.table(file="trait.txt", header=TRUE)
weather <- read.table(file="weather.txt", header=TRUE)

intervals <- 365 # number of time intervals per year (e.g. 365 days, 52 weeks)
period <- 365 # intervals in the past that can be included in weather index.
index <- array(data=NA, dim=c(length(weather$Precipitation)-period), dimnames=NULL)
# Note that weather data spans a longer period than traitdata, because we correlate
# trait expression data with past weather (e.g. we look back 365 days / 52 weeks maximum).

# first run a model without weather index that can be used as a reference
null_model <- coxph(Surv(trait$Time-1,trait$Time,trait$Event) ~ factor(trait$Age)
	+ frailty(trait$IndividualID), data=trait)

#############################################################
#    WEIBULL method, see eqn. 5 in van de Pol & Cockburn    #
#############################################################

### OBJECT DESCRIPTIONS
#
# intervals = number of time intervals per year (e.g. 365 days, 52 weeks)
# period = intervals in the past that can be included in weather index
# theta, delta, gamma = parameters of Weibull distribution
# y = binary variable indicating when event occured
# traitdata = dataframe giving individual ID, year, time, age
# weatherdata = daily weather variable (vector)
# j = normalised index of days from min to max of period
# weight = distribution to weight climate effects
# index = summed weighted weather values, starting period at each possible day
# (index[1] starts at day 1, index[2] starts at day 2, etc)


# define a function that returns the -LogLikelihood of the extended Cox Model
ModelLogLikelihood<-function(par, y, traitdata, weatherdata) {
  theta <- par[1]     # shape parameter
  gamma <- par[2]     # scale parameter
  delta <- par[3]     # location parameter

  j <- seq(1:period)/period # scale j to interval [0,1]

  # calculate weights using the Weibull probability distribution function
  weight <- theta*gamma*(j-delta)^(theta-1)*exp(-gamma*(j-delta)^theta)

  # calculate index from weather data
  for (a in 1:(length(weatherdata)-period)){
		# beginning at start and  moving period forward until 
		# ...last weather value sampled is at end of weather data
	index[a] <- sum(weatherdata[(a+1):(a+intervals)]*weight) 
		# weighted index of weather over period a+1 to a+year
	}
  index <- (index-mean(index))/sd(index) 
	# standardize index to mean=0 & standard deviation=1
	# (to help with convergence?)
  
  # plot the weight function and corresponding weather index being evaluated
  par(mfrow=c(3,1))
  plot(weatherdata[(1+intervals):(1+4*intervals)],type="s", ylab="precipitation", xlab="days")
  plot(index[1:(3*intervals)],type="s", ylab="index", xlab="days")
  plot((weight/sum(weight)),type="l", ylab="weight", xlab="days")

  # run regression model using weather index
  model<-coxph(Surv(traitdata[,3]-1, traitdata[,3], y) ~
	index[ traitdata[,3] + ((traitdata[,2]-min(traitdata[,2]))*intervals) ] 
	+ factor(traitdata[,5]) 
	+ frailty(traitdata[,1]), 
	data=traitdata
	)
	# Explanation:
	# traitdata[,2]-min(traitdata[,2] gives years starting at 0 from first year
	# ...so this multiplied by intervals gives number of days since start of first year
	# traitdata[,3] gives time step within year (day)
	# ...so added all together calculates weighted climate for each indiv at each time step
	# factor(traitdata[,5]) includes age as a factor
	# frailty term gives individual effects
  logl<-model$loglik[2]-null_model$loglik[2]
  print(c(theta,gamma,delta,logl))
  return(logl)

  }

# Now use optim() routine to find parameter values of Weibull weight function that
# maximize the likelihood of the above model. Depending on the data structure, the
# complexity of the regression model and sample size, the likelihood landscape can
# vary in ruggedness which may affect the probability of getting stuck on local
# optima. Consequently, it is important to experiment with different initial
# conditions (i.e. par) and tune the optimization routine (e.g. ndeps, parscale).
# For more details, see the R-documentation of the optim() function.

result<-optim(par=c(1,1,-0.01), fn=ModelLogLikelihood, y=trait$Event, traitdata=trait, weatherdata=weather$Precipitation, control=list(ndeps=c(0.001,0.001,0.001), parscale=c(100,100,1)), method="L-BFGS-B", lower=c(0.0001,0.0001,-Inf), upper=c(Inf,Inf,0))

# once the optimization routine is finished, we output the final model
theta<-result$par[1]
gamma<-result$par[2]
delta<-result$par[3]
j<-seq(1:period)/period
weight<-theta*gamma*(j-delta)^(theta-1)*exp(-gamma*(j-delta)^theta)
for (a in 1:(length(weather$Precipitation)-period)) {index[a]<-sum(weather$Precipitation[(a+1):(a+intervals)]*weight) }
index<-(index-mean(index))/sd(index)
model<-coxph(Surv(trait$Time-1, trait$Time, trait$Event)~index[trait$Time +((trait$Year-min(trait$Year))*intervals)]+factor(trait$Age)+frailty(trait$IndividualID), data=trait)
summary(model)
par(mfrow=c(3,1))
plot(weather$Precipitation[(1+intervals):(4*intervals)],type="s", ylab="precipitation", xlab="days")
plot(index[1:(3*intervals)],type="s", ylab="index", xlab="days")
plot((weight/sum(weight)),type="l", ylab="weight", xlab="days")
print(c(theta,gamma,delta, model$loglik[2]-null_model$loglik[2]))
print("finished")

#############################################################
#    GEV method, see eqn. 6 in van de Pol & Cockburn        #
#############################################################
# first run a model without weather index that can be used as a reference
null_model<-coxph(Surv(trait$Time-1,trait$Time,trait$Event) ~factor(trait$Age)+frailty(trait$IndividualID), data=trait)

# write a function that returns the -LogLikelihood of the extended Cox Model
ModelLogLikelihood<-function(par, y, traitdata, weatherdata) {
  kappa<-par[1]     # shape parameter
  lambda<-par[2]    # scale parameter
  mu<-par[3]        # location parameter
  range<-10         # rescaling parameter to adjust the range of j over which the GEV function is calculated
  j<-seq(-range,range,by=(2*range/period))
  # calculate weights using the GEV probability distribution function
  weight<-dgev(j[1:period], loc=mu, scale=lambda, shape=kappa, log = FALSE)
  # the GEV function produces "NA" for some positive or negative values of j if
  # the parameter constraint on kappa, lambda and mu (see main text below eqn. 6)
  # is not satisfied. We put such values to zero.
  weight[is.na(weight)]<-0
  if(sum(weight)==0) weight<-weight+1
  # calculate index from weather data
  for (a in 1:(length(weatherdata)-period)) { index[a]<-sum(weatherdata[(a+1):(a+intervals)]*weight) }
  index<-(index-mean(index))/sd(index) # standardize index to mean=0 mean & standard deviation=1
  # plot the weight function and corresponding weather index being evaluated
  par(mfrow=c(3,1))
  plot(weatherdata[(1+intervals):(4*intervals)],type="s", ylab="precipitation", xlab="days")
  plot(index[1:(3*intervals)],type="s",  ylab="index", xlab="days")
  plot((weight/sum(weight)),type="l", ylab="weight", xlab="days")
  # run regression model using weather index
  model<-coxph(Surv(traitdata[,3]-1, traitdata[,3], y)~index[traitdata[,3] +((traitdata[,2]-min(traitdata[,2]))*intervals)] +factor(traitdata[,5]) +frailty(traitdata[,1]), data=traitdata)
  logl<-model$loglik[2]-null_model$loglik[2]
  print(c(kappa,lambda,mu,logl))
  return(-logl)
}
# Now use optim() routine to find parameter values of the GEV weight function that
# maximizes the likelihood of the above model. Finding the global maximum for the
# GEV is typically harder than for the Weibull because of the asymptotic behaviour
# of the GEV probability distribution function when close to kappa=0.
# Therefore, it is advisable to try a wide range of initial conditions (i.e. par).
# The results from the Weibull method should ideally be used to choose similarly shaped
# GEV functions as initial conditions.
result<-optim(par=c(1,1,0), fn=ModelLogLikelihood, y=trait$Event, traitdata=trait, weatherdata=weatherdata$Precipitation, control=list(ndeps=c(0.001,0.001,0.005), parscale=c(1,1,1)),  method = "L-BFGS-B", lower=c(-Inf,0.001,-Inf), upper=c(Inf,Inf,Inf))

# once the optimization routine is finished interpret the final model
kappa<-result$par[1]
lambda<-result$par[2]
mu<-result$par[3]
range<-10
j<-seq(-range,range,by=(2*range/period))
weight<-dgev(j[1:period], loc=mu, scale=lambda, shape=kappa, log = FALSE)
weight[is.na(weight)]<-0
if(sum(weight)==0) weight<-weight+1
for (a in 1:(length(weather$Precipitation)-period)) { index[a]<-sum(weather$Precipitation[(a+1):(a+intervals)]*weight) }
index<-(index-mean(index))/sd(index)
model<-coxph(Surv(trait$Time-1, trait$Time, trait$Event) ~index[trait$Time+((trait$Year-min(trait$Year))*intervals)] +factor(trait$Age)+frailty(trait$IndividualID), data=trait)
summary(model)
par(mfrow=c(3,1))
plot(weather$Precipitation[(1+intervals):(4*intervals)],type="s", ylab="precipitation", xlab="days")
plot(index[1:(3*intervals)],type="s",  ylab="index", xlab="days")
plot((weight/sum(weight)),type="l", ylab="weight", xlab="days")
print(c(kappa,lambda,mu, model$loglik[2]-null_model$loglik[2]))
print("finished")