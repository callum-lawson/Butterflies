####################################################################################
# Simulations of second-order discrete population dynamics based on BMS parameters #
####################################################################################

m1coefs <- read.csv("Output/secondorder_m1coefs_27Jun2017.csv")
m2coefs <- read.csv("Output/secondorder_m2coefs_27Jun2017.csv")

m1 <- subset(m1coefs,warn==0)
m2 <- subset(m2coefs,warn==0)
ns <- nrow(m1)
nm <- 2 # number of models

# Constant environment ----------------------------------------------------

nt <- 1000

N <- array(dim=c(nt,nm,ns))
N[1,,] <- 0.1
N[2,,] <- 0.1
for(t in 3:nt){
  N[t,1,] <- with(m1, N[t-1,1,]*exp(b0 + b1*log(N[t-1,1,])))
  N[t,2,] <- with(m2, N[t-1,2,]*exp(b0 + b1*log(N[t-1,2,]) + b2*log(N[t-2,2,])))
}

mcols <- c("blue","red")
par(mfrow=c(5,4),mar=c(1,1,1,1),bty="l")
for(i in 1:ns){
  matplot(I(1:nt),log(N[,,i]),type="l",lty=1,col=mcols,main=m1$species[i])
}
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",title="order",legend=c("1","2"),col=mcols,lty=1,bty="n")

# Variable environment ----------------------------------------------------

nt <- 1000
eps <- rnorm(nt,0,0.33)

N <- array(dim=c(nt,nm,ns))
N[1,,] <- 0.1
N[2,,] <- 0.1
for(t in 3:nt){
  N[t,1,] <- with(m1, N[t-1,1,]*exp(b0 + b1*log(N[t-1,1,]) + eps[t]))
  N[t,2,] <- with(m2, N[t-1,2,]*exp(b0 + b1*log(N[t-1,2,]) + b2*log(N[t-2,2,]) + eps[t]))
}

mcols <- c("blue","red")
par(mfrow=c(5,4),mar=c(1,1,1,1),bty="l")
for(i in 1:ns){
  matplot(I(1:nt),log(N[,,i]),type="l",lty=1,col=mcols,main=m1$species[i])
}
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",title="order",legend=c("1","2"),col=mcols,lty=1,bty="n")

par(mfrow=c(1,1))
i <- 6
matplot(I(1:nt),log(N[,,i]),type="l",lty=1,col=mcols,main=m1$species[i])

plot(t(apply(log(N),c(2,3),mean)))
abline(0,1)
plot(t(apply(log(N),c(2,3),sd)))
abline(0,1)
  # means quite similar but variance almost always higher when second-order included
exp(apply(t(apply(log(N),c(2,3),sd)),1,diff))
  # but not *that* much higher - only by around 5-20%
