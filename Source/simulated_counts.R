
nt <- 20
Tbar <- mean(curdat$temp8)
T <- rep(5.5,nt)
tseq <- 1:nt

### Temperature interaction

alpha <- params[1]
gamma1 <- params[2]
beta2 <- params[3]
beta3 <- params[4]

N <- numeric(nt)
lambda <- numeric(nt)
N[1] <- mean(curdat$count)
Lbar <- mean(curdat$tlen)

for (t in 1:(nt-1)) {
	#alpha <- rep(alpha,nt)
	lambda[t] <- exp(
		alpha
		+ gamma1*log(N[t]/Tbar)			# is it correct to divide by Tbar here?
		+ beta2*T[t]
		+ beta3*log(N[t]/Tbar)*T[t]
		)
	N[t+1] <- rpois(1,lambda=N[t]*lambda[t])
	}

plot(N~tseq)

### Random model


nt <- 20
Tbar <- mean(curdat$temp8)
T <- rep(5.5,nt)
tseq <- 1:nt

### Temperature interaction

params <- coef(mod1g)
alpha <- params[1]
gamma1 <- params[2]
beta2 <- params[3]
beta3 <- params[4]

N <- numeric(nt)
lambda <- numeric(nt)
N[1] <- mean(curdat$count)
Lbar <- mean(curdat$tlen)

for (t in 1:(nt-1)) {
	#alpha <- rep(alpha,nt)
	lambda[t] <- exp(
		alpha
		+ gamma1*log(N[t]/Tbar)			# is it correct to divide by Tbar here?
		+ beta2*T[t]
		+ beta3*log(N[t]/Tbar)*T[t]
		)
	N[t+1] <- rpois(1,lambda=N[t]*lambda[t])
	}

