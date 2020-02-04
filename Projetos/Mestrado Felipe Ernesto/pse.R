library("pse")

dtnorm <- function(x, mean = 0, sd = 1, a = -Inf, b = Inf)
{
	if(length(x) == 1)
	{
		if(x < a | x > b)
			return(0)
		else
			return(dnorm(x = x, mean = mean, sd = sd)/(pnorm(q = b, mean = mean, sd = sd) - pnorm(q = a, mean = mean, sd = sd)))
	}
	else
	{
		return(sapply(x, dtnorm, mean = mean, sd = sd, a = a, b = b))
	}
}

ptnorm <- function(q, mean = 0, sd = 1, a = -Inf, b = Inf)
{
	if(length(q) == 1)
	{
		if(q < a)
			return(0)
		else if (q >= a & q <= b)
			return((pnorm(q = q, mean = mean, sd = sd) - pnorm(q = a, mean = mean, sd = sd))/(pnorm(q = b, mean = mean, sd = sd) - pnorm(q = a, mean = mean, sd = sd)))
		else if (q > b)
			return(1)
	}
	else
	{
		return(sapply(q, ptnorm, mean = mean, sd = sd, a = a, b = b))
	}
}

qtnorm <- function(p, mean = 0, sd = 1, a = -Inf, b = Inf)
{
	if(length(p) == 1)
	{
		if(p < 0 | p > 1)
			warning("p não pode ser menor que zero nem maior que 1.")
		else
			return(qnorm(p = p*(pnorm(q = b, mean = mean, sd = sd) - pnorm(q = a, mean = mean, sd = sd)) + pnorm(q = a, mean = mean, sd = sd), mean = mean, sd = sd))
	}
	else
	{
		return(sapply(p, qtnorm, mean = mean, sd = sd, a = a, b = b))
	}
}

A <- 0.6 #km2
mulamb0 <- 1.26
sigmalamb0 <- 0.44
mut <- 1.89
sigmat <- 0.21
muK <- 3.24
sigmaK <- 0.86
mund <- 1.31
sigmand <- 1.92
maxnd <- log(3)
rholamb0t <- -0.11
rholamb0K <- 0
rholamb0nd <- -0.18
rhotK <- 0.3
rhotnd <- -0.254
rhoKnd <- 0.79

#factors <- c("R", "s", "delta", "lamb0", "mu0", "b", "d", "dm")
factors <- c("lamb0", "t", "K", "nd", "angvis", "u", "v", "w", "dm", "mm")
q <- c("qnorm", "qnorm", "qnorm", "qtnorm", "qunif", "qunif", "qunif", "qunif", "qunif", "qunif")
q.arg <- list(list(mean=mulamb0, sd=sigmalamb0), list(mean=mut, sd=sigmat), list(mean=muK, sd=sigmaK), list(mean=mund, sd=sigmand, b=maxnd), list(min=15, max=360), list(min=1, max=5), list(min=1, max=50), list(min=1, max=50), list(min=1, max=10), list(min=1, max=10))

COR <- matrix(c(1, rholamb0t, rholamb0K, rholamb0nd, 0, 0, 0, 0, 0, 0,
								rholamb0t, 1, rhotK, rhotnd, 0, 0, 0, 0, 0, 0,
								rholamb0K, rhotK, 1, rhoKnd, 0, 0, 0, 0, 0, 0,
								rholamb0nd, rhotnd, rhoKnd, 1, 0, 0, 0, 0, 0, 0,
								0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
								0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
								0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
								0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
								0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
								0, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow = 10, ncol = 10)

opts <- list(COR)
myLHS <- LHS(model=NULL, factors, N=50, q=q, q.arg=q.arg, opts = opts, nboot=0)
LHSdata <- get.data(myLHS)

t <- exp(LHSdata$t)
K <- exp(LHSdata$K)
nd <- exp(LHSdata$nd)
u <- LHSdata$u
v <- LHSdata$v
w <- LHSdata$w

R <- sqrt(A/(pi*K))
s <- u*R
delta <- 4*nd^2/(s^2*t*pi)
lamb0 <- exp(LHSdata$lamb0)
mu0 <- 1/t
b <- lamb0*pi*R^2/v
d <- mu0*pi*R^2/w
angvis <- LHSdata$angvis
dm <- LHSdata$dm
mm <- LHSdata$mm

LHSTWoLife <- data.frame(R=R, s=s, delta=delta, lamb0=lamb0, mu0=mu0, b=b, d=d, angvis=angvis, dm=dm, mm=mm)
write.csv(LHSTWoLife, file="Hipercubo.csv")
