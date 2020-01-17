library("pse")

A <- 0.6 #km2
mulamb0 <- 1.26
sigmalamb0 <- 0.44
mut <- 1.89
sigmat <- 0.21
muK <- 3.24
sigmaK <- 0.86
mund <- 1.31
sigmand <- 1.92
rholamb0t <- -0.11
rholamb0K <- 0
rholamb0nd <- -0.18
rhotK <- 0.3
rhotnd <- -0.254
rhoKnd <- 0.79

#factors <- c("R", "s", "delta", "lamb0", "mu0", "b", "d", "dm")
factors <- c("lamb0", "t", "K", "nd", "angvis", "u", "v", "w", "dm")
q <- c("qnorm", "qnorm", "qnorm", "qnorm", "qunif", "qunif", "qunif", "qunif", "qunif")
q.arg <- list(list(mean=mulamb0, sd=sigmalamb0), list(mean=mut, sd=sigmat), list(mean=muK, sd=sigmaK), list(mean=mund, sd=sigmand), list(min=15, max=360), list(min=1, max=5), list(min=1, max=50), list(min=1, max=50), list(min=1, max=10))

COR <- matrix(c(1, rholamb0t, rholamb0K, rholamb0nd, 0, 0, 0, 0, 0,
								rholamb0t, 1, rhotK, rhotnd, 0, 0, 0, 0, 0,
								rholamb0K, rhotK, 1, rhoKnd, 0, 0, 0, 0, 0,
								rholamb0nd, rhotnd, rhoKnd, 1, 0, 0, 0, 0, 0,
								0, 0, 0, 0, 1, 0, 0, 0, 0,
								0, 0, 0, 0, 0, 1, 0, 0, 0,
								0, 0, 0, 0, 0, 0, 1, 0, 0,
								0, 0, 0, 0, 0, 0, 0, 1, 0,
								0, 0, 0, 0, 0, 0, 0, 0, 1), nrow = 9, ncol = 9)

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

LHSTWoLife <- data.frame(R=R, s=s, delta=delta, lamb0=lamb0, mu0=mu0, b=b, d=d, angvis=angvis, dm=dm)
write.csv(LHSTWoLife, file="Hipercubo.csv")
