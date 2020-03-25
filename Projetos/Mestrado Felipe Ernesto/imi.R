library(doMC)
registerDoMC(10)
library("bbmle")
#LHS <- "LHS-1"
HabProp <- seq(0, 100, 5)
Conf <- 1:10
f=0.8

taxas <- read.csv(paste(LHS,"/taxas/taxas.csv", sep=""))
taxas$X<-NULL
taxas <- taxas[taxas$N>0,]
abund <- read.csv(paste(LHS,"/taxas/abund.csv", sep=""))
abund$X<-NULL
dists <- list()
areas <- list()
emi <- data.frame(h=c(), conf=c(), p=c(), EN=c())
for(i in HabProp)
{
  for(l in Conf)
  {
    h = formatC(i, width = 3, format = "d", flag = "0")
    conf = formatC(l, width = 2, format = "d", flag = "0")
    arquivo <- paste("LHS-1/output/output-0.80_", h, "_", conf, "_1.txt", sep ="")
    dados = file(arquivo, "r")
  	dpaisagem = readLines(dados, n=9)
  	nfrag = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])
    close(dados)
    conf = formatC(l, width = 5, format = "d", flag = "0")
    ds <- scan(paste("Distancias/dists_0.80_", h,"_", conf, ".txt", sep=""))
    ds <- matrix(ds, nrow=nfrag, ncol=nfrag, byrow=TRUE)
    dists[[paste(i, "-", l, sep="")]] <- ds
    
    conf = formatC(l, width = 2, format = "d", flag = "0")
    arqarea <- file(paste("LHS-1/output/output-0.80_", h, "_", conf, "_1.txt", sep=""))
    as <- readLines(arqarea, n=5)[5]
    as <- as.double(unlist(strsplit(as, " "))[-c(1,2)])
    close(arqarea)
    areas[[paste(i, "-", l, sep="")]] <- as  
    
    tax <- taxas[taxas$h==i & taxas$conf==l,]
    ab <- abund[abund$h==i & abund$conf==l,]
    frags <- unique(tax$p)
    for(p in frags)
    {
      E <- mean(tax$e[tax$p==p]/tax$N[tax$p==p])
      N <- mean(ab$meanN[ab$p==p])
      emi <- rbind(emi, data.frame(h=c(i), conf=c(l), p=c(p), EN=c(E*N)))
    }
  }
}


IF1single <- function(linha, alpha, b)
{
  di <- dists[[paste(linha[1], "-", linha[3], sep="")]][linha[5],-(linha[5])]
  ar <- areas[[paste(linha[1], "-", linha[3], sep="")]][-(linha[5])]
  sum(exp(-exp(alpha)*di)*ar^b)
}

IF1 <- function(alpha, b, conf)
{
  if1 <- apply(taxas[taxas$conf==conf,], 1, IF1single, alpha, b)
  return(if1)
}

LL1 <- function(lalpha, b, conf){
    if1 <- IF1(lalpha, b, conf)
    -sum(dpois(taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf], lambda = if1*taxas$deltat[taxas$conf==conf], log=TRUE))
}


LL1k <- function(lalpha, b, k, conf){
    if1k <- k*IF1(lalpha, b, conf)
    -sum(dpois(taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf], lambda = if1k*taxas$deltat[taxas$conf==conf], log=TRUE))
}


IF2single <- function(linha, alpha, b, c)
{
  di <- dists[[paste(linha[1], "-", linha[3], sep="")]][linha[5],-(linha[5])]
  ar <- areas[[paste(linha[1], "-", linha[3], sep="")]][-(linha[5])]
  ai <- areas[[paste(linha[1], "-", linha[3], sep="")]][linha[5]]
  ai^c*sum(exp(-exp(alpha)*di)*ar^b)
}

IF2 <- function(alpha, b, c, conf)
{
  
  if2 <- apply(taxas[taxas$conf==conf,], 1, IF2single, alpha, b, c)
  return(if2)
}

LL2 <- function(lalpha, b, c, conf){
    if2 <- IF2(lalpha, b, c, conf)
    -sum(dpois(taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf], lambda = if2*taxas$deltat[taxas$conf==conf], log=TRUE))
}


IFEsingle <- function(linha, alpha)
{
  EN <- c()
  di <- c()
  ar <- c()
  frags <- unique(emi$p[emi$h==linha[1] & emi$conf==linha[3]])
  for(p in frags)
  {
    if(p!=linha[5])
    {
      EN <- c(EN, emi$EN[emi$h==linha[1] & emi$conf==linha[3] & emi$p==p])
      di <- c(di, dists[[paste(linha[1], "-", linha[3], sep="")]][linha[5], p])
      ar <- c(ar, areas[[paste(linha[1], "-", linha[3], sep="")]][p])
    }
  }
  if(length(EN)>0)
    return(sum(exp(-exp(alpha)*di)*EN))
  else
    return(NA)
}

IFE <- function(alpha, conf)
{
  ife <- apply(taxas[taxas$conf==conf,], 1, IFEsingle, alpha)
  return(ife)
}

LLE <- function(lalpha, conf){
    ife <- IFE(lalpha, conf)
    x <- taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf]
    lambda = ife*taxas$deltat[taxas$conf==conf]
    -sum(dpois(x[is.na(lambda)==FALSE], lambda = lambda[is.na(lambda)==FALSE], log=TRUE))
}

LLEK <- function(lalpha, k, conf){
    ifek <- k*IFE(lalpha, conf)
    x <- taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf]
    lambda <- ifek*taxas$deltat[taxas$conf==conf]
    -sum(dpois(x[is.na(lambda)==FALSE], lambda = lambda[is.na(lambda)==FALSE], log=TRUE))
}

IFE2single <- function(linha, alpha, c)
{
  EN <- c()
  di <- c()
  ar <- c()
  frags <- unique(emi$p[emi$h==linha[1] & emi$conf==linha[3]])
  for(p in frags)
  {
    if(p!=linha[5])
    {
      EN <- c(EN, emi$EN[emi$h==linha[1] & emi$conf==linha[3] & emi$p==p])
      di <- c(di, dists[[paste(linha[1], "-", linha[3], sep="")]][linha[5], p])
      ar <- c(ar, areas[[paste(linha[1], "-", linha[3], sep="")]][p])
    }
  }
  ai <- areas[[paste(linha[1], "-", linha[3], sep="")]][linha[5]]
  if(length(EN)>0)
    return(ai^c*sum(exp(-exp(alpha)*di)*EN))
  else
    return(NA)
}

IFE2 <- function(alpha, c, conf)
{
  ife2 <- apply(taxas[taxas$conf==conf,], 1, IFE2single, alpha, c)
  return(ife2)
}

LLE2 <- function(lalpha, c, conf){
    ife2 <- IFE2(lalpha, c, conf)
    x <- taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf]
    lambda = ife2*taxas$deltat[taxas$conf==conf]
    -sum(dpois(x[is.na(lambda)==FALSE], lambda = lambda[is.na(lambda)==FALSE], log=TRUE))
}

saida <- foreach(conf = 1:2) %dopar%
{
  fit1 <- mle2(LL1, start = list(lalpha=log(1), b=1), fixed = list(conf=conf))
  fit1k <- mle2(LL1k, start = list(lalpha=log(1), b=1, k=2), fixed = list(conf=conf))
  fit2 <- mle2(LL2, start = list(lalpha=log(1), b=1, c=1), fixed=list(conf=conf))
  fite <- mle2(LLE, start = list(lalpha=log(1)), fixed = list(conf=conf))
  fitek <- mle2(LLEK, start = list(lalpha=log(1), k=2), fixed = list(conf=conf))
  fite2 <- mle2(LLE2, start = list(lalpha=log(1), c=1), fixed=list(conf=conf))
  metrica1 <- IF1(coef(fit1)[1], coef(fit1)[2], conf)
  metrica1k <- coef(fit1k)[3]*IF1(coef(fit1k)[1], coef(fit1k)[2], conf)
  metrica2 <- IF2(coef(fit2)[1], coef(fit2)[2], coef(fit2)[3], conf)
  metricaE <- IFE(coef(fite)[1], conf)
  metricaEk <- coef(fitek)[2]*IFE(coef(fitek)[1], conf)
  metricaE2 <- IFE2(coef(fite2)[1], coef(fite2)[2], conf)
  
  tam <- length(taxas$i[taxas$conf==conf])
  # fitimi <- rbind(fitimi, data.frame(f=rep(f, tam), conf=rep(conf, tam), nobs=taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf], metrica1=metrica1*taxas$deltat[taxas$conf==conf], metrica1k=metrica1k*taxas$deltat[taxas$conf==conf], metrica2=metrica2*taxas$deltat[taxas$conf==conf]))
  fitimi <- data.frame(f=rep(f, tam), conf=rep(conf, tam), nobs=taxas$i[taxas$conf==conf]*taxas$deltat[taxas$conf==conf], metrica1=metrica1*taxas$deltat[taxas$conf==conf], metrica1k=metrica1k*taxas$deltat[taxas$conf==conf], metrica2=metrica2*taxas$deltat[taxas$conf==conf],  metricaE=metricaE*taxas$deltat[taxas$conf==conf], metricaEk=metricaEk*taxas$deltat[taxas$conf==conf], metricaE2=metricaE2*taxas$deltat[taxas$conf==conf])
  modelsimi <- list(metrica1 = fit1, metrica1k = fit1k, metrica2=fit2, metricaE=fite, metricaEk=fitek, metricaE2=fite2)
  #aicimi <- rbind(aicimi, data.frame(f=c(f), conf=c(conf), m1trica1=c(AIC(metrica1)), metrica1k=c(AIC(metrica1k)), metrica2=c(AIC(metrica2))))
  aicimi <- AICtab(fit1, fit1k, fit2, fite, fitek, fite2)
  list(fitimi = fitimi, modelsimi = modelsimi, aicimi = aicimi)
}

fitimi <- data.frame(f=c(), conf=c(), nobs=c(), metrica1=c(), metrica1k=c(), metrica2=c(), metricaE=c(), metricaEk=c())
aicimi <- list()
modelsimi <- list()
for(conf in 1:2)
{
  fitimi <- rbind(fitimi, saida[[conf]]$fitimi)
  modelsimi[[toString(conf)]] <- saida[[conf]]$modelsimi
  aicimi[[toString(conf)]] <- saida[[conf]]$aicimi
}
write.csv(fitimi, file=paste(LHS,"/fit/fitimi.csv", sep=""))
save(aicimi, file=paste(LHS,"/fit/aicimi.rdata", sep=""))
save(modelsimi, file=paste(LHS,"/fit/modelsimi.rdata", sep=""))
