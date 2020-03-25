#LHS <- "LHS-1"
taxas <- read.csv(paste(LHS,"/taxas/taxas.csv", sep=""))

HabProp <- seq(0,100,5)
Conf <- 1:10
f=0.8
fitdem <- data.frame(h=c(), conf=c(), p=c(), b0=c(), b1=c(), d0=c(), d1=c(), K=c(), r=c())

for(i in HabProp)
{
  for(l in Conf)
  {    
    h = formatC(i, width = 3, format = "d", flag = "0")
    conf = formatC(l, width = 2, format = "d", flag = "0")
    arqarea <- file(paste("LHS-1/output/output-0.80_", h, "_", conf, "_1.txt", sep=""))
    area <- readLines(arqarea, n=5)[5]
    area <- as.double(unlist(strsplit(area, " "))[-c(1,2)])
    #area <- c(512*512*0.03*0.03 - sum(area), area)
    close(arqarea)
    
    
    taxash <- taxas[taxas$h==i & taxas$conf==l,]
    frags <- unique(taxash$p)
    for(p in frags)
    {
      taxasp <- taxash[taxash$p==p,]
      lmb <- lm(taxasp$b/taxasp$N ~ taxasp$N)
      lmd <- lm(taxasp$d/taxasp$N ~ taxasp$N)
    
      b0 = coef(lmb)[1]
      b1 = -coef(lmb)[2]
      d0 = coef(lmd)[1]
      d1 = coef(lmd)[2]
      
      r = b0 - d0
      K = (b0-d0)/(b1+d1)
        
      fitdem <- rbind(fitdem, data.frame(h=c(i), f=c(f), conf=c(l), p=c(p), b0=c(b0), b1=c(b1), d0=c(d0), d1=c(d1), K=c(K), r=c(r), area=c(area[p])))
    }
    rm(taxash)
  }
}

aicdem <- data.frame(conf=c(), r=c(), K=c())
modelsdem <- list()

lmK <- lm(fitdem$K ~ fitdem$area)
lmr <- lm(fitdem$r ~ fitdem$area)
modelsdem <- list(r = lmr, K = lmK)
aicdem <- rbind(aicdem, data.frame(f=c(f), r=c(AIC(lmr)), K=c(AIC(lmK))))

write.csv(fitdem, file=paste(LHS,"/fit/fitdem.csv", sep=""))
write.csv(aicdem, file=paste(LHS,"/fit/aicdem.csv", sep=""))
save(modelsdem, file=paste(LHS,"/fit/modelsdem.rdata", sep=""))