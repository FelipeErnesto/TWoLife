# Talvez não precise fitar nenhuma função, talvez só usar os dados da
# simulação pra pegar os valores de E e I baseado em atributos do fragmento
# seja o suficiente para validar o modelo, e dps ter a cara das métricas
# às quais N* responde

#LHS <- "LHS-1"
library("bbmle")
library("lme4")

HabProp <- seq(0, 100, 5)
Conf <- 1:10
f = 0.8
areas <- list()
for(i in HabProp)
{
  for(l in Conf)
  {
      h = formatC(i, width = 3, format = "d", flag = "0")
      conf = formatC(l, width = 2, format = "d", flag = "0")
      arqarea <- file(paste("LHS-1/output/output-0.80_", h, "_", conf,"_1.txt", sep=""))
      as <- readLines(arqarea, n=5)[5]
      as <- as.double(unlist(strsplit(as, " "))[-c(1,2)])
      close(arqarea)
      areas[[paste(i, "-", l, sep="")]] <- as
  }
}

taxas <- read.csv(paste(LHS,"/taxas/taxas.csv", sep=""))
taxas$X<-NULL
taxas <- taxas[taxas$N>0,]
areaemi <- rep(-1, length(taxas$h))

for(i in 1:length(taxas$h))
{
  line <-  taxas[i,]
  areaemi[i] <- areas[[paste(line$h, "-", line$conf, sep="")]][line$p]
}

fitemi <- data.frame(f=taxas$f, conf=taxas$conf, E = taxas$e/taxas$N, area = areaemi)

#aicemi <- data.frame(conf=c(), modl=c(), modh=c(), mode=c(), modlog=c())
modelsemi <- list()

modl <- lm(fitemi$E ~ fitemi$area)
modlog <- lm(fitemi$E ~ log(fitemi$area))
modh <- nls(E ~ a*area^b, start = list(a=1,b=1), data=fitemi)
mode <- nls(E ~ a*exp(-b*area), start = list(a=1,b=1), data=fitemi)
#modelsemi[[toString(conf)]] <- list(modl = modl, modh = modh, mode=mode, modlog=modlog)
modelsemi <- list(modl = modl, modh = modh, mode=mode, modlog=modlog)
#aicemi <- rbind(aicemi, data.frame(f=c(f), conf=c(conf), modl=c(AIC(modl)), modh=c(AIC(modh)), mode=c(AIC(mode))))
aicemi <- AICtab(modl, modlog, modh, mode)
class(aicemi) <- "data.frame"

# mode2 <- glm(E ~ area, data=fitemi, family=Gamma(link="log"))
# modlog2 <- glm(E ~ log(area), data=fitemi, family=Gamma(link="identity"))
# modh2 <- glm(E ~ log(area), data=fitemi, family=Gamma(link="log"))
# modl2 <- glm(E ~ area, data=fitemi, family=Gamma(link="identity"))
# modelsemi2 <- list(modl2 = modl2, modh2 = modh2, mode2=mode2, modlog2=modlog2)
# aicemi2 <- AICtab(modl2, modlog2, modh2, mode2)
# class(aicemi2) <- "data.frame"

write.csv(fitemi, file=paste(LHS,"/fit/fitemi.csv", sep=""))
write.csv(aicemi, file=paste(LHS,"/fit/aicemi.csv", sep=""))
save(modelsemi, file=paste(LHS,"/fit/modelsemi.rdata", sep=""))

# write.csv(aicemi2, file=paste(LHS,"/fit/aicemi2.csv", sep=""))
# save(modelsemi2, file=paste(LHS,"/fit/modelsemi2.rdata", sep=""))