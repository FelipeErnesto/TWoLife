LHS="LHS-1"
HabProp <- seq(0, 100, 5)
HabFrag <- c(0.80)
Conf <- 1:10
Replica <- 1:5

for(j in HabFrag)
{
	taxas <- data.frame(h = c(), f=c(), conf= c(), rep = c(), p =c (), N = c(), deltat = c(), b = c(), d = c(), i = c(), e = c())
	abund <- data.frame(h = c(), f=c(), conf= c(), rep = c(), p = c(), meanN = c())
	for(i in HabProp)
	{
		for(k in Conf)
		{
			for(l in Replica)
			{
	  		h = formatC(i, width = 3, format = "d", flag = "0")
				f = format(j, nsmall=2)
				conf = formatC(k, width = 2, format = "d", flag = "0")
	  
	      arquivo <- paste(LHS,"/eventos/eventos-output-", f, "_", h, "_", conf, "_", l, ".txt", sep="")
	      x <- read.csv(arquivo)
				
	      frags <- rev(as.integer(names(sort(table(x$Patch[x$Patch>0])))))
							
				
				for(p in frags)
	      {
	        F <- x[x$Patch==p,]
					if(nrow(F)<=1)
					{
						next
					}
	        F["diff"] <- c(-1, diff(F$t))
					F <- F[-1,]
					
					tN <- tapply(F$diff, F$N, sum)
					
					#plot(sort(unique(F$N)), tN) # pra ver a distribuição de tempos entre os Ns
					me <- sum(sort(unique(F$N))*tN/sum(tN))
					dev <- sqrt(sum((tN/sum(tN))*(sort(unique(F$N))-me)^2))
					
					abund <- rbind(abund, data.frame(h = c(i), f=c(j), conf = c(k), rep = c(l), p = c(p), meanN = c(me)))
					
					F <- F[me-1.5*dev <= F$N & F$N <= me+1.5*dev,]
					
					tN <- tapply(F$diff, F$N, sum)
					eventosN <- table(F$N, F$Event)

					
	        ab <- sort(unique(F$N))
					ab <- ab[tN>0]
					
					if("b" %in% levels(F$Event))
						taxab = eventosN[,"b"][tN>0]/tN[tN>0]
					else
						taxab = 0
					if("d" %in% levels(F$Event))
						taxad = eventosN[,"d"][tN>0]/tN[tN>0]
					else
						taxad = 0
					if("i" %in% levels(F$Event))
						taxai = eventosN[,"i"][tN>0]/tN[tN>0]
					else
						taxai = 0
					if("e" %in% levels(F$Event))
						taxae = eventosN[,"e"][tN>0]/tN[tN>0]
					else
						taxae = 0
					
					tN <- tN[tN>0]
					
					minN <- 10
					if(length(ab) >= minN)
					{
						df <- data.frame(h=rep(i, length(ab)), f=rep(j, length(ab)), conf=rep(k,length(ab)), rep=rep(l,length(ab)), p=rep(p, length(ab)), N = ab, deltat = tN, b = taxab, d = taxad, i = taxai, e = taxae )
						taxas <- rbind(taxas, df)
					}
	      }
			}
    }
  }
}
dir.create(paste(LHS, "/taxas", sep=""))
write.csv(file=paste(LHS,"/taxas/taxas.csv", sep=""), x=taxas)
write.csv(file=paste(LHS,"/taxas/abund.csv", sep=""), x=abund)
