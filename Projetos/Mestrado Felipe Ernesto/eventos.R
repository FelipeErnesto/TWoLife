library(doMC)
registerDoMC(8)

eventos<-function(arquivo, npop){
	nlines <- as.integer(unlist(strsplit(system(paste("wc -l", arquivo), intern=TRUE), split=" "))[1])
	dados = file(arquivo, "r")
	dpaisagem = readLines(dados, n=9)
	npatches = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])

	IDmax <- 1
	patches1 <- rep(-1, nlines)
	patches2 <- rep(-1, nlines)
	patches3 <- rep(-1, nlines)
	eventost <- as.double(rep(-1, nlines))
	eventosID <- as.integer(rep(-1, nlines))
	eventosEvent <- as.integer(rep(-1, nlines))
	eventosPatch <- as.integer(rep(-1, nlines))
	eventosN <- as.integer(rep(-1, nlines))
	patchPop <- rep(0, npatches+1)
	saidas <- data.frame(t=c(), Patch=c(), ID=c(), nlin=c())

	for(i in 1:npop)
	{
		lin = readLines(dados, n=1)
		lin<-strsplit(lin, " ")
		line <- unlist(lin)
		patches1[IDmax] <- strtoi(line[3])
		patches2[IDmax] <- -1
		patches3[IDmax] <- -1
		IDmax <- IDmax + 1
		patchPop[strtoi(line[3])+1] <- patchPop[strtoi(line[3])+1]+1
	}
	
	nlin <- 1
	lin = readLines(dados, n=1)
	while(lin != "EOF")
	{
		lin<-strsplit(lin, " ")
		line <- unlist(lin)
		acao <-strtoi(line[2])
		p <- strtoi(line[4])
		ID <- strtoi(line[3])
		t <- as.double(line[1])
		ind <- ID
		
		if(acao == 0)
		{
			eventost[nlin] <- t
			eventosID[nlin] <- ID
			eventosEvent[nlin] <- 2
			eventosPatch[nlin] <- p
			eventosN[nlin] <- patchPop[p+1]
			nlin <- nlin + 1

			if(patches1[ind]==0 & patches2[ind]>0)
			{
				saiu <- saidas[saidas$Patch==patches2[ind] & saidas$ID==ID,]
				saidas <- saidas[saidas$Patch!=patches2[ind] | saidas$ID!=ID,]
				eventosEvent[saiu$nlin] <- 4
			}
			
			patchPop[p+1] <- patchPop[p+1] - 1
			
		}

		if(acao == 1)
		{
			eventost[nlin] <- t
			eventosID[nlin] <- ID
			eventosEvent[nlin] <- 1
			eventosPatch[nlin] <- p
			eventosN[nlin] <- patchPop[p+1]
			nlin <- nlin + 1
			patchPop[p+1] <- patchPop[p+1] + 1
			patches1[IDmax] <- p
			patches2[IDmax] <- -1
			patches3[IDmax] <- -1
			IDmax <- IDmax + 1
		}
		if(acao == 2)
		{
			if( p!=patches1[ind]  )
			{
				patches3[ind] <- patches2[ind]
				patches2[ind] <- patches1[ind]
				patches1[ind] <- p

				if(patches1[ind]>0)
				{
					if(patches2[ind] > 0)
					{
						eventost[nlin] <- t
						eventosID[nlin] <- ID
						eventosEvent[nlin] <- 4
						eventosPatch[nlin] <- patches2[ind]
						eventosN[nlin] <- patchPop[patches2[ind]+1]
						nlin <- nlin + 1
						eventost[nlin] <- t
						eventosID[nlin] <- ID
						eventosEvent[nlin] <- 3
						eventosPatch[nlin] <- patches1[ind]
						eventosN[nlin] <- patchPop[patches1[ind]+1]
						nlin <- nlin + 1
					}
					
					else if(patches3[ind] > -1 & patches3[ind]!=patches1[ind])
					{
						saiu <- saidas[saidas$Patch==patches3[ind] & saidas$ID==ID,]
						saidas <- saidas[saidas$Patch!=patches3[ind] | saidas$ID!=ID,]
						eventost[nlin] <- t
						eventosID[nlin] <- ID
						eventosEvent[nlin] <- 3
						eventosPatch[nlin] <- patches1[ind]
						eventosN[nlin] <- patchPop[patches1[ind]+1]
						nlin <- nlin + 1
						eventosEvent[saiu$nlin] <- 4
					}
					else if(patches3[ind]==patches1[ind])
					{
						saidas <- saidas[saidas$Patch!=patches3[ind] | saidas$ID!=ID,]
						eventost[nlin] <- t
						eventosID[nlin] <- ID
						eventosEvent[nlin] <- 6
						eventosPatch[nlin] <- patches1[ind]
						eventosN[nlin] <- patchPop[patches1[ind]+1]
						nlin <- nlin + 1
					}
				}
				else if(patches1[ind]==0)
				{
					saidas <- rbind(saidas, data.frame(t=c(t), Patch=c(patches2[ind]), ID=c(ID), nlin=c(nlin)))
					eventost[nlin] <- t
					eventosID[nlin] <- ID
					eventosEvent[nlin] <- 5
					eventosPatch[nlin] <- patches2[ind]
					eventosN[nlin] <- patchPop[patches2[ind]+1]
					nlin <- nlin + 1
				}
				
				patchPop[patches1[ind]+1] <- patchPop[patches1[ind]+1] + 1
				patchPop[patches2[ind]+1] <- patchPop[patches2[ind]+1] - 1
				
			}
		}

		if(acao == 3)
		{
			if(patches1[ind]>0)
			{
				eventost[nlin] <- t
				eventosID[nlin] <- ID
				eventosEvent[nlin] <- 4
				eventosPatch[nlin] <- patches1[ind]
				eventosN[nlin] <- patchPop[patches1[ind]+1]
				nlin <- nlin + 1
			}
			else if(patches1[ind]==0 & patches2[ind]>0)
			{
				saiu <- saidas[saidas$Patch==patches2[ind] & saidas$ID==ID,]
				saidas <- saidas[saidas$Patch!=patches2[ind] | saidas$ID!=ID,]
				eventosEvent[saiu$nlin] <- 4
			}
			
			patchPop[p+1] <- patchPop[p+1] - 1

		}

		lin = readLines(dados, n=1)
	}
	out <- strsplit(arquivo, split="/")[[1]][2]
	arqab <- paste0("eventos/abundancia-", out)
	arqtaxas <- paste0("eventos/eventos-", out)

	eventosEvent <- ordered(eventosEvent, levels=c(-1,1,2,3,4,5,6))
	levels(eventosEvent) <- c("z", "b", "d", "i", "e", "saiu", "entrou")
	eventos <- data.frame(t = eventost, ID = eventosID, Event = eventosEvent, Patch = eventosPatch, N = eventosN)
	write.csv(eventos[eventos$t >=5,], arqtaxas, row.names=FALSE)
	
	file.create(arqab)
	dfab <- data.frame(Patch=0:npatches, N=patchPop)
	write.csv(dfab, arqab, row.names=FALSE)

	close(dados)
}

for(j in 1:50)
{
	setwd(paste("LHS-", j, sep =""))
	files <- list.files(pattern="output-*", path="output", full.names=T)
	dir.create("eventos")
	foreach(i = files) %dopar%
	{
		eventos(i, 300)
	}
	setwd("..")
}
