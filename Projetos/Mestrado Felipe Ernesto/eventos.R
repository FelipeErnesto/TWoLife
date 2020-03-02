library(doMC)
registerDoMC(8)

eventos<-function(arquivo, npop){
	dados = file(arquivo, "r")
	dpaisagem = readLines(dados, n=9)
	npatches = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])

	IDlist<-rep(0, npop)
	patches <- matrix( rep(0, 3*npop), nrow = npop, ncol = 3 )
	eventos <- data.frame(t=double(), ID=integer(), Event = factor(), Patch = integer(), N = integer()) # Data frame que armezana os eventos relacionados às quatro taxas
	levels(eventos$Event) <- c("b", "d", "i", "e", "saiu", "entrou")
	patchPop <- rep(0, npatches+1)
	saidas <- data.frame(t=c(), Patch=c(), ID=c())

	for(i in 1:npop)
	{
		lin = readLines(dados, n=1)
		lin<-strsplit(lin, " ")
		line <- unlist(lin)
		patches[i,1]<-strtoi(line[3])
		patches[i,2] <- -1
		patches[i,3] <- -1
		IDlist[i] <- strtoi(line[2])
		patchPop[strtoi(line[3])+1] <- patchPop[strtoi(line[3])+1]+1
	}
	maxID = npop

	lin = readLines(dados, n=1)
	while(lin != "EOF")
	{
		lin<-strsplit(lin, " ")
		line <- unlist(lin)
		acao <-strtoi(line[2])
		p <- strtoi(line[4])
		ID <- strtoi(line[3])
		t <- as.double(line[1])
		ind <- match(ID, IDlist)
		
		if(acao == 0)
		{	
			eventos <- rbind(eventos, data.frame(t=c(t), ID=c(ID), Event = factor(c("d"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(p), N = c(patchPop[p+1])))	
			
			if(patches[ind,1]==0 & patches[ind,2]>0)
			{
				saiu <- saidas[saidas$Patch==patches[ind,2] & saidas$ID==ID,]
				saidas <- saidas[saidas$Patch!=patches[ind,2] | saidas$ID!=ID,]
				eventos[eventos$t == saiu$t & eventos$ID == ID & eventos$Event == "saiu" & eventos$Patch == patches[ind,2], ]$Event <- "e"
			}
			
			patchPop[p+1] <- patchPop[p+1] - 1
			
			if(length(IDlist)>2)
			{
				patches<-patches[-ind,]
				IDlist<-IDlist[-ind]
			}
			
			else
			{
					patches<-patches[-ind,]
					patches<-matrix(patches, nrow=1, ncol=3)
					IDlist<-IDlist[-ind]
			}
		}

		if(acao == 1)
		{
			maxID<-maxID+1
			eventos <- rbind(eventos, data.frame(t=c(t), ID=c(ID), Event = factor(c("b"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(p), N = c(patchPop[p+1])))
			patchPop[p+1] <- patchPop[p+1] + 1
			patches<-rbind(patches, c(p, -1, -1))
			IDlist<- c(IDlist, maxID)
		}
		if(acao == 2)
		{
			if( p!=patches[ind, 1]  )
			{
				patches[ind,3] <- patches[ind,2]
				patches[ind,2] <- patches[ind,1]
				patches[ind,1] <- p

				if(patches[ind,1]>0)
				{
					if(patches[ind,2] > 0)
					{
						eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("e"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,2]), N = c(patchPop[patches[ind,2]+1])))
						eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("i"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,1]), N = c(patchPop[patches[ind,1]+1])))
					}
					
					else if(patches[ind,3] > -1 & patches[ind,3]!=patches[ind,1])
					{
						saiu <- saidas[saidas$Patch==patches[ind,3] & saidas$ID==ID,]
						saidas <- saidas[saidas$Patch!=patches[ind,3] | saidas$ID!=ID,]
						eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("i"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,1]), N = c(patchPop[patches[ind,1]+1])))
						eventos[eventos$t == saiu$t & eventos$ID == ID & eventos$Event == "saiu" & eventos$Patch == patches[ind,3], ]$Event <- "e"
					}
					else if(patches[ind,3]==patches[ind,1])
					{
						saidas <- saidas[saidas$Patch!=patches[ind,3] | saidas$ID!=ID,]
						eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("entrou"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,1]), N = c(patchPop[patches[ind,1]+1])))
					}
				}
				else if(patches[ind,1]==0)
				{
					saidas <- rbind(saidas, data.frame(t=c(t), Patch=c(patches[ind,2]), ID=c(ID)))
					eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("saiu"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,2]), N = c(patchPop[patches[ind,2]+1])))
				}
				
				patchPop[patches[ind,1]+1] <- patchPop[patches[ind,1]+1] + 1
				patchPop[patches[ind,2]+1] <- patchPop[patches[ind,2]+1] - 1
				
			}
		}

		if(acao == 3)
		{
			if(patches[ind,1]>0)
			{
				eventos <- rbind(eventos, data.frame(t=c(t), ID = c(ID), Event = factor(c("e"), levels=c("b", "d", "i", "e", "saiu", "entrou")), Patch = c(patches[ind,1]), N = c(patchPop[patches[ind,1]+1])))
			}
			else if(patches[ind,1]==0 & patches[ind,2]>0)
			{
				saiu <- saidas[saidas$Patch==patches[ind,2] & saidas$ID==ID,]
				saidas <- saidas[saidas$Patch!=patches[ind,2] | saidas$ID!=ID,]
				eventos[eventos$t == saiu$t & eventos$ID == ID & eventos$Event == "saiu" & eventos$Patch == patches[ind,2], ]$Event <- "e"
			}
			
			patchPop[p+1] <- patchPop[p+1] - 1

			if(length(IDlist)>2)
			{
				patches<-patches[-ind,]
				IDlist<-IDlist[-ind]
			}
			else
			{
					patches<-patches[-ind,]
					patches<-matrix(patches, nrow=1, ncol=3)
					IDlist<-IDlist[-ind]
			}
		}

		lin = readLines(dados, n=1)
	}
	out <- strsplit(arquivo, split="/")[[1]][2]
	arqab <- paste0("eventos/abundancia-", out)
	arqtaxas <- paste0("eventos/eventos-", out)

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
