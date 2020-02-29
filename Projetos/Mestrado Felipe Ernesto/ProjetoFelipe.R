library(doMC)
registerDoMC(8)

# Função "wrapper" para a chamada em C:
#
# Os passos abaixo foram adaptados de http://users.stat.umn.edu/~geyer/rc/

Sys.setenv("PKG_CPPFLAGS" = "-fopenmp -DPARALLEL") # liga biblioteca de paralelismo
system("rm ../../TWoLife.so") #limpa sources velhos
system("rm ../../TWoLife.o") #limpa sources velhos
system ("R CMD SHLIB ../../TWoLife.cpp") ## compila no R
dyn.load("../../TWoLife.so") ## carrega os source resultantes como biblioteca dinamica no R

# Generates the landscape with specified conditions.
# numb.cells represents both the lenght AND width of the landscape, so numb.cells=100 creates a 100x100 landscape
# Land.shape can be 0 = XXX or 1 = XXX.
# Bound.condition can be 0 = XXX or 1 = XXX.
readLandscape <- function (h, f, r, np)
{
	land <- rep(-1, np*np)
	nome <- paste("../Paisagens/", f, "_", h, "_", r, ".txt", sep="")
	paisagem <- file(nome, "r")
	for(i in seq(0, np-1))
	{
		lin = readLines(paisagem, n=1)
		lin <- strsplit(lin, " ")
		line <- unlist(lin)
		for(j in seq(0, np-1))
		{
			land[1+np*j+i] <- strtoi(line[j+1])
		}
	}
	close(paisagem)
	return(land)
}

TWoLife <- function (
					 raio=0.1,
					 N=80,
					 AngVis=360,
					 passo=5,
					 taxa.move=0.5,
					 taxa.basal=0.6,
					 taxa.morte=0.1,
					 incl.birth=0.5/0.01,
					 incl.death=0,
					 density.type=0,
					 death.mat=7,
					 move.mat=7,
					 landscape,
					 tempo=20,
           ini.config=0,
           out.code="1")
{
	if(class(landscape) != "landscape") {
		stop("Error in function TWoLife: you must provide a valid landscape. See ?Landscape")
	}
  if(raio>landscape$numb.cells*landscape$cell.size/2)
  {stop("Error in function TWoLife: the radius must be lower than or equal to the half of landscape side (radius <= numb.cells*cell.size/2)")}

  saida.C <- .C("TWoLife",
              as.double(raio),# 1
              as.integer(N),# 2
              as.double(AngVis),# 3
              as.double(passo),# 4
              as.double(taxa.move),# 5
              as.double(taxa.basal),# 6
              as.double(taxa.morte),# 7
              as.double(incl.birth),# 8
              as.double(incl.death),# 9
              as.integer(landscape$numb.cells),# 10
              as.double(landscape$cell.size),# 11
              as.integer(landscape$land.shape),# 12
              as.integer(density.type),# 13
              as.double(death.mat), # 14
							as.double(move.mat), # 14
              as.integer(ini.config), #15
              as.integer(landscape$bound.condition), #16
              as.integer(landscape$scape), #17
              as.double(tempo), #18
              as.integer(0), # 19
              as.double(rep(0, 5000)), # 20
              as.double(rep(0,5000)), # 21
              as.character(out.code)
              ## verificar se precisa definir o tamanho e se isto nao dará problemas (dois ultimos argumentos)
				  )
	n <- saida.C[[19]]
	x <- saida.C[[20]]
	y <- saida.C[[21]]
	x <- x[1:n]; y <- y[1:n]
	return(data.frame(x=x,y=y))
}

projetoFelipe<-function()
{
#Diretório que receberá os outputs para esta combinação de parâmetros
#Sempre que for começar a rodar lembrar de colocar "run 1" em iteration.txt
ite <- file("iteration.txt", "r")
linha <- readLines(ite, n=1)
nite <- strtoi(unlist(strsplit(linha, " "))[2])
dir <- paste("LHS-", nite, sep = "")
dir.create(dir)
cat(paste("run ", nite+1, sep = ""), file = "iteration.txt", append = FALSE, sep = "\n")
close(ite)

#Leitura dos parâmetros de "Hipercubo.csv"
LHS <- read.csv("Hipercubo.csv")

setwd(dir)

R <- LHS$R[nite]
s <- LHS$s[nite]
delta <- LHS$delta[nite]
lamb0 <- LHS$lamb0[nite]
mu0 <- LHS$mu0[nite]
b <- LHS$b[nite]
d <- LHS$d[nite]
angvis <- LHS$angvis[nite]
dm <- LHS$dm[nite]
mm <- LHS$mm[nite]

#Arquivo metadata com os valores dos parâmetros
file.create("METADATA.txt")
cat(date(), file = "METADATA.txt", append = TRUE, sep = "\n")
cat("", file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("R:", R), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("s:", s), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("delta:", delta), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("lamb0:", lamb0), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("mu0:", mu0), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("b:", b), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("d:", d), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("angvis:", angvis), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("dm:", dm), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("mm:", mm), file = "METADATA.txt", append = TRUE, sep = "\n")


HabProp <- seq(0, 100, 5)
HabFrag <- c(0.8)
Rep <- data.frame(Conf = rep(1:10, each=5), Rep = rep(1:5, 10))
for(i in HabProp)
{
	for(j in HabFrag)
	{
		foreach(k = iter(Rep, by="row")) %dopar%
		{

				h = formatC(i, width = 3, format = "d", flag = "0")
				f = format(j, nsmall=2)
				c = formatC(k$Conf, width = 5, format = "d", flag = "0")
				scape <- readLandscape(h, f, c, 128)
				land <- list(numb.cells = 128, cell.size = 0.03, bound.condition = 0, land.shape = 1, scape = scape )
				class(land) <- "landscape"
				c <- formatC(k$Conf, width = 2, format = "d", flag = "0")
				o.c <- paste("output-",f,"_",h, "_", c,"_", k$Rep, ".txt", sep="")
				TWoLife(raio=R, N=300, AngVis=angvis, passo=s, taxa.move=delta, taxa.basal=lamb0, taxa.morte=mu0, incl.birth=b, incl.death=d, density.type=1, death.mat=dm, move.mat=mm, landscape = land, tempo=50, ini.config=1, out.code=o.c)

		}
	}
}

setwd("..")

}

projetoFelipe()