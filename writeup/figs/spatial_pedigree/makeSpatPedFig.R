################################################################
################################################################
#	spatial pedigree conceputal figure
################################################################
################################################################

################################
#	source in useful functions
#	and load required libraries
################################

source("../useful_slim_functions.R")
library(plot3D)


################################
#	run spatial slimulation 
#		with pedigree tracking
################################

call <- "slim 2d_ped.slim"
system(call)


################################
#	make figure
################################

# pick a focal individual
i <- 4
# pick a range of generations to plot
gens <- c(100,97)
# specify the most recent generation simulated
maxGen <- 100

# need to plot everything to get it into the composite 
#	plot list, then can discard this file and
#	replot cleanly with plotdev()
pdf(file="tmp.pdf")
	ancs <- getAncestors(focalInd=i,gen1=gens[1],stopGen=gens[2])
	xlim <- range(ancs[,3]) + c(-0.02,0.02)
	ylim <- range(ancs[,4]) + c(-0.02,0.02)
	plotPed(gens=gens[1]:gens[2],xlim=xlim,ylim=ylim,zlim=NULL,maxGen=maxGen)
	points3D(x=ancs[,3],
			 y=ancs[,4],
			 z=maxGen-ancs[,2],
			 pch=19,
			 cex=1.2,
			 col=adjustcolor(2,0.3),
			 add=TRUE)
	plotPedLines(ancs,maxGen=maxGen,lwd=1)
dev.off()

pdf(file="../spatial_pedigree.pdf",width=8,height=8)
	plotdev()
dev.off()

file.remove("tmp.pdf")
