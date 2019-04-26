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

if(!file.exists("gen_output/gen_100")){
	call <- "slim 2d_ped.slim"
	system(call)
}


################################
#	make figure
################################

# pick a range of generations to plot
gens <- c(100,97)

# pick a focal individual
# goodOnes <- sapply(1:100,
				# function(i){
					# ancs <- getAncestors(focalInd=i,gen1=gens[1],stopGen=gens[2])
					# length(unique(ancs[which(ancs[,2]==gens[2]),1])) == 2^(gens[1]-gens[2])
				# })

i <- 5	#sample(which(goodOnes),1)

# specify the most recent generation simulated
maxGen <- 100

# need to plot everything to get it into the composite 
#	plot list, then can discard this file and
#	replot cleanly with plotdev()


p.cols <- color.func((gens[2]+1):gens[1],c("darkblue","blue","deepskyblue"))
m.cols <- color.func((gens[2]+1):gens[1],c("darkred","red","coral1"))

getAncCols <- function(ancs,p.cols,m.cols){
	ancGens <- max(ancs[,2]) - ancs[,2]
	anc.cols <- numeric(nrow(ancs))
	anc.cols[1] <- "darkorchid4"
	for(i in 1:max(ancGens)){
		anc.cols[2^i:(2^(i+1)-1)] <- c(rep(p.cols[i],2^(i-1)),
									   rep(m.cols[i],2^(i-1)))
	}
	return(anc.cols)
}

pdf(file="tmp.pdf")
	ancs <- getAncestors(focalInd=i,gen1=gens[1],stopGen=gens[2])
	anc.cols <- getAncCols(ancs,p.cols,m.cols)
	buffer <- 0.02 * c(-1,1)
	xlim <- range(ancs[,3]) + buffer
	ylim <- range(ancs[,4]) + buffer
	plotPed(gens=gens[1]:gens[2],xlim=xlim,ylim=ylim,zlim=NULL,maxGen=maxGen,theta=20,phi=9,r=1000)
	plotPedLines(ancs,maxGen=maxGen,lwd=3,p.cols=p.cols,m.cols=m.cols)
	points3D(x=ancs[,3],
			 y=ancs[,4],
			 z=maxGen-ancs[,2],
			 pch=19,
			 cex=2.4,
			 col=anc.cols,
			 colkey=list(plot=FALSE),
			 colvar=NULL,
			 add=TRUE)
dev.off()


pdf(file="../spatial_pedigree.pdf",width=16,height=8)
	par(mfrow=c(1,2))
	plotNonSpPed(ancs,mar=c(4,0,4,0),segLwd=3,pt.cex=2.5)
	plotdev()
dev.off()


# fix time slices perspective
# fix time slice labeling
# fix color discrepancy btwn sp and nsp
# fix time slice color and other inds color
#	fix associated legend

file.remove("tmp.pdf")