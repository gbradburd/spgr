# reads in the slim spatial pedigree output
#	and converts to a matrix
getSlimGenOutput <- function(outputFile){
	x <- scan(outputFile,what=character(),sep="\n",quiet=TRUE)
	x <- Reduce(rbind,strsplit(x," "))
	colnames(x) <- c("ID","Par1","Par2","GPar1","GPar2","GPar3","GPar4","X","Y")
	return(x)
}

# reads spatial pedigree output for a generation
#	and converts X/Y coordinates to a numeric matrix
getCoords <- function(x){
	coords <- cbind(as.numeric(x[,8]),as.numeric(x[,9]))
	row.names(coords) <- x[,1]
	return(coords)
}

# plots the spatial coordinates for a given generation
#	with z coordinates determined by the difference 
#	between the maximum number of generations run 
#	(aka the present) and the focal generation
plotGen <- function(coords,z,add=FALSE,xlim=NULL,ylim=NULL,zlim=NULL,maxGen,theta=20,phi=-1,r=1){
	coords <- coords[-c(which(coords[,1] < xlim[1] | coords[,1] > xlim[2]),
						which(coords[,2] < ylim[1] | coords[,2] > ylim[2])),]
	z <- maxGen - z
	rect3D(x0 = xlim[1],x1 = xlim[2],
		   y0 = ylim[1],y1 = ylim[2],
		   z0 = z,
		   col=adjustcolor("gray",0.3),
		   bty='n',theta=theta,phi=phi,r=r,
		   colkey=FALSE,xlab="",ylab="",zlab="",
		   xlim=xlim,ylim=ylim,zlim=zlim,
		   add=add)
	scatter3D(x=coords[,1],
			  y=coords[,2],
			  z=rep(z,length(coords[,1])),
			  pch=19,col=adjustcolor("darkgray",0.5),
			  cex=1.2,
			  add=TRUE)
	return(invisible("plotted"))
}

# plots the spatial coordinates for a set of generations
#	with z coordinates determined by the difference 
#	between the maximum number of generations run 
#	(aka the present) and the focal generation
plotPed <- function(gens,xlim=NULL,ylim=NULL,zlim=NULL,maxGen,theta=20,phi=-1,r=1){
	par(mar=c(0,0,0,0))
	if(is.null(xlim)){
		xlim <- c(0,1)
	}
	if(is.null(ylim)){
		ylim <- c(0,1)
	}
	if(is.null(zlim)){
		zlim <- range(maxGen-gens) + c(-0.3,0.3)
	}
	x <- getSlimGenOutput(sprintf("gen_output/gen_%s",gens[1]))
	coords <- getCoords(x)
	plotGen(coords,z=gens[1],add=FALSE,xlim=xlim,ylim=ylim,zlim=zlim,maxGen=maxGen,theta=theta,phi=phi,r=r)
	lapply(gens[2:length(gens)],
		function(g){
			x <- getSlimGenOutput(sprintf("gen_output/gen_%s",g))
			coords <- getCoords(x)
			plotGen(coords,z=g,add=TRUE,xlim=xlim,ylim=ylim,zlim=NULL,maxGen=maxGen)
		})
	text3D(x=mean(xlim)-diff(range(xlim))/6,y=min(ylim)+diff(range(ylim))/9,z=min(zlim),labels=c("longitude"),add=TRUE,cex=1.5)
#	text3D(x=min(xlim)-diff(range(xlim))/3,y=mean(ylim),z=min(zlim)+diff(range(zlim))/7,labels=c("latitude"),add=TRUE,cex=1.5)
	text3D(x=max(xlim)+diff(range(xlim))/4.3,y=min(ylim),z=min(zlim)+diff(range(zlim))/6,labels=c("latitude"),add=TRUE,cex=1.5)
#	text3D(x=max(xlim)+diff(range(xlim))/10,y=min(ylim),z=mean(zlim),labels=c("time"),add=TRUE,cex=1.5)
#	text3D(x=max(xlim)+diff(range(xlim))/20,y=min(ylim),z=min(zlim)+diff(range(zlim))/5,labels=c("present"),add=TRUE,cex=1.5)
#	text3D(x=max(xlim)+diff(range(xlim))/20,y=min(ylim),z=max(zlim)-diff(range(zlim))/5,labels=c("past"),add=TRUE,cex=1.5)
	# arrows3D(x0=max(xlim)+diff(range(xlim))/12,
			 # y0=min(ylim),
			 # z0=min(zlim)+diff(range(zlim))/4,
			 # x1=max(xlim)+diff(range(xlim))/12,
			 # y1=min(ylim),
			 # z1=max(zlim)-diff(range(zlim))/4.25,
			 # add=TRUE)
}


# get the indices of a set of parent IDs 
#	in a set of pedigree IDs
getParentIndices <- function(pedIDs,parentIDs){
	parentIndices <- unlist(
						lapply(parentIDs,
							function(x){
								which(pedIDs==x)}))
	return(parentIndices)
}

# get the space/time locations of ancestors 
#	of a particular focal individual 
#	within a specified list of generations
getAncestors <- function(focalInd,gen1,stopGen){
#	recover()
	nAncs <- sum(2^(0:(gen1-stopGen)))
	ancsInSpace <- matrix(rep(NA,2*4*nAncs),nrow=nAncs,ncol=4)
	gen <- gen1
	x <- getSlimGenOutput(sprintf("gen_output/gen_%s",gen))
	ancsInSpace[1,] <- c(x[focalInd,1],gen,as.numeric(x[focalInd,c(8,9)]))
	parents <- x[focalInd,2:3]
	while(gen > stopGen){
		gen <- gen - 1
		x <- getSlimGenOutput(sprintf("gen_output/gen_%s",gen))
		parIndex <- getParentIndices(x[,1],parents)
		ancIndex <- max(which(!is.na(ancsInSpace),arr.ind=TRUE)[,1]) + 
					1:length(parIndex)
		ancsInSpace[ancIndex,] <- matrix(c(x[parIndex,1],
										   rep(gen,length(parIndex)),
										   x[parIndex,c(8,9)]),ncol=4)
		parents <- unlist(lapply(parIndex,function(i){x[i,2:3]}))
	}
	ancsInSpace <- matrix(as.numeric(ancsInSpace),ncol=4)
	return(ancsInSpace)
}

plotPedLines <- function(ancs,maxGen,lwd,p.cols=1,m.cols=1){
	#recover()
	gens <- maxGen-unique(ancs[,2])
	ancs[,2] <- maxGen-ancs[,2]
	if(m.cols==1 && p.cols==1){
		p.cols <- rep(1,2^(max(gens)+1)-1)
		m.cols <- rep(1,2^(max(gens)+1)-1)
	}
	invisible(
		lapply(max(gens):1,function(g){
			if(length(which(ancs[,2] == g)) > 1){
				if(g > 1){
					nPairs <- length(which(ancs[,2] == g))/2
					lapply(1:nPairs,
						function(i){
							seg.col <- ifelse(i <= nPairs/2,p.cols[g],m.cols[g])
							pair <- ancs[ancs[,2]==g,][(i*2-1):(i*2),]
							child <- ancs[ancs[,2]==(g-1),,drop=FALSE][i,]
							segments3D(x0=pair[1,3],y0=pair[1,4],z0=g,
										x1=pair[2,3],y1=pair[2,4],z1=g,
										lwd=lwd,
										lty=2,add=TRUE,col=seg.col)
							segments3D(x0=mean(pair[,3]),y0=mean(pair[,4]),z0=g,
										x1=child[3],y1=child[4],z1=g-1,
										lwd=lwd,col=seg.col,
										lty=1,add=TRUE)
					})
				} else if(g==1) {
						pair <- ancs[ancs[,2]==g,]
						child <- ancs[ancs[,2]==(g-1),,drop=FALSE]
						segments3D(x0=pair[1,3],y0=pair[1,4],z0=g,
									x1=mean(pair[,3]),
									y1=mean(pair[,4]),z1=g,
									lwd=lwd,
									lty=2,add=TRUE,col=p.cols[g])
						segments3D(x0=pair[2,3],y0=pair[2,4],z0=g,
									x1=mean(pair[,3]),
									y1=mean(pair[,4]),z1=g,
									lwd=lwd,
									lty=2,add=TRUE,col=m.cols[g])
						segments3D(x0=mean(pair[,3])-0.0003,y0=mean(pair[,4])-0.0003,z0=g,
									x1=child[3]-0.0002,y1=child[4]-0.0001,z1=g-1,
									lwd=lwd,col=p.cols[g],
									lty=1,add=TRUE)
						segments3D(x0=mean(pair[,3])+0.0003,y0=mean(pair[,4])+0.0003,z0=g,
									x1=child[3]+0.0002,y1=child[4]+0.0002,z1=g-1,
									lwd=lwd,col=m.cols[g],
									lty=1,add=TRUE)
				}
			}
		})
	)
}

color.func <- function(x,cols){
	cols <- colorRampPalette(cols)(length(x))[findInterval(x,sort(x))]
	return(cols)
}

plotNonSpPed <- function(ancs,mar=c(3,0,3,0),anc.cols=NULL,segLwd=1.5,pt.cex=2){
	#recover()
	nGen <- length(unique(ancs[,2]))
	if(is.null(anc.cols)){
		p.cols <- color.func(1:nGen,c("darkblue","blue","deepskyblue"))
		m.cols <- color.func(1:nGen,c("darkred","red","coral1"))
	}
	par(mar=mar)
	plot(0,type='n',xlim=c(0,1.4),ylim=c(1,nGen+0.2),xlab="",ylab="",main="",xaxt='n',yaxt='n',bty='n')
	for(n in nGen:2){
		if(n == nGen){
			x <- cumsum(rep((1/(2^(nGen-1)+1)),2^(nGen-1)))
			y <- rep(nGen,2^(nGen-1))
		} else {
			x <- (x[seq(1,2^n,by=2)] + x[seq(2,2^n,by=2)])/2
			y <- rep(n,2^(n-1))
		}
		points(x=x,y=y,pch=19,col=c(rep(p.cols[n],(2^n)/4),rep(m.cols[n],(2^n)/4)),cex=pt.cex)
		if(n > 2){
			nPairs <- 2^(n-1)/2
			lapply(1:nPairs,
				function(i){
					seg.col <- ifelse(i <= nPairs/2,p.cols[n],m.cols[n])
					segments(x0=x[i*2-1],x1=x[i*2],y0=n,y1=n,lty=2,col=seg.col,lwd=segLwd)
					segments(x0=mean(c(x[i*2-1],x[i*2])),
							 x1=mean(c(x[i*2-1],x[i*2])),
							 y0=n,
							 y1=n-1,lty=1,
							 col=seg.col,lwd=segLwd)
				})
		}
		if(n == 2){
			segments(x0=x[1],x1=mean(x),y0=n,y1=n,lty=2,col=p.cols[n],lwd=segLwd)
			segments(x0=x[2],x1=mean(x),y0=n,y1=n,lty=2,col=m.cols[n],lwd=segLwd)
			segments(x0=mean(x)-0.002,
					 x1=mean(x)-0.002,
					 y0=n,
					 y1=n-1,lty=1,
					 col=p.cols[n],lwd=segLwd)
			segments(x0=mean(x)+0.002,
					 x1=mean(x)+0.002,
					 y0=n,
					 y1=n-1,lty=1,
					 col=m.cols[n],lwd=segLwd)
		}
	}
	points(x=0.5,y=1,col="darkorchid4",cex=pt.cex,pch=19)
	text(x=1.2,y=1,labels="now",cex=2)
	text(x=1.2,y=2,labels="parents",cex=2)
	text(x=1.2,y=3,labels="grandparents",cex=2)
	text(x=1.2,y=4,labels="great\ngrandparents",cex=2)
	return(invisible("plotted"))
}