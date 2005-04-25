

imageplot.bma <- function (bicreg.out,color="default", ...) 
{ 
	keep.mar<- par()$mar
	par(mar = c(5,6,4,2) + .1)

	# Plot the which matrix of bicreg as an image. 
	which.out <- bicreg.out$which 
	nvar <- ncol (which.out) 
	nmodel <- nrow (which.out) 
	par (las=1) 
	if (color=="default") 
	{ 
		image (c(0,cumsum(bicreg.out$postprob)),1:nvar,-which.out[1:nmodel,nvar:1, drop=FALSE], 
			xlab="Model #",ylab="",xaxt="n",yaxt="n", xlim=c(0,1),
			main="Models selected by BMA", ...) 
	} 
	if (color=="blackandwhite") 
	{ 
		image (cumsum(bicreg.out$postprob),1:nvar,-which.out[1:nmodel,nvar:1, drop=FALSE], 
			xlab="Model #",ylab="",xaxt="n",yaxt="n",col=c("black","white"), 
			main="Models selected by BMA", ...) 
	} 
	xat<- (cumsum(bicreg.out$postprob) + c(0,cumsum(bicreg.out$postprob[-nmodel])) )/2
	axis (1,at=xat,labels=1:nmodel, ...) 
	axis (2,at=1:nvar, labels=rev(bicreg.out$namesx), ...) 

	par(mar = keep.mar)
} 


