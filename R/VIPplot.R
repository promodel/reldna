VIPplot<-function#Function to plot VIP-data
###
###
(VIPdata, ##<< matrix with VIP-values
ncomp, ##<< integer, the number of component to plot VIP-values 
col='darkgrey', ##<< color of plot body 
drawlines=TRUE, ##<< logical, if TRUE: lines of cutoffs are drawn
cutoff=seq(0.7, 1.1, 0.1), ##<< numeric vector, specifies cutoffs to be drawn (optional)
cutoffcol=rainbow(length(cutoff)), ##<< character vector, specifies colors for cutoffs (optional) 
numlabels=TRUE, ##<< logical, if TRUE: numbers of variables higher or equal to cutoffs' values are shown 
main='VIP plot', ##<< character string, the title of the plot
ylab=paste('VIP[', ncomp, ']', sep=''), ##<< character string, the y label
xlab='Angstrom', ##<< character string, the x label 
lwd=2, ##<< integer, if drawlines=TRUE: specifies width of lines (optional)
cex=1, ##<< integer, if drawlines=TRUE: specifies symbol size (optional)
save=FALSE, ##<< logical, to save plot into file or not? 
filetype, ##<< if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
filename='VIP', ##<< character string, if save=TRUE specifies a name of the file
show=TRUE, ##<< logical, to show plot or not? 
... ##<< arguments passed to plot
){
  vp<-function(){
    par(mar=c(4.5, 5, 2, 1))
    plot(VIPdata[ , ncomp], type='h', col=col, main=main, xlab=xlab, ylab=ylab, ...)
	if (drawlines==TRUE){  
		abline(h=cutoff, col=cutoffcol, lwd=lwd)
	}
	if (drawlines==TRUE & numlabels==TRUE){
		num<-apply(as.array(cutoff), 1, function(x) as.character(length(which(VIPdata[,ncomp]>=x))))
		edge<-dim(VIPdata)[1]
		text(x=rep(edge-edge/10, length(cutoff)), y=cutoff+0.025, labels=num, col=cutoffcol, cex=cex)
	}
	if (drawlines==FALSE & numlabels==TRUE){ stop('Cannot put numlabels without drawlines=TRUE') }
   }
	if (show==TRUE) { vp() }
	if (save==TRUE) {filesave(vp(), filetype, filename)}
}
