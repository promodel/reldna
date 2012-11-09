VIP_cutoff<-function#Function to determine VIP-values higher or equal to cutoffs
###
###
(VIPdata, ##<< matrix with VIP-values 
ncomp, ##<< integer, the number of component to plot VIP-values   
cutoffs, ##<< numeric vector, specifies cutoffs 
save=FALSE, ##<< logical, to save results into Rdata-file or not? 
filename=paste('VIP_cutoff_', paste(cutoffs, collapse='_'), '.Rdata', sep='')##<< character string, if save=TRUE specifies a name of the file
){
	if(length(cutoff)==1){
		result<-names(VIPdata[ , ncomp])[which(VIPdata[ , ncomp]>=cutoffs)]
	}
	if(length(cutoff)>1){
		result<-lapply(cutoffs, function(x) names(VIPdata[ , ncomp])[which(VIPdata[ , ncomp]>=x)])
	}
	if (save==TRUE){ 
		save(result, file=filename) 
	}
	return(invisible(result))
} 
