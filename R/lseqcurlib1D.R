lseqcurlib1D <-function(
### main function to calculate electrostatic profile of DNA
s, ##<< DNA sequence
bound, ##<< define a fragment of interest
width=1, ##<< smoothing window width
ref, ##<< reference position 
filename=NA #<< name of the file to save data in. Empty string or NA value means file would not be saved
#,name ##<< name of the library
){
	library(seqinr)
	width<-floor(width/2)
	s<-tolower(gsub(' ', '', s))
	l<-getLength(s)
	nseq<-s2n(s2c(s), levels=c("a", "t", "g", "c"), base4=FALSE)

#	load('rise_twist.Rdata')
#	data(rise_twist,qqs,w)
	data<-rise_twist

	di<-array(0, dim=c(16, l-1))
	for (k in 1:16) {
		i<-gregexpr(data$dilet[k], s, ignore.case=TRUE)
		i<-as.vector(i[[1]])
		di[k,i]<-array(1, dim=c(1, length(i)))
	}
	for (m in 1:(l-1)){ 
		if(sum(di[, m])!=1){ 
			di[, m]<-di[,m-1]
	}}

	geom<-array(0, dim=c(1,l))
	geom[1, 2:l]<-data$rise%*%di  

	risem<-cumsum(geom)

	rm(data, geom, di,  s)

	risem<-risem-risem[ref]

	risef<-floor(risem)
	i<-which(abs(risem-risef-1)<1e-10)
	risef[i]<-risef[i]+1

	risec<-ceiling(risem)
	rspar<-risem-risef
	risef<-floor(risem-min(risem))+1
	risec<-ceiling(risem-min(risem))+1
#	libname<-paste(name, '.Rdata', sep='')
#	data(libname)
	zlib<-dim(qqs)[2]
	lz<-max(risec)-min(risec)+zlib
	pot<-array(0, dim=c(1,lz))
	qqm<-apply(qqs, c(2,3), FUN=mean)
	for (i in 1:l){
		q<-qqm[,nseq[i]]
		pot[risef[i]:(risef[i]+zlib-1)]<-pot[risef[i]:(risef[i]+zlib-1)]+(1-rspar[i])*q
		pot[risec[i]:(risec[i]+zlib-1)]<-pot[risec[i]:(risec[i]+zlib-1)]+rspar[i]*q
	}

	rm(qqm, rspar)

	i1<-(zlib/2-width):(max(risef)+zlib/2-width)

	mpot<-pot
	mpot<-mpot[(risef[bound[1]]-8+zlib/2):(risef[bound[2]]+8+zlib/2)]
	risem<-risem[bound[1]:bound[2]]-risem[ref]
	i1<-risef[bound[1]:bound[2]]-risef[bound[1]]+9
  
  if(!is.na(filename)&!nchar(gsub('^ +','',gsub(' +$','',filename)))>0){
    filename<-gsub(' +','_',gsub('^ +','',gsub(' +$','',filename)))
    save(mpot, risef, i1, file=paste(filename,'.lseqcurlib_data.Rdata',sep=''))
  }

	elstatlist<-list(mpot=mpot, risef=risef, i1=i1)
	return (elstatlist)
### list containing a few electrostatic data: pot, mpot, risef, i1
}
