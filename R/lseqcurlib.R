lseqcurlib <-function(
### main function to calculate electrostatic profile of DNA
s, ##<< DNA sequence
bound, ##<< define a fragment of interest
width, ##<< smoothing window width
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

	geom<-array(0, dim=c(2,l))
	geom[1,2:l]<-data$rise%*%di  
	geom[2,2:l]<-data$twist%*%di          
	twistm<-cumsum(geom[2,])
	risem<-cumsum(geom[1,])

	rm(data, geom, di,  s)

	twistm<-twistm-twistm[ref]
	risem<-risem-risem[ref]

	twistf<-floor(twistm)
	i<-which(abs(twistm-twistf-1)<1e-10)
	twistf[i]<-twistf[i]+1
	risef<-floor(risem)
	i<-which(abs(risem-risef-1)<1e-10)
	risef[i]<-risef[i]+1

	twistc<-ceiling(twistm)
	twpar<-twistm-twistf
	risec<-ceiling(risem)
	rspar<-risem-risef
	twistf<-(twistf%%360)+1
	twistc<-(twistc%%360)+1
	i<-which(twistf==360)
	twistf[i]<-359
	i<-which(twistc==360)
	twistc[i]<-359
	risef<-floor(risem-min(risem))+1
	risec<-ceiling(risem-min(risem))+1
#	libname<-paste(name, '.Rdata', sep='')
#	load(libname)
	zlib<-dim(qqs)[2]
	lz<-max(risec)-min(risec)+zlib
	pot<-array(0, dim=c(360,lz))
	for (i in 1:l){
		q<-qqs[,,nseq[i]]
		par<-twpar[i]
		f<-twistf[i]
		c<-twistc[i]
		qq<-(1-par)*rbind(q[(361-f):360,], q[1:(360-f),])+par*rbind(q[(361-c):360,], q[1:(360-c),])
		pot[,risef[i]:(risef[i]+zlib-1)]<-pot[,risef[i]:(risef[i]+zlib-1)]+(1-rspar[i])*qq
		pot[,risec[i]:(risec[i]+zlib-1)]<-pot[,risec[i]:(risec[i]+zlib-1)]+rspar[i]*qq
	}

	rm(twistm, twistf, twpar, qqs, rspar, qq, par, c, f)

	i1<-(zlib/2-width):(max(risef)+zlib/2-width)
	i2<-(zlib/2+width):(max(risef)+zlib/2+width)
	i1<-(i1-1)*360+1
	i2<-i2*360
	mpot<-array(0, dim=c(1, length(i1)))
	for (i in 1:length(i1)){
   		mpot[i]<-mean(pot[i1[i]:i2[i]])
	}

	pot<-pot[,(risef[bound[1]]-8+zlib/2):(risef[bound[2]]+8+zlib/2)]
	mpot<-mpot[(risef[bound[1]]-8+zlib/2):(risef[bound[2]]+8+zlib/2)]
	risem<-risef[bound[1]:bound[2]]-risef[ref]
	i1<-risef[bound[1]:bound[2]]-risef[bound[1]]+9

	if(!is.na(filename)&!nchar(gsub('^ +','',gsub(' +$','',filename)))>0){
	  filename<-gsub(' +','_',gsub('^ +','',gsub(' +$','',filename)))
	  save(pot, mpot, risef, i1, file=paste(filename,'.lseqcurlib_data.Rdata',sep=''))
	}
	elstatlist<-list(pot=pot, mpot=mpot, risef=risef, i1=i1)
	return (elstatlist)
### list containing a few electrostatic data: pot, mpot, risef, i1
}
