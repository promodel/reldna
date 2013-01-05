lseqcurlib1D <-function(
### main function to calculate electrostatic profile of DNA
s, ##<< DNA sequence
bound, ##<< define a fragment of interest
width=1, ##<< smoothing window width
ref, ##<< reference position 
filename=NA #<< name of the file to save data in. Empty string or NA value means file would not be saved
#,name ##<< name of the library
){
  width<-floor(width/2)
  geom<-dnaGeom(s)
	risem<-geom$risem-geom$risem[ref]

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
	for (i in 1:geom$l){
		q<-qqm[,geom$nseq[i]]
		pot[risef[i]:(risef[i]+zlib-1)]<-pot[risef[i]:(risef[i]+zlib-1)]+(1-rspar[i])*q
		pot[risec[i]:(risec[i]+zlib-1)]<-pot[risec[i]:(risec[i]+zlib-1)]+rspar[i]*q
	}

	rm(qqm, rspar)

	i1<-(zlib/2-width):(max(risef)+zlib/2-width)

	mpot<-pot
	mpot<-mpot[(risef[bound[1]]-8+zlib/2):(risef[bound[2]]+8+zlib/2)]
	i1<-risef[bound[1]:bound[2]]-risef[ref]+9
  zout<-floor(risem[bound[1]]-9):(risem[bound[2]]+9)
  x<-(risem[bound[1]]-8):(risem[bound[2]]+8)
	risem<-risem[bound[1]:bound[2]]-risem[ref]
	if(!is.na(filename)&!nchar(gsub('^ +','',gsub(' +$','',filename)))>0){
    filename<-gsub(' +','_',gsub('^ +','',gsub(' +$','',filename)))
    save(mpot, risef, i1, file=paste(filename,'.lseqcurlib_data.Rdata',sep=''))
  }
  sgz<-data.frame(pos=1:length(s),risem=risem,
                  minZ=(risem-zlib/2),
                  maxZ=(risem+zlib/2-1))
  sgz$minZ[sgz$minZ<min(zout)]<-min(zout)
  sgz$maxZ[sgz$maxZ>max(zout)]<-max(zout)
  sgz$minI<-floor(sgz$minZ)-min(zout)+1
  sgz$maxI<-ceiling(sgz$maxZ)-min(zout)+1  
  
	elstatlist<-list(mpot=mpot, risef=risem, i1=i1,x=x,seq=s,bound=bound,ref=ref,zmap=sgz[sgz$minZ<sgz$maxZ,])
  class(elstatlist)<-'elDNA1d'
  ##value<< list with eight components:
  ##\item{mpot}{1D profile of electrostatic potential along Z axis of DNA;} 
  ##\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
  ##\item{i1}{index of the base pair geometrical center nearest mesh point ;} 
  ##\item{x}{Z coordinates;} 
  ##\item{seq}{DNA sequence used to calculate profile;} 
  ##\item{bound}{boundaries of the part of interest within the sequence;} 
  ##\item{ref}{index of the base pair that suppose to be placed at the origin;} 
  ##\item{zmap}{data frame of geometical properties of base pairs like index, coordinate of the center, part of profile influenced by its charges.} 
  return (elstatlist)
}
