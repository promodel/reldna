dnaGeom<-function(
  ### Function to calculate geometry of the DNA double helix.   
  s##<< sequence of the DNA span to be build
  ){
  if(!require(seqinr)){
    stop('Required library "seqinr" is not installed.')
  }
  s<-tolower(gsub('\\s', '', s,perl=TRUE))
  l<-getLength(s)
  nseq<-s2n(s2c(s), levels=c("a", "t", "g", "c"), base4=FALSE)
  istarts <- seq(from = 1, to = l-1, by = 1)
  oligos <- s2c(s)[istarts]
  oligos <- paste(oligos, s2c(s)[istarts + 1], sep = "")

  #  load('rise_twist.Rdata')
  #	data(rise_twist,qqs,w)
  data<-rise_twist
  oligos.levels<-tolower(data$dilet)
  din<-factor(oligos, levels = oligos.levels)

  di<-array(0, dim=c(16, l-1))
#   for (k in 1:16) {
#     i<-gregexpr(data$dilet[k], s, ignore.case=TRUE)
#     i<-as.vector(i[[1]])
#     di[k,i]<-array(1, dim=c(1, length(i)))
#   }
#   for (m in 1:(l-1)){ 
#     if(sum(di[, m])!=1){ 
#       di[, m]<-di[,m-1]
#     }}
  for(m in 1:l-1){
    di[din[m],m]<-1
  }

  geom<-array(0, dim=c(2,l))
  phi<-cumsum(c(0,data$twist%*%di/180*pi))
  risem<-cumsum(c(0,data$rise%*%di))
  
  rm(data, geom, di)
  geom<-list(risem=risem,twist=phi,nseq=nseq,l=l)
  class(geom)<-'DNAgeom'
  ##value<< object of class 'DNAgeom', list with four slots:
  ##\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
  ##\item{twist}{that gives the orientation of the X axis of the base pair;} 
  ##\item{nseq}{simple numerical encoding of a DNA sequence by \code{seqinr::s2n} ;} 
  ##\item{l}{length of the sequence.} 
  return(geom)
  ##references<< http://chem.rutgers.edu/~xiangjun/3DNA/images/bp_step_hel.gif
}
