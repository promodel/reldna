dnaGeom<-function(
  ### Function to calculate geometry of the DNA double helix.   
  s##<< sequence of the DNA span to be build
  ){
  if(!require(seqinr)){
    stop('Required library "seqinr" is not installed.')
  }
  s<-tolower(gsub(' ', '', s))
  l<-getLength(s)
  nseq<-s2n(s2c(s), levels=c("a", "t", "g", "c"), base4=FALSE)
  
  #  load('rise_twist.Rdata')
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
  phi<-cumsum(data$twist%*%di/180*pi)
  risem<-cumsum(data$rise%*%di)
  
  rm(data, geom, di)
  geom<-list(risem=risem,twist=phi,nseq=nseq,l=l)
  class(geom)<-'DNAgeom'
  return(geom)
  ### object of class 'DNAgeom' with two slots: rise, which gives Z coordinate of the corresponding base pair, and 
  ### twist, which gives orientation of the X axis of the base pair.
  ##references<< http://chem.rutgers.edu/~xiangjun/3DNA/images/bp_step_hel.gif
}
