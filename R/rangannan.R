freeG<-function(  
  ### function to calulate the profile of DNA free energy for long sequences
  ### according to A. Kanhere, M. Bansal, BMC Bioinformatics 6, 1 (2005).
  ### =======================================================
  s, ##<< DNA sequence
  bound, ##<< define a fragment of interest
  width=1 ##<< smoothing window width
){
  if (!require(zoo)) {
    stop('Required library "zoo" is not installed.')
  }
  oligos<-getOligos(s)
  l<-length(oligos)
  if(length(bound)>2){
    warning(paste('Length of "bound" is',length(bound),'when 2 is expected. First two values of bound are used.'))
    bound<-bound[1:2]
  }else if(length(bound)==1){
    bound<-c(bound,l-bound)
  }
  data<-rangannan
  oligos.levels<-tolower(data$Step)
  din<-factor(oligos, levels = oligos.levels)
  di<-array(0, dim=c(16, l))
  for(m in 1:l){
    di[din[m],m]<-1
  }
  fenergy<-rollsum(t(data[,2]%*%di),14,align = 'left')
  cfg<-cumsum(fenergy)
l<-length(cfg)
e1<-cfg[100:l]-c(0,cfg[1:(l-100)])
d<-e1[1:(l-250)]-e1[151:(l-100)]
res<-list(dG=fenergy[bound[1]:bound[2]],E100=e1[bound[1]:bound[2]]/100.0,d=d[bound[1]:bound[2]]/100.0  ,seq=s,bound=bound)
class(res)<-'rangannan'
  return(res)
}

