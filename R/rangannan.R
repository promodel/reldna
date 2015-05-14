freeG<-function(  
  ### function to calulate the profile of DNA free energy for long sequences
  ### according to A. Kanhere, M. Bansal, BMC Bioinformatics 6, 1 (2005).
  ### =======================================================
  s, ##<< DNA sequence
  bound, ##<< define a fragment of interest
  width=1 ##<< smoothing window width
){
  oligos<-getOligos(s)
  data<-rangannan
  oligos.levels<-tolower(data$Step)
  din<-factor(oligos, levels = oligos.levels)
  
  di<-array(0, dim=c(16, l-1))
  for(m in 1:l-1){
    di[din[m],m]<-1
  }
  
  fenergy<-t(data[,2]%*%di)
  return(fenergy)
}

