`trecase1irwsh_b0` = function(trc, asn, asnp, thp, Xmatr, est.inits=TRUE,inits=NULL)
{
  nInd = length(trc)
  nAse = length(asn)
  nGenes = 1
  nPar = 3+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trecase1irwsh_b0", as.integer(nInd), as.integer(nAse), as.integer(nGenes), as.integer(thp), as.double(trc),as.double(asn),as.double(asnp),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}

`trecase1irwsh_b1` = function(trc, asn, asnp, thp, Xmatr, est.inits=TRUE,inits=NULL, b0=0)
{
  nInd = length(trc)
  nAse = length(asn)
  nGenes=length(b0)
  nPar = 3+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trecase1irwsh_b1", as.integer(nInd), as.integer(nAse), as.integer(nGenes), as.integer(thp), as.double(trc),as.double(asn),as.double(asnp),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits),as.double(b0), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}

`trecase1irwsh_b01` = function(trc, asn, asnp, thp, Xmatr, est.inits=TRUE,inits=NULL)
{
  nInd = length(trc)
  nAse = length(asn)
  nGenes=1
  nPar = 4+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trecase1irwsh_b01", as.integer(nInd), as.integer(nAse), as.integer(nGenes), as.integer(thp), as.double(trc),as.double(asn),as.double(asnp),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}
