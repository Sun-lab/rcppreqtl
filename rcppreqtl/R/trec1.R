`trec1irw0` = function(trc, thp, Xmatr, est.inits=TRUE,inits=NULL)
{
  nInd = length(trc)
  nGenes=1
  nPar = 1+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trec1irw0", as.integer(nInd), as.integer(nGenes), as.integer(thp),as.double(trc),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}

`trec1irw_b0` = function(trc, thp, Xmatr, est.inits=TRUE,inits=NULL)
{
  nInd = length(trc)
  nGenes=1
  nPar = 2+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trec1irw_b0", as.integer(nInd), as.integer(nGenes), as.integer(thp),as.double(trc),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}

`trec1irw_b1` = function(trc, thp, Xmatr, est.inits=TRUE,inits=NULL, b0=1)
{
  nInd = length(trc)
  nGenes=length(b0)
  nPar = 2+ncol(Xmatr)
  
  if(is.null(inits))inits = rep(0,nPar+1)  
  if(length(inits)!=nPar){
    inits = rep(0,nPar+1)
  }else{
    inits = c(inits,0)
  }
      
  Z = .C("max_trec1irw_b1", as.integer(nInd), as.integer(nGenes), as.integer(thp),as.double(trc),
  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits),as.double(b0), PACKAGE="rcppreqtl")
  #cat("returned: ",Z$optval,"\n")
  return(Z$optval)
}
