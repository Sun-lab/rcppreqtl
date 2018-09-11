fit = function(subset=NULL, data, traceit=FALSE){#ensure all the input parameters are there
#data = dat; subset=1:100; traceit=FALSE
  #list(haplotype=haplotypef, haplotype4=haplotype, 
  #  params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt), 
  #  X=Xmatr, trc=trcm, asn=asnf, asnp=asnpf, asn4=asnm, asnp4=asnpm)

  Xmatr = data$X
  nbeta = ncol(data$X)
  llfsh = nbeta + 5  
  ll0sh = nbeta + 4  
  ll0tr = nbeta + 3
  llsh = nbeta + 2
  if(is.null(subset))subset=1:nrow(data$trc)
  nsubit = length(subset)
  trec_shC = matrix(NA,nrow=nsubit,ncol=llsh)
  trec_b0C = trec_b1C = matrix(NA,nrow=nsubit,ncol=ll0tr)
  trecase_b1C = matrix(NA,nrow=nsubit,ncol=ll0sh)
  trecase_b01C = matrix(NA,nrow=nsubit,ncol=llfsh)
#  if(traceit)message("trace 1")
  rownames(trec_shC) = rownames(trec_b0C) = rownames(trecase_b1C) = rownames(trecase_b01C) = rownames(data$trc)[subset]
  trecase_b1C0 = trecase_b0C = trecase_b1C
#  if(traceit)message("trace 1: ", length(subset))

  trcm = data$trc
  asnm = data$asn
  asnpm = data$asnp
  thp = data$haplotype
#  if(traceit)message("trace 2")

  for(j in 1:nsubit){   
  #j = 1
    i = subset[j]
    inits = c(-2,rep(0,nbeta))
    trc = trcm[i,]
    asn = asnm[i,]
    asnp = asnpm[i,]
#    if(traceit)message("trace 2a", length(inits))

    initsirw0 = trec1irw0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits,est.inits=T)
#    message(paste(initsirw0, collapse=" "))
#    if(traceit)message("trace 2b", length(initsirw0))
    if(any(!is.finite(initsirw0))){
      initsirw0 = trec1irw0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits,est.inits=F)
    }  
#    message(paste(initsirw0, collapse=" "))
#    if(traceit)message("trace 2c ", length(initsirw0), " ", ncol(trec_shC), " ", nrow(trec_shC), " ", i, " ", j)
    trec_shC[i,]=initsirw0
#  if(traceit)message("trace 2d")
    b0i = log(median(trc[thp==3]))-log(median(trc[thp==0]))
    b0i2 = log((median(trc[thp==1])+median(trc[thp==2]))/2)-log(median(trc[thp==0]))
    if(!is.finite(b0i)){
      b0i = b0i2
    }
    if(!is.finite((b0i+b0i2)/2)){
      b0i = 0
    }
    inits0 = c(initsirw0[1],b0i,initsirw0[-c(1,length(initsirw0))])
    initsirw = trec1irw_b0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits0,est.inits=F)
    trec_b0C[i,] = initsirw
    
    b0i = initsirw[2];inits1 = initsirw[-(nbeta+3)];inits1[2] = 0
    b1i = sign(log(mean(trc[thp==1]))-log(mean(trc[thp==2])))
    if(!is.finite(b1i)){
      b1i = 0
    }
#  if(traceit)message("trace 2e")
  
    inits1[2] = b1i
    initsirw = trec1irw_b1(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits1,est.inits=F,b0=b0i)
    trec_b1C[i,] = initsirw
    trec_b1C[i,]
#  if(traceit)message("trace 3")
    
    #short
    #trecase  
    #additive
    ini_b0ta = c(trec_b0C[i,1],trec_b0C[i,1],trec_b0C[i,-c(1,ll0tr)])
    if(abs(ini_b0ta[3])>8)ini_b0ta[3]=sign(ini_b0ta[3])*8
    if(abs(ini_b0ta[4])>8)ini_b0ta[4]=sign(ini_b0ta[4])*8
    est_b0ta = trecase1irw_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
    trecase_b0C[i,] = est_b0ta#fit b0 assuming no b1
    trecase_b0C[i,]
  
    #poo est 1 (pre joined)
    fit_b0ta = trecase_b0C[i,3]
    ini_b1ta = trecase_b0C[i,-8];ini_b1ta[3]=0
    ini_b1ta[2] = -ini_b1ta[1]
    est_b1ta = trecase1irw_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=fit_b0ta)
    trecase_b1C[i,] = est_b1ta
    trecase_b1C[i,]
    
    #joined
    ini_b01ta = c(trecase_b1C[i,1:2],fit_b0ta,trecase_b1C[i,3:(ll0sh-1)])
    est_b01ta = trecase1irw_b01(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b01ta,est.inits=F)
    trecase_b01C[i,] = est_b01ta
    trecase_b01C[i,]
  
    #poo, at b0=0
    ini_b1ta = trecase_b01C[i,-c(3,llfsh)];#ini_b1ta[3]=0
    est_b1ta = trecase1irw_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=0)
    trecase_b1C0[i,] = est_b1ta#fit b1 assuming no b0
    trecase_b1C0[i,]
  
    #additive redo, at b1=0
    ini_b0ta = trecase_b01C[i,-c(4,llfsh)];#ini_b0ta[3]=0
    if(abs(ini_b0ta[3])>8)ini_b0ta[3]=sign(ini_b0ta[3])*8
    if(abs(ini_b0ta[4])>8)ini_b0ta[4]=sign(ini_b0ta[4])*8
    est_b0ta = trecase1irw_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
    delt0 = trecase_b0C[i,ll0sh]-est_b0ta[ll0sh]
    if((delt0> 1e-4) & traceit){
      message(i, " upd H0(b0=0)comm = ", format(delt0, scientific=T, digits=2))
      trecase_b0C[i,] = est_b0ta    
    }
#  if(traceit)message("trace 4")
  
    #f
    #joined
    ini_b01ta = trecase_b01C[i,-(nbeta+5)]
    trecase_b0C[i,nbeta+4]
    trecase_b1C0[i,nbeta+4]
    trecase_b01C[i, nbeta+5]
    if(trecase_b01C[i, llfsh]>trecase_b1C0[i,ll0sh]){
      ini_b01ta[1:2] = trecase_b1C0[i, 1:2]
      ini_b01ta[3] = 0
      ini_b01ta[4:(llfsh-1)] = trecase_b1C0[i, 3:(ll0sh-1)]      
    }
    if(trecase_b01C[i, llfsh]>trecase_b0C[i,ll0sh]){
      ini_b01ta[1:3] = trecase_b1C0[i, 1:3]
      ini_b01ta[4] = 0
      ini_b01ta[5:(llfsh-1)] = trecase_b1C0[i, 4:(ll0sh-1)]      
    }
    est_b01taf = trecase1irw_b01(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b01ta,est.inits=F)
    rerun = F
    if(est_b01taf[llfsh] < (trecase_b01C[i,llfsh])){
      rerun = T
      trecase_b01C[i,] = est_b01taf
    }
#  if(traceit)message("trace 5")

    if(rerun){  
      #poo, at b0=0
      ini_b1ta = trecase_b01C[i,-c(3,llfsh)]
      est_b1ta = trecase1irw_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=0)
      if(est_b1ta[ll0sh] < (trecase_b1C0[i,ll0sh])){
        trecase_b1C0[i,] = est_b1ta#fit b1 assuming no b0
      }  
      #additive
      ini_b0ta = trecase_b01C[i,-c(4,llfsh)]
      est_b0ta = trecase1irw_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
      if(est_b0ta[ll0sh] < (trecase_b0C[i,ll0sh])){
        trecase_b0C[i,] = est_b0ta#fit b1 assuming no b0
      }  
    }
    if(i%%10==0)cat(i,"'th iteration\n")
  }
#  if(traceit)message("trace 6")
  apply(trecase_b01C,2,median)
  
  pval_b0 = 2*(trecase_b0C[, ll0sh]-trecase_b01C[,llfsh])
  pval_b1 = 2*(trecase_b1C0[, ll0sh]-trecase_b01C[,llfsh])
  if(traceit){
    message("clear fails for full likelihood: ", sum(pval_b0< -.1 | pval_b1< -.1))
  }
  fixb0 = pval_b0<0
  fixb1 = pval_b1<0
  pval_b0[fixb0] = 0
  pval_b1[fixb1] = 0
  pval_b0 = pchisq(pval_b0, df=1, lower.tail=F)
  pval_b1 = pchisq(pval_b1, df=1, lower.tail=F)
  
  list(full=cbind(trecase_b01C, pval_b0, pval_b1), testadd=trecase_b0C, testpoo=trecase_b1C0)
}


#common od assumption
fitsh = function(subset=NULL, data, traceit=FALSE){#ensure all the input parameters are there
#data = dat; subset=1:100; traceit=FALSE
  #list(haplotype=haplotypef, haplotype4=haplotype, 
  #  params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt), 
  #  X=Xmatr, trc=trcm, asn=asnf, asnp=asnpf, asn4=asnm, asnp4=asnpm)

  Xmatr = data$X
  nbeta = ncol(data$X)
  llfsh = nbeta + 5  
  ll0sh = nbeta + 4  
  ll0tr = nbeta + 3
  llsh = nbeta + 2
  if(is.null(subset))subset=1:nrow(data$trc)
  nsubit = length(subset)
  trec_shC = matrix(NA,nrow=nsubit,ncol=llsh)
  trec_b0C = trec_b1C = matrix(NA,nrow=nsubit,ncol=ll0tr)
  trecase_b1C = matrix(NA,nrow=nsubit,ncol=ll0sh)
  trecase_b01C = matrix(NA,nrow=nsubit,ncol=llfsh)
  rownames(trec_shC) = rownames(trec_b0C) = rownames(trecase_b1C) = rownames(trecase_b01C) = rownames(data$trc)[subset]
  trecase_b1C0 = trecase_b0C = trecase_b1C
  
  trcm = data$trc
  asnm = data$asn
  asnpm = data$asnp
  thp = data$haplotype

  for(j in 1:nsubit){   
  #j = 1
    i = subset[j]
    inits = c(-2,rep(0,nbeta))
    trc = trcm[i,]
    asn = asnm[i,]
    asnp = asnpm[i,]
    
    initsirw0 = trec1irw0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits,est.inits=T)
    if(any(!is.finite(initsirw0))){
      initsirw0 = trec1irw0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits,est.inits=F)
    }  
    trec_shC[i,]=initsirw0
    b0i = log(median(trc[thp==3]))-log(median(trc[thp==0]))
    b0i2 = log((median(trc[thp==1])+median(trc[thp==2]))/2)-log(median(trc[thp==0]))
    if(!is.finite(b0i)){
      b0i = b0i2
    }
    if(!is.finite((b0i+b0i2)/2)){
      b0i = 0
    }
    inits0 = c(initsirw0[1],b0i,initsirw0[-c(1,length(initsirw0))])
    initsirw = trec1irw_b0(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits0,est.inits=F)
    trec_b0C[i,] = initsirw
    
    b0i = initsirw[2];inits1 = initsirw[-(nbeta+3)];inits1[2] = 0
    b1i = sign(log(mean(trc[thp==1]))-log(mean(trc[thp==2])))
    if(!is.finite(b1i)){
      b1i = 0
    }
  
    inits1[2] = b1i
    initsirw = trec1irw_b1(trc=trc,thp=thp,Xmatr=Xmatr,inits=inits1,est.inits=F,b0=b0i)
    trec_b1C[i,] = initsirw
    trec_b1C[i,]
    
    #short
    #trecase  
    #additive
    ini_b0ta = c(trec_b0C[i,1],trec_b0C[i,1],trec_b0C[i,-c(1,ll0tr)])
    if(abs(ini_b0ta[3])>8)ini_b0ta[3]=sign(ini_b0ta[3])*8
    if(abs(ini_b0ta[4])>8)ini_b0ta[4]=sign(ini_b0ta[4])*8
    est_b0ta = trecase1irwsh_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
    trecase_b0C[i,] = est_b0ta#fit b0 assuming no b1
    trecase_b0C[i,]
  
    #poo est 1 (pre joined)
    fit_b0ta = trecase_b0C[i,3]
    ini_b1ta = trecase_b0C[i,-8];ini_b1ta[3]=0
    ini_b1ta[2] = -ini_b1ta[1]
    est_b1ta = trecase1irwsh_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=fit_b0ta)
    trecase_b1C[i,] = est_b1ta
    trecase_b1C[i,]
    
    #joined
    ini_b01ta = c(trecase_b1C[i,1:2],fit_b0ta,trecase_b1C[i,3:(ll0sh-1)])
    est_b01ta = trecase1irwsh_b01(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b01ta,est.inits=F)
    trecase_b01C[i,] = est_b01ta
    trecase_b01C[i,]
  
    #poo, at b0=0
    ini_b1ta = trecase_b01C[i,-c(3,llfsh)];#ini_b1ta[3]=0
    est_b1ta = trecase1irwsh_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=0)
    trecase_b1C0[i,] = est_b1ta#fit b1 assuming no b0
    trecase_b1C0[i,]
  
    #additive redo, at b1=0
    ini_b0ta = trecase_b01C[i,-c(4,llfsh)];#ini_b0ta[3]=0
    if(abs(ini_b0ta[3])>8)ini_b0ta[3]=sign(ini_b0ta[3])*8
    if(abs(ini_b0ta[4])>8)ini_b0ta[4]=sign(ini_b0ta[4])*8
    est_b0ta = trecase1irwsh_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
    delt0 = trecase_b0C[i,ll0sh]-est_b0ta[ll0sh]
    if((delt0> 1e-4) & traceit){
      message(i, " upd H0(b0=0)comm = ", format(delt0, scientific=T, digits=2))
      trecase_b0C[i,] = est_b0ta    
    }
  
    #f
    #joined
    ini_b01ta = trecase_b01C[i,-(nbeta+5)]
    trecase_b0C[i,nbeta+4]
    trecase_b1C0[i,nbeta+4]
    trecase_b01C[i, nbeta+5]
    if(trecase_b01C[i, llfsh]>trecase_b1C0[i,ll0sh]){
      ini_b01ta[1:2] = trecase_b1C0[i, 1:2]
      ini_b01ta[3] = 0
      ini_b01ta[4:(llfsh-1)] = trecase_b1C0[i, 3:(ll0sh-1)]      
    }
    if(trecase_b01C[i, llfsh]>trecase_b0C[i,ll0sh]){
      ini_b01ta[1:3] = trecase_b1C0[i, 1:3]
      ini_b01ta[4] = 0
      ini_b01ta[5:(llfsh-1)] = trecase_b1C0[i, 4:(ll0sh-1)]      
    }
    est_b01taf = trecase1irwsh_b01(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b01ta,est.inits=F)
    rerun = F
    if(est_b01taf[llfsh] < (trecase_b01C[i,llfsh])){
      rerun = T
      trecase_b01C[i,] = est_b01taf
    }

    if(rerun){  
      #poo, at b0=0
      ini_b1ta = trecase_b01C[i,-c(3,llfsh)]
      est_b1ta = trecase1irwsh_b1(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b1ta,est.inits=F,b0=0)
      if(est_b1ta[ll0sh] < (trecase_b1C0[i,ll0sh])){
        trecase_b1C0[i,] = est_b1ta#fit b1 assuming no b0
      }  
      #additive
      ini_b0ta = trecase_b01C[i,-c(4,llfsh)]
      est_b0ta = trecase1irwsh_b0(trc=trc,asn=asn,asnp=asnp,thp=thp,Xmatr=Xmatr,inits=ini_b0ta,est.inits=F)
      if(est_b0ta[ll0sh] < (trecase_b0C[i,ll0sh])){
        trecase_b0C[i,] = est_b0ta#fit b1 assuming no b0
      }  
    }
    if(i%%10==0)cat(i,"'th iteration\n")
  }
  apply(trecase_b01C,2,median)
  
  pval_b0 = 2*(trecase_b0C[, ll0sh]-trecase_b01C[,llfsh])
  pval_b1 = 2*(trecase_b1C0[, ll0sh]-trecase_b01C[,llfsh])
  if(traceit){
    message("clear fails for full likelihood: ", sum(pval_b0< -.1 | pval_b1< -.1))
  }
  fixb0 = pval_b0<0
  fixb1 = pval_b1<0
  pval_b0[fixb0] = 0
  pval_b1[fixb1] = 0
  pval_b0 = pchisq(pval_b0, df=1, lower.tail=F)
  pval_b1 = pchisq(pval_b1, df=1, lower.tail=F)
  
  list(full=cbind(trecase_b01C, pval_b0, pval_b1), testadd=trecase_b0C, testpoo=trecase_b1C0)
}
