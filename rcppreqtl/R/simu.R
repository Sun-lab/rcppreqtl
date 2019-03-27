
makeXmatr = function(ss){
  thp = rep(0:3,8)
  nind = length(thp)
  Xmatr = matrix(c(rep(1,nind),rep(0:1,each=nind/2),rep(0:3,each=nind/4),rnorm(nind)),ncol=4)
  
  if(ss>1){
    for(j in 2:ss){
      Xmatr2 = Xmatr;Xmatr2[,4]=rnorm(nind)
      Xmatr = rbind(Xmatr,Xmatr2)
      thp = c(thp,thp)
      nind = nrow(Xmatr)
    }
  }
  list(Xmatr=Xmatr,  thp=thp)
}


simu2 = function(num, Xmatr, haplotype, totmean, percase=0.1, dblcnt=0, phiNB=1, phiBB=0.5, b0=0, b1=0, betas=rep(1,4)){
  #library(VGAM)
  #library(MASS)
  
  rho = phiBB/(1+phiBB)
  
  nind = length(haplotype)
  nind2 = nrow(Xmatr)
  if(nind!=nind2)stop("haplotype and Xmatr dimension mismatch")
  nbeta = ncol(Xmatr)
  nbeta2 = length(betas)
  if(nbeta!=nbeta2)stop("betas and Xmatr dimension mismatch")

  #total read count structure
  etas <- rep(0, nind) #base will be eta3 <- 0 - AxA
  etas[haplotype==1] = log1p(exp(b0+b1)) - log1p(exp(b1))#AxB
  etas[haplotype==2] = log1p(exp(b0-b1)) - log1p(exp(-b1))#AxB
  etas[haplotype==3] = b0 #BxB
  lmu <- Xmatr%*%betas + etas
  mu = exp(lmu[,1])
  mu = mu/mean(mu)*totmean
  #allele specific count structure
  prb = rep(0, nind)
  prb[haplotype==0] = b1
  prb[haplotype==1] = b0 + b1
  prb[haplotype==2] = -b0 + b1
  prb[haplotype==3] = b1
  logiti = function(x){1/(1+exp(-x))}
  prb = logiti(prb)

  asn10 = asn20 = asn30 = asn40 = asnp10 = asnp20 = asnp30 = asnp40 =
  asn4m = asnp4m = asn3m = asnp3m = asn2m = asnp2m = asn1m = asnp1m = trcm = matrix(NA,nrow=num,ncol=nind)
  #having
  # np/(alpha+beta)= n p and
  #(alpha+beta+n)/(alpha+beta+1)=[1+(n-1)rho]
  #gives
  #alpha = p(1-rho)/rho
  #beta = (1-p)(1-rho)/rho
  if(rho>0){
    alphasBB = prb*(1-rho)/rho
    betasBB = (1-prb)*(1-rho)/rho
  }

  for(i in 1:num){
    if(rho>0){
      prb1 = rbeta(nind, alphasBB, betasBB)
    }else{
      prb1 = prb
    }
  
    trc = rnbinom(n=nind,size=1/phiNB,mu=mu)
    asn0 = round(trc * percase)
        
    asn1 = round(.5*asn0)
    asn2 = asn0 - asn1
    
    asn4d = asn3d = asn2d = asn1d =
    asnp4d = asnp3d = asnp2d = asnp1d =
    asnp4 = asnp3  = asnp1 = asnp2  = rep(0, nind)    
    if(rho>0){
      prb1 = rbeta(nind, alphasBB, betasBB)
    }else{
      prb1 = prb
    }
    asnp1 = rbinom(nind, size=asn1,prob=prb1)
    asnp2 = rbinom(nind, size=asn2,prob=prb1)
    

    trcm[i,] = trc

    asn1m[i,] = asn10[i,] = asn1
    asn2m[i,] = asn20[i,] = asn2

    asnp1m[i,] = asnp10[i,] = asnp1
    asnp2m[i,] = asnp20[i,] = asnp2

    if(dblcnt>0){
      to.flip = rbinom(n=nind, size=asn1, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp1d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp1/asn1)[kp])
      asn1d = to.flip
      }
    
      to.flip = rbinom(n=nind, size=asn2, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp2d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp2/asn2)[kp])
      asn2d = to.flip
      }
    }

    asn1m[i,] = asn1m[i,] + asn2d
    asnp1m[i,] = asnp1m[i,] + asnp2d
    asn2m[i,] = asn2m[i,] + asn1d
    asnp2m[i,] = asnp2m[i,] + asnp1d
  }

  #one gene with no double-counting
  asnf = asn10 + asn20
  asnpf = asnp10 + asnp20
  haplotypef = haplotype

  asnm = cbind(asn1m, asn2m)
  asnpm = cbind(asnp1m, asnp2m)
  haplotyped = c(haplotype, haplotype)
  rownames(asnf) = rownames(asnpf) = rownames(asnm) = rownames(asnpm) = rownames(trcm) = sprintf("gene%s", 1:num)

  haplotyped = matrix(rep(haplotyped, each=num), nrow=num)#multisnp
  haplotypef = matrix(rep(haplotypef, each=num), nrow=num)

  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, 
             haplotype2S=haplotyped, asn2S=asnm, asnp2S=asnpm,
             haplotype4S=NULL, asn4S=NULL, asnp4S=NULL,
             haplotype8S=NULL, asn8S=NULL, asnp8S=NULL,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
#  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
#             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}


simu4= function(num, Xmatr, haplotype, totmean, percase=0.1, dblcnt=0, phiNB=1, phiBB=0.5, b0=0, b1=0, betas=rep(1,4)){
  #library(VGAM)
  #library(MASS)
  
  rho = phiBB/(1+phiBB)
  
  nind = length(haplotype)
  nind2 = nrow(Xmatr)
  if(nind!=nind2)stop("haplotype and Xmatr dimension mismatch")
  nbeta = ncol(Xmatr)
  nbeta2 = length(betas)
  if(nbeta!=nbeta2)stop("betas and Xmatr dimension mismatch")

  #total read count structure
  etas <- rep(0, nind) #base will be eta3 <- 0 - AxA
  etas[haplotype==1] = log1p(exp(b0+b1)) - log1p(exp(b1))#AxB
  etas[haplotype==2] = log1p(exp(b0-b1)) - log1p(exp(-b1))#AxB
  etas[haplotype==3] = b0 #BxB
  lmu <- Xmatr%*%betas + etas
  mu = exp(lmu[,1])
  mu = mu/mean(mu)*totmean
  #allele specific count structure
  prb = rep(0, nind)
  prb[haplotype==0] = b1
  prb[haplotype==1] = b0 + b1
  prb[haplotype==2] = -b0 + b1
  prb[haplotype==3] = b1
  logiti = function(x){1/(1+exp(-x))}
  prb = logiti(prb)

  asn10 = asn20 = asn30 = asn40 = 
  asnp10 = asnp20 = asnp30 = asnp40 =
  asn4m = asnp4m = asn3m = asnp3m = asn2m = asnp2m = asn1m = asnp1m = trcm = matrix(NA,nrow=num,ncol=nind)
  #having
  # np/(alpha+beta)= n p and
  #(alpha+beta+n)/(alpha+beta+1)=[1+(n-1)rho]
  #gives
  #alpha = p(1-rho)/rho
  #beta = (1-p)(1-rho)/rho
  if(rho>0){
    alphasBB = prb*(1-rho)/rho
    betasBB = (1-prb)*(1-rho)/rho
  }

  for(i in 1:num){
    if(rho>0){
      prb1 = rbeta(nind, alphasBB, betasBB)
    }else{
      prb1 = prb
    }
  
    trc = rnbinom(n=nind,size=1/phiNB,mu=mu)
    asn0 = round(trc * percase)
        
    asn1 = round(.5*asn0)
    asn2 = asn0 - asn1
    
    asn3 = round(.5*asn1)
    asn4 = round(.5*asn2)
    asn1 = asn1-asn3
    asn2 = asn2-asn4
    asn4d = asn3d = asn2d = asn1d =
    asnp4d = asnp3d = asnp2d = asnp1d =
    asnp4 = asnp3  = asnp1 = asnp2  = rep(0, nind)    
    if(rho>0){
      prb1 = rbeta(nind, alphasBB, betasBB)
    }else{
      prb1 = prb
    }
    asnp1 = rbinom(nind, size=asn1,prob=prb1)
    asnp2 = rbinom(nind, size=asn2,prob=prb1)
    asnp3 = rbinom(nind, size=asn3,prob=prb1)
    asnp4 = rbinom(nind, size=asn4,prob=prb1)    

    trcm[i,] = trc

    asn1m[i,] = asn10[i,] = asn1
    asn2m[i,] = asn20[i,] = asn2
    asn3m[i,] = asn30[i,] = asn3
    asn4m[i,] = asn40[i,] = asn4

    asnp1m[i,] = asnp10[i,] = asnp1
    asnp2m[i,] = asnp20[i,] = asnp2
    asnp3m[i,] = asnp30[i,] = asnp3
    asnp4m[i,] = asnp40[i,] = asnp4

    if(dblcnt>0){
      to.flip = rbinom(n=nind, size=asn1, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp1d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp1/asn1)[kp])
      asn1d = to.flip
      }
    
      to.flip = rbinom(n=nind, size=asn2, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp2d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp2/asn2)[kp])
      asn2d = to.flip
      }
    
      to.flip = rbinom(n=nind, size=asn3, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp3d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp3/asn3)[kp])
        asn3d = to.flip
      }
  
      to.flip = rbinom(n=nind, size=asn4, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp4d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp4/asn4)[kp])
        asn4d = to.flip
      }  
    }

    asn1m[i,] = asn1m[i,] + asn4d
    asnp1m[i,] = asnp1m[i,] + asnp4d
    asn2m[i,] = asn2m[i,] + asn1d
    asnp2m[i,] = asnp2m[i,] + asnp1d
    asn3m[i,] = asn3m[i,] + asn2d
    asnp3m[i,] = asnp3m[i,] + asnp2d
    asn4m[i,] = asn4m[i,] + asn3d
    asnp4m[i,] = asnp4m[i,] + asnp3d
  }

  #one gene with no double-counting
  asnf = asn10 + asn20 + asn30 + asn40
  asnpf = asnp10 + asnp20 + asnp30 + asnp40
  haplotypef = haplotype

  asnm = cbind(asn1m, asn2m, asn3m, asn4m)
  asnpm = cbind(asnp1m, asnp2m, asnp3m, asnp4m)
  haplotyped = c(haplotype, haplotype, haplotype, haplotype)
  
  rownames(asnf) = rownames(asnpf) = rownames(asnm) = rownames(asnpm) = rownames(trcm) = sprintf("gene%s", 1:num)
  haplotype1S = matrix(rep(haplotypef, each=num), nrow=num)
  haplotype4S = matrix(rep(haplotyped, each=num), nrow=num)

  readCounts(haplotype=haplotype1S, trc=trcm, asn=asnf, asnp=asnpf, 
             haplotype2S=NULL, asn2S=NULL, asnp2S=NULL,
             haplotype4S=haplotype4S, asn4S=asnm, asnp4S=asnpm,
             haplotype8S=NULL, asn8S=NULL, asnp8S=NULL,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
#  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
#             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}



simu8 = function(num, Xmatr, haplotype, totmean, percase=0.1, dblcnt=0, phiNB=1, phiBB=0.5, b0=0, b1=0, betas=rep(1,4)){
#num=niter;Xmatr=Xm$Xmatr; haplotype=Xm$thp; totmean=100; percase=0.2; dblcnt=0.2; phiNB=1; phiBB=0.5; b0=0; b1=0; betas=rep(1,4)
#dblcnt=0
  #library(VGAM)
  #library(MASS)
  
  rho = phiBB/(1+phiBB)
  
  nind = length(haplotype)
  nind2 = nrow(Xmatr)
  if(nind!=nind2)stop("haplotype and Xmatr dimension mismatch")
  nbeta = ncol(Xmatr)
  nbeta2 = length(betas)
  if(nbeta!=nbeta2)stop("betas and Xmatr dimension mismatch")

  #total read count structure
  etas <- rep(0, nind) #base will be eta3 <- 0 - AxA
  etas[haplotype==1] = log1p(exp(b0+b1)) - log1p(exp(b1))#AxB
  etas[haplotype==2] = log1p(exp(b0-b1)) - log1p(exp(-b1))#AxB
  etas[haplotype==3] = b0 #BxB
  lmu <- Xmatr%*%betas + etas
  mu = exp(lmu[,1])
  mu = mu/mean(mu)*totmean
  #allele specific count structure
  prb = rep(0, nind)
  prb[haplotype==0] = b1
  prb[haplotype==1] = b0 + b1
  prb[haplotype==2] = -b0 + b1
  prb[haplotype==3] = b1
  logiti = function(x){1/(1+exp(-x))}
  prb = logiti(prb)

  asn10 = asn20 = asn30 = asn40 =  asn50 = asn60 = asn70 = asn80 = 
  asnp10 = asnp20 = asnp30 = asnp40 = asnp50 = asnp60 = asnp70 = asnp80 =
 
  asn8dm = asnp8dm = asn7dm = asnp7dm = asn6dm = asnp6dm = asn5dm = asnp5dm = 
  asn4dm = asnp4dm = asn3dm = asnp3dm = asn2dm = asnp2dm = asn1dm = asnp1dm = 
  trcm = matrix(0,nrow=num,ncol=nind)
  #having
  # np/(alpha+beta)= n p and
  #(alpha+beta+n)/(alpha+beta+1)=[1+(n-1)rho]
  #gives
  #alpha = p(1-rho)/rho
  #beta = (1-p)(1-rho)/rho
  if(rho>0){
    alphasBB = prb*(1-rho)/rho
    betasBB = (1-prb)*(1-rho)/rho
  }

  for(i in 1:num){
    if(rho>0){
      prb1 = rbeta(nind, alphasBB, betasBB)
    }else{
      prb1 = prb
    }
  
    trc = rnbinom(n=nind,size=1/phiNB,mu=mu)
    asn0 = round(trc * percase)
        
    asn1 = round(.5*asn0)
    asn2 = asn0 - asn1
    
    asn3 = round(.5*asn1)
    asn4 = round(.5*asn2)

    asn1 = asn1-asn3
    asn2 = asn2-asn4

    asn5 = round(.5*asn1)
    asn1 = asn1-asn5
    asn6 = round(.5*asn2)
    asn2 = asn2-asn6
    asn7 = round(.5*asn3)
    asn3 = asn3-asn7
    asn8 = round(.5*asn4)
    asn4 = asn4-asn8
    
    asnp8 = asnp7 = asnp6 = asnp5 = asnp4 = asnp3 = asnp2 = asnp1 = rep(NA, nind)    
    
    asn8d = asn7d = asn6d = asn5d = asn4d = asn3d = asn2d = asn1d = 
    asnp4d = asnp3d = asnp2d = asnp1d = asnp5d = asnp6d = asnp7d = asnp8d = rep(0, nind)    

    asnp1 = rbinom(nind, size=asn1,prob=prb1)
    asnp2 = rbinom(nind, size=asn2,prob=prb1)
    asnp3 = rbinom(nind, size=asn3,prob=prb1)
    asnp4 = rbinom(nind, size=asn4,prob=prb1)
    asnp5 = rbinom(nind, size=asn5,prob=prb1)
    asnp6 = rbinom(nind, size=asn6,prob=prb1)
    asnp7 = rbinom(nind, size=asn7,prob=prb1)
    asnp8 = rbinom(nind, size=asn8,prob=prb1)

    trcm[i,] = trc

    asn10[i,] = asn1
    asn20[i,] = asn2
    asn30[i,] = asn3
    asn40[i,] = asn4
    asn50[i,] = asn5
    asn60[i,] = asn6
    asn70[i,] = asn7
    asn80[i,] = asn8

    asnp10[i,] = asnp1
    asnp20[i,] = asnp2
    asnp30[i,] = asnp3
    asnp40[i,] = asnp4
    asnp50[i,] = asnp5
    asnp60[i,] = asnp6
    asnp70[i,] = asnp7
    asnp80[i,] = asnp8
  }

  #4 SNP with no double-counting
  asn4S1m = asn10 + asn20
  asn4S2m = asn30 + asn40
  asn4S3m = asn50 + asn60
  asn4S4m = asn70 + asn80

  asnp4S1m = asnp10 + asnp20
  asnp4S2m = asnp30 + asnp40
  asnp4S3m = asnp50 + asnp60
  asnp4S4m = asnp70 + asnp80
  
  #2 SNP with no double-counting
  asn2S1m = asn4S1m + asn4S2m
  asn2S2m = asn4S3m + asn4S4m

  asnp2S1m = asnp4S1m + asnp4S2m
  asnp2S2m = asnp4S3m + asnp4S4m

  #2 SNP with no double-counting
  asn2S1m = asn4S1m + asn4S2m
  asn2S2m = asn4S3m + asn4S4m

  asnp2S1m = asnp4S1m + asnp4S2m
  asnp2S2m = asnp4S3m + asnp4S4m

  asn1S1m = asn2S1m + asn2S2m
  asnp1S1m = asnp2S1m + asnp2S2m
  
  #add double-counting
  for(i in 1:num){
    if(dblcnt>0){
      to.flip = rbinom(n=nind, size=asn1, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp1d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp1/asn1)[kp])
      asn1d = to.flip
      }
      asn1dm[i,] = asn1d
      asnp1dm[i,] = asnp1d
    
      to.flip = rbinom(n=nind, size=asn2, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp2d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp2/asn2)[kp])
      asn2d = to.flip
      }
      asn2dm[i,] = asn2d
      asnp2dm[i,] = asnp2d
    
      to.flip = rbinom(n=nind, size=asn3, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp3d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp3/asn3)[kp])
        asn3d = to.flip
      }
      asn3dm[i,] = asn3d
      asnp3dm[i,] = asnp3d
  
      to.flip = rbinom(n=nind, size=asn4, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp4d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp4/asn4)[kp])
        asn4d = to.flip
      }  
      asn4dm[i,] = asn4d
      asnp4dm[i,] = asnp4d

      to.flip = rbinom(n=nind, size=asn5, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp5d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp5/asn5)[kp])
        asn5d = to.flip
      }  
      asn5dm[i,] = asn5d
      asnp5dm[i,] = asnp5d
      
      to.flip = rbinom(n=nind, size=asn6, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp6d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp6/asn6)[kp])
        asn6d = to.flip
      }  
      asn6dm[i,] = asn6d
      asnp6dm[i,] = asnp6d

      to.flip = rbinom(n=nind, size=asn7, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp7d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp7/asn7)[kp])
        asn7d = to.flip
      }  
      asn7dm[i,] = asn7d
      asnp7dm[i,] = asnp7d

      to.flip = rbinom(n=nind, size=asn8, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp8d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp8/asn8)[kp])
        asn8d = to.flip
      }  
      asn8dm[i,] = asn8d
      asnp8dm[i,] = asnp8d

    }
  }
  #add double-counting 8SNP
  asn1m = asn10 + asn8dm
  asnp1m = asnp10 + asnp8dm
  asn2m = asn20 + asn1dm
  asnp2m = asnp20 + asnp1dm
  asn3m = asn30 + asn2dm
  asnp3m = asnp30 + asnp2dm
  asn4m = asn40 + asn3dm
  asnp4m = asnp40 + asnp3dm
  asn5m = asn50 + asn4dm
  asnp5m = asnp50 + asnp4dm
  asn6m = asn60 + asn5dm
  asnp6m = asnp60 + asnp5dm
  asn7m = asn70 + asn6dm
  asnp7m = asnp70 + asnp6dm
  asn8m = asn80 + asn7dm
  asnp8m = asnp80 + asnp7dm
  sum(asn1m)/sum(asn10)
  sum(asn2m)/sum(asn20)
  sum(asn3m)/sum(asn30)
  sum(asn4m)/sum(asn40)
  sum(asn5m)/sum(asn50)
  sum(asn6m)/sum(asn60)
  sum(asn7m)/sum(asn70)
  sum(asn8m)/sum(asn80)

  #add double-counting 4SNP
  asn4S1d = asn4S1m + asn7dm + asn8dm
  asnp4S1d = asnp4S1m + asnp7dm + asnp8dm
  asn4S2d = asn4S2m + asn1dm + asn2dm
  asnp4S2d = asnp4S2m + asnp1dm + asnp2dm
  asn4S3d = asn4S3m + asn3dm + asn4dm
  asnp4S3d = asnp4S3m + asnp3dm + asnp4dm
  asn4S4d = asn4S4m + asn5dm + asn6dm
  asnp4S4d = asnp4S4m + asnp5dm + asnp6dm
  sum(asn4S1d)/sum(asn4S1m)
  sum(asn4S2d)/sum(asn4S2m)
  sum(asn4S3d)/sum(asn4S3m)
  sum(asn4S4d)/sum(asn4S4m)

  #add double-counting 2SNP
  asn2S1d = asn2S1m + asn5dm + asn6dm + asn7dm + asn8dm
  asnp2S1d = asnp2S1m + asnp5dm + asnp6dm + asnp7dm + asnp8dm
  asn2S2d = asn2S2m + asn1dm + asn2dm + asn3dm + asn4dm
  asnp2S2d = asnp2S2m + asnp1dm + asnp2dm + asnp3dm + asnp4dm
  sum(asn2S1d)/sum(asn2S1m)
  sum(asn2S2d)/sum(asn2S2m)
  
  haplotype1S = haplotype
  haplotype2S = c(haplotype1S, haplotype1S)
  haplotype4S = c(haplotype2S, haplotype2S)
  haplotype8S = c(haplotype4S, haplotype4S)

  asnm8 = cbind(asn1m, asn2m, asn3m, asn4m, asn5m, asn6m, asn7m, asn8m)
  asnpm8 = cbind(asnp1m, asnp2m, asnp3m, asnp4m, asnp5m, asnp6m, asnp7m, asnp8m)

  #asnm4 = cbind(asn4S1m, asn4S2m, asn4S3m, asn4S4m)
  #asnpm4 = cbind(asnp4S1m, asnp4S2m, asnp4S3m, asnp4S4m)
  asnm4 = cbind(asn4S1d, asn4S1d, asn4S1d, asn4S1d)
  asnpm4 = cbind(asnp4S1d, asnp4S2d, asnp4S3d, asnp4S4d)

  #asnm2 = cbind(asn2S1m, asn2S2m)
  #asnpm2 = cbind(asnp2S1m, asnp2S2m)
  asnm2 = cbind(asn2S1d, asn2S2d)
  asnpm2 = cbind(asnp2S1d, asnp2S2d)
  
  asnm1 = asn1S1m
  asnpm1 = asnp2S1m

  rownames(asnm1) = rownames(asnpm1) = 
  rownames(asnm2) = rownames(asnpm2) = 
  rownames(asnm4) = rownames(asnpm4) = 
  rownames(asnm8) = rownames(asnpm8) =  rownames(trcm) = sprintf("gene%s", 1:num)
  haplotype1S = matrix(rep(haplotype1S, each=num), nrow=num)
  haplotype2S = matrix(rep(haplotype2S, each=num), nrow=num)
  haplotype4S = matrix(rep(haplotype4S, each=num), nrow=num)
  haplotype8S = matrix(rep(haplotype8S, each=num), nrow=num)

  readCounts(haplotype=haplotype1S, trc=trcm, asn=asnm1, asnp=asnpm1, 
             haplotype2S=haplotype2S, asn2S=asnm2, asnp2S=asnpm2,
             haplotype4S=haplotype4S, asn4S=asnm4, asnp4S=asnpm4,
             haplotype8S=haplotype8S, asn8S=asnm8, asnp8S=asnpm8,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}

