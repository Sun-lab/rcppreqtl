
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

  haplotyped = matrix(rep(haplotyped, each=num), nrow=num)
  haplotypef = matrix(rep(haplotypef, each=num), nrow=num)

  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
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
    asnp3m[i,] = asnp3m[i,] + asnp1d
    asn4m[i,] = asn4m[i,] + asn4d
    asnp4m[i,] = asnp4m[i,] + asnp4d
  }

  #one gene with no double-counting
  asnf = asn10 + asn20 + asn30 + asn40
  asnpf = asnp10 + asnp20 + asnp30 + asnp40
  haplotypef = haplotype

  asnm = cbind(asn1m, asn2m, asn3m, asn4m)
  asnpm = cbind(asnp1m, asnp2m, asnp3m, asnp4m)
  haplotyped = c(haplotype, haplotype, haplotype, haplotype)
  
  rownames(asnf) = rownames(asnpf) = rownames(asnm) = rownames(asnpm) = rownames(trcm) = sprintf("gene%s", 1:num)
  haplotyped = matrix(rep(haplotyped, each=num), nrow=num)
  haplotypef = matrix(rep(haplotypef, each=num), nrow=num)

  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}



simu8 = function(num, Xmatr, haplotype, totmean, percase=0.1, dblcnt=0, phiNB=1, phiBB=0.5, b0=0, b1=0, betas=rep(1,4)){
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
  asn8m = asnp8m = asn7m = asnp7m = asn6m = asnp6m = asn5m = asnp5m = 
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

    asn1m[i,] = asn10[i,] = asn1
    asn2m[i,] = asn20[i,] = asn2
    asn3m[i,] = asn30[i,] = asn3
    asn4m[i,] = asn40[i,] = asn4
    asn5m[i,] = asn50[i,] = asn5
    asn6m[i,] = asn60[i,] = asn6
    asn7m[i,] = asn70[i,] = asn7
    asn8m[i,] = asn80[i,] = asn8

    asnp1m[i,] = asnp10[i,] = asnp1
    asnp2m[i,] = asnp20[i,] = asnp2
    asnp3m[i,] = asnp30[i,] = asnp3
    asnp4m[i,] = asnp40[i,] = asnp4
    asnp5m[i,] = asnp50[i,] = asnp5
    asnp6m[i,] = asnp60[i,] = asnp6
    asnp7m[i,] = asnp70[i,] = asnp7
    asnp8m[i,] = asnp80[i,] = asnp8

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

      to.flip = rbinom(n=nind, size=asn5, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp5d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp5/asn5)[kp])
        asn5d = to.flip
      }  
      
      to.flip = rbinom(n=nind, size=asn6, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp6d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp6/asn6)[kp])
        asn6d = to.flip
      }  
      to.flip = rbinom(n=nind, size=asn7, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp7d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp7/asn7)[kp])
        asn7d = to.flip
      }  
      to.flip = rbinom(n=nind, size=asn8, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp8d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp8/asn8)[kp])
        asn8d = to.flip
      }  

    }

    asn1m[i,] = asn1m[i,] + asn4d
    asnp1m[i,] = asnp1m[i,] + asnp4d
    asn2m[i,] = asn2m[i,] + asn1d
    asnp2m[i,] = asnp2m[i,] + asnp1d
    asn3m[i,] = asn3m[i,] + asn2d
    asnp3m[i,] = asnp3m[i,] + asnp1d
    asn4m[i,] = asn4m[i,] + asn4d
    asnp4m[i,] = asnp4m[i,] + asnp4d

    asn5m[i,] = asn5m[i,] + asn5d
    asnp5m[i,] = asnp5m[i,] + asnp5d
    asn6m[i,] = asn6m[i,] + asn6d
    asnp6m[i,] = asnp6m[i,] + asnp6d
    asn7m[i,] = asn7m[i,] + asn7d
    asnp7m[i,] = asnp7m[i,] + asnp7d
    asn8m[i,] = asn8m[i,] + asn8d
    asnp8m[i,] = asnp8m[i,] + asnp8d
  }

  #one gene with no double-counting
  asnf = asn10 + asn20 + asn30 + asn40 + asn50 + asn60 + asn70 + asn80
  asnpf = asnp10 + asnp20 + asnp30 + asnp40 + asnp50 + asnp60 + asnp70 + asnp80
  haplotypef = haplotype

  asnm = cbind(asn1m, asn2m, asn3m, asn4m, asn5m, asn6m, asn7m, asn8m)
  asnpm = cbind(asnp1m, asnp2m, asnp3m, asnp4m, asnp5m, asnp6m, asnp7m, asnp8m)
  haplotyped = c(haplotype, haplotype, haplotype, haplotype, haplotype, haplotype, haplotype, haplotype)
  
  rownames(asnf) = rownames(asnpf) = rownames(asnm) = rownames(asnpm) = rownames(trcm) = sprintf("gene%s", 1:num)
  haplotyped = matrix(rep(haplotyped, each=num), nrow=num)
  haplotypef = matrix(rep(haplotypef, each=num), nrow=num)

  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}

