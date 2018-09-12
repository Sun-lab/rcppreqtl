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
  etas <- rep(0, length(haplotype)) #base will be eta3 <- 0 - AxA
  etas[haplotype==1] = log1p(exp(b0+b1)) - log1p(exp(b1))#AxB
  etas[haplotype==2] = log1p(exp(b0-b1)) - log1p(exp(-b1))#AxB
  etas[haplotype==3] = b0 #BxB
  lmu <- Xmatr%*%betas + etas
  mu = exp(lmu[,1])
  mu = mu/mean(mu)*totmean
  #allele specific count structure
  prb = rep(0, length(haplotype))
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
      prb1 = rbeta(length(prb), alphasBB, betasBB)
    }else{
      prb1 = prb
    }
  
    trc = rnbinom(n=length(lmu),size=1/phiNB,mu=mu)
    asn0 = round(trc * percase)
        
    asn1 = round(.5*asn0)
    asn2 = asn0 - asn1
    
    asn3 = round(.5*asn1)
    asn4 = round(.5*asn2)
    asn1 = asn1-asn3
    asn2 = asn2-asn4
    asn4d = asn3d = asn2d = asn1d =
    asnp4d = asnp3d = asnp2d = asnp1d =
    asnp4 = asnp3  = asnp1 = asnp2  = rep(0, length(asn1))    
    if(rho>0){
      asnp1 = rbinom(length(asnp1), size=asn1,prob=prb1)
      asnp2 = rbinom(length(asnp2), size=asn2,prob=prb1)
      asnp3 = rbinom(length(asnp3), size=asn3,prob=prb1)
      asnp4 = rbinom(length(asnp4), size=asn4,prob=prb1)    
    }else{
      asnp1 = rbinom(length(asnp1), size=asn1, prob=prb)
      asnp2 = rbinom(length(asnp2), size=asn2, prob=prb)
      asnp3 = rbinom(length(asnp3), size=asn3, prob=prb)
      asnp4 = rbinom(length(asnp4), size=asn4, prob=prb)
    }

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
      to.flip = rbinom(n=length(asnp1), size=asn1, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp1d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp1/asn1)[kp])
      asn1d = to.flip
      }
    
      to.flip = rbinom(n=length(asnp2), size=asn2, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp2d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp2/asn2)[kp])
      asn2d = to.flip
      }
    
      to.flip = rbinom(n=length(asnp3), size=asn3, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
        asnp3d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp3/asn3)[kp])
        asn3d = to.flip
      }
  
      to.flip = rbinom(n=length(asnp4), size=asn4, prob=dblcnt)
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
  haplotyped = matrix(rep(haplotyped, each=num), nrow=num)
  haplotypef = matrix(rep(haplotypef, each=num), nrow=num)

  readCounts(haplotype=haplotypef, trc=trcm, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
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
  etas <- rep(0, length(haplotype)) #base will be eta3 <- 0 - AxA
  etas[haplotype==1] = log1p(exp(b0+b1)) - log1p(exp(b1))#AxB
  etas[haplotype==2] = log1p(exp(b0-b1)) - log1p(exp(-b1))#AxB
  etas[haplotype==3] = b0 #BxB
  lmu <- Xmatr%*%betas + etas
  mu = exp(lmu[,1])
  mu = mu/mean(mu)*totmean
  #allele specific count structure
  prb = rep(0, length(haplotype))
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
      prb1 = rbeta(length(prb), alphasBB, betasBB)
    }else{
      prb1 = prb
    }
  
    trc = rnbinom(n=length(lmu),size=1/phiNB,mu=mu)
    asn0 = round(trc * percase)
        
    asn1 = round(.5*asn0)
    asn2 = asn0 - asn1
    
    asn4d = asn3d = asn2d = asn1d =
    asnp4d = asnp3d = asnp2d = asnp1d =
    asnp4 = asnp3  = asnp1 = asnp2  = rep(0, length(asn1))    
    if(rho>0){
      asnp1 = rbinom(length(asnp1), size=asn1,prob=prb1)
      asnp2 = rbinom(length(asnp2), size=asn2,prob=prb1)
    }else{
      asnp1 = rbinom(length(asnp1), size=asn1, prob=prb)
      asnp2 = rbinom(length(asnp2), size=asn2, prob=prb)
    }

    trcm[i,] = trc

    asn1m[i,] = asn10[i,] = asn1
    asn2m[i,] = asn20[i,] = asn2

    asnp1m[i,] = asnp10[i,] = asnp1
    asnp2m[i,] = asnp20[i,] = asnp2

    if(dblcnt>0){
      to.flip = rbinom(n=length(asnp1), size=asn1, prob=dblcnt)
      kp = which(to.flip>0)
      if(length(kp)>0){
      asnp1d[kp] = rbinom(n=length(kp), size=to.flip[kp], prob=(asnp1/asn1)[kp])
      asn1d = to.flip
      }
    
      to.flip = rbinom(n=length(asnp2), size=asn2, prob=dblcnt)
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