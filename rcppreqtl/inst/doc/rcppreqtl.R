### R code from vignette source 'rcppreqtl.Rnw'

###################################################
### code chunk number 1: initialize
###################################################
library(rcppreqtl, lib.loc="/nas02/home/z/h/zhabotyn/research/package9/")


###################################################
### code chunk number 2: rcppreqtl.Rnw:31-32
###################################################
options(width = 80)


###################################################
### code chunk number 3: rcppreqtl.Rnw:63-76
###################################################
percase = 0.1
dblcnt = 0.2
mn = 100
b0 = 0;b1 = 0;th = .5;dv=4;niter = 100;betas = c(3,.2,.05,.5);ss=2
set.seed(12345)
library(VGAM)
library(MASS)
phiNB = th
phiBB = th/dv
dep = makeXmatr(ss)
dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, 
            percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, 
            b0=b0, b1=b1, betas=betas)


###################################################
### code chunk number 4: rcppreqtl.Rnw:81-83
###################################################
#fit trecase autosome genes:
fullest = fit(subset=1:2, data=dat, traceit=FALSE)


