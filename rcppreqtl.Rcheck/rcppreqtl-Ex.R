pkgname <- "rcppreqtl"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('rcppreqtl')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("fit")
### * fit

flush(stderr()); flush(stdout())

### Name: fit
### Title: Optimization wrapper, maximizing the joint model of total (TReC)
###   and allele specific (ASE) counts for autosomes
### Aliases: fit
### Keywords: methods

### ** Examples
## Not run: 
##D # fitting autosome data for a full model with allele-specific counts collected on gene level:
##D percase = 0.1
##D dblcnt = 0.2
##D mn = 100
##D b0 = 0
##D b1 = 0
##D phiNB = .5;
##D phiBB=phiNB/4
##D niter = 10
##D betas = c(3,.2,.05,.5)
##D ss=2
##D 
##D dep = makeXmatr(ss)
##D dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
##D fit(subset=NULL, data=dat, traceit=FALSE)
## End(Not run)


cleanEx()
nameEx("fitsh")
### * fitsh

flush(stderr()); flush(stdout())

### Name: fitsh
### Title: Optimization wrapper, maximizing the joint model of total (TReC)
###   and allele specific (ASE) counts for autosomes
### Aliases: fitsh
### Keywords: methods

### ** Examples
## Not run: 
##D # fitting autosome data for a full model with allele-specific counts collected on gene level:
##D percase = 0.1
##D dblcnt = 0.2
##D mn = 100
##D b0 = 0
##D b1 = 0
##D phiNB = .5;
##D phiBB=phiNB/4
##D niter = 10
##D betas = c(3,.2,.05,.5)
##D ss=2
##D 
##D dep = makeXmatr(ss)
##D dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
##D fitsh(subset=NULL, data=dat, traceit=FALSE)
## End(Not run)


cleanEx()
nameEx("makeXmatr")
### * makeXmatr

flush(stderr()); flush(stdout())

### Name: makeXmatr
### Title: Create example design matrix for simulations
### Aliases: makeXmatr
### Keywords: methods

### ** Examples
## Not run: 
##D # fitting autosome data for a full model with allele-specific counts collected on gene level:
##D percase = 0.1
##D dblcnt = 0.2
##D mn = 100
##D b0 = 0
##D b1 = 0
##D phiNB = .5;
##D phiBB=phiNB/4
##D niter = 100
##D betas = c(3,.2,.05,.5)
##D ss=2
##D 
##D dep = makeXmatr(ss)
## End(Not run)


cleanEx()
nameEx("readCounts")
### * readCounts

flush(stderr()); flush(stdout())

### Name: readCounts
### Title: A list object that should be used as input to optimization fit
###   or fitsh function.
### Aliases: readCounts
### Keywords: utilities

### ** Examples
## Not run: 
##D # see total read counts (TReC) for first 2 X chromosome genes of a data example:
##D rc = readCounts(haplotype=haplotypef, trc=trc, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
##D              X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
## End(Not run)


cleanEx()
nameEx("simu2")
### * simu2

flush(stderr()); flush(stdout())

### Name: simu2
### Title: Simulate a dataset in a format acceptable by a fit function
### Aliases: simu2
### Keywords: methods

### ** Examples
## Not run: 
##D # fitting autosome data for a full model with allele-specific counts collected on gene level:
##D percase = 0.1
##D dblcnt = 0.2
##D mn = 100
##D b0 = 0
##D b1 = 0
##D phiNB = .5;
##D phiBB=phiNB/4
##D niter = 100
##D betas = c(3,.2,.05,.5)
##D ss=2
##D 
##D dep = makeXmatr(ss)
##D dat = simu2(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
## End(Not run)


cleanEx()
nameEx("simu4")
### * simu4

flush(stderr()); flush(stdout())

### Name: simu4
### Title: Simulate a dataset in a format acceptable by a fit function
### Aliases: simu4
### Keywords: methods

### ** Examples
## Not run: 
##D # fitting autosome data for a full model with allele-specific counts collected on gene level:
##D percase = 0.1
##D dblcnt = 0.2
##D mn = 100
##D b0 = 0
##D b1 = 0
##D phiNB = .5;
##D phiBB=phiNB/4
##D niter = 100
##D betas = c(3,.2,.05,.5)
##D ss=2
##D 
##D dep = makeXmatr(ss)
##D dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
## End(Not run)


### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
