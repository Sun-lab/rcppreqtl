\name{fitsh}
\alias{fitsh}
\title{Optimization wrapper, maximizing the joint model of total (TReC) and allele specific (ASE) counts for autosomes}
\description{
  Performs optimization of joint TReC and ASE for autosome and tests with lrt test for two hypotheses hypotheses:
  additive and parent of origin. This function modification assumes that phiNB and phiBB are common for the gene of interest. 
  Otherwise the model is exactly the same.
}
\usage{
fit(subset=NULL, data, traceit=FALSE)}
\arguments{
  \item{subset}{a subset of entries to be tested, but default is set to NULL which leads to fitting all the genes in the table}
  \item{data}{an object of class readCounts including read counts and necessary supporting information}
  \item{traceit}{include more debug output, by default set to FALSE}
}
\value{
  a list of following matrices:
  \item{full}{matrix with columns: log(phiNB) and -log(phiBB), coefficients, -log(likelihood) for the full model fit as well as appended two p-value tests for additive and parent of origin effect}  
  \item{testadd}{matrix with columns: log(phiNB) and -log(phiBB), coefficients, -log(likelihood) for the restricted to b0=0 model fit}
  \item{testpoo}{matrix with columns: log(phiNB) and -log(phiBB), coefficients, -log(likelihood) for the restricted to b1=0 model fit}
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{fitsh}}, \code{\link{data}}, \code{\link{simu}}, \code{\link{readCounts}}.}
\keyword{ methods }

\examples{\dontrun{
# fitting autosome data for a full model with allele-specific counts collected on gene level:
percase = 0.1
dblcnt = 0.2
mn = 100
b0 = 0
b1 = 0
phiNB = .5;
phiBB=phiNB/4
niter = 100
betas = c(3,.2,.05,.5)
ss=2

dep = makeXmatr(ss)
dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
commest = fitsh(subset=NULL, data=dat, traceit=FALSE)
mean(commest$full[,"pval_b0"]<.05)
}}
