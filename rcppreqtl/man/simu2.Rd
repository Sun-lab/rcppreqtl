\name{simu2}
\alias{simu2}
\title{Simulate a dataset in a format acceptable by a fit function}
\description{
  Creates an object of a class readCounts with data simulated according to a provided setup
}
\usage{
  simu2(num, Xmatr, haplotype, totmean, percase=0.1, dblcnt=0, phiNB=1, phiBB=0.5, b0=0, b1=0, betas=rep(1,4))
fullest = fit(subset=NULL, data=dat, traceit=FALSE)
}
\arguments{
  \item{num}{number of iterations}
  \item{Xmatr}{design matrix for total read counts covariates}
  \item{haplotype}{classes of haplotypes 0 - AA, 1 - AB, 2 - BA, 3 - BB where first letter represents paternal allele}
  \item{totmean}{average total gene expression}
  \item{percase}{percentage of reads classified as allele-specific, default value 10\%}
  \item{dblcnt}{optional output considering a simulation split into 2 SNPs with double-counting. Default value 0.}
  \item{phiNB}{over-dispersion for Negative-Binomial distribution, default value 1}
  \item{phiBB}{over-dispersion for Beta-Binomial distribution, default value 0.5}
  \item{b0}{additive eQTL, default value 0}
  \item{b1}{parent of origin effect, default value 0}
  \item{betas}{covariates for design matrix Xmatr}

}
\value{
  \item{an object of class readCounts}{simulated dataset that can be used to fit the model}  
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{fitsh}}, \code{\link{data}}, \code{\link{simu4}}, \code{\link{simu2}}, \code{\link{readCounts}}.}
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
dat = simu2(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
}}
