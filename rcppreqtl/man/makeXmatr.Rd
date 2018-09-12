\name{makeXmatr}
\alias{makeXmatr}
\title{Create example design matrix for simulations}
\description{
  Produces a design matrix of several sample sizes which can be used to generate simulated dataset
}
\usage{
makeXmatr(ss)
}
\arguments{
  \item{ss}{Sample size class: ss=1 implies dample size 32, ss=2 implies sample size 64, etc}
}
\value{
  a design matrix of 4 variables
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
}}
