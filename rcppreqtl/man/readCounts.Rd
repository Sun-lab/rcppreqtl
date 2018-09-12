\name{readCounts}
\alias{readCounts}
\title{A list object that should be used as input to optimization fit or fitsh function.}
\description{
  It should contain at least total read counts (TReC), overall allele-specific counts (ASE), paternal allele-specific counts and haplotype classification 0 to 4. Also should include X matrix with covariates such as intercept, library depth, principal components.
}
\value{  
  \item{haplotype}{a matrix defining the haplotype status of the individual for each gene - AB=0, BA=1, AA=2, BB=3 coded as PaternalMaternal haplotypes. Each row - one gene}
  \item{trc}{matrix of TReC counts. Each row - one gene}
  \item{asn}{matrix of ASE counts for subject (column) for corresponding genes (row).}  
  \item{asnp}{matrix of ASE counts belonging to paternal allele, for correponding subjects and genes as in \bold{asn}.}  
  \item{haplotypeA}{haplotype modification for a setup with allele-specific reads collected for multiple SNPs in a gene. This block by default is NULL, but can be used by simulation function to study scenarios when reads are collected not on gene level, but on SNP level}
  \item{asnA}{allele-specific count corresponding to \bold{haplotypeA}}
  \item{asnpA}{paternal allele-specific count corresponding to \bold{haplotypeA}}
  \item{X}{design matrix for total read counts - a place to include intercept, library depth, other covariates such as batch effect}
  \item{params}{if data is produced by simulation function parameter values can be stored here for further comparisons}
  \item{settings}{other important settings for the simulation such as mean of total read counts, percentage of allele specific counts, percentage of reads double-counted with neighboring SNPs can be stored here}
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\keyword{ utilities }
\seealso{\code{\link{fit}}, \code{\link{fitsh}}, \code{\link{simu2}}, \code{\link{simu4}}.}

\examples{\dontrun{
# see total read counts (TReC) for first 2 X chromosome genes of a data example:
rc = readCounts(haplotype=haplotypef, trc=trc, asn=asnf, asnp=asnpf, haplotypeA=haplotyped, asnA=asnm, asnpA=asnpm,
             X=Xmatr, params=c(phiNB, phiBB, b0, b1, betas), settings=c(totmean, percase, dblcnt))
}}