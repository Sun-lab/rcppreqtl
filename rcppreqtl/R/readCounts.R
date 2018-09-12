readCounts = function(haplotype, trc, asn, asnp, haplotypeA=NULL, asnA=NULL, asnpA=NULL,
             X, params=NULL, settings=NULL)
{
  if(ncol(trc) != ncol(haplotype)){
    stop("number of columns of y should match length of kappas")
  }
        rc = list(
                haplotype=haplotype,
                trc = trc,
                asn = asn,
                asnp = asnp,
                haplotypeA = haplotypeA,
                asnA = asnA,
                asnpA = asnpA,
                X = X,
                params=params,
                settings=settings
       )

        ## Set the name for the class
        class(rc) = append(class(rc),"ReadCounts")
        return(rc)
}