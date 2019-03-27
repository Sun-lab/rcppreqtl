readCounts = function(haplotype, trc, asn, asnp, 
                      haplotype2S=NULL, asn2S=NULL, asnp2S=NULL,
                      haplotype4S=NULL, asn4S=NULL, asnp4S=NULL,
                      haplotype8S=NULL, asn8S=NULL, asnp8S=NULL,
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
                haplotype2S = haplotype2S,
                asn2S = asn2S,
                asnp2S = asn2S,
                haplotype4S = haplotype4S,
                asn4S = asn4S,
                asnp4S = asn4S,
                haplotype8S = haplotype8S,
                asn8S = asn8S,
                asnp8S = asn8S,
                X = X,
                params=params,
                settings=settings
       )

        ## Set the name for the class
        class(rc) = append(class(rc),"ReadCounts")
        return(rc)
}