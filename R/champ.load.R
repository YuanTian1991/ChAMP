if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","snp.hit","multi.hit","probe.features","probe.features.epic"))

champ.load <- function(directory = getwd(),
                       methValue="B",
                       filterDetP=TRUE,
                       detPcut=0.01,
                       removeDetP = 0,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=TRUE,
                       filterSNPs=TRUE,
                       filterMultiHit=TRUE,
                       filterXY=TRUE,
                       arraytype="450K")
{
    message("[===========================]")
    message("[<<<< ChAMP.LOAD START >>>>>]")
    message("-----------------------------")
	message("Loading data from ", directory)
	
    myDir <- directory
    suppressWarnings(targets <- read.metharray.sheet(myDir))
    rgSet <- read.metharray.exp(targets = targets,extended=TRUE)
    if(arraytype=="EPIC") rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilmn10.hg19")

	sampleNames(rgSet)=rgSet[[1]]
	pd <- pData(rgSet)
	mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
	
    message("<< Read DataSet Success. >>\n")

    message("The fraction of failed positions per sample: ")
    numfail <- matrix(colMeans(detP>detPcut))
    rownames(numfail) <- colnames(detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    
    if(filterDetP)
    {
        mset.f = mset[rowSums(detP >= detPcut) <= removeDetP*ncol(detP),]
        
        if(removeDetP==0)
        {
            message("Filtering probes with a detection p-value above ",detPcut," in one or more samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.")
        }else{
            message("Filtering probes with a detection p-value above ",detPcut," in at least ",removeDetP*100,"% of samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.")
        }
        mset=mset.f
    }
    
    message("<< Filter DetP Done. >>\n")

    if(filterBeads)
    {
        bc=beadcount(rgSet)
        bc2=bc[rowSums(is.na(bc))< beadCutoff*(ncol(bc)),]
        mset.f2 = mset[featureNames(mset) %in% row.names(bc2),]
        message("Filtering probes with a beadcount <3 in at least ",beadCutoff*100,"% of samples, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    message("<< Filter Beads Done. >>\n")

    if(filterNoCG)
    {
        mset.f2=mset
        dropMethylationLoci(mset,dropCH=T)
        message("Filtering non-cg probes, has removed ",dim(mset.f2)[1]-dim(mset)[1]," from the analysis.")
    }
    message("<< Filter NoCG Done. >>\n")

    if(filterSNPs)
    {
        data(snp.hit)
        mset.f2=mset[!featureNames(mset) %in% snp.hit$TargetID,]
        message("Filtering probes with SNPs as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    message("<< Filter SNP Done. >>\n")
    
    if(filterMultiHit)
    {
        data(multi.hit)
        mset.f2=mset[!featureNames(mset) %in% multi.hit$TargetID,]
        message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    message("<< Filter MultiHit Done. >>\n")

    if(filterXY)
	{
        if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
        mset.f2=mset[featureNames(mset) %in% row.names(autosomes),]
        message("Filtering probes on the X or Y chromosome has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
	}
    message("<< Filter XY chromosome Done. >>\n")

    
    if(methValue=="B")beta.raw = getBeta(mset, "Illumina") else beta.raw = getM(mset)
    message(paste(if(methValue=="B") "[Beta" else "[M","value is selected as output.]\n"))

    intensity <-  minfi::getMeth(mset) + minfi::getUnmeth(mset)
    detP <- detP[which(row.names(detP) %in% row.names(beta.raw)),]

    if(min(beta.raw, na.rm=TRUE)==0) beta.raw[beta.raw==0] <- 0.000001
    message("Zeros in your dataset have been replaced with 0.000001\n")

    message("The analysis will be proceed with ", dim(beta.raw)[1], " probes and ",dim(beta.raw)[2], " samples.\n")
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
	return(list(mset=mset,rgSet=rgSet,pd=pd,intensity=intensity,beta=beta.raw,detP=detP))
}
