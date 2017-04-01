if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","EPIC.manifest.hg38","EPIC.manifest.pop.hg38","hm450.manifest.hg38","hm450.manifest.pop.hg38","multi.hit","probe.features","probe.features.epic"))

champ.load <- function(directory = getwd(),
                       methValue="B",
                       filterDetP=TRUE,
                       detSamplecut=0.1,
                       detPcut=0.01,
                       removeDetP=0,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=TRUE,
                       filterSNPs=TRUE,
                       population=NULL,
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

    message("The fraction of failed positions per sample\n 
            (You may need to delete samples with high proportion of failed probes\n): ")
    numfail <- matrix(colMeans(detP>detPcut))
    rownames(numfail) <- colnames(detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < detSamplecut)

    if(any(numfail >= detSamplecut))
    message("The detSamplecut parameter is : ",detSamplecut, "\nSamples : ",
            paste(rownames(numfail)[which(numfail >= detSamplecut)],collapse=",")," will be deleted.\n",
            "There are ",length(RemainSample)," samples left for analysis.\n")
    
    rgSet <- rgSet[,RemainSample]
    detP <- detP[,RemainSample]
    mset <- mset[,RemainSample]
    pd <- pd[RemainSample,]

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
        mset.f2=dropMethylationLoci(mset,dropCH=T)
        message("Filtering non-cg probes, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset <- mset.f2
    }
    message("<< Filter NoCG Done. >>\n")

    if(filterSNPs)
    {
        if(arraytype=="450K")
        {
            if(is.null(population))
            {
                message("Using general 450K SNP list for filtering.")
                data(hm450.manifest.hg38)
                maskname <- rownames(hm450.manifest.hg38)[which(hm450.manifest.hg38$MASK.general==TRUE)]
            }else if(!population %in% c("AFR","EAS","EUR","SAS","AMR","GWD","YRI","TSI",
                                        "IBS","CHS","PUR","JPT","GIH","CHB","STU","ITU",
                                        "LWK","KHV","FIN","ESN","CEU","PJL","ACB","CLM",
                                        "CDX","GBR","BEB","PEL","MSL","MXL","ASW"))
            {
                message("Using general 450K SNP list for filtering.")
                data(hm450.manifest.hg38)
                maskname <- rownames(hm450.manifest.hg38)[which(hm450.manifest.hg38$MASK.general==TRUE)]
            }else
            {
                message("Using ",population," specific 450K SNP list for filtering.")
                data(hm450.manifest.pop.hg38)
                maskname <- rownames(hm450.manifest.pop.hg38)[which(hm450.manifest.pop.hg38[,paste("MASK.general",population,sep=".")]==TRUE)]
            }
        }else
        {
            if(is.null(population))
            {
                message("Using general EPIC SNP list for filtering.")
                data(EPIC.manifest.hg38)
                maskname <- rownames(EPIC.manifest.hg38)[which(EPIC.manifest.hg38$MASK.general==TRUE)]
            }else if(!population %in% c("AFR","EAS","EUR","SAS","AMR","GWD","YRI","TSI",
                                        "IBS","CHS","PUR","JPT","GIH","CHB","STU","ITU",
                                        "LWK","KHV","FIN","ESN","CEU","PJL","ACB","CLM",
                                        "CDX","GBR","BEB","PEL","MSL","MXL","ASW"))
            {
                message("Using general EPIC SNP list for filtering.")
                data(EPIC.manifest.hg38)
                maskname <- rownames(EPIC.manifest.hg38)[which(EPIC.manifest.hg38$MASK.general==TRUE)]
            }else
            {
                message("Using ",population," specific EPIC SNP list for filtering.")
                data(EPIC.manifest.pop.hg38)
                maskname <- rownames(EPIC.manifest.pop.hg38)[which(EPIC.manifest.pop.hg38[,paste("MASK.general",population,sep=".")]==TRUE)]
            }
        }
        mset.f2=mset[!featureNames(mset) %in% maskname,]
        message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
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

    if(min(beta.raw, na.rm=TRUE)<=0) beta.raw[beta.raw<=0] <- min(beta.raw[beta.raw > 0])
    message("Zeros in your dataset have been replaced with smallest positive value.\n")

    if(max(beta.raw, na.rm=TRUE)>=0) beta.raw[beta.raw>=1] <- max(beta.raw[beta.raw < 1])
    message("One in your dataset have been replaced with largest value below 1.\n")

    message("The analysis will be proceed with ", dim(beta.raw)[1], " probes and ",dim(beta.raw)[2], " samples.\n")
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
	return(list(mset=mset,rgSet=rgSet,pd=pd,intensity=intensity,beta=beta.raw,detP=detP))
}
