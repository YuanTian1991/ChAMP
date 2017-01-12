if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","EPIC.manifest.hg38","EPIC.manifest.pop.hg38","hm450.manifest.hg38","hm450.manifest.pop.hg38","multi.hit","probe.features","probe.features.epic"))


champ.filter <- function(beta,
                         detP=NULL,
                         pd,
                         filterDetP=TRUE,
                         detSamplecut=0.1,
                         detPcut=0.01,
                         removeDetP=0,
                         filterNoCG = TRUE,
                         filterSNPs = TRUE, 
                         population = NULL,
                         filterMultiHit = TRUE,
                         filterXY = TRUE, 
                         arraytype = "450K")
{
    message("[===========================]")
    message("[<<<< ChAMP.FILTER START >>>>>]")
    message("-----------------------------")


    message("This function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix, you also need to input a detected P matrix.")

    if(is.null(beta) || is.null(pd) || !class(beta) %in% c("matrix","data.frame"))
        stop("One beta matrix and one pd file must be provided.")

    mset <- beta

    if(filterDetP && is.null(detP))
    {
        message("If you want to do filtering on detected p value for probes, detP Matrix MUST be provided. filterDetP parameter will be set FALSE.")
        filterDetP <- FALSE
    }

    if(filterDetP)
    {
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

        detP <- detP[,RemainSample]
        mset <- mset[,RemainSample]
        pd <- pd[RemainSample,]

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

    if(filterNoCG)
    {
        mset.f2=mset[which(substr(rownames(mset),1,2)=="cg"),]
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
                message("Seems your population name is wrong. Using general 450K SNP list for filtering.")
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
                message("Seems your population name is wrong. Using general EPIC SNP list for filtering.")
                data(EPIC.manifest.hg38)
                maskname <- rownames(EPIC.manifest.hg38)[which(EPIC.manifest.hg38$MASK.general==TRUE)]
            }else
            {
                message("Using ",population," specific EPIC SNP list for filtering.")
                data(EPIC.manifest.pop.hg38)
                maskname <- rownames(EPIC.manifest.pop.hg38)[which(EPIC.manifest.pop.hg38[,paste("MASK.general",population,sep=".")]==TRUE)]
            }
        }
        mset.f2=mset[!rownames(mset) %in% maskname,]
        message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    message("<< Filter SNP Done. >>\n")
    
    if(filterMultiHit)
    {
        data(multi.hit)
        mset.f2=mset[!rownames(mset) %in% multi.hit$TargetID,]
        message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    message("<< Filter MultiHit Done. >>\n")

    if(filterXY)
	{
        if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
        mset.f2=mset[rownames(mset) %in% row.names(autosomes),]
        message("Filtering probes on the X or Y chromosome has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
	}
    message("<< Filter XY chromosome Done. >>\n")

    
    if(min(mset, na.rm=TRUE)==0) mset[mset==0] <- 0.000001
    message("Zeros in your dataset have been replaced with 0.000001\n")

    message("The analysis will be proceed with ", dim(mset)[1], " probes and ",dim(mset)[2], " samples.\n")
    message("[<<<<< ChAMP.FILTER END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
	return(list(beta=mset,pd=pd))
}

