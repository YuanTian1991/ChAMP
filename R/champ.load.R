if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","EPIC.manifest.hg19","EPIC.manifest.pop.hg19","hm450.manifest.hg19","hm450.manifest.pop.hg19","multi.hit","probe.features","probe.features.epic"))

champ.load <- function(directory = getwd(),
                       method = "ChAMP",
                       methValue="B",
                       autoimpute=TRUE,
                       filterDetP=TRUE,
                       ProbeCutoff=0,
                       SampleCutoff=0.1,
                       detPcut=0.01,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=TRUE,
                       filterSNPs=TRUE,
                       population=NULL,
                       filterMultiHit=TRUE,
                       filterXY=TRUE,
                       force=FALSE,
                       arraytype="450K")
{
    message("[===========================]")
    message("[<<<< ChAMP.LOAD START >>>>>]")
    message("-----------------------------")

    mybeadcount <- function(x)
    {
        #select out bead count dataframe
        getNBeads(x) -> nb
        #match rownames of beadcount dataframe to addresses
        getProbeInfo(x,type="I")->typeIadd
        match(typeIadd$AddressA,rownames(nb))->typeImatchA
        match(typeIadd$AddressB,rownames(nb))->typeImatchB

        #match rownames of beadcount dataframe to addresses
        getProbeInfo(x,type="II")->typeIIadd
        match(typeIIadd$Address,rownames(nb))->typeIImatch

        nb->nbcg

        locusNames <- getManifestInfo(x, "locusNames")
        bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
        dimnames = list(locusNames, sampleNames(x)))
    
        TypeII.Name <- getProbeInfo(x, type = "II")$Name
        bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,]
       
        TypeI <- getProbeInfo(x, type = "I")

        bc_temp->bcB
        bc_temp->bcA        
    
        bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB,]
        bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA,]

        which(bcB<3) -> bcB3
        which(bcA<3) -> bcA3
        bcA->bcA2
        bcB->bcB2
        bcA2[bcA3]<-NA
        bcA2[bcB3]<-NA
    
        data.frame(bcA2)->bc
        bc
    }

    if(method=="minfi")
    {
        message("\n[ Loading Data with Minfi Method ]")
        message("----------------------------------")
    
	message("Loading data from ", directory)
	
    myDir <- directory
    suppressWarnings(targets <- read.metharray.sheet(myDir))
    rgSet <- read.metharray.exp(targets = targets,extended=TRUE,force=force)
    if(arraytype=="EPIC") rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b3.hg19")

	sampleNames(rgSet)=rgSet[[1]]
	pd <- pData(rgSet)
	mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
    message("<< Read DataSet Success. >>\n")

    if(methValue=="B") tmp = getBeta(mset, "Illumina") else tmp = getM(mset)
    tmp[detP >= detPcut] <- NA 
    message("The fraction of failed positions per sample\n 
            (You may need to delete samples with high proportion of failed probes\n): ")

    numfail <- matrix(colMeans(is.na(tmp)))
    rownames(numfail) <- colnames(detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < SampleCutoff)

    if(any(numfail >= SampleCutoff))
    message("The detSamplecut parameter is : ",SampleCutoff, "\nSamples : ",
            paste(rownames(numfail)[which(numfail >= SampleCutoff)],collapse=",")," will be deleted.\n",
            "There are ",length(RemainSample)," samples left for analysis.\n")
    
    rgSet <- rgSet[,RemainSample]
    detP <- detP[,RemainSample]
    mset <- mset[,RemainSample]
    pd <- pd[RemainSample,]
    tmp <- tmp[,RemainSample]


    if(filterDetP)
    {
        mset.f = mset[rowSums(is.na(tmp)) <= ProbeCutoff*ncol(detP),]
        
        if(ProbeCutoff==0)
        {
            message("Filtering probes with a detection p-value above ",detPcut," in one or more samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.")
        }else{
            message("Filtering probes with a detection p-value above ",detPcut," in at least ",ProbeCutoff*100,"% of samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample file to identify potentially bad samples.")
        }
        mset=mset.f
        tmp <- tmp[rowSums(is.na(tmp)) <= ProbeCutoff*ncol(detP),]
        message("<< Filter DetP Done. >>\n")
    }

    if(sum(is.na(tmp))==0){
       message("\nThere is no NA values in your matrix, there is no need to imputation.\n")
    }else
    {
        message("\nThere are ",sum(is.na(tmp))," NA remain in filtered Data Set. Impute can be done for remain NAs, but not suitable for small number samples. For small Data Set (like only 20 samples), we suggest you set parameter ProbeCutoff as 0 in champ.load() here, which would remove all NA involved probe no matter how many samples of those probes are NA.\n")
    }

    if(autoimpute & sum(is.na(tmp)) > 0){
        message("Impute will be conducted here for remain ",sum(is.na(tmp)),"  NAs. Note that if you don't do this, NA values would be kept in your data set. You may use champ.impute() function to do more complex imputation as well.")
        # Open a file to send messages to save lot's of information from impute.knn
        message("\nImpute function is working now, it may need couple minutes...")
        zz <- file("ImputeMessage.Rout", open="wt")
        sink(zz)
        sink(zz, type="message")
        tmp <- impute.knn(tmp,k=5)$data
        sink(type="message")
        sink()
        message("<< Imputation Done. >>\n")
    }
    

    if(filterBeads)
    {
        bc=mybeadcount(rgSet)
        bc2=bc[rowSums(is.na(bc)) < beadCutoff*(ncol(bc)),]

        mset.f2 = mset[featureNames(mset) %in% row.names(bc2),]
        tmp <- tmp[rownames(tmp) %in% row.names(bc2),]
        message("Filtering probes with a beadcount <3 in at least ",beadCutoff*100,"% of samples, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
        message("<< Filter Beads Done. >>\n")
    }

    if(filterNoCG)
    {
        mset.f2=dropMethylationLoci(mset,dropCH=T)
        tmp <- tmp[rownames(tmp) %in% featureNames(mset.f2),]
        message("Filtering non-cg probes, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset <- mset.f2
        message("<< Filter NoCG Done. >>\n")
    }

    if(filterSNPs)
    {
        if(arraytype=="450K")
        {
            if(is.null(population))
            {
                message("Using general 450K SNP list for filtering.")
                data(hm450.manifest.hg19)
                maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general==TRUE)]
            }else if(!population %in% c("AFR","EAS","EUR","SAS","AMR","GWD","YRI","TSI",
                                        "IBS","CHS","PUR","JPT","GIH","CHB","STU","ITU",
                                        "LWK","KHV","FIN","ESN","CEU","PJL","ACB","CLM",
                                        "CDX","GBR","BEB","PEL","MSL","MXL","ASW"))
            {
                message("Using general 450K SNP list for filtering.")
                data(hm450.manifest.hg19)
                maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general==TRUE)]
            }else
            {
                message("Using ",population," specific 450K SNP list for filtering.")
                data(hm450.manifest.pop.hg19)
                maskname <- rownames(hm450.manifest.pop.hg19)[which(hm450.manifest.pop.hg19[,paste("MASK_general_",population,sep="")]==TRUE)]
            }
        }else
        {
            if(is.null(population))
            {
                message("Using general EPIC SNP list for filtering.")
                data(EPIC.manifest.hg19)
                maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general==TRUE)]
            }else if(!population %in% c("AFR","EAS","EUR","SAS","AMR","GWD","YRI","TSI",
                                        "IBS","CHS","PUR","JPT","GIH","CHB","STU","ITU",
                                        "LWK","KHV","FIN","ESN","CEU","PJL","ACB","CLM",
                                        "CDX","GBR","BEB","PEL","MSL","MXL","ASW"))
            {
                message("Using general EPIC SNP list for filtering.")
                data(EPIC.manifest.hg19)
                maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general==TRUE)]
            }else
            {
                message("Using ",population," specific EPIC SNP list for filtering.")
                data(EPIC.manifest.pop.hg19)
                maskname <- rownames(EPIC.manifest.pop.hg19)[which(EPIC.manifest.pop.hg19[,paste("MASK_general_",population,sep="")]==TRUE)]
            }
        }
        mset.f2=mset[!featureNames(mset) %in% maskname,]
        tmp <- tmp[!rownames(tmp) %in% maskname,]
        message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
        message("<< Filter SNP Done. >>\n")
    }
    
    if(filterMultiHit)
    {
        data(multi.hit)
        mset.f2=mset[!featureNames(mset) %in% multi.hit$TargetID,]
        tmp <- tmp[!rownames(tmp) %in% multi.hit$TargetID,]
        message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
        message("<< Filter MultiHit Done. >>\n")
    }

    if(filterXY)
	{
        if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
        mset.f2=mset[featureNames(mset) %in% row.names(autosomes),]
        tmp <- tmp[rownames(tmp) %in% row.names(autosomes),]
        message("Filtering probes on the X or Y chromosome has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
        message("<< Filter XY chromosome Done. >>\n")
	}

    
    message(paste(if(methValue=="B") "[Beta" else "[M","value is selected as output.]\n"))
    beta.raw <- tmp

    intensity <-  minfi::getMeth(mset) + minfi::getUnmeth(mset)
    detP <- detP[which(row.names(detP) %in% row.names(beta.raw)),]

    if(min(beta.raw, na.rm=TRUE)<=0) beta.raw[beta.raw<=0] <- min(beta.raw[beta.raw > 0])
    message("Zeros in your dataset have been replaced with smallest positive value.\n")

    if(max(beta.raw, na.rm=TRUE)>=0) beta.raw[beta.raw>=1] <- max(beta.raw[beta.raw < 1])
    message("One in your dataset have been replaced with largest value below 1.\n")

    message("The analysis will be proceed with ", dim(beta.raw)[1], " probes and ",dim(beta.raw)[2], " samples.\n")
    message("Current Data Set contains ",sum(is.na(beta.raw))," NA in ", if(methValue=="B") "[Beta]" else "[M]"," Matrix.\n")

    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
	return(list(mset=mset,rgSet=rgSet,pd=pd,intensity=intensity,beta=beta.raw,detP=detP))

    } else {
        message("\n[ Loading Data with ChAMP Method ]")
        message("----------------------------------")

        message("Note that ChAMP method will NOT return rgSet or mset, they object defined by minfi. Which means, if you use ChAMP method to load data, you can not use SWAN or FunctionNormliazation method in champ.norm() (you can use BMIQ or PBC still). But All other function should not be influenced.\n")


        myImport <- champ.import(directory,arraytype=arraytype)
        
        if(methValue=="B")
            myLoad <- champ.filter(beta=myImport$beta,
                                   M=NULL,
                                   pd=myImport$pd,
                                   intensity=myImport$intensity,
                                   Meth=NULL,
                                   UnMeth=NULL,
                                   detP=myImport$detP,
                                   beadcount=myImport$beadcount,
                                   autoimpute=autoimpute,
                                   filterDetP=filterDetP,
                                   ProbeCutoff=ProbeCutoff,
                                   SampleCutoff=SampleCutoff,
                                   detPcut=detPcut,
                                   filterBeads=filterBeads,
                                   beadCutoff=beadCutoff,
                                   filterNoCG=filterNoCG,
                                   filterSNPs=filterSNPs,
                                   population=population,
                                   filterMultiHit=filterMultiHit,
                                   filterXY=filterXY,
                                   arraytype=arraytype)
        else
            myLoad <- champ.filter(beta=NULL,
                                   M=myImport$M,
                                   pd=myImport$pd,
                                   intensity=myImport$intensity,
                                   Meth=NULL,
                                   UnMeth=NULL,
                                   detP=myImport$detP,
                                   beadcount=myImport$beadcount,
                                   autoimpute=autoimpute,
                                   filterDetP=filterDetP,
                                   ProbeCutoff=ProbeCutoff,
                                   SampleCutoff=SampleCutoff,
                                   detPcut=detPcut,
                                   filterBeads=filterBeads,
                                   beadCutoff=beadCutoff,
                                   filterNoCG=filterNoCG,
                                   filterSNPs=filterSNPs,
                                   population=population,
                                   filterMultiHit=filterMultiHit,
                                   filterXY=filterXY,
                                   arraytype=arraytype)

    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
	return(myLoad)
    }
}
