if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","snp.hit","multi.hit","probe.features"))

champ.load <- function(directory = getwd(),
                       methValue="B",
                       resultsDir=paste(getwd(), "resultsChamp",sep="/"),
                       filterXY=TRUE,
                       QCimages=TRUE,
                       filterDetP=TRUE,
                       detPcut=0.01,
                       removeDetP = 0,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=FALSE,
                       filterSNPs=TRUE,
                       filterMultiHit=TRUE,
                       arraytype="450K")
{
	read.450k.sheet<-NA
	rm(read.450k.sheet)
	read.450k.exp<-NA
	rm(read.450k.exp)
	pData<-NA
	rm(pData)
	getGreen<-NA
	rm(getGreen)
	getRed<-NA
	rm(getRed)
	preprocessRaw<-NA
	rm(preprocessRaw)	
	getMeth<-NA
	rm(getMeth)
	getUnmeth<-NA
	rm(getUnmeth)		
	getBeta<-NA
	rm(getBeta)
	getM<-NA
	rm(getM)
	plotBetasByType<-NA
	rm(plotBetasByType)
	mdsPlot<-NA
	rm(mdsPlot)
    detectionP<-NA
    rm(detectionP)

	
	message("Loading data from ", directory)

	if(!file.exists(resultsDir))
	{
		dir.create(resultsDir)
		message("Creating results directory. Results will be saved in ", resultsDir)
	}			
	
    myDir= directory
    suppressWarnings(targets <- read.metharray.sheet(myDir))
	rgSet <- read.metharray.exp(targets = targets, extended=TRUE)
    if(arraytype=="EPIC")
        rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilmn10.hg19")
	sampleNames(rgSet)=rgSet[[1]]
	pd<-pData(rgSet)
	green=getGreen(rgSet)
	red=getRed(rgSet)
	mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
	
    message("Read DataSet Success.\n")

    failed <- detP>detPcut
    numfail = colMeans(failed) # Fraction of failed positions per sample
    message("The fraction of failed positions per sample: ")
    print(numfail)
    fileName=paste(resultsDir,"/failedSample",".txt",sep="")
    write.table(numfail, fileName, row.names=T, col.names=paste("Sample_Name","Fraction_Failed_Probes",sep="\t"), quote=F,sep="\t")
    
    if(filterDetP)
    {
        mset.f = mset[rowSums(detP >= detPcut) <= removeDetP*ncol(detP),]
        
        if(removeDetP==0)
        {
            message("Filtering probes with a detection p-value above ",detPcut," in one or more samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.")
        }else{
            message("Filtering probes with a detection p-value above ",detPcut," in at least ",removeDetP*100,"% of samples has removed ",dim(mset)[1]-dim(mset.f)[1]," probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.")
        }

        mset=mset.f
    }
    
    message("Filter DetP Done.\n")

    if(filterBeads)
    {
        bc=beadcount(rgSet)
        bc2=bc[rowSums(is.na(bc))< beadCutoff*(ncol(bc)),]
        mset.f2 = mset[featureNames(mset) %in% row.names(bc2),]
        message("Filtering probes with a beadcount <3 in at least ",beadCutoff*100,"% of samples, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }

    message("Filter Beads Done.\n")
    
    if(filterNoCG)
    {
        mset.f2=mset
        dropMethylationLoci(mset,dropCH=T)
        message("Filtering non-cg probes, has removed ",dim(mset.f2)[1]-dim(mset)[1]," from the analysis.")
    }
    
    message("Filter NoCG Done.\n")

    if(filterSNPs)
    {
        data(snp.hit)
        mset.f2=mset[!featureNames(mset) %in% snp.hit$TargetID,]
        message("Filtering probes with SNPs as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }

    message("Filter SNP Done.\n")
    
    if(filterMultiHit)
    {
        data(multi.hit)
        mset.f2=mset[!featureNames(mset) %in% multi.hit$TargetID,]
        message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
    }
    
    message("Filter MultiHit Done.\n")

    intensity=getMeth(mset)+getUnmeth(mset)
    
    if(methValue=="B")
    {
        beta.raw = getBeta(mset, "Illumina")
    }else{
        beta.raw = getM(mset)
        }
        
    head(beta.raw)    
    
    detP=detP[which(row.names(detP) %in% row.names(beta.raw)),]
    
    if(filterXY)
	{
        if(arraytype=="EPIC")
        {
            data(probe.features.epic)
        }else{
		data(probe.features)
        }
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
        mset.f2=mset[featureNames(mset) %in% row.names(autosomes),]
		beta.raw=beta.raw[row.names(beta.raw) %in% row.names(autosomes), ]
        detP=detP[row.names(detP) %in% row.names(autosomes), ]
        intensity=intensity[row.names(intensity) %in% row.names(autosomes), ]
        message("Filtering probes on the X or Y chromosome has removed ",dim(mset)[1]-dim(mset.f2)[1]," from the analysis.")
        mset=mset.f2
	}

    message("Filter XY chromosome Done.\n")
    
    totalProbes=dim(beta.raw)[1]

	#cluster
	if(QCimages)
	{
		if(ncol(beta.raw) > 0)
		{
            if(min(beta.raw, na.rm=TRUE)==0)
            {
                message("Zeros in your dataset have been replaced with 0.000001")
                beta.raw[beta.raw==0]<-0.000001
            }
            
		#save images...
        imageName1=paste(resultsDir,"raw_mdsPlot.pdf",sep="/")
        imageName2=paste(resultsDir,"raw_densityPlot.pdf",sep="/")
        imageName3=paste(resultsDir,"raw_SampleCluster.jpg",sep="/")
        main1=paste("Density plot of raw data (",totalProbes," probes)",sep="")
        main2=paste("All samples before normalization (",totalProbes, " probes)",sep="")
            
		pdf(imageName1,width=6,height=4)
		mdsPlot(beta.raw,numPositions=1000,sampGroups=pd$Sample_Group,sampNames=pd$Sample_Name)
		dev.off()
				
		pdf(imageName2,width=6,height=4)
		densityPlot(rgSet,sampGroups=pd$Sample_Group,main=main1,xlab="Beta")
		dev.off()
		
		if(ncol(beta.raw) < 60)
		{	
			jpeg(imageName3)
			betar_d<-dist(t(beta.raw))
			plot(hclust(betar_d),main=main2,cex=0.8)
			dev.off()
		}else{
			message("Cluster image is not saved when the number of samples exceeds 60.")
			}
		}
	}
    message("The analysis will proceed with ", dim(beta.raw)[1], " probes and ",dim(beta.raw)[2], " samples.")
	return(list(mset=mset,rgSet=rgSet,pd=pd,intensity=intensity,beta=beta.raw,detP=detP));
}
