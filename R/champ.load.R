champ.load <-
function(directory = getwd(), methValue="B", resultsDir=paste(getwd(), "resultsChamp",sep="/"), filterXY=TRUE, QCimages=TRUE, filter=TRUE, detPcut=0.01)
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
    suppressWarnings(targets <- read.450k.sheet(myDir))
	rgSet <- read.450k.exp(base = myDir, targets = targets)
	sampleNames(rgSet)=rgSet[[1]]
	pd<-pData(rgSet)
	green=getGreen(rgSet)
	red=getRed(rgSet)
	mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
	
    if(filter)
    {
        failed <- detP>detPcut
        numfail = colMeans(failed) # Fraction of failed positions per sample
        message("The fraction of failed positions per sample: ")
        print(numfail)
        fileName=paste(resultsDir,"/failedSample",".txt",sep="")
        write.table(numfail, fileName, row.names=T, col.names=paste("Sample_Name","Fraction_Failed_Probes",sep="\t"), quote=F,sep="\t")
        mset.f = mset[rowSums(detP <= detPcut) == ncol(detP),]
        message("By filtering with a detection p-value of ",detPcut," a total of ",dim(mset)[1]-dim(mset.f)[1]," probes have been removed from the analysis.")
        mset=mset.f
    }
        
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
		data(probe.features)
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
		beta.raw=beta.raw[row.names(beta.raw) %in% row.names(autosomes), ]
        detP=detP[row.names(detP) %in% row.names(autosomes), ]
        intensity=intensity[row.names(intensity) %in% row.names(autosomes), ]
	}
    
    totalProbes=dim(beta.raw)[1]

	#cluster
	if(QCimages)
	{
		if(ncol(beta.raw) > 0)
		{
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
