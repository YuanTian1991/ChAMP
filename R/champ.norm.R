champ.norm <-
function(beta=myLoad$beta, rgSet=myLoad$rgSet, pd=myLoad$pd,mset=myLoad$mset, sampleSheet="sampleSheet.txt", resultsDir=paste(getwd(), "resultsChamp", sep="/"), methValue="B", fromIDAT=TRUE, norm="BMIQ", fromFile = FALSE, betaFile, filter=TRUE, filterXY=TRUE, QCimages=TRUE, plotBMIQ=TRUE)
{
	detectionP<-NA
	rm(detectionP)
	preprocessSWAN<-NA
	rm(preprocessSWAN)
	getBeta<-NA
	rm(getBeta)	
	getM<-NA
	rm(getM)
	cwd=getwd()
	data(probe.features)	

	message("Normalizing data with ",norm)
	if(fromIDAT==F)
	{
		if(norm=="SWAN")
		{
			message("Swan normalization is not available from this data. You need to use champ.load() with raw IDAT files.")
			return()
			
		}else if(fromFile==F)
		{
			beta.p=betaFile
		}
		else
		{
			#check if genome studio file
			beta.p=champ.read(betaFile,sampleSheet)
            #add in filter for detP if it is a genome studio file
		}
        if(filterXY)
        {
            autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
            beta.p=beta.p[row.names(beta.p) %in% row.names(autosomes), ]
        }
	}else
	{
        if(norm == "SWAN")
		{
			mset <-preprocessSWAN(rgSet, mset)
			if(methValue=="B")
			{	
				beta.p = getBeta(mset, "Illumina")

			}else{beta.p = getM(mset)}
		}else{
        beta.p=beta
        }
        if(filterXY)
		{
			autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
            beta.p=beta.p[row.names(beta.p) %in% row.names(autosomes), ]	
		}
	}	
	if(norm=="BMIQ" | norm == "PBC")
	{
		### create design.v file
		data(probeInfoALL.lv)
		match(rownames(beta.p),probeInfoALL.lv[[5]]) -> map.idx
		mapto <- function(tmp.v){ return(tmp.v[map.idx])}
 		probeInfo.lv  <- lapply(probeInfoALL.lv,mapto)
		design.v <- probeInfo.lv[[2]]
		
		if(min(beta.p, na.rm=TRUE)==0)
		{
			message("Zeros in your dataset have been replaced with 0.000001")
			beta.p[beta.p==0]<-0.000001
		}

		if(norm == "BMIQ")
 		{
			newDir=paste(resultsDir,"Normalization",sep="/")
 			if(plotBMIQ){if(!file.exists(resultsDir)){dir.create(resultsDir)}
 			if(!file.exists(newDir))
			{
				dir.create(newDir)
			}}
 			design.v<-as.numeric(design.v)
 			bmiq=beta.p
			hf.v <- vector();

			if(plotBMIQ){setwd(newDir)}
 			for(s in 1:ncol(beta.p))
 			{
				sID=colnames(beta.p)[s]
				beta.v <- beta.p[,s];


				bmiq.o <- BMIQ(beta.v,design.v,doH=TRUE,nL=3,nfit=5000,niter=10,plots=plotBMIQ,sampleID=sID);
				bmiq[,s] <- bmiq.o$nbeta;
				hf.v[s] <- bmiq.o$hf;

			}
			setwd(cwd)
			beta.p=bmiq
		}else if(norm == "PBC")
		{
			beta.p=DoPBC(beta.p,design.v)
		}
	}else if(norm == "NONE")
	{
		
	}else{
		if(norm!="SWAN"){
		message("You have not selected a valid normalization method. No normalization will be preformed.")}
		}
	#else if(norm == "SubQ")
	#{
		#beta.norm=DoSubQ(beta.raw)
		#add normalization to colname _subQ

	#}
    
	totalProbes=dim(beta.p)[1]
	
	if(QCimages)
	{
		imageName1=paste(resultsDir,"norm_mdsPlot.pdf",sep="/")
        imageName2=paste(resultsDir,"norm_densityPlot.pdf",sep="/")
        imageName3=paste(resultsDir,"norm_SampleCluster.jpg",sep="/")
        main1=paste("Density plot of normalized data (",totalProbes," probes)",sep="")
        main2=paste("All samples after normalization (",totalProbes, " probes)",sep="")
		
		#save images...
		pdf(imageName1,width=6,height=4)
		mdsPlot(beta.p,numPositions=1000,sampGroups=pd$Sample_Group,sampNames=pd$Sample_Name)
		dev.off()
		
		pdf(imageName2,width=6,height=4)
		densityPlot(beta.p,sampGroups=pd$Sample_Group,main=main1,xlab=methValue)
		dev.off()
		
        if(ncol(beta.p) < 60)
        {
            #cluster
            jpeg(imageName3)
            betar_d<-dist(t(beta.p))
            plot(hclust(betar_d),main=main2,cex=0.8)
            dev.off()
    }

	}
	return(list(beta=beta.p))
}
