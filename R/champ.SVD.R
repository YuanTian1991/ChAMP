champ.SVD <-
function(beta=myNorm$beta, rgSet=myLoad$rgSet,detP=myLoad$detP,pd=myLoad$pd, loadFile=FALSE, betaFile="beta.txt", sampleSheet="sampleSheet.txt", methProfile=FALSE, methFile="MethylationProbeProfile.txt", controlProfile=FALSE, controlFile="ControlProbeProfile.txt", studyInfo=FALSE,studyInfoFile="studyInfo.txt", infoFactor=c(),resultsDir=paste(getwd(),"resultsChamp",sep="/"))
{
    impute.knn<-NA
	rm(impute.knn)
    
	message("Performing SVD")
	if(!is.null(rgSet))
	{
		if(is.null(detP)){detP <- detectionP(rgSet)}
		if(is.null(beta))
		{
			mset <- preprocessRaw(rgSet)
			beta=getBeta(mset,"Illumina")
			
		}
	}else{
		if(is.null(detP))
		{
			message("You need a matrix of detection p-values to run SVD. Please rerun with parameter detP")
			return()
		}
	}
	
	if(loadFile)
	{
		beta.p=champ.read(betaFile,sampleSheet)
		beta.m=beta.p$beta
		pd=beta.p$pd
	}else{
		beta.m<-beta		
	}
	if(methProfile)
	{
		if(!file.exists(methFile)){message("You don't have a MethylationProbeProfile.txt file in this directory"); return()}
		tmp.df<- read.table(methFile,header=TRUE,sep="\t",fill=T);
		tmpTarget <- as.vector(tmp.df$TargetID);
		tmp.df<- tmp.df[,c(grep("X(.*)\\.",names(tmp.df)))]

		tmpB.idx <- grep("AVG_Beta",colnames(tmp.df));
		tmpI.idx <- grep("Intensity",colnames(tmp.df));
		tmpU.idx <- grep("Signal_A",colnames(tmp.df));
		tmpM.idx <- grep("Signal_B",colnames(tmp.df));
		tmpPV.idx <- grep("Detection",colnames(tmp.df));

		data.l <- list(B=as.matrix(tmp.df[,tmpB.idx]),I=as.matrix(tmp.df[,tmpI.idx]),U=as.matrix(tmp.df[,tmpU.idx]),M=as.matrix(tmp.df[,tmpM.idx]),pv=as.matrix(tmp.df[,tmpPV.idx]));
		for(i in 1:length(data.l))
		{
  			rownames(data.l[[i]]) <- tmpTarget;
		}
		rm(tmp.df,tmpB.idx,tmpI.idx,tmpU.idx,tmpM.idx,tmpPV.idx,tmpTarget);

		for(i in 1:length(data.l))
		{
  			tmp.v <- gsub(".AVG_Beta","",colnames(data.l[[i]]));
  			tmp.v <- gsub(".Detection.Pval","",tmp.v);
 			tmp.v <- gsub(".Intensity","",tmp.v);
  			tmp.v <- gsub(".Signal_A","",tmp.v);
  			tmp.v <- gsub(".Signal_B","",tmp.v);
  			colnames(data.l[[i]]) <- tmp.v;
		}
		beta.m <- data.l$B
		rm(tmp.v)
	}

	if(!is.null(rgSet))
	{
		bc1=getControlAddress(rgSet,controlType=c("BISULFITE CONVERSION I"))
		bc2=getControlAddress(rgSet,controlType=c("BISULFITE CONVERSION II"))
		ext=getControlAddress(rgSet,controlType=c("EXTENSION"))
		tr=getControlAddress(rgSet,controlType=c("TARGET REMOVAL"))
		hyb=getControlAddress(rgSet,controlType=c("HYBRIDIZATION"))
			
		#columns are samples
		CPP=rbind(getGreen(rgSet)[bc1[1:3],],getRed(rgSet)[bc1[7:9],],getRed(rgSet)[bc2[1:4],],getGreen(rgSet)[tr[1:2],],getGreen(rgSet)[hyb[1:3],],getRed(rgSet)[ext[1:2],],getGreen(rgSet)[ext[3:4],])
		controlNames <- c("BSC-I C1 Grn","BSC-I C2 Grn","BSC-I C3 Grn","BSC-I C4 Red","BSC-I C5 Red","BSC-I C6 Red","BSC-II C1 Red","BSC-II C2 Red","BSC-II C3 Red","BSC-II C4 Red","Target Removal 1 Grn","Target Removal 2 Grn","Hyb (Low) Grn","Hyb (Medium) Grn","Hyb (High) Grn","Extension (A) Red","Extension (T) Red","Extension (C) Grn","Extension (G) Grn");
		rownames(CPP)=controlNames
			
		#log2
		dataC2.m <- matrix(nrow=length(controlNames),ncol=ncol(beta.m));
		colnames(dataC2.m) <- colnames(beta.m);
		rownames(dataC2.m) <- controlNames;
		for(r in 1:nrow(dataC2.m))
		{
  			dataC2.m[r,] <- log2(as.numeric(CPP[r,]));
		}

	}else{	
		if(controlProfile)
		{	
			
		ControlProbeProfiles<-read.delim(controlFile)
		ControlProbeProfiles_s <- ControlProbeProfiles[,grep("Signal",colnames(ControlProbeProfiles))];
		ControlProbeProfiles_t <- ControlProbeProfiles[,grep("ID",colnames(ControlProbeProfiles))];
		ControlProbeProfiles <- cbind(ControlProbeProfiles_s,ControlProbeProfiles_t);

		dataC.m <- ControlProbeProfiles;

		dataC.m<- dataC.m[,c(grep("X(.*)\\.",names(dataC.m)))]
		dataC.l <- list();
		dataC.l[[1]] <- dataC.m[,grep("Signal_Grn",colnames(dataC.m))];
		dataC.l[[2]] <- dataC.m[,grep("Signal_Red",colnames(dataC.m))];
		names(dataC.l) <- c("Grn","Red");
		colnames(dataC.l$Grn) <- gsub("X","",gsub("_","",gsub(".Signal_Grn","",colnames(dataC.l$Grn))));
		colnames(dataC.l$Red) <- gsub("X","",gsub("_","",gsub(".Signal_Red","",colnames(dataC.l$Red))));

		rownames(dataC.l$Grn) <- ControlProbeProfiles$ProbeID
		rownames(dataC.l$Red) <- ControlProbeProfiles$ProbeID

		### build probe control data matrix
		seltypeC.v <- c("BSC-I C1 Grn", "BSC-I C2 Grn", "BSC-I C3 Grn", "BSC-I C4 Red", "BSC-I C5 Red", "BSC-I C6 Red", "BSC-II C1 Red", "BSC-II C2 Red", "BSC-II C3 Red", "BSC-II C4 Red", "Target Removal 1 Grn", "Target Removal 2 Grn", "Hyb (Low) Grn", "Hyb (Medium) Grn", "Hyb (High) Grn", "Extension (A) Red", "Extension (T) Red", "Extension (C) Grn", "Extension (G) Grn");
		seltypeC.idx <- c(1:3,7:9,13:16,833:834,21:23,17:20); #this chooses the specific control columns...
		grn.idx <- grep("Grn",seltypeC.v);
		red.idx <- grep("Red",seltypeC.v);

		## select specific controls
		dataC.m <- matrix(nrow=length(seltypeC.v),ncol=ncol(beta.m));
		colnames(dataC.m) <- colnames(beta.m);
		rownames(dataC.m) <- seltypeC.v;

		for(i in 1:ncol(dataC.m))
		{
			dataC.m[grn.idx,][,i] <- dataC.l$Grn[seltypeC.idx[grn.idx],][,i];
			dataC.m[red.idx,][,i] <- dataC.l$Red[seltypeC.idx[red.idx],][,i];
		}

		#log2
		dataC2.m <- matrix(nrow=length(seltypeC.v),ncol=ncol(beta.m));
		colnames(dataC2.m) <- colnames(beta.m);
		rownames(dataC2.m) <- seltypeC.v;
		for(r in 1:nrow(dataC2.m))
		{
  			dataC2.m[r,] <- log2(as.numeric(dataC.m[r,]));
		}
		dataC2.m=t(dataC2.m)
		}
	}

	PhenoTypes.lv <- list();

	################Customise Phenotype Data########################
	pd=pd[order(pd$Sample_Name),]
	well=unique(pd$Sample_Well)
	plates=unique(pd$Sample_Plate)
	groups=unique(pd$Sample_Group)
	pool=unique(pd$Pool_ID)
	chips=unique(pd$Slide)
	arrays=unique(pd$Array)


	allP=list(c(well,plates,groups,pool,chips,array))

	p=1
	a=1
	###Well
	if(length(well)>1)
	{
		tmp.v <- as.vector(pd$Sample_Well);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)); 
		for(i in 2:length(well))
		{
			PhenoTypes.lv[[p]][grep(well[i],tmp.v)] <- i; 
		}
		names(PhenoTypes.lv)[p]<-"Sample_Well"
		p=p+1
		a=a+1
	}
	
	##Plate
	if(length(plates)>1)
	{
		tmp.v <- as.vector(pd$Sample_Plate);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)) #1
		for(i in 2:length(plates))
		{
			PhenoTypes.lv[[p]][grep(plates[i],tmp.v)] <- 2
		}
		names(PhenoTypes.lv)[p]<-"Sample_Plate"
		p=p+1
		a=a+1
	}

	##Group
	if(length(groups)>1)
	{
		tmp.v <- as.vector(pd$Sample_Group);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		for(i in 2:length(groups))
		{
			PhenoTypes.lv[[p]][grep(groups[i],tmp.v)] <- i; #R
		}
		names(PhenoTypes.lv)[p]<-"Sample_Group"
		p=p+1
		a=a+1
	}
	
	##Pool_ID
	if(length(pool)>1)
	{
		tmp.v <- as.vector(pd$Pool_ID);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		for(i in 2:length(pool))
		{
			PhenoTypes.lv[[p]][grep(pool[i],tmp.v)] <- i;
		}
		names(PhenoTypes.lv)[p]<-"Pool_ID"
		p=p+1
		a=a+1
	}

	### Slide
	if(length(chips)>1)
	{
		tmp.v <- as.vector(pd$Slide);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		#check for just 1
		for(i in 2:length(chips))
		{
			PhenoTypes.lv[[p]][grep(as.numeric(chips[i]),tmp.v)] <- i; 
		}
		names(PhenoTypes.lv)[p]<-"Slide"
		p=p+1
		a=a+1
	}

	### Sentrix_Position
	if(length(arrays)>1)
	{
		tmp.v <- as.vector(pd$Array);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)); ## R01C01
		for(i in 2:length(arrays))
		{
			PhenoTypes.lv[[p]][grep(arrays[i],tmp.v)] <- i; ## R01C02
		}
		names(PhenoTypes.lv)[p]<-"Array"
		p=p+1
		a=a+1
	}
	########add in studyInfo
	
	b=1
	if(studyInfo)
	{
		tmp.info<- read.table(studyInfoFile,header=TRUE,sep="\t",fill=T,stringsAsFactors=FALSE,as.is=T)
		tmp.info=tmp.info[which(tmp.info$Sample_Name %in% pd$Sample_Name),]
		if(!is.null(tmp.info))
		{
			tmp.info=tmp.info[order(tmp.info$Sample_Name),]
			tmp.info=as.matrix(tmp.info)
			#match tmp.info[1] with pd=$Sample_Names
			for(n in 2:ncol(tmp.info))
			{	
				tmp.v <- as.vector(tmp.info[,n]);
				if(length(unique(tmp.v))>1)
				{
					PhenoTypes.lv[[p]] <- rep(1,length(tmp.v))
					phenos=unique(tmp.info[,n])
					for(i in 2:length(phenos))
					{
						PhenoTypes.lv[[p]][grep(phenos[i],tmp.v)]<-i
					}
					names(PhenoTypes.lv)[p]<-colnames(tmp.info)[n]
				}else{
					p=p-1
					b=b-1}
				p=p+1
				b=b+1
			}
		}else{ message("Your sample info file doesn't have the same samples as the sample sheet")}
	}
	p=p-1
	a=a-1
	b=b-1
	print(a)
	print(b)
	print(p)
	
	if(is.null(infoFactor))
	{
		factor.log <- c(rep(TRUE,p))
	}else{
		factor.log <- c(rep(TRUE,a))
		factor.log<- c(factor.log,infoFactor)
		}
	if(!is.null(rgSet) | controlProfile==T)
	{
		factor.log <- c(factor.log,rep(FALSE,19))
	}

	##DoQC.R
	common.idx <- 1:nrow(beta.m);
	selCL.idx <- 1:ncol(beta.m);

	### find coverage per sample
    ndet.v <- vector();
	det.li <- list();

	#detection p-value
    for(s in 1:ncol(beta.m))
	{
		det.li[[s]] <- which(detP[,s]<0.05);
	  	ndet.v[s] <- length(det.li[[s]]);
	}
	covS.v <- ndet.v/nrow(beta.m);
	
    ### global coverage
	for(s in selCL.idx)
	{
		common.idx <- intersect(common.idx,det.li[[s]]);
	}

	### impute remaining missing values

    capture.output(irawdataCL.m <- impute.knn(beta.m[common.idx,selCL.idx],k=5)$data)
	
    ### DoSVD.R
    tmp.m<- irawdataCL.m-rowMeans(irawdataCL.m)
	rmt.o <- EstDimRMTv2(tmp.m);
	svd.o <- svd(tmp.m);
    if(rmt.o$dim >6)
    {
        topPCA <- 6;
    }else{topPCA <- rmt.o$dim}
        
	#Check coloumn idx & change to appropriate coloumns from dataC2.m(below includes control profile, sentrix and clin.info)
    selcat.idx <- c(1:length(PhenoTypes.lv));
	if(!is.null(rgSet) | controlProfile)
	{
		svdPV.m <- matrix(nrow=topPCA,ncol=length(selcat.idx)+nrow(dataC2.m));
		colnames(svdPV.m) <- c(names(PhenoTypes.lv)[selcat.idx],rownames(dataC2.m));
		print(cbind(factor.log[1:length(selcat.idx)],names(PhenoTypes.lv)[selcat.idx]));
	}else
	{
		svdPV.m <- matrix(nrow=topPCA,ncol=length(selcat.idx));
		colnames(svdPV.m) <- c(names(PhenoTypes.lv)[selcat.idx]);
		print(cbind(factor.log[1:length(selcat.idx)],names(PhenoTypes.lv)[selcat.idx]));
	}
	for(c in 1:topPCA)
	{
		tmp.v <- svd.o$v[,c];
	  	for(f in 1:length(selcat.idx))
	  	{
    		if(factor.log[f])
    		{
    	  		svdPV.m[c,f] <- kruskal.test(tmp.v ~ as.factor(PhenoTypes.lv[[selcat.idx[f]]]))$p.value;
    		}else 
    		{
      			svdPV.m[c,f] <- summary(lm(tmp.v ~ PhenoTypes.lv[[f]]))$coeff[2,4];
    		}
  		}
  		if(!is.null(rgSet) | controlProfile)
  		{
  			for(f in (length(selcat.idx)+1):ncol(svdPV.m))
  			{
    			svdPV.m[c,f] <- summary(lm(tmp.v ~ dataC2.m[f-length(selcat.idx),selCL.idx]))$coeff[2,4];
  			}
  		}
	}
	### image heatmap
	myPalette <- c("darkred","red","orange","pink","white");
	breaks.v <- c(-200,-10,-5,-2,log10(0.05),0);
	imageName=paste(resultsDir,"SVDsummary.pdf",sep="/")
	pdf(imageName,width=8,height=8);
	par(mar=c(5,15,2,1),xpd=TRUE);
	image(x=1:nrow(svdPV.m), y=1:ncol(svdPV.m), z=log10(svdPV.m), col=myPalette, breaks=breaks.v, xlab="", ylab="", axes=FALSE, main= "Singular Value Decomposition Analysis (SVD)");
	axis(1,at=1:nrow(svdPV.m),labels=paste("PC-",1:nrow(svdPV.m),sep=""),las=2);
	suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
    legend(x=-3, y=5, legend=c(expression("p < 1x"~10^{-10}), expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), fill=c("darkred","red","orange","pink","white"));    
	dev.off();

	#imageName=paste(resultsDir,"ScatterPC12.pdf",sep="/")
	#pdf(imageName,width=10,height=10);
	#plot(svd.o$v[,1],svd.o$v[,2],col=PhenoTypes.lv$Sample_Group,pch=c(15,15,15)[PhenoTypes.lv$Group],xlab="PC1",ylab="PC2",main="Plot of Priciple Components 1 & 2 showing groups");
	#legend('topright',c(groups[1],groups[2]),fill=c("black","red"),bty="n");
	#dev.off();
}
