if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","probe.features","probe.features.epic","bloodCtl"))

champ.CNA <- function(intensity=myLoad$intensity,
                      pheno=myLoad$pd$Sample_Group,
                      control=TRUE,
                      controlGroup="champCtls",
                      sampleCNA=TRUE,
                      groupFreqPlots=TRUE,
                      Rplot=FALSE,
                      PDFplot=TRUE,
                      freqThreshold=0.3,
                      resultsDir="./CHAMP_CNA",
                      arraytype="450K")
{
    message("[===========================]")
    message("[<<<<< ChAMP.CNA START >>>>>]")
    message("-----------------------------")

    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.CNA Results will be saved in ",resultsDir," .\n")

    if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)

    message("ChaMP.CNA does not provide batch Correct on intensity data now, but you can use champ.runCombat to correct slides batch yourself.")
	
	if(control)
	{
        message("<< Create Control Data >>")
        if(controlGroup != "champCtls" & !(controlGroup %in% pheno))
        {
        	message("You have chosen ", controlGroup, " as the reference and this does not exist in your sample sheet (column Sample_Group). The analysis will run with ChAMP blood controls.")
        	controlGroup="champCtls"
        }
        if(controlGroup == "champCtls")
        {
            message("<< Combining champ bloodCtl dataset into your intensity dataset as control >>")
            data(champBloodCtls)
            ctlIntensity=bloodCtl$intensity
            intensity <- cbind(intensity,ctlIntensity[rownames(intensity),])
            pheno <- c(pheno,bloodCtl$pd$Sample_Group)
            message("champ bloodCtl dataset contains only two samples, they will be used as control groups.")
        }
	}

	#Extracts names of samples 
	names <- colnames(intensity)
	#Quantile normalises intensities	
	intsqn <- normalize.quantiles(as.matrix(intensity))
	colnames(intsqn)<-names
	#Calculates Log2
	intsqnlog<-log2(intsqn)
    	
	if(control)
	{
        message("<< Calculate mean value difference between each sample to mean control samples >>")
        intsqnlogratio <- apply(intsqnlog[,which(!pheno %in% controlGroup)],2,function(x) x - rowMeans(as.data.frame(intsqnlog[,which(pheno %in% controlGroup)])))
	}else
	{
        message("<< Calculate mean value difference between each sample to mean all samples >>")
        intsqnlogratio <- apply(intsqnlog[,which(!pheno %in% controlGroup)],2,function(x) x - rowMeans(intsqnlog))
	}

    ints <- data.frame(intensity,probe.features[rownames(intensity),c("MAPINFO","CHR")])
	ints$MAPINFO <- as.numeric(ints$MAPINFO)

	#Replaces Chr X and Y with 23 and 24
	levels(ints$CHR)[levels(ints$CHR)=='X'] <- '23'
	levels(ints$CHR)[levels(ints$CHR)=='Y'] <- '24'

    message("<< Generate CHR and MAPINFO information >>")
	CHR <- ints$CHR
    MAPINFO <- ints$MAPINFO


	#Runs CNA and generates individual DNA Copy Number profiles
    sampleResult <- list()
	if(sampleCNA)
	{
	    message("<< Processing Samples >>")
		for(i in 1:ncol(intsqnlogratio))
		{
			CNA.object <- CNA(cbind(intsqnlogratio[,i]), CHR, MAPINFO ,data.type = "logratio", sampleid = paste(colnames(intsqnlogratio)[i],"qn"))
			smoothed.CNA.object <- smooth.CNA(CNA.object)
			segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)
			seg <- segment.smoothed.CNA.object$output
            sampleResult[[colnames(intsqnlogratio)[i]]] <- seg

            if(Rplot)
                plot(segment.smoothed.CNA.object, plot.type = "w", ylim=c(-6,6))
			if(PDFplot)
			{
				imageName <- paste(colnames(intsqnlogratio)[i],"qn.pdf",sep="")
				imageName=paste(resultsDir,imageName,sep="/")
				pdf(imageName)
				plot(segment.smoothed.CNA.object, plot.type = "w", ylim=c(-6,6))
				dev.off()
			}
		}
	}
	

    groupResult <- list()
	if(groupFreqPlots)
	{
        message("<< Processing Groups >>")
        if(control) groups <- setdiff(unique(pheno),controlGroup) else groups <- unique(pheno)	
        tmp_pheno <- pheno[pheno %in% groups]
		for(g in 1:length(groups))
		{
		
			data_group=intsqnlogratio[,which(tmp_pheno == groups[g])]
			ints_group=ints[,which(tmp_pheno == groups[g])]
			row.names(ints_group)=row.names(ints)
	
			group.CNA.object <- CNA(data_group, CHR, MAPINFO,data.type = "logratio", sampleid = paste(paste(colnames(data_group),tmp_pheno[which(tmp_pheno == groups[g])]),"qn"))
			group.smoothed.CNA.object <- smooth.CNA(group.CNA.object)
			group.segment.smoothed.CNA.object <- segment(group.smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)
            seg <- group.segment.smoothed.CNA.object$output
            groupResult[[groups[g]]] <- seg
		
			group.freq = glFrequency(group.segment.smoothed.CNA.object,freqThreshold)		
	
            innerplot <- function(ints,group.freq,g)
            {
                ints = ints[order(ints$CHR,ints$MAPINFO),]
                labels_chr <- data.matrix(summary(as.factor(ints$CHR)))
                test1<- data.frame(labels_chr,row.names(labels_chr) )
                test <- data.frame(unique(CHR))
                colnames(test) = c("chr")
                colnames(test1) = c("count","chr")
                F1 <- merge(test,test1, by="chr", sort=T)
                for(i in 2:length(row.names(F1))){F1[i,2] = F1[i-1,2] + F1[i,2] ; }
                F1$label <- NULL ; F1[1,3] <- F1[1,2] / 2 ;	
                for (i in 2:length(row.names(F1))){ F1[i,3] <- (F1[i,2]+F1[i-1,2])/2; }
	
                y1=group.freq$gain
                y2=group.freq$loss

                graphTitle = paste("Frequency Plot of ",groups[g]," Samples",sep="")
                #plot gain
                plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = graphTitle , ylim=range(-1, 1), xlab='Chromosome Number',  ylab=paste('Fraction of Samples with Gain or Loss (n=',dim(data_group)[2],")",sep=""),xaxs = "i", yaxs = "i")

                #plot loss
                points(y2, type='h', col="red")

                #label for chromosomes
                x= c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
                y= c(1:length(F1[,2]))
                axis(1, at = c(F1[,2]), labels =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
                axis(1, at = c(F1[,3]), labels =F1$chr, tick = FALSE );
                axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
            }
			#begin plot
            if(Rplot) innerplot(ints,group.freq,g)
            if(PDFplot)
            {
                imageName1=paste(groups[g],"_","FreqPlot.pdf",sep="")
                imageName1=paste(resultsDir,imageName1,sep="/")
                pdf(imageName1, width = 10.0, height = 9.0)
                innerplot(ints,group.freq,g)
                dev.off()
            }
		}
	}

    message("[<<<<<< ChAMP.CNA END >>>>>>]")
    message("[===========================]")
    return(list(sampleResult=sampleResult,groupResult=groupResult))
}
