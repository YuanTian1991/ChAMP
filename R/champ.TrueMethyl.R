if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","probe.features"))

champ.TrueMethyl <- function(beta.norm = myNorm$beta,
                             pd=myLoad$pd,
                             adjPVal=0.05,
                             adjust.method="BH",
                             compare.group=c("oxBS","BS"),
                             resultsDir=paste(getwd(),"resultsChamp",sep="/"),
                             bedFile=TRUE,
                             arraytype="450K")
{	
	makeContrasts<-NA
	rm(makeContrasts)
	lmFit<-NA
	rm(lmFit)
	contrasts.fit<-NA
	rm(contrasts.fit)
	eBayes<-NA
	rm(eBayes)
	topTable<-NA
	rm(topTable)

    if(arraytype=="EPIC"){
        data(probe.features.epic)
    }else{
        data(probe.features)
    }
	
    groupLabel=compare.group
    
    checkLabel=unique(pd$Sample_Group)
    
    if(any(!(groupLabel %in% checkLabel))){
        message("The group labels that have been defined ",groupLabel[1]," and ", groupLabel[2]," do not exist in your sample sheet. Please edit the Sample_Group column or the compare.group parameter. ChAMP will use information in your Sample_Group columnn.")
        groupLabel = checkLabel
    }
    
	if(length(groupLabel)>2)
	{
        message("Your dataset has more than two groups. ChAMP will compare the first two groups.")
    }
	
        oxBS=pd[which(pd$Sample_Group==groupLabel[1]),]
        BS=pd[which(pd$Sample_Group==groupLabel[2]),]
        all=c(oxBS$Sample_Name,BS$Sample_Name)
        data=matrix(NA,length(row.names(beta.norm)),length(all))
        row.names(data)=row.names(beta.norm)
        colnames(data)=all
        message(groupLabel[1]," will be treated as your oxBS sample and ", groupLabel[2], " will be treated as your BS sample.")
	
	#quicker way to do this?
	if(length(beta.norm)==0)
	{
		message("Your dataset is empty and MVP list cannot be produced")
		return()
	}
    #change this loop...
	for(i in 1:length(all))
	{
		done=F
		j=1
		while(done==F)
		{
			if(colnames(data)[i]==colnames(beta.norm)[j])
			{
				data[,i]=beta.norm[,j] 
				done=T
                
			}else{
				if(j<ncol(beta.norm))
				{
					j=j+1
				}else{

					print("There is an error in this dataset")
				}
			}
		}
		
	}

	numsamples=length(oxBS$Sample_Name)
	numtest=length(BS$Sample_Name)
	design <- model.matrix(~0 + factor(c(rep("oxBS",numsamples), rep("BS",numtest))))
	print(paste("contrast ",groupLabel[1]," (oxBS sample) ", "against ", groupLabel[2]," (BS sample).",sep=""))
	colnames(design) <- c("oxBS", "BS")
	contrast.matrix <- makeContrasts(oxBS-BS, levels=design)

	fit <- lmFit(data, design)
	fit2 <- contrasts.fit(fit, contrast.matrix)

	ok=TRUE
	tryCatch(fit3 <- eBayes(fit2),
      warning=function(w) 
      {
      	cat("No sample variance.\n")
      	ok <- FALSE
      })
      if(ok)
      {
          results<- topTable(fit3, coef=1, number=dim(data)[1], adjust.method=adjust.method,p.value=adjPVal)

          message("You have found ", dim(results)[1], " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal)
          if(dim(results)[1]==0)
          {
              message("No bedfile will be generated for tophits but a full MVP list with all p-values is being saved")
          }else{
              
          if(colnames(results[1])=="ID")
          {
                  row.names(results)=results$ID
                  colnames(results)[1]<-"probeID"
                  
              }else{
                  results$probeID=row.names(results)
                  
              }
          
          resList=data[which(row.names(data) %in% results$probeID),]
	
          #join with annotation
          resList_anno<-data.frame(probe.features[match(row.names(resList),row.names(probe.features)),],resList)
          resList_anno$probeID=row.names(resList_anno)        
          
          if(bedFile)
          {
              
              if(!is.null(resList_anno))
              {
              	bedfile=champ.bedfile(resList_anno)
              	#add pvalue to file name
              	fileName1=paste(resultsDir,"/MVP_",adjPVal,"_",groupLabel[1],"vs",groupLabel[2],"_",adjust.method,"adjust_",dim(resList)[1],".bed",sep="")
              	write.table(bedfile,fileName1,row.names=F,col.names=F,quote=F, sep = "\t")
	         }
          }
       }
          
          resultsALL<- topTable(fit3, coef=1, number=dim(data)[1], adjust.method=adjust.method,p.value=1)
          if(colnames(resultsALL[1])=="ID")
          {
              row.names(resultsALL)=resultsALL$ID
              colnames(resultsALL)[1]<-"probeID"
          }else{
              resultsALL$probeID=row.names(resultsALL)
              
          }
          resultsALL_anno<-data.frame(resultsALL,probe.features[match(row.names(resultsALL),row.names(probe.features)),])

          oxBS.data=data[,which(colnames(data) %in% oxBS$Sample_Name),drop=FALSE]
          BS.data=data[,which(colnames(data) %in% BS$Sample_Name),drop=FALSE]
          data=data.frame(data)
          data$oxBS_AVG=rowMeans(oxBS.data)
          data$BS_AVG=rowMeans(BS.data)
          data$hmC<-data$BS_AVG-data$oxBS_AVG
          
          negative<- data[which(data$hmC < 0),]
          sig.neg=negative[row.names(negative) %in% row.names(results),]
          percent=(dim(sig.neg)[1]/(dim(results)[1]))*100
          percent=round(percent,digits=2)
          message("Your data has a total of ", dim(sig.neg)[1], " negative values. This represents ", percent, "% of your results. These will be filtered out and considered false positives but please consider whether this number is out of the acceptable range.")
          positive=data[which(data$hmC > 0),]
          positive=positive[which(rownames(positive) %in% rownames(results)),]
           message("After filtering for negative values you have found ", dim(positive)[1], " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal)
           
          colnames(data)[length(data)-2]=paste(groupLabel[1],"_oxBS_AVG",sep="")
          colnames(data)[length(data)-1]=paste(groupLabel[2],"_BS_AVG",sep="")
          data1=data[(length(data)-2):length(data)]
          resultsALL_anno<-data.frame(resultsALL_anno,data1[match(row.names(resultsALL_anno),row.names(data1)),])
          
          fileName2=paste(resultsDir,"/MVP_ALL_",groupLabel[1],"vs",groupLabel[2],"_",adjust.method,"adjust.txt",sep="")
          write.table(resultsALL_anno, fileName2 ,quote=F,sep="\t",row.names=F)
          return(resultsALL_anno)
          
      }else(print("There are no significant MVPs"))

}
