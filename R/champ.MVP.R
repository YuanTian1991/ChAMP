champ.MVP <-
function(beta.norm = myNorm$beta, pd=myLoad$pd, adjPVal=0.05, adjust.method="BH", compare.group=c("C","T"),resultsDir=paste(getwd(),"resultsChamp",sep="/"),bedFile=TRUE)
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
	data(probe.features)	
	
    groupLabel=compare.group
    
    checkLabel=unique(pd$Sample_Group)
    
    if(!(groupLabel %in% checkLabel)){
        message("The group labels that have been defined ",groupLabel," do not exist in your sample sheet. Please edit the Sample_Group column or the compare.group parameter. ChAMP will use information in your Sample_group columnn.")
        groupLabel = checkLabel
    }
    
	if(length(groupLabel)>2)
	{
        message("Your dataset has more than two groups. ChAMP will compare the first two groups.")
    }
	
        controls=pd[which(pd$Sample_Group==groupLabel[1]),]
        test=pd[which(pd$Sample_Group==groupLabel[2]),]
        all=c(controls$Sample_Name,test$Sample_Name)
        data=matrix(NA,length(row.names(beta.norm)),length(all))
        row.names(data)=row.names(beta.norm)
        colnames(data)=all
	
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

	numsamples=length(controls$Sample_Name)
	numtest=length(test$Sample_Name)
	design <- model.matrix(~0 + factor(c(rep("C",numsamples), rep("T",numtest))))
	print(paste("contrast",groupLabel[1],groupLabel[2],sep=" "))
	colnames(design) <- c("C", "T")
	contrast.matrix <- makeContrasts(T-C, levels=design)

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

          control.data=data[,which(colnames(data) %in% controls$Sample_Name)]
          test.data=data[,which(colnames(data) %in% test$Sample_Name)]
          data=data.frame(data)
          data$C_AVG=rowMeans(control.data)
          data$T_AVG=rowMeans(test.data)
          data$deltaBeta<-data$T_AVG-data$C_AVG

          colnames(data)[length(data)-2]=paste(groupLabel[1],"_AVG",sep="")
          colnames(data)[length(data)-1]=paste(groupLabel[2],"_AVG",sep="")
          data1=data[(length(data)-2):length(data)]
          resultsALL_anno<-data.frame(resultsALL_anno,data1[match(row.names(resultsALL_anno),row.names(data1)),])
          
          fileName2=paste(resultsDir,"/MVP_ALL_",groupLabel[1],"vs",groupLabel[2],"_",adjust.method,"adjust.txt",sep="")
          write.table(resultsALL_anno, fileName2 ,quote=F,sep="\t",row.names=F)
          return(results.file=resultsALL_anno)
          
      }else(print("There are no significant MVPs"))

}
