if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad"))

champ.QC <- function(beta = myLoad$beta,
                     pheno=myLoad$pd$Sample_Group,
                     mdsPlot=TRUE,
                     densityPlot=TRUE,
                     dendrogram=TRUE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     Feature.sel="None",
                     resultsDir="./CHAMP_QCimages/")
{
    message("[===========================]")
    message("[<<<<< ChAMP.QC START >>>>>>]")
    message("-----------------------------")
    ### Prepare Checking ###
    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.QC Results will be saved in ",resultsDir)
    message("[QC plots will be proceed with ",dim(beta)[1], " probes and ",dim(beta)[2], " samples.]\n")
    if(min(beta,na.rm=TRUE)==0)
    {
        beta[beta==0] <- 0.000001
        message("[",length(which(beta==0))," Zeros dectect in your dataset, will be replaced with 0.000001]\n")
    }
    if(ncol(beta)!=length(pheno)) stop("Dimension of DataSet Samples, pheno and name must be the same. Please check your input.")
    message("<< Prepare Data Over. >>")

    if(mdsPlot)
    {
        if(Rplot) mdsPlot(beta,numPositions=1000,sampGroups=pheno,colnames(beta))
        if(PDFplot){
            pdf(paste(resultsDir,"raw_mdsPlot.pdf",sep="/"),width=6,height=4)
            mdsPlot(beta,numPositions=1000,sampGroups=pheno,colnames(beta))
            dev.off()
        }
        message("<< plot mdsPlot Done. >>\n")
    }

    if(densityPlot)
    {
        if(Rplot) densityPlot(beta,sampGroups=pheno,main=paste("Density plot of raw data (",nrow(beta)," probes)",sep=""),xlab="Beta")
        if(PDFplot){
            pdf(paste(resultsDir,"raw_densityPlot.pdf",sep="/"),width=6,height=4)
            densityPlot(beta,sampGroups=pheno,main=paste("Density plot of raw data (",nrow(beta)," probes)",sep=""),xlab="Beta")
            dev.off()
        }
        message("<< Plot densityPlot Done. >>\n")
    }

    if(dendrogram)
    {
        if(Feature.sel=="None")
        {
            message("< Dendrogram Plot Feature Selection Method >: No Selection, directly use all CpGs to calculate distance matrix.")
            hc <- hclust(dist(t(beta)))
        }else if(Feature.sel=="SVD")
        {
            message("< Dendrogram Feature Selection Method >: Use top SVD CpGs to calculate distance matrix.")
            SVD <- svd(beta)
            rmt.o <- EstDimRMT(beta - rowMeans(beta))
            M <- SVD$v[,1:rmt.o$dim]
            rownames(M) <- colnames(beta)
            colnames(M) <- paste("Component",c(1:rmt.o$dim))
            hc <- hclust(dist(M))
        }
        dend <- as.dendrogram(hc)
        MyColor <- rainbow(length(table(pheno)))
        names(MyColor) <- names(table(pheno))
        labels_colors(dend) <- MyColor[pheno[order.dendrogram(dend)]]
        dend <- dend %>% set("labels_cex",0.8)
        dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex",0.6) %>% set("leaves_col",MyColor[pheno[order.dendrogram(dend)]])

        if(Rplot)
        {
            plot(dend,center=TRUE,main=paste("All samples before normalization (",nrow(beta), " probes)",sep=""))
            legend("topright",fill=MyColor,legend=names(MyColor))
        }
        if(PDFplot)
        {
            pdf(paste(resultsDir,"raw_SampleCluster.pdf",sep="/"),width=floor(log(ncol(beta))*3),height=floor(log(ncol(beta))*2))
            plot(dend,center=TRUE,main=paste("All samples before normalization (",nrow(beta), " probes)",sep=""))
            legend("topright",fill=MyColor,legend=names(MyColor))
            dev.off()
        }
        message("<< Plot dendrogram Done. >>\n")
    }

    message("[<<<<<< ChAMP.QC END >>>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.norm() next.]\n")
}
