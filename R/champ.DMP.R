if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","probe.features.epic","probe.features"))

champ.DMP <- function(beta = myNorm,
                      pheno = myLoad$pd$Sample_Group,
                      compare.group = NULL,
                      adjPVal = 0.05,
                      adjust.method = "BH",
                      arraytype = "450K",
                      ###The following options are parsed down to qqman package for Manhattan and Q-Q plots
                      resultsDir = "./CHAMP_DMPimages/",
                      man.plot = T,
                      probes = F,
                      chr = F,
                      dotsize = 1,
                      sug.line = F,
                      gen.line = T,
                      gen.alpha = 0.05,
                      q.plot = T)
{
    message("[===========================]")
    message("[<<<<< ChAMP.DMP START >>>>>]")
    message("-----------------------------")

    CalculateDMP <- function(beta,pheno,tmp_compare,adjPVal=adjPVal,adjust.method=adjust.method)
    {
        message("  -----------------------------")
        message("  Start to Compare : ",tmp_compare[1],", ",tmp_compare[2])
        p <- pheno[which(pheno %in% tmp_compare)]
        tmpbeta <- beta[,which(pheno %in% tmp_compare)]
        design <- model.matrix( ~ 0 + p)
        contrast.matrix <- makeContrasts(contrasts=paste(colnames(design)[2:1],collapse="-"), levels=colnames(design))
        message("  Contrast Matrix")
        print(contrast.matrix)
        fit <- lmFit(tmpbeta, design)
        fit2 <- contrasts.fit(fit,contrast.matrix)
        tryCatch(fit3 <- eBayes(fit2),
                 warning=function(w) 
                {stop("limma failed, No sample variance.\n")})
        DMP <- topTable(fit3,coef=1,number=nrow(tmpbeta),adjust.method=adjust.method,p.value=adjPVal)
        message("  You have found ",sum(DMP$adj.P.Val <= adjPVal), " significant DMPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
        message("  Calculate DMP for ",tmp_compare[1]," and ",tmp_compare[2]," done.")
        return(DMP)
    }

    message("!!! Important !!! New Modification has been made on champ.DMP(): \n")
    message("    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.\n")
    message("    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted \"pheno\" parameter is \"numeric\" type.\n")

    Compare <- NULL

    message("--------------------------------")

    if(is.null(pheno) | length(unique(pheno))<=1)
    {
        stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
    }else
    {
        if(length(pheno)!=ncol(beta)) stop("Your Pheno's length is not in accord with your beta value's ncol.")
        message("\n[ Section 1:  Check Input Pheno Start ]\n")
        if(class(pheno)=="numeric")
        {
            message("  You pheno is numeric type.")
            message("    pheno parameter contains :",length(pheno)," values.")
            message("    pheno parameter ranges from ",min(pheno)," to ",max(pheno))
            
        } else {
            message("  You pheno is ",class(pheno)," type.")
            message("    Your pheno information contains following groups. >>")
            sapply(unique(pheno),function(x) message("    <",x,">:",sum(pheno==x)," samples."))
            message("    [The power of statistics analysis on groups contain very few samples may not strong.]")

            if(length(unique(pheno)) == 2) 
            {
                message("    pheno contains only 2 phenotypes")
                if(is.null(compare.group))
                {
                    message("    compare.group parameter is NULL, two pheno types will be added into Compare List.")
                    Compare <- list(x1=unique(pheno))
                } else if(sum(compare.group %in% unique(pheno))==2)
                {
                    message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
                    Compare <- list(x1=unique(pheno))
                } else {
                    stop(" You have specified compare.group, but it's not in accord with your pheno parameter. Please recheck your compare.group or pheno.")
                }
            } else if(length(unique(pheno)) > 2) {
                message("    pheno contains ",length(unique(pheno))," phenotypes")

                if(is.null(compare.group))
                {
                    message("    compare.group parameter is NULL, EACH PAIR of phenotypes will be added into Compare List.")
                    Compare <- as.list(data.frame(combn(unique(pheno),2)))
                } else if(sum(compare.group %in% unique(pheno))==2) {
                    message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
                    Compare <- list(x1=sort(compare.group))
                } else {
                    stop("    Your pheno parameter contains multiple phenotypes, but values in your compare.group parameter are not all find in them. Please recheck your compare.group or pheno.")
                }
            } else {
                stop("    !!! Something wrong with your pheno. Please check your input.")
            }

            tmpnamelist <- vector()
            for(i in 1:length(Compare))
            {
                tmpname <- paste(Compare[[i]][1],Compare[[i]][2],sep="_to_")
                message("    ",tmpname," compare group : ",Compare[[i]][1],", ",Compare[[i]][2])
                tmpnamelist <- c(tmpnamelist,tmpname)
            }
            names(Compare) <- tmpnamelist
        }
        message("\n[ Section 1:  Check Input Pheno Done ]\n")
    }


    DMPs <- list()
    if(is.null(Compare)){
        message("\n[ Section 2:  Find Numeric Covariates Linear Regression CpGs Start ]\n")
      if(man.plot | q.plot) message("  P-values for all probes will be calculated for plotting.")
        df <- data.frame(pheno=pheno)
        model.matrix <- model.matrix(~ pheno, data=df)
        fit1 <- lmFit(beta,model.matrix)
        fit2 <- eBayes(fit1)
        DMP <- topTable(fit2,coef=2,number=nrow(beta),adjust.method=adjust.method,p.value=adjPVal)
        message("  You have found ",sum(DMP$adj.P.Val <= adjPVal), " significant DMPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
      if(man.plot | q.plot) {  
      DMPs[["plotDMP"]] <- suppressMessages(topTable(fit2,coef=2,number=nrow(beta),adjust.method=adjust.method,p.value=1))
        }
        if(sum(DMP$adj.P.Val <= adjPVal)!=0)
            DMPs[["NumericVariable"]] <- DMP
    } else {
        message("\n[ Section 2:  Find Differential Methylated CpGs Start ]\n")
      if(man.plot | q.plot) message("  P-values for all probes will be calculated for plotting, but only for first comparison pair, ",Compare[[1]][1]," and ",Compare[[1]][2],".")
        for(i in names(Compare))
        {
            DMP <- CalculateDMP(beta,pheno,Compare[[i]],adjPVal,adjust.method)
            if(sum(DMP$adj.P.Val <= adjPVal)!=0)
                DMPs[[i]] <- DMP
        }
      if(man.plot | q.plot) {  
      DMPs[["plotDMP"]] <- suppressMessages(CalculateDMP(beta,pheno,Compare[[1]],1,adjust.method))
        }
    }
    message("\n[ Section 2:  Find Numeric Vector Related CpGs Done ]\n")

    if(length(DMPs)==0) stop("ChAMP.DMP Have not detected even one significant CpGs. You may try other threshold.")

    message("\n[ Section 3:  Match Annotation Start ]\n")
    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    if(man.plot | q.plot) Compare[["plotDMP"]] <- Compare[[1]]
    for(i in names(DMPs))
    {
        com.idx <- intersect(rownames(DMPs[[i]]),rownames(probe.features))
        if(!is.null(Compare))
        {
            avg <-  cbind(rowMeans(beta[com.idx,which(pheno==Compare[[i]][1])]),rowMeans(beta[com.idx,which(pheno==Compare[[i]][2])]))
            avg <- cbind(avg,avg[,2]-avg[,1])
            colnames(avg) <- c(paste(Compare[[i]],"AVG",sep="_"),"deltaBeta")
            DMPs[[i]] <- data.frame(DMPs[[i]][com.idx,],avg,probe.features[com.idx,])
        } else {
            DMPs[[i]] <- data.frame(DMPs[[i]][com.idx,],probe.features[com.idx,])
        }
    }
    message("\n[ Section 3:  Match Annotation Done ]\n")
    
    if(man.plot==T | q.plot==T) {
    message("\n[ Section 4:  Plotting Start ]\n")
      
    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message(" champ.DMP plots will be saved in ",resultsDir)
   
    if(probes!=F & !is.character(probes)){
      message(" 'probes' argument must be FALSE or a character string. Setting to FALSE.")
      probes = F
      }
    
      if(chr!=F & !is.numeric(chr)) {
      message(" 'chr' argument must be FALSE or a numeric string. Setting to FALSE.")
      chr = F
      }
    
      if(chr==F) chr = as.numeric(levels(DMPs[["plotDMP"]][,"CHR"])[levels(DMPs[["plotDMP"]][,"CHR"]) %in% unique(DMPs[["plotDMP"]][,"CHR"])])
    
    ##Transfer necessary data, set coloumn names
    DMP <- DMPs[["plotDMP"]]
    plot <- data.frame(rownames(DMP),(DMP["CHR"]),"",DMP["P.Value"],(DMP["MAPINFO"]))
    colnames(plot) <- c("SNP","CHR","BP","P","MAPINFO")

    if(man.plot) {

    ##Check if probes exist in data
    if(probes==F){
      message("No probes will be highlighted")
        } else if(probes!=F & isTRUE(probes %in% plot$SNP[plot$CHR %in% chr]==T)) {
      message("Probe list OK, highlighting ",length(probes)," probe(s)")
          } else {
      errP <- matrix(setdiff(probes,plot$SNP[plot$CHR %in% chr]))
      colnames(errP) <- c("")
      message("The following probe(s) for highlighting were not in the selected data:")
      print(errP[,1])
      stop("Please check for spelling errors, or if probes are on the selected chromosome(s)")
      }

    ##Prepare data
    plot$CHR <- as.integer(as.character(factor(plot$CHR)))
    plot$SNP <- as.character(plot$SNP)
    plot$MAPINFO <- as.numeric(plot$MAPINFO)
    plot <- plot[order(plot$CHR,plot$MAPINFO),]

    ##Count probes on CHR
    ch <- data.frame(table(plot$CHR))

    ##Number probes
    BP = c()
    for(i in unique(ch[,1])){BP <- c(BP,1:ch[i,2])}
    plot$BP <- BP

    ##Make plot
    message(" Drawing Manhattan plot")
    if(gen.line) {
    gen.line <- gen.alpha/dim(plot)[1]
    message(" Genome-wide line drawn at p=",formatC(gen.line,format="e", digits=2))
    }
      
    if(sug.line!=F)  message(" Suggestive line drawn at p=",formatC(sug.line,format="e", digits=2))
      
    if(is.numeric(gen.line) && is.numeric(sug.line)) {
    tiff(paste(resultsDir,"Manhattan.tiff",sep=""), width=1024, height=425)
    suppressWarnings(qqman::manhattan(subset(plot, CHR %in% chr), main="Manhattan plot", cex=dotsize, suggestiveline=-log10(sug.line), genomewideline=-log10(gen.line), highlight=probes))
    dev.off()
    } else if(is.numeric(gen.line) && !is.numeric(sug.line)) {
    tiff(paste(resultsDir,"Manhattan.tiff",sep=""), width=1024, height=425)
    suppressWarnings(qqman::manhattan(subset(plot, CHR %in% chr), main="Manhattan plot", cex=dotsize, suggestiveline=F, genomewideline=-log10(gen.line), highlight=probes))
    dev.off()
    } else if(!is.numeric(gen.line) && is.numeric(sug.line)) {
    tiff(paste(resultsDir,"Manhattan.tiff",sep=""), width=1024, height=425)
    suppressWarnings(qqman::manhattan(subset(plot, CHR %in% chr), main="Manhattan plot", cex=dotsize, suggestiveline=-log10(sug.line), genomewideline=F, highlight=probes))
    dev.off()
    } else {
    tiff(paste(resultsDir,"Manhattan.tiff",sep=""), width=1024, height=425)
    suppressWarnings(qqman::manhattan(subset(plot, CHR %in% chr), main="Manhattan plot", cex=dotsize, suggestiveline=F, genomewideline=F, highlight=probes))
    dev.off()
    }
    }

    if(q.plot) {
    #Q-Q plot
    message(" Drawing Q-Q plot")
    tiff(paste(resultsDir,"QQplot.tiff",sep=""), width=512, height=512)
    qqman::qq(plot$P, main="Q-Q Plot of EWAS p-values")
    dev.off()
    }
      
    message("\n[ Section 4:  Plotting Done ]\n")
    }
    message("[<<<<<< ChAMP.DMP END >>>>>>]")
    message("[===========================]")
    message("[You may want to process DMP.GUI() or champ.GSEA() next.]\n")
    return(DMPs)
}
