if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","mylimma","detectCores","probe.features","illumina450Gr","seqlevels<-","seqlevels"))

champ.DMR <- function(beta=myNorm,
                      pheno=myLoad$pd$Sample_Group,
                      arraytype="450K",
                      method = "Bumphunter",
                      minProbes=7,
                      adjPvalDmr=0.05,
                      cores=3,
                      ## following parameters are specifically for Bumphunter method.
                      maxGap=300,
                      cutoff=NULL,
                      pickCutoff=TRUE,
                      smooth=TRUE,
                      smoothFunction=loessByCluster,
                      useWeights=FALSE,
                      permutations=NULL,
                      B=250,
                      nullMethod="bootstrap",
                      ## following parameters are specifically for probe ProbeLasso method.
                      DMP=myDMP,
                      meanLassoRadius=375,
                      minDmrSep=1000,
                      minDmrSize=50,
                      adjPvalProbe=0.05,
                      Rplot=T,
                      PDFplot=T,
                      resultsDir="./CHAMP_ProbeLasso/",
                      ## following parameters are specifically for DMRcate method.
                      rmSNPCH=T,
                      dist=2,
                      mafcut=0.05,
                      lambda=1000,
                      C=2)
{
    message("[===========================]")
    message("[<<<<< ChAMP.DMR START >>>>>]")
    message("-----------------------------")


     if(length(which(is.na(beta)))>0) message(length(which(is.na(beta)))," NA are detected in your beta Data Set, which may cause fail or uncorrect of SVD analysis. You may want to impute NA with champ.impute() function first.")


    #if(arraytype=="EPIC"){
    #    data(probe.features.epic)
    #}else{
    #    data(probe.features)
    #}

    if(arraytype=="EPIC"){
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
    }else{
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
    }
    probe.features <- getAnnotation(RSobject)

    if(!class(pheno) %in% c("character","factor","numeric")) stop("pheno parameter must be a category vector, could be character, factor or numeric.")
    if(cores > detectCores()) cores <- detectCores()

    if(method=="Bumphunter")
    {
        message("<< Find DMR with Bumphunter Method >>")

        message(cores," cores will be used to do parallel BMIQ computing.")
        registerDoParallel(cores = cores)
        

        cpg.idx <- intersect(rownames(beta),rownames(probe.features))
        Anno <- probe.features[cpg.idx,]
        Anno <- Anno[order(Anno$chr,Anno$pos),]
        cpg.idx <- rownames(Anno)

        cl <- clusterMaker(Anno$chr,Anno$pos,maxGap=maxGap)
        names(cl) <- cpg.idx
        bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl)>minProbes)))]

        message("According to your data set, champ.DMR() detected ",sum(table(cl)>minProbes)," clusters contains MORE THAN ",minProbes," probes within",maxGap," maxGap. These clusters will be used to find DMR.\n")

        X <- cbind(1,(as.numeric(as.factor(pheno))-1))
        Beta <- beta[bumphunter.idx,]
        Beta <- replace(Beta,which(Beta <= 0.001),0.001)
        Beta <- replace(Beta,which(Beta >= 0.999),0.999)
        Y <- log((Beta/(1-Beta)),2)
        
        Bumps <- bumphunter(Y,
                            design=X,
                            chr=Anno[bumphunter.idx,]$chr,
                            pos=Anno[bumphunter.idx,]$pos,
                            cluster=cl[bumphunter.idx],
                            cutoff=cutoff,
                            pickCutoff=pickCutoff,
                            smooth=smooth,
                            smoothFunction=smoothFunction,
                            useWeights=useWeights,
                            permutations=permutations,
                            verbose=TRUE,
                            B=B,
                            nullMethod=nullMethod)

        message("<< Calculate DMR success. >>")
        DMR <- Bumps$table[which(Bumps$table$p.valueArea <= adjPvalDmr),]
        message("Bumphunter detected ",nrow(DMR)," DMRs with P value <= ",adjPvalDmr,".")

        if(nrow(DMR) == 0) stop("No DMR detected.")

        rownames(DMR) <- paste("DMR",1:nrow(DMR),sep="_")
        #DMRProbes <- apply(DMR,1,function(x) Anno[which(Anno$chr==x[1] & Anno$pos>= as.numeric(x[2]) & Anno$pos<= as.numeric(x[3])),])
        DMR <- data.frame(DMR[,1:3],width=DMR[,3]-DMR[,2],strand="*",DMR[,4:14])
        colnames(DMR)[1:3] <- c("seqnames","start","end") 

        #OutputDMR <- list(BumphunterDMR=DMR,BumphunterDMRProbes=DMRProbes)
        OutputDMR <- list(BumphunterDMR=DMR)

    } else if(method == "ProbeLasso")
    {
        if (!file.exists(resultsDir)) dir.create(resultsDir)
        message("champ.DMR Results will be saved in ",resultsDir)

        message("<< Find DMR with ProbeLasso Method >>")
        gc()
        if(arraytype=="EPIC") data(illuminaEPICGr) else data(illumina450Gr)
        if(length(which(DMP$adj.P.Val < adjPvalProbe))==0) stop("There is no probe show significant difference from champ.DMP() function.")


        myResultsGr <- illumina450Gr[match(rownames(DMP), names(illumina450Gr))]
        myResultsGr$P.Value <- DMP$P.Value[match(names(myResultsGr), rownames(DMP))];
        myResultsGr$adj.P.Val <- DMP$adj.P.Val[match(names(myResultsGr), rownames(DMP))]
        seqlevels(myResultsGr) <- sort(seqlevels(myResultsGr));
        myResultsGr <- sort(myResultsGr,ignore.strand=T) # sort for later
        ### readjust pValues after masking
        myResultsGr$adj.P.Val <- p.adjust(mcols(myResultsGr)$P.Value, method = "BH")
        ### Probe spacing and quantile derivation

        message("<< Get closestProbe for each Probe >>")
        closestProbe <- as.data.frame(distanceToNearest(myResultsGr,ignore.strand=T))$distance
        closestProbeSp <- split(closestProbe, mcols(myResultsGr)$featureCgi); rm(closestProbe)

        message("<< Get lassoQuantileThreshold for each featureCgi >>")
        lassoQuantileDeviation <- abs(meanLassoRadius - rowMeans(as.data.frame(lapply(closestProbeSp,function(x) quantile(x,(1:1000)/1000)))))
        lassoQuantileThreshold <- which.min(lassoQuantileDeviation) / 1000;
        lassoSizes <- lapply(closestProbeSp, function(x) quantile(x, lassoQuantileThreshold, na.rm = T))

        message("<< Get expend ranges for each probe >>")
        myResultsGrSp <- split(myResultsGr, myResultsGr$featureCgi) # splits myResultsGr by 'featureCgi'; length = 28
        lassoGr <- mapply(function(x, y) promoters(x, upstream = y, downstream = y, ignore.strand = TRUE), x = myResultsGrSp, y = lassoSizes)
        lassoGr <- unlist(GRangesList(lassoGr)); rm(myResultsGrSp)
        myResultsSigGr <- myResultsGr[which(mcols(myResultsGr)$adj.P.Val < adjPvalProbe)]
        lassoProbeCountOverlap <- countOverlaps(lassoGr, myResultsSigGr, ignore.strand = T);rm(myResultsSigGr)

        message("<< Get DMR from overlapped probes >>")
        dmrGr <- reduce(lassoGr[which(lassoProbeCountOverlap >= minProbes)], min.gapwidth = minDmrSep, ignore.strand=TRUE);
        rm(lassoProbeCountOverlap, lassoGr) # lassos capturing 'minSigProbesLasso', merged 
        strand(dmrGr) <- '*'
        dmrGr <- dmrGr[which(width(dmrGr) > minDmrSize)] # remove DMRs < minDmrSize
        probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
        pValuesGr <- myResultsGr[probeIndex$subjectHits, "P.Value"]
        myBetas <- beta[match(names(pValuesGr), rownames(beta)), ]
        myBetas <- split(as.data.frame(myBetas), probeIndex$queryHits)
            
        message("<< Get adjusted P value for DMR >>")
        correl <- lapply(myBetas, function(x) cor(t(x)))
        weights <- lapply(correl, function(x) 1/apply(x^2,1,sum)); rm(correl)
        dmrQP <- qnorm(mcols(pValuesGr)$P.Value); dmrQP <- split(dmrQP, probeIndex$queryHits)
        dmrQPW <- mapply("*", dmrQP, weights); rm(dmrQP)
        if(class(dmrQPW) == "matrix") dmrStat <- sum(dmrQPW) else dmrStat <- lapply(dmrQPW, sum)
        rm(dmrQPW)
        dmrSd <- lapply(weights, function(x) sqrt(sum(x^2))); rm(weights)
        dmrP <- mapply(function(x,y) pnorm(x,0, sd=y), dmrStat, dmrSd); rm(dmrStat, dmrSd)
        dmrP <- p.adjust(dmrP, method = "BH")
        goodDmr <- which(dmrP < adjPvalDmr)
        dmrGr <- dmrGr[goodDmr] 
        dmrP <- dmrP[goodDmr]
        dmrpRank <- rank(dmrP, ties.method="min"); rm(goodDmr)



        ### get pvalues and betas for GOOD DMRs
        message("<< Get Start-End Ranges for each DMR >>")
        probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
        dmrProbesGr <- myResultsGr[probeIndex$subjectHits]
        myBetas <- beta[match(names(dmrProbesGr), rownames(beta)), ]; myBetas <- as.data.frame(myBetas)
        dmrCoreStart <- start(dmrProbesGr)
        dmrCoreEnd <- end(dmrProbesGr)
        myBetas <- split(myBetas, probeIndex$queryHits)
        dmrCoreStart <- split(dmrCoreStart, probeIndex$queryHits); dmrCoreStart <- sapply(dmrCoreStart, min)
        dmrCoreEnd <- split(dmrCoreEnd, probeIndex$queryHits); dmrCoreEnd <- sapply(dmrCoreEnd, max)
        ### calculate methylation scores for each DMR in each Sample_Group

        message("<< Calculate Methylation Scores for each DMR >>")
        groupIndex <- pheno
        dmrGroupMeans <- do.call(rbind, lapply(myBetas, function(x) sapply(split(t(x), groupIndex), mean)))
        colnames(dmrGroupMeans) <- paste("betaAv", colnames(dmrGroupMeans), sep = "_")
        probeGroupMeans <- lapply(myBetas, function(x) split(as.data.frame(t(x)), groupIndex)); rm(groupIndex, myBetas)
        probeGroupMeans <- lapply(probeGroupMeans, function(x) lapply(x, colMeans))
        probeGroupMeans <- do.call(rbind, lapply(probeGroupMeans, function(x) t(do.call(rbind, x))))
        colnames(probeGroupMeans) <- paste("betaAv", colnames(probeGroupMeans), sep = "_")

        ### probe-level data and DMR metadata

        message("<< Generate Probe-level Data >>")
        myDmrProbesGr <- myResultsGr[probeIndex$subjectHits]
        myDmrProbesGr <- as(cbind(as.data.frame(myDmrProbesGr), probeGroupMeans), "GRanges");
        rm(probeGroupMeans)
        myDmrProbesGr$dmrNo <- probeIndex$queryHits
        myDmrProbesGr$dmrP <- dmrP[probeIndex$queryHits]
        myDmrProbesGr$dmrpRank <- dmrpRank[probeIndex$queryHits]
        myDmrProbesGr$dmrChrom <- seqnames(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrStart <- start(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrEnd <- end(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrSize <- width(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrCoreStart <- dmrCoreStart[probeIndex$queryHits]
        myDmrProbesGr$dmrCoreEnd <- dmrCoreEnd[probeIndex$queryHits]    
        myDmrProbesGr$dmrCoreSize <- myDmrProbesGr$dmrCoreEnd - myDmrProbesGr$dmrCoreStart + 1

        ### DMR metadata
        message("<< Generate DMR metadata >>")
        myDmrGr <- dmrGr
        myDmrGr$dmrNo <- unique(probeIndex$queryHits)
        myDmrGr$dmrP <- dmrP; rm(dmrP)
        myDmrGr$dmrpRank <- dmrpRank; rm(dmrpRank)
        myDmrGr$dmrChrom <- seqnames(dmrGr) 
        myDmrGr$dmrStart <- start(dmrGr)
        myDmrGr$dmrEnd <- end(dmrGr)
        myDmrGr$dmrSize <- width(dmrGr); rm(dmrGr)
        myDmrGr$dmrCoreStart <- dmrCoreStart
        myDmrGr$dmrCoreEnd <- dmrCoreEnd
        myDmrGr$dmrCoreSize <- myDmrGr$dmrCoreEnd - myDmrGr$dmrCoreStart + 1
        genes <- split(as.data.frame(myResultsGr)[probeIndex$subjectHits, c("ensemblID", "geneSymbol")], probeIndex$queryHits); rm(probeIndex)
        myDmrGr$ensemblID <- sapply(genes, function(x) paste(unique(unlist(strsplit(x$ensemblID, ";"))), collapse = ";"))
        myDmrGr$geneSymbol <- sapply(genes, function(x) paste(unique(unlist(strsplit(x$geneSymbol, ";"))), collapse = ";")); rm(genes)
        myDmrGr <- as(cbind(as.data.frame(myDmrGr), dmrGroupMeans), "GRanges"); rm(dmrGroupMeans) 

	    # plot of lassos
        interplot <- function(illumina450Gr,lassoSizes)
        {
            cgiNames <- levels(illumina450Gr$cgi)
            featureNames <- levels(illumina450Gr$feature)
            sfo <- c(6, 7, 3, 1, 4, 2, 5)
            scgio <- c(1, 4, 3, 2)
            sfcgio <- rep((sfo - 1) *4, each = 4) + rep(scgio, 7)
            lassoSizes <- round(unlist(lassoSizes))
            #imageName <- paste(getwd(),"myLassos.pdf",sep="/")
            #pdf(imageName,width=9,height=9)

            par(mar = c(7, 4, 4, 3)+0.5)
            plot(c(1,28), y=c(range(0.3*sqrt(lassoSizes))[1]*0.8, range(0.3*sqrt(lassoSizes))[2]*1.2), type="n", xaxt="n", xlab="", yaxt="n", ylab="lasso radius [bp]", main=paste("lasso quantile = ", round(lassoQuantileThreshold,2), "\nmean lasso radius = ", meanLassoRadius, "bp", sep = ""), bty="n")
            segments(1:28, rep(0,28), 1:28, 0.3*sqrt(lassoSizes[sfcgio]), lty=3, col="grey")
            points(1:28, 0.3*sqrt(lassoSizes[sfcgio]), pch=16, cex=0.3*sqrt(lassoSizes[sfcgio]), col=rep(rainbow(7,alpha=0.5)[sfo], each=4))
            text(1:28, 0.3*sqrt(lassoSizes[sfcgio]), lassoSizes[sfcgio], pos=3, cex=0.8)
            axis(1, at = 1:28, labels = rep(cgiNames[scgio], 7), las = 2)
            par(xpd = T)
            segments(seq(1, 28, 4), rep(-2.5, 7), seq(4, 28, 4), rep(-2.5, 7))
            mtext(text = featureNames[sfo], side = 1, at = seq(2.5, 28, 4), line = 5.5, las = 1, cex.axis = 1)
            axis(2, at=c(0,max(0.3*sqrt(lassoSizes))), labels=F)
        }
        if(Rplot) interplot(illumina450Gr,lassoSizes)
        if(PDFplot)
        {
            pdf(paste(resultsDir,"myLassos.pdf",sep="/"),width=9,height=9)
            interplot(illumina450Gr,lassoSizes)
            dev.off()
        }
        DMRProbes <- as.data.frame(myDmrProbesGr)
        DMRProbes <- data.frame(probe.features[rownames(DMRProbes),],DMRProbes[,which(colnames(DMRProbes)=="P.Value"):which(colnames(DMRProbes)=="dmrNo")])
        DMRProbes <- split(DMRProbes,DMRProbes$dmrNo)
        DMR <- as.data.frame(myDmrGr)

        message("ProbeLasso detected ",nrow(DMR)," DMRs with P value <= ",adjPvalDmr,".")
        if(nrow(DMR) == 0) stop("No DMR detected.")

        rownames(DMR) <- paste("DMR",DMR$dmrNo,sep="_")
        names(DMRProbes) <- rownames(DMR)

        #OutputDMR <- list(ProbeLassoDMR=DMR,ProbeLassoDMRProbes=DMRProbes)
        OutputDMR <- list(ProbeLassoDMR=DMR)

    }else if(method=="DMRcate")
    {

        message(cores," cores will be used to do parallel DMRcate computing.")

        message("<< Find DMR with DMRcate Method >>")
        myMs <- logit2(beta)
        if(rmSNPCH) myMs <- rmSNPandCH(myMs, dist=dist, mafcut=mafcut)
        design <- model.matrix(~ pheno)

        if(arraytype=="450K")
        {
            myannotation <- cpg.annotate(datatype="array", myMs,design=design,coef=ncol(design), analysis.type="differential",annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),what="M")
        }else
        {
            myannotation <- cpg.annotate(datatype="array", myMs,design=design,coef=ncol(design), analysis.type="differential",annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19"),what="M")
        }
        M <- do.call("cbind", lapply(myannotation, as.data.frame))
        colnames(M) <- names(myannotation)

        dmrcoutput <- dmrcate(myannotation, min.cpgs = minProbes, lambda=lambda, C=C,mc.cores = cores)
        data(dmrcatedata)
        DMR <- as.data.frame(extractRanges(dmrcoutput, genome = "hg19"))

        message("Bumphunter detected ",nrow(DMR)," DMRs with mafcut as= ",adjPvalDmr,".")
        if(nrow(DMR) == 0) stop("No DMR detected.")

        if(nrow(DMR)!=0)
        {
            DMRProbes <- apply(DMR,1,function(x) M[which(M[,3]==x[1] & M[,4]>= as.numeric(x[2]) & M[,4]<= as.numeric(x[3])),])
            rownames(DMR) <- paste("DMR",1:nrow(DMR),sep="_")
            names(DMRProbes) <- rownames(DMR)
            for(i in names(DMRProbes)) rownames(DMRProbes[[i]]) <- DMRProbes[[i]]$ID
            X <- lapply(DMRProbes,as.data.frame)
            DMRProbes <- lapply(X,function(x) cbind(ID=x[,1],probe.features[rownames(x),],x[,c(2,5,6,7)]))
            
            OutputDMR <- list(DMRcateDMR=DMR)
        }else
        {
            OutputDMR <- NULL
        }
        #OutputDMR <- list(DMRcateDMR=DMR,DMRcateDMRProbes=DMRProbes)

    } else {
        stop("Please assign correct DMR method: 'Bumphunter' or 'ProbeLasso'")
    }
    message("[<<<<<< ChAMP.DMR END >>>>>>]")
    message("[===========================]")
    message("[You may want to process DMR.GUI() or champ.GSEA() next.]\n")
    return(OutputDMR)
}
