if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","mylimma","detectCores","probe.features","illumina450Gr","seqlevels<-","seqlevels"))

champ.DMR <- function(betaNorm=myNorm$beta,
                      design=myLoad$pd$Sample_Group,
                      maxGap=300,
                      cutoff=0.5,
                      minProbes=7,
                      smooth=TRUE,
                      smoothFunction=loessByCluster,
                      useWeights=FALSE,
                      permutations=NULL,
                      B=250,
                      pickCutoff=FALSE,
                      pickCutoffQ=0.99,
                      nullMethod="bootstrap",
                      verbose=TRUE,
                      cores=3,
                      arraytype="450K",
                      method = "Bumphunter",
                      ## following parameter are for probe lasso method.
                      resultsFile=mylimma,
                      meanLassoRadius=375,
                      minSigProbesLasso=3,
                      minDmrSep=1000,
                      minDmrSize=50,
                      adjPvalProbe=0.05,
                      adjPvalDmr=0.05,
                      pData=myLoad$pd)
{
    if(length(which(is.na(betaNorm)))!=0)
    {
        stop("NA value is not allowed in betaNorm matrix. Please impute betaNorm matrix first.")
    }

    if(method=="Bumphunter")
    {
        if(arraytype=="EPIC"){
            data(probe.features.epic)
        }else{
            data(probe.features)
        }
        if(!class(design) %in% c("character","factor","numeric"))
        {
            stop("design parameter must be a vector, could be character, factor or numeric.")
        }

        if(is.null(cores)){
            registerDoParallel(cores = detectCores())
        }else{
            registerDoParallel(cores = cores)
        }

        cpg.idx <- intersect(rownames(betaNorm),rownames(probe.features))
        Anno <- probe.features[cpg.idx,]
        Anno <- Anno[order(Anno$CHR,Anno$MAPINFO),]
        cpg.idx <- rownames(Anno)

        cl <- clusterMaker(Anno$CHR,Anno$MAPINFO,maxGap=maxGap)
        names(cl) <- cpg.idx
        bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl)>minProbes)))]

        X <- cbind(rep(1,length(design)),(as.numeric(as.factor(design))-1))
        Beta <- betaNorm[bumphunter.idx,]
        Beta <- replace(Beta,which(Beta <= 0.001),0.001)
        Beta <- replace(Beta,which(Beta >= 0.999),0.999)
        Y <- log((Beta/(1-Beta)),2)
        
        Bumps <- bumphunter(Y,
                            design=X,
                            chr=Anno[bumphunter.idx,]$CHR,
                            pos=Anno[bumphunter.idx,]$MAPINFO,
                            cluster=cl[bumphunter.idx],
                            cutoff=cutoff,
                            pickCutoff=pickCutoff,
                            pickCutoffQ=pickCutoffQ,
                            smooth=smooth,
                            smoothFunction=smoothFunction,
                            useWeights=useWeights,
                            permutations=permutations,
                            B=B,
                            verbose=verbose,
                            nullMethod=nullMethod)
        return(list(myDmrs=Bumps,myDmrProbes=Anno[bumphunter.idx,]))
    } else if(method == "ProbeLasso"){
        gc()
        if(arraytype=="EPIC"){
            data(illuminaEPICGr)
        }else{
            data(illumina450Gr)
        }
        #data(illumina450Gr)
        if(length(which(resultsFile$adj.P.Val < adjPvalProbe))==0){
        stop("There is no probe show significant difference from champ.MVP() function.")
    }
        ### create annotation file with study p.values ###
        myResultsGr <- illumina450Gr[match(rownames(resultsFile), names(illumina450Gr))]
        myResultsGr$P.Value <- resultsFile$P.Value[match(names(myResultsGr), rownames(resultsFile))];
        myResultsGr$adj.P.Val <- resultsFile$adj.P.Val[match(names(myResultsGr), rownames(resultsFile))]
        seqlevels(myResultsGr) <- sort(seqlevels(myResultsGr));

        myResultsGr <- sort(myResultsGr,ignore.strand=T) # sort for later

        ### readjust pValues after masking
        myResultsGr$adj.P.Val <- p.adjust(mcols(myResultsGr)$P.Value, method = "BH")
	
        ### Probe spacing and quantile derivation
        closestProbe <- as.data.frame(distanceToNearest(myResultsGr,ignore.strand=T))$distance
	
        ### calculate lasso sizes for each featureCgi
        closestProbeSp <- split(closestProbe, mcols(myResultsGr)$featureCgi); rm(closestProbe)

        lassoQuantileDeviation <- numeric(1000)
        for (i in 1:1000)
        {
            lassoQuantileDeviation[i] <- abs(meanLassoRadius - mean(sapply(closestProbeSp, function(x) quantile(x, i/1000))))
        }
        lassoQuantileThreshold <- which(lassoQuantileDeviation == min(lassoQuantileDeviation)) / 1000; rm(lassoQuantileDeviation) 
        lassoSizes <- lapply(closestProbeSp, function(x) quantile(x, lassoQuantileThreshold, na.rm = T)); rm(closestProbeSp)

        ### create GRanges object of lassos
        myResultsGrSp <- split(myResultsGr, myResultsGr$featureCgi) # splits myResultsGr by 'featureCgi'; length = 28
        lassoGr <- mapply(function(x, y) promoters(x, upstream = y, downstream = y, ignore.strand = TRUE), x = myResultsGrSp, y = lassoSizes)
        lassoGr <- unlist(GRangesList(lassoGr)); rm(myResultsGrSp)

        ### get probes capturing at least 'minSigProbesLasso'
        myResultsSigGr <- myResultsGr[which(mcols(myResultsGr)$adj.P.Val < adjPvalProbe)]
        lassoProbeCountOverlap <- countOverlaps(lassoGr, myResultsSigGr, ignore.strand = T);
        rm(myResultsSigGr)

        ### create DMR GRanges object
        dmrGr <- reduce(lassoGr[which(lassoProbeCountOverlap >=minSigProbesLasso)], min.gapwidth = minDmrSep, ignore.strand=TRUE);
        rm(lassoProbeCountOverlap, lassoGr) # lassos capturing 'minSigProbesLasso', merged 
        strand(dmrGr) <- '*'
        dmrGr <- dmrGr[which(width(dmrGr) > minDmrSize)] # remove DMRs < minDmrSize

        ### get pvalues and betas for probes in DMRs
        probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
        pValuesGr <- myResultsGr[probeIndex$subjectHits, "P.Value"]
        myBetas <- betaNorm[match(names(pValuesGr), rownames(betaNorm)), ]

        myBetas <- split(as.data.frame(myBetas), probeIndex$queryHits)
			
        ### calculate weights and adjust pValues for Dmr
        correl <- lapply(myBetas, function(x) cor(t(x)))
        weights <- lapply(correl, function(x) 1/apply(x^2,1,sum)); rm(correl)
        dmrQP <- qnorm(mcols(pValuesGr)$P.Value); dmrQP <- split(dmrQP, probeIndex$queryHits)
        dmrQPW <- mapply("*", dmrQP, weights); rm(dmrQP)

        if(class(dmrQPW) == "matrix") # in the case of a single DMR, a single matrix will be returned
        {
            dmrStat <- sum(dmrQPW)
        }else
        {
            dmrStat <- lapply(dmrQPW, sum)#
        }; 
        rm(dmrQPW)
        dmrSd <- lapply(weights, function(x) sqrt(sum(x^2))); rm(weights)
        dmrP <- mapply(function(x,y) pnorm(x,0, sd=y), dmrStat, dmrSd); rm(dmrStat, dmrSd)
        dmrP <- p.adjust(dmrP, method = "BH")
        goodDmr <- which(dmrP < adjPvalDmr)
        dmrGr <- dmrGr[goodDmr] 
        dmrP <- dmrP[goodDmr]
        dmrpRank <- rank(dmrP, ties.method="min"); rm(goodDmr)

        ### get pvalues and betas for GOOD DMRs
        probeIndex <- as.data.frame(findOverlaps(dmrGr, myResultsGr))
        dmrProbesGr <- myResultsGr[probeIndex$subjectHits]
        myBetas <- betaNorm[match(names(dmrProbesGr), rownames(betaNorm)), ]; myBetas <- as.data.frame(myBetas)
        dmrCoreStart <- start(dmrProbesGr)
        dmrCoreEnd <- end(dmrProbesGr)

        myBetas <- split(myBetas, probeIndex$queryHits)
        dmrCoreStart <- split(dmrCoreStart, probeIndex$queryHits); dmrCoreStart <- sapply(dmrCoreStart, min)
        dmrCoreEnd <- split(dmrCoreEnd, probeIndex$queryHits); dmrCoreEnd <- sapply(dmrCoreEnd, max)

        ### calculate methylation scores for each DMR in each Sample_Group
        groupIndex <- pData$Sample_Group[match(colnames(betaNorm), rownames(pData))]
        dmrGroupMeans <- do.call(rbind, lapply(myBetas, function(x) sapply(split(t(x), groupIndex), mean)))
        colnames(dmrGroupMeans) <- paste("betaAv", colnames(dmrGroupMeans), sep = "_")
        probeGroupMeans <- lapply(myBetas, function(x) split(as.data.frame(t(x)), groupIndex)); rm(groupIndex, myBetas)
        probeGroupMeans <- lapply(probeGroupMeans, function(x) lapply(x, colMeans))
        probeGroupMeans <- do.call(rbind, lapply(probeGroupMeans, function(x) t(do.call(rbind, x))))
        colnames(probeGroupMeans) <- paste("betaAv", colnames(probeGroupMeans), sep = "_")

        ### probe-level data and DMR metadata
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
	    cgiNames <- levels(mcols(illumina450Gr)$cgi)
	    featureNames <- levels(mcols(illumina450Gr)$feature)
	    sfo <- c(6, 7, 3, 1, 4, 2, 5)
	    scgio <- c(1, 4, 3, 2)
	    sfcgio <- rep((sfo - 1) *4, each = 4) + rep(scgio, 7)
	    lassoSizes <- round(unlist(lassoSizes))
	    imageName <- paste(getwd(),"myLassos.pdf",sep="/")
	    pdf(imageName,width=9,height=9)

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
	    dev.off()
	    rm(lassoSizes, lassoQuantileThreshold, cgiNames, featureNames, sfo, scgio, sfcgio)
        return(list(myDmrProbes= myDmrProbesGr,myDmrs= myDmrGr))
    }else{
        stop("Please assign correct DMR method: 'Bumphunter' or 'ProbeLasso'")
    }
}
