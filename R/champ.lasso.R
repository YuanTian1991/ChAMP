champ.lasso <-
function(fromFile=FALSE, uploadResults=FALSE, uploadFile="limma.txt", limma, beta.norm=myNorm$beta, pd=myLoad$pd, filterXY=TRUE, image=TRUE, mafPol.lower=0, mafPol.upper=0.05, popPol="eur",lassoStyle="max", lassoRadius=2000, minSigProbesLasso=3, minDmrSep=1000, minDmrSize=0, adjPVal=0.05, adjust.method="BH", resultsDir=paste(getwd(),"resultsChamp",sep="/"), bedFile=TRUE, DMRpval=0.05, batchDone=FALSE, normSave)
{
	data(probe.features)
	data(probe.450K.VCs.af)
    count<-NA
    rm(count)
    
	message("Run Probe Lasso DMR Hunter")
	
    if(uploadResults)
	{
		resultsFile <- read.table(uploadFile, header=T)
	}else if(fromFile){
		resultsFile<-limma
		
	}else{
		resultsFile=champ.MVP(beta.norm=beta.norm, pd=pd, adjPVal=adjPVal, bedFile=bedFile, adjust.method=adjust.method, resultsDir=resultsDir)
	}
    
    #myResults <- data.frame(resultsFile[,grep("adj.P.Val", colnames(resultsFile))])
    #rownames(myResults) <- rownames(resultsFile)
    #colnames(myResults) <- c("adj.P.Val")
    #myResults <- data.frame(myResults, probe.features[match(rownames(myResults), rownames(probe.features)),c(1:4,12:14)])
    
    myResults <- data.frame("P.Value" = resultsFile$P.Value, "adj.P.Val" = resultsFile$adj.P.Val, probe.features[match(rownames(resultsFile), rownames(probe.features)), ], row.names = rownames(resultsFile))
    
    if(filterXY) # converts $CHR into numeric for sorting purposes
	{
		myResults <- myResults[which(!(myResults$CHR=="X" | myResults$CHR=="Y")),]
		myResults$CHR <- as.numeric(as.character(myResults$CHR))
        
    }else{
		myResults$CHR <- as.character(myResults$CHR)
		myResults$CHR <- replace(myResults$CHR, which(myResults$CHR == "X"), "23")
		myResults$CHR <- replace(myResults$CHR, which(myResults$CHR == "Y"), "24")
		myResults$CHR <- as.numeric(myResults$CHR)
	}
    emptyDMRList=as.data.frame(setNames(replicate(21,numeric(0),simplify=F),c("probeID","adj.P.Val","CHR","MAPINFO","arm","gene.1","feature","cgi","feat.cgi","pol.af.f","pol.af.r","lasso.radius","dmr.no","dmr.start","dmr.end","dmr.size","dmr.core,start","dmr.core.end","dmr.core.size","dmr.p","deltaBeta")))
    
    if(dim(myResults)[1]==0)
    {
		message("Your dataset is empty and champ.lasso() cannot proceed")
		return(emptyDMRList)
        
	}else {
        if(min(myResults$adj.P.Val)>adjPVal)
        {
            if(batchDone)
            {
                message("The adusted p-values in your dataset exceed the cutoff you have chosen champ.lasso() cannot proceed. You might like to rerun champ.lasso with the normalised beta values without batchCorrect.")
                return(emptyDMRList)
                #champ.lasso(beta.norm=normSave,pd=pd,resultsDir=resultsDir,bedFile=bedFile,batchDone=F)
                
            }else{
                message("The adusted p-values in your dataset exceed the cutoff you have chosen and champ.lasso() cannot proceed.")
                return(emptyDMRList)
            }
            
        }else
        if(count(myResults$adj.P.Val<adjPVal)[which(count(myResults$adj.P.Val<adjPVal)$x==TRUE),][1,2] < 3)
        {
            message("There are not enough MVPs in your dataset for champ.lasso() to proceed.")
            return(emptyDMRList)
        }
        
        pop.ref <- data.frame("popPol" = c("asn", "amr", "afr", "eur"), "col.f"=c(9, 13, 17, 21), "col.r"=c(10,14,18,22))
        vc <- data.frame("ID" = rownames(probe.450K.VCs.af), "pol.af.f" = probe.450K.VCs.af[,pop.ref$col.f[match(popPol, pop.ref$popPol)]], "pol.af.r" = probe.450K.VCs.af[,pop.ref$col.r[match(popPol, pop.ref$popPol)]])
        myPol <- cbind(vc$pol.af.f[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.f[match(rownames(myResults), vc$ID)], vc$pol.af.r[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.r[match(rownames(myResults), vc$ID)])
        myPol <- cbind(apply(myPol[, 1:2], 1, min), apply(myPol[, 3:4], 1, min))
        myResults <- data.frame(myResults, "pol.af.f" = myPol[,1], "pol.af.r" = myPol[,2])
        rm(pop.ref, vc, myPol)
        
        if(mafPol.lower == 0 & mafPol.upper == 0 )
        {
            myResults <- myResults[myResults$pol.af.f == 0 & myResults$pol.af.r == 0,]
        }else{
            
            myResults <-myResults[myResults$pol.af.f >= mafPol.lower &
            myResults$pol.af.f <= mafPol.upper & #changed to and
            myResults$pol.af.r >= mafPol.lower &
            myResults$pol.af.r <= mafPol.upper,]
        }
       	myResults$adj.P.Val <- p.adjust(myResults$P.Value, method = "BH")
        myResults <- myResults[order(myResults$CHR, myResults$MAPINFO),]
        
        ###Probe spacing and quantile derivation
        
        myResults.split <- split(myResults, paste(myResults$CHR, myResults$arm), drop=T) #split by chromosome & arm
        
        
        probe.spacing <- lapply(myResults.split, function(x) apply(cbind(c(diff(x$MAPINFO), tail(diff(x$MAPINFO), n = 1)), c(head(diff(x$MAPINFO), n = 1), diff(x$MAPINFO))), 1, min))	# nearest probe calculations
        myResults <- data.frame(do.call(rbind, myResults.split), "nrst.probe" = unlist(probe.spacing))
        rm(myResults.split, probe.spacing)
        
        
        lasso.quantiles <- do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.cgi), function(x) ecdf(x)(lassoRadius)))
        if(lassoStyle == "max")
        {
            value.lasso.quantile <- min(lasso.quantiles)
        }else{
            value.lasso.quantile <- max(lasso.quantiles)
        }
        rm(lasso.quantiles)
        lasso.radius <- round(do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.cgi), function(x) quantile(x, value.lasso.quantile, na.rm=T))))
        myResults$lasso.radius <- lasso.radius[match(myResults$feat.cgi, rownames(lasso.radius))]
        
        myResults.sig <- myResults[myResults$adj.P.Val < adjPVal,]
        
        lasso.gr <- GRanges(seqnames=paste("chr", myResults.sig$CHR, sep = ""), ranges=IRanges(start=myResults.sig$MAPINFO - myResults.sig$lasso.radius, end=myResults.sig$MAPINFO + myResults.sig$lasso.radius))
        probe.gr <- GRanges(seqnames=paste("chr", myResults.sig$CHR, sep = ""), ranges=IRanges(start=myResults.sig$MAPINFO, end=myResults.sig$MAPINFO))
        lasso.probe.countOverlap <- countOverlaps(lasso.gr, probe.gr)
        
        myResults.sig <- data.frame(myResults.sig, lasso.probe.countOverlap)
        probeKeepers <- which(lasso.probe.countOverlap >= minSigProbesLasso)
        
        
        if(length(probeKeepers)==0)
        {
            return(emptyDMRList)
        }else{
            
        myResults.sig.cap <- myResults.sig[probeKeepers,] # retains rows (probes) that capture no less than "no.probes"
        myResults.sig.cap <- myResults.sig.cap[order(myResults.sig.cap$CHR, myResults.sig.cap$MAPINFO),] # orders object by chromosome then position
        rm(lasso.gr, probe.gr, lasso.probe.countOverlap, myResults.sig, probeKeepers)
        
        #big code replacement before this point...
        
        lasso.coord <- data.frame("CHR" = myResults.sig.cap$CHR, "arm" = myResults.sig.cap$arm, "lasso.start"=myResults.sig.cap$MAPINFO - myResults.sig.cap$lasso.radius, "lasso.end"=myResults.sig.cap$MAPINFO + myResults.sig.cap$lasso.radius)
        lasso.coord <- lasso.coord[order(lasso.coord$CHR, lasso.coord$lasso.start),]
        lasso.seq <- split(lasso.coord, paste(lasso.coord$CHR, lasso.coord$arm))
        rm(lasso.coord)
        lasso.bp <- vector("list", length(lasso.seq)) # genomic coordinates of every base within all lassos, by chromosome
        names(lasso.bp) <- names(lasso.seq)
        for (k in 1:length(lasso.bp)) # k = CHR
        {
            dd <- vector("list", nrow(lasso.seq[[k]]))
            for (i in 1:length(dd)) # i = lasso bounds
            {
                dd[[i]] <- seq(lasso.seq[[k]][i, 3], lasso.seq[[k]][i, 4])
            }
            lasso.bp[[k]] <- sort(do.call(c, dd))
        }
        
        rm(lasso.seq)
        lasso.bp <- lapply(lasso.bp, unique) # leaves unique genome coordinates. NB names are factors-based
        diffs <- unlist(lapply(lasso.bp, function(x) c(FALSE, diff(x) <= minDmrSep)))
        chr <- sapply(strsplit(names(lasso.bp), " "), "[[", 1)
        chr.un <- unlist(chr)
        chr.un.nu <- as.numeric(chr.un)
        rm(chr, chr.un)
        len1 <- lapply(lasso.bp, length)
        len.un <- unlist(len1)
        len.un.vec <- as.vector(len.un)
        rm(len1, len.un)
        chr <- rep(chr.un.nu, len.un.vec)
        rm(chr.un.nu)
        bp <- unlist(lasso.bp)
        
        dmr.start.index <- which(diffs == FALSE)
        dmr.end.index <- c(dmr.start.index[-1] - 1, length(bp))
        dmrs <- data.frame("CHR" = chr[dmr.start.index],
        "dmr.start" = bp[dmr.start.index],
        "dmr.end" = bp[dmr.end.index])
        dmrs <- dmrs[order(dmrs$CHR, dmrs$dmr.start), ]
        rm(chr, diffs, bp, dmr.start.index, dmr.end.index)
        dmrs$dmr.size <- (dmrs$dmr.end - dmrs$dmr.start) + 1
        dmrs <- dmrs[dmrs$dmr.size >= minDmrSize,]
        rownames(dmrs) <- 1:nrow(dmrs) # renames rows after removing small DMRs
        dmrs$dmr.no <- 1:nrow(dmrs) # renames DMRs after removing small DMRs
        
        dmrs.gr <- GRanges(seqnames=paste("chr", dmrs$CHR, sep = ""), ranges=IRanges(start=dmrs$dmr.start, end=dmrs$dmr.end))
        probe.gr <- GRanges(seqnames=paste("chr", myResults$CHR, sep = ""), ranges=IRanges(start=myResults$MAPINFO, end=myResults$MAPINFO))
        dmr.probes.gr <- findOverlaps(dmrs.gr, probe.gr)
        dmr.probes <- as.data.frame(dmr.probes.gr)
        rm(dmrs.gr, probe.gr, dmr.probes.gr)
        
        toMatch <- c("dmr.no", "dmr.start", "dmr.end", "dmr.size")
        dmr.col.keep <- match(toMatch, colnames(dmrs))
        myDf <- data.frame(myResults[dmr.probes$subjectHits, ], dmrs[match(dmr.probes$queryHits, rownames(dmrs)), dmr.col.keep ])
        
        core.start <- lapply(split(myDf$MAPINFO, myDf$dmr.no), min)
        core.end <- lapply(split(myDf$MAPINFO, myDf$dmr.no), max)
        core.start.un <- unlist(core.start)
        core.end.un <- unlist(core.end)
        len <- as.vector(table(as.factor(myDf$dmr.no)))
        myDf$dmr.core.start <- rep(core.start.un, len)
        myDf$dmr.core.end <- rep(core.end.un, len)
        myDf$dmr.core.size <- (myDf$dmr.core.end - myDf$dmr.core.start) + 1
        rm(core.start, core.start.un, core.end, core.end.un, len)


#don't think this works
        toMatch <- paste(c("p.", "q."), collapse = "|")
        rownames(myDf) <- sapply(strsplit(rownames(myDf), toMatch), "[[", 2)
        
        rm(dmrs, toMatch)
        
        if(filterXY) # retro-converts $CHR into a factor
        {
            myDf$CHR <- as.factor(as.character(myDf$CHR))
        }else
        {
            myDf$CHR <- as.character(myDf$CHR)
            myDf$CHR <- replace(myDf$CHR, which(myDf$CHR == 23), "X")
            myDf$CHR <- replace(myDf$CHR, which(myDf$CHR == 24), "Y")
            myDf$CHR <- as.factor(myDf$CHR)
        }
        
        dmr.beta.means=resultsFile[match(rownames(myDf),rownames(resultsFile)),]
        myDf <- data.frame(myDf, dmr.beta.means[,17:19,])
        rm(dmr.beta.means)
        
        dmr.beta <- split(as.data.frame(beta.norm[match(rownames(myDf), rownames(beta.norm)),]), myDf$dmr.no)
        corel <- lapply(dmr.beta, function(x) cor(t(x)))
        weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
        dmr.ind.p <- split(myDf$adj.P.Val, myDf$dmr.no)
        dmr.qp <- lapply(dmr.ind.p, qnorm)
        dmr.qp.w <- mapply("*", dmr.qp, weights)
        
        if(class(dmr.qp.w) == "matrix")
        {
            dmr.stat <- sum(dmr.qp.w)
        }else
        {
            dmr.stat <- lapply(dmr.qp.w, sum)
        }
        
        dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
        dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
        dmr.rep <- as.numeric(summary(dmr.ind.p)[,1])
        
        #changed this
        dmr.p <- p.adjust(dmr.p, method = "BH")
        dmr.p.rank <- rank(dmr.p, ties.method="min")
        
        myDf$dmr.p <- rep(dmr.p, dmr.rep)
        
        myDf$dmr.p.rank <- rep(dmr.p.rank, dmr.rep)
        colKeep1 <- grep("AVG", colnames(myDf))
        dmr.means <- do.call(rbind, lapply(split(myDf[, colKeep1], myDf$dmr.no), colMeans))
        deltaBeta=dmr.means[,2]-dmr.means[,1]
        toMatch<-c("CHR","arm")
        if("gene" %in% colnames(myDf))
        {
            toMatch[length(toMatch)+1]="gene"
            
        }else{
            toMatch[length(toMatch)+1]="gene.1"
        }
        toMatch[4:12]<-c("dmr.no", "dmr.start", "dmr.end", "dmr.size", "dmr.core.start", "dmr.core.end", "dmr.core.size","dmr.p", "dmr.p.rank")
        colKeep2 <- c(match(toMatch, colnames(myDf)))
        myDfTrunc <- do.call(rbind, lapply(split(myDf[, colKeep2], myDf$dmr.no), function(x) x[1, ]))
        index <- match(c("dmr.p", "dmr.p.rank"), colnames(myDfTrunc))
        myDfTrunc <- data.frame(myDfTrunc[, 1:(index[1]-1)], dmr.means, deltaBeta, myDfTrunc[, index])
        rm(dmr.beta, corel, weights, dmr.ind.p, dmr.qp, dmr.qp.w, dmr.stat, dmr.sd, dmr.p, dmr.rep, dmr.p.rank, colKeep1, dmr.means, toMatch, colKeep2, index)
        
        if(image)
        {
            plotThreshold(myResults,lasso.radius,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
        }
        
        message("You have found ",max(myDf$dmr.no)," DMRs.")
        
        dmrs <- myDfTrunc[order(myDfTrunc$dmr.no),]
        
        sig.dmrs=myDfTrunc[which(myDfTrunc$dmr.p < DMRpval),]
        message("You have found ",length(unique(myDfTrunc$dmr.no))," significant DMRs with a dmr.p < ",DMRpval,".")
        
        
        fileName=paste(resultsDir,"/DMRs_",DMRpval,"_",max(myDfTrunc$dmr.no),".txt",sep="")
        write.table(myDfTrunc, fileName ,quote=F,sep="\t",row.names=F)
        
        fileName2=paste(resultsDir,"/DMRprobes_","_",dim(myDf)[1],".txt",sep="")
        write.table(myDf, fileName2 ,quote=F,sep="\t",row.names=F)
        
        ########make bedfile
        if(bedFile)
        {
            bedfile=data.frame(matrix(NA,ncol=4,nrow=nrow(myDfTrunc)))
            colnames(bedfile)[1]="chr"
            colnames(bedfile)[2]="start"
            colnames(bedfile)[3]="end"
            colnames(bedfile)[4]="gene"
            bedfile$chr<-myDfTrunc$CHR
            bedfile$start<-myDfTrunc$dmr.start
            bedfile$end<-myDfTrunc$dmr.end
            bedfile$gene<-myDfTrunc$gene
            bedfile<-bedfile[order(bedfile$chr,bedfile$start),]

            #add pvalue to file name
            fileName1=paste(resultsDir,"/DMR_",DMRpval,"_",max(myDfTrunc$dmr.no),".bed",sep="")
            write.table(bedfile,fileName1,row.names=F,col.names=F,quote=F, sep = "\t")
        }
    }
    return(list("dmr.probes" = myDf, "dmrs" = myDfTrunc))
}}
