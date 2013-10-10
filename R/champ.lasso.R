champ.lasso <-
function(fromFile=FALSE, uploadResults=FALSE, uploadFile="limma.txt", limma, beta.norm=myNorm$beta, pd=myLoad$pd, filterXY=TRUE, image=TRUE, mafPol.lower=0, mafPol.upper=0.05, popPol="eur",lassoStyle="max", lassoRadius=2000, minSigProbesLasso=3, minDmrSep=1000, minDmrSize=0, adjPVal=0.05, adjust.method="BH", resultsDir=paste(getwd(),"resultsChamp",sep="/"), bedFile=TRUE, DMRpval=0.05, batchDone=FALSE, normSave)
{
	data(probe.features)
	data(probe.450K.VCs.af)

	message("Run Probe Lasso DMR Hunter")
    
    go = T
	
    if(uploadResults)
	{
		resultsFile <- read.table(uploadFile, header=T)
	}else if(fromFile){
		resultsFile<-limma
		
	}else{
		resultsFile=champ.MVP(beta.norm=beta.norm, pd=pd, adjPVal=adjPVal, bedFile=bedFile, adjust.method=adjust.method, resultsDir=resultsDir)
	}	

    myResults <- data.frame(resultsFile[,grep("adj.P.Val", colnames(resultsFile))])
    rownames(myResults) <- resultsFile$probeID
    colnames(myResults) <- c("adj.P.Val")
    myResults <- data.frame(myResults, probe.features[match(rownames(myResults), rownames(probe.features)),c(1:4,20)])

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
    
    if(dim(myResults)[1]==0)
    {
		message("Your dataset is empty and champ.lasso() cannot proceed")
		return(beta.norm)
        
	}else 
        if(min(myResults$adj.P.Val)>adjPVal)
        {
            if(batchDone)
            {
                message("The adusted p-values in your dataset exceed the cutoff you have chosen champ.lasso() cannot proceed. You might like to rerun champ.lasso with the normalised beta values without batchCorrect.")
                return(beta.norm)
                #champ.lasso(beta.norm=normSave,pd=pd,resultsDir=resultsDir,bedFile=bedFile,batchDone=F)
            
            }else{
                message("The adusted p-values in your dataset exceed the cutoff you have chosen and champ.lasso() cannot proceed.")
                return(beta.norm)
                }
                
        }else 
            if(count(myResults$adj.P.Val<adjPVal)[which(count(myResults$adj.P.Val<adjPVal)$x==TRUE),][1,2] < 3)
            {
                message("There are not enough MVPs in your dataset for champ.lasso() to proceed.")
                return(beta.norm)
            }
    #else{}
    
        
        pop.ref <- data.frame("popPol" = c("asn", "amr", "afr", "eur"), "col.f"=c(9, 13, 17, 21), "col.r"=c(10,14,18,22))
        vc <- data.frame("ID" = rownames(probe.450K.VCs.af), "pol.af.f" = probe.450K.VCs.af[,pop.ref$col.f[match(popPol, pop.ref$popPol)]], "pol.af.r" = probe.450K.VCs.af[,pop.ref$col.r[match(popPol, pop.ref$popPol)]])
        myPol <- cbind(vc$pol.af.f[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.f[match(rownames(myResults), vc$ID)], vc$pol.af.r[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.r[match(rownames(myResults), vc$ID)])
        myPol <- cbind(apply(myPol[, 1:2], 1, min), apply(myPol[, 3:4], 1, min))
        myResults <- data.frame(myResults, "pol.af.f" = myPol[,1], "pol.af.r" = myPol[,2])

        if(mafPol.lower == 0 & mafPol.upper == 0 )
        {
            myResults <- myResults[myResults$pol.af.f == 0 & myResults$pol.af.r == 0,]
        }else{

            myResults <-myResults[myResults$pol.af.f >= mafPol.lower &
                        myResults$pol.af.f <= mafPol.upper |
                        myResults$pol.af.r >= mafPol.lower &
                        myResults$pol.af.r <= mafPol.upper,]
            }
        myResults <- myResults[order(myResults$CHR, myResults$MAPINFO),]
    
        ###Probe spacing and quantile derivation
    
        myResults.split <- split(myResults$MAPINFO, paste(myResults$CHR, myResults$arm), drop=T) #split by chromosome & arm
    
        myResults.split.rownames <- split(rownames(myResults), paste(myResults$CHR, myResults$arm), drop=T) # rownames split by chromosome
        probe.spacing <- lapply(myResults.split, diff)	# nearest probe calculations
        up.down <- data.frame("upstream" = unlist(lapply(probe.spacing, function(x) c(0, x))),
        "downstream" = unlist(lapply(probe.spacing, function(x) c(x,0))), row.names = unlist(myResults.split.rownames))
        up.down$nrst.probe <- apply(up.down, 1, min)
        rpl <- apply(up.down[which(up.down$nrst.probe == 0),],1,max)
        up.down$nrst.probe <- replace(up.down$nrst.probe, which(up.down$nrst.probe == 0), rpl)

        myResults <- data.frame(myResults, "nrst.probe" = up.down$nrst.probe[match(rownames(myResults), rownames(up.down))])

    
        lasso.quantiles <- do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.rel), function(x) ecdf(x)(lassoRadius)))
        if(lassoStyle == "max")
        {
            value.lasso.quantile <- min(lasso.quantiles)
        }else{
            value.lasso.quantile <- max(lasso.quantiles)
        }
    
        lasso.radii <- round(do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.rel), function(x) quantile(x, value.lasso.quantile, na.rm=T))))
        myResults$lasso.radii <- lasso.radii[match(myResults$feat.rel, rownames(lasso.radii))]
	
        myResults.sig <- myResults[myResults$adj.P.Val < adjPVal,]
        
    
        myResults.sig.splitarm <- split(myResults.sig[,c(3,ncol(myResults.sig))], paste(myResults.sig$CHR, myResults.sig$arm), drop=T) #split MAPINFO & lasso.radii by chromosome & arm
        
    
        myResults.sig.rownames <- split(rownames(myResults.sig), paste(myResults.sig$CHR, myResults.sig$arm), drop=T) # rownames split by chromosome

        no.pr <- vector("list", length(myResults.sig.splitarm)) # counts no. of probes captured by lasso
        for (k in 1:length(myResults.sig.splitarm))
        {
            for (i in 1:nrow(myResults.sig.splitarm[[k]]))
            {
                no.pr[[k]][i] <- 	length(which(myResults.sig.splitarm[[k]][,1] > myResults.sig.splitarm[[k]][i,1] - myResults.sig.splitarm[[k]][i,2] &
                myResults.sig.splitarm[[k]][,1] < myResults.sig.splitarm[[k]][i,1] + myResults.sig.splitarm[[k]][i,2]))
            }
        }
    
        
        tmp1 <- data.frame("lassoRadius" = do.call("rbind", myResults.sig.splitarm)[,2], "lasso.captured.probes"=unlist(no.pr), row.names = unlist(myResults.sig.rownames))
    
        if(minSigProbesLasso > max(tmp1$lasso.captured.probes))
        {
            message("Your lassos have not managed to capture at least ", minSigProbesLasso, " probes. No DMRs can be found")
            return(beta.norm)
        }else{
    
        tmp2 <- data.frame(myResults.sig[match(rownames(tmp1), rownames(myResults.sig)), ], "lasso.captured.probes"=tmp1$lasso.captured.probes)
        tmp3 <- tmp2[tmp2$lasso.captured.probes >= minSigProbesLasso,] # retains rows (probes) that capture no less than "no.probes"
        tmp3 <- tmp3[order(tmp3$CHR, tmp3$MAPINFO),] # orders file by chromosome then position
	
        lasso.coord <- data.frame("CHR" = tmp3$CHR, "lasso.start"=tmp3$MAPINFO - tmp3$lasso.radii, "lasso.end"=tmp3$MAPINFO + tmp3$lasso.radii)
        lasso.coord <- lasso.coord[order(lasso.coord$CHR, lasso.coord$lasso.start),]
	
        lasso.diffs <- numeric(nrow(lasso.coord)) # looking for overlapping lassos
        for (i in 2:length(lasso.diffs))
        {
            lasso.diffs[i] <- lasso.coord$lasso.start[i] - lasso.coord$lasso.end[i-1]
        }
        
        dmr.marker <- numeric(length(lasso.diffs))
        for (i in 2:length(dmr.marker))
        {
            dmr.marker[i] <- ifelse(tmp3$CHR[i] != tmp3$CHR[i-1] | lasso.diffs[i] > minDmrSep, 0, 1)
        }
	
        dmr.no <- numeric(length(dmr.marker))
        dmr.no[1] <- 1
        for (i in 2:length(dmr.marker))
        {
            dmr.no[i] <- ifelse(dmr.marker[i] == 0, dmr.no[i-1]+1, dmr.no[i-1])
        }
	
        dmrBylasso <- data.frame(lasso.coord, dmr.no)
    
        dmrs <- data.frame("chr"=unlist(lapply(split(dmrBylasso$CHR,dmrBylasso$dmr.no),unique)), "dmr.no"=unique(dmrBylasso$dmr.no), "dmr.start"=unlist(lapply(split(dmrBylasso[,2:3], dmrBylasso$dmr.no),min)), "dmr.end"=unlist(lapply(split(dmrBylasso[,2:3], dmrBylasso$dmr.no),max)))
        dmrs$dmr.size <- dmrs$dmr.end - dmrs$dmr.start +1
        dmrs <- dmrs[dmrs$dmr.size >= minDmrSize,]
        rownames(dmrs) <- 1:nrow(dmrs) # renames rows after removing small DMRs
        dmrs$dmr.no <- 1:nrow(dmrs) # renames DMRs after removing small DMRs
        
        dmr.probes <- vector("list", nrow(dmrs)) # collects all probes (sig and nonsig) within DMR boundaries
        names(dmr.probes) <- dmrs$dmr.no
        for (i in 1:length(dmr.probes))
        {
            dmr.probes[[i]] <- which(myResults$CHR == dmrs$chr[i] & myResults$MAPINFO >= dmrs$dmr.start[i] & myResults$MAPINFO <= dmrs$dmr.end[i])
        }
        dmr.probes.lengths <- unlist(lapply(dmr.probes, length))
        dmr.probes <- data.frame("dmr.no"=rep(as.numeric(names(dmr.probes)), dmr.probes.lengths), "probe.no"=unlist(dmr.probes))
        
        tmp4 <- data.frame("ID" = rownames(myResults)[dmr.probes$probe.no], myResults[dmr.probes$probe.no,], dmrs[dmr.probes$dmr.no,-1], row.names=1:nrow(dmr.probes))
        
        if(filterXY) # retro-converts $CHR into a factor
        {
            tmp4$CHR <- as.factor(as.character(tmp4$CHR))
        }else{
            tmp4$CHR <- as.character(tmp4$CHR)
            tmp4$CHR <- replace(tmp4$CHR, which(tmp4$CHR == 23), "X")
            tmp4$CHR <- replace(tmp4$CHR, which(tmp4$CHR == 24), "Y")
            tmp4$CHR <- as.factor(tmp4$CHR)
        }

	if(is.null(beta.norm))
	{
		message("beta values missing...cannot compute dmr.p")
	}else{	
    
        dmr.beta <- split(as.data.frame(beta.norm[match(tmp4$ID, rownames(beta.norm)),]), tmp4$dmr.no)
		corel <- lapply(dmr.beta, function(x) cor(t(x)))
		weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
		dmr.ind.p <- split(tmp4$adj.P.Val, tmp4$dmr.no)
		dmr.qp <- lapply(dmr.ind.p, qnorm)
		dmr.qp.w <- mapply("*", dmr.qp, weights)
		dmr.stat <- lapply(dmr.qp.w, sum)
		dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
		dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
		dmr.rep <- as.numeric(summary(dmr.ind.p)[,1])
		dmr.p <- dmr.p*length(dmr.p)/rank(dmr.p)
        dmr.p <- rep(dmr.p, dmr.rep)
      }
    
        if(image)
        {
            plotThreshold(myResults,lasso.radii,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
        }
    
        message("You have found ",max(dmr.probes$dmr.no)," significant DMRs.")
        resultsFile=resultsFile[c(1,30)]
        colnames(resultsFile)[1]="ID"
        tmp4=merge(resultsFile,tmp4,by="ID")
 
        tmp4 <- tmp4[order(tmp4$dmr.no),]
	if(!is.null(dmr.p))
	{
		
        dmrList=data.frame(tmp4,"dmr.p" = dmr.p)
        dmrList=dmrList[which(dmrList$dmr.p < DMRpval),]
	}
        colnames(dmrList)[1]="probeID"

        fileName=paste(resultsDir,"/DMR_",DMRpval,"_",max(dmrList$dmr.no),".txt",sep="")
        write.table(dmrList, fileName ,quote=F,sep="\t",row.names=F)
	
        ########make bedfile
        if(bedFile)
        {
            bedfile=champ.bedfile(dmrList)
            #add pvalue to file name
            fileName1=paste(resultsDir,"/DMR_",DMRpval,"_",max(dmrList$dmr.no),".bed",sep="")
            write.table(bedfile,fileName1,row.names=F,col.names=F,quote=F, sep = "\t")
        }

    }
    return(dmrList)
    
}
