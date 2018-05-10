if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad"))

champ.Block <- function(beta=myNorm,
                        pheno=myLoad$pd$Sample_Group,
                        arraytype="450K",
                        maxClusterGap=250000,
                        B=500,
                        bpSpan=250000,
                        minNum=10,
                        cores=3)
{
    message("[===========================]")
    message("[<<<< ChAMP.Block START >>>>]")
    message("-----------------------------")


    #############  ConstResgion.R  ################
    if(arraytype=="EPIC"){
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19"))
    }else{
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
    }
    if(cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores = cores)

    RSanno <- getAnnotation(RSobject)[,c("chr","pos","Relation_to_Island")]
    RSanno <- RSanno[order(RSanno$chr,RSanno$pos),]
    sbeta.m <- beta[rownames(RSanno),]
    sregion.v <- RSanno$Relation_to_Island

    message("<< Load Annotation Successfully >>")

    ### define open sea clusters
    opensea.idx <- which(sregion.v=="OpenSea")
    openseaCLID.v <- boundedClusterMaker(RSanno$chr[opensea.idx],RSanno$pos[opensea.idx],assumeSorted=TRUE,maxGap=500,maxClusterWidth=1500)
    names(openseaCLID.v) <- rownames(RSanno)[opensea.idx]

    ### define island clusters
    island.idx <- which(sregion.v=="Island")
    cpgiCLID.v <- boundedClusterMaker(RSanno$chr[island.idx],RSanno$pos[island.idx],assumeSorted=TRUE,maxGap=300,maxClusterWidth=1500)
    names(cpgiCLID.v) <- rownames(RSanno)[island.idx]

    ### define shore/shelve clusters
    sh.idx <- setdiff(1:length(sregion.v),c(opensea.idx,island.idx))
    shCLID.v <- boundedClusterMaker(RSanno$chr[sh.idx],RSanno$pos[sh.idx],assumeSorted=TRUE,maxGap=500,maxClusterWidth=1500)
    names(shCLID.v) <- rownames(RSanno)[sh.idx]

    allCLID.v <- rep(NA,length(sregion.v))
    allCLID.v[opensea.idx] <- paste("OS",openseaCLID.v,sep="")
    allCLID.v[island.idx] <- paste("CPGI",cpgiCLID.v,sep="")
    allCLID.v[sh.idx] <- paste("SH",shCLID.v,sep="")
    names(allCLID.v) <- rownames(RSanno)

    nf <- length(unique(allCLID.v))
    npg.v <- summary(factor(allCLID.v),maxsum=nf)

    message("<< Get Clusters by cgi-info Successfully >>")

    ### average profiles
    avbetaCL.m <- rowsum(sbeta.m,allCLID.v,reorder=FALSE)
    nf.v <- npg.v[match(rownames(avbetaCL.m),names(npg.v))]
    avbetaCL.m <- avbetaCL.m/nf.v
    rownames(sbeta.m) -> sortedCpGsbeta.v;
    rm(sbeta.m);

    message("<< Calculate Average Beta Value Successfully >>")
    ### find midpoint, start, end pos of each region
    tmpPOS.v <- RSanno$pos[match(names(allCLID.v),rownames(RSanno))]
    tmpCHR.v <- substr(RSanno$chr[match(names(allCLID.v),rownames(RSanno))],4,100)
    tmpCHR.v[which(tmpCHR.v=="X")] <- "23"
    tmpCHR.v[which(tmpCHR.v=="Y")] <- "24"
    tmpCHR.v <- as.numeric(tmpCHR.v)
    
    posCL.m <- rowsum(cbind(tmpPOS.v,tmpCHR.v),allCLID.v,reorder=FALSE);
    nf.v <- npg.v[match(rownames(posCL.m),names(npg.v))]
    posCL.m <- posCL.m/nf.v;
    colnames(posCL.m) <- c("Pos","Chr");

    message("<< Generate Block Position Successfully >>")
    os.idx <- grep("OS",rownames(avbetaCL.m));
    avbetaOS.m <- avbetaCL.m[os.idx,];
    posOS.m <- posCL.m[os.idx,];

    design.m <- data.frame(1,Sample_Group=as.numeric(as.factor(pheno))-1)
    selAUT.idx <- which(posOS.m[,2]<=22);
    blocks <- clusterMaker(chr=posOS.m[selAUT.idx,2],pos=posOS.m[selAUT.idx,1],assumeSorted = TRUE, maxGap = maxClusterGap);
    blockID.v <- levels(as.factor(blocks));
    message("<< New Clusters are generated for blocks >>")

    blockPROP.m <- data.frame(CHR=aggregate(posOS.m[,2],by=list(blocks),function(x) x[1])[,2],
                              aggregate(posOS.m[,1],by=list(blocks),function(x) c(mean(x),min(x),max(x),max(x)-min(x),length(x)))[,2])
    colnames(blockPROP.m) <- c("CHR","AvPos","Start","End","Size","NumberOS");

    message("<< Generate information for New Clusters >>")

    bh2.o <- bumphunter(avbetaOS.m[selAUT.idx,],design.m,chr=posOS.m[selAUT.idx,2],pos=posOS.m[selAUT.idx,1],cluster=blocks,coef=2, pickCutoff=TRUE, pickCutoffQ=0.75,smooth=TRUE,smoothFunction=loessByCluster, useWeights=FALSE, B=B,bpSpan=bpSpan,minNum=minNum);
    message("<< Run Bumphunter Successfully >>")

    Block <- bh2.o$tab
    rownames(Block) <- paste("Block",1:nrow(Block),sep="_")


    message("[<<<<< ChAMP.BLOCK END >>>>>]")
    message("[===========================]")
    message("[You may want to process Block.GUI() next.]\n")
    return(list(Block=Block,clusterInfo=blockPROP.m,allCLID.v=allCLID.v,avbetaCL.m=avbetaCL.m,posCL.m=posCL.m))
}
