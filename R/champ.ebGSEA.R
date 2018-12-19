if(getRversion() >= "3.1.0") utils::globalVariables("PathwayList")

champ.ebGSEA <- function(beta=myNorm, pheno=myLoad$pd$Sample_Group, minN=5, adjPval=0.05, arraytype="450K")
{
    #library(Hmisc);
    #library(globaltest);
  message("ebGSEA requires no NA in your beta and pheno parameter.")

    mapEIDtoCpG <- function(beta,arraytype) {
        if(arraytype=="EPIC"){
            message("  Extracting annotation from IlluminaHumanMethylationEPICilm10b2.hg19.")
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
        }else{
            message("  Extracting annotation from IlluminaHumanMethylation450kilmn12.hg19.")
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
        }
        RSanno <- as.data.frame(getAnnotation(RSobject))
        message("  Removing Non-CG probes out of annotation.")
        ann.keep <- RSanno[substr(rownames(RSanno),1,2)=="cg" & RSanno$UCSC_RefGene_Name != "",]

        message("  Flat all genes on each CpG.")
        geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
        names(geneslist)<-rownames(ann.keep)

        grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
        names(grouplist)<-rownames(ann.keep)

        flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
        flat$symbol<-as.character(flat$symbol)
        flat$group <- as.character(flat$group)

        flat$cpg<- substr(rownames(flat),1,10)

        flat$alias <- alias2SymbolTable(flat$symbol)

        #eg <- toTable(org.Hs.egSYMBOL2EG)
        #m <- match(flat$alias,eg$symbol)
        #flat$entrezid <- eg$gene_id[m]        #transform symbol to entrezid
        #flat <- flat[!is.na(flat$entrezid),]  #get rid of NA entrezid 

        message("  Removing all duplicated CpG-genes.")
        id<-paste(flat$cpg,flat$alias,sep=".")
        d <- duplicated(id)
        flat.u <- flat[!d,]

        message("  Annotation Prepared.")
        mapEIDtoCpG <- split(flat.u$cpg, flat.u$alias)
        return(mapEIDtoCpG)
    }

    gseaWTfn <- function(termEID.v,rankEID.v,minN=minN){
        commonEID.v <- intersect(termEID.v,rankEID.v);
        nrep <- length(commonEID.v);
        if(length(commonEID.v)>=minN){
            otherEID.v <- setdiff(rankEID.v,termEID.v);
            match(commonEID.v,rankEID.v) -> rank1.idx;
            match(otherEID.v,rankEID.v) -> rank2.idx;
            wilcox <- wilcox.test(rank1.idx,rank2.idx,alt="less")
            pv <- wilcox$p.value
            n1n2 <- nrep * length(otherEID.v)
            auc <- (n1n2-wilcox$statistic)/(n1n2)

            ## add kpmt here.
            pop.v <- 1:length(rankEID.v);
            names(pop.v) <- rankEID.v;
            obs.v <- commonEID.v
            pvKPMT <- kpmt::kpmt(pop=pop.v,obs=obs.v,tail="lower")[[4]]

            out.v <- c(nrep,auc,pv,pvKPMT);
        } else {
            out.v <- c(nrep,0,1,1);
        }
        return(out.v);
    }


    message("\n[ Section 1:  Generate Annotation Start ]\n")
    mapEID <- mapEIDtoCpG(beta,arraytype)
    message("\n[ Section 1:  Generate Annotation Done ]\n")


    message("\n[ Section 2:  Running Global Test Start ]\n")
    if(class(pheno)=="numeric")
    {
        message("  Applying Linear Model on Global Test. It could be very slow...")
        gt.o <- gt(response=pheno, alternative=t(beta),model="linear",directional = FALSE,
                   standardize = FALSE, permutations = 0, subsets=mapEID,trace=FALSE);
    } else {
        message("  Applying Binary Model on Global Test. It could be very slow...")
        gt.o <- gt(response=(as.numeric(as.factor(pheno))-1), alternative=t(beta),model="logistic",directional = FALSE,
                   standardize = FALSE, permutations = 0, subsets=mapEID,trace=FALSE);
    }
    resGT.m <- result(gt.o);
    tmp.s <- sort(resGT.m[,1],index.return=TRUE);
    sresGT.m <- resGT.m[tmp.s$ix,];
    message("\n[ Section 2:  Running Global Test Done ]\n")

    message("\n[ Section 3:  GSEA on Pathway Start ]\n")
    message("  Loading MsigDB PathwayList information.")
    data(PathwayList)

    message("  Doing wilcox test and Known Population Median Test, it could be slow here.")
    gseaWT.m <- do.call(rbind,lapply(PathwayList,function(x) gseaWTfn(termEID.v=x, rankEID.v=rownames(sresGT.m),minN=minN)))
    colnames(gseaWT.m) <- c("nREP","AUC","P","P(KPMT)")
    rownames(gseaWT.m) <- names(PathwayList);

    message("  Adjusting Pathway P value with BH method.")
    tmp.s <- sort(gseaWT.m[,3],decreasing=FALSE,index.return=TRUE);
    sgseaWT.m <- gseaWT.m[tmp.s$ix,];
    padj.v <- p.adjust(sgseaWT.m[,3],method="BH");
    message("  Your adjPval is ",adjPval, ", only pathway below BH adjusted P value 0.05 would be returned.")
    sel.idx <- which(padj.v <= adjPval);

    message("  Forming up final result.")
    if(length(sel.idx) > 1) {
        topGSEAwt.m <- cbind(sgseaWT.m[sel.idx,],padj.v[sel.idx]);
        colnames(topGSEAwt.m) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");

        topGSEAwt.lm <- list();
        topGSEAwt.lm[[1]] <- topGSEAwt.m;
        tmp.s <- sort(sgseaWT.m[,2],decreasing=TRUE,index.return=TRUE);
        topGSEAwt.lm[[2]] <- sgseaWT.m[tmp.s$ix,];
        names(topGSEAwt.lm) <- c("Rank(P)","Rank(AUC)");
    } else if (length(sel.idx)==1) {
        topGSEAwt.v <- as.vector(c(sgseaWT.m[sel.idx,],padj.v[sel.idx]));
        names(topGSEAwt.v) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");
        topGSEAwt.lm <- list("Rank(P)"=topGSEAwt.v,"Rank(AUC)"=topGSEAwt.v,"POI"=rownames(sgseaWT.m)[sel.idx]);
    } else {
        message("No significant pathway detected. You may try relax the threshold like adjPval.")
        topGSEAwt.lm <- list();
    }

    message("\n[ Section 3:  GSEA on Pathway Done ]\n")

    return(topGSEAwt.lm)
} 
