if(getRversion() >= "3.1.0") utils::globalVariables(c("myDMP","myDMR","myNorm","PathwayList"))

champ.GSEA <- function(beta=myNorm,
                       DMP=myDMP[[1]],
                       DMR=myDMR,
                       CpGlist=NULL,
                       Genelist=NULL,
                       pheno=myLoad$pd$Sample_Group,
                       method="fisher",
                       arraytype="450K",
                       Rplot=TRUE,
                       adjPval=0.05)
{
    message("[===========================]")
    message("[<<<< ChAMP.GSEA START >>>>>]")
    message("-----------------------------")

    ##### Defind a inter function #####
    PasteVector <- function(v){
        vt <- v[1];
        if(length(v) > 1){
            for(g in 2:length(v)){
                vt <- paste(vt,v[g],sep=" ")
            }
        }
        vt <- paste(vt," EnD",sep="");
        out.v <- sub(" EnD","",vt);
        out.v <- sub("NA , ","",out.v);
        out.v <- sub(" , NA","",out.v);
        out.v <- sub(" , NA , "," , ",out.v);
        return(out.v);
    }


    if(method=="fisher")
    {
        data(PathwayList)
        if(arraytype=="EPIC"){
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
        }else{
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
        }
        RSanno <- getAnnotation(RSobject)[,c("chr","pos","Name","UCSC_RefGene_Name")]
        ueid.v <- unique(unlist(sapply(RSanno$UCSC_RefGene_Name,function(x) strsplit(x,split=";")[[1]])))
        bias.Data <- table(unlist(sapply(RSanno$UCSC_RefGene_Name,function(x) unique(strsplit(x,split=";")[[1]]))))
        bias.Data <- as.numeric(bias.Data[ueid.v])

        loi.lv <- list()
        if(!is.null(DMP))
        {
            cpg.idx <- rownames(DMP)[which(DMP$adj.P.Val <= adjPval)]
            loi.lv[["DMP"]] <- unique(unlist(sapply(RSanno[cpg.idx,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]])))
        }
        if(!is.null(DMR))
        {
            #DMR[[1]]$seqnames <- as.character(DMR[[1]]$seqnames)
            cpg.idx <- unique(unlist(apply(DMR[[1]],1,function(x) rownames(RSanno)[which(RSanno$chr==x[1] & RSanno$pos >= as.numeric(x[2]) & RSanno$pos <= as.numeric(x[3]))])))
            loi.lv[["DMR"]] <- unique(unlist(sapply(RSanno[cpg.idx,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]])))
        }
        if(!is.null(CpGlist))
        {
            if(class(CpGlist)=="character")
            {
                loi.lv[["CpG"]] <- unique(unlist(sapply(RSanno[CpGlist,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]])))
            }else if(class(CpGlist)=="list")
            {
                if(names(CpGlist) %in% c("DMP","DMR","CpG","Gene")) stop("Your CpG list can not contain names as follows: \"DMP\",\"DMR\",\"CpG\",\"Gene\".")
                if(names(CpGlist)==NULL | sum(is.na(names(CpGlist)))!=0) stop("Please specify your CpG list names.")
                cpg.tmp <- lapply(CpGlist,function(x) unique(unlist(sapply(RSanno[x,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]]))))
                loi.lv <- c(loi.lv,cpg.tmp)
            }else
            {
                stop("CpGlist Format is not correct, it must be character vector or character list.")
            }
        }
        if(!is.null(Genelist))
        {
            if(class(Genelist)=="character")
            {
                loi.lv[["Gene"]] <- Genelist
            }else if(class(Genelist)=="list")
            {
                if(names(Genelist) %in% c("DMP","DMR","CpG","Gene")) stop("Your Gene list can not contain names as follows: \"DMP\",\"DMR\",\"CpG\",\"Gene\".")
                if(names(Genelist)==NULL | sum(is.na(names(Genelist)))!=0) stop("Please specify your Gene list names.")
                loi.lv <- c(loi.lv,Genelist)
            }else
            {
                stop("Genelist Format is not correct, it must be character vector or character list.")
            }
        }
        message("<< Prepare Gene List Ready  >>")
        
        ######################  Do GSEA   ###########################

        selGEID.lv <- loi.lv;
        
        ### How many of the MSigDB list genes are on array?
        listnames.v <- names(PathwayList);
        nrepGEID.m <- matrix(nrow=length(listnames.v),ncol=3);
        rownames(nrepGEID.m) <- listnames.v;
        colnames(nrepGEID.m) <- c("nList","nRep","fRep");

        nrepGEID.m[,"nList"] <- unlist(lapply(PathwayList,function(x) length(x)))
        nrepGEID.m[,"nRep"] <- unlist(lapply(PathwayList,function(x) length(intersect(x,ueid.v))))
        nrepGEID.m[,"fRep"] <- nrepGEID.m[,"nRep"]/nrepGEID.m[,"nList"]


        ### Do enrichment analysis
        listsummary.lm <- list();

        message("<< Start Do GSEA on each Gene List  >>")

        for(i in 1:length(selGEID.lv))
        {
            message("<< Do GSEA on Gene list ",names(loi.lv)[i],">>")
            listPV.v <- vector();
            listOR.v <- vector();
            fisher.lm <- list();
            novlap.v <- vector();
            ovlapG.lv <- list();

            selEID.v <- selGEID.lv[[i]];

            # Blow is pale fisher exact test method.
            # ==================
            if(method == "fisher")
            {
                message("<< Pale Fisher Exact Test will be used to do GSEA >>")
                message(" << The category information is downloaded from MsigDB, and only simple Fisher Exact Test will be used to calculate GSEA. This method is suitable if your genes has equalivalent probability to be enriched. If you are using CpGs mapping genes, gometh method is recommended.>> ")
                InMSigDBlist_InList <- unlist(lapply(PathwayList,function(x) length(intersect(x,selEID.v))))
                NotInMSigDBlist_InList <- length(selEID.v)-InMSigDBlist_InList
                InMSigDBlist_NotInList <- unlist(lapply(PathwayList,function(x) length(intersect(x,ueid.v))-length(intersect(x,selEID.v)))) 
                NotInMSigDBlist_NotInList <- length(ueid.v)-InMSigDBlist_InList-NotInMSigDBlist_InList-InMSigDBlist_NotInList 

                ovlapG.lv_2 <- lapply(PathwayList,function(x) intersect(x,selEID.v))
                fisher.lm_2 <- cbind(InMSigDBlist_InList,NotInMSigDBlist_InList,InMSigDBlist_NotInList,NotInMSigDBlist_NotInList) 
                listPV.v_2 <- t(apply(fisher.lm_2,1,function(x) unlist(fisher.test(matrix(x,2,2),alternative="greater")[c(1,3)])))
                info <- cbind(fisher.lm_2,listPV.v_2)

                selL.idx <- which(nrepGEID.m[,3]>0.6); ### remove lists with less than 60% representation on array
                listsummary.m_2 <- data.frame("Gene_List"=listnames.v[selL.idx],
                                              nrepGEID.m[selL.idx,],
                                              "nOVLAP"=fisher.lm_2[,1][selL.idx],
                                              "OR"=listPV.v_2[selL.idx,2],
                                              "P-value"=listPV.v_2[selL.idx,1],
                                              "adjPval"=p.adjust(listPV.v_2[selL.idx,1],method="BH"),
                                              "Genes"=unlist(lapply(ovlapG.lv_2[selL.idx],function(x) PasteVector(x))))
                listsummary.m_2 <- listsummary.m_2[order(listsummary.m_2$adjPval),]
                pvals <- listsummary.m_2[which(listsummary.m_2$adjPval <= adjPval),]
            }else if(method == "goseq")
            {
                message("goseq method will be used to do GSEA")
                message("This method is supported by goseq package, which is developed to address inequalivalent issue between number of CpGs and genes.")
                # Below is goseq method
                # goseq method is incoperated here because not all genes contain same length of CpGs.
                # goseq would generated a probability weight function between significant genes and 
                # number of CpGs contained by each genes, one plot will be plotted on that.
                # Then goseq function would use corrected gene list to do GSEA.
                # ===================
                DE.genes <- as.integer(ueid.v %in% selEID.v)
                names(DE.genes) <- ueid.v
                pwf <- nullp(DE.genes,"hg19","geneSymbol",bias.data=bias.Data,plot.fit=Rplot)
                pvals <- goseq(pwf,'hg19','geneSymbol')
                over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
                pvals <- cbind(pvals,over_represented_adjPvalue)
            }
            listsummary.lm[[i]] <- pvals;
            message("<< Done for Gene list ",names(selGEID.lv)[i]," >>");
        }
        names(listsummary.lm) <- names(selGEID.lv);
    } else if(method=="gometh")
    {
        if(arraytype=="EPIC"){
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
        }else{
            RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
        }

        RSanno <- getAnnotation(RSobject)[,c("chr","pos","Name","UCSC_RefGene_Name")]
        ueid.v <- rownames(beta)
        loi.lv <- list()
        if(!is.null(DMP))
        {
            cpg.idx <- rownames(DMP)
            loi.lv[["DMP"]] <- cpg.idx
        }
        if(!is.null(DMR))
        {
            cpg.idx <- unique(unlist(apply(DMR[[1]],1,function(x) rownames(RSanno)[which(RSanno$chr==x[1] & RSanno$pos >= as.numeric(x[2]) & RSanno$pos <= as.numeric(x[3]))])))
            loi.lv[["DMR"]] <- cpg.idx
        }
        if(!is.null(CpGlist))
        {   
            loi.lv[["CpG"]] <- CpGlist
        }
        if(!is.null(Genelist))
        {
            message("gometh function only take CpGs to calculate GSEA, so your Gene List can not be used in this function. Please use \"fisher\"method.")
        }

        message("<< Prepare CpG List Ready  >>")
        listsummary.lm <- list()
        for(i in names(loi.lv))
        {
            message("  Calculating GSEA with gometh method on ",i," CpG list")
            message("  Note that gometh method would count numbers of CpGs in each genes and correct this bias.")
            tmp <- gometh(loi.lv[[i]], all.cpg=ueid.v, collection=c("GO"), array.type = arraytype, plot.bias=Rplot, prior.prob=TRUE)
            A <- tmp[tmp$FDR <= adjPval,]
            listsummary.lm[[i]] <- A[order(A$FDR),]
        }
    } else if (method=="ebayes") {
        message("<< ebayes GSEA method  >>")
        message("ebayes GSEA method would do GSEA based on global test on all CpG sites, thus you don't need to provide any DMP or DMR here, but instead, please assign pheno parameter.")
        if(is.null(pheno)) stop("Pheno parameter is required for ebayes GSEA method.")

        if(class(pheno)=="numeric") {
            message("  pheno parameter is continues variable, linear method would be used when calculating ebayes GSEA.")
        } else if (class(pheno) %in% c("character","factor")) {
            if(length(table(pheno) == 2)) message("  pheno parameter is category variable, logistic method would be ued when calculating ebayes GSEA.")
            else stop("If your pheno parameter is category variable, only two phenotypes are allowed, say Cancer vs Normal, please modify your data set.")
        }
        listsummary.lm <- champ.ebGSEA(beta=beta,pheno=pheno,adjPval=adjPval,minN=5,arraytype=arraytype)
    } else {
        stop(" You must assign method parameter as fisher, gometh or ebayes.")
    }

    message("[<<<<< ChAMP.GSEA END >>>>>>]")
    message("[===========================]")
    return(listsummary.lm)
} 
