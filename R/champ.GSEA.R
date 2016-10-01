if(getRversion() >= "3.1.0") utils::globalVariables(c("myMVP","myDMR","myNorm","PathwayList"))

champ.GSEA <- function(beta=myNorm,
                       MVP=myMVP,
                       DMR=myDMR,
                       CpGlist=NULL,
                       Genelist=NULL,
                       arraytype="450K",
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
    data(PathwayList)

    if(arraytype=="EPIC"){
         RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b2.hg19"))
    }else{
        RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylation450k",annotation = "ilmn12.hg19"))
    }
    RSanno <- getAnnotation(RSobject)[,c("chr","pos","Name","UCSC_RefGene_Name")]
    ueid.v <- unique(unlist(sapply(RSanno$UCSC_RefGene_Name,function(x) strsplit(x,split=";")[[1]])))

    loi.lv <- list()
    if(!is.null(MVP))
    {
        cpg.idx <- rownames(MVP)[which(MVP$adj.P.Val <= adjPval)]
        loi.lv[["MVP"]] <- unique(unlist(sapply(RSanno[cpg.idx,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]])))
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
            loi.lv[["CpG"]] <- unique(unlist(sapply(RSanno[cpg.idx,"UCSC_RefGene_Name"],function(x) strsplit(x,split=";")[[1]])))
        }else if(class(CpGlist)=="list")
        {
            if(names(CpGlist) %in% c("MVP","DMR","CpG","Gene")) stop("Your CpG list can not contain names as follows: \"MVP\",\"DMR\",\"CpG\",\"Gene\".")
            if(names(CpGlist)==NULL | length(is.na(names(CpGlist)))!=0) stop("Please specify your CpG list names.")
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
            if(names(Genelist) %in% c("MVP","DMR","CpG","Gene")) stop("Your Gene list can not contain names as follows: \"MVP\",\"DMR\",\"CpG\",\"Gene\".")
            if(names(Genelist)==NULL | length(is.na(names(Genelist)))!=0) stop("Please specify your Gene list names.")
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

    message("<< Do GSEA on each Gene List  >>")

    for(i in 1:length(selGEID.lv))
    {
        message("<< Do GSEA on Gene list ",names(loi.lv)[i],">>")
        listPV.v <- vector();
        listOR.v <- vector();
        fisher.lm <- list();
        novlap.v <- vector();
        ovlapG.lv <- list();

        selEID.v <- selGEID.lv[[i]];

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
        listsummary.m_2 <- listsummary.m_2[which(listsummary.m_2$adjPval <= adjPval),]
        listsummary.lm[[i]] <- listsummary.m_2;

        message("<< Done for Gene list ",names(selGEID.lv)[i]," >>");
    }
    names(listsummary.lm) <- names(selGEID.lv);

    message("[<<<<< ChAMP.GSEA END >>>>>>]")
    message("[===========================]")
    return(list(GSEA=listsummary.lm,GeneList=loi.lv))
} 
