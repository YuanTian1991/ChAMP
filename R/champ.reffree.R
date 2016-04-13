if(getRversion() >= "3.1.0") utils::globalVariables("myLoad")

champ.reffree <- function(beta=myLoad$beta,design=myLoad$pd$Sample_Group,K=NULL,nboot=10)
{
    if(is.null(K))
    {
        tmp.m <- na.omit(beta - rowMeans(beta));
        rmt.o <- EstDimRMT(tmp.m);
        k <- rmt.o$dim;
    }else{
        k <- K
    }

    if(!is.numeric(design))
    {
        if(class(design)=="matrix"){
            design <- apply(design,2,function(x) (as.numeric(as.factor(x))-1))
            denDf <- dim(design)[1]
        }else{
            design <- (as.numeric(as.factor(design))-1)
            denDf <- length(design)
        }
    }

    if(nboot <= 0)
    {
        message("nboot parameter must be a positive value.")
        return
    }

    rf.o <- RefFreeEwasModel(beta,cbind(1,design),K=k)
    rfB.o <- BootRefFreeEwasModel(rf.o,nboot);

    seBeta <- apply(rfB.o[,,"B",], 1:2, sd)
    seBstar <- apply(rfB.o[,,"B*",], 1:2, sd)

    pvBeta <- 2*pt(abs(rf.o$Beta)/seBeta,denDf,lower.tail=FALSE)
#    pvBstar <- 2*pt(abs(rf.o$Bstar)/seBstar,denDf,lower.tail=FALSE)
    qvBeta <- apply(pvBeta,2,function(x) qvalue(x)$qvalue)
#    qvBstar <- apply(pvBstar,2,function(x) qvalue(x)$qvalue)
    
    return(list(RefFreeEWASModel=rf.o,pvBeta=pvBeta,qvBeta=qvBeta))
}
