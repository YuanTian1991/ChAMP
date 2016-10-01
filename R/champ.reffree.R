if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad"))

champ.reffree <- function(beta=myNorm,
                          pheno=myLoad$pd$Sample_Group,
                          K=NULL,
                          nboot=50)
{
    message("[===========================]")
    message("[<<< ChAMP.REFFREE START >>>]")
    message("-----------------------------")

    if(is.null(K))
    {
        tmp.m <- na.omit(beta - rowMeans(beta));
        rmt.o <- EstDimRMT(tmp.m);
        k <- rmt.o$dim;
    }else{
        k <- K
    }

    message("<< Measure numbers of latent variables success >>")
    message("champ.reffree will proceed with ",k," components.")

    if(!is.numeric(pheno))
    {
        if(class(pheno)=="matrix"){
            pheno <- apply(pheno,2,function(x) (as.numeric(as.factor(x))-1))
            denDf <- dim(pheno)[1]
        }else{
            pheno <- (as.numeric(as.factor(pheno))-1)
            denDf <- length(pheno)
        }
    }

    if(nboot <= 0)
    {
        message("nboot parameter must be a positive value.")
        return
    }

    rf.o <- RefFreeEwasModel(beta,cbind(1,pheno),K=k)
    rfB.o <- BootRefFreeEwasModel(rf.o,nboot);

    message("<< Calculate RefFreeEWASModel Success >>")

    seBeta <- apply(rfB.o[,,"B",], 1:2, sd)
    seBstar <- apply(rfB.o[,,"B*",], 1:2, sd)

    pvBeta <- 2*pt(abs(rf.o$Beta)/seBeta,denDf,lower.tail=FALSE)
#    pvBstar <- 2*pt(abs(rf.o$Bstar)/seBstar,denDf,lower.tail=FALSE)
    qvBeta <- apply(pvBeta,2,function(x) qvalue(x)$qvalue)
#    qvBstar <- apply(pvBstar,2,function(x) qvalue(x)$qvalue)

    message("Generate p value and q value success.")
    

    message("[<<<< ChAMP.REFBASE END >>>>]")
    message("[===========================]")
    return(list(RefFreeEWASModel=rf.o,pvBeta=pvBeta,qvBeta=qvBeta))
}
