if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","probe.features.epic","probe.features"))

champ.MVP <- function(beta = myNorm,
                      pheno = myLoad$pd$Sample_Group,
                      adjPVal = 0.05,
                      adjust.method = "BH",
                      compare.group = NULL,
                      arraytype = "450K")
{
    message("[===========================]")
    message("[<<<<< ChAMP.MVP START >>>>>]")
    message("-----------------------------")


    if(is.null(pheno) | length(unique(pheno))<=1)
    {
        stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
    }else
    {
        message("<< Your pheno information contains following groups. >>")
        sapply(unique(pheno),function(x) message("<",x,">:",sum(pheno==x)," samples."))
        message("[The power of statistics analysis on groups contain very few samples may not strong.]")
    }
	
    if(is.null(compare.group))
    {
        message("You did not assign compare groups. The first two groups: <",unique(pheno)[1],"> and <",unique(pheno)[2],">, will be compared automatically.")
        compare.group <- unique(pheno)[1:2]
    }else if(sum(compare.group %in% unique(pheno))==2)
    {
        message("As you assigned, champ.svd will compare ",compare.group[1]," and ",compare.group[2],".")
    }else
    {
        message("Seems you did not assign correst compare groups. The first two groups: <",unique(pheno)[1],"> and <",unique(pheno)[2],">, will be compared automatically.")
        compare.group <- unique(pheno)[1:2]
    }
    pheno <- pheno[which(pheno %in% compare.group)]
    beta <- beta[,which(pheno %in% compare.group)]
	design <- model.matrix( ~ 0 + pheno)
    contrast.matrix <- makeContrasts(contrasts=paste(colnames(design)[2:1],collapse="-"), levels=colnames(design))
    message("Test")
    message("\n<< Contrast Matrix >>")
    print(contrast.matrix)

    message("\n<< All beta, pheno and model are prepared successfully. >>")
	
	fit <- lmFit(beta, design)
	fit2 <- contrasts.fit(fit,contrast.matrix)
	tryCatch(fit3 <- eBayes(fit2),
      warning=function(w) 
      {
      	stop("limma failed, No sample variance.\n")
      })
    MVP <- topTable(fit3,coef=1,number=nrow(beta),adjust.method=adjust.method,p.value=adjPVal)
    message("You have found ",sum(MVP$adj.P.Val <= adjPVal), " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
    message("\n<< Calculate MVP successfully. >>")

    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    com.idx <- intersect(rownames(MVP),rownames(probe.features))
    avg <-  cbind(rowMeans(beta[com.idx,which(pheno==compare.group[1])]),rowMeans(beta[com.idx,which(pheno==compare.group[2])]))
    avg <- cbind(avg,avg[,2]-avg[,1])
    colnames(avg) <- c(paste(compare.group,"AVG",sep="_"),"deltaBeta")
    MVP <- data.frame(MVP[com.idx,],avg,probe.features[com.idx,])

    message("[<<<<<< ChAMP.MVP END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.DMR() next.]\n")
    return(MVP)
}
