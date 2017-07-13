if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","probe.features.epic","probe.features"))

champ.DMP <- function(beta = myNorm,
                      pheno = myLoad$pd$Sample_Group,
                      compare.group = NULL,
                      adjPVal = 0.05,
                      adjust.method = "BH",
                      arraytype = "450K")
{
    message("[===========================]")
    message("[<<<<< ChAMP.DMP START >>>>>]")
    message("-----------------------------")

    CalculateDMP <- function(beta,pheno,tmp_compare,adjPVal=adjPVal,adjust.method=adjust.method)
    {
        message("  -----------------------------")
        message("  Start to Compare : ",tmp_compare[1],", ",tmp_compare[2])
        p <- pheno[which(pheno %in% tmp_compare)]
        tmpbeta <- beta[,which(pheno %in% tmp_compare)]
        design <- model.matrix( ~ 0 + p)
        contrast.matrix <- makeContrasts(contrasts=paste(colnames(design)[2:1],collapse="-"), levels=colnames(design))
        message("  Contrast Matrix")
        print(contrast.matrix)
        fit <- lmFit(tmpbeta, design)
        fit2 <- contrasts.fit(fit,contrast.matrix)
        tryCatch(fit3 <- eBayes(fit2),
                 warning=function(w) 
                {stop("limma failed, No sample variance.\n")})
        DMP <- topTable(fit3,coef=1,number=nrow(tmpbeta),adjust.method=adjust.method,p.value=adjPVal)
        message("  You have found ",sum(DMP$adj.P.Val <= adjPVal), " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
        message("  Calculate DMP for ",tmp_compare[1]," and ",tmp_compare[2]," done.")
        return(DMP)
    }

    message("!!! Important !!! New Modification has been made on champ.DMP(): \n")
    message("    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.\n")
    message("    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted \"pheno\" parameter is \"numeric\" type.\n")

    Compare <- NULL

    message("--------------------------------")

    if(is.null(pheno) | length(unique(pheno))<=1)
    {
        stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
    }else
    {
        if(length(pheno)!=ncol(beta)) stop("Your Pheno's length is not in accord with your beta value's ncol.")
        message("\n[ Section 1:  Check Input Pheno Start ]\n")
        if(class(pheno)=="numeric")
        {
            message("  You pheno is numeric type.")
            message("    pheno parameter contains :",length(pheno)," values.")
            message("    pheno parameter ranges from ",min(pheno)," to ",max(pheno))
            
        } else {
            message("  You pheno is ",class(pheno)," type.")
            message("    Your pheno information contains following groups. >>")
            sapply(unique(pheno),function(x) message("    <",x,">:",sum(pheno==x)," samples."))
            message("    [The power of statistics analysis on groups contain very few samples may not strong.]")

            if(length(unique(pheno)) == 2) 
            {
                message("    pheno contains only 2 phenotypes")
                if(is.null(compare.group))
                {
                    message("    compare.group parameter is NULL, two pheno types will be added into Compare List.")
                    Compare <- list(x1=unique(pheno))
                } else if(sum(compare.group %in% unique(pheno))==2)
                {
                    message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
                    Compare <- list(x1=unique(pheno))
                } else {
                    stop(" You have specified compare.group, but it's not in accord with your pheno parameter. Please recheck your compare.group or pheno.")
                }
            } else if(length(unique(pheno)) > 2) {
                message("    pheno contains ",length(unique(pheno))," phenotypes")

                if(is.null(compare.group))
                {
                    message("    compare.group parameter is NULL, EACH PAIR of phenotypes will be added into Compare List.")
                    Compare <- as.list(data.frame(combn(unique(pheno),2)))
                } else if(sum(compare.group %in% unique(pheno))==2) {
                    message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
                    Compare <- list(x1=sort(compare.group))
                } else {
                    stop("    Your pheno parameter contains multiple phenotypes, but values in your compare.group parameter are not all find in them. Please recheck your compare.group or pheno.")
                }
            } else {
                stop("    !!! Something wrong with your pheno. Please check your input.")
            }

            tmpnamelist <- vector()
            for(i in 1:length(Compare))
            {
                tmpname <- paste(Compare[[i]][1],Compare[[i]][2],sep="_to_")
                message("    ",tmpname," compare group : ",Compare[[i]][1],", ",Compare[[i]][2])
                tmpnamelist <- c(tmpnamelist,tmpname)
            }
            names(Compare) <- tmpnamelist
        }
        message("\n[ Section 1:  Check Input Pheno Done ]\n")
    }


    DMPs <- list()
    if(is.null(Compare)){
        message("\n[ Section 2:  Find Numeric Covariates Linear Regression CpGs Start ]\n")
        df <- data.frame(pheno=pheno)
        model.matrix <- model.matrix(~ pheno, data=df)
        fit1 <- lmFit(beta,model.matrix)
        fit2 <- eBayes(fit1)
        DMP <- topTable(fit2,coef=2,number=nrow(beta),adjust.method=adjust.method,p.value=adjPVal)
        message("  You have found ",sum(DMP$adj.P.Val <= adjPVal), " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
        if(sum(DMP$adj.P.Val <= adjPVal)!=0)
            DMPs[["NumericVariable"]] <- DMP
    } else {
        message("\n[ Section 2:  Find Differential Methylated CpGs Start ]\n")
        for(i in names(Compare))
        {
            DMP <- CalculateDMP(beta,pheno,Compare[[i]],adjPVal,adjust.method)
            if(sum(DMP$adj.P.Val <= adjPVal)!=0)
                DMPs[[i]] <- DMP
        }
    }
    message("\n[ Section 2:  Find Numeric Vector Related CpGs Done ]\n")

    if(length(DMPs)==0) stop("ChAMP.DMP Have not detected even one significant CpGs. You may try other threshold.")

    message("\n[ Section 3:  Match Annotation Start ]\n")
    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)

    for(i in names(DMPs))
    {
        com.idx <- intersect(rownames(DMPs[[i]]),rownames(probe.features))
        if(!is.null(Compare))
        {
            avg <-  cbind(rowMeans(beta[com.idx,which(pheno==Compare[[i]][1])]),rowMeans(beta[com.idx,which(pheno==Compare[[i]][2])]))
            avg <- cbind(avg,avg[,2]-avg[,1])
            colnames(avg) <- c(paste(Compare[[i]],"AVG",sep="_"),"deltaBeta")
            DMPs[[i]] <- data.frame(DMPs[[i]][com.idx,],avg,probe.features[com.idx,])
        } else {
            DMPs[[i]] <- data.frame(DMPs[[i]][com.idx,],probe.features[com.idx,])
        }
    }
    message("\n[ Section 3:  Match Annotation Done ]\n")

    message("[<<<<<< ChAMP.DMP END >>>>>>]")
    message("[===========================]")
    message("[You may want to process DMP.GUI() or champ.GSEA() next.]\n")
    return(DMPs)
}
