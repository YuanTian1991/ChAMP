if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","probeInfoALL.lv","%dopar%","makeCluster","foreach"))

champ.norm <- function(beta=myLoad$beta,
                       rgSet=myLoad$rgSet,
                       mset=myLoad$mset,
                       resultsDir="./CHAMP_Normalization/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=3)
{
    message("[===========================]")
    message("[>>>>> ChAMP.NORM START <<<<<<]")
    message("-----------------------------")

    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.norm Results will be saved in ",resultsDir)


    message("[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]\n")

    if(method=="SWAN")
    {
        beta.p = getBeta(preprocessSWAN(rgSet,mset),"Illumina")
    }else if(method == "BMIQ")
    {
        message("<< Normalizing data with BMIQ Method >>")
        message("Note that,BMIQ function may fail for bad quality samples (Samples did not even show beta distribution).")
        if(arraytype=="EPIC") data(probeInfoALL.epic.lv) else data(probeInfoALL.lv)
        design.v <- as.numeric(lapply(probeInfoALL.lv,function(x) x)$Design[match(rownames(beta),probeInfoALL.lv$probeID)])
        if(min(beta,na.rm=TRUE)==0)
        {
            beta[beta==0] <- 0.000001
            message("Zeros in your dataset have been replaced with 0.000001\n")
        }
        current_cwd <- getwd()
        if(plotBMIQ)
            setwd(resultsDir)

        if(cores>1)
        {
            if(cores > detectCores()) cores <- detectCores()
            message(cores," cores will be used to do parallel BMIQ computing.")
            registerDoParallel(makeCluster(cores))
            beta.p <- foreach(x = 1:ncol(beta),.combine = cbind,.export=c("champ.BMIQ","blc")) %dopar% champ.BMIQ(beta[,x],design.v,sampleID=colnames(beta)[x],plots=plotBMIQ)$nbeta
        }else
        {
            beta.p <- sapply(1:ncol(beta),function(x) champ.BMIQ(beta[,x],design.v,sampleID=colnames(beta)[x],plots=plotBMIQ)$nbeta)
        }
        setwd(current_cwd)
    }else if(method == "PBC")
	{
        message("<< Normalizing data with PBC Method >>")
        if(arraytype=="EPIC") data(probeInfoALL.epic.lv) else data(probeInfoALL.lv)
        design.v <- as.numeric(lapply(probeInfoALL.lv,function(x) x)$Design[match(rownames(beta),probeInfoALL.lv$probeID)])
        if(min(beta,na.rm=TRUE)==0)
        {
            beta[beta==0] <- 0.000001
            message("Zeros in your dataset have been replaced with 0.000001\n")
        }
        beta.p=DoPBC(beta,design.v)
	}else if(method=="FunctionalNormalization")
    {
        if(is.null(rgSet)) stop("rgSet not found, it is required for FunctionalNormalization")
        if(arraytype=="EPIC") rgSet@annotation[2] <- "ilm10b3.hg19"
        beta.p <- getBeta(preprocessFunnorm(rgSet))[rownames(beta),]
    }else
    {
        stop("Please Select Normalization Method from: BMIQ,PBC, FunctionalNormalization and SWAN.")
    }
    rownames(beta.p) <- rownames(beta)
    colnames(beta.p) <- colnames(beta)
    message("[>>>>> ChAMP.NORM END <<<<<<]")
    message("[===========================]")
    message("[You may want to process champ.SVD() next.]\n")
	return(beta.p)
}
