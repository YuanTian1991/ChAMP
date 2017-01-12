if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad"))

champ.runCombat <- function(beta=myNorm,
                            pd=myLoad$pd,
                            batchname=c("Slide"),
                            logitTrans=TRUE)
{


    message("[===========================]")
    message("[<< CHAMP.RUNCOMBAT START >>]")
    message("-----------------------------")

    message("<< Preparing files for ComBat >>")
    
   if(length(which(is.na(beta)))>0) message(length(which(is.na(beta)))," NA are detected in your beta Data Set, which may cause fail or uncorrect of runCombat analysis. You may want to impute NA with champ.impute() function first.")

    message("[Combat correction will be proceed with ",dim(beta)[1], " probes and ",dim(beta)[2], " samples.]\n")

	################ Customise Phenotype Data ########################

    if(is.null(pd) | class(pd)=="list") stop("pd parameter in Data Frame or Matrix is necessary And must contain at least tow factors. If your pd is a list, please change its Format.")
    if(class(pd)=="matrix") pd <- as.data.frame(pd)

    valid.idx <- which(!colnames(pd) %in% c("Sample_Name","Sample_Group","Project","filenames","Basename") &
                       apply(pd,2,function(x) length(unique(x)))!=1 &
                       apply(pd,2,function(x) all(table(x)>=2)))
    if(length(valid.idx)==0) stop("There is not valid factor can be corrected. Factor can be corrected must contian at least two phenotypes, each of them must contain at least two samples. Please check if your covariates fulfill these requirement.")
    PhenoTypes.lv_tmp <- as.data.frame(pd[,valid.idx])
    colnames(PhenoTypes.lv_tmp) <- colnames(pd)[valid.idx]
    PhenoTypes.lv <- as.data.frame(apply(PhenoTypes.lv_tmp,2,function(x) if(class(x)!="numeric") as.factor(as.numeric(as.factor(x)))))
    if(!is.null(rownames(pd))) rownames(PhenoTypes.lv) <- rownames(pd)

    if(ncol(PhenoTypes.lv)>=1)
    {
        message("<< Following Factors in your pd(sample_sheet.csv) could be applied to Combat: >>")
        sapply(colnames(PhenoTypes.lv_tmp),function(x) message("<",x,">(",class(PhenoTypes.lv[[x]]),"):",paste(unique(PhenoTypes.lv_tmp[,x]),collapse=", ")))
        message("[champ.runCombat have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv).]")
    }else
    {
        stop("You don't have even one factor with at least two value to be analysis. Maybe your factors contains only one value, no variation at all...")
    }

    if(ncol(pd) > ncol(PhenoTypes.lv))
    {
        message("\n<< Following Factors in your pd(sample_sheet.csv) can not be corrected: >>")
        sapply(setdiff(colnames(pd),colnames(PhenoTypes.lv)),function(x) message("<",x,">"))
        message("[Factors are ignored because they only indicate Name or Project, or they contain ONLY ONE value across all Samples, or some phenotype contains less than 2 Samples.]")
    }

    if(all(batchname %in% colnames(PhenoTypes.lv)))
    {
        message("As your assigned in batchname parameter: ",paste(batchname,collapse=",")," will be corrected by Combat function.")
    }else
    {
        stop(setdiff(batchname,colnames(PhenoTypes.lv))," factors is not valid to run Combat, please recheck your dataset.")
    }
 	
    if(min(beta)==0)
    {
        message("Zeros in your dataset have been replaced with 0.000001")
        beta[beta==0]<-0.000001
    }

    beta_2 <- beta 
    innercombat <- function(beta,batch,formdf)
    {
        if(logitTrans) beta <- logit2(beta)
        print(formdf)
        mod <- model.matrix(formdf,data=pd)
        message("Generate mod success. Started to run ComBat, which is quite slow...")
        combat <- ComBat(dat=beta,batch=batch,mod=mod,par.prior=TRUE)
        if(logitTrans)combat=ilogit2(combat)
        return(combat)
    }
    for(i in 1:length(batchname))
    {
        message("\n<< Start Correcting ",batchname[i]," >>")
        if(i+1 <= length(batchname))
            formdf <- as.formula(paste(" ~ ",paste(c("Sample_Group",batchname[(i+1):length(batchname)]),collapse=" + "),sep=""))
        else
            formdf <- as.formula(" ~ Sample_Group")
        beta <- innercombat(beta,pd[,batchname[i]],formdf)
    }

    if(is.null(beta))
    {
        message("Sorry, champ.runCombat failed. Your old dataset will be returned.")
        return(beta_2) 
    }else
    {
        message("champ.runCombat success. Corrected dataset will be returned.")
        return(beta)
    }

    message("[<<< CHAMP.RUNCOMBAT END >>>]")
    message("[===========================]")
    message("[You may want to process champ.SVD() next.]\n")
}
