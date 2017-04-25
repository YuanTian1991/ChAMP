if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","impute.knn"))

champ.impute <- function(beta=myLoad$beta,
                         pd=myLoad$pd,
                         method="Combine",
                         k=5,
                         ProbeCutoff=0.2,
                         SampleCutoff=0.1)
{
    message("[===========================]")
    message("[<<< ChAMP.IMPUTE START >>>>]")
    message("-----------------------------")

    if(method=="Combine")
    {
        message("<method>:Combine. (Suitable for Large Data Set)")
        ## Eliminate probes and samples with too many NA values
        colNA <- apply(beta,2,function(x) length(which(is.na(x)))/length(x))
        colValid <- which(colNA < SampleCutoff)
        message("(1): ",ncol(beta)-length(colValid)," Samples contain ",SampleCutoff," or above NA will be removed.")

        tmp_beta <- beta[,colValid]

        rowNA <- apply(tmp_beta,1,function(x) length(which(is.na(x)))/length(x))
        rowValid <- which(rowNA < ProbeCutoff)
        message("(2): ",nrow(tmp_beta)-length(rowValid)," Probes contain ",ProbeCutoff," or above NA will be removed.")

        beta <- tmp_beta[rowValid,]
        if(!is.null(pd)) pd <- pd[colValid,]
        data.m <- impute.knn(beta)$data
        message("(3): The rest NA are imputed by KNN method which parameter k as ",k,".")
    }else if(method=="Delete")
    {
        message("<method>:Delete. (Suitable for Small Data Set)")
        colNA <- apply(beta,2,function(x) length(which(is.na(x)))/length(x))
        colValid <- which(colNA < SampleCutoff)
        message("(1): ",ncol(beta)-length(colValid)," Samples contain ",SampleCutoff," or above NA will be removed.")
        beta <- beta[,colValid]
        if(!is.null(pd)) pd <- pd[colValid,]
        data.m <- beta[-unique(which(is.na(beta),arr.ind=T)[,1]),]
        message("(2): ",nrow(beta)-nrow(data.m)," Probe are removed. Even they contain only one NA.")
    }else if(method=="KNN")
    {
        message("<method>:KNN. (Suitable for Large Data Set)")
        data.m <- impute.knn(beta,k=k)$data
        message("All NA are imputed by KNN method with parameter k as ",k,".")
    }
    message("[<<<<< ChAMP.IMPUTE END >>>>]")
    message("[===========================]")
    if(!is.null(pd)) return(list(beta=data.m,pd=pd)) else return(beta=data.m)
}
