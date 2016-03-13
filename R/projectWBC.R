projectWBC <- function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE, lessThanOne=FALSE){

    if(is.null(contrastWBC)) Xmat = coefWBC
    else Xmat = coefWBC %*% t(contrastWBC)

    nCol = dim(Xmat)[2]
    nSubj = dim(Y)[2]
    
    mixCoef = matrix(0, nSubj, nCol)
    rownames(mixCoef) = colnames(Y)
    colnames(mixCoef) = colnames(Xmat)

    if(nonnegative){
        #library(quadprog)

        if(lessThanOne){
            Amat = cbind(rep(-1,nCol), diag(nCol))
            b0vec = c(-1,rep(0,nCol))
        }
        else{
            Amat = diag(nCol)
            b0vec = rep(0,nCol)
        }
        for(i in 1:nSubj){
            obs = which(!is.na(Y[,i]))
            Dmat = t(Xmat[obs,])%*%Xmat[obs,]
            mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
        }
    }
    else{
        for(i in 1:nSubj){
            obs = which(!is.na(Y[,i]))
            Dmat = t(Xmat[obs,])%*%Xmat[obs,]
            mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
    }
    return(mixCoef)
}
