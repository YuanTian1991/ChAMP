if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","CellTypeMeans27K","CellTypeMeans450K"))

champ.refbase <- function(beta=myLoad$beta,arraytype="450K")
{
    if(arraytype == "27K")
    {
        data(CellTypeMeans27K)     # From Accomando et al. (27K)
        DMRs27K <- intersect(rownames(CellTypeMeans27K), rownames(beta))
        cellFrac <- projectWBC(beta[DMRs27K,], CellTypeMeans27K[DMRs27K,], lessThanOne=TRUE)
    }else
    {
        data(CellTypeMeans450K)    # From Reinius et al. (450K)
        DMRs450K <- intersect(rownames(CellTypeMeans450K), rownames(beta))
        cellFrac <- projectWBC(beta[DMRs450K,],CellTypeMeans450K[DMRs450K,],lessThanOne=TRUE)
    }

    lm.o <- lm(t(beta) ~ cellFrac[,1]+cellFrac[,2]+cellFrac[,3]+cellFrac[,4]+cellFrac[,5]+cellFrac[,6])
    tmp.m <- t(lm.o$res)+rowMeans(beta);

	return(list(CorrectedBeta=tmp.m,CellFraction=cellFrac))
}
