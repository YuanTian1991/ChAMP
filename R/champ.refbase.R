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
    message("Mean value for each estimated Cell Proportion:")
    print(colMeans(cellFrac))
    message(names(which.min(colMeans(cellFrac)))," has smallest cell proportion, all other cell proportions will be corrected by linear regression method.")

    lm.o <- lm(t(beta) ~ cellFrac[,-1*which.min(colMeans(cellFrac))])
    tmp.m <- t(lm.o$res)+rowMeans(beta);

	return(list(CorrectedBeta=tmp.m,CellFraction=cellFrac))
}
