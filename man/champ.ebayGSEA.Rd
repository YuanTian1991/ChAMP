\name{champ.ebayGSEA}
\alias{champ.ebayGSEA}
\title{
    Empirical Bayes GSEA method.
}
\description{
    This is a newly created method for conduct bias-free GSEA from 450K or EPIC data set. This method use global test to detect significance of genes from DNA methylation data sets directly, instead of simply select genes mapped my DMPs pr DMRs. By applying this method, users could find GSEA without bias from inequality number of CpGs of genes, and detect some marginal significant genes for GSEA process. After global test, Empirical Bayes method would use wilcox test to enrich genes to pathways. Note that you can directly use champ.GSEA() to use this method, just need to set "method" parameter in champ.GSEA() to run this method.
}

\usage{
champ.ebayGSEA(beta=myNorm, pheno=myLoad$pd$Sample_Group, adjPval=0.05, arraytype="450K")
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
    \item{pheno}{If use ebayes method, user needs to provide phenotype information to conduct global test. (default = myLoad$pd$Sample_Group)}
    \item{adjPval}{Adjusted p value cutoff for all calculated GSEA result. (default = 0.05)}
    \item{arraytype}{Which kind of array your data set is? (default = "450K")}
}
\value{
    \item{nREP}{Number of genes enriched in this pathway.}
    \item{AUC}{Area under curve from wilcox test.}
    \item{P}{P value detected for each pathway.}
    \item{adjP}{Adjusted P value for each pathway, using BH method.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myGSEA.ebayes <- champ.ebayGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450K")
        }
}

\author{
Yuan Tian, Danyue Dong
}
\keyword{ GSEA }
