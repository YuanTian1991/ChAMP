\name{champ.GSEA}
\alias{champ.GSEA}
\title{
    Do GSEA for MVP and DMR results.
}
\description{
    This function would do GSEA on the results of champ functions like MVP and DMR. However users may also add individual CpGs and genes in it.
}

\usage{
    champ.GSEA(beta=myNorm,
               MVP=myMVP,
               DMR=myDMR,
               CpGlist=NULL,
               Genelist=NULL,
               arraytype="450K",
               adjPval=0.05)
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{MVP}{Results from champ.MVP() function. (default = myMVP)}
    \item{DMR}{Results from champ.DMR() function. (default = myDMR)}
    \item{CpGlist}{Apart from previous parameters, if you have any other CpGs list want to do GSEA, you can input them here as a list. (default = NULL)}
    \item{Genelist}{Apart from previous parameters, if you have any other Gene list want to do GSEA. you can inpute them here as a list. (default = NULL)}
    \item{arraytype}{Which kind of array your data set is? (default = "450K")}
    \item{adjPval}{Adjusted p value cutoff for all calculated GSEA result. (default = 0.05)}
}
\value{
    \item{Block}{A data.frame contains all detected Blocks, with colnames as chr, start, end, value, area, cluster, indexStart, indexEnd, L, clusterL, p.value, fwer, p.valueArea, fwerArea. The result format is actually the same as Bumphunter, you may refer to Bumphunter packages to get more explaination about the result.}
    \item{clusterInfo}{When champ.Block() detection significant Blocks, a group of candidate Blocks would be detected out at first, this is the data frame of all candidate Blocks. The "TRUE" Blocks in above value are located in these candidate Blocks.} 
    \item{allCLID.v}{The first step of detectiong methylation Blocks is to get each probes into a cluster(region). This value is the clustering result of each probes.}
    \item{avbetaCL.m}{The beta matrix for each cluster. The value is calculated by taking mean value of all probes located in each cluster.}
    \item{posCL.m}{Position of each cluster, which is calculated by average all probes' position in each cluster.}
    \item{adjPval}{P value cutoff for calculated GSEA results. (default = 0.05)}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myMVP <- champ.MVP()
        myDMR <- champ.DMR()
        myGSEA <- champ.GSEA()
        }
}

\author{
Yuan Tian
}
\keyword{ GSEA }
