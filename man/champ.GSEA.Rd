\name{champ.GSEA}
\alias{champ.GSEA}
\title{
    Do GSEA for DMP, DMR and other methylation data related results.
}
\description{
    New Modification: We adaptted gometh method as default function to do GSEA while correcting bias caused by inequality number of CpGs on each genes. Previously, we used goseq package to do this. gometh function for missMethyl package also used the idea of goseq but implemented totally different from goseq function. The idea of gometh is firstly fit a line between probability that a gene would be choosed as significant differential genes, and numbers of CpGs enriched in these genes. Then the probability would be fixed while doing GSEA. For more details user may refer missMethyl package.
    This function would do GSEA on the results of champ functions like DMP and DMR. However users may also add individual CpGs and genes in it. There are two method are incoporated into champ.GSEA function here. One is old Fisher Exact Test method, which will used information downloaded from MSigDB and do fisher exact test to calculated the enrichment status for each pathways. And the other is "gometh" method, which will use missMethyl package to correct the inequality between number of genes and number of CpGs, then do GSEA. User may assign parameter "method" as "gometh" or "fisher" to choose which method they want to use.
}

\usage{
    champ.GSEA(beta=myNorm,
               DMP=myDMP[[1]],
               DMR=myDMR,
               CpGlist=NULL,
               Genelist=NULL,
               method="fisher",
               arraytype="450K",
               Rplot=TRUE,
               adjPval=0.05)
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{DMP}{Results from champ.DMP() function. (default = myDMP)}
    \item{DMR}{Results from champ.DMR() function. (default = myDMR)}
    \item{CpGlist}{Apart from previous parameters, if you have any other CpGs list want to do GSEA, you can input them here as a list. (default = NULL)}
    \item{Genelist}{Apart from previous parameters, if you have any other Gene list want to do GSEA. you can inpute them here as a list. (default = NULL)}
    \item{method}{Which method would be used to do GSEA?"gometh" or "fisher". (default = "gometh")}
    \item{arraytype}{Which kind of array your data set is? (default = "450K")}
    \item{Rplot}{If gometh method was chosen, should Probability Weight plot will be plotted. More information please check gometh package. (default = TRUE)}
    \item{adjPval}{Adjusted p value cutoff for all calculated GSEA result. (default = 0.05)}
}
\value{
    For fisher Method:
    \item{Block}{A data.frame contains all detected Blocks, with colnames as chr, start, end, value, area, cluster, indexStart, indexEnd, L, clusterL, p.value, fwer, p.valueArea, fwerArea. The result format is actually the same as Bumphunter, you may refer to Bumphunter packages to get more explaination about the result.}
    \item{clusterInfo}{When champ.Block() detection significant Blocks, a group of candidate Blocks would be detected out at first, this is the data frame of all candidate Blocks. The "TRUE" Blocks in above value are located in these candidate Blocks.} 
    \item{allCLID.v}{The first step of detectiong methylation Blocks is to get each probes into a cluster(region). This value is the clustering result of each probes.}
    \item{avbetaCL.m}{The beta matrix for each cluster. The value is calculated by taking mean value of all probes located in each cluster.}
    \item{posCL.m}{Position of each cluster, which is calculated by average all probes' position in each cluster.}
    \item{adjPval}{P value cutoff for calculated GSEA results. (default = 0.05)}

    For gometh method, the returned value are:
    \item{category}{GO pathway's index.}
    \item{over_represented_pvalue}{The p value for genes' over representing in this pathway.}
    \item{under_represented_pvalue}{The p value for genes' under representing in this pathway.(Not likely to be used)}
    \item{numDEInCat}{Numbers of Different Methylation Genes in this pathway.}
    \item{numInCat}{Numbers of all genes related to this pathway.}
    \item{term}{The short explaination for this pathway.}
    \item{ontology over_represented_adjPvalue}{The ajusted over representing p value with "BH" method. User may used this one to select qualitied Pathways.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myDMP <- champ.DMP()
        myDMR <- champ.DMR()
        myGSEA <- champ.GSEA()
        }
}

\author{
Yuan Tian
}
\keyword{ GSEA }
