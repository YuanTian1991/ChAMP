\name{champ.GSEA}
\alias{champ.GSEA}
\title{
    Do GSEA for DMP, DMR and other methylation data related results.
}
\description{
    This function would do GSEA on the results of champ functions like DMP and DMR. However users may also add individual CpGs and genes in it. There are three method are incoporated into champ.GSEA function here. One is old Fisher Exact Test method, which will used information downloaded from MSigDB and do fisher exact test to calculated the enrichment status for each pathways. And another method is "gometh" method, which will use missMethyl package to correct the inequality between number of genes and number of CpGs, then do GSEA. The third and newest method is Empirical Bayes (ebayes) method, which does not need DMP or DMR information, but would directly calculate global test across all CpGs then do GSEA. User may assign parameter "method" as "ebayes", "gometh" or "fisher" to choose which method they want to use.
}

\usage{
    champ.GSEA(beta=myNorm,
               DMP=myDMP[[1]],
               DMR=myDMR,
               CpGlist=NULL,
               Genelist=NULL,
               pheno=myLoad$pd$Sample_Group,
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
    \item{pheno}{If use ebayes method, user needs to provide phenotype information to conduct global test. (default = myLoad$pd$Sample_Group)}
    \item{method}{Which method would be used to do GSEA?"gometh","fisher", or"ebayes". "ebayes" is our new unbias GSEA method, you could refer to champ.ebGSEA() function to know more. (default = "fisher")}
    \item{arraytype}{Which kind of array your data set is? (default = "450K")}
    \item{Rplot}{If gometh method was chosen, should Probability Weight plot will be plotted. More information please check gometh package. (default = TRUE)}
    \item{adjPval}{Adjusted p value cutoff for all calculated GSEA result. (default = 0.05)}
}
\value{
    For fisher Method:
    \item{Genelist}{List of pathway we get by enriching genes onto annotation database.}
    \item{nOVLAP}{Number of genes overlapped in your significant gene list and annotated pathways.}
    \item{OR}{Odds Ratio calculated for each enrichment.}
    \item{P-value}{Significance calculated from fisher exact test.}
    \item{adjPval}{Adjusted P value from "BH" method.}
    \item{Genes}{Name of genes enriched in each pathway.}

    For gometh method, the returned value are:
    \item{category}{GO pathway's index.}
    \item{over_represented_pvalue}{The p value for genes' over representing in this pathway.}
    \item{under_represented_pvalue}{The p value for genes' under representing in this pathway.(Not likely to be used)}
    \item{numDEInCat}{Numbers of Different Methylation Genes in this pathway.}
    \item{numInCat}{Numbers of all genes related to this pathway.}
    \item{term}{The short explaination for this pathway.}
    \item{ontology over_represented_adjPvalue}{The ajusted over representing p value with "BH" method. User may used this one to select qualitied Pathways.}

    For ebayes method:
    \item{nREP}{Number of genes enriched in this pathway.}
    \item{AUC}{Area under curve from wilcox test.}
    \item{P(WT)}{P value detected for each pathway from Wilxos Test.}
    \item{P(KPMT)}{P value from Known Population Median Test}
    \item{adjP}{Adjusted P value for each pathway, using BH method.}
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
Yuan Tian, Danyue Dong
}
\keyword{ GSEA }
