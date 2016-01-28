\name{champ.lasso}
\alias{champ.lasso}
\title{
Probe Lasso DMR Hunter
}
\description{
A method for identifying DMRs (differentially methylated regions) using a feature based dynamic window. Also offers the option to filter
SNPs based on data from the 1000 Genomes Project.
}
\usage{
champ.lasso(fromFile = FALSE, uploadResults = FALSE, uploadFile = "limma.txt", limma, 
beta.norm = myNorm$beta, pd = myLoad$pd, filterXY = TRUE, image = TRUE, mafPol.lower = 0, 
mafPol.upper = 0.05, popPol = "eur", lassoStyle = "max", lassoRadius = 2000, 
minSigProbesLasso = 3, minDmrSep = 1000, minDmrSize = 0, adjPVal = 0.05, 
adjust.method = "BH", resultsDir = paste(getwd(), "resultsChamp", sep = "/"), 
bedFile = TRUE, DMRpval = 0.05, batchDone = FALSE, normSave)
}

\arguments{
	\item{fromFile}{if}
	\item{uploadResults}{Set uploadResults=TRUE if you haven't loaded data from .idat files and need to upload the limma file} 
	\item{uploadFile}{If uploadResults=TRUE this is the file name} 
	\item{limma}{If}
	\item{beta.norm}{A matrix of values representing the methylation scores for each sample (M or B). The default assumes you ran champ.norm and saved the output to "norm".}
	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{filterXY}{If filterXY=T, probes from the X and Y chromosomes are removed.}
	\item{image}{If image=T, images will be saved as a pdf file in the resultsDir.}
	\item{mafPol.lower}{The lower limit for the minor allele frequencies of included polymorphisms}
	\item{mafPol.upper}{The upper limit for the minor allele frequencies of included polymorphisms}
	\item{popPol}{Indicates the population on which to base the polymorphic frequency Asian (asn), American (amr), African (afr) or Northern European (eur)}
	\item{lassoStyle}{Determines whether lassoRadius is the minimum (min) or maximum (max) lasso size, default = "max"}
	\item{lassoRadius}{The lasso size, default = 2000}
	\item{minSigProbesLasso}{The minimum number of significant probes to be captured in lasso, default = 3}
	\item{minDmrSep}{The minimum seperation (bp) between neighbouring DMRs, default = 1000}
	\item{minDmrSize}{The minimum DMR size (bp), default = 0}
	\item{adjPVal}{The minimum threshold of significance for probes to be includede in DMRs, default = 0.05}
	\item{adjust.method}{The p-value adjustment method to be used for the limma analyis, default= "BH" (Bonferroni-Hochberg)}
	\item{resultsDir}{Directory where results will be saved. Default is to create a folder called "resultsChamp"in the current working directory.}
	\item{bedFile}{If bedFile=TRUE, the DMRs will be saved in bedfile format for downstream analysis. Default is TRUE.}
	\item{DMRpval}{This is the significance threshold for including DMRs in the final DMR list.}
	\item{batchDone}{Internal variable to indicate if combat batch correction was performed.}
  	\item{normSave}{Internal variable to store normalized, not-batch corrected beta values.}
}

\examples{
	data(testDataSet)
	myLoad=testDataSet
	myNorm=champ.norm(norm="NONE")
}

\value{
\item{dmrList}{A matrix of DMRs is returned containing columns for probeID, deltaBeta, adjusted p-value, chromosome, map info, 
	chromosome arm, nearest feature, SNP allele frequency on forward strand, SNP allele frequence on reverse strand, distance 
	of nearest probe, radius of lasso that captured DMR, DMR number, DMR start, DMR end, DMR size, p-value for DMR}

}
\author{
Stirling, L
}

\keyword{ DMR }
