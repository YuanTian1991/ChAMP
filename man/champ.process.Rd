\name{champ.process}
\alias{champ.process}

\title{
Process function to run all methods in ChAMP pipeline.
}
\description{
This function allows the user to run the entire pipeline in one function. Arguments allow user to select functions if desired.
}

\usage{
champ.process(fromIDAT = TRUE, fromFile = FALSE, directory = getwd(), resultsDir = 
paste(getwd(), "resultsChamp", sep = "/"), methValue = "B", filterDetP = TRUE, 
detPcut = 0.01, filterXY = TRUE, removeDetP = 0, filterBeads = TRUE, beadCutoff = 
0.05, filterNoCG = FALSE, QCimages = TRUE, batchCorrect = TRUE, runSVD = 
TRUE, studyInfo = FALSE, infoFactor = c(), norm = "BMIQ", adjust.method = "BH", 
adjPVal = 0.05, runDMR = TRUE,  runCNA = TRUE, plotBMIQ = FALSE, DMRpval = 0.05,
sampleCNA=TRUE,plotSample = TRUE,groupFreqPlots=TRUE,freqThreshold=0.3, bedFile 
= FALSE, methProfile = FALSE, controlProfile = FALSE,arraytype="450K")
}

\arguments{
	\item{fromIDAT}{If fromIDAT=TRUE, data is imported from .idat files with an associated sample sheet (.csv). If rawdata=FALSE then data is uploaded from a text file (saved as "beta.txt". Default is TRUE.)} 
	\item{fromFile}{The}
	\item{directory}{The directory where the .idat files and sample sheet are located, default is current working directory.}
	\item{resultsDir}{Directory where results will be saved. Default is to create a folder called "resultsChamp"in the current working directory.}
	\item{methValue}{Indicates whether you prefer the methylation scores to be calculated as m-values (M) or beta-values (B). Default is B.}
	\item{filterDetP}{If filter=TRUE, probes that have a detection p-value below detPcut are filtered out. Default is TRUE.}
	\item{detPcut}{If filter=TRUE, this value with be used as the significance threshold for filtering out probes based on the detection p-value. Default=0.01.}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed. Default is TRUE.}
	\item{QCimages}{If QCimages=TRUE, then quality control images are saved to the resultsDir. Default is TRUE.}
	\item{removeDetP}{The removeDetP parameter represents the fraction of samples that can contain a detection p-value above the detPcut. Default is 0.}
  	\item{filterBeads}{If filterBeads=TRUE, probes with a beadcount less than 3 will be removed depending on the beadCutoff value. Default is TRUE.}
  	\item{beadCutoff}{The beadCutoff represents the fraction of samples that must have a beadcount less than 3 before the probe is removed. Default is 0.05 or 5 percent of samples.}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed. Default is FALSE.}
	\item{batchCorrect}{If batchCorrect=TRUE, then the ComBat batch correction will be performed on batch effects related to bead chip. Default is TRUE.}
	\item{runSVD}{If runSVD=TRUE, SVD analysis for identifying batch effects will be performed. Default is TRUE.}
	\item{studyInfo}{If runSVD = TRUE, additional study covariate information can be included in the SVD analysis. Default is FALSE.}
	\item{infoFactor}{This}
	\item{norm}{This specifies which normalization method will be used. Values can be BMIQ (by default), PBC, SWAN or NONE.}
	\item{adjPVal}{The minimum threshold of significance for probes to be includede in DMRs, default = 0.05}
	\item{adjust.method}{The p-value adjustment method to be used for the limma analyis, default= BH (Bonferroni-Hochberg)}
	\item{runDMR}{If runDMR=TRUE, runs the probe lasso method for finding DMRs. This will result in an MVP list with p-values and a DMR list with p-values. Default is TRUE.}
	\item{runCNA}{If runCNA=TRUE, copy number abberation analysis will be performed. Default is TRUE.}
	\item{plotBMIQ}{If plotBMIQ=TRUE and norm="BMIQ", BMIQ plots will be saved. Default is TRUE.}
	\item{DMRpval}{If runDMR=TRUE, this value will be used as the cutoff for the DMR p-value. Default is 0.05.}
	\item{sampleCNA}{If sampleCNA=TRUE, then . Default is TRUE.}
	\item{plotSample}{If plotSample=TRUE, CNA plots will be saved. Default is TRUE.}
	\item{groupFreqPlots}{If groupFreqPlots=T, then}
	\item{freqThreshold}{If groupFreqPlots=T, then freqThreshold will be used as the cutoff for calling a gain or loss. Default is 0.03.}
	\item{bedFile}{if bedFile = TRUE. MVP list will be saved as an additional file in bedfile format for downstream analysis. Defaults is TRUE.}
	\item{methProfile}{If methProfile=TRUE then the beta values will be uploaded using the MethylationProbeProfile file from Genome Studio. Default is FALSE.}
	\item{controlProfile}{If rawdata = FALSE and runSVD = TRUE, then it is useful to have a control probe profile file exported from Genome Studio so that internal control probes can be included in the SVD analyis. Default is FALSE.}
    \item{arraytype}{Choose microarray type is 450K or EPIC.}
}

\examples{
\donttest{
	directory=system.file("extdata",package="ChAMPdata")
	champ.process(directory=directory)
}
}

\author{
Morris, T
}

\keyword{ 450K data }
