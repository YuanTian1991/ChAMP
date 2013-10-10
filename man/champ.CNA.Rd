\name{champ.CNA}
\alias{champ.CNA}
\title{
Inference of Copy Number Abberrations from intensity values.}
\description{
This function enables CNA profiles to be built using methylation data from Illumina HumanMethylation450 BeadChips. 
}
\usage{
champ.CNA(intensity = myLoad$intensity, pd = myLoad$pd, loadFile = FALSE, batchCorrect = TRUE, 
file = "intensity.txt", resultsDir = paste(getwd(), "resultsChamp", sep = "/"), 
sampleCNA=TRUE, plotSample=TRUE, filterXY = TRUE, groupFreqPlots=TRUE,freqThreshold=0.3)
}

\arguments{
	\item{intensity}{A matrix of intensity values for each sample. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{loadFile}{If loadFile=TRUE, intensity data will be loaded from a separate file. Default is FALSE.}
	\item{batchCorrect}{If batchCorrect=TRUE ComBat will be run on the data to correct for batch effects due to sentrixID/slide number. Default is TRUE.}
	\item{file}{If loadFile=T this is the name of the file with the intensity values. Default is "intensity.txt".} 
	\item{resultsDir}{Directory where results will be saved. Default is a folder in the current working directory called "resultsChamp".}
	\item{sampleCNA}{If sampleCNA=TRUE, then . Default is TRUE.}	
	\item{plotSample}{If sampleCNA=TRUE and plotSample=TRUE, then CNA plots will be saved for each sample. Default is TRUE.}
	\item{filterXY}{Probes from X and Y chromosomes are removed. Default is TRUE.}
	\item{groupFreqPlots}{If groupFreqPlots=T, then}
	\item{freqThreshold}{If groupFreqPlots=T, then freqThreshold will be used as the cutoff for calling a gain or loss. Default is 0.03.}
}
\examples{
	data(testDataSet)
	myLoad=testDataSet
	champ.CNA(batchCorrect=FALSE)
}
\references{
Feber, A et. al. (2013). CNA profiling using high density DNA methylation arrays, under review. Genome Biology.

}
\author{
Feber, A\cr
adapted by Morris, T
}

\keyword{ copynumber }
