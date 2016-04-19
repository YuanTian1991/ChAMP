\name{champ.norm}
\alias{champ.norm}
\title{
Normalization of HumanMethylation450 data
}
\description{
Option to normalize data with a selection of normalization methods.
}

\usage{
champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset, 
sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp", 
sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile, 
filter = TRUE, filterXY = TRUE, QCimages = FALSE, plotBMIQ = FALSE,arraytype="450K")
}

\arguments{
    \item{arraytype}{Choose microarray type is 450K or EPIC.}
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{rgSet}{An rgSet object that was created when data was loaded the data from the .idat files. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{mset}{Loads an mset object that was created when data was loaded from the .idat files. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{sampleSheet}{If the data has not been loaded from .idat files and fromFile=TRUE then this points to the required sampleSheet. Default is "sampleSheet.txt".}
	\item{resultsDir}{Directory where results will be saved. Default is a folder in the current working directory called "resultsChamp".}
	\item{methValue}{Indicates whether you prefer the methylation scores to be calculated as m-values (M) or beta-values (B). Default is B.}
	\item{fromIDAT}{If fromIDAT=T,}
	\item{norm}{This specifies which normalization method will be used. Values can be BMIQ (by default), PBC, SWAN or NONE.}
  	\item{fromFile}{If loadFile=TRUE, then the beta values and sample sheet need to be  uploaded.}
	\item{betaFile}{If }
	\item{filter}{Not yet implemented. If fromFile=T and this is from a genome studio file, probes that have a detection p-value below detPcut are filtered out. Default is TRUE.}
	\item{filterXY}{If fromFile=True, probes from X and Y chromosomes are removed. Default is TRUE.}
	\item{QCimages}{If QCimages=TRUE, then quality control images are saved to the resultsDir. Default is TRUE.}
	\item{plotBMIQ}{If plotBMIQ=TRUE and norm="BMIQ", BMIQ plots will be saved. Default is TRUE.}
}
\value{
	\item{beta}{A matrix of normalised methylation scores (M or beta values) for all probes and all samples.}
}
\references{

Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J, Gomez-Cabrero D,
Beck S. A beta-mixture quantile normalization method for correcting probe design 
bias in Illumina Infinium 450k DNA methylation data. Bioinformatics. 2013 Jan
15;29(2):189-96.

Dedeurwaerder S, Defrance M, Calonne E, Denis H, Sotiriou C, Fuks F.Evaluation of the Infinium Methylation 450K technology. Epigenomics. 2011,Dec;3(6):771-84.

Touleimat N, Tost J. Complete pipeline for Infinium Human Methylation 450K 
BeadChip data processing using subset quantile normalization for accurate DNA
methylation estimation. Epigenomics. 2012 Jun;4(3):325-41.
}

\author{
Morris, T. wrote the wrappers
}

\keyword{ normalization }
