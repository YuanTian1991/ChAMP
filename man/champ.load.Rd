\name{champ.load}
\alias{champ.load}

\title{
Upload of raw HumanMethylation450 data from IDAT files.

}
\description{
Function that loads data from IDAT files to calculate intensity and produce quality control images.
}
\usage{
champ.load(directory = getwd(), methValue = "B", resultsDir = paste(getwd(),
 "resultsChamp", sep = "/"), filterXY = TRUE, QCimages = TRUE, filter = TRUE, 
 detPcut = 0.01)
}

\arguments{
	
	\item{directory}{Location of IDAT files, default is current working directory.}
	\item{methValue}{Indicates whether you prefer m-values M or beta-values B.}
	\item{resultsDir}{Directory where results will be saved.}
	\item{QCimages}{If QCimages=T, then images will be saved.}
	\item{filter}{If filter = T, then probes above the detPcut will be filtered out.}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed. Default is TRUE.}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. Default is 0.01}
}

\value{
	\item{mset}{mset object}
	\item{rgSet}{rgset object}
	\item{pd}{pd file of all sample information from Sample Sheet}
	\item{intensity}{A matrix of intensity values for all probes and all samples.}
	\item{beta}{A matrix of methylation scores (M or beta values) for all probes and all samples.}
	\item{detP}{A matrix of detection p-values for all probes and all samples.}
}
\examples{
myLoad=champ.load(directory=system.file("extdata",package="ChAMPdata"))
}
\author{
Morris, T
}
