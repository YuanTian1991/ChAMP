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
 "resultsChamp", sep = "/"), filterXY = TRUE, QCimages = TRUE, filterDetP = TRUE, 
 detPcut = 0.01,  removeDetP = 0, filterBeads=FALSE, beadCutoff=0.05, filterNoCG=FALSE)
}

\arguments{
	
	\item{directory}{Location of IDAT files, default is current working directory.}
	\item{methValue}{Indicates whether you prefer m-values M or beta-values B.}
	\item{resultsDir}{Directory where results will be saved.}
	\item{QCimages}{If QCimages=T, then images will be saved.}
	\item{filterDetP}{If filter = T, then probes above the detPcut will be filtered out.}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed. Default is TRUE.}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. Default is 0.01}
  	\item{removeDetP}{The removeDetP parameter represents the fraction of samples that can contain a detection p-value above the detPcut. Default is 0.}
  	\item{filterBeads}{If filterBeads=TRUE, probes with a beadcount less than 3 will be removed depending on the beadCutoff value. Default is FALSE until the next Bioconductor release.}
  	\item{beadCutoff}{The beadCutoff represents the fraction of samples that must have a beadcount less than 3 before the probe is removed. Default is 0.05 or 5\% of samples.}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed. Default is FALSE.}
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
myLoad=champ.load(directory=system.file("extdata",package="ChAMPdata"),filterBeads=TRUE)
}
\author{
Morris, T
}
