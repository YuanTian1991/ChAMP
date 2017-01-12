\name{champ.filter}
\alias{champ.filter}

\title{
    Do filtering on beta or M matrix. So that user who have no IDAT file can also do filtering.
}
\description{
    Function provided a series of filtering on oone solo beta matrix, which contained most filtering on champ.load(). For user may have no IDAT file but only beta matrix or detected P matrix, they can use this function to do filtering, then do analysis.
}
\usage{
champ.filter(beta,
             detP=NULL,
             pd,
             filterDetP=TRUE,
             detSamplecut=0.1,
             detPcut=0.01,
             removeDetP = 0,
             filterNoCG=TRUE,
             filterSNPs=TRUE,
             population=NULL,
             filterMultiHit=TRUE,
             filterXY=TRUE,
             arraytype="450K")
}

\arguments{
	
	\item{beta}{One single beta matrix to do filtering, this matrix can be beta matrix or M value matrix.}
	\item{detP}{Detected P value matrix for corresponding beta matrix, it MUST be 100\% corresponding, which can be ignored if you don't have.(default = NULL)}
    \item{pd}{pd file related to this beta matrix, must be provided, because maybe filtering would be on pd file.}
	\item{filterDetP}{If filter = TRUE, then probes above the detPcut will be filtered out.(default = TRUE)}
    \item{detSamplecut}{The detection p value threshhold for samples. Samples with above proportion of failed p value will be removed.}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. (default = 0.01)}
  	\item{removeDetP}{The removeDetP parameter represents the fraction of samples that can contain a detection p-value above the detPcut.(default = 0)}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed.(default = TRUE)}
  	\item{filterSNPs}{If filterSNPs=TRUE, probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed.(default = TRUE)}
    \item{population}{If you want to do filtering on specifical populations you may assign this parameter as one of "AFR","EAS"... The full list of population is in http://www.internationalgenome.org/category/population/. (default = TRUE)}
  	\item{filterMultiHit}{If filterMultiHit=TRUE, probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed.(default = TRUE)}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed.(default = TRUE)}
	\item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
}

\value{
	\item{beta}{Filtered matrix of methylation scores (M or beta values) for all probes and all samples.}
	\item{pd}{pd file of all sample information from Sample Sheet, which would be very frequently by following functions as DEFAULT input, thus it's not very necessarily, please don't modify it.}
}
\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"),
                             filterDetP=F,
                             filterBeads=F,
                             filterSNPs=F,
                             filterXY=F,
                             filterNoCG=F,
                             filterMultiHit=F)
        myfilter <- champ.filter(myLoad$beta,myLoad$detP,myLoad$pd)
    }
}

\references{
    Zhou W, Laird PW and Shen H: Comprehensive characterization, annotation and innovative use of Infinium DNA Methylation BeadChip probes. Nucleic Acids Research 2016
}

\author{
Yuan Tian
}
