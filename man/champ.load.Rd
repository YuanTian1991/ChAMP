\name{champ.load}
\alias{champ.load}

\title{
    Upload of raw HumanMethylation450K or HumanMethylationEPIC data from IDAT files.
}
\description{
    Function that loads data from IDAT files to calculate intensity. Some kinds of filtering will be conducted as well such as unqualied CpGs, SNP, multihit sites, and XY chromosomes related CpGs.
}
\usage{
champ.load(directory = getwd(),
           methValue="B",
           filterDetP=TRUE,
           detPcut=0.01,
           removeDetP = 0,
           filterBeads=TRUE,
           beadCutoff=0.05,
           filterNoCG=TRUE,
           filterSNPs=TRUE,
           filterMultiHit=TRUE,
           filterXY=TRUE,
           arraytype="450K")
}

\arguments{
	
	\item{directory}{Location of IDAT files, default is current working directory.(default = getwd())}
	\item{methValue}{Indicates whether you prefer m-values M or beta-values B. (default = "B")}
	\item{filterDetP}{If filter = TRUE, then probes above the detPcut will be filtered out.(default = TRUE)}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. (default = 0.01)}
  	\item{removeDetP}{The removeDetP parameter represents the fraction of samples that can contain a detection p-value above the detPcut.(default = 0)}
  	\item{filterBeads}{If filterBeads=TRUE, probes with a beadcount less than 3 will be removed depending on the beadCutoff value.(default = TRUE)}
  	\item{beadCutoff}{The beadCutoff represents the fraction of samples that must have a beadcount less than 3 before the probe is removed.(default = 0.05)}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed.(default = TRUE)}
  	\item{filterSNPs}{If filterSNPs=TRUE, probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed.(default = TRUE)}
  	\item{filterMultiHit}{If filterMultiHit=TRUE, probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed.(default = TRUE)}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed.(default = TRUE)}
	\item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
}

\value{
	\item{mset}{mset object from minfi package, with filtering CpGs discarded.}
	\item{rgSet}{rgset object from minfi package function read.metharray.exp(), contains all information of a .idat methylation dataset. If you want to do more analysis than functions provided by ChAMP, you can take this as a start point.}
	\item{pd}{pd file of all sample information from Sample Sheet, which would be very frequently by following functions as DEFAULT input, thus it's not very necessarily, please don't modify it.}
	\item{intensity}{A matrix of intensity values for all probes and all samples, the information would be used in champ.CNA() function. CpGs has been filtered as well.}
	\item{beta}{A matrix of methylation scores (M or beta values) for all probes and all samples.}
	\item{detP}{A matrix of detection p-values for all probes and all samples.}
}
\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
    }
}
\author{
Yuan Tian
}
