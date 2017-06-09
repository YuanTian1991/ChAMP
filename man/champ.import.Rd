\name{champ.import}
\alias{champ.import}

\title{
    Read data from IDAT file.
}
\description{
    Function provided by ChAMP to extract value from IDAT file, and mapping between CpGs and Probes on Chip. The older version of ChAMP used minfi to load data, this is a version provided by ChAMP. The function would read data from one directory, which contains IDAT files and phenotype data csv files. Then champ.import() would firstly read the csv file, mapping between each sample and IDAT file. Then champ.import() would read IDAT file for each sample. After reading Green and Red Channal, Meth Matrix, UnMeth Matrixn beta value, intensity, detect P value, bead count would be calculated. Above are matrix would be used in champ.filter(). Note that, champ.import() would NOT do batch correction. And data read by champ.import() can not be used for SWAN normalization and FunctionNormalization in champ.norm() function. If user want to use SWAN, you may still consider champ.load() function, but remember to set "method" parameter as "minfi", which is "ChAMP" in default.
}
\usage{
champ.import(directory = getwd(),
             offset=100,
             arraytype="450K")
}

\arguments{
	
	\item{directory}{Location of IDAT files, default is current working directory.(default = getwd())}
	\item{offset}{offset is set to make sure no inf value would be returned.(default = 100)}
	\item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
}

\value{
	\item{beta}{A matrix of beta methylation scores for all probes and all samples (No filtering has been don).}
	\item{M}{A matrix of M methylation scores for all probes and all samples (No filtering has been done).}
	\item{pd}{pd file of all sample information from Sample Sheet, which would be very frequently by following functions as DEFAULT input, thus it's not very necessarily, please don't modify it.}
	\item{intensity}{A matrix of intensity values for all probes and all samples, the information would be used in champ.CNA() function. It has not been filtered. Actually, intensity are the sum of Meth Matrix and UnMeth Matrix.}
	\item{detP}{A matrix of detection p-values for all probes and all samples.}
	\item{beadcount}{A matrix beads for each probe on each sample. Value less then 3 has been set NA.}
	\item{Meth}{Methylated Matrix for all probe and all samples.}
	\item{UnMeth}{UnMethylated Matrix for all probe and all samples.}
}
\examples{
    \dontrun{
        myimport <- champ.import(directory=system.file("extdata",package="ChAMPdata"))
    }
}

\author{
Yuan Tian
}
