\name{champ.load}
\alias{champ.load}

\title{
    Upload of raw HumanMethylation450K or HumanMethylationEPIC data from IDAT files.
}
\description{
    Function that loads data from IDAT files to calculate intensity. Some kinds of filtering will be conducted as well such as unqualied CpGs, SNP, multihit sites, and XY chromosomes related CpGs The old champ.load() function used minfis read.metharray.exp as API to read IDAT file, but now we also provide another approch that use that use "ChAMP" method to read data. The two method should return the same result. However, based on my checking, "ChAMP" method and "Minfi" method would reach exactly the same result at last on 450K array, but for EPIC array, since "ChAMP" is using latest B4 version manifest, so there are some mino difference between the result of these two method. Basically, "ChAMP" method is a recode version of minfi's function, except not using any Objects or GenomeRange class, which make ChAMP method's result much easier to read, however, minfi's API provided more functions, and two normalization function "SWAN" and "FunctionNormalization" both need Object of minfi's method. So for user want to use minfi method to read IDAT file, you can set method parameter as "minfi". For user who want to use "ChAMP" method, actually "ChAMP" method is a combination of champ.import() function and champ.filter() function. You may use champ.import() function to get more result, like detectP value, beadcount, Meth and Unmeth Matrix and so one. We will keep working and make "ChAMP" method more comprehensive in the future.
}
\usage{
champ.load(directory = getwd(),
           method = "ChAMP",
           methValue="B",
           autoimpute=TRUE,
           filterDetP=TRUE,
           ProbeCutoff=0,
           SampleCutoff=0.1,
           detPcut=0.01,
           filterBeads=TRUE,
           beadCutoff=0.05,
           filterNoCG=TRUE,
           filterSNPs=TRUE,
           population=NULL,
           filterMultiHit=TRUE,
           filterXY=TRUE,
           arraytype="450K")
}

\arguments{
	
	\item{directory}{Location of IDAT files, default is current working directory.(default = getwd())}
    \item{method}{which way you want to use to load data? "ChAMP" is a new reading method provided by ChAMP author, would "minfi" means using minfi's loading function. (default = "ChAMP")}
	\item{methValue}{Indicates whether you prefer m-values M or beta-values B. (default = "B")}
    \item{autoimpute}{If after filtering (or not do filtering) there are NA values in it, should impute.knn(k=3) should be done for the rest NA?}
	\item{filterDetP}{If filter = TRUE, then probes above the detPcut will be filtered out.(default = TRUE)}
    \item{ProbeCutoff}{The NA ratio threshhold for probes. Probes with above proportion of NA will be removed.}
    \item{SampleCutoff}{The failed p value (or NA) threshhold for samples. Samples with above proportion of failed p value (NA) will be removed.}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. (default = 0.01)}
  	\item{filterBeads}{If filterBeads=TRUE, probes with a beadcount less than 3 will be removed depending on the beadCutoff value.(default = TRUE)}
  	\item{beadCutoff}{The beadCutoff represents the fraction of samples that must have a beadcount less than 3 before the probe is removed.(default = 0.05)}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed.(default = TRUE)}
  	\item{filterSNPs}{If filterSNPs=TRUE, probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed.(default = TRUE)}
    \item{population}{If you want to do filtering on specifical populations you may assign this parameter as one of "AFR","EAS"... The full list of population is in http://www.internationalgenome.org/category/population/. (default = TRUE)}
  	\item{filterMultiHit}{If filterMultiHit=TRUE, probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed.(default = TRUE)}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed.(default = TRUE)}
	\item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
}

\value{
    For ChAMP method:
    \item{beta}{A matrix of beta methylation scores for all probes and all samples (Filtered by champ.filter()).}
    \item{M}{A matrix of M methylation scores for all probes and all samples (Filtered by champ.filter()).}
    \item{pd}{pd file of all sample information from Sample Sheet, which would be very frequently by following functions as DEFAULT input, thus it's not very necessarily, please don't modify it. (filtered bty champ.filter())}
    \item{intensity}{A matrix of intensity values for all probes and all samples, the information would be used in champ.CNA() function. It has not been filtered. Actually, intensity are the sum of Meth Matrix and UnMeth Matrix. (Filtered by champ.filter())}


    For minfi method:
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

\references{
    Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD and Irizarry RA (2014). Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays. Bioinformatics, 30(10), pp. 1363-1369. doi: 10.1093/bioinformatics/btu049.
    
    Jean-Philippe Fortin, Timothy Triche, Kasper Hansen. Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array. bioRxiv 065490; doi: https://doi.org/10.1101/065490

    Zhou W, Laird PW and Shen H: Comprehensive characterization, annotation and innovative use of Infinium DNA Methylation BeadChip probes. Nucleic Acids Research 2016
}

\author{
Yuan Tian
}
