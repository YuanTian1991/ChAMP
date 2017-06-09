\name{champ.filter}
\alias{champ.filter}

\title{
    Do filtering on beta, M, intensity, Meth and UnMeth matrix.
}
\description{
    Do filtering on beta, M, intensity, Meth and UnMeth matrix. So that user who have no IDAT file can also do filtering. This function has been totally recoded, firstly it is not take result from champ.import() as input and do filtering on that. So that user can use champ.import() + champ.filter() to generate data set. The other purpose of doing that is user can do any filtering on any of above 5 matrix, as long as they have a single matrix. Not that some accessory data sets are required for some methods, like you MUST provide detect P value matrix so that champ.filtering can do filtering on detect P value. Also you have to provide beadcount information so that champ.filtering can do filtering on beadcount. Also remember that, if you want to keep pd file in accord with your data matrix, you can surely input pd, but make sure Sample_Name of your pd file is EXACTLY the same as your data matrix's colnames. Also keep in mind that, if you want to do filtering on multiple data matrix, you MUST make sure they have EXACTLY the same rownames and colnames. The function would filtering all matrix at the same time, so keeping the two names same make sure champ.filter() is not doing wrong filtering on different data sets.
}
\usage{
champ.filter(beta=myImport$beta,
             M=NULL,
             pd=myImport$pd,
             intensity=NULL,
             Meth=NULL,
             UnMeth=NULL,
             detP=NULL,
             beadcount=NULL,
             autoimpute=TRUE,
             filterDetP=TRUE,
             ProbeCutoff=0,
             SampleCutoff=0.1,
             detPcut=0.01,
             filterBeads=TRUE,
             beadCutoff=0.05,
             filterNoCG = TRUE,
             filterSNPs = TRUE,
             population = NULL,
             filterMultiHit = TRUE,
             filterXY = TRUE,
             fixOutlier = TRUE,
             arraytype = "450K")
}

\arguments{
	
	\item{beta}{One single beta matrix to do filtering. (default = myImport$beta).}
	\item{M}{One single M matrix to do filtering. (default = NULL).}
    \item{pd}{pd file related to this beta matrix, suggest provided, because maybe filtering would be on pd file. (default = myImport$pd)}
	\item{intensity}{intensity matrix. (default = NULL).}
	\item{Meth}{Methylated matrix. (default = NULL).}
	\item{UnMeth}{UnMethylated matrix. (default = NULL).}
	\item{detP}{Detected P value matrix for corresponding beta matrix, it MUST be 100\% corresponding, which can be ignored if you don't have.(default = NULL)}
	\item{beadcount}{Beadcount information for Green and Red Channal, need for filterBeads.(default = NULL)}
    \item{autoimpute}{If after detect P filtering, some NA are still exist in your data set (Only beta or M matrix), should imputation be done one them. Should only be done on big data set. Before do imputation, checking process would be done ahead to make sure Detect P, ProbeCutoff, beta or M valule are exist. (default = TRUE)}
	\item{filterDetP}{If filter = TRUE, then probes above the detPcut will be filtered out.(default = TRUE)}
    \item{SampleCutoff}{The detection p value threshhold for samples. Samples with above proportion of failed p value will be removed. (default = 0.1)}
    \item{ProbeCutoff}{The detection p value threshhold for Probe. After removing failed Samples(controled by SampleCutoff parameter), probes with above proportion of failed p value will be removed.(default = 0)}
  	\item{detPcut}{The detection p-value threshhold. Probes about this cutoff will be filtered out. (default = 0.01)}
    \item{filterBeads}{probes with less then 3 beads would be set NA. If for one probe, number of NAs above certian ratio, filtering would be conducter on that probe. (default = TRUE)}
    \item{beadCutoff}{Ratio threshhold that a probe should be removed for failed in beadcount check (default = 0.05).}
  	\item{filterNoCG}{If filterNoCG=TRUE, non-cg probes are removed.(default = TRUE)}
  	\item{filterSNPs}{If filterSNPs=TRUE, probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed.(default = TRUE)}
    \item{population}{If you want to do filtering on specifical populations you may assign this parameter as one of "AFR","EAS"... The full list of population is in http://www.internationalgenome.org/category/population/. (default = TRUE)}
  	\item{filterMultiHit}{If filterMultiHit=TRUE, probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed.(default = TRUE)}
	\item{filterXY}{If filterXY=TRUE, probes from X and Y chromosomes are removed.(default = TRUE)}
	\item{fixOutlier}{If fixOutlier=TRUE, in beta matrix only, value below 0 would be replaced as minium positive value, would value above 1 would be replaced as maxium value below 1.(default = TRUE)}
	\item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
}

\value{
	\item{Objects}{A list of data sets you want to filtering and inputted into this function.}
}
\examples{
    \dontrun{
        myimport <- champ.import(directory=system.file("extdata",package="ChAMPdata"))
        myfilter <- champ.filter(beta=myImport$beta,pd=myImport$pd,detP=myImport$detP,beadcount=myImport$beadcount)
    }
}

\references{
    Zhou W, Laird PW and Shen H: Comprehensive characterization, annotation and innovative use of Infinium DNA Methylation BeadChip probes. Nucleic Acids Research 2016
}

\author{
Yuan Tian
}
