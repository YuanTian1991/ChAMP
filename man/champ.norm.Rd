\name{champ.norm}
\alias{champ.norm}
\title{
Normalization for HumanMethylation450 or HumanMethylationEPIC data
}
\description{
Option to normalize data with a selection of normalization methods. There are four functions could be selected: "PBC","BMIQ","SWAN" and "FunctionalNormalize". SWAN method call for BOTH rgSet and mset input, FunctionNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. BMIQ method is the default function, which would also return normalised density plots in PDF format in results Dir. FunctionalNormalize is provided in minfi package, which ONLY support 450K data yet. Not that BMIQ function might fail if you sample's beta value distribution is not beta distribution, which occationally happen when too many CpGs are deleted while loading .idat files with champ.load() function. Also multi-cores parallel is conductable for BMIQ function, if your server or computer is good enought with more than one cores, you may assign more cores like 10 to accelerate the process. No matter what method you selected, they all will return the same result: Normalize beta matrix with effect of Type-I and Type-II probes corrected.
}

\usage{
    champ.norm(beta=myLoad$beta,
               rgSet=myLoad$rgSet,
               mset=myLoad$mset,
               resultsDir="./CHAMP_Normalization/",
               method="BMIQ",
               plotBMIQ=FALSE,
               arraytype="450K",
               cores=3)
}

\arguments{
	\item{beta}{Original beta matrix waiting to be normalized. NA value are not recommended, thus you may want to use champ.impute to impute data first. colname of each sample MUST be marked. (default = myLoad$beta)}
	\item{rgSet}{Original full information matrix from champ.load(), which is required by "SWAN" and "FunctionNormalization" method. (default = myLoad$rgSet)}
	\item{mset}{mset object from minfi package, with filtering CpGs discarded, which is required by "SWAN" method. (default = myLoad$mset)}
    \item{resultsDir}{The folder where champ.norm()'s PDF file should be saved. (default = "./CHAMP_Normalization/")}
    \item{method}{Method to do normalization: "PBC","BMIQ","SWAN" and "FunctionalNormalize". (default = "BMIQ")}
    \item{plotBMIQ}{If "BMIQ" method is choosen, should champ.norm() plot normalized plot in PDF and save it in resultsDir. (default = FALSE)}
    \item{arraytype}{Choose microarray type is "450K" or "EPIC".(default = "450K")}
    \item{cores}{If "BMIQ" method is choosen, how many cores shall be used to run parallel. (default = 3)}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        }
}

\value{
	\item{beta.p}{A matrix of normalised methylation scores (M or beta values) for all probes and all samples.}
}
\references{
    Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J, Gomez-Cabrero D, Beck S. A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450k DNA methylation data. Bioinformatics. 2013 Jan 15;29(2):189-96.

    Dedeurwaerder S, Defrance M, Calonne E, Denis H, Sotiriou C, Fuks F.Evaluation of the Infinium Methylation 450K technology. Epigenomics. 2011,Dec;3(6):771-84.

    Touleimat N, Tost J. Complete pipeline for Infinium Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation. Epigenomics. 2012 Jun;4(3):325-41.

    Fortin J. P. et al. Functional normalization of 450k methylation array data improves replication in large cancer studies. Genome Biol. 15, 503 (2014). 
}

\author{
Yuan Tian wrote the wrappers
}

\keyword{ normalization }
