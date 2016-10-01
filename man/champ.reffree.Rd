\name{champ.reffree}
\alias{champ.reffree}
\title{
    Applying RefFreeEWAS Methold to beta valued methylation data.
}
\description{
    Applying RefFreeEWAS method to beta valued methylation data. This method does not rely on puritied cell reference, thus can be easily used on tissue data set, while RefbaseEWAS can only be used to whole blood samples. Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types. This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption. Returns mixture-adjusted Beta and unadjusted Bstar, as well as estimates of various latent quantities.
}

\usage{
    champ.reffree(beta=myNorm,
                  pheno=myLoad$pd$Sample_Group,
                  K=NULL,
                  nboot=50)
}
\arguments{
    \item{beta}{Methylation beta valueed dataset user want to do RefFreeEWAS. (default = myNorm)}
	\item{pheno}{Phenotype user want to find differential methylation result. This parameter MUST be a vecter or a matrix. Though Characters are allowed because inside the function, character covariates will be transformed into numeric, we still recommend user input numeric deigned covariates matrix or vector.(default = myLoad$pd$Sample_Group)}
    \item{K}{Number of latent variable. If this value was ignored, function will use Random Matrix Theory from isva pacakge to estimate latent variables. (default = NULL)}
    \item{nboot}{Number for Bootstrap on result of RefFreeEWAS. (default=50)}
}
\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))	
        myNorm <- champ.norm()
        myRefFree <- champ.reffree()
    }
}

\value{
	\item{RefFreeEWASModel}{RefFreeEWASModel S4 Object from RefFreeEWAS pacakge, contains adjusted beta value and unadjusted beta value (Bstar).}
	\item{pvBeta}{p value of each covariates, calculated from cell type mixture corrected Beta value.}
	\item{qvBeta}{q value of each covariates, calculated from cell type mixture corrected Beta value.}
}
\references{
    Houseman EA, Kile ML, et al., Reference-free deconvolution of DNA methylation data and mediation by cell composition effects (2016). http://biorxiv.org/content/early/2016/01/23/037671.
}

\author{
    Houseman EA, Yuan Tian, Andrew Teschendorff
}

\keyword{ RefFreeEWAS }
