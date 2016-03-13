\name{champ.reffree}
\alias{champ.reffree}
\title{
    Applying RefFreeEWAS Methold to beta valued methylation data.
}
\description{
    Applying RefFreeEWAS method to beta valued methylation data. This method does not rely on puritied cell reference, thus can be easily used on tissue data set, while RefbaseEWAS can only be used to whole blood samples. Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types. This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption. Returns mixture-adjusted Beta and unadjusted Bstar, as well as estimates of various latent quantities.
}

\usage{
champ.reffree(beta=myLoad$beta,design=myLoad$pd$Sample_Group,K=NULL,nboot=10)
}
\arguments{
    \item{beta}{Methylation beta valueed dataset user want to do RefFreeEWAS.}
	\item{design}{Design matrix (subjects x covariates). This parameter MUST be a vecter or a matrix. Though Characters are allowed because inside the function, character covariates will be transformed into numeric, we still recommend user input numeric deigned covariates matrix or vector.}
    \item{K}{Number of latent variable. If this value was ignored, function will use Random Matrix Theory from isva pacakge to estimate latent variables.}
    \item{nboot}{Number for Bootstrap on result of RefFreeEWAS.}
}
\examples{
	myLoad=champ.load(directory=system.file("extdata",package="ChAMPdata"))	
    myRefFree=champ.reffree()
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

\keyword{ celltype }
\keyword{ RefFreeEWAS }
