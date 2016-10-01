\name{champ.refbase}
\alias{champ.refbase}
\title{
    Applying References-Base Methold to beta valued methylation data.
}
\description{
    Applying References-Based Methold to correct cell-proportion in a methylation dataset. Reference-based method use purified whole blood cell-type specific methylation value to correct beta value dataset. Cell Proportions for each cell-type will be detected, and lm function will be used to correct beta value for 5 largest cell types. Cell type with smallest cell proportion will not be corrected.
}
\usage{
    champ.refbase(beta=myNorm,
                  arraytype="450K")
}
\arguments{
    \item{beta}{whole blood beta methylation dataset user want to correct. (default = myNorm)}
	\item{arraytype}{There are two types of purified cell-type specific references can be chosen, "450K" and "27K". By default, 450K value will be used, but user may choose 27K as well. (default = myNorm)}
}
\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myRefbase <- champ.refbase()
    }
}

\value{
	\item{CorrectedBea}{A beta valued matrix, with all value get corrected with RefBaseEWAS method. Be aware, champ.refbase will only correct top 5 cell types with largest mean cell proportions, and leave the cell with smallest mean cell proportion. User may check CellFraction result to find out which cell types are get corrected.}
	\item{CellFraction}{Proportion for each cell type.}
}
\references{
    Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, et al. (2012) DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics 13: 86. doi: 10.1186/1471-2105-13-86. pmid:22568884
}

\author{
    Houseman EA, Yuan Tian, Andrew Teschendorff
}

\keyword{ RefbaseEWAS }
