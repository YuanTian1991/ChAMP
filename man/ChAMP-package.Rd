\name{ChAMP-package}
\alias{ChAMP-package}
\alias{ChAMP}
\docType{package}
\title{
ChAMP-Chip Analysis Methylation Pipeline
}
\description{
A pipeline that enables pre-processing of 450k data, a selection of normalization methods and novel methods for downstream analysis including Probe Lasso DMR Hunter and Copy Number Aberration analysis.
}
\details{
\tabular{ll}{
Package: \tab ChAMP\cr
Type: \tab Package\cr
Version: \tab 1.1.7\cr
Date: \tab 2014-03-05\cr
License: \tab GPL-3\cr
}
The full analysis pipeline can be run with all defaults using champ.process()

Alternatively, it can be run in steps using all funtions separately.
}

\author{
Tiffany Morris, Lee Butcher, Andy Feber, Andrew Teschendorff, Ankur Chakravarthy, Stephen Beck\cr
\bold{Maintainer}: Tiffany Morris <tiffany.morris@ucl.ac.uk>
}
\keyword{ package }
\keyword{ methylation}
\keyword{ HumanMethylation450 }
\keyword{ 450k }
\keyword{ array }
\keyword{ Beadchip}

\examples{
\donttest{
	directory=system.file("extdata",package="ChAMPdata")
	champ.process(directory=directory)
	myLoad=champ.load()
	myNorm=champ.norm()
	champ.SVD()
	batchNorm=champ.runCombat()
	limma=champ.MVP()
	lasso=champ.lasso()
	champ.CNA
	}
}