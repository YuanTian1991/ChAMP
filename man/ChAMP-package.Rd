\name{ChAMP-package}
\alias{ChAMP-package}
\alias{ChAMP}
\docType{package}
\title{
ChAMP-Chip Analysis Methylation Pipeline
}
\description{
A pipeline that enables pre-processing of 450K or EPIC data, a selection of normalization methods and a bundle of analysis method including SVD checking, Batch effect correction, DMP, DMR, Block detection, Cell proportion detection, GSEA pathway detection, EpiMod module detection, and copy number variance detection. ChAMP provided a very comprehensive analysis pipeline for EPIC or 450K data set.
}
\details{
\tabular{ll}{
Package: \tab ChAMP\cr
Type: \tab Package\cr
Version: \tab 2.8.6\cr
Date: \tab 2017-07-19\cr
License: \tab GPL-3\cr
}
The full analysis pipeline can be run with all defaults using champ.process()
Alternatively, it can be run in steps using all functions separately.
}

\author{
Yuan Tian, Tiffany Morris, Lee Stirling, Andy Feber, Andrew Teschendorff, Ankur Chakravarthy, Stephen Beck\cr
\bold{Maintainer}: Yuan Tian  <champ450k@gmail.com>
}
\keyword{ package }
\keyword{ methylation}
\keyword{ HumanMethylation450 }
\keyword{ 450k }
\keyword{ array }
\keyword{ Beadchip }
\keyword{ DNAMethylation }
\keyword{ EPIC }
\examples{
\donttest{
	directory=system.file('extdata',package='ChAMPdata')
	champ.process(directory=directory)
    ### run champ functions separately.
    myLoad <- champ.load(directory)
    myImpute <- champ.impute()
    champ.QC()
    myNorm <- champ.norm()
    champ.SVD()
    myCombat <- champ.runCombat()
    myDMP <- champ.DMP()
    myDMR <- champ.DMR()
    myBlock <- champ.Block()
    myGSEA <- champ.GSEA()
    myEpiMod <- champ.EpiMod()
    myCNA <- champ.CNA()
    myRefbase <- champ.refbase() ### for blood sample only

    CpG.GUI()
    QC.GUI()
    DMP.GUI()
    DMR.GUI()
    Block.GUI()
	}
}
