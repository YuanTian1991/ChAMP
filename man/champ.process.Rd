\name{champ.process}
\alias{champ.process}

\title{
Process function to run all methods in ChAMP pipeline.
}
\description{
This function allows the user to run the entire pipeline in one function. Arguments allow user to select functions if desired. Note that it maybe run during champ.process() if there is any problem during the process, thus run ChAMP functions one by one is actually recommended.
}

\usage{
    champ.process(runload=TRUE,
                  directory = getwd(),
                  filters=c("XY","DetP","Beads","NoCG","SNP","MultiHit"),
                  #---champ.impute parameters below---#
                  runimpute=TRUE,
                  imputemethod="Combine",
                  #---champ.QC parameters below---#
                  runQC=TRUE,
                  QCplots=c("mdsPlot","densityPlot","dendrogram"),
                  #---champ.norm parameters below---#
                  runnorm=TRUE,
                  normalizationmethod="BMIQ",
                  #---champ.SVD parameters below---#
                  runSVD=TRUE,
                  RGEffect=FALSE,
                  #---champ.runCombat parameters below---#
                  runCombat=TRUE,
                  batchname=c("Slide"),
                  #---champ.DMP parameters below---#
                  runDMP=TRUE,
                  #---champ.DMR parameters below---#
                  runDMR=TRUE,
                  DMRmethod="Bumphunter",
                  #---champ.Block parameters below---#
                  runBlock=TRUE,
                  #---champ.GSEA parameters below---#
                  runGSEA=TRUE,
                  #---champ.EpiMod parameters below---#
                  runEpiMod=TRUE,
                  #---champ.CNA parameters below---#
                  runCNA=TRUE,
                  control=TRUE,
                  controlGroup="champCtls",
                  #---champ.refbase parameters below---#
                  runRefBase=FALSE,
                  #---universal settings---#
                  compare.group=NULL,
                  adjPVal=0.05,
                  resultsDir="./CHAMP_RESULT/",
                  arraytype="450K",
                  PDFplot=TRUE,
                  Rplot=TRUE,
                  cores=3,
                  saveStepresults=TRUE)

}

\arguments{
	\item{runload}{If champ.load() should be run? (default = TRUE)}
    \item{directory}{The folder directory of .idat files. (default = getwd())}
    \item{filters}{A character vector indicates filters should be done if load data from .idat files. You can remove some of the filters in it if you don't need that much. (default = c("XY","DetP","Beads","NoCG","SNP","MultiHit"))}
    \item{runimpute}{If champ.impute() should be run? Note that if your data contains too many NA, champ.impute() may remove not only CpGs, but also samples. (default = TRUE)}
    \item{imputemethod}{Which imputation method should be applied into champ.impute().}
    \item{runQC}{If champ.QC() should be run? (default = TRUE)}
    \item{QCplots}{A character vector indicates plots should be drawn by champ.QC(). You can remove some plots in it if you don't need them. (default = c("mdsPlot","densityPlot","dendrogram"))}
    \item{runnorm}{If champ.norm() should be run? (default = TRUE)}
    \item{normalizationmethod}{Which normalization method should be selected by champ.norm().}
    \item{runSVD}{If champ.SVD() should be run? (default = TRUE)}
    \item{RGEffect}{If Red Gree color Effect should be calculated in champ.SVD(). (default = FALSE)}
    \item{runCombat}{If champ.runCombat() should be run? (default = TRUE)}
    \item{batchname}{A character vector indicates what factors should be corrected by champ.runCombat(). (default = c("Slide"))}
    \item{runDMP}{If champ.DMP() should be run? (default = TRUE)}
    \item{runDMR}{If champ.DMR() should be run? (default = TRUE)}
    \item{DMRmethod}{Which DMR method should be applied by champ.DMR()? (default = TRUE)}
    \item{runBlock}{If champ.Block() should be run? (default = TRUE)}
    \item{runGSEA}{If champ.GSEA() should be run? (default = TRUE)}
    \item{runEpiMod}{If champ.EpiMod() should be run? (default = TRUE)}
    \item{runCNA}{If champ.CNA() should be run? (default = TRUE)}
    \item{control}{If champ.CNA() should be calculate copy number variance between case and control? (The other option for champ.CNA() is calculate copy number variance for each sample to the averaged value). (default = TRUE)}
    \item{controlGroup}{Which pheno should be treated as control group while running champ.CNA().(default = "champCtls")}
    \item{runRefBase}{If champ.refbase() should be run? (default = TRUE)}
    \item{compare.group}{Which two phenos should be compared in champ.DMP()?}
    \item{adjPVal}{The adjusted p value for each function's significant cutoff.}
    \item{resultsDir}{The directory where result should be stored. (default = "./CHAMP_RESULT/")}
    \item{arraytype}{If the data set under analysis is "450K" or "EPIC"? (default = "450K")}
    \item{PDFplot}{If PDF files should be plotted during running? (default = TRUE)}
    \item{Rplot}{If R plots should be plotted during running? (default = TRUE)}
    \item{cores}{How many cores should be used for parallel running during champ.process()? (default = 3)}
    \item{saveStepresults}{If result of each steps should be saved as .rd file into resultsDir folder? (default = TRUE)}
}

\value{
    \item{CHAMP_RESULT}{A list contains all results from each champ.method.}
}

\examples{
    \dontrun{
        directory=system.file("extdata",package="ChAMPdata")
        champ.process(directory=directory)
        }
}

\author{
    Yuan Tian
}
