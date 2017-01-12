\name{champ.CNA}
\alias{champ.CNA}
\title{
    Inference of Copy Number Abberrations from intensity values.
}
\description{
    This function enables CNA profiles to be built using methylation data from Illumina HumanMethylation450K and HumanMethylationEPIC BeadChips. This function provide options to find Copy Number Abberrations between two phenotype (.e.g. Cancer & Normal), or the function would take the average value of your dataset as control and detect if some value are out of average status. For user want to detect abberrations between phenotypes, they can specify controlGroup in parameter, or they can simply used packaged dataset as control. Two kinds of plot would be returned, the abberrations of each sample, and the abberrations of each phenotype. The older version of ChAMP provide batchcorrect for intensity dataset, but it's nolonger provided here, user may use champ.runCombat() function to correct batch effect just like they correct beta matrix.
}
\usage{
    champ.CNA(intensity=myLoad$intensity,
              pheno=myLoad$pd$Sample_Group,
              control=TRUE,
              controlGroup="champCtls",
              sampleCNA=TRUE,
              groupFreqPlots=TRUE,
              Rplot=FALSE,
              PDFplot=TRUE,
              freqThreshold=0.3,
              resultsDir="./CHAMP_CNA",
              arraytype="450K")
}

\arguments{
	\item{intensity}{A matrix of intensity values for each sample. (default = myLoad$intensity)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
	\item{control}{If champ.CNA() should calculate the difference between groups(controls and case) of not(with average). (default = TRUE)} 
	\item{controlGroup}{which phenotype in your pheno parameter shall be treated as control type is you want to comparision between two groups. If this value was missing or invalid, the function would automatically use packaged Blood sample(champCtls) as control. (default = "champCtls")}
	\item{sampleCNA}{If sampleCNA=TRUE, then each sample's Copy Number Abberrations would be calculated and plotted. (default = TRUE)}	
	\item{groupFreqPlots}{If groupFreqPlots=TRUE, then each group's Copy Number  Abberrations Frequence would be calculated and plotted. (default = TRUE)}
	\item{freqThreshold}{If groupFreqPlots=T, then freqThreshold will be used as the cutoff for calling a gain or loss. (default = 0.3)}
    \item{PDFplot}{If PDFplot would be generated and save in resultsDir. (default = TRUE)}
    \item{Rplot}{If Rplot would be generated and save in resultsDir. Note if you are doing analysis on a server remotely, please make sure the server could connect your local graph applications. (For example X11 for linux.) (default = TRUE)}
    \item{arraytype}{Choose microarray type is 450K or EPIC.}
    \item{resultsDir}{The directory where PDF files would be saved. (default = "./CHAMP_CNA/")}
}

\value{
    \item{sampleResult}{The Copy Number Abberrations result calculated and ploted for each Sample.}
    \item{groupResult}{The Copy Number Abberrations result calculated and ploted for each Group.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myCNA <- champ.CNA()
    }
}
\references{
Feber, A et. al. (2014). CNA profiling using high density DNA methylation arrays. Genome Biology.
}
\author{
Feber, A\cr
adapted by Yuan Tian
}

\keyword{ copynumber }
