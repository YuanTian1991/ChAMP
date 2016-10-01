\name{champ.QC}
\alias{champ.QC}

\title{
    Plot quality control plot, mdsplot, densityPlot, dendrogram for a data set.
}
\description{
    champ.QC() function would plot some summary plot for a dataset, including mdsplot, densityPlot, dendrogram. You may use QC.GUI() function to see even more plot interactively, like heatmap, Type-I and Type-II probes plot. Note that the dendrogram would do it's best to modify plot size automatically, but if you have too many samples like 1000+, the speed would be slow and the plot might be hard to read.
}
\usage{
    champ.QC(beta = myLoad$beta,
             pheno=myLoad$pd$Sample_Group,
             mdsPlot=TRUE,
             densityPlot=TRUE,
             dendrogram=TRUE,
             PDFplot=TRUE,
             Rplot=TRUE,
             Feature.sel="None",
             resultsDir="./CHAMP_QCimages/")
}

\arguments{
	\item{beta}{beta matrix want to be analysed. NA value are not recommended, thus you may want to use champ.impute to impute data first. colname of each sample MUST be marked. (default = myLoad$beta)}
    \item{pheno}{one Phenotype categorical vector for your dataset. NO list or dataframe or numeric. (default = myLoad$pd$Sample_Group)}
    \item{mdsPlot}{If mdsPlot would be plotted. (default = TRUE)}
    \item{densityPlot}{If densityPlot would be plotted. (default = TRUE)}
    \item{dendrogram}{If dendrogram would be plotted. (default = TRUE)}
    \item{PDFplot}{If PDFplot would be generated and save in resultsDir. (default = TRUE)}
    \item{Rplot}{If Rplot would be generated and save in resultsDir. Note if you are doing analysis on a server remotely, please make sure the server could connect your local graph applications. (For example X11 for linux.) (default = TRUE)}
    \item{Feature.sel}{Featrue Selection method when champ.QC() calculate dendrogram. Two options are provided, "None" means no featrue selection would be done, all probes would be used to calculate distance between each sample. "SVD" method means champ.QC() would firstly do SVD deconvolution on beta dataset, then use Random Theory Matrix mathod in "isva" package to calculated numbers of latent variable, and the "distance" between samples would be calcuated by top components of SVD result (similar to PCA). (default = "None")}
    \item{resultsDir}{The directory where PDF files would be saved. (default = "./CHAMP_QCimages/")}
}

\note{
    You can try to use QC.GUI() to do similar but interactively analysis.
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        champ.QC()
    }
}
\author{
Yuan Tian
}
