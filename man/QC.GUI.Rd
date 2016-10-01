\name{QC.GUI}
\alias{QC.GUI}
\title{
    Generate interactive plot for Quality Control of a dataset
}
\description{
    A Shiny, Plotly and Web Brower based analysis interface. QC.GUI() is aimed to provide  mdsplot, densityPlot, Type-I&Type-II densityplot, dendrogram(no interactable) and heatmap for top 1000 variale CpGs. 
    In the first tab,mdsplot are plotted based on the distance calculated by top 1000 variable CpGs. For dendrogram, if there are only less than 10 samples, the distance between samples are calculated by all CpGs, if there are more than 10 samples, QC.GUI() would apply SVD doconvolution on the dataset first then extract top significant components as latent variabls and calculate distance between samples. For the heatmap, if your dataset contains less than 1000 CpGs, all CpGs would be plotted, but if your dataset contains more than 1000 CpGs, the top 1000 variable CpGs would be selected and plot.
}

\usage{
    QC.GUI(beta=myLoad$beta,
           pheno=myLoad$pd$Sample_Group,
           arraytype="450K")
}

\arguments{
    \item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
     \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}

\value{
    Totally five tabs would be generated on opened webpage.
    \item{mdsplot}{A mdsplot used to see the clustering result and similarity between sampels.}
    \item{TypeDensity}{A two-line density Plot indicate Type-I CpGs and Type-II CpGs.}
    \item{QCplot}{Beta distribution of each sample. You may use it to check samples with low qualities.}
    \item{Dendrogram}{Dendrogram of all samples. If there are only less than 10 samples, the distance between samples are calculated by all CpGs, if there are more than 10 samples, QC.GUI() would apply SVD doconvolution on the dataset first then extract top significant components as latent variabls and calculate distance between samples.}
    \item{heatmap}{Heatmap for top 1000 variale CpGs.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        QC.GUI()
        }
}

\note{
    Please make sure you are running R locally or connected with local graph software(X11) remotely.
}


\author{
    Yuan Tian
}
\keyword{ shiny }
\keyword{ plotly }
