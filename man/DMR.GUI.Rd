\name{DMR.GUI}
\alias{DMR.GUI}
\title{
    Generate interactive plot for the result of champ.DMR() function,
}
\description{
    A Shiny, Plotly and Web Brower based analysis interface. DMR.GUI() is aimed to provide a comprehensive interactive analysis platform for the result of champ.DMR(). The last panel indicate parameters user may be used to select significant DMRs, here I only provided minprobes and p value as two threshold cutoff. After opening this web page, user may select their cutoff, then press submit, the webpage would calculate the result automatically. 
    User could check the DMRtable in first tab easily, users can rank and select certain genes in the table, the content of the table might be changed based on the cutoff you selected in left panel. 
    The second tab is the CpGtable, which extract all CpGs involved in selected CpGs. Note that maybe not all CpGs are DMPs.
    The thrid tab provide the plot of the DMR, just like gene plot in DMP.GUI(). Above the plot are CpGs information involved in this DMR.
    The fourth panel provide a heatmap of all CpGs involved in significant DMRs, and a gene enrichment plot. Both plot maybe not very clear to look, but user may zoom in for these two plots. Again be careful if you have a very big dataset.
    Note that the runDMP parameters will indicate if DMR.GUI() shall calculated DMP for all CpGs, which may cause slight different in the CpG table and the gene enrichment plot.
    And though there are three ways to calculate DMR, all three results from champ.DMR() are applicatable for this function. The title would changed automatically for different result.
}

\usage{
    DMR.GUI(DMR=myDMR,
            beta=myNorm,
            pheno=myLoad$pd$Sample_Group,
            runDMP=TRUE,
            compare.group=NULL,
            arraytype="450K")
}

\arguments{
    \item{DMR}{The result from champ.DMR(), all three DMR methods' result are supported. (default = myDMR)}
    \item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{runDMP}{If DMP result sould be calculated and combined into the result of CpGs annotation.}
    \item{compare.group}{compare.group is a parameter to assign which two phenotypes you wish to analysis, if your pheno contains only 2 phenotyes you can leave it as NULL, but if your pheno contains multiple phenotypes, you MUST specify compare.group. (default = NULL)} 
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}

\value{
    Totally four tabs would be generated on opened webpage.
    \item{DMRtable}{The DMR list of all significant DMR you selected by cutoff in left panel.}
    \item{CpGtable}{A CpGs annotation (with p value and t value if runDMP=TRUE) of all CpGs related with selected DMRs in tab 1.}
    \item{DMRPlot}{Dots and lines of all significant CpGs involved in one DMR, the distance between CpGs are equal, and the feature and Cgi information are marked down the plot. Above the plot, is the CpGs list involved in this DMR.}
    \item{Summary}{CpG enrichment gene barplot, hyper CpGs and hyper CpGs may be marked if runDMP=TRUE. Below is the heatmap for all significant DMRs related CpGs. Both plots maybe not that clear but zoomable.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myDMR <- champ.DMR() # All three methods supported.
        DMR.GUI()
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
