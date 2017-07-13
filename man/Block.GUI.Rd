\name{Block.GUI}
\alias{Block.GUI}
\title{
    Generate interactive plot for the result of champ.Block() function,
}
\description{
    A Shiny, Plotly and Web Brower based analysis interface. Block.GUI() is aimed to provide a comprehensive interactive analysis platform for the result of champ.Block(). The left panel indicate parameters user may be used to select significant Block, here I only provided minium number of clusters and p value as two threshold cutoff. After opening this web page, user may select their cutoff, then press submit, the webpage would calculate the result automatically.
    User could check the Blocktable in first tab easily, users can rank and select certain genes in the table, the content of the table might be changed based on the cutoff you selected in left panel.
    The second tab provide the mapping information from CpGs to Blocks, which will makes your easier to find connection between CpGs to clusters then Blocks. 
    The third tab is the plot of Block and the clusters' differential methylation information, you may search the Block you want to check by left panel, note that if there is only one significant cluster in the Block you selected, the plot might not be show properly.
}

\usage{
    Block.GUI(Block=myBlock,
              beta=myNorm,
              pheno=myLoad$pd$Sample_Group,
              runDMP=TRUE,
              compare.group=NULL,
              arraytype="450K")
}

\arguments{
    \item{Block}{The result from champ.Block(). (default = myBlock)}
    \item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{runDMP}{If DMP result sould be calculated and combined into the result of CpGs annotation.}
    \item{compare.group}{compare.group is a parameter to assign which two phenotypes you wish to analysis, if your pheno contains only 2 phenotyes you can leave it as NULL, but if your pheno contains multiple phenotypes, you MUST specify compare.group. (default = NULL)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}

\value{
    Totally three tabs would be generated on opened webpage.
    \item{Blocktable}{The Block list of all significant Blocks selected by cutoff in left panel.}
    \item{CpGtable}{Information of all significant CpGs selected by cutoff in left panel. More importantly, it also contains mapping information each between CpG ID, Cluster ID and Block ID.}
    \item{BlockPlot}{Dots and lines of all clusters involved in one Block, the xaix is based on real Map information of clusters. Above the plot, is the differential methylation information of clusters contained in this Block.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myBlock <- champ.Block()
        Block.GUI()
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
