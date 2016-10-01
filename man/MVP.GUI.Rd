\name{MVP.GUI}
\alias{MVP.GUI}
\title{
    Generate interactive plot for the result of champ.MVP() function,
}
\description{
    A Shiny, Plotly and Web Brower based analysis interface. MVP.GUI() is aimed to provide a comprehensive interactive analysis platform for the result of champ.MVP(). The MVP panel indicate parameters user may be used to select significant CpGs, here I only provided abslogFC and p value as two threshold cutoff. After opening this web page, user may select their cutoff, then press submit, the webpage would calculate the result automatically.
    User could check the MVPtable in first tab easily, users can rank and select certain genes in the table, the content of the table might be changed based on the cutoff you selected in left panel.
    The second tab provide the heatmap of all significant CpGs you selected, be careful that if there are too many CpGs, the memory consumption might be large. 
    The third tab provide barplots of proportions of feature and CpGs in for your selected CpGs. 
    The fourth tab is the plot of gene and the wikigene information of certain gene, you may search the gene you want to check by left panel, note that if there is only one significant CpG in the gene you selected, the plot might not be show properly. 
    The last panel provide a boxplot of CpGs and a gene enrichment plot, you may use this gene enrichment plot to find interesting genes.
}

\usage{
    MVP.GUI(MVP=myMVP,
            beta=myNorm,
            pheno=myLoad$pd$Sample_Group)
}

\arguments{
    \item{MVP}{The result from champ.MVP(). (default = myMVP)}
    \item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
}

\value{
    Totally five tabs would be generated on opened webpage.
    \item{MVPtable}{The DMP list of all significant CpGs selected by cutoff in left panel.}
    \item{Heatmap}{Heatmap of all significant CpGs selected by cutoff in left panel.}
    \item{Feature&CpG}{Barplot of feature and Cgi information for all significant CpGs selected by cutoff in left panel.}
    \item{Gene}{Dots and lines of all significant CpGs involved in one gene, the distance between CpGs are equal, and the feature and Cgi information are marked down the plot. Below the plot, is the wikigene information extracted from website.}
    \item{CpG}{Boxplot for CpGs you want to check, you can search CpGs based on the left panel. Below is the gene enrichment plot, hyper CpGs and hyper CpGs are separated.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myMVP <- champ.MVP()
        MVP.GUI()
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
