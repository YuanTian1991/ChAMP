\name{DMP.GUI}
\alias{DMP.GUI}
\title{
    Generate interactive plot for the result of champ.DMP() function,
}
\description{
    New Modification: We now makes DMP.GUI() would detect numeric variables, which means, if you use champ.DMP() calculated some numeric variable related CpGs, you can continuely use this DMP.GUI() function to draw nice plot for these CpGs and genes even. For CpG plot, while for categorical variables, boxplot would be plotted, we will plot scatter plot for numeric variable (like age) now. For gene plot, we will firstly devide your covariates into couple groups (default is 4), then treat it as a categorical variable. By doing this, you may see that your CpGs should significantly difference lines for difference phenotypes.
    Also, since now champ.DMP() would calculate pairwise comparision for covariate contains more than 2 phenotypes. All result of DMP would be stored into a list, no longer directly myDMP again, so if you have multiple result from champ.DMP(), please inpute each of them here into DMP.GUI(), like DMP.GUI(myDMP[[1]]...), DMP.GUI(myDMP[[2]]...), DMP.GUI(myDMP[[3]]...)
    A Shiny, Plotly and Web Brower based analysis interface. DMP.GUI() is aimed to provide a comprehensive interactive analysis platform for the result of champ.DMP(). The left panel indicate parameters user may be used to select significant CpGs, here I only provided abslogFC and p value as two threshold cutoff. After opening this web page, user may select their cutoff, then press submit, the webpage would calculate the result automatically.
    User could check the DMPtable in first tab easily, users can rank and select certain genes in the table, the content of the table might be changed based on the cutoff you selected in left panel.
    The second tab provide the heatmap of all significant CpGs you selected, be careful that if there are too many CpGs, the memory consumption might be large. 
    The third tab provide barplots of proportions of feature and CpGs in for your selected CpGs. 
    The fourth tab is the plot of gene and the wikigene information of certain gene, you may search the gene you want to check by left panel, note that if there is only one significant CpG in the gene you selected, the plot might not be show properly. 
    The last panel provide a boxplot of CpGs and a gene enrichment plot, you may use this gene enrichment plot to find interesting genes.
}

\usage{
    DMP.GUI(DMP=myDMP[[1]],
            beta=myNorm,
            pheno=myLoad$pd$Sample_Group,
            cutgroupnumber=4)
}

\arguments{
    \item{DMP}{The result from champ.DMP(). (default = myDMP)}
    \item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. In our new upgrading work, DMP.GUI() also accept numeric variables. (default = myLoad$pd$Sample_Group)}
    \item{cutgroupnumber}{This parameters only works if your pheno parameter is a numeric variable, when DMP.GUI() plot gene plot, we will automatically devide your phenotype into couple groups, then treat it as a categorical variable. You may modify this parameter here to tell DMP.GUI how many groups should be devide. Note that this parameter should be setted based on number of value in your pheno parameter. (default = 4)}
}

\value{
    Totally five tabs would be generated on opened webpage.
    \item{DMPtable}{The DMP list of all significant CpGs selected by cutoff in left panel.}
    \item{Heatmap}{Heatmap of all significant CpGs selected by cutoff in left panel.}
    \item{Feature&CpG}{Barplot of feature and Cgi information for all significant CpGs selected by cutoff in left panel.}
    \item{Gene}{Dots and lines of all significant CpGs involved in one gene, the distance between CpGs are equal, and the feature and Cgi information are marked down the plot. Below the plot, is the wikigene information extracted from website.}
    \item{CpG}{Boxplot for CpGs you want to check, you can search CpGs based on the left panel. Below is the gene enrichment plot, hyper CpGs and hyper CpGs are separated.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myDMP <- champ.DMP()
        DMP.GUI()
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
