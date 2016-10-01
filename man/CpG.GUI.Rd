\name{CpG.GUI}
\alias{CpG.GUI}
\title{
    Generate interactive plot for summary information of a list of CpGs
}
\description{
    A Shiny, Plotly and Web Brower based analysis interface. CpG.GUI() is aimed to generate summary of a list of CpGs. Feature distribution, CpG island distribution .e.g. It's call for X11 similar graph software locally if you are doing analysis on server. Also the RAM memory might be large if you have a very big dataset. This function can be used anytime you have a list of CpGs from any analysis, you simply need to imput the CpGs and specify the array type, a web brower interactive interface would be generated automatically. The plots are interactive thus you can make easier and better analysis on your data, and also download them at any size (jpg only).
}

\usage{
    CpG.GUI(CpG=rownames(myLoad$beta),
            arraytype="450K")
}

\arguments{
	\item{CpG}{A list of CpG you want to do plot summary. MUST be a vector with CpG ID. (default = rownames(myLoad$beta))}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}

\value{
    Totally four plots would be generated on opened webpage.
    \item{chromosome_barplot}{A chromosome barplot for the CpG list}
    \item{feature_barplot}{A feature barplot for the CpG list}
    \item{cgi_barplot}{A cgi barplot for the CpG list}
    \item{type_barplot}{A type-I and type-II barplot for the CpG list}
}


\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        CpG.GUI()
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
