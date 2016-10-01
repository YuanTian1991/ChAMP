\name{champ.EpiMod}
\alias{champ.EpiMod}
\title{
    champ.EpiMod() infer differential methylation interactome hotspots. 
}
\description{
    champ.EpiMod() used FEM package to identify interactome hotspots of differential promoter methylation. By "interactome hotspot" we mean a connected subnetwork of the protein interaction network (PIN) with an exceptionally large average edge-weight density in relation to the rest of the network. The weight edges are constructed from the statistics of association of DNA methylation with the phenotype of interest. Thus, the EpiMod algorithm can be viewed as a functional supervised algorithm, which uses a network of relations between genes (in our case a PPI network), to identify subnetworks where a significant number of genes are associated with a phenotype of interest (POI). We call these "hotspots" also Functional Epigenetic Modules (FEMs). You can get more detailed information in FEM package.
}

\usage{
    champ.EpiMod(beta=myNorm,
                 pheno=myLoad$pd$Sample_Group,
                 nseeds=100,
                 gamma=0.5,
                 nMC=1000,
                 sizeR.v=c(1,100),
                 minsizeOUT=10,
                 resultsDir="./CHAMP_EpiMod/",
                 PDFplot=TRUE,
                 arraytype="450K")
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{nseeds}{An integer specifying the number of seeds and therefore modules to search for. (default = 100)}
    \item{gamma}{A parameter of the spin-glass algorithm, which determines the average module size. Default value generally leads to modules in the desired size range (10-100 genes). (default = 0.5)}
    \item{nMC}{Number of Monte Carlo runs for establishing statistical significance of modularity values under randomisation of the molecular profiles on the network. (default = 1000)}
    \item{sizeR.v}{Desired size range for modules. (default = c(1,100))}
    \item{minsizeOUT}{ Minimum size of modules to report as interesting. (default = 10)}
    \item{resultsDir}{The directory where PDF files would be saved. (default = "./CHAMP_QCimages/")}
    \item{PDFplot}{If PDFplot would be generated and save in resultsDir. (default = TRUE)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}
\value{
\item{EpiMod.o}{A data frame of all probes with an adjusted p-value for significance of differential methylation containing columns for logFC, AveExpr, t, P.Value, adj.P.Val, B, C_AVG, T_AVG, deltaBeta, CHR, MAPINFO, Strand, Type, gene, feature, cgi, feat.cgi, UCSC_CpG_Islands_Name, DHS, Enhancer, Phantom, Probe_SNPs, Probe_SNPs_10, you can turn to FEM package for more informations.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myEpiMod <- champ.EpiMod()
        }
}

\references{
    Y J, M W and AE T (2014). A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control. Bioinformatics, 30(16), pp. 2360-2366.
}

\author{
    Teschendorff, A\cr
    adapted by Yuan Tian
}
\keyword{ EpiMod }
