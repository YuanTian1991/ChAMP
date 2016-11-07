\name{champ.DMP}
\alias{champ.DMP}
\title{
Identify Differential Methylation Positions (DMP) in Illumina HumanMethylation450 or HumanMethylationEPIC data.
}
\description{
This function would use limma package to calculate differential methylation probes between two phenotypes. It's easy to use but do remember to specify which two phenotypes you want to calculate in compare.group parameter. Note that, if the compare.group parameter is NULL, or the factor in it are not find in pheno, the first two phenotypes would be analysed automatically. Note that the result of champ.DMR() would be used as inpute of champ.GSEA() function, thus we suggest user not change the internal structure of the result of champ.DMR() function.
}
\usage{
    champ.DMP(beta = myNorm,
              pheno = myLoad$pd$Sample_Group,
              adjPVal = 0.05,
              adjust.method = "BH",
              compare.group = NULL,
              arraytype = "450K")
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
  	\item{adjPVal}{The minimum threshold of significance for probes to be considered an DMP. (default = 0.05)}
	\item{adjust.method}{The p-value adjustment method to be used for the limma analyis, (default= BH (Benjamini-Hochberg))}
	\item{compare.group}{compare.group is a parameter to assign which two phenotypes you wish to analysis, if it's missed(NULL) or can not fulfill the condition of the dataset, the first two phenotypes in your pheno would be selected as compare.group automatically. (default = NULL)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
}
\value{
\item{DMP}{A data frame of all probes with an adjusted p-value for significance of differential methylation containing columns for logFC, AveExpr, t, P.Value, adj.P.Val, B, C_AVG, T_AVG, deltaBeta, CHR, MAPINFO, Strand, Type, gene, feature, cgi, feat.cgi, UCSC_CpG_Islands_Name, DHS, Enhancer, Phantom, Probe_SNPs, Probe_SNPs_10}
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
    The internal structure of the result of champ.DMP() function should not be modified if it's not necessary caused it would be assigned as inpute for some other functions like DMP.GUI(), champ.DMR() or champ.GSEA().
    You can try to use DMP.GUI() to do interactively analysis on the result of champ.DMP().
}


\author{
Yuan Tian
}
\keyword{ limma }
