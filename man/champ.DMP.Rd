\name{champ.DMP}
\alias{champ.DMP}
\title{
Identify Differential Methylation Positions (DMP) and Numeric Covariate related CpGs in Illumina HumanMethylation450 or HumanMethylationEPIC data.
}
\description{
    New modification: champ.DMP() can now find numeric variable related CpGs, and do pairwise comparision between more than 2 phenotypes' covariate.
This function would use limma package to calculate differential methylation probes between two phenotypes. Or use linear regression model to calcuate CpGs related with certain variables. Now in new version champ.DMP() we still have compare.group parameter, but if compare.group is NULL, and user's pheno variable contains more than 2 phenotypes, champ.DMP() would calculate pairwise DMP between each pair of them. Note that the result of champ.DMP() would be used as inpute of champ.GSEA() and DMP.GUI() function, thus we suggest user not change the internal structure of the result of champ.DMP() function.
}
\usage{
    champ.DMP(beta = myNorm,
              pheno = myLoad$pd$Sample_Group,
              compare.group = NULL,
              adjPVal = 0.05,
              adjust.method = "BH",
              arraytype = "450K",
	      resultsDir = "./CHAMP_DMPimages/")
              man.plot = T,
              probes = c(),
              chr = c(),
              dotsize = 1,
              sug.line = F,
              gen.line = T,
              gen.alpha = 0.05,
              q.plot = T))
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{pheno}{Covariate that you want to do analysis, it might be a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. Or it can be a numeric variable like age. (default = myLoad$pd$Sample_Group)}
    \item{compare.group}{If your pheno is categorical variable, you may specify this parameter to ask champ.DMP() only compare certain two phenotypes. If your pheno contains more than 2 phenotypes, and compare.group is NULL, pairwise comparision would be done between each two phenotypes. You may set the value as compare.group=c("C","T"), it must be a vector contains only two charactor element. (default = NULL)}
  	\item{adjPVal}{The minimum threshold of significance for probes to be considered an DMP. (default = 0.05)}
	\item{adjust.method}{The p-value adjustment method to be used for the limma analyis, (default= BH (Benjamini-Hochberg))}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
	\item{resultsDir}{Choose where plots are saved (default = "./CHAMP_DMPimages/")}
	\item{man.plot}{Should Manhattan plot be produced (default = TRUE)}
	\item{probes}{For Manhattan plot, choose probes to highlight, e.g., c("cg26090107", "cg05970496") (default = empty)}
	\item{chr}{For Manhattan plot, choose chromosomes to plot. Input must be numerical, e.g., c(1,4,6) (default = all chromosomes)}
	\item{dotsize}{For Manhattan plot, choose dotsize (between 0 and 1) (default = 1)}
	\item{sug.line}{For Manhattan plot, print a suggestive threshold line. Input must be a numerical, e.g., 6e-08, or FALSE (default = FALSE)}
	\item{gen.line}{For Manhattan plot, print a genome-wide threshold line (default = TRUE)}
	\item{gen.alpha}{For Manhattan plot, choose genome-wide significance level (default = 0.05)}
	\item{q.plot}{Should Q-Q plot be produced (default = TRUE)}
}
\value{
\item{DMP}{A list DMP results. Each element in the list is a data frame of all probes with an adjusted p-value for significance of differential methylation containing columns for logFC, AveExpr, t, P.Value, adj.P.Val, B, C_AVG, T_AVG, deltaBeta, CHR, MAPINFO, Strand, Type, gene, feature, cgi, feat.cgi, UCSC_CpG_Islands_Name, DHS, Enhancer, Phantom, Probe_SNPs, Probe_SNPs_10. These values are directly calculated from limma package, user may read limma manual for more information. deltaBeta is the same as logFC, we kept it here cause maybe old users would stil using it. XXX_AVG is mean value of XXX pheno type in your pheno parameter. Note for numeric variables, the returned result will be named as "NumericVariable", it contains most features as output for categorical covariates except for XXX_AVG and deltaBeta. Furthermore, function produces Manhattan plot as known from GWAS studies, and Q-Q plot to check for inflated p-values.}
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

\references{
    Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47

    Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016). Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. Annals of Applied Statistics 10(2), 946-963.

}

\author{
Yuan Tian
}
\keyword{ limma }
