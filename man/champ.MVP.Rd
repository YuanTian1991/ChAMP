\name{champ.MVP}
\alias{champ.MVP}
\title{
Identify Most Variable Positions in Illumina HumanMethylation450 data.
}
\description{
This function 
}
\usage{
champ.MVP(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.05, adjust.method = "BH", 
compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"), 
bedFile = FALSE)
}

\arguments{
	\item{beta.norm}{A matrix of values representing the methylation scores for each sample (M or B). The default assumes you ran champ.norm and saved the output to "norm"".}
	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
  	\item{adjPVal}{The minimum threshold of significance for probes to be considered an MVP, default = 0.05}
	\item{adjust.method}{The p-value adjustment method to be used for the limma analyis, default= BH (Bonferroni-Hochberg)}
	\item{compare.group}{Not yet implemented}
	\item{resultsDir}{Directory where results will be saved. Default is a folder in the current working directory called "resultsChamp".}
	\item{bedFile}{If bedFile=TRUE, the MVPs will be saved in bedfile format for downstream analysis.}
}
\examples{
	data(testDataSet)
	myLoad=testDataSet
	myNorm=champ.norm(norm="NONE")
}
\value{
\item{results.file}{A matrix of all probes with an adjusted p-value for significance of differential methylation containing columns for probeID, logFC, AveExpr, t, P.Value, adjusted p-value, B, chromosome, map info, chromosome arm, closest gene.1, gene.2, gene.3, gene.4, closest feature.1, feature.2, feature.3, feature.4, UCSC_CpG_ISALNDS_NAME, Relation to UCSC CpG Island, Phantom, DMR, Enhancer, HMM_Island, regulatory feature name, regulatory feature group, feature relation, average of first sample group, average of second sample group, delta beta}
}
\author{
Morris, T
}
\keyword{ limma }
