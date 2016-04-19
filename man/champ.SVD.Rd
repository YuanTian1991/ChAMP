\name{champ.SVD}
\alias{champ.SVD}

\title{
Singular Value Decomposition analysis for batch effects prediciton in HumanMethylation450 data
}

\description{
Runs Singular Value Decomposition on a dataset to estimate the impact of batch effects. 
}
\usage{
champ.SVD(beta = myNorm$beta, rgSet = myLoad$rgSet, detP = myLoad$detP, pd = myLoad$pd, 
loadFile = FALSE, betaFile = "beta.txt", sampleSheet = "sampleSheet.txt", methProfile = FALSE, 
methFile = "MethylationProbeProfile.txt", controlProfile = FALSE, 
controlFile = "ControlProbeProfile.txt", studyInfo = FALSE, studyInfoFile = "studyInfo.txt", 
infoFactor = c(), resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). The default assumes you ran champ.norm and saved the output to "myNorm".}
	\item{rgSet}{An rgSet object that was created when data was loaded the data from the .idat files. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{detP}{A matrix of detection p-values for each sample. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
	\item{loadFile}{If loadFile=TRUE, then the beta values and sample sheet need to be  uploaded}  
	\item{betaFile}{If loadFile=T, }
	\item{sampleSheet}{If the data has not been loaded from .idat files and fromFile=TRUE then this points to the required sampleSheet. Default is "sampleSheet.txt"}
	\item{methProfile}{If methprofile=TRUE then the beta values will be uploaded using the MethylationProbeProfile file from Genome Studio}
	\item{methFile}{If methProfile=TRUE then the beta values will be uploaded using the MethylationProbeProfile from Genome Studio. This is the name of the file. Default is "MethylationProbeProfile.txt"}
	\item{controlProfile}{If rawdata = FALSE and runSVD = TRUE, then it is useful to have a control probe profile file exported from Genome Studio so that internal control probes can be included in the SVD analyis. Default is FALSE.}
	\item{controlFile}{If controlProfile = TRUE then the control probe values will be uploaded using the ControlProbeProfile from Genome Studio. This is the name of the file. Default is "ControlProbeProfile.txt"}
	\item{studyInfo}{If studyInfo=TRUE, additional study covariate information can be included in the SVD analysis. Default is FALSE.}
	\item{infoFactor}{This.}	
	\item{studyInfoFile}{If studyInfo =T, this file will include the additional study information. Default is "studyInfo.txt".}
	\item{resultsDir}{Directory where results will be saved. Default is to create a folder called "resultsChamp"in the current working directory.}
}
\references{
Teschendorff, A. E., Menon, U., Gentry-Maharaj, A., Ramus, S. J., Gayther, S. A., Apostolidou, S., Jones, A., Lechner, M., Beck, S., Jacobs, I. J., and Widschwendter, M. (2009). An epigenetic signature in peripheral blood predicts active ovarian cancer. PLoS One, 4(12), e8274

}
\author{
Teschendorff, A \cr
adapted by Morris, T
}

\keyword{ batch effects }
