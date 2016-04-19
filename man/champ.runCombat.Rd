\name{champ.runCombat}
\alias{champ.runCombat}

\title{
Function that uses ComBat to correct for batch effects related to slide/BeadChip.
}
\description{
This function formats data to run through ComBat batch correction. If beta values are used the data is first logit transformed. 
}

\usage{
champ.runCombat(beta.c = myNorm$beta, pd = myLoad$pd, logitTrans = TRUE)
}

\arguments{
  	\item{beta.c}{A matrix of values representing the methylation scores for each sample (M or B). The default assumes you ran champ.norm and saved the output to "norm".}
  	\item{pd}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
    \item{logitTrans}{If logitTrans=T then your data will be logit transformed before the Combat correction and inverse logit transformed after correction. This is T by default for Beta values but if you have selected M values it will revert to False. It is also False when used with CNA as those are intensity values that don't need to be transformed.
    	}
}
\value{
	\item{beta}{The matrix of values represeting the methylation scores for each sample after ComBat batch correction.}
}
\author{
T. Morris
}

\keyword{ ComBat }

