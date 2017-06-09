\name{champ.runCombat}
\alias{champ.runCombat}

\title{
Function that uses ComBat to correct for batch effects. Multiple batch effect correction is allowed.
}
\description{
This function formats data to run through ComBat batch correction. If beta values are used the data is first logit transformed. Them Combat function from "sva" package would be used to do batch correction. Note that multi-batch correction is supported, user just need to assign name of batch need to be corrected. Note Combat function is a little bit critical to dataset, thus you have futher question or higher lever of application of Combat, you may turn to "sva" pacakge for help. After inputing pd file, champ.runCombat() would automatically detect all correctable factors and list them below, if your assigned batchname is correct, champ.runCombat() would start to do batch correction. Note that in new version champ.runCombat() function, we will check if user's variable and batch confound with each other.
}

\usage{
    champ.runCombat(beta=myNorm,
                    pd=myLoad$pd,
                    variablename="Sample_Group",
                    batchname=c("Slide"),
                    logitTrans=TRUE)
}

\arguments{
  	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). (default = myNorm).}
  	\item{pd}{This data.frame includes the information from the sample sheet. (default = myLoad$pd).}
  	\item{variablename}{Variable name which batch should be corrected for, in previous version of ChAMP, variablename was "Sample_Group". (default = "Sample_Group").}
    \item{batchname}{A character vector of name indicates which batch factors shall be corrected. (default = c("Slide"))}
    \item{logitTrans}{If logitTrans=T then your data will be logit transformed before the Combat correction and inverse logit transformed after correction. This is T by default for Beta values but if you have selected M values, it should be FALSE. It is also FALSE when used with CNA as those are intensity values that don't need to be transformed.
    }
}

\value{
	\item{beta}{The matrix of values represeting the methylation scores for each sample after ComBat batch correction.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        champ.SVD()
        myCombat <- champ.runCombat()
        }
}

\references{
    Johnson WE et a. Adjusting batch effects in microarray expression data using empirical bayes methods. Biostatistics. 2007;8(1):118-127.
}

\author{
    Yuan Tian
}

\keyword{ ComBat }
