\name{champ.impute}
\alias{champ.impute}

\title{
    Conduct imputation for NA value on beta matrix and corresponding pd(Sample_sheet.csv) file.
}
\description{
    champ.impute will conduct imputation on beta matrix contains missing value. This function can be used for any beta dataset, along with their corresponding pd files. If you loaded this file with champ.load(), champ.impute() function will automatically loaded myLoad$beta as inputted beta matrix, while take myLoad$pd as pd input. There are totally three method provided in champ.impute() function. "Delete" is simply remove all NA related CpGs and Samples contain certian proportion of missing value, which is suitable for Small DataSets. "KNN" method use impute.knn() function from "impute" to do imputation on all missing value, which is rather popular but would cause trouble if DataSets contains few samples, no CpGs or samples woule be deleted. "Combine" method would remove all Samples and CpGs with certian proportions of missing value, then do KNN imputation for the rest (Default).
}
\usage{
    champ.impute(beta=myLoad$beta,
                 pd=myLoad$pd,
                 method="Combine",
                 k=5,
                 ProbeCutoff=0.2,
                 SampleCutoff=0.1)
}

\arguments{
	\item{beta}{Data matrix want to be imputed, user can input M matrix or intensity matrix even. (default = myLoad$beta)}
    \item{pd}{Phenotype file for your data set. It's optional for this function, but if during imputation some samples contain too many NA values dicarded, your old pd file might not be able to work for imputed data properly any more. (default = myLoad$pd)}
    \item{method}{Imputation method optional, only "Combine","KNN","Delete" are feasible. (default = "Combine").}
    \item{k}{Number of neighbors to be used in the imputation (default = 5)}
    \item{ProbeCutoff}{Proportion of for probes shall be removed. Any probes with NA value proportion above this parameter will be removed. (default = 0.2)}
    \item{SampleCutoff}{Proportion of for Sample shall be removed. Any Sample with NA value proportion above this parameter will be removed. (default = 0.1)}
}

\value{
	\item{beta}{The matrix get imputed}
    \item{pd}{The pd file corresponding to imputed matrix, if provided.}
}
\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myImpute <- champ.impute()
    }
}
\author{
Yuan Tian
}
