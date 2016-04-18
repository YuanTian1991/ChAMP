\name{champ.DMR}
\alias{champ.DMR}
\title{
    Applying Bumphunter or ProbeLasso Algorithms to detect Different Methylation Regions in a beta valued Methylation Dataset.
}
\description{
    Applying Bumphunter or ProbeLasso Algorithms to Estimate regions for which a genomic profile deviates from its baseline value. Originally implemented to detect differentially methylated genomic regions between two populations. By default, we recommend user do champ.DMR on normalized beta value on two populations, like case to control. The function will return detected DMR and estimated p value. The two algorithms specified in this function is different, while Bumphunter calcuated averaged candidate bumps methylation value between case and control, ProbeLasso need Different Methylated Probes (DMP) from champ.MVP as input parameter and find DMRs around those DMPs. Thus parameters is different for two algorithms.
}

\usage{
    champ.DMR(betaNorm=myNorm$beta,design=myLoad$pd$Sample_Group,maxGap=300,
            cutoff=0.5,minProbes=7,smooth=TRUE,smoothFunction=loessByCluster,
            useWeights=FALSE,permutations=NULL,B=250,pickCutoff=FALSE,
            pickCutoffQ=0.99,nullMethod="bootstrap",verbose=TRUE,cores=3,arraytype="450K",
            method = "Bumphunter",resultsFile=mylimma,meanLassoRadius=375,minSigProbesLasso=3,minDmrSep=1000,
            minDmrSize=50,adjPvalProbe=0.05,adjPvalDmr=0.05,pData=myLoad$pd)

}
\arguments{
    Since there are two function incoporated to detect DMRs, user may specify each function to do DMR detection, Bumphunter or ProbeLasso. Both methods are available for both 450K and EPIC beadarray. But they are controled by different parameters, thus users shall be careful when they specify parameters for corresponding algorithm.
    Parameters shared by two algorithms:
    \item{betaNorm}{Methylation beta valueed dataset user want to do RefFreeEWAS. We recommend to use normalized beta value. In champ.DMR function, beta value will be transformed to M value. NA value is NOT allowed into this function, thus user may need to do some imputation work beforehead. This parameter is essential for both two algorithms.}
    \item{arraytype}{Assign types of manifest shall be used to do DMR detection, "EPIC" or "450K" are provided. Annotation of probes is playing an important role in clustering probes across genome. This parameter is both important and available for both Bumphunter and ProbeLasso functions.}
    \item{method}{Specify the method users want to use to do DMR detection. There are two options: "Bumphunter" or "ProbeLasso". The default value is "Bumphunter".}

    Parameters specified for Bumphunter algorithm:
	\item{design}{Design vector (subjects x covariates). This parameter MUST be a vecter. Though characters, factor and numeric are all allowed because inside the function, character covariates will be transformed into numeric, we still recommend user input numeric deigned covariates vector. This parameter is specified for Bumphunter algorithm.}
    \item{cutoff}{A numeric value. Values of the estimate of the genomic profile above the cutoff or below the negative of the cutoff will be used as candidate regions. It is possible to give two separate values (upper and lower bounds). If one value is given, the lower bound is minus the value. This parameter is specified for Bumphunter.}
    \item{pickCutoff}{Should bumphunter attempt to pick a cutoff using the permutation distribution? This parameter is specified for Bumphunter.}
    \item{pickCutoffQ}{The quantile used for picking the cutoff using the permutation distribution. This parameter is specified for Bumphunter.}
    \item{maxGap}{If cluster is not provided this maximum location gap will be used to define cluster via the clusterMaker function. This parameter is specified for Bumphunter.}
    \item{minProbes}{Threshold to filtering clusters with too few probes in it. After region detection, champ.DMR will only select probes in clusters contain more than minsample probes to continue the program. This parameter is specified for Bumphunter.}
    \item{nullMethod}{Method used to generate null candidate regions, must be one of \sQuote{bootstrap} or \sQuote{permutation} (defaults to \sQuote{permutation}). However, if covariates in addition to the outcome of interest are included in the design matrix (ncol(design)>2), the \sQuote{permutation} approach is not recommended. See vignette and original paper for more information. This parameter is specified for Bumphunter.}
    \item{smooth}{A logical value. If TRUE the estimated profile will be smoothed with the smoother defined by \code{smoothFunction}, This parameter is specified for Bumphunter.}
    \item{smoothFunction}{A function to be used for smoothing the estimate of the genomic profile. Two functions are provided by the package: \code{loessByCluster} and \code{runmedByCluster}. This parameter is specified for Bumphunter.}
    \item{useWeights}{A logical value. If \code{TRUE} then the standard errors of the point-wise estimates of the profile function will be used as weights in the loess smoother \code{loessByCluster}. If the \code{runmedByCluster} smoother is used this argument is ignored. This parameter is specified for Bumphunter.}
    \item{B}{An integer denoting the number of resamples to use when computing null distributions. This defaults to 0. If \code{permutations} is supplied that defines the number of permutations/bootstraps and \code{B} is ignored. The default value is 250, This parameter is specified for Bumphunter.}
    \item{permutations}{ is a matrix with columns providing indexes to be used to scramble the data and create a null distribution when \code{nullMethod} is set to  permutations. If the bootstrap approach is used this argument is ignored. If this matrix is not supplied and \code{B}>0 then these indexes are created using the function \code{sample}. This parameter is specified for Bumphunter.}
    \item{verbose}{logical value. If \code{TRUE}, it writes out some messages indicating progress. If \code{FALSE} nothing should be printed. This parameter is specified for Bumphunter.}
    \item{cores}{The embeded DMR detection function, bumphunter, could automatically use more parallel to accelerate the program. User may assgin number of cores could be used on users's computer. The default value is 3. User may use \code{detectCore()} function to detect number of cores in total. This parameter is specified for Bumphunter.}

    Parameters specified for ProbeLasso algorithm:
    \item{resultsFile}{Different Methylated Probes (DMP) detected from champ.MVP() function, which used limma function to find all CpGs show significant different methylation value. It's a MUST provided parameter for ProbeLasso algorithm.}
    \item{meanLassoRadius}{Radius around each DMP to detect DMR, the default value is 375.}
    \item{minSigProbesLasso}{The minimum number of significant probes to be captured in lasso, default = 3.}
    \item{minDmrSep}{The minimum seperation (bp) between neighbouring DMRs, default = 1000.}
    \item{minDmrSize}{The minimum DMR size (bp), default = 50.}
    \item{adjPvalProbe}{The minimum threshold of significance for probes to be includede in DMRs, default = 0.05}
    \item{adjPvalDmr}{This is the significance threshold for including DMRs in the final DMR list.}
    \item{pData}{This data.frame includes the information from the sample sheet. The default assumes you ran champ.load and saved the output to "myLoad".}
}
\examples{
	myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))	
    myNorm <- champ.norm()
    myDMR <- champ.DMR(B=10)
}

\value{
	\item{myDmrs}{Different Methylation Regions detected by champ.DMR. For different algorithms, myDmrs would be in different structure.}
	\item{myDmrProbes}{Different Methylated Probes detected in DMRs.}
}
\references{
    Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD and Irizarry RA (2014). "Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays." Bioinformatics, 30(10), pp. 1363-1369. http://doi.org/10.1093/bioinformatics/btu049.
    Butcher LM and Beck S (2015). "Probe Lasso: A novel method to rope in differentially methylated regions with 450K DNA methylation data." Methods, 72, pp. 21-28. http://doi.org/10.1016\%2Fj.ymeth.2014.10.036.
}


\author{
    Butcher, L,Aryee MJ, Irizarry RA, Andrew Teschendorff, Yuan Tian
}
\keyword{ DMR }
\keyword{ Bumphunter }
\keyword{ ProbeLasso }
