\name{champ.DMR}
\alias{champ.DMR}
\title{
    Applying Bumphunter, DMRcate or ProbeLasso Algorithms to detect Different Methylation Regions in a beta valued Methylation Dataset.
}
\description{
    Applying Bumphunter, DMRcate or ProbeLasso Algorithms to Estimate regions for which a genomic profile deviates from its baseline value. Originally implemented to detect differentially methylated genomic regions between two populations. By default, we recommend user do champ.DMR on normalized beta value on two populations, like case to control. The function will return detected DMR and estimated p value. The three algorithms specified in this function is different, while Bumphunter and DMRcate calcuated averaged candidate bumps methylation value between case and control. Thus parameters is different for three algorithms. Note that the result of champ.DMR() would be used as inpute of champ.GSEA() function, thus we suggest user not change the internal structure of the result of champ.DMR() function. 
}

\usage{
    champ.DMR(beta=myNorm,
              pheno=myLoad$pd$Sample_Group,
              compare.group=NULL,
              arraytype="450K",
              method = "Bumphunter",
              minProbes=7,
              adjPvalDmr=0.05,
              cores=3,
              ## following parameters are specifically for Bumphunter method.
              maxGap=300,
              cutoff=NULL,
              pickCutoff=TRUE,
              smooth=TRUE,
              smoothFunction=loessByCluster,
              useWeights=FALSE,
              permutations=NULL,
              B=250,
              nullMethod="bootstrap",
              ## following parameters are specifically for probe ProbeLasso       method.
              meanLassoRadius=375,
              minDmrSep=1000,
              minDmrSize=50,
              adjPvalProbe=0.05,
              Rplot=T,
              PDFplot=T,
              resultsDir="./CHAMP_ProbeLasso/",
              ## following parameters are specifically for DMRcate method.
              rmSNPCH=T,
              fdr=0.05,
              dist=2,
              mafcut=0.05,
              lambda=1000,
              C=2)
}
\arguments{
    Since there are three methods incoporated to detect DMRs, user may specify which function to do DMR detection, Bumphunter DMRcate or ProbeLasso. All three methods are available for both 450K and EPIC beadarray. But they are controled by different parameters, thus users shall be careful when they specify parameters for corresponding algorithm.
    Parameters shared by three algorithms:
    \item{beta}{Methylation beta valueed dataset user want to detect DMR. We recommend to use normalized beta value. In Bumphunter method, beta value will be transformed to M value. NA value is NOT allowed into this function, thus user may need to do some imputation work beforehead. This parameter is essential for both two algorithms. (default = myNorm)}
    \item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for  example "Cancer", "Normal"... Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{compare.group}{ProbeLasso Method does not allow pheno contains more than 2 phenotypes, so if your want use ProbeLasso method, but pheno parameter contains more than 2 phenotypes, you MUST specify compare.group as "compare.group=c("A","B")" to make sure ProbeLasso only works on ONLY two phenotypes. If your pheno parameter contains only 2 phenotypes, you can leave it as NULL. (default=NULL)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
    \item{method}{Specify the method users want to use to do DMR detection. There are three options: "Bumphunter", "DMRcate" or "ProbeLasso". (default = "Bumphunter").}
    \item{minProbes}{Threshold to filtering clusters with too few probes in it. After region detection, champ.DMR will only select DMRs contain more than minProbes to continue the program. (default = 7)}
    \item{adjPvalDmr}{This is the significance threshold for including DMRs in the final DMR list. (default = 0.05)}
    \item{cores}{The embeded DMR detection function, bumphunter and DMRcate, could automatically use more parallel to accelerate the program. User may assgin number of cores could be used on users's computer. User may use \code{detectCore()} function to detect number of cores in total. (default = 3)}

    Parameters specific for Bumphunter algorithm:
    \item{maxGap}{The maximum length for a DMR should be detected, regions longer then this would be discarded. (default = 300)}
    \item{cutoff}{A numeric value. Values of the estimate of the genomic profile above the cutoff or below the negative of the cutoff will be used as candidate regions. It is possible to give two separate values (upper and lower bounds). If one value is given, the lower bound is minus the value. (default = NULL)}
    \item{pickCutoff}{A bool value to indicate if bumphunter algorithm will automatically select the threshold of DMRs. If the value is TRUE, bumphunter will automatically generated 0.99 cutoff from permutation. If user think this threshold is not suitable, user may set their own cutoff here. (default = TRUE)}
    \item{smooth}{A logical value. If TRUE the estimated profile will be smoothed with the smoother defined by \code{smoothFunction}. (default = TRUE)}
    \item{smoothFunction}{A function to be used for smoothing the estimate of the genomic profile. Two functions are provided by the package: \code{loessByCluster} and \code{runmedByCluster}. (default = loessByCluster)}
    \item{useWeights}{A logical value. If \code{TRUE} then the standard errors of the point-wise estimates of the profile function will be used as weights in the loess smoother \code{loessByCluster}. If the \code{runmedByCluster} smoother is used this argument is ignored. (default = FALSE)}
    \item{permutations}{ is a matrix with columns providing indexes to be used to scramble the data and create a null distribution when \code{nullMethod} is set to  permutations. If the bootstrap approach is used this argument is ignored. If this matrix is not supplied and \code{B}>0 then these indexes are created using the function \code{sample}. (default = NULL)}
    \item{B}{An integer denoting the number of resamples to use when computing null distributions. If \code{permutations} is supplied that defines the number of permutations/bootstraps and \code{B} is ignored. (default = 250)}
    \item{nullMethod}{Method used to generate null candidate regions, must be one of \sQuote{bootstrap} or \sQuote{permutation} (defaults to \sQuote{permutation}). However, if covariates in addition to the outcome of interest are included in the design matrix (ncol(design)>2), the \sQuote{permutation} approach is not recommended. See vignette and original paper for more information. (default = "bootstrap")}

    Parameters specific for ProbeLasso algorithm:
    \item{meanLassoRadius}{Radius around each DMP to detect DMR. (default = 375)}
    \item{minDmrSep}{The minimum seperation (bp) between neighbouring DMRs. (default = 1000.)}
    \item{minDmrSize}{The minimum DMR size (bp). (default = 50)}
    \item{adjPvalProbe}{The minimum threshold of significance for probes to be includede in DMRs. (default = 0.05)}
    \item{PDFplot}{If PDFplot would be generated and save in resultsDir. (default = TRUE)}
    \item{Rplot}{If Rplot would be generated and save in resultsDir. Note if you are doing analysis on a server remotely, please make sure the server could connect your local graph applications. (For example X11 for linux.) (default = TRUE)}
    \item{resultsDir}{The directory where PDF files would be saved. (default = "./CHAMP_ProbeLasso/")}

    Parameters specific for Dmrcate algorithm:
    \item{rmSNPCH}{Filters a matrix of M-values (or beta values) by distance to SNP. Also (optionally) removes crosshybridising probes and sex-chromosome probes. (default = TRUE)}
    \item{fdr}{FDR cutoff (Benjamini-Hochberg) for which CpG sites are individually called as significant. Used to index default thresholding in dmrcate(). Highly recommended as the primary thresholding parameter for calling DMRs.}
    \item{dist}{Maximum distance (from CpG to SNP) of probes to be filtered out. See details for when Illumina occasionally lists a CpG-to-SNP distance as being < 0. (default = 2)}
    \item{mafcut}{Minimum minor allele frequency of probes to be filtered out. (default = 0.05)}
    \item{lambda}{Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR bookend definition; gaps >= lambda between significant CpG sites will be in separate DMRs. Support is truncated at 5*lambda. See DMRcate package for further info. (default = 1000)}
    \item{C}{Scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma. Empirical testing shows that when lambda=1000, near-optimal prediction of sequencing-derived DMRs is obtained when C is approximately 2, i.e. 1 standard deviation of Gaussian kernel = 500 base pairs. Cannot be < 0.2. (default = 2)}
}

\value{
	\item{myDmrs}{A data.frame in a list contains Different Methylation Regions detected by champ.DMR. For different algorithms, myDmrs would be in different structure and named as "BumphunterDMR", "DMRcateDMR" and "ProbeLassoDMR". They may contain some different informations, caused by their method. However all three kinds of result are already suitable for champ.GSEA() analysis, so please don't modify the stucture if it's not necessary.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))	
        myNorm <- champ.norm()
        myDMR <- champ.DMR()
        DMR.GUI()
    }
}

\note{
    The internal structure of the result of champ.DMR() function should not be modified if it's not necessary      caused it would be assigned as inpute for some other functions like champ.GSEA().
    You can try to use DMR.GUI() to do interactively analysis on the result of champ.DMR().
}

\note{
    The internal structure of the result of champ.DMR() function should not be modified if it's not necessary caused it would be assigned as inpute for some other functions like DMR.GUI() and champ.GSEA().
    You can try to use DMR.GUI() to do interactively analysis on the result of champ.DMR().
}

\references{
    Jaffe AE et a. Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies. Int J Epidemiol. 2012;41(1):200-209.

    Butcher LM, Beck S. Probe lasso: A novel method to rope in differentially methylated regions with 450K dna methylation data. Methods. 2015;72:21-28.

    Peters TJ, Buckley MJ, Statham AL, et al. De novo identification of differentially methylated regions in the human genome. Epigenetics & Chromatin. 2015;8(1):1-16.
}

\author{
    Butcher, L,Aryee MJ, Irizarry RA, Andrew Teschendorff, Yuan Tian
}
\keyword{ DMR }
