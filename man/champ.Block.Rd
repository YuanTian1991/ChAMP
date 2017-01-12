\name{champ.Block}
\alias{champ.Block}
\title{
Identify Differential Methylation Blocks in Illumina HumanMethylation450 or HumanMethylationEPIC data.
}
\description{
    This function would detect all methylation Blocks exist in your dataset, methylation Block should be calculated based on the average value of clusters across whole genome. Firstly champ.Block would calculate all clusters in the dataset with clustermaker() function provided by Bumphunter package. Then, only OpenSea Clusters would be picked out to calculate Block. Block can be seen as "large clusters" generated from all small OpenSea Clusters. The algrorithms is similar to the normal DMR-detection one. We will firstly collapse all OpenSea Clusters (or to say regions) into one dot on genome, using average beta value to represent their beta value, and using average position to represent their position. Then we do clustering on these collapsed regions with Bumphunter algrorithms but bigger ranges.
}

\usage{
    champ.Block(beta=myNorm,
                pheno=myLoad$pd$Sample_Group,
                arraytype="450K",
                maxClusterGap=250000,
                B=500,
                bpSpan=250000,
                minNum=10,
                cores=3)
}

\arguments{
	\item{beta}{A matrix of values representing the methylation scores for each sample (M or B). Better to be imputed and normalized data. (default = myNorm)}
	\item{pheno}{This is a categorical vector representing phenotype of factor wish to be analysed, for example "Cancer", "Normal"...  Tow or even more phenotypes are allowed. (default = myLoad$pd$Sample_Group)}
    \item{arraytype}{Choose microarray type is 450K or EPIC. (default = "450K")}
    \item{maxClusterGap}{ Max gap between clusters when calculating region at first step. (default = 250000)}
    \item{B}{An integer denoting the number of resamples to use when computing null distributions. If \code{permutations} is supplied that defines the number of permutations/bootstraps and \code{B} is ignored. (default = 250)}
    \item{bpSpan}{The maximum length for a Block should be detected, regions longer then this would be discarded. (default = 250000)}
    \item{minNum}{Threshold to filtering Blocks with too few probes in it. After region detection, champ.Block will only select Blocks contain more than minNum clusters(OpenSea Regions) to continue the program. (default = 10)}
    \item{cores}{The embeded DMR detection function, bumphunter, could automatically use more parallel to accelerate the program. User may    assgin number of cores could be used on users's computer. User may use \code{detectCore()} function to detect number of cores in total.       (default = 3)}
}
\value{
    \item{Block}{A data.frame contains all detected Blocks, with colnames as chr, start, end, value, area, cluster, indexStart, indexEnd, L, clusterL, p.value, fwer, p.valueArea, fwerArea. The result format is actually the same as Bumphunter, you may refer to Bumphunter packages to get more explaination about the result.}
    \item{clusterInfo}{When champ.Block() detection significant Blocks, a group of candidate Blocks would be detected out at first, this is the data frame of all candidate Blocks. The "TRUE" Blocks in above value are located in these candidate Blocks.} 
    \item{allCLID.v}{The first step of detectiong methylation Blocks is to get each probes into a cluster(region). This value is the clustering result of each probes.}
    \item{avbetaCL.m}{The beta matrix for each cluster. The value is calculated by taking mean value of all probes located in each cluster.}
    \item{posCL.m}{Position of each cluster, which is calculated by average all probes' position in each cluster.}
}

\examples{
    \dontrun{
        myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
        myNorm <- champ.norm()
        myBlock <- champ.Block()
        Block.GUI()
        }
}

\note{
    The internal structure of the result of champ.Block() function should not be modified if it's not necessary caused it would be assigned as inpute for some other functions like Block.GUI().
    You can try to use Block.GUI() to do interactively analysis on the result of champ.Block().
}

\references{
    Hansen KD, Timp W, Bravo HC, et al. Increased methylation variation in epigenetic domains across cancer types. Nat Genet. 2011;43(8):768-775.
}

\author{
Yuan Tian
}
\keyword{ Block }
