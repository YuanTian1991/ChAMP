## ----eval=FALSE----------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest","limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend","IlluminaHumanMethylationEPICmanifest","FEM","matrixStats"))

## ----eval=TRUE,message=FALSE, warning=FALSE------------------------------
library("ChAMP")

## ----eval=FALSE----------------------------------------------------------
#  testDir=system.file("extdata",package="ChAMPdata")
#  myLoad <- champ.load(testDir,arraytype="450K")

## ----eval=FALSE----------------------------------------------------------
#  data(EPICSimData)

## ---- out.width = 800, fig.retina = NULL,echo=F--------------------------
knitr::include_graphics("Figure/ChAMP_Pipeline.png")

## ----eval=FALSE----------------------------------------------------------
#  champ.process(directory = testDir)

## ----eval=FALSE----------------------------------------------------------
#  myLoad <- cham.load(testDir)
#  CpG.GUI()
#  champ.QC() # Alternatively: QC.GUI()
#  myNorm <- champ.norm()
#  champ.SVD()
#  # If Batch detected, run champ.runCombat() here.
#  myDMP <- champ.DMP()
#  DMP.GUI()
#  myDMR <- champ.DMR()
#  DMR.GUI()
#  myBlock <- champ.Block()
#  Block.GUI()
#  myGSEA <- champ.GSEA()
#  myEpiMod <- champ.EpiMod()
#  myCNA <- champ.CNA()
#  myRefFree <- champ.reffree()
#  # If DataSet is Blood samples, run champ.refbase() here.

## ----eval=FALSE----------------------------------------------------------
#  # myLoad <- champ.load(directory = testDir,arraytype="EPIC")
#  # We simulated EPIC data from beta value instead of .idat file,
#  # but user may use above code to read .idat files directly.
#  # Here we we started with myLoad.
#  
#  data(EPICSimData)
#  CpG.GUI(arraytype="EPIC")
#  champ.QC() # Alternatively QC.GUI(arraytype="EPIC")
#  myNorm <- champ.norm(arraytype="EPIC")
#  champ.SVD()
#  # If Batch detected, run champ.runCombat() here.This data is not suitable.
#  myDMP <- champ.DMP(arraytype="EPIC")
#  DMP.GUI()
#  myDMR <- champ.DMR()
#  DMR.GUI()
#  myDMR <- champ.DMR(arraytype="EPIC")
#  DMR.GUI(arraytype="EPIC")
#  myBlock <- champ.Block(arraytype="EPIC")
#  Block.GUI(arraytype="EPIC") # For this simulation data, not Differential Methylation Block is detected.
#  myGSEA <- champ.GSEA(arraytype="EPIC")
#  myEpiMod <- champ.EpiMod(arraytype="EPIC")
#  myRefFree <- champ.reffree()
#  
#  # champ.CNA(arraytype="EPIC")
#  # champ.CNA() function call for intensity data, which is not included in our Simulation data.

## ----eval=FALSE----------------------------------------------------------
#  library("doParallel")
#  detectCores()

## ----eval=FALSE----------------------------------------------------------
#  myLoad <- champ.load(testDir)
#  ## We are not running this code here because it cost about 1 minute.

## ----eval=TRUE-----------------------------------------------------------
data(testDataSet)

## ----eval=TRUE-----------------------------------------------------------
myLoad$pd

## ----eval=FALSE----------------------------------------------------------
#  CpG.GUI(CpG=rownames(myLoad$beta),arraytype="450K")

## ---- out.width = 800, fig.retina = NULL,echo=F--------------------------
knitr::include_graphics("Figure/CpGGUI.png")

## ----eval=TRUE,dpi=100,fig.width=7,fig.height=4,message=FALSE------------
champ.QC()

## ----eval=FALSE----------------------------------------------------------
#  QC.GUI(CpG=rownames(myLoad$beta),arraytype="450K")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/QCGUI.jpg")

## ----eval=FALSE----------------------------------------------------------
#  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/BMIQ.jpg")

## ----eval=TRUE,dpi=100,fig.width=8,fig.height=8,message=FALSE,warning=FALSE----
champ.SVD(beta=myNorm,pd=myLoad$pd)

## ----eval=FALSE----------------------------------------------------------
#  myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))

## ----eval=TRUE,warning=FALSE,message=FALSE-------------------------------
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group)

## ----eval=TRUE-----------------------------------------------------------
head(myDMP)

## ----eval=FALSE----------------------------------------------------------
#  DMP.GUI(DMP=myDMP,beta=myNorm,pheno=myLoad$pd$Sample_Group)

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMP-1.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMP-2.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMP-3.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMP-4.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMP-5.png")

## ----eval=FALSE,message=FALSE,warning=TRUE-------------------------------
#  myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter")

## ----eval=TRUE-----------------------------------------------------------
head(myDMR$DMRcateDMR)

## ----eval=FALSE----------------------------------------------------------
#  DMR.GUI(DMR=myDMR)
#  # It might be a little bit slow to open DMR.GUI() because function need to extract annotation for CpGs from DMR. Might take 30 seconds.

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMR-1.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMR-2.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMR-3.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/DMR-4.png")

## ----eval=FALSE----------------------------------------------------------
#  myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450K")

## ----eval=TRUE-----------------------------------------------------------
head(myBlock$Block)

## ----eval=FALSE----------------------------------------------------------
#  Block.GUI(Block=myBlock,beta=myNorm,pheno=myLoad$pd$Sample_Group,runDMP=TRUE,compare.group=NULL,arraytype="450K")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/Block-1.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/Block-2.png")

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/Block-3.png")

## ----eval=FALSE----------------------------------------------------------
#  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP,DMR=myDMR,arraytype="450K",adjPval=0.05)
#  # myDMP and myDMR could (not must) be used directly.

## ----eval=TRUE-----------------------------------------------------------
head(myGSEA$DMP)
# Above is the GSEA result for differential methylation probes.
head(myGSEA$DMR)
# Above is the GSEA result for differential methylation regions.
# Too many information may be printed, so we are not going to show the result here.

## ----eval=FALSE----------------------------------------------------------
#  myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/EpiMod.jpg")

## ----eval=FALSE----------------------------------------------------------
#  myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)

## ---- out.width = 800, fig.retina = NULL,echo=FALSE----------------------
knitr::include_graphics("Figure/CNAGroupPlot.jpg")

## ----eval=FALSE----------------------------------------------------------
#  myRefFree <- champ.reffree(beta=myNorm,pheno=myLoad$pd$Sample_Group)

## ----eval=TRUE-----------------------------------------------------------
head(myRefFree$qvBeta)

## ----eval=FALSE----------------------------------------------------------
#  myRefBase <- champ.refbase(beta=myNorm,arraytype="450K")
#  # Our test data set is not blood.

