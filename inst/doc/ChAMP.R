### R code from vignette source 'ChAMP.Rnw'

###################################################
### code chunk number 1: InstallLibraries (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c('minfi', 'DNAcopy', 'impute', 'marray', 'limma',
## 'preprocessCore', 'RPMM', 'sva', 'IlluminaHumanMethylation450kmanifest',
## 'wateRmelon','isva','quadprog','bumphunter','doParallel',
## 'qvalue','RefFreeEWAS','GenomicRanges','plyr'))


###################################################
### code chunk number 2: loadChAMPLibrary
###################################################
library(ChAMP)


###################################################
### code chunk number 3: loadTest
###################################################
testDir=system.file("extdata",package="ChAMPdata")


###################################################
### code chunk number 4: loadTest2
###################################################
data(testDataSet)
myLoad=testDataSet


###################################################
### code chunk number 5: processFunctionIDAT (eval = FALSE)
###################################################
## champ.process(directory = testDir)


###################################################
### code chunk number 6: processSAVE (eval = FALSE)
###################################################
## save(myLoad,file="currentStudyloadedData.RData")
## load("currentStudyloadedData.RData")
## champ.process(fromIDAT=FALSE)


###################################################
### code chunk number 7: processSTEPbySTEP (eval = FALSE)
###################################################
## myLoad <- champ.load(directory = testDir)
## myNorm <- champ.norm()
## champ.SVD()
## batchNorm <- champ.runCombat()
## limma <- champ.MVP()
## myDMR <- champ.DMR()
## myRefBase <- champ.refbase()
## myRefFree <- champ.reffree()
## champ.CNA()


###################################################
### code chunk number 8: EPICprocessSTEPbySTEP (eval = FALSE)
###################################################
## # myLoad <- champ.load(directory = testDir,arraytype="EPIC")
## # We simulated EPIC data from beta value instead of .idat file,
## # but user may use above code to read .idat files directly. 
## # Here we we started with myLoad.
## data(EPICSimData)
## 
## myNorm <- champ.norm(arraytype="EPIC")
## champ.SVD()
## batchNorm <- champ.runCombat()
## myrefbase <- champ.refbase()
## myreffree <- champ.reffree()
## limma <- champ.MVP(arraytype="EPIC")
## myDMR <- champ.DMR(arraytype="EPIC")
## 
## # champ.CNA(arraytype="EPIC")
## # champ.CNA() function call for intensity data, which is not included in 
## # out Simulation data.


###################################################
### code chunk number 9: SampleSheet
###################################################
myLoad$pd


###################################################
### code chunk number 10: load (eval = FALSE)
###################################################
## myLoad=champ.load(directory = testDir, filterBeads=TRUE)


###################################################
### code chunk number 11: loadFunction
###################################################
myLoad = champ.load(directory=testDir)


###################################################
### code chunk number 12: normFunction
###################################################
myNorm=champ.norm()


###################################################
### code chunk number 13: svdFunction
###################################################
champ.SVD()


###################################################
### code chunk number 14: combatFunction (eval = FALSE)
###################################################
## batchNorm=champ.runCombat()


###################################################
### code chunk number 15: mvpFunction (eval = FALSE)
###################################################
## limma=champ.MVP()
## head(limma)


###################################################
### code chunk number 16: lassoFunction (eval = FALSE)
###################################################
## lasso <- champ.DMR(resultFiles=limma,method="ProbeLasso",arraytype="450K")
## bump <- champ.DMR(method="Bumphunter",arraytype="450K")
## if(!is.null(lasso))
## {
## head(lasso)
## }


###################################################
### code chunk number 17: cnaFunction (eval = FALSE)
###################################################
## CNA=champ.CNA()


###################################################
### code chunk number 18: RefFunction (eval = FALSE)
###################################################
## refbase <- champ.refbase()
## reffree <- champ.reffree()


