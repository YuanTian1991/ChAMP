champ.process <-
function(fromIDAT=TRUE, fromFile=FALSE,directory = getwd(), resultsDir=paste(getwd(),"resultsChamp",sep="/"), methValue="B", filterDetP = TRUE, detPcut=0.01, filterXY = TRUE, removeDetP = 0, filterBeads = TRUE, beadCutoff = 0.05, filterNoCG= FALSE, QCimages = TRUE, batchCorrect=TRUE, runSVD = TRUE, studyInfo=FALSE, infoFactor=c(), norm = "BMIQ", adjust.method="BH", adjPVal=0.05, runDMR=TRUE, runCNA=TRUE,plotBMIQ=FALSE,DMRpval=0.05, sampleCNA=TRUE,plotSample=TRUE,groupFreqPlots=TRUE,freqThreshold=0.3,bedFile=FALSE, methProfile=FALSE,controlProfile=FALSE,arraytype="450K")
{
    batchDone=FALSE

	if(!file.exists(resultsDir))
	{
        dir.create(resultsDir)
		message("Creating results directory. Results will be saved in ", resultsDir)
        #log<-print("Creating results directory. Results will be saved in ", resultsDir)
    }
	if(fromIDAT)
    {
        myLoad=champ.load(directory=directory, methValue=methValue, resultsDir=resultsDir, QCimages=QCimages, filterDetP=filterDetP, detPcut=detPcut,filterXY=filterXY,filterBeads=filterBeads,beadCutoff=beadCutoff,filterNoCG=filterNoCG,arraytype=arraytype)
	}else if(fromFile){
		myLoad=champ.read()
	}else{}

    #####NORMALIZATION########
    #options are BMIQ, SWAN, PBC or NONE
    
    myNorm=champ.norm(beta=myLoad$beta, rgSet=myLoad$rgSet, pd=myLoad$pd, mset=myLoad$mset, methValue=methValue, plotBMIQ=plotBMIQ, QCimages=QCimages, norm=norm, resultsDir=resultsDir,arraytype=arraytype)
	
	if(runSVD)
	{
        champ.SVD(beta=myNorm$beta, pd=myLoad$pd, rgSet=myLoad$rgSet, detP=myLoad$detP, studyInfo=studyInfo, infoFactor=infoFactor,resultsDir=resultsDir)
	}
    
	if(batchCorrect)
	{
        if(methValue=="B")
        {
            batchNorm=champ.runCombat(beta.c=myNorm$beta, pd=myLoad$pd, logitTrans=TRUE)
        }else{
            batchNorm=champ.runCombat(beta.c=myNorm$beta, pd=myLoad$pd, logitTrans=FALSE)
        }
        batchDone=TRUE
	}
	
	if(runDMR)
	{
        dmr <- champ.DMR(arraytype=arraytype)
    }
    
	if(runCNA)
	{
		champ.CNA(intensity=myLoad$intensity,pd=myLoad$pd,batchCorrect=batchCorrect,resultsDir=resultsDir,sampleCNA=sampleCNA,plotSample=plotSample,groupFreqPlots=groupFreqPlots,freqThreshold=freqThreshold,arraytype=arraytype)
	}
	message("ChAMP complete.")
	
}
