champ.process <-
function(fromIDAT=TRUE, fromFile=FALSE,directory = getwd(), resultsDir=paste(getwd(),"resultsChamp",sep="/"), methValue="B", filterDetP = TRUE, detPcut=0.01, filterXY = TRUE, removeDetP = 0, filterBeads = TRUE, beadCutoff = 0.05, filterNoCG= FALSE, QCimages = TRUE, batchCorrect=TRUE, runSVD = TRUE, studyInfo=FALSE, infoFactor=c(), norm = "BMIQ", adjust.method="BH", adjPVal=0.05, runDMR=TRUE, runCNA=TRUE,plotBMIQ=FALSE,DMRpval=0.05, sampleCNA=TRUE,plotSample=TRUE,groupFreqPlots=TRUE,freqThreshold=0.3,bedFile=FALSE, methProfile=FALSE,controlProfile=FALSE)
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
        myLoad=champ.load(directory=directory, methValue=methValue, resultsDir=resultsDir, QCimages=QCimages, filterDetP=filterDetP, detPcut=detPcut,filterXY=filterXY,filterBeads=filterBeads,beadCutoff=beadCutoff,filterNoCG=filterNoCG)
	}else if(fromFile){
		myLoad=champ.read()
	}else{}

    #####NORMALIZATION########
    #options are BMIQ, SWAN, PBC or NONE
    
    myNorm=champ.norm(beta=myLoad$beta, rgSet=myLoad$rgSet, pd=myLoad$pd, mset=myLoad$mset, methValue=methValue, plotBMIQ=plotBMIQ, QCimages=QCimages, norm=norm, resultsDir=resultsDir)
	
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
		if(batchDone==F)
		{
			dmr=champ.lasso(pd=myLoad$pd,beta.norm=myNorm$beta,resultsDir=resultsDir,bedFile=bedFile, adjust.method = adjust.method, adjPVal=adjPVal, DMRpval=DMRpval, batchDone=batchDone)
		}else{
			dmr=champ.lasso(pd=myLoad$pd, beta.norm=batchNorm$beta, resultsDir=resultsDir, bedFile=bedFile, batchDone=batchDone, adjust.method=adjust.method, adjPVal=adjPVal, DMRpval=DMRpval,normSave=myNorm$beta)
		}	
    }
    
	if(runCNA)
	{
		champ.CNA(intensity=myLoad$intensity,pd=myLoad$pd,batchCorrect=batchCorrect,resultsDir=resultsDir,sampleCNA=sampleCNA,plotSample=plotSample,groupFreqPlots=groupFreqPlots,freqThreshold=freqThreshold)
	}
	message("ChAMP complete.")
	
}
