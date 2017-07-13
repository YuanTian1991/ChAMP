champ.process <- function(runload=TRUE,
                          directory = getwd(),
                          filters=c("XY","DetP","Beads","NoCG","SNP","MultiHit"),
                          #---champ.impute parameters below---#
                          runimpute=TRUE,
                          imputemethod="Combine",
                          #---champ.QC parameters below---#
                          runQC=TRUE,
                          QCplots=c("mdsPlot","densityPlot","dendrogram"),
                          #---champ.norm parameters below---#
                          runnorm=TRUE,
                          normalizationmethod="BMIQ",
                          #---champ.SVD parameters below---#
                          runSVD=TRUE,
                          RGEffect=FALSE,
                          #---champ.runCombat parameters below---#
                          runCombat=TRUE,
                          batchname=c("Slide"),
                          #---champ.DMP parameters below---#
                          runDMP=TRUE,
                          #---champ.DMR parameters below---#
                          runDMR=TRUE,
                          DMRmethod="Bumphunter",
                          #---champ.Block parameters below---#
                          runBlock=TRUE,
                          #---champ.GSEA parameters below---#
                          runGSEA=TRUE,
                          #---champ.EpiMod parameters below---#
                          runEpiMod=TRUE,
                          #---champ.CNA parameters below---#
                          runCNA=TRUE,
                          control=TRUE,
                          controlGroup="champCtls",
                          #---champ.refbase parameters below---#
                          runRefBase=FALSE,
                          #---universal settings---#
                          compare.group=NULL,
                          adjPVal=0.05,
                          resultsDir="./CHAMP_RESULT/",
                          arraytype="450K",
                          PDFplot=TRUE,
                          Rplot=TRUE,
                          cores=3,
                          saveStepresults=TRUE)
{
    message("[===========================]")
    message("[<<< ChAMP.PROCESS START >>>]")
    message("-----------------------------")

    message("This is champ.process() function, which would run ALL analysis ChAMP package can be done on your package. But sometimes it's not always very smoothly to conduct all function because of mistake in dataset or unsignificant results. This if you encounter problem in champ.proces(), you shall try run champ step by step.^_^")


    if(!file.exists(resultsDir))
    {
        dir.create(resultsDir)
        message("Creating results directory. All Results will be saved in ", resultsDir)
    }

    CHAMP.RESULT <- list()

    ### Applying champ.load() function.
    if(runload)
    {
        myfilter <- rep(FALSE,6)
        names(myfilter) <- c("XY","DetP","Beads","NoCG","SNP","MultiHit")
        myfilter[filters] <- TRUE

        message("\nRunning champ.load()...")
        myLoad <- champ.load(directory=directory,
                             methValue="B",
                             autoimpute=TRUE,
                             filterXY=myfilter["XY"],
                             filterDetP=myfilter["DetP"],
                             ProbeCutoff=0,
                             SampleCutoff=0.1,
                             detPcut=0.01,
                             filterBeads=myfilter["Beads"],
                             beadCutoff=0.05,
                             filterNoCG=myfilter["NoCG"],
                             filterSNPs=myfilter["SNP"],
                             filterMultiHit=myfilter["MultiHit"],
                             arraytype=arraytype)
        if(saveStepresults)
        {
            save(myLoad,file=paste(resultsDir,"/myLoad.rda",sep=""))
            message("champ.load()'s result \"myLoad\" has been saved in ",resultsDir," as \"myLoad.rda.\"")
        }
        message("Run champ.load() Over!\n")
        gc()
        CHAMP.RESULT[["champ.load"]] <- myLoad
    }

    tmpbeta <- myLoad$beta
    tmppd <- myLoad$pd

    ### Applying champ.impute() function.
    if(runimpute)
    {
        message("\nRunning champ.impute()...")
        myImpute <- champ.impute(beta=myLoad$beta,
                                 pd=myLoad$pd,
                                 method=imputemethod,
                                 k=5,
                                 ProbeCutoff=0.2,
                                 SampleCutoff=0.1)
        if(saveStepresults)
        {
            save(myImpute,file=paste(resultsDir,"/myImpute.rda",sep=""))
            message("champ.impute()'s result \"myImpute\" has been saved in ",resultsDir," as \"myImpute.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.impute"]] <- myImpute
        message("Run champ.impute() Over!\n")
    }

    tmppd <- myImpute$pd

    ### Applying champ.QC() function.
    if(runQC)
    {
        message("\nRunning champ.QC()...")
        myQCplots <- rep(FALSE,3)
        names(myQCplots) <- c("mdsPlot","densityPlot","dendrogram")
        myQCplots[QCplots] <- TRUE

        champ.QC(beta = tmpbeta,
                 pheno=tmppd$Sample_Group,
                 mdsPlot=myQCplots["mdsPlot"],
                 densityPlot=myQCplots["densityPlot"],
                 dendrogram=myQCplots["dendrogram"],
                 PDFplot=PDFplot,
                 Rplot=Rplot,
                 Feature.sel="None",
                 resultsDir=paste(resultsDir,"/CHAMP_QCimages/",sep=""))

        if(PDFplot)
        {
            message("Plots of champ.QC() has been saved in ",paste(resultsDir,"/CHAMP_QCimages/",sep=""))
        }
        gc()
        message("Run champ.QC() Over!\n")
    }


    ### Applying champ.norm() function.
    if(runnorm)
    {
        message("\nRunning champ.norm()...")
        myNorm <- champ.norm(beta=tmpbeta,
                             rgSet=myLoad$rgSet,
                             mset=myLoad$mset,
                             resultsDir=paste(resultsDir,"/CHAMP_Normalization/",sep=""),
                             method=normalizationmethod,
                             plotBMIQ=PDFplot,
                             arraytype=arraytype,
                             cores=cores)
        if(saveStepresults)
        {
            save(myNorm,file=paste(resultsDir,"/myNorm.rda",sep=""))
            message("champ.norm()'s result \"myNorm\" has been saved in ",resultsDir," as \"myNorm.rda.\"")
            if(normalizationmethod=="BMIQ" & PDFplot==TRUE)
                message("Plots of champ.norm() has been saved in ",paste(resultsDir,"/CHAMP_Normalization/",sep=""))
        }
        gc()
        CHAMP.RESULT[["champ.norm"]] <- myNorm
        message("Run champ.norm() Over!\n")
    }

    ### Applying champ.SVD() function.
    if(runSVD)
    {
        message("\nRunning champ.SVD()...")
        champ.SVD(beta=myNorm,
                  rgSet=myLoad$rgSet,
                  pd=tmppd,
                  RGEffect=RGEffect,
                  PDFplot=PDFplot,
                  Rplot=Rplot,
                  resultsDir=paste(resultsDir,"/CHAMP_SVDimages/",sep=""))

        if(PDFplot)
        {
            message("Plots of champ.SVD() has been saved in ",paste(resultsDir,"/CHAMP_SVDimages/",sep=""))
        }
        gc()
        message("Run champ.SVD() Over!\n")
    }


    ### Applying champ.runCombat() function.
    if(runCombat)
    {
        message("\nRunning champ.runCombat()...")
        myCombat <- champ.runCombat(beta=myNorm,
                                    pd=tmppd,
                                    batchname=batchname,
                                    logitTrans=TRUE)

        if(saveStepresults)
        {
            save(myCombat,file=paste(resultsDir,"/myCombat.rda",sep=""))
            message("champ.runCombat()'s result \"myCombat\" has been saved in ",resultsDir," as \"myCombat.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.runCombat"]] <- myCombat
        message("Run champ.runCombat() Over!\n")
    }

    ### Applying champ.DMP() function.
    if(runDMP)
    {
        message("\nRunning champ.DMP()...")
        myDMP <- champ.DMP(beta = myNorm,
                           pheno = tmppd$Sample_Group,
                           adjPVal = adjPVal,
                           adjust.method = "BH",
                           compare.group = compare.group,
                           arraytype = arraytype)

        if(saveStepresults)
        {
            save(myDMP,file=paste(resultsDir,"/myDMP.rda",sep=""))
            message("champ.DMP()'s result \"myDMP\" has been saved in ",resultsDir," as \"myDMP.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.DMP"]] <- myDMP
        message("Run champ.DMP() Over!\n")
    }

    ### Applying champ.DMR() function.
    if(runDMR)
    {
        message("\nRunning champ.DMR()...")
        
        myDMR <- champ.DMR(beta=myNorm,
                           pheno=tmppd$Sample_Group,
                           arraytype=arraytype,
                           method = DMRmethod,
                           minProbes=7,
                           adjPvalDmr=adjPVal,
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
                           cores=cores,
                           ## following parameters are specifically for probe ProbeLasso method.
                           meanLassoRadius=375,
                           minDmrSep=1000,
                           minDmrSize=50,
                           adjPvalProbe=adjPVal,
                           Rplot=Rplot,
                           PDFplot=PDFplot,
                           resultsDir=paste(resultsDir,"/CHAMP_ProbeLasso/",sep=""),
                           ## following parameters are specifically for DMRcate method.
                           rmSNPCH=T,
                           dist=2,
                           mafcut=0.05)

        if(saveStepresults)
        {
            save(myDMR,file=paste(resultsDir,"/myDMR.rda",sep=""))
            message("champ.DMR()'s result \"myDMR\" has been saved in ",resultsDir," as \"myDMR.rda.\"")
            if(DMRmethod=="ProbeLasso" & PDFplot==TRUE)
                message("Plots of champ.DMR() have been saved in ",paste(resultsDir,"/CHAMP_ProbeLasso/",sep=""))
        }
        gc()
        CHAMP.RESULT[["champ.DMR"]] <- myDMR
        message("Run champ.DMR() Over!\n")
    }

    ### Applying champ.Block() function.
    if(runBlock)
    {
        message("\nRunning champ.Block()...")
        myBlock <- champ.Block(beta=myNorm,
                               pheno=tmppd$Sample_Group,
                               arraytype=arraytype,
                               maxClusterGap=250000,
                               B=500,
                               bpSpan=250000,
                               minNum=10,
                               cores=cores)


        if(saveStepresults)
        {
            save(myBlock,file=paste(resultsDir,"/myBlock.rda",sep=""))
            message("champ.Block()'s result \"myBlock\" has been saved in ",resultsDir," as \"myBlock.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.Block"]] <- myBlock
        message("Run champ.Block() Over!\n")
    }

    ### Applying champ.GSEA() function.
    if(runGSEA)
    {
        message("\nRunning champ.GSEA()...")
        myGSEA <- champ.GSEA(beta=myNorm,
                             DMP=myDMP,
                             DMR=myDMR,
                             CpGlist=NULL,
                             Genelist=NULL,
                             arraytype=arraytype,
                             adjPval=adjPVal)

        if(saveStepresults)
        {
            save(myGSEA,file=paste(resultsDir,"/myGSEA.rda",sep=""))
            message("champ.GSEA()'s result \"myGSEA\" has been saved in ",resultsDir," as \"myGSEA.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.GSEA"]] <- myGSEA
        message("Run champ.GSEA() Over!\n")
    }

    ### Applying champ.EpiMod() function.
    if(runEpiMod)
    {
        message("\nRunning champ.EpiMod()...")

        myEpiMod <- champ.EpiMod(beta=myNorm,
                                 pheno=tmppd$Sample_Group,
                                 nseeds=100,
                                 gamma=0.5,
                                 nMC=1000,
                                 sizeR.v=c(1,100),
                                 minsizeOUT=10,
                                 resultsDir=paste(resultsDir,"/CHAMP_EpiMod/",sep=""),
                                 PDFplot=PDFplot,
                                 arraytype=arraytype)


        if(saveStepresults)
        {
            save(myEpiMod,file=paste(resultsDir,"/myEpiMod.rda",sep=""))
            message("champ.EpiMod()'s result \"myEpiMod\" has been saved in ",resultsDir," as \"myEpiMod.rda.\"")
            if(PDFplot==TRUE)
                message("Plots of champ.EpiMod() have been saved in ",paste(resultsDir,"/CHAMP_EpiMod/",sep=""))
        }
        gc()
        CHAMP.RESULT[["champ.EpiMod"]] <- myEpiMod
        message("Run champ.EpiMod() Over!\n")
    }

    ### Applying champ.CNA() function.
    if(runCNA)
    {
        message("\nRunning champ.CNA()...")

        myCNA <- champ.CNA(intensity=myLoad$intensity,
                           pheno=myLoad$pd$Sample_Group,
                           resultsDir="./CHAMP_CNA",
                           control=TRUE,
                           controlGroup="champCtls",
                           sampleCNA=TRUE,
                           groupFreqPlots=TRUE,
                           Rplot=Rplot,
                           PDFplot=PDFplot,
                           freqThreshold=0.3,
                           arraytype=arraytype)


        if(saveStepresults)
        {
            save(myCNA,file=paste(resultsDir,"/myCNA.rda",sep=""))
            message("champ.CNA()'s result \"myCNA\" has been saved in ",resultsDir," as \"myCNA.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.CNA"]] <- myCNA
        message("Run champ.CNA() Over!\n")
    }

    if(runRefBase)
    {
        message("\nRunning champ.refbase()...")
        myRefbase <- champ.refbase(beta=myNorm,
                                   arraytype=arraytype)

        if(saveStepresults)
        {
            save(myRefbase,file=paste(resultsDir,"/myRefbase.rda",sep=""))
            message("champ.refbase()'s result \"myRefbase\" has been saved in ",resultsDir," as \"myRefbase.rda.\"")
        }
        gc()
        CHAMP.RESULT[["champ.refbase"]] <- myRefbase
        message("Run champ.refbase() Over!\n")
    }

    message("[<<<< ChAMP.PROCESS END >>>>]")
    message("[===========================]")
    return(CHAMP.RESULT)
}
