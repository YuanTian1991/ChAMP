if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","hprdAsigH.m"))

champ.EpiMod <- function(beta=myNorm,
                         pheno=myLoad$pd$Sample_Group,
                         nseeds=100,
                         gamma=0.5,
                         nMC=1000,
                         sizeR.v=c(1,100),
                         minsizeOUT=10,
                         resultsDir="./CHAMP_EpiMod/",
                         PDFplot=TRUE,
                         arraytype="450K")
{
    message("[===========================]")
    message("[<<< ChAMP.EpiMod START >>>>]")
    message("-----------------------------")

    ### Prepare Checking ###
    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.EpiMod Results will be saved in ",resultsDir)

    data(hprdAsigH)
    message("<< Load PPI network hprdAsigH >>")

    if(arraytype=="EPIC")
        statM.o <- GenStatM(beta,pheno,arraytype)
    else
        statM.o <- GenStatM(beta,pheno,"450k")
    message("<< Generate statM.o >>")
        
    intEpi.o=DoIntEpi450k(statM.o,hprdAsigH.m,c=1)

    message("<< Calculate EpiMod.o >>")
    EpiMod.o=DoEpiMod(intEpi.o,
                      nseeds=nseeds,
                      gamma=gamma,
                      nMC=nMC,
                      sizeR.v=sizeR.v,
                      minsizeOUT=minsizeOUT,
                      writeOUT=TRUE,
                      ew.v=NULL);

    if(PDFplot)
    {
        message("<< Draw All top significant module plot in PDF >>")
        tmpdir <- getwd()
        setwd(resultsDir)
        for(i in names(EpiMod.o$topmod)) FemModShow(EpiMod.o$topmod[[i]],name=i,EpiMod.o,mode="Epi")
        setwd(tmpdir)
    }

    message("[<<<< ChAMP.EpiMod END >>>>>]")
    message("[===========================]")
    return(EpiMod.o)
}
