if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm","myLoad","ControlProbes"))

champ.SVD <- function(beta=myNorm,
                      rgSet=NULL,
                      pd=myLoad$pd,
                      RGEffect=FALSE,
                      PDFplot=TRUE,
                      Rplot=TRUE,
                      resultsDir="./CHAMP_SVDimages/")
{
    message("[===========================]")
    message("[<<<<< ChAMP.SVD START >>>>>]")
    message("-----------------------------")

    ### Some inner functions
   GenPlot <- function(thdens.o,estdens.o,evalues.v){
       minx <- min(min(thdens.o$lambda),min(evalues.v));
       maxx <- max(max(thdens.o$lambda),max(evalues.v));
       miny <- min(min(thdens.o$dens),min(estdens.o$y));
       maxy <- max(max(thdens.o$dens),max(estdens.o$y));
   }
 
    EstDimRMTv2 <- function(data.m)
    {    
        ### standardise matrix
        M <- data.m;
        for(c in 1:ncol(M)) M[,c] <- (data.m[,c]-mean(data.m[,c]))/sqrt(var(data.m[,c]));
        sigma2 <- var(as.vector(M));
        Q <- nrow(data.m)/ncol(data.m);
        thdens.o <- thdens(Q,sigma2,ncol(data.m));
        C <- 1/nrow(M) * t(M) %*% M;

        eigen.o <- eigen(C,symmetric=TRUE);
        estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);

        GenPlot(thdens.o,estdens.o,eigen.o$values);
        intdim <- length(which(eigen.o$values > thdens.o$max));
        return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o));
    }
    thdens <- function(Q,sigma2,ns)
    {
        lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
        lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
        delta <- lambdaMAX - lambdaMIN;#  print(delta);
        roundN <- 3;
        step <- round(delta/ns,roundN);
        while(step==0){
            roundN <- roundN+1;
            step <- round(delta/ns,roundN);
        }
        lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
        dens.v <- vector();
        ii <- 1;
        for(i in lambda.v){
            dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
            ii <- ii+1;
        }
        return(list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v));
    }


    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.SVD Results will be saved in ",resultsDir," .\n")

    if(length(which(is.na(beta)))>0) message(length(which(is.na(beta)))," NA are detected in your beta Data Set, which may cause fail or uncorrect of SVD analysis. You may want to impute NA with champ.impute() function first.")

    message("[SVD analysis will be proceed with ",dim(beta)[1], " probes and ",dim(beta)[2], " samples.]\n")

    message("\n[ champ.SVD() will only check the dimensions between data and pd, instead if checking if Sample_Names are correctly matched (because some user may have no Sample_Names in their pd file),thus please make sure your pd file is in accord with your data sets (beta) and (rgSet).]\n")


	################ Customise Phenotype Data ########################

    if(is.null(pd) | class(pd)=="list") stop("pd parameter in Data Frame or Matrix is necessary And must contain at least tow factors. If your pd is a list, please change its Format.")
    if(class(pd)=="matrix") pd <- as.data.frame(pd)

    PhenoTypes.lv_tmp <- pd[,!colnames(pd) %in% c("Sample_Name","Project","filenames","Basename") & apply(pd,2,function(x) length(unique(x)))!=1]
    PhenoTypes.lv <- PhenoTypes.lv_tmp

    if(!is.null(rownames(pd))) rownames(PhenoTypes.lv) <- rownames(pd)

    if(ncol(PhenoTypes.lv)>=2)
    {
        message("<< Following Factors in your pd(sample_sheet.csv) will be analysised: >>")
        sapply(colnames(PhenoTypes.lv_tmp),function(x) message("<",x,">(",class(PhenoTypes.lv[[x]]),"):",paste(unique(PhenoTypes.lv_tmp[,x]),collapse=", ")))
        message("[champ.SVD have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv), if you don't want to analysis some of them, please remove them manually from your pd variable then retry champ.SVD().]")
    }else
    {
        stop("You don't have even one factor with at least two value to be analysis. Maybe your factors contains only one value, no variation at all...")
    }

    if(ncol(pd) > ncol(PhenoTypes.lv))
    {
        message("\n<< Following Factors in your pd(sample_sheet.csv) will not be analysis: >>")
        sapply(setdiff(colnames(pd),colnames(PhenoTypes.lv)),function(x) message("<",x,">"))
        message("[Factors are ignored because they only indicate Name or Project, or they contain ONLY ONE value across all Samples.]")
    }
    
    #### PhenoTypes.lv prepare ready.
    if(RGEffect==TRUE & is.null(rgSet)) message("If you want to check Effect of Control Probes, you MUST provide rgSet parameter. Now champ.SVD can only analysis factors in pd.")
	if(!is.null(rgSet) & RGEffect)
	{
        if(rgSet@annotation[1]=="IlluminaHumanMethylation450k") data(ControlProbes450K) else data(ControlProbesEPIC)
        dataC2.m <- as.data.frame(log2(apply(ControlProbes,1,function(x) if(x[3]=="Grn") getGreen(rgSet)[x[2],] else getRed(rgSet)[x[2],])))
        PhenoTypes.lv <- cbind(PhenoTypes.lv,dataC2.m)
        message("\n<< Following rgSet information have been added to PhenoTypes.lv. >>")
        sapply(colnames(dataC2.m),function(x) message("<",x,">"))
        #### RG Information extract!
    }	

    if(nrow(PhenoTypes.lv)==ncol(beta)) message("\n<< PhenoTypes.lv generated successfully. >>") else stop("Dimension of your pd file (and rgSet information) is not equal to your beta matrix.")

	######################## Do SVD #############################

    tmp.m <- beta-rowMeans(beta)
	rmt.o <- EstDimRMTv2(tmp.m);
	svd.o <- svd(tmp.m);
    if(rmt.o$dim > 20) topPCA <- 20  else topPCA <- rmt.o$dim
        
    svdPV.m <- matrix(nrow=topPCA,ncol=ncol(PhenoTypes.lv));
    colnames(svdPV.m) <- colnames(PhenoTypes.lv);

	for(c in 1:topPCA)
        for(f in 1:ncol(PhenoTypes.lv))
            if(class(PhenoTypes.lv[,f])!="numeric")
                svdPV.m[c,f] <- kruskal.test(svd.o$v[,c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
            else
                svdPV.m[c,f] <- summary(lm(svd.o$v[,c] ~ PhenoTypes.lv[[f]]))$coeff[2,4];

    message("<< Calculate SVD matrix successfully. >>")

	######################## Plot SVD Image #############################
	myPalette <- c("darkred","red","orange","pink","white");
	breaks.v <- c(-200,-10,-5,-2,log10(0.05),0);

    if(Rplot)
    {
        par(mar=c(5,15,2,1));
        image(x=1:nrow(svdPV.m), y=1:ncol(svdPV.m), z=log10(svdPV.m), col=myPalette, breaks=breaks.v, xlab="", ylab="", axes=FALSE, main= "Singular Value Decomposition Analysis (SVD)");
        axis(1,at=1:nrow(svdPV.m),labels=paste("PC-",1:nrow(svdPV.m),sep=""),las=2);
        suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
        legend(x=-(topPCA/2.5),y=3,legend=c(expression("p < 1x"~10^{-10}),expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), fill=c("darkred","red","orange","pink","white"),par('usr')[2], par('usr')[4], xpd=NA);    
    }
    if(PDFplot)
    {
        pdf(paste(resultsDir,"SVDsummary.pdf",sep=""),width=8,height=8);
        par(mar=c(5,15,2,1),xpd=TRUE);
        image(x=1:nrow(svdPV.m), y=1:ncol(svdPV.m), z=log10(svdPV.m), col=myPalette, breaks=breaks.v, xlab="", ylab="", axes=FALSE, main= "Singular Value Decomposition Analysis (SVD)");
        axis(1,at=1:nrow(svdPV.m),labels=paste("PC-",1:nrow(svdPV.m),sep=""),las=2);
        suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
        legend(x=-(topPCA/2.5),y=3,legend=c(expression("p < 1x"~10^{-10}),expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), fill=c("darkred","red","orange","pink","white"),par('usr')[2], par('usr')[4], xpd=NA,);    
        dev.off();
    }

    message("<< Plot SVD matrix successfully. >>")

    message("[<<<<<< ChAMP.SVD END >>>>>>]")
    message("[===========================]")
    message("[If the batch effect is not significant, you may want to process champ.DMP() or champ.DMR() or champ.BlockFinder() next, otherwise, you may want to run champ.runCombat() to eliminat batch effect, then rerun champ.SVD() to check corrected result.]\n")
}
