if(getRversion() >= "3.1.0") utils::globalVariables(c("sampleNames<-","EPIC.manifest.hg19","EPIC.manifest.pop.hg19","hm450.manifest.hg19","hm450.manifest.pop.hg19","multi.hit","probe.features","probe.features.epic"))


champ.filter <- function(beta=myImport$beta,
                         M=NULL,
                         pd=myImport$pd,
                         intensity=NULL,
                         Meth=NULL,
                         UnMeth=NULL,
                         detP=NULL,
                         beadcount=NULL,
                         autoimpute=TRUE,
                         filterDetP=TRUE,
                         ProbeCutoff=0,
                         SampleCutoff=0.1,
                         detPcut=0.01,
                         filterBeads=TRUE,
                         beadCutoff=0.05,
                         filterNoCG = TRUE,
                         filterSNPs = TRUE, 
                         population = NULL,
                         filterMultiHit = TRUE,
                         filterXY = TRUE, 
                         fixOutlier = TRUE,
                         arraytype = "450K")
{
  message("[===========================]")
  message("[<<<< ChAMP.FILTER START >>>>>]")
  message("-----------------------------")
  
  message("\nIn New version ChAMP, champ.filter() function has been set to do filtering on the result of champ.import(). You can use champ.import() + champ.filter() to do Data Loading, or set \"method\" parameter in champ.load() as \"ChAMP\" to get the same effect.")
  
  message("\nThis function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix and Bead Count, you also need to input a detected P matrix and Bead Count information.")
  
  message("\nNote that if you want to filter more data matrix, say beta, M, intensity... please make sure they have exactly the same rownames and colnames.")
  
  message("\n\n[ Section 1:  Check Input Start ]")
  
  Objects <- list("beta"=beta,
                  "M"=M,
                  "intensity"=intensity,
                  "Meth"=Meth,
                  "UnMeth"=UnMeth)
  
  Objects <- Objects[which(lapply(Objects,FUN=is.null)==FALSE)]
  if(length(Objects)==0) stop("  At least one Data Set needed.")
  message("  You have inputed ",paste(names(Objects),collapse=",")," for Analysis.")
  if(length(unique(lapply(Objects,FUN=rownames))) != 1) stop("  !!!  You objects have different rownames. Please make sure your Matrix are in accord on Rows.")
  if(length(unique(lapply(Objects,FUN=colnames))) != 1) stop("  !!!  You objects have different colnames. Please make sure your Matrix are in accord on Rows.")
  
  Accessory <- list("detP"=detP,
                    "beadcount"=beadcount)
  
  FilterOption <- list("filterDetP"=filterDetP,
                       "autoimpute"=autoimpute,
                       "filterBeads"=filterBeads,
                       "filterMultiHit"=filterMultiHit,
                       "filterSNPs"=filterSNPs,
                       "filterNoCG"=filterNoCG,
                       "filterXY"=filterXY)
  
  ### Checking pd file
  if(!is.null(pd)) {
    message("\n  pd file provided, checking if it's in accord with Data Matrix...")
    if(nrow(pd) == ncol(Objects[[1]])) {
      if("Sample_Name" %in% names(pd)){
        if(identical(as.character(pd$Sample_Name),colnames(Objects[[1]]))) 
          message("    pd file check success.")
        else
          stop("    Your pd file's Sample_Name is different from your Data Matrix colnames.")
      } else {
        message("    !!! Your pd file does not have Sample_Name column, we can not check your Sample_Name, please make sure the pd file is correct.")
      }
    } else {
      stop("    pd file and Data matrix have different dimentions.")
    }
  }
  
  ### Checking Detect P value
  if(FilterOption$filterDetP == TRUE)
  {
    message("\n  Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...")
    if(!is.null(Accessory$detP)) {
      if(identical(rownames(Accessory$detP),rownames(Objects[[1]])) & identical(colnames(Accessory$detP),colnames(Objects[[1]])))
      {
        message("    detP check success.")
      }
      else {
        message("    !!! Your detP does not have the EXACT same rowname and colname as Data Matrix.")
        message("    !!! filterDetP can not be done. filterDetP is reset FALSE now.")
        FilterOption$filterDetP <- FALSE
        Accessory$detP <- NULL
      }  
    } else {
      message("    !!! Parameter detP is not found, filterDetP is reset FALSE now.")
      FilterOption$filterDetP <- FALSE
      Accessory$detP <- NULL
    }
  }
  
  ### Checking Beadcount value
  if(FilterOption$filterBeads == TRUE)
  {
    message("\n  Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...")
    if(!is.null(Accessory$beadcount)) {
      if(identical(rownames(Accessory$beadcount),rownames(Objects[[1]])) & identical(colnames(Accessory$beadcount),colnames(Objects[[1]])))
      {
        message("    beadcount check success.")
      }
      else {
        message("    !!! Your beadcount does not have the EXACT same rowname and colname as Data Matrix.")
        message("    !!! filterBeads can not be done. filterBeads is reset FALSE now.")
        FilterOption$filterBeads <- FALSE
        Accessory$beadcount <- NULL
      }  
    } else {
      message("    !!! Parameter beadcount is not found, filterBeads is reset FALSE now.")
      FilterOption$filterBeads <- FALSE
      Accessory$beadcount <- NULL
    }
  }
  
  if(FilterOption$autoimpute == TRUE)
  {
    message("\n  parameter autoimpute is TRUE. Checking if the conditions are fulfilled...")
    if("beta" %in% names(Objects) | "M" %in% names(Objects)){
      if(ProbeCutoff > 0 & !is.null(detP)) {
        message("    autoimpute check success.")
      } else {
        message("    !!! ProbeCutoff is 0, which means you have no needs to do imputation. autoimpute has been reset FALSE.")
        FilterOption$autoimpute <- FALSE
      }
    } else {
      message("    !!! beta matrix or M matrix are required for impute. autoimpute has been reset FALSE.")
      FilterOption$autoimpute <- FALSE
    }
  }
  
  if(sum(FilterOption==TRUE)>0) {
    message("\n  Checking Finished :" , paste(names(which(FilterOption==TRUE)),collapse=",")," would be done on ", paste(names(Objects),collapse=","),".")
  } else {
    stop("  No Filtering would be done. Please check your input parameters.")
  }
  if(length(Accessory) > 0) {
    message("  You also provided :", paste(names(Accessory),collapse=",")," .")
  }
  message("[ Section 1: Check Input Done ]")
  
  ### Start Filtering Here
  message("\n\n[ Section 2: Filtering Start >>")
  
  if(FilterOption$filterDetP == TRUE)
  {
    message("\n  Filtering Detect P value Start")
    message("    The fraction of failed positions per sample")
    message("    You may need to delete samples with high proportion of failed probes:\n")
    
    numfail <- matrix(colMeans(Accessory$detP >= detPcut))
    rownames(numfail) <- colnames(Accessory$detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < SampleCutoff)
    
    if(any(numfail >= SampleCutoff))
    {
      message("\n    The detSamplecut parameter is : ",SampleCutoff)
      message("    Samples : ", paste(rownames(numfail)[which(numfail >= SampleCutoff)],collapse=",")," will be deleted.")
      message("    There are ",length(RemainSample)," samples remained for analysis.")
      Objects <- lapply(Objects,function(x) x[,RemainSample])
      Accessory <- lapply(Accessory,function(x) x[,RemainSample])
    }
    
    RemainProbe <- rowSums(Accessory$detP > detPcut) <= ProbeCutoff * length(RemainSample)
    if(ProbeCutoff==0)
    {
      message("\n    Filtering probes with a detection p-value above ",detPcut,".")
      message("    Removing ",sum(RemainProbe==FALSE)," probes.")
      message("    If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples")
    } else {
      message("\n    Filtering probes with a detection p-value above ",detPcut," in at least", ProbeCutoff*100,"% Samples.")
      message("    Removing ",sum(RemainProbe==FALSE)," probes.")
      message("    If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples")
    }
    
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
    if(sum(Accessory$detP > detPcut) > 0) message("    There are still ",sum(Accessory$detP > detPcut), " failed probes exist in your data set, imputation is recommended.")
  }
  
  
  if(FilterOption$autoimpute == TRUE)
  {
    message("\n  Autoimpute Start") 
    if(sum(Accessory$detP > detPcut) == 0) {
      message("    No NAs (failed probes) exist in your data set any more, you don't need to do any imputation.")
    } else
    {
      message("    There are ",sum(Accessory$detP > detPcut), " NAs (failed probes) exists in your data set.")
      message("    Impute.knn will be conducted for remain NAs. (NOT suitable for small data sets)\n")
      
      ### Using sink function to remove messages
      zz <- file("ImputeMessage.Rout", open="wt")
      sink(zz)
      sink(zz, type="message")
      
      if("beta" %in% names(Objects))
      {
        message("    Doing imputation on beta matrix.")
        Objects$beta[Accessory$detP > detPcut] <- NA
        Objects$beta <- suppressMessages(impute.knn(Objects$beta, k=5)$data)
      }
      if("M" %in% names(Objects))
      {
        message("    Doing imputation on M matrix.")
        Objects$M[Accessory$detP > detPcut] <- NA
        Objects$M <- suppressMessages(impute.knn(Objects$M,k=5)$data)
      }
      sink(type="message")
      sink()
    }
  }
  
  if(FilterOption$filterBeads == TRUE)
  {
    message("\n  Filtering BeadCount Start")
    RemainProbe <- rowSums(is.na(Accessory$beadcount)) < beadCutoff*(ncol(Accessory$beadcount))
    message("    Filtering probes with a beadcount <3 in at least ",beadCutoff*100,"% of samples.")
    message("    Removing ",sum(RemainProbe == FALSE)," probes")
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
  }
  
  if(FilterOption$filterNoCG == TRUE)
  {
    message("\n  Filtering NoCG Start")
    RemainProbe <- substr(rownames(Objects[[1]]),1,2) == "cg"
    message("    Only Keep CpGs, removing ", sum(RemainProbe == FALSE) ," probes from the analysis.")
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
  }
  
  if(FilterOption$filterSNPs== TRUE)
  {
    message("\n  Filtering SNPs Start")
    if(arraytype %in% c("EPIC", "EPICv2")) {
      message("\n    !!! Important, since version 2.29.1, ChAMP set default `EPIC` arraytype as EPIC version 2. ",
              "\n        You can set 'EPIC' or 'EPICv2' to use version 2 EPIC annotation",
              "\n        If you want to use the old version (v1), please specify arraytype parameter as `EPICv1`. ",
              "\n        For 450K array, still use `450K`\n")
      data("EPICv2_Mask")
    } else if(arraytype == "EPICv1") {
      data("EPICv1_Mask")
    } else { 
      data("Hm450K_Mask")
    }
    
    if(is.null(population)) {
      message("    Using general mask options")
      maskname <- rownames(Mask)[Mask$MASK_general]
    } else if (paste0("MASK_general_", population) %in% colnames(Mask)) {
      message("    Using ", population, " specific SNP option.")
      maskname <- rownames(Mask)[Mask[, paste0("MASK_general_", population)]]
    } else {
      message("    Seems the population does not exist in this array annotation, using general options.")
      maskname <- rownames(Mask)[Mask$MASK_general]
    }
    
    RemainProbe <- !rownames(Objects[[1]]) %in% maskname
    message("    Removing ", sum(RemainProbe == FALSE) ," probes from the analysis.")
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
  }
  
  
  
  if(FilterOption$filterMultiHit == TRUE)
  {
    message("\n  Filtering MultiHit Start")
    data(multi.hit)
    RemainProbe <- !rownames(Objects[[1]]) %in% multi.hit$TargetID
    message("    Filtering probes that align to multiple locations as identified in Nordlund et al")
    message("    Removing ", sum(RemainProbe == FALSE) ," probes from the analysis.")
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
  }
  
  
  if(FilterOption$filterXY == TRUE)
  {
    message("\n  Filtering XY Start")
    if(arraytype %in% c("EPIC", "EPICv2")) {
      data("probe.features.epicv2")
    } else if(arraytype == "EPICv1") {
      data("probe.features.epicv1")
    } else if(arraytype == "450K") { 
      data("probe.features")
    } else (
      stop("arraytype must be `EPICv2`, `EPICv1`, `450K`")
    )
    
    RemainProbe <- rownames(Objects[[1]]) %in% (rownames(probe.features)[!probe.features$CHR %in% c("X","Y", "chrX", "chrY")])
    message("    Filtering probes located on X,Y chromosome, removing ", sum(RemainProbe == FALSE) ," probes from the analysis.")
    Objects <- lapply(Objects,function(x) x[RemainProbe,])
    Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
  }
  
  if(!is.null(pd))
  {
    message("\n  Updating PD file")
    if(FilterOption$filterDetP==FALSE)
    {
      message("    filterDetP parameter is FALSE, so no Sample Would be removed.")
      RemainSample <- 1:nrow(pd)
    }
    pd <- pd[RemainSample,]
    Objects <- append(Objects,list(pd=pd))
  }
  
  if(fixOutlier & "beta" %in% names(Objects))
  {
    message("\n  Fixing Outliers Start")
    message("    Replacing all value smaller/equal to 0 with smallest positive value.")
    Objects$beta[Objects$beta <= 0] <- min(Objects$beta[which(Objects$beta > 0)])
    message("    Replacing all value greater/equal to 1 with largest value below 1..")
    Objects$beta[Objects$beta >= 1] <- max(Objects$beta[which(Objects$beta < 1)])
  }
  message("[ Section 2: Filtering Done ]")
  
  message("\n All filterings are Done, now you have ", nrow(Objects[[1]]), " probes and ",ncol(Objects[[1]]), " samples.\n")
  
  message("[<<<<< ChAMP.FILTER END >>>>>>]")
  message("[===========================]")
  message("[You may want to process champ.QC() next.]\n")
  return(Objects)
}

