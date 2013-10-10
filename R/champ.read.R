champ.read <-
function(betaFile="beta.txt",sampleSheet="sampleSheet.txt",resultsDir)
{
	beta <-read.table(betaFile,row.names =T, sep = "\t")
	pd <- read.table(sampleSheet) #need to check
	if(file.exists(resultsDir))
	{
		message("The directory ",resultsDir," already exists. To avoid overwriting please rename in the arguments or delete this folder.")
		
	}else{ 
    
		dir.create(resultsDir)	
	}	
	return(list(beta=beta,pd=pd))
}
