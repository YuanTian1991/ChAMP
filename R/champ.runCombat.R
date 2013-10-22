champ.runCombat <-
function(beta.c=myNorm$beta,pd=myLoad$pd,logitTrans=TRUE)
{
 	
    message("Preparing files for ComBat")
    
    batch = pd$Slide
 	
    if(length(unique(batch))<2)
	{
		message("You must have at least two different slides for Combat to correct for batch effects")
		return(list(beta=beta.c,pd=pd))
        
	}else if(any(count(batch)[2]<2)){
        
        message("You must have at least two samples on each slide for Combat to correct for batch effects")
        return(list(beta=beta.c,pd=pd))
        
    }else{
        
        if(min(beta.c)==0)
        {
            message("Zeros in your dataset have been replaced with 0.000001")
            beta.c[beta.c==0]<-0.000001
        }
        if(logitTrans)
        {
            message("Your data is being logit transformed before batch correction")
            mod = model.matrix(~as.factor(Sample_Group), data=pd)
            log=logit2(beta.c)
            message("Beginning batch correction")
            combat=champ.ComBat(dat=log,batch=batch,mod=mod,par.prior=TRUE)
            message("Your data is being inverse logit transformed")
            if(!is.null(combat)){combat=ilogit2(combat)}
            
        }else{
            message("Beginning batch correction without logit transformation.")
            mod = model.matrix(~as.factor(Sample_Group), data=pd)
            combat=champ.ComBat(dat=beta.c,batch=batch,mod=mod,par.prior=TRUE)
            
        }
        if(is.null(combat))
        {
            return(list(beta=beta.c))
        }else{
            return(list(beta=combat))
        }
        
    }
    
    
}
