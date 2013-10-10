champ.colname <-
function(beta,pd)
{
	c=colnames(beta)
	for( i in 1:length(c)){c[i]= pd$Sample_Name[i]}
	colnames(beta)=c
	return(beta)
}
