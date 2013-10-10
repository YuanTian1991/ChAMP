champ.bedfile <-
function(info)
{
	bedfile=data.frame(matrix(NA,ncol=4,nrow=nrow(info)))

	#check for MAPINFO column
	colnames(bedfile)[1]="chr"
	colnames(bedfile)[2]="start"
	colnames(bedfile)[3]="end"
	colnames(bedfile)[4]="probe"
	bedfile$probe=info$probeID
	bedfile$chr<-paste("chr",info$CHR,sep="")
	bedfile$start<-info$MAPINFO-1
	bedfile$end<-info$MAPINFO

	bedfile<-bedfile[order(bedfile$chr,bedfile$start),]
	return(bedfile)
}
