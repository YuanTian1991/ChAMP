plotThreshold <-
function(myResults,lasso.radii,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
{
data(probe.features)

	no.feature <- length(summary(as.factor(probe.features$feature.1)))
	no.relation <- length(summary(as.factor(probe.features$RELATION_TO_UCSC_CPG_ISLAND))) - 2 # merges North- & South- shore & shelf relations
	no.feat.rel <- length(summary(as.factor(probe.features$feat.rel)))
    
	spc.thrd.qt <- matrix(nrow=99, ncol=no.feat.rel)
	for (i in 1:99)
	{
		spc.thrd.qt[i,] <- as.numeric(tapply(myResults$nrst.probe, myResults$feat.rel, function(x) quantile(x,i/100, na.rm=T))) # calculates nearest probe split by feature
	}
    
    imageName <- paste(resultsDir,"myLassos.pdf",sep="/")
    pdf(imageName,width=9,height=9)
	par(bty="n", adj=0.5, mar=c(8,4,7,2)+0.1)
	boxplot(myResults$nrst.probe~myResults$feat.rel, outline=F, ylab="nearest neighbour [bp]", col=rep(rainbow(no.feature),each=no.relation), las=2, lwd=0.2)
    axis(side=3, at=c(1:no.feat.rel), labels=paste("N =", format(as.numeric(tapply(myResults$nrst, myResults$feat.rel, length)), big.mark=",", scientific=F), sep=""), tick=FALSE, cex.axis=0.8, las=2)
    title(main="a", adj=0, cex.main=2)
	plot(log10(spc.thrd.qt[,1]), ylim=log10(range(spc.thrd.qt)), xlab="quantile", ylab="lasso radius [bp]", lty=1, lwd=2, col=rainbow(no.feature)[1], type="l", xaxt="n", yaxt="n")
    axis(side=1, at=seq(0,100,10), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,5,1), labels=c(0,10,100,1000,10000,100000), las=2)
    for (i in 2:no.feat.rel)
    {
        lines(log10(spc.thrd.qt[,i]), lty=rep(c(1,2,3,6),no.feature)[i], lwd=2, col=rep(rainbow(no.feature),each=no.relation)[i])
    }
    legend(0,5,legend=c(unique(sapply(strsplit(names(summary(as.factor(probe.features$feat.rel))), "_"), "[[",1)),unique(sapply(strsplit(names(summary(as.factor(probe.features$feat.rel))), "_"), "[[",2))), col=c(rainbow(7),rep("black",4)), lty=c(rep(1,7),c(1,2,3,6)), lwd=2, bty="n", cex=0.8)
    if(lassoStyle == "max")
    {
        segments(c(0,value.lasso.quantile*100), rep(log10(lassoRadius), 2), rep(value.lasso.quantile*100,2), c(log10(lassoRadius),0), lty=2)
    }else
    {
        segments(c(0,value.lasso.quantile*100), c(log10(lassoRadius), log10(max(lasso.radii))), rep(value.lasso.quantile*100,2), c(log10(lassoRadius),0), lty=2)
    }
    title(main="b", adj=0, cex.main=2)
	plot(c(1,28), y=c(range(0.3*sqrt(lasso.radii))[1]*0.8, range(0.3*sqrt(lasso.radii))[2]*1.2), type="n", xaxt="n", xlab="", yaxt="n", ylab="lasso radius [bp]", main=paste("lasso quantile = ", round(value.lasso.quantile,2), sep=""), bty="n")
    segments(1:28, rep(0,28), 1:28, 0.3*sqrt(lasso.radii), lty=3, col="grey")
    points(1:28, 0.3*sqrt(lasso.radii), pch=16, cex=0.3*sqrt(lasso.radii), col=rep(rainbow(7,alpha=0.5), each=4))
    text(1:28, 0.3*sqrt(lasso.radii), lasso.radii, pos=3, cex=0.8)
    axis(1, at=1:28,rownames(lasso.radii), las=2, cex.axis=0.8, tick=F)
    axis(2, at=c(0,max(0.3*sqrt(lasso.radii))), labels=F)
    title(main="c", adj=0, cex.main=2)
	dev.off()
}
