plotThreshold <-
function(myResults,lasso.radius,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
{
    data(probe.features)
    
    no.feature <- length(summary(as.factor(probe.features$feature)))
    names.feature <- names(summary(as.factor(probe.features$feature)))
    no.cgi <- length(summary(as.factor(probe.features$cgi)))
    names.cgi <- names(summary(as.factor(probe.features$cgi)))
    no.feat.cgi <- length(summary(as.factor(probe.features$feat.cgi)))
    names.feat.cgi <- names(summary(as.factor(probe.features$feat.cgi)))
    
	spc.thrd.qt <- matrix(nrow=99, ncol=no.feat.cgi)
	for (i in 1:99)
	{
		spc.thrd.qt[i,] <- as.numeric(tapply(myResults$nrst.probe, myResults$feat.cgi, function(x) quantile(x,i/100, na.rm=T))) # calculates nearest probe split by feature
	}
    
    sfo <- c(6, 7, 3, 1, 4, 2, 5)
	scgio <- c(1, 4, 3, 2)
	sfcgio <- rep((sfo - 1) *4, each = 4) + rep(scgio, 7)
    
	tmp1 <- split(myResults$nrst.probe, myResults$feat.cgi)
    
    imageName <- paste(getwd(),"myLassos.pdf",sep="/")
    pdf(imageName,width=9,height=9)
    ## plot 1 of 3
    par(bty="n", adj=0.5, mar=c(8,4,7,2)+0.1)
    plot(c(0.5, 28.5), c(1, max(unlist(lapply(tmp1, function(x) quantile(x, .95))))*1.1), type = "n", xaxt = "n", xlab = "", ylab = "nearest neighbour [bp]", log = "y")
    par(cex = 1)
    for (i in 1:28)
    {
        boxplot(tmp1[[sfcgio[i]]], col = rep(rainbow(7), each = 4)[sfcgio[i]], at = i, add = T, outline = F, yaxt = "n", cex.axis  = 1)
    }
    axis(1, at = 1:28, labels = rep(names.cgi[scgio], 7), las = 2)
    axis(3, at = 1:28, labels=paste("N =", format(unlist(lapply(tmp1, length))[sfcgio], big.mark=",", scientific=F), sep=""), tick=FALSE, cex.axis=0.8, las=2)
    8.
    segments(seq(4.5, 28, 4), rep(0.05, 7), seq(4.5, 28, 4), rep(max(unlist(lapply(tmp1, function(x) quantile(x, .95))))*1.1, 7), lty = 2, col = "lightgrey")
    
    # plot 2 of 3
    par(bty="n", adj=0.5, mar=c(8,4,4,2)+0.1)
    plot(c(0, 100), y = range(spc.thrd.qt), ylim=range(spc.thrd.qt), xlab="quantile", ylab="lasso radius [bp]", xaxt="n", log = "y", type = "n")
    axis(side=1, at=seq(0,100,10), labels=seq(0,1,0.1))
    for (i in 1:28)
    {
        lines(spc.thrd.qt[,sfcgio[i]], lty=rep(c(1,2,3,6),no.feature)[sfcgio[i]], lwd=2, col=rep(rainbow(no.feature),each=no.cgi)[sfcgio[i]])
    }
    legend("topleft",legend=names.feature[sfo], col=rainbow(7)[sfo], lty=rep(1,7), lwd=2, bty="n", cex=0.8)
    legend(20, 1e5, legend=names.cgi[scgio], col=rep("black",4), lty=c(1,2,3,6)[scgio], lwd=2, bty="n", cex=0.8)
    if(lassoStyle == "max")
    {
        segments(c(-5,value.lasso.quantile*100), rep(lassoRadius, 2), rep(value.lasso.quantile*100,2), c(lassoRadius,1), lty=2)
    }else
    {
        segments(c(-5,value.lasso.quantile*100), c(lassoRadius, max(lasso.radius)), rep(value.lasso.quantile*100,2), c(lassoRadius,0), lty=2)
    }
    
    # plot 3 of 3
    plot(c(1,28), y=c(range(0.3*sqrt(lasso.radius))[1]*0.8, range(0.3*sqrt(lasso.radius))[2]*1.2), type="n", xaxt="n", xlab="", yaxt="n", ylab="lasso radius [bp]", main=paste("lasso quantile = ", round(value.lasso.quantile,2), sep=""), bty="n")
    segments(1:28, rep(0,28), 1:28, 0.3*sqrt(lasso.radius[sfcgio]), lty=3, col="grey")
    points(1:28, 0.3*sqrt(lasso.radius[sfcgio]), pch=16, cex=0.3*sqrt(lasso.radius[sfcgio]), col=rep(rainbow(7,alpha=0.5)[sfo], each=4))
    text(1:28, 0.3*sqrt(lasso.radius[sfcgio]), lasso.radius[sfcgio], pos=3, cex=0.8)
    axis(1, at = 1:28, labels = rep(names.cgi[scgio], 7), las = 2)
    par(xpd = T)
    segments(seq(1, 28, 4), rep(-3.5, 7), seq(4, 28, 4), rep(-3.5, 7))
    mtext(text = names.feature[sfo], side = 1, at = seq(2.5, 28, 4), line = 5.5, las = 1, cex.axis = 1)
    axis(2, at=c(0,max(0.3*sqrt(lasso.radius))), labels=F)
    dev.off()
    rm(lasso.radius)
    #return(myDf)
    #gc()
}