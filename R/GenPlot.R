GenPlot <-
function(thdens.o,estdens.o,evalues.v){


minx <- min(min(thdens.o$lambda),min(evalues.v));
maxx <- max(max(thdens.o$lambda),max(evalues.v));
miny <- min(min(thdens.o$dens),min(estdens.o$y));
maxy <- max(max(thdens.o$dens),max(estdens.o$y));


#plot(thdens.o$lambda,thdens.o$dens,xlim=c(minx,maxx),ylim=c(miny,maxy),type="b",col="green");
#i <- min(which(estdens.o$x > min(evalues.v)));
#f <- max(which(estdens.o$x < max(evalues.v)));
#points(x=estdens.o$x[i:f],y=estdens.o$y[i:f],type="b",col="red");

}
