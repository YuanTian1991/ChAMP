EstDimRMTv2 <-
function(data.m){

 ### standardise matrix
 M <- data.m;
 for(c in 1:ncol(M)){
  M[,c] <- (data.m[,c]-mean(data.m[,c]))/sqrt(var(data.m[,c]));
 }
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
