EstDimRMTv2 <-
function(data.m){

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

thdens <-
function(Q,sigma2,ns){

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
