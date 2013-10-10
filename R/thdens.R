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
