CheckBMIQ <-
function(beta.v,design.v,pnbeta.v){### pnbeta is BMIQ normalised profile

type1.idx <- which(design.v==1);
type2.idx <- which(design.v==2);

beta1.v <- beta.v[type1.idx];
beta2.v <- beta.v[type2.idx];
pnbeta2.v <- pnbeta.v[type2.idx];
  

}
