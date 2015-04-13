DoPBC <-
function(beta.m,design.v){

  mval.m <- log2(beta.m/(1-beta.m));
  type1.idx <- which(design.v==1)
  type2.idx <- which(design.v==2)  
  mvalT.m <- mval.m;
  for(s in 1:ncol(beta.m)){

    neg.idx <- which(mval.m[,s]<0);
    pos.idx <- which(mval.m[,s]>0);

    neg1.idx <- intersect(neg.idx,type1.idx);
    neg2.idx <- intersect(neg.idx,type2.idx);

    pos1.idx <- intersect(pos.idx,type1.idx);
    pos2.idx <- intersect(pos.idx,type2.idx);

    d.o <- density(mval.m[neg2.idx,s],kernel="gaussian",bw=0.5);
    peakU2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos2.idx,s],kernel="gaussian",bw=0.5);
    peakM2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[neg1.idx,s],kernel="gaussian",bw=0.5);
    peakU1 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos1.idx,s],kernel="gaussian",bw=0.5);
    peakM1 <- abs(d.o$x[which.max(d.o$y)]);

    mvalT.m[neg2.idx,s] <- (mval.m[neg2.idx,s]/peakU2)*peakU1;
    mvalT.m[pos2.idx,s] <- (mval.m[pos2.idx,s]/peakM2)*peakM1;
    print(paste("Done for sample ",s,sep=""));
  }

  betaT.m <- 2^mvalT.m/(2^mvalT.m+1);
  return(betaT.m);
}
