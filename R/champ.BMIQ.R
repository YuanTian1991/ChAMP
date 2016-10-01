### BMIQ.R & CheckBMIQ.R
### This function adjusts for the type-2 bias in Illumina Infinium 450k data.
### Author: Andrew Teschendorff
### Date v_1.1: Nov 2012
### Date v_1.2: 6th Apr 2013
### Date v_1.3: 29th May 2013
### Date v_1.4: 23rd Sep 2014(corrected minor inconsequential bug)
### Date v_1.5: 19th Oct 2015 (corrected minor inconsequential bug)

### SUMMARY
### BMIQ is an intra-sample normalisation procedure, adjusting for the bias in type-2 probe values, using a 3-step procedure published in Teschendorff AE et al "A Beta-Mixture Quantile Normalisation method for correcting probe design bias in Illumina Infinium 450k DNA methylation data", Bioinformatics 2012 Nov 21.


### INPUT:
### beta.v: vector consisting of beta-values for a given sample. NAs are not allowed, so these must be removed or imputed prior to running BMIQ. Beta-values that are exactly 0 or 1 will be replaced by the minimum positive above zero or maximum value below 1, respectively.
### design.v: corresponding vector specifying probe design type (1=type1,2=type2). This must be of the same length as beta.v and in the same order.
### doH: perform normalisation for hemimethylated type2 probes. By default TRUE.
### nfit: number of probes of a given design to use for the fitting. Default is 10000. Smaller values will make BMIQ run faster at the expense of a small loss in accuracy. For most applications, even 5000 is ok.
### nL: number of states in beta mixture model. 3 by default. At present BMIQ only works for nL=3.
### th1.v: thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type1 probes hemi-methylated and methylated, and will be refined by the EM algorithm. Default values work well in most cases.
### th2.v: thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type2 probes hemi-methylated and methylated, and will be refined by the EM algorithm. By default this is null, and the thresholds are estimated based on th1.v and a modified PBC correction method.
### niter: maximum number of EM iterations to do. By default 5.
### tol: tolerance threshold for EM algorithm. By default 0.001.
### plots: logical specifying whether to plot the fits and normalised profiles out. By default TRUE.
### sampleID: the ID of the sample being normalised.

### OUTPUT
### A list with the following elements:
### nbeta: the normalised beta-profile for the sample
### class1: the assigned methylation state of type1 probes
### class2: the assigned methylation state of type2 probes
### av1: mean beta-values for the nL classes for type1 probes.
### av2: mean beta-values for the nL classes for type2 probes.
### hf: the "Hubble" dilation factor
### th1: estimated thresholds used for type1 probes
### th2: estimated thresholds used for type2 probes

require(RPMM);

champ.BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=10000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){

type1.idx <- which(design.v==1);
type2.idx <- which(design.v==2);

beta1.v <- beta.v[type1.idx];
beta2.v <- beta.v[type2.idx];

### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
if(min(beta1.v)==0){
  beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
}
if(min(beta2.v)==0){
  beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
}
if(max(beta1.v)==1){
  beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
}
if(max(beta2.v)==1){
  beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
}

### estimate initial weight matrix from type1 distribution
w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
w0.m[which(beta1.v > th1.v[2]),3] <- 1;

### fit type1
print("Fitting EM beta mixture to type1 probes");
rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
em1.o <- blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
subsetclass1.v <- apply(em1.o$w,1,which.max);
subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]]))),mean(c(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]]))));
class1.v <- rep(2,length(beta1.v));
class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
nth1.v <- subsetth1.v;
print("Done");

### generate plot from estimated mixture
if(plots){
print("Check");
tmpL.v <- as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
tmpB.v <- vector();
for(l in 1:nL){
  tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
}

pdf(paste("Type1fit-",sampleID,".pdf",sep=""),width=6,height=4);
plot(density(beta1.v));
d.o <- density(tmpB.v);
points(d.o$x,d.o$y,col="green",type="l")
legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
dev.off();
}



### Estimate Modes 
d1U.o <- density(beta1.v[class1.v==1])
d1M.o <- density(beta1.v[class1.v==3])
mod1U <- d1U.o$x[which.max(d1U.o$y)]
mod1M <- d1M.o$x[which.max(d1M.o$y)]
d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
d2M.o <- density(beta2.v[which(beta2.v>0.6)]);
mod2U <- d2U.o$x[which.max(d2U.o$y)]
mod2M <- d2M.o$x[which.max(d2M.o$y)]


### now deal with type2 fit
th2.v <- vector();
th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
th2.v[2] <- nth1.v[2] + (mod2M-mod1M);

### estimate initial weight matrix 
w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
w0.m[which(beta2.v > th2.v[2]),3] <- 1;

print("Fitting EM beta mixture to type2 probes");
rand.idx <- sample(1:length(beta2.v),nfit,replace=FALSE)
em2.o <- blc(matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
print("Done");

### for type II probes assign to state (unmethylated, hemi or full methylation)
subsetclass2.v <- apply(em2.o$w,1,which.max);
subsetth2.v <- c(mean(c(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]]))),mean(c(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]]))));
class2.v <- rep(2,length(beta2.v));
class2.v[which(beta2.v < subsetth2.v[1])] <- 1;
class2.v[which(beta2.v > subsetth2.v[2])] <- 3;


### generate plot
if(plots){
tmpL.v <- as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
tmpB.v <- vector();
for(lt in 1:nL){
  tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
}
pdf(paste("Type2fit-",sampleID,".pdf",sep=""),width=6,height=4);
plot(density(beta2.v));
d.o <- density(tmpB.v);
points(d.o$x,d.o$y,col="green",type="l")
legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
dev.off();
}

classAV1.v <- vector();classAV2.v <- vector();
for(l in 1:nL){
  classAV1.v[l] <-  em1.o$mu[l,1];
  classAV2.v[l] <-  em2.o$mu[l,1];
}

### start normalising type2 probes
print("Start normalising type 2 probes");
nbeta2.v <- beta2.v;
### select U probes
lt <- 1;
selU.idx <- which(class2.v==lt);
selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
### find prob according to typeII distribution
p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
nbeta2.v[selUR.idx] <- q.v;
p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
nbeta2.v[selUL.idx] <- q.v;

### select M probes
lt <- 3;
selM.idx <- which(class2.v==lt);
selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
### find prob according to typeII distribution
p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
nbeta2.v[selMR.idx] <- q.v;


if(doH){ ### if TRUE also correct type2 hemimethylated probes
### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
lt <- 2;
selH.idx <- c(which(class2.v==lt),selML.idx);
minH <- min(beta2.v[selH.idx])
maxH <- max(beta2.v[selH.idx])
deltaH <- maxH - minH;
#### need to do some patching
deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])

## new maximum of H probes should be
nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
## new minimum of H probes should be
nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
ndeltaH <- nmaxH - nminH;

### perform conformal transformation (shift+dilation)
## new_beta_H(i) = a + hf*(beta_H(i)-minH);
hf <- ndeltaH/deltaH ;
### fix lower point first
nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);

}

pnbeta.v <- beta.v;
pnbeta.v[type1.idx] <- beta1.v;
pnbeta.v[type2.idx] <- nbeta2.v;

### generate final plot to check normalisation
if(plots){
 print("Generating final plot");
 d1.o <- density(beta1.v);
 d2.o <- density(beta2.v);
 d2n.o <- density(nbeta2.v);
 ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
 pdf(paste("CheckBMIQ-",sampleID,".pdf",sep=""),width=6,height=4)
 plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
 points(d1.o$x,d1.o$y,col="red",type="l");
 points(d2n.o$x,d2n.o$y,col="blue",type="l");
 legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
 dev.off();
}

print(paste("Finished for sample ",sampleID,sep=""));

return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));

}



CheckBMIQ <- function(beta.v,design.v,pnbeta.v){### pnbeta is BMIQ normalised profile

type1.idx <- which(design.v==1);
type2.idx <- which(design.v==2);

beta1.v <- beta.v[type1.idx];
beta2.v <- beta.v[type2.idx];
pnbeta2.v <- pnbeta.v[type2.idx];
  

}
