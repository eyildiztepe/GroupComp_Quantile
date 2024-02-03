source(file=file.choose()) #funs.R 
source(file=file.choose()) #Rallfunv43.txt
#Rallfunv43.txt can be downloaded from https://osf.io/xhe8u/

if(!require(Hmisc))install.packages("Hmisc")
library(Hmisc)

lab1 <- c("Type-7", "HD","THD", "NO")
qi <- c(0.025,0.5,0.975)
n1 <- 20
n2 <- 20
r <- 1000
nboot <- 2000
sig.level <- matrix(0,nrow=r, ncol=4)
colnames(sig.level) <- lab1

for(k in 1:r){
  x<-ghdist(n=n1,g=0.2,h=0)
  y<-ghdist(n=n2,g=0.2,h=0)
  dif1 <- dif2 <- dif3 <- dif4 <- NULL
  dif1<-quantile(x,qi)-quantile(y,qi) #Type-7
  dif2<-hdquantile(x,qi)-hdquantile(y,qi) #Harrell-Davis
  dif3<-thdquantile(x,qi)-thdquantile(y,qi) #Trimmed Harrell-Davis
  dif4<-noquantile(x,qi)-noquantile(y,qi) # NO
  bvec1 <- bvec2 <- bvec3 <- bvec4 <-matrix(NA,nrow=nboot,ncol=length(qi))

  for(i in 1:nboot){
    xi<-sample(x,size=length(x),replace=T)
    yi<-sample(y,size=length(y),replace=T)
    bvec1[i,]<-(quantile(xi,qi)-quantile(yi,qi))
    bvec2[i,]<-(hdquantile(xi,qi)-hdquantile(yi,qi))
    bvec3[i,]<-(thdquantile(xi,qi)-thdquantile(yi,qi))
    bvec4[i,]<-(noquantile(xi,qi)-noquantile(yi,qi))
  }

  vecz<-rep(0,length(qi))
  
  smat1 <- var(bvec1)
  smat2 <- var(bvec2)
  smat3 <- var(bvec3)
  smat4 <- var(bvec4)
  
  bvec1 <- rbind(bvec1,vecz)
  bvec2 <- rbind(bvec2,vecz)
  bvec3 <- rbind(bvec3,vecz)
  bvec4 <- rbind(bvec4,vecz)

  db1<-sqrt(mahalanobis(bvec1,dif1,smat1))
  db2<-sqrt(mahalanobis(bvec2,dif2,smat2))
  db3<-sqrt(mahalanobis(bvec3,dif3,smat3))
  db4<-sqrt(mahalanobis(bvec4,dif4,smat4))
  
  sig.level[k,1] <- sum(db1[nboot+1] < db1[1:nboot]) / nboot
  sig.level[k,2] <- sum(db2[nboot+1] < db2[1:nboot]) / nboot
  sig.level[k,3] <- sum(db3[nboot+1] < db3[1:nboot]) / nboot
  sig.level[k,4] <- sum(db4[nboot+1] < db4[1:nboot]) / nboot
}
apply(sig.level,2,function(x){sum(x<0.05)})/r
