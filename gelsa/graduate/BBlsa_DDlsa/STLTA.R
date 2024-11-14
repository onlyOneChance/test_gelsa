library(qvalue)
library(forecast)
library(igraph)
library(DMwR)

### calcualte LS score
LSA<-function(x,y,D){

  n<-length(x)
  x1<-scale(qnorm(rank(x)/(n+1)))
  y1<-scale(qnorm(rank(y)/(n+1)))
  p<-matrix(0,n+1,n+1)
  N<-matrix(0,n+1,n+1)
  for (i in 1:n){
    for (j in max(1,i-D):min(i+D,n)){
      p[i+1,j+1]<-max(0,p[i,j]+x1[i]*y1[j])
      N[i+1,j+1]<-max(0,N[i,j]-x1[i]*y1[j])
    }
  }
  p<-p[-1,]
  p<-p[,-1]
  N<-N[-1,]
  N<-N[,-1]
  pm<-max(p)
  Nm<-max(N)
  s<-max(pm,Nm)
  return(s)
}
oLSA<-function(x,y,D){
  n<-length(x)
  p<-matrix(0,n+1,n+1)
  N<-matrix(0,n+1,n+1)
  for (i in 1:n){
    for (j in max(1,i-D):min(i+D,n)){
      p[i+1,j+1]<-max(0,p[i,j]+x[i]*y[j])
      N[i+1,j+1]<-max(0,N[i,j]-x[i]*y[j])
    }
  }
  p<-p[-1,]
  p<-p[,-1]
  N<-N[-1,]
  N<-N[,-1]
  pm<-max(p)
  Nm<-max(N)
  s<-max(pm,Nm)
  return(s)
}
#### generate AR(1) model
ARM<-function(n,rho){
    inter<-rep(0,n+100)
    epsilon<-rnorm(n+100,0,1)
    for (j in 1:(n+99))
    {
      inter[j+1]<-epsilon[j+1]+rho*inter[j]
    }
    return(inter[101:(n+100)])
}
library(MASS)
CAR<-function(n,rho1,rho2,rho){
  r1<-rho
  ri<-(1-rho1*rho2)*rho
  epsilon<-matrix(0,n,2)
  epsilon[1,]<-mvrnorm(1,c(0,0),matrix(c(1,rho,rho,1),2,2))
  epsilon[2:n,]<-mvrnorm(n-1,mu = c(0,0),
                         Sigma = matrix(c(1-rho1^2,ri,ri,1-rho2^2),2,2))
  inter<-matrix(0,n,2)
  inter[1,]<-epsilon[1,]
  for (j in 1:(n-1))
  {
    inter[j+1,1]<-epsilon[j+1,1]+rho1*inter[j,1]
    inter[j+1,2]<-epsilon[j+1,2]+rho2*inter[j,2]
  }
  return(inter)
}

SCAR<-function(s,n,rho1,rho2,rho){
  r1<-rho
  ri<-(1-rho1*rho2)*rho
  epsilon<-matrix(0,n+1,2)
  epsilon[1,]<-mvrnorm(1,c(0,0),matrix(c(1,rho,rho,1),2,2))
  epsilon[2:n,]<-mvrnorm(n,mu = c(0,0),
                         Sigma = matrix(c(1-rho1^2,ri,ri,1-rho2^2),2,2))
  inter<-matrix(0,n+1,2)
  inter[1,]<-s
  for (j in 1:n)
  {
    inter[j+1,1]<-epsilon[j+1,1]+rho1*inter[j,1]
    inter[j+1,2]<-epsilon[j+1,2]+rho2*inter[j,2]
  }
  return(inter[2:(n+1)])
}
###generate local correlation model
locn<-function(n,p,rho,rho1,rho2){
  nn<-ceiling(n*p)
  a<-floor((n-nn)/2)
  b<-n-nn-a
  if (a == 0) {
    x<-CAR(nn,rho1,rho2,rho)
    return(x)
  }
  else {
  x1<-CAR(a,rho1,rho2,0)
  x<-SCAR(x1[a,],nn,rho1,rho2,rho)
  x2<-SCAR(x[nn,],b,rho1,rho2,0)
  return(rbind(x1,x,x2))
  }
}
locn2<-function(n,p,rho,rho1,rho2){
  nn<-ceiling(n*p)
  a<-floor((n-nn)/2)
  b<-n-nn-a
  if (a == 0) {
    x<-CAR(nn,rho1,rho2,rho)
    return(x)
  }
  else {
  x1<-mvrnorm(a,c(0,0),matrix(c(1,0,0,1),2,2))
  x2<-mvrnorm(b,c(0,0),matrix(c(1,0,0,1),2,2))
  x<-CAR(nn,rho1,rho2,rho)
  return(rbind(x1,x,x2))
  }
}
#### calculate asympetic p-value
ld<-function(d,x){
  if (x==0) {
    ld<-1} else{
    s<-rep(0,1000)
    s[1]<-(1/x^2+1/(pi^2))*exp(-pi^2/(2*x^2))
    ps<-1
    i<-2
    while (ps>1e-5){
      ps<-(1/x^2+1/((2*i-1)^2*pi^2))*exp(-(2*i-1)^2*pi^2/(2*x^2))
      s[i]<-s[i-1]+ps
      i<-i+1
    }
  ld<-1-8^(2*d+1)*max(s)^(2*d+1)
  }
  return(ld)
}

### omega estimate
omega<-function(data){
  n<-length(data)
  autocov<-acf(data,type="covariance",plot=FALSE,lag.max=n-1)[[1]][,1,1]
  demeandata<-data-mean(data)
  rhohat<-sum(demeandata[2:n]*demeandata[1:(n-1)])/sum((demeandata[1:(n-1)])^2)
  alphahat<-4*rhohat^2/(1-rhohat^2)^2
  lnDDB<-ceiling(1.1447*(alphahat*n)^(1/3))
  if (lnDDB>1) sigmasqDDB<-sum(c(autocov[lnDDB:1],autocov[2:lnDDB])*(1-abs(seq(-lnDDB+1,lnDDB-1)/lnDDB))) else  sigmasqDDB<-autocov[1]
  return(sigmasqDDB)
}


#### bandwidth select
bw<-function(data){
  n<-length(data)
  demeandata<-data-mean(data)
  rhohat<-sum(demeandata[2:n]*demeandata[1:(n-1)])/sum((demeandata[1:(n-1)])^2)
  alphahat<-4*rhohat^2/(1-rhohat^2)^2
  l<-ceiling(1.1447*(alphahat*n)^(1/3))
  return(l)
}

#### adjusted omega estimate
omega2<-function(x,y){
  n<-length(x)
  b<-bw(x*y)
  cov1<-acf(x,type="covariance",plot=FALSE,lag.max=n-1)[[1]][,1,1]
  cov2<-acf(y,type="covariance",plot=FALSE,lag.max=n-1)[[1]][,1,1]
  cov3<-cov1*cov2
  if (!is.na(b)){
  if (b>1) sigma<-sum(c(cov3[b:1],cov3[2:b])*(1-abs(seq(-b+1,b-1)/b)))
  else sigma<-cov3[1]
  }
  else sigma<-cov3[1]
  return(sigma)
}

#### all in
aLSA<-function(x,y,D){
  n<-length(x)
  x1<-scale(x)
  y1<-scale(y)
  score<-oLSA(x1,y1,D)
  ome<-omega2(x1,y1)
  if (is.na(ome)) app<-ld(D,score/sqrt(n)) else {
  ap<-score/sqrt(ome*n)
  app<-ld(D,ap)
  }
  return(app)
}

###size
IARM<-function(n,alpha,rho){
    inter<-rep(0,n+100)
    epsilon<-rnorm(n+100,0,1)
    for (j in 1:(n+99))
    {
      inter[j+1]<-epsilon[j+1]+rho*inter[j]+alpha
    }
    return(inter[101:(n+100)])
}


library(markovchain)
LSA<-function(x, y, D){
  n<-length(x)
  P<-matrix(0, n+1, n+1)
  N<-matrix(0, n+1, n+1)
  for (i in 1:n){
    for (j in max(1, i-D):min(i+D, n)){
      P[i+1, j+1] <- max(0, P[i, j] + x[i] * y[j])
      N[i+1, j+1] <- max(0, N[i, j] - x[i] * y[j])
    }
  }
  P <- P[-1,]
  P <- P[,-1]
  N <- N[-1,]
  N <- N[,-1]
  Pm <- max(P)
  Nm <- max(N)
  s <- max(Pm, Nm)
  return(s)
}
#### calculate asympetic p-value
AsypeticSignificance <- function(d, x){
  if (x == 0) {
    AsypeticSignificance <- 1
    } else{
      s <- rep(0, 1000)
      s[1] <- (1 / x^2 + 1 / (pi^2)) * exp(-pi^2 / (2 * x^2))
      ps <- 1
      i <- 2
      while (ps > 1e-5){
        ps <- (1 / x^2 + 1 / ((2 * i-1)^2 * pi^2)) * exp(-(2 * i - 1)^2 * pi^2 / (2 * x^2))
        s[i] <- s[i-1] + ps
        i <- i + 1
      }
      AsypeticSignificance <- 1 - 8^(2 * d + 1) * max(s)^(2 * d + 1)
    }
  return(AsypeticSignificance)
}

bw <- function(data){
  n <- length(data)
  demeandata <- data - mean(data)
  rhohat <- sum(demeandata[2:n] * demeandata[1:(n-1)])/sum((demeandata[1:(n-1)])^2)
  alphahat <- 4 * rhohat^2 / (1 - rhohat^2)^2
  l <- ceiling(1.1447 * (alphahat * n)^(1/3))
  return(l)
}

omega<-function(x,y){
  n<-length(x)
  b<-bw(x*y)
  cov1<-acf(x,type="covariance",plot=FALSE,lag.max=n-1)[[1]][,1,1]
  cov2<-acf(y,type="covariance",plot=FALSE,lag.max=n-1)[[1]][,1,1]
  cov3<-cov1*cov2
  if (!is.na(b)){
    if (b>1) {
      b <- min(b, n)
      sigma<-sum(c(cov3[b:1],cov3[2:b])*(1-abs(seq(-b+1,b-1)/b)))
    }else sigma<-cov3[1]
  }
  else sigma<-cov3[1]
  return(sigma)
}

DDLTA<-function(x,y,D){
  n<-length(x)
  score<-LSA(x,y,D)
  ome<-omega(x,y)
  if (is.na(ome)) app<-AsypeticSignificance(D,score/sqrt(n)) else {
    ap<-score/sqrt(ome*n)
    app<-AsypeticSignificance(D,ap)
  }
  return(app)
}

discretization <- function(x,t){
  n <- length(x)
  b <- diff(x)/abs(x)[-n]
  b[b>=t] <- 1
  b[b<= (-1)*t] <- -1
  b[b>-t & b< t] <- 0
  for (i in which(!complete.cases(b))){
    b[i] <- sign(x[i+1])
  }
  return(b)
}

convergence_matrix <- function(A){
  delta <- 1
  B <- A
  while (delta > 10e-10){
    A <- A%*%A
    delta <- sum(abs(A - B))
    B <- A
  }
  return(A)
}

markovchain_sigma_2d <- function(ax, ay){
  return(1 + 2 * (2 * ax - 1) * (2 * ay - 1) / (1 - (2 * ax - 1) * (2 * ay - 1)))
}

markovchain_sigma_3d <- function(bx, dx, cx, by, dy, cy) {
  return(4*(dx / (1 - bx - cx + 2 * dx)) * (dy / (1 - by - cy + 2 * dy)) 
         * (1 + 2 * (bx - cx) * (by - cy) / (1 - (bx - cx) * (by - cy))))
}

markovchain_sigma_2d_3d <- function(ax, by, dy, cy) {
  return(2 * (dy / (1 - by - cy + 2 * dy)) 
         * (1 + 2 * (2 * ax - 1) * (by - cy) / (1 - (2 * ax - 1) * (by - cy))))
}

Sigma <- function(x,y){
  T_X <- markovchainFit(data = x)$estimate@transitionMatrix
  T_Y <- markovchainFit(data = y)$estimate@transitionMatrix
  if (dim(T_X)[1] == 2 & dim(T_Y)[1] == 2){
    ax <- (T_X[1,1] + T_X[2,2]) / 2
    ay <- (T_Y[1,1] + T_Y[2,2]) / 2
    return(markovchain_sigma_2d(ax, ay))
  }else{
    if (dim(T_X)[1] == 3 & dim(T_Y)[1] == 3){
      bx <- (T_X[1,1] + T_X[3,3]) / 2
      dx <- (T_X[2,1] + T_X[2,3]) / 2
      cx <- (T_X[3,1] + T_X[1,3]) / 2
      by <- (T_Y[1,1] + T_Y[3,3]) / 2
      dy <- (T_Y[2,1] + T_Y[2,3]) / 2
      cy <- (T_Y[3,1] + T_Y[1,3]) / 2
      return(markovchain_sigma_3d(bx, dx, cx, by, dy, cy))
    }else{
      if (dim(T_X)[1] == 3 & dim(T_Y)[1] == 2){
        bx <- (T_X[1,1] + T_X[3,3]) / 2
        dx <- (T_X[2,1] + T_X[2,3]) / 2
        cx <- (T_X[3,1] + T_X[1,3]) / 2
        ay <- (T_Y[1,1] + T_Y[2,2]) / 2
        return(markovchain_sigma_2d_3d(ay, bx, dx, cx))
      }
      if (dim(T_X)[1] == 2 & dim(T_Y)[1] == 3){
        ax <- (T_X[1,1] + T_X[2,2]) / 2
        by <- (T_Y[1,1] + T_Y[3,3]) / 2
        dy <- (T_Y[2,1] + T_Y[2,3]) / 2
        cy <- (T_Y[3,1] + T_Y[1,3]) / 2
        return(markovchain_sigma_2d_3d(ax, by, dy, cy))
      }
    }
  }
}


per <- function(x, y, D, t){
  n <- length(x)
  dx <- discretization(x, t)
  dy <- discretization(y, t)
  z <- LSA(dx, dy, D)
  a <- rep(0, 1000)
  for (k in 1:1000){
    a[k] <- LSA(discretization(sample(x,n), t), 
                discretization(sample(y,n), t), D)
  }
  l <- length(a[a > z]) / 1000
  return(l)
}

M3fece <- read.csv('M3_feces.csv',row.names = 1)
M3fece <- as.matrix(M3fece[apply(M3fece>0,1,sum)>=0.6*dim(M3fece)[2],])

M3_fece <- t(apply(M3fece, 1, function(x) x / apply(M3fece, 2, sum)))
n<-dim(M3_fece)[2]
sam_num<-dim(M3_fece)[1]

acfP<-rep(0,sam_num)
for (i in 1:sam_num){
  acfP[i]<-acf(M3_fece[i,],plot=FALSE)[[1]][2]
}
sum(acfP>1.96/sqrt(n))

BJtest = rep(0, sam_num)
for (i in 1:sam_num){
  BJtest[i] = Box.test(M3_fece[i,])$p.value
}
sum(BJtest<0.05)
which(BJtest<0.05)

t <- 0.5
Perm_Pvalue<-matrix(0,sam_num,sam_num)
TLTA_Pvalue<-matrix(0,sam_num,sam_num)
TDLTA_Pvalue<-matrix(0,sam_num,sam_num)
DD<-matrix(0,sam_num,sam_num)
for (i in 2:sam_num){
  for (j in 1:(i-1)){
    x <- M3_fece[i,] - mean(M3_fece[i,])
    y <- M3_fece[j,] - mean(M3_fece[j,])
    dx <- discretization(x, t)
    dy <- discretization(y, t)
    TLTA_Pvalue[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * 1.25))
    TDLTA_Pvalue[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * Sigma(dx, dy)))
    Perm_Pvalue[i,j] <- per(x, y, 3, t)
  }
}
TLTA_qvalue<-qvalue(TLTA_Pvalue[lower.tri(TLTA_Pvalue)])
TDLTA_qvalue<-qvalue(TDLTA_Pvalue[lower.tri(TDLTA_Pvalue)])
Perm_qvalue<-qvalue(Perm_Pvalue[lower.tri(Perm_Pvalue)])
summary(TLTA_qvalue)
summary(TDLTA_qvalue)
summary(Perm_qvalue)

n1<-sum(Perm_qvalue$qvalues<0.05)
n2<-sum(TLTA_qvalue$qvalues<0.05)
n3<-sum(TDLTA_qvalue$qvalues<0.05)
n12<-sum(Perm_qvalue$qvalues<0.05 & TLTA_qvalue$qvalues<0.05)
n13<-sum(Perm_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)
n23<-sum(TLTA_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)
n123<-sum(Perm_qvalue$qvalues<0.05 & TLTA_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)

library(VennDiagram)
draw.triple.venn(area1=n1,area2=n2,area3=n3,n12=n12,n23=n23,n13=n13,n123=n123,
                 category=c("Permutation","TLTA","TDLTA"),
                 fill = c('green','blue','red'),cex.prop=TRUE)

nn1<-sum(Perm_qvalue$qvalues<0.01)
nn2<-sum(TLTA_qvalue$qvalues<0.01)
nn3<-sum(TDLTA_qvalue$qvalues<0.01)
nn12<-sum(Perm_qvalue$qvalues<0.01&TLTA_qvalue$qvalues<0.01)
nn13<-sum(Perm_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)
nn23<-sum(TLTA_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)
nn123<-sum(Perm_qvalue$qvalues<0.01&TLTA_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)

library(VennDiagram)
draw.triple.venn(area1=nn1,area2=nn2,area3=nn3,n12=nn12,n23=nn23,n13=nn13,n123=nn123,
                 category=c("Permutation","TLTA","TDLTA"),
                 fill = c('green','blue','red'))


t <- 0
Perm_Pvalue2<-matrix(0,sam_num,sam_num)
TLTA_Pvalue2<-matrix(0,sam_num,sam_num)
TDLTA_Pvalue2<-matrix(0,sam_num,sam_num)
for (i in 2:sam_num){
  for (j in 1:(i-1)){
    x <- M3_fece[i,] - mean(M3_fece[i,])
    y <- M3_fece[j,] - mean(M3_fece[j,])
    dx <- discretization(x, t)
    dy <- discretization(y, t)
    TLTA_Pvalue2[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * 1.25))
    TDLTA_Pvalue2[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * Sigma(dx, dy)))
    # DD2[i,j]<-aLSA(M3_fece[i,],M3_fece[j,],3)
    Perm_Pvalue2[i,j] <- per(x, y, 3, t)
  }
}
TLTA_qvalue2<-qvalue(TLTA_Pvalue2[lower.tri(TLTA_Pvalue2)])
TDLTA_qvalue2<-qvalue(TDLTA_Pvalue2[lower.tri(TDLTA_Pvalue2)])
Perm_qvalue2<-qvalue(Perm_Pvalue2[lower.tri(Perm_Pvalue2)])
summary(TLTA_qvalue2)
summary(TDLTA_qvalue2)
summary(Perm_qvalue2)

n1<-sum(Perm_qvalue2$qvalues<0.05)
n2<-sum(TLTA_qvalue2$qvalues<0.05)
n3<-sum(TDLTA_qvalue2$qvalues<0.05)
n12<-sum(Perm_qvalue2$qvalues<0.05 & TLTA_qvalue2$qvalues<0.05)
n13<-sum(Perm_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)
n23<-sum(TLTA_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)
n123<-sum(Perm_qvalue2$qvalues<0.05 & TLTA_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)

library(VennDiagram)
draw.triple.venn(area1=n1,area2=n2,area3=n3,n12=n12,n23=n23,n13=n13,n123=n123,
                 category=c("Permutation","TLTA","TDLTA"),
                 fill = c('green','blue','red'),cex.prop=TRUE)

env<-read.csv('env.csv')
envf<-env[1]
env<-as.matrix(env[-1])
anyNA(env)
index <- which(is.na(env),arr.ind = TRUE)
for (i in 1:dim(index)[1]){
  env[index[i,1],index[i,2]] = sum(env[index[i,1],index[i,2]-1],env[index[i,1],index[i,2]+1])/2
}

bac<-read.csv('PML.csv')
OTU<-bac[1]
colnames(envf)<-colnames(OTU)
fac<-rbind(envf,OTU)
P_matrix<-as.matrix(bac[-1])
All <- rbind(env, t(apply(P_matrix, 1, function(x) x / apply(P_matrix, 2, sum))))
n<-dim(All)[2]
sam_num<-dim(All)[1]
rownames(All) <- lapply(fac, as.character)$taxa

acfP<-rep(0,sam_num)
for (i in 1:sam_num){
  acfP[i]<-acf(All[i,],plot=FALSE)[[1]][2]
}
sum(acfP>1.96/sqrt(n))

BJtest = rep(0, sam_num)
for (i in 1:sam_num){
  BJtest[i] = Box.test(All[i,])$p.value
}
sum(BJtest<0.05)
which(BJtest<0.05)

t <- 0.5
Perm_Pvalue<-matrix(0,sam_num,sam_num)
TLTA_Pvalue<-matrix(0,sam_num,sam_num)
TDLTA_Pvalue<-matrix(0,sam_num,sam_num)
for (i in 2:sam_num){
  for (j in 1:(i-1)){
    x <- All[i,] - mean(All[i,])
    y <- All[j,] - mean(All[j,])
    dx <- discretization(x, t)
    dy <- discretization(y, t)
    TLTA_Pvalue[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * 1.25))
    TDLTA_Pvalue[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * Sigma(dx, dy)))
    Perm_Pvalue[i,j] <- per(x, y, 3, t)
  }
}
TLTA_qvalue<-qvalue(TLTA_Pvalue[lower.tri(TLTA_Pvalue)])
TDLTA_qvalue<-qvalue(TDLTA_Pvalue[lower.tri(TDLTA_Pvalue)])
Perm_qvalue<-qvalue(Perm_Pvalue[lower.tri(Perm_Pvalue)])
summary(TLTA_qvalue)
summary(TDLTA_qvalue)
summary(Perm_qvalue)

n1<-sum(Perm_qvalue$qvalues<0.05)
n2<-sum(TLTA_qvalue$qvalues<0.05)
n3<-sum(TDLTA_qvalue$qvalues<0.05)
n12<-sum(Perm_qvalue$qvalues<0.05 & TLTA_qvalue$qvalues<0.05)
n13<-sum(Perm_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)
n23<-sum(TLTA_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)
n123<-sum(Perm_qvalue$qvalues<0.05 & TLTA_qvalue$qvalues<0.05 & TDLTA_qvalue$qvalues<0.05)

library(VennDiagram)
draw.triple.venn(area1=n1,area2=n2,area3=n3,n12=n12,n23=n23,n13=n13,n123=n123,
                 category=c("Permutation","TLTA","TDLTA"),
                 fill = c('green','blue','red'),cex.prop=TRUE)

nn1<-sum(Perm_qvalue$qvalues<0.01)
nn2<-sum(TLTA_qvalue$qvalues<0.01)
nn3<-sum(TDLTA_qvalue$qvalues<0.01)
nn12<-sum(Perm_qvalue$qvalues<0.01&TLTA_qvalue$qvalues<0.01)
nn13<-sum(Perm_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)
nn23<-sum(TLTA_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)
nn123<-sum(Perm_qvalue$qvalues<0.01&TLTA_qvalue$qvalues<0.01&TDLTA_qvalue$qvalues<0.01)

library(VennDiagram)
draw.triple.venn(area1=nn1,area2=nn2,area3=nn3,n12=nn12,n23=nn23,n13=nn13,n123=nn123,
                 category=c("Permutation","TLTA","TDLTA"))


t <- 0
Perm_Pvalue2<-matrix(0,sam_num,sam_num)
TLTA_Pvalue2<-matrix(0,sam_num,sam_num)
TDLTA_Pvalue2<-matrix(0,sam_num,sam_num)
for (i in 2:sam_num){
  for (j in 1:(i-1)){
    x <- All[i,] - mean(All[i,])
    y <- All[j,] - mean(All[j,])
    dx <- discretization(x, t)
    dy <- discretization(y, t)
    TLTA_Pvalue2[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * 1.25))
    TDLTA_Pvalue2[i,j] <- AsypeticSignificance(3, LSA(dx, dy, 3) / sqrt(n * Sigma(dx, dy)))
    #Perm_Pvalue2[i,j] <- per(x, y, 3, t)
  }
}
TLTA_qvalue2<-qvalue(TLTA_Pvalue2[lower.tri(TLTA_Pvalue2)])
TDLTA_qvalue2<-qvalue(TDLTA_Pvalue2[lower.tri(TDLTA_Pvalue2)])
Perm_qvalue2<-qvalue(Perm_Pvalue2[lower.tri(Perm_Pvalue2)])
summary(TLTA_qvalue2)
summary(TDLTA_qvalue2)
summary(Perm_qvalue2)

n1<-sum(Perm_qvalue2$qvalues<0.05)
n2<-sum(TLTA_qvalue2$qvalues<0.05)
n3<-sum(TDLTA_qvalue2$qvalues<0.05)
n12<-sum(Perm_qvalue2$qvalues<0.05 & TLTA_qvalue2$qvalues<0.05)
n13<-sum(Perm_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)
n23<-sum(TLTA_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)
n123<-sum(Perm_qvalue2$qvalues<0.05 & TLTA_qvalue2$qvalues<0.05 & TDLTA_qvalue2$qvalues<0.05)

library(VennDiagram)
draw.triple.venn(area1=n1,area2=n2,area3=n3,n12=n12,n23=n23,n13=n13,n123=n123,
                 category=c("Permutation","TLTA","TDLTA"),
                 fill = c('green','blue','red'),cex.prop=TRUE)