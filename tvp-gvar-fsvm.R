setwd("!tvp-gvar")

require(stochvol)
require(mvnfast)
require(MASS)
require(Rcpp)
require(GIGrvg)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("threshold_functions.cpp")
source("aux_func.R")
load("data_raw.rda")

Yraw <- Yraw[,grepl(paste0(c("CPI","UNR","INP","LIR","EQP","SIR","FCI"),collapse="|"),colnames(Yraw))]
Yraw <- window(Yraw,start=c(1995,1),end=c(2019,12))
Ysd <- apply(Yraw,2,sd)
Yraw <- apply(Yraw,2,function(x){(x-mean(x))/sd(x)})

# select some variables
run <- as.integer(Sys.getenv("SGE_TASK_ID"))
ids <- c(1:2)
lags <- c(3) # c(1,2,3)
facs <- c(10,15)#c(1,2,6,10,15,20) # seq(4,8,by=2)
vola <- 0.2
pr_imp <- 1
consc <- FALSE
fsvrw <- TRUE

combi <- expand.grid("id"=ids,"p"=lags,"fac"=facs,"vol"=vola,"pri"=pr_imp,"consc"=consc,"fsvrw"=fsvrw,stringsAsFactors = F)
sl.run <- combi[run,]

start <- "1995-01-01"
end <- "2019-12-01"
dates <- seq(from=as.Date(start),to=as.Date(end),by="month")

cons.cons <- sl.run[["consc"]]
fsv.rw <- sl.run[["fsvrw"]]

dir.create("results/",showWarnings = FALSE)
pname <- paste0("GVAR_p",sl.run[["p"]],"_fac",sl.run[["fac"]],"_id",sl.run[["id"]])

# -----------------------------------------------------------------------------------------------
# MCMC setup
check.time <- TRUE # indicate progress bar
p <- sl.run[["p"]] # lags
q <- sl.run[["fac"]] # factors
cons <- TRUE # include constant
dates <- dates[(p+1):length(dates)] # adjust for lag-length

nburn <- 2000
nsave <- 8000
thinfac <- 1/4
nthin <- round(thinfac * nsave)

ntot <- nburn + nsave
thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
nhor <- 37
in.thin <- 0

tvpar <- TRUE # time varying parameter model
shrink.lambda <- TRUE # shrinkage on the factor loadings
svmean.exp <- FALSE # include the volatilities in the mean as exp() or not

toplot <- seq(1,ntot,by=10)

# -----------------------------------------------------------------------------------------------
# data transformations
Y <- Yraw[(p+1):nrow(Yraw),]
m <- ncol(Y)
T <- nrow(Y)

cN <- unique(gsub("\\..*","",colnames(Yraw))) # names of the countries
N <- length(cN)
vN <- unique(gsub(".*\\.","",colnames(Yraw))) # names of the variables
G <- length(vN) # number of variables per country
k <- 2*G*p+cons
id.unc <- k+1

ident.countries <- matrix(seq(1,m),G,N)
ident.label <- as.vector(ident.countries[1,])

X <- mlag(Yraw,p)
if (cons){
  X <- cbind(X[(p+1):nrow(X),],1,matrix(rnorm(T),T,1))
}else{
  X <- cbind(X[(p+1):nrow(X),],matrix(rnorm(T),T,1))
}
K <- ncol(X)

lab.names <- NULL
for (jj in 1:p) lab.names <- c(lab.names,paste0(colnames(Y),"lag_",jj))
if (cons){
  colnames(X) <- c(lab.names,"cons","UNC")
}else{
  colnames(X) <- c(lab.names,"UNC")
}

ident.matrix <- matrix(NA,K,m)
rownames(ident.matrix) <- colnames(X)
colnames(ident.matrix) <- colnames(Y)

for (jj in 1:m){
  c.select <- substr(colnames(ident.matrix)[[jj]],1,2)
  var.select <- substr(colnames(ident.matrix)[[jj]],4,6)
  sl.all <- grepl(c.select,rownames(ident.matrix))
  
  ident.matrix[sl.all,jj] <- var.select
}
ident.country <- t(matrix(seq(1,m),G,N))

# -----------------------------------------------------------------------------------------------
# prior setup
A_prior <- At_prior <- matrix(0,k+1,m) # GVAR representation coefficients
A_draw <- array(0,dim=c(k+1,m))
At_draw <- array(0,dim=c(T,k+1,m))

theta  <- matrix(1,k+1,m)
omega_mat <- matrix(0.01,k+1,m)
omega_pr <- array(1,dim=c(k+1,m))

# -----------------------------------------------------------------------------------------------
# OLS quantities and initial values
A_OLS <- ginv(crossprod(X))%*%crossprod(X,Y)
SIGMA_OLS <- crossprod(Y-X%*%A_OLS)/(T-K)
epsilon <- eta <-  Y-X%*%A_OLS

# full system coefficients
PHI_draw_t <- array(0,c(T,K,m))
dimnames(PHI_draw_t) <- list(as.character(dates),colnames(X),colnames(Y))

sv_draw <- list()
sv_latent <- list()
sv_priors <- list()

for (mm in seq_len(m)){
  sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  sv_latent[[mm]] <- rep(0,T)
  sv_priors[[mm]] <- specify_priors(
    mu = sv_normal(mean = 0, sd = 10),
    phi = sv_beta(50,1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 0.5),
    nu = sv_infinity(),
    rho = sv_constant(0)
  )
}

# factor structure in the error terms
Lambda <- matrix(0,m,q) # to store factor loadings
theta.lambda <- matrix(1,m,q)
id.lambda <- rbind(lower.tri(matrix(1,q,q)),matrix(TRUE,m-q,q))
q.free <- sum(id.lambda)

F_PC <- princomp(epsilon)
F_draw <- F_PC$scores[,1:q,drop=FALSE]
ident.fac <- matrix(1:q,q,N)

svtemp <- matrix(NA,T,q)
for(qq in 1:q){
  svtemp[,qq] <- as.numeric(apply(latent(svsample(F_draw[,qq],draws=1000,burnin=1000,priorphi=c(50,1.5))),2,median))
}

# volatilities in the error term
Sig.t <- matrix(1,T,m) # idiosyncratic variances
H <- matrix(exp(apply(svtemp,1,mean)),T,1)
H.acc <- H*0
sig.H <- sl.run[["vol"]] # variance of the SV process of the factor

sv_var <- rep(1,m)
fsv_para <- c(0,1,0.1)

# Horseshoe for VAR coefficients
A.lambda <- A_draw^0
A.nu <- A_draw^0
A.tau <- 1
A.zeta <- 1

# Horseshoe for state innovation variances
O.lambda <- A_draw^0
O.nu <- A_draw^0
O.tau <- 1
O.zeta <- 1

# Horseshoe for factor loadings
L.lambda <- Lambda^0
L.nu <- Lambda^0
L.tau <- 1
L.zeta <- 1

# -----------------------------------------------------------------------------------------------
# storage objects
PHI_store <- array(NA,dim=c(nthin,T,K,m))
Sigma_store <- array(NA,dim=c(nthin,T,m))
H_store <- array(NA,dim=c(nthin,T,1))

# -----------------------------------------------------------------------------------------------
# Gibbs sampling
if(check.time == TRUE){
  t.start <- Sys.time()
  pb <- txtProgressBar(min = 1, max = ntot, style = 3)
}

for (irep in 1:ntot){
  # Sample GVAR and full system coefficients
  Y_ <- Y - tcrossprod(F_draw,Lambda)
  for(mm in 1:m){
    sl.mm <- which(apply(ident.countries==mm,2,sum)==1)
    
    # draw time-varying coefficients
    W_i <- create.W_lag(sl.C=sl.mm,G=G,N=N,p=p,cons=cons,ex=1,weights=W[sl.mm,-sl.mm])
    X_iraw <- tcrossprod(X,W_i)
    Y_i <- Y_[,mm] - (X_iraw %*% A_draw)[,mm]
    X_i <- X_iraw %*% diag(omega_mat[,mm])
    
    if(tvpar){
      At_draw[,,mm] <- t(KF_fast(t(as.matrix(Y_i)),X_i,as.matrix(exp(Sig.t[,mm])),t(matrix(1,k+1,T)),k+1,1,T,matrix(0,k+1,1),diag(k+1)*1/10^15))
    }else{
      At_draw[,,mm] <- matrix(0,T,ncol(X_i))
    }
    
    # draw constant coefficients and shrinkage parameters
    normalizer <- exp(-0.5*Sig.t[,mm])
    Y_i <- Y_[,mm]*normalizer
    X_i <- cbind(X_iraw,X_iraw*At_draw[,,mm])*normalizer
    A_pr <- c(A_prior[,mm],At_prior[,mm])
    V_pr <- c(theta[,mm],omega_pr[,mm])
    draw <- get.A(y=Y_i,x=X_i,A_pr=A_pr,V_pr=V_pr)
    
    A_draw[,mm] <- draw[1:(k+1)]
    sqrtomega <- draw[(k+1+1):(2*(k+1))]
    omega_mat[,mm] <- sqrtomega
    
    for(tt in 1:T){
      Atvp_draw <- A_draw[,mm] + omega_mat[,mm]*At_draw[tt,,mm]
      PHI_draw_t[tt,,mm] <- crossprod(W_i,Atvp_draw)
      epsilon[tt,mm] <- Y[tt,mm] - X_iraw[tt,] %*% Atvp_draw
      eta[tt,mm] <- Y_[tt,mm] - X_iraw[tt,] %*% Atvp_draw
    }
  }

  # shrinkage
  A_mat <- A_draw
  O_mat <- omega_mat
  
  hs_draw <- get.hs(bdraw=as.numeric(A_mat)-as.numeric(A_prior),lambda.hs=as.numeric(A.lambda),nu.hs=as.numeric(A.nu),tau.hs=A.tau,zeta.hs=A.zeta)
  theta <- matrix(hs_draw$psi,k+1,m)
  A.lambda <- matrix(hs_draw$lambda,k+1,m)
  A.nu <- matrix(hs_draw$nu,k+1,m)
  A.tau <- hs_draw$tau
  A.zeta <- hs_draw$zeta
  
  hs_draw <- get.hs(bdraw=as.numeric(O_mat),lambda.hs=as.numeric(O.lambda),nu.hs=as.numeric(O.nu),tau.hs=O.tau,zeta.hs=O.zeta)
  omega_pr <- matrix(hs_draw$psi,k+1,m)
  O.lambda <- matrix(hs_draw$lambda,k+1,m)
  O.nu <- matrix(hs_draw$nu,k+1,m)
  O.tau <- hs_draw$tau
  O.zeta <- hs_draw$zeta
  
  # idiosyncratic stochastic volatilities
  for(mm in 1:m){
    svdraw_mm <- svsample_fast_cpp(eta[,mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
    sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
    sv_tmp <- svdraw_mm$latent
    sv_tmp[exp(sv_tmp)<1e-4] <- log(1e-4)
    Sig.t[,mm] <- sv_latent[[mm]] <- sv_tmp
  }

  # draw factors and loadings
  F_draw <- get.factors(e=epsilon,S=Sig.t,H=H,L=Lambda,q=q,t=T)
  Lambda <- get.Lambda(eps=epsilon,fac=F_draw,S=Sig.t,pr=theta.lambda,m=m,q=q)
  if(shrink.lambda){
    hs_draw <- get.hs(bdraw=as.numeric(Lambda[id.lambda]),lambda.hs=as.numeric(L.lambda[id.lambda]),nu.hs=as.numeric(L.nu[id.lambda]),tau.hs=L.tau,zeta.hs=L.zeta)
    theta.lambda[id.lambda] <- hs_draw$psi
    L.lambda[id.lambda] <- hs_draw$lambda
    L.nu[id.lambda] <- hs_draw$nu
    L.tau <- hs_draw$tau
    L.zeta <- hs_draw$zeta
  }
  
  epsilon.hat <- matrix(0,T,m)
  for (tt in 1:T){
    epsilon.hat[tt,] <- Y[tt,] - X[tt,-K] %*% PHI_draw_t[tt,-K,]
  }
  
  if(fsv.rw){
    sv.draw <- SV_JPR(mu_h0=0,sig_h0=10,c_h=0,rho_h=1,sig_h=fsv_para[3],Em=epsilon.hat,hv=log(H),L=Lambda,gamma=PHI_draw_t[,K,],Sigma=exp(Sig.t),vnr=m,mean.exp=svmean.exp,scale_h=1.2)
    H.acc <- H.acc+(H!=exp(sv.draw$hv)*1)
    H <- exp(sv.draw$hv)
    
    ht_eps <- diff(as.numeric(c(sv.draw$h0,sv.draw$hv)))
    fsv_para[3] <- get.svpara_rw(ht_eps,0.2,0.1)
  }else{
    sv.draw <- SV_JPR(mu_h0=0,sig_h0=10,c_h=fsv_para[1],rho_h=fsv_para[2],sig_h=fsv_para[3],Em=epsilon.hat,hv=log(H),L=Lambda,gamma=PHI_draw_t[,K,],Sigma=exp(Sig.t),vnr=m,mean.exp=svmean.exp,scale_h=0.1)
    H <- exp(sv.draw$hv)
    
    ht_y <- as.numeric(c(sv.draw$h0,sv.draw$hv))
    yy <- ht_y[2:length(ht_y)]
    xx <- cbind(1,ht_y[1:(length(ht_y)-1)])
    fsv_para <- get.svpara(y=yy,X=xx,b0=c(0,1),V0=diag(2),mm=0.2,vv=0.1)
  }
  
  if(svmean.exp){
    X[,K] <- H
  }else{
    X[,K] <- log(H)
  }

  # store quantities
  if (irep %in% thin.set){
    in.thin <- in.thin+1
    PHI_store[in.thin,,,] <- PHI_draw_t
    Sigma_store[in.thin,,] <- Sig.t
    H_store[in.thin,,] <- H
  }

  if(irep %in% toplot){
    ts.plot(ts(log(H),start=c(1995,1+p),frequency=12),
            ylab="Uncertainty",main=paste0("MCMC: ",irep),
            ylim=c(0,max(H)))
    # for(mm in 1:m){
    #   pdf(file=paste0("tmp_plot/",colnames(Y)[mm],".pdf"),width=5,height=5)
    #   ts.plot(ts(PHI_draw_t[,"UNC",colnames(Y)[mm]],start=c(1991,1+p),frequency=12),ylab=paste0(colnames(Y)[mm],", MCMC: ",irep))
    #   dev.off()
    # }
    message("\n H-acceptance rate: ",paste0(c("Low="," Med="," High="),round(100*quantile(H.acc/irep,probs=c(0.25,0.50,0.75)),digits=1)))
  }
  
  if(check.time == TRUE){
    setTxtProgressBar(pb, irep)
  }
}

# -------------------------------------------------------------------------------------------------------
# plot structural inference
colfunc <- colorRampPalette(c("red", "yellow"))
irf_store <- array(NA,dim=c(nthin,T,m,nhor))
dimnames(irf_store) <- list(NULL,as.character(dates),colnames(Y),paste0("hor",0:(nhor-1)))

for(irep in 1:nthin){
  PHI_t <- PHI_store[irep,,,]
  H_t <- H_store[irep,,]
  A_t <- PHI_t[,1:(K-(cons+1)),]
  imp_t <- PHI_t[,K,]
  
  irf <- array(NA,dim=c(T,m*p,nhor))
  for(tt in 1:T){
    AA <- get.companion(A_t[tt,,],c(m,0,p))$MM
    eigen.crit <- abs(max(Re(eigen(AA)$value)))
    if(eigen.crit>1.05) next
    irf[tt,,1] <- c(Ysd*imp_t[tt,],rep(0,m*(p-1)))
    for(hh in 2:nhor){
      irf[tt,,hh] <- AA %*% irf[tt,,hh-1]
    }
  }
  irf_store[irep,,,] <- irf[,1:m,]
}

irf_post <- apply(irf_store,c(2,3,4),quantile,probs=c(0.16,0.25,0.50,0.75,0.84),na.rm=TRUE)
irf_cum0 <- apply(apply(irf_store[,,,1],c(1,2,3),sum),c(2,3),quantile,probs=c(0.16,0.25,0.50,0.75,0.84),na.rm=TRUE)
irf_cum1 <- apply(apply(irf_store[,,,1:12],c(1,2,3),sum),c(2,3),quantile,probs=c(0.16,0.25,0.50,0.75,0.84),na.rm=TRUE)
irf_cum3 <- apply(apply(irf_store,c(1,2,3),sum),c(2,3),quantile,probs=c(0.16,0.25,0.50,0.75,0.84),na.rm=TRUE)

# save output file
dir.create("mcmc",showWarnings = FALSE)
ret.obj <- list("H_mcmc"=H_store,"irf_post"=irf_post,"irf_cum0"=irf_cum0,"irf_cum1"=irf_cum1,"irf_cum3"=irf_cum3)
save(ret.obj,file=paste0("mcmc/",pname,".rda"))

# plot of the uncertainty measure
par(mar=c(3,4,1,1))
pdf(file=paste0("results/unc-",pname,".pdf"),width=10,height=6)
ts.plot(ts(cbind(0,t(apply(sqrt(H_store),c(2),quantile,probs=c(0.16,0.50,0.84)))),start=c(1995,1+p),frequency=12),
        ylab=paste0("uncertainty"),col=c("red",rep("black",3)),lwd=c(2,1,2,1))
dev.off()

for(cc in 1:length(cN)){
  sl.vars <- grepl(cN[cc],colnames(Y))
  sl.labs <- colnames(Y)[sl.vars]
  tmp_post <- irf_post[,,sl.vars,]
  tmp_cum <- irf_cum[,,sl.vars]
  
  # posterior median response over time
  pdf(file=paste0("results/irft-",pname,"_",cN[cc],".pdf"),width=10,height=6)
  par(mfrow=c(ceiling(sqrt(length(sl.labs))),ceiling(sqrt(length(sl.labs)))),mar=c(3,4,1,1))
  for(j in 1:length(sl.labs)){
    ts.plot(cbind(t(tmp_post[2,,j,]),0),ylab=sl.labs[j],col=c(colfunc(T),"black"))
  }
  dev.off()
  
  pdf(file=paste0("results/irfc-",pname,"_",cN[cc],".pdf"),width=10,height=6)
  par(mfrow=c(ceiling(sqrt(length(sl.labs))),ceiling(sqrt(length(sl.labs)))),mar=c(3,4,1,1))
  for(j in 1:length(sl.labs)){
    ts.plot(ts(cbind(0,t(tmp_cum[,,j])),start=c(1995,1+p),frequency=12),col=c("red",rep("black",3)),lwd=c(2,1,2,1),ylab=sl.labs[j])
  }
  dev.off()
}


