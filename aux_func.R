get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# -----------------------------------------------------------------------------------------------
# function to remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# -----------------------------------------------------------------------------------------------
# function to lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

# -----------------------------------------------------------------------------------------------
# function to create weights matrix
create.W_lag <- function(sl.C, G, N, p, cons, ex=NULL, weights=NULL){
  sl.ID <- matrix(seq(1,G*N),G,N)
  W_i <- matrix(0, 2*G, G*N)
  W_i[1:G,sl.ID[,sl.C]] <- diag(G)
  
  if (is.null(weights)){
    for (jj in 1:G){
      W_i[(G+jj), sl.ID[jj,-sl.C]] <- 1  
    }
  }else{
    for (jj in 1:G){
      W_i[(G+jj), sl.ID[jj,-sl.C]] <- weights 
    }
  }
  
  W_i <- t(apply(W_i, 1, function(x) x/sum(x)))
  if(p>1){
    W_i <- kronecker(diag(p),W_i)
  }
  if(cons){
    W_i <- rbind(cbind(W_i,0),0)
    W_i[nrow(W_i),ncol(W_i)] <- 1
  }
  if(!is.null(ex)){
    for(ii in 1:ex){
      W_i <- rbind(cbind(W_i,0),0)
      W_i[nrow(W_i),ncol(W_i)] <- 1
    }
  }
  return(W_i)
}

create.W_ctp <- function(sl.C, G, N, weights=NULL){
  sl.ID <- matrix(seq(1,G*N),G,N)
  W_c <- matrix(0, G, G*N)
  
  if (is.null(weights)){
    for (jj in 1:G){
      W_c[jj, sl.ID[jj,-sl.C]] <- 1  
    }
  }else{
    for (jj in 1:G){
      W_c[jj, sl.ID[jj,-sl.C]] <- weights 
    }
  }
  W_c <- t(apply(W_c, 1, function(x) x/sum(x)))
  return(W_c)
}

# -----------------------------------------------------------------------------------------------
# function to draw VAR coefficients
get.A <- function(y,x,A_pr,V_pr){
  V_post <- try(chol2inv(chol(crossprod(x)+diag(1/V_pr))),silent=TRUE)
  if (is(V_post,"try-error")) V_post <- ginv(crossprod(x)+diag(1/V_pr))
  A_post <- V_post%*%(crossprod(x,y)+diag(1/V_pr)%*%A_pr)
  
  A_draw <- try(A_post+t(chol(V_post))%*%rnorm(ncol(x)),silent=TRUE)
  if (is(A_draw,"try-error")) A_draw <- mvrnorm(1,A_post,V_post)
  return(A_draw)
}

# -----------------------------------------------------------------------------------------------
# function to draw the factor loadings (basic linear regression)
get.facload <- function(yy,xx,l_sd){
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

get.Lambda <- function(eps,fac,S,pr,m,q){
  L <- matrix(0,m,q)
  for(jj in 1:m){
    if (jj<=q){
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj]-fac[,jj])/normalizer
      xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
      if (jj>1){
        l_sd <- pr[jj,1:(jj-1)]
        lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      }else{
        lambda0 <- 1
      }
      
      if (jj>1){
        L[jj,1:(jj-1)] <- lambda0
        L[jj,jj] <- 1
      }else if (jj==1){
        L[jj,jj] <- 1
      }
    }else{
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

# sample the latent factors
get.factors <- function(e,S,H,L,q,t){
  F_raw <- matrix(0,t,q)
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    F_temp <- try(fac.mean + t(chol(fac.Sigma)) %*% rnorm(q),silent=TRUE)
    if (is(F_temp,"try-error")) F_temp <- fac.mean + t(chol(fac.Sigma+diag(q)*1e-6)) %*% rnorm(q)
    F_raw[tt,] <- F_temp
  }
  return(F_raw)
}

# -----------------------------------------------------------------------------------------------
# get standard regression coefficients
get.reg <- function(y,x,sig,b0,V0){
  k <- ncol(x)

  Vpo <- solve(crossprod(x)/sig + diag(k)*(1/V0))
  bpo <- Vpo %*% (crossprod(x,y)/sig + (diag(k)*(1/V0)) %*% matrix(b0,k,1))
  
  bdraw <- bpo + t(chol(Vpo)) %*% rnorm(k)
  return(bdraw)
}

# -----------------------------------------------------------------------------------------------
# obtain companion matrix
get.companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  if(nd){
    MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+1)),c(matrix(0,1,(nn*nl)),1))
  }else{
    MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
  }
  
  return(list(MM=MM,Jm=Jm))
}

# -----------------------------------------------------------------------------------------------
# JPR-algorithm
SV_JPR <- function(mu_h0,sig_h0,c_h=0,rho_h,sig_h,Em,hv,L,gamma,Sigma,vnr,scale_h,mean.exp=TRUE){
  require(mvtnorm)
  LL <- tcrossprod(L)
  for (it in 0:T){
    if (it==0){
      ht_sig <- 1/(mu_h0/sig_h0 + rho_h^2/sig_h)
      ht_mu <- ht_sig*(mu_h0/sig_h0 + rho_h*(hv[1]-c_h)/sig_h)
      h0d <- ht_mu + sqrt(ht_sig) * rnorm(1)
    }else{
      # define prior mean and variance t-by-t
      if (it==1){
        h_mut <- ((1-rho_h)*c_h + rho_h*(h0d + hv[2]))/(1+rho_h^2)
        h_sig <- sig_h/(1+rho_h^2)
      }else if (it==T){
        h_mut <-c_h + rho_h*hv[T-1]
        h_sig <- sig_h
      }else{
        h_mut <- ((1-rho_h)*c_h + rho_h*(hv[it-1] + hv[it+1]))/(1+rho_h^2)
        h_sig <- sig_h/(1+rho_h^2)
      }
      
      h_old <- hv[it]
      h_prop <- h_old + rnorm(1,0,sqrt(scale_h))
      if(mean.exp){
        e.0 <- Em[it,]-gamma[it,]*exp(h_prop)
        e.1 <- Em[it,]-gamma[it,]*exp(h_old)
      }else{
        e.0 <- Em[it,]-gamma[it,]*h_prop
        e.1 <- Em[it,]-gamma[it,]*h_old
      }
      v.0 <- exp(h_prop)*LL+diag(Sigma[it,])
      v.1 <- exp(h_old)*LL+diag(Sigma[it,])
      
      alphanum <- mvnfast::dmvn(e.0,mu=rep(0,vnr),sigma=v.0,log=TRUE) + dnorm(h_prop,h_mut,sqrt(h_sig),log=TRUE)
      alphaden <- mvnfast::dmvn(e.1,mu=rep(0,vnr),sigma=v.1,log=TRUE) + dnorm(h_old,h_mut,sqrt(h_sig),log=TRUE)
      u <- log(runif(1,0,1))<(alphanum-alphaden)
      hv[it] <- u*h_prop+(1-u)*h_old
    }
  }
  return(list(h0=h0d,hv=hv))
}

# draw conjugate regression coefs
get.svpara <- function(y,X,b0,V0,mm,vv){
  # prior specified in terms of mean "mm" and var "vv" of IG(a0,a1)
  a0 <- (2+((mm^2)/vv))
  a1 <- (mm + ((mm^3)/vv))
  
  n <- length(y)
  k <- ncol(X)
  
  par1 <- (a0+n)/2
  var <- solve(crossprod(X)+solve(V0))
  mean <- matrix(var%*%(crossprod(X,y)+crossprod(solve(V0),b0)),k,1)
  par2 <- a0*a1 + sum((y-crossprod(t(X),mean))^2)
  par2 <- (par2 + crossprod(t(crossprod(mean-b0,V0)),mean-b0))/2
  
  sig2 <- 1/rgamma(1,par1,par2)
  var <- var*sig2
  mean <- mean + crossprod(t(chol(var)),rnorm(k))
  return(c(mean,sig2))
}

get.svpara_rw <- function(eps,mm,vv){
  # prior specified in terms of mean "mm" and var "vv" of IG(a0,a1)
  a0 <- (2+((mm^2)/vv))
  a1 <- (mm + ((mm^3)/vv))
  
  n <- length(eps)
  return(1/rgamma(1,a0+n/2,a1+sum(eps^2)/2))
}
