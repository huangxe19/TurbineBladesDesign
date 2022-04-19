library(tidyverse)
library(dplyr)
library(plgp)
library(lhs)
library(laGP)
library(maximin)
library(nloptr)


# generate maximin design
mymaximin <- function(n, m, T=100000, Xorig=NULL) {   
  ## return a maximin design given dimension and size
  X <- matrix(runif(n*m), ncol=m)     ## initial design
  d <- plgp::distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  if(!is.null(Xorig)) {               ## new code
    md2 <- min(distance(X, Xorig))
    if(md2 < md) md <- md2
  }
  
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,]                   ## random row selection
    X[row,] <- runif(m)               ## random new row
    d <- plgp::distance(X)
    d <- d[upper.tri(d)]
    mdprime <- min(d)
    if(!is.null(Xorig)) {             ## new code
      mdprime2 <- min(plgp::distance(X, Xorig))
      if(mdprime2 < mdprime) mdprime <- mdprime2
    }
    if(mdprime > md) { md <- mdprime  ## accept
    } else { X[row,] <- xold }        ## reject
  }
  
  return(X)
}

# generate beta-distribution design
BetaDistD <- function(n, m, a = 2,b = 5, T=100000) {  
  ## return a beta-distribution design (Beta(2,5)) given dimension and size
  X <- matrix(runif(n*m), ncol=m)    
  D_init <- distance(X)
  D_init <- D_init[upper.tri(D_init)]
  KSD <- ks.test(D_init,"pbeta",a,b)$statistic
  
  for(t in 1:T) {
    row <- sample(1:n, 1)             
    Xold <- X                         
    X[row,] <- runif(m)             
    D_s <- distance(X)
    D_s <- D_s[upper.tri(D_s)]
    KSDprime <- ks.test(D_s,"pbeta",a,b)$statistic
    if(KSDprime < KSD) { KSD <- KSDprime 
    } else { X <- Xold }              
  }
  
  return(X)
}


# convert design to actual scale
convert_scale <- function(OldMax=1, OldMin=0, OldValue, NewMax, NewMin){
  ## convert a vector into given scale
  OldRange = (OldMax - OldMin)  
  NewRange = (NewMax - NewMin)  
  NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
}
convert_vars <- function(oldmatrix){
  ## convert the six variables into actual scale
  young <- mapply(convert_scale, OldValue = oldmatrix[,1], NewMax = 300e9, NewMin=200e9)
  poisson <- mapply(convert_scale, OldValue = oldmatrix[,2], NewMax = 0.1, NewMin=0.49)
  CTE <- mapply(convert_scale, OldValue = oldmatrix[,3], NewMax = 5e-6, NewMin=1.5e-5)
  thermal <- mapply(convert_scale, OldValue = oldmatrix[,4], NewMax = 5, NewMin=15)
  temp <- mapply(convert_scale, OldValue = oldmatrix[,5], NewMax = 50, NewMin=360)
  pressure <- mapply(convert_scale, OldValue = oldmatrix[,6], NewMax = 1e5, NewMin=4.8e5)
  newmatrix <- matrix(c(young, poisson, CTE, thermal, temp, pressure),
                      byrow=F, ncol=6)
  return(newmatrix)
}


# interface simulator between matlab and r
simulator <- function(input){
  ## return output of simulator.p given a vector of input
  ymod = input[1]; prat = input[2]; cte = input[3]; 
  therm = input[4]; ctemp = input[5]; press = input[6]
  
  # set a variable in R and send to MATLB
  setVariable(matlab,
              ymod = ymod,
              prat = prat,
              cte = cte,
              therm = therm,
              ctemp = ctemp,
              press = press)
  
  evaluate(matlab, "[stress,displ] = simulator(ymod,prat,cte,therm,ctemp,press)")
  
  return(c(getVariable(matlab, "stress")$stress,getVariable(matlab, "displ")$displ))
}


# Set optimization options
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 160000,
              "local_opts" = local_opts,
              "print_level" = 0 )

# objective function
eval_f <- function(x, gpi_stress, gpi_disp, Xref,n){
  return(- sqrt(alcGPsep(gpi_stress, matrix(x, nrow=1), Xref)))
}

# constraint
eval_g <- function(x, gpi_stress, gpi_disp, Xref, n) {
  p_d <- predGPsep(gpi_disp, matrix(x, nrow=1), lite=TRUE)
  return(1.3e-3 - (p_d$mean-n*sqrt(p_d$s2)))
}

# optimize ALC objective under constraint
optim_alc <-function(x_start, eval_f, eval_g, opts, gpi_stress, gpi_disp, xref, n) {
  ## return new point that optimized ALC objective given constraint
  ## take optimization options, start point, objective, constraint, gp models, and reference set as input
  xnew <- nloptr(x0=x_start,
                 eval_f=eval_f,
                 lb = c(200e9, 0.1, 5e-6, 5, 50, 1e5),
                 ub = c(300e9, 0.49, 1.5e-5, 15, 350, 4.8e5),
                 eval_g_ineq = eval_g,
                 opts = opts,
                 gpi_stress = gpi_stress,
                 gpi_disp =gpi_disp,
                 Xref = xref,
                 n = n)
  return(xnew$solution)
}


# ALC
ALC <- function(X_raw, y, XX_raw, yy, niter) {
  ## take train and test set(in range(0,1)) and number of sequential runs as input
  ## return final design, y, and rmse
  
  num_para <- ncol(X_raw)
  X <- convert_vars(X_raw)
  XX <- convert_vars(XX_raw)
  
  names(X) <- c("x1", "x2","x3", "x4", "x5", "x6")
  names(XX) <- c("x1", "x2","x3", "x4", "x5", "x6")
  
  y_s <- y[,1]
  y_d <- y[,2]
  yy_s <- yy[,1]
  
  # GP for stress
  g <- garg(list(mle=TRUE), y_s)
  d <- darg(list(mle=TRUE), X)
  gpi <- newGPsep(X, y_s, d=rep(d$start,6), g=1/1000, dK=TRUE)
  mle <- jmleGPsep(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)
  p <- predGPsep(gpi, XX, lite=TRUE)
  rmse.alc <- sqrt(mean((yy_s - p$mean)^2))
  
  # GP for displacement
  g_disp <- garg(list(mle=TRUE,min=0), y_d)
  d_disp <- darg(list(mle=TRUE), X_raw)
  gpi_disp <- newGPsep(X_raw, y_d, d=rep(d_disp$start, 6), g=1/1000, dK=TRUE)
  mle_disp <- jmleGPsep(gpi_disp, c(d_disp$min, d_disp$max), c(g_disp$min, g_disp$max), d_disp$ab, g_disp$ab)
  
  # search for Xnew
  Xref <- convert_vars(randomLHS(100, num_para))
  x0_raw <- matrix(runif(6),ncol=6)
  x0 <- convert_vars(x0_raw) %>% as.vector()
  xnew <- optim_alc(x_start    =x0, 
                    eval_f     =eval_f, 
                    eval_g     =eval_g, 
                    opts       =opts, 
                    gpi_stress =gpi, 
                    gpi_disp   =gpi_disp,
                    xref       =Xref,
                    n          =1.96)
  
  xnew <-matrix(xnew, nrow=1,byrow=T)
  X <- rbind(X, xnew)
  y_out <- simulator(xnew)
  y_s <- c(y_s, y_out[1])
  y_d <- c(y_d, y_out[2])
  
  # update GP for stress
  updateGPsep(gpi, xnew, y_s[length(y_s)])
  mle <- rbind(mle, jmleGPsep(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab))
  p <- predGPsep(gpi, XX, lite=TRUE)
  rmse.alc <- c(rmse.alc, sqrt(mean((yy_s - p$mean)^2)))
  
  # update GP for displacement
  updateGPsep(gpi_disp, xnew, y_d[length(y_d)])
  mle_disp <- rbind(mle_disp, 
                    jmleGPsep(gpi_disp, c(d_disp$min, d_disp$max), c(g_disp$min, g_disp$max), d_disp$ab, g_disp$ab))
  p_disp <- predGPsep(gpi_disp, XX, lite=TRUE)
  
  d <- darg(list(mle=TRUE), X)
  d_disp <- darg(list(mle=TRUE), X_raw)
  

  for(i in 1:(niter-1)) {
    # search for Xnew
    Xref <- convert_vars(randomLHS(100, num_para))
    x0_raw <- matrix(runif(6),ncol=6)
    x0 <- convert_vars(x0_raw) %>% as.vector()
    print(x0)
    
    if (i<=10) {n <- 1.96} 
      else if (10 < i & i>= 20) {n <- 1.64}
      else {n <- 0}
    
    xnew <- optim_alc(x_start    =x0, 
                      eval_f     =eval_f, 
                      eval_g     =eval_g, 
                      opts       =opts, 
                      gpi_stress =gpi, 
                      gpi_disp   =gpi_disp,
                      xref       =Xref,
                      n = n)
    
    xnew <-matrix(xnew, nrow=1,byrow=T)
    X <- rbind(X, xnew)
    y_out <- simulator(xnew)
    y_s <- c(y_s, y_out[1])
    y_d <- c(y_d, y_out[2])
    
    # update GP for stress
    updateGPsep(gpi, xnew, y_s[length(y_s)])
    mle <- rbind(mle, jmleGPsep(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab))
    p <- predGPsep(gpi, XX, lite=TRUE)
    rmse.alc <- c(rmse.alc, sqrt(mean((yy_s - p$mean)^2)))
    
    # update GP for displacement
    updateGPsep(gpi_disp, xnew, y_d[length(y_d)])
    mle_disp <- rbind(mle_disp, 
                      jmleGPsep(gpi_disp, c(d_disp$min, d_disp$max), c(g_disp$min, g_disp$max), d_disp$ab, g_disp$ab))
    p_disp <- predGPsep(gpi_disp, XX, lite=TRUE)
  }
  
  return(list(X = X, y = cbind(y_s, y_d), rmse.alc = rmse.alc, mle = mle, gpi = gpi))
}

# EI function
EI <- function(gpi, x, fmin, pred=predGPsep) {
  # returns EI value
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  p <- pred(gpi, x, lite=TRUE)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d/sigma
  ei <- d*pnorm(dn) + sigma*dnorm(dn)
  return(ei)
}

# EI objective
obj.EI <- function(x, fmin, gpi, pred=predGPsep)
  - EI(gpi, x, fmin, pred)


eps <- sqrt(.Machine$double.eps) 

# EI search
EI.search <- function(X, y, gpi, pred=predGPsep, multi.start=5, tol=eps) {
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(X[m,], nrow=1)
  if(multi.start > 1)
    start <- rbind(start, convert_vars(randomLHS(multi.start - 1, ncol(X))))
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    if(EI(gpi, start[i,], fmin) <= tol) { out <- list(value=-Inf); next }
    out <- optim(start[i,], obj.EI, method="L-BFGS-B",
                 lower=c(200E9,0.1 ,5E-6 ,5  ,50 ,1E5), 
                 upper=c(300E9,0.49,1.5E-5,15,350,4.8E5), 
                 gpi=gpi, pred=pred, fmin=fmin)
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "s3", "s4","s5", "s6",
                    "x1", "x2", "x3", "x4","x5", "x6","val")
  solns <- solns[solns$val > tol,]
  return(solns)
}

# wrapper for EI optimization
optim.EI <- function(f, X, y, end) {
  ## take function to optimize, initial design X and y, number of iterations as inputs
  
  d <- darg(list(mle=TRUE), X)
  gpi <- newGPsep(X, y, d=rep(d$start, 6), g=1e-6, dK=TRUE)
  mle <- mleGPsep(gpi, param='d',tmin= d$min, tmax=d$max, ab=d$ab)

  ## optimization loop of sequential acquisitions
  maxei <- c()
  for(i in 1:(end-1)) {
    solns <- EI.search(X, y, gpi)
    m <- which.max(solns$val)
    maxei <- c(maxei, solns$val[m])
    xnew <- as.matrix(solns[m,7:12])%>% unlist()
    ynew <- f(xnew)[2]
    updateGPsep(gpi, matrix(xnew,nrow = 1), ynew)
    mle <- rbind(mle, mle <- mleGPsep(gpi, param='d',tmin= d$min, tmax=d$max, ab=d$ab))
    X <- rbind(X, xnew)
    y <- c(y, ynew)
  }
  ## clean up and return
  deleteGPsep(gpi)
  return(list(X=X, y=y, maxei=maxei))
}