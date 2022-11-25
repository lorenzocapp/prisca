#' Prisca Algorithm 1 implementation
#'
#' @param y vector of observations
#' @param L number of single effects (upper bound on changes in variances)
#' @param a0 hyperparameter on the Gamma distribution ()
#' @param eps stopping rule of the ELBO maximization (default is .01)
#' @param limit maximum number of iterations for ELBO convengece (default is 10000)
#' @param grid default is FALSE, which corresponds to one observations per time point. If grid=TRUE, it bins the data into a number of intervals specified by n.grid (see below)
#' @param AR default is FALSE. If AR=TRUE, it fits AR-PRISCA with r=1 (see paper for details)
#' @param tf default is FALSE, If AR=TRUE, it fits TF-PRISCA with using the trend filter algorithm (see paper for details)
#' @param x default is NULL. When grid=TRUE, one can inputs at which time instance is y has been collected
#' @param n.grid default is NULL. When grid=TRUE, one can inputs how many intervals to define
#' @param b0 default is equal to a0. It is the rate parameter of the Gaussian

#'
#' @return  \code{ELBO}: vector of all ELBO values, 
#'           \code{alpha.mat}: It is a L x T matrix. Each row is a posterior of an effect \gamma_l (\alpha_l)
#'           \code{phihat}: If AR=TRUE, it outputs the AR coefficient
#'           \code{meanhat}: If tf=TRUE, it outputs the tf estimated mean
#' @export


prisca <- function(y,L,a0,eps=.01,limit=10000,
                   grid=FALSE,AR=FALSE,tf=FALSE,
                   x=NULL,
                   n.grid=NULL,
                   b0=NULL,
                   ...){
  
  "Function to do variance change point detection
  It outputs:
  -ELBO: vector of all ELBO values
  -alpha.mat: T x L matrix with probabilities of change point locations
  -if AR=TRUE, it outputs also the estimated AR coefficient"
  
  if (grid==FALSE){
    "This is the case where you have one observations per time point"
    T <- length(y)
    Sy2 <- y^2
    Sn <- revcumsum(rep(1,T))
  } else if (grid==TRUE){
    sample <- data.frame(x,y)
    init<-initialize(sample,n.grid)
    n <- init$n.vec
    Sy2 <- init$Sy2
    Sn <- revcumsum(n)
    T <- length(n)
  }
  
  if (length(b0)==0){b0=a0}
  #Initialization
  a.mat <- matrix(rep(a0,T*L),ncol=T) #posterior a matrix never changes
  b.mat <- matrix(rep(b0,T*L),ncol=T) #where I am going to store posterior b
  tau2.mat <- a.mat/b.mat
  expec.mat <- tau2.mat #expected precision
  alpha.mat <- matrix(rep(0,T*L),ncol=T) #This is the prior on the weights
  pi <- rep(1/T,T)
  dELBO.s <- c(-10000)
  A <- rep(0,L)
  C <- rep(0,L)
  #Posterior computation
  a.mat <- matrix(rep(a0+Sn/2,L),ncol=T,byrow = T) 
  dig.a <- digamma(a0+Sn/2) #trigamma function of vector a posterior
  a <- a0+Sn/2 # vector of a posterior
  loggamma.a <- lgamma(a) #log \Gamma (a)
  
  
  
  i <- 1
  delta <- 1
  while (abs(delta)>eps){
    
    if (AR==TRUE){
      "AR part. We deal with unknown phi, since, known phi it is enough a transformation"
      expecprec <- matrixStats::colProds(expec.mat)[-1]
      phihat <- lm (y[2:T]~ -1 +y[1:(T-1)],weights=expecprec)$coefficients #we multiply and not 1/w since it is the precision
      resy <- c(y[1],y[2:T]-phihat*y[1:(T-1)])
      Sy2 <- resy^2
    } else if (tf==TRUE){
      expecprec <- matrixStats::colProds(expec.mat)
      mean <- trend_filter_bic(y,k,prec=expecprec)
      resy <- y-mean
      Sy2 <- resy^2
    }
    "Loops for the variance starts below"
    for (l in 1:L){
      
      
      #Think about an alternative to the one below. Overall I am doing many more calculations here. 
      resSy2 <- scale_effect(tau2.mat[-l,,drop=FALSE],alpha.mat[-l,,drop=FALSE],Sy2,L-1) #This gotta be such I am just removing one, not all?
      #vec_log_prob.for <- vec_log_marginal_prob(a0,b0,resSy2,Sn)
      vec_log_prob.for <- logmarginal(a0,b0,resSy2,Sn)
      #plot(vec_log_prob.for)
      nc.oneway <- matrixStats::logSumExp(vec_log_prob.for)
      incl.prob <- exp(vec_log_prob.for-nc.oneway)
      alpha.mat[l,] <- incl.prob
      alpha.mat[l,alpha.mat[l,]<1e-300] <- 1e-300
      # alpha.mat[l,] <- posterior_alpha(a,b,resSy2,Sn)
      b.mat[l,] <- b0 + revcumsum(resSy2)/2
      tau2.mat[l,] <- a.mat[l,]/b.mat[l,]
      #tau2.mat[l,]<- posterior_precision_oneway(a,b,resSy2,n.plus.for)
      
      #ELBO
      #Term A
      A.l <- sum(alpha.mat[l,]*(log(pi/alpha.mat[l,])-(a*log(b.mat[l,])-loggamma.a)+
                                  (a0-a)*(dig.a-log(b.mat[l,]))-(b0-b.mat[l,])*tau2.mat[l,]))
      A[l]<-A.l
      #Term B
      expec.mat[l,] <- cumsum(alpha.mat[l,]*tau2.mat[l,])+(1-cumsum(alpha.mat[l,]))
      #Term C
      C.l <- sum(alpha.mat[l,]*Sn/2*(dig.a-log(b.mat[l,])))
      C[l] <- C.l
    }
    #Complete ELBO calculations and compute the change  
    B <- sum(matrixStats::colProds(expec.mat)*Sy2/2)
    ELBO <- sum(A) -B +sum(C)
    dELBO.s <- c(dELBO.s,ELBO)
    delta <- ELBO-dELBO.s[i]
    i <- i+1 
    if (i>limit){break}
  }
  
  #Output
  out <- list(ELBO=dELBO.s,alpha.mat=alpha.mat,y=y)
  if (AR==TRUE){out$phihat <- phihat}
  if (tf==TRUE){out$meanhat <- mean}
  
  return(out=out)
}