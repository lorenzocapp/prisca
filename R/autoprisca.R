#' Autoprisca: identical to prisca Algorithm 1 except that L is determined through a heuristic (see paper for details)
#'
#' @param y vector of observations
#' @param a0 hyperparameter on the Gamma distribution ()
#' @param eps stopping rule of the ELBO maximization (default is .01)
#' @param grid default is FALSE, which corresponds to one observations per time point. If grid=TRUE, it bins the data into a number of intervals specified by n.grid (see below)
#' @param x default is NULL. When grid=TRUE, one can inputs at which time instance is y has been collected
#' @param n.grid default is NULL. When grid=TRUE, one can inputs how many intervals to define
#' @param b0 default is equal to a0. It is the rate parameter of the Gaussian

#'
#' @return  \code{ELBO}: vector of all ELBO values, 
#'           \code{alpha.mat}: It is a L X T matrix. Each row is a posterior of an effect \gamma_l (\alpha_l)
#' @export


autoprisca <- function(y,a0,eps=.01,grid=FALSE,x=NULL,
                       n.grid=NULL, b0=NULL,...){
  
  "Function to do variance change point detection
  This is an automated version, where you don't need to specify L. It will be probably be merged 
  to the main function.
  - It still misses the opportunity to do AR detection
  
  It outputs:
  -ELBO: vector of all ELBO values
  -alpha.mat: T x L matrix with probabilities of change point locations"
  
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
  #Initialization autoPRISCA
  L <- 2
  detectedA <- 1 #This is just in the initialization
  stableL <- FALSE
  
  while (stableL==FALSE){
    
    #Initialization
    if (L==2){#meaning it is the first cycle
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
    } else { #add one line
      a.mat <- rbind(rep(a0,T),a.mat)
      b.mat <- rbind(rep(b0,T),b.mat)
      tau2.mat <- rbind(rep(a0,T)/rep(b0,T),tau2.mat)
      alpha.mat <- rbind(rep(0,T),alpha.mat)
      expec.mat <- rbind(tau2.mat[1,],expec.mat) 
      a.mat[1,] <- a0+Sn/2
    }
    
    
    i <- 1
    delta <- 1
    while (abs(delta)>eps){
      #  #Scale the y's (afterwards I need to unscale them)
      #  #Wrong below! It is not how you scale them!! In order to have a AR
      # res <- scale_variance(tau2.mat,alpha.mat,y)
      # phi.hat <- ar(res,order=1)$ar
      # phi.hat <- 0.4
      # #Here it can be pretty much everything. I guess that one would need to scale 
      # #the y's to remove the scaling component
      # Sy2 <- c(y[1]^2,(y[2:T]-phi.hat*y[1:(T-1)])^2)
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
    }
    
    #This part is the autoupdate (L=3 is the minimum number for which it can stop)
    cs <- credible_set(alpha.mat,p=0.9,delta=3)
    if (L==2){
      detectedB <- length(cs) #For L=2 we complete the logical condition
      L <- L+1
    } else {
      detectedC <- length(cs)
      if (detectedA==detectedB & detectedA==detectedC & detectedB==detectedC){
        stableL <- TRUE
      } else {
        detectedA <- detectedB
        detectedB <- detectedC
        L <- L+1
      }
    }
    
  }
  return(list(ELBO=dELBO.s,alpha.mat=alpha.mat,y=y))
}