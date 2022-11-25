
revcumsum <- function(x){
  return(rev(cumsum(rev(x))))
}


initialize <- function(sample,n.grid){
  #Defines a grid
  #Bin observations into the grid points
  #Compute local statistics within each interval that are needed to compute posteriors
  
  
  
  sample <- sample[order(sample$x),]
  n <- length(sample$x)
  grid.b<-seq(0,1,1/n.grid)
  grid <- grid.b[-n.grid-1]+(grid.b[2]-grid.b[1])/2
  #bracketing
  for (j in 1:n){
    i.low <- min(which(sample$x[j]<=grid.b))-1
    i.low <- max(1,i.low) #to include if there is a x=0
    sample$x[j] <- grid[i.low]
  }
  #binning
  count.grid <- table(sample$x) #how many obs per grid point
  x.observed <-names(count.grid)
  idn<-which(as.character(grid) %in% x.observed)
  n.vec<-rep(0,n.grid)
  n.vec[idn]<-count.grid
  sample$ysq <- sample$y^2
  S_ysq0 <-tapply(sample$ysq,sample$x, FUN=sum)
  Sy2<-rep(0,n.grid)
  Sy2[idn]<-S_ysq0
  return(list(Sy2=Sy2,n.vec=n.vec,
              sample.grid=sample,grid.bound=grid.b,grid.mid=grid))
}







logmarginal <- function(a,b,Sy2,Sn){
  "
  NEW version (computation boost)
  Compute the vector of marginal probability one way.
  Trying to reducee the computational burden of the initial function
  "
  T <- length(Sn)
  null <- cumsum(c(0,Sy2))[-(T+1)]/2
  alternative <- lgamma(a+Sn/2)-log(b+revcumsum(Sy2)/2)*(a+Sn/2)
  vec_logprob <- -null+alternative
  return(vec_logprob)
}

posterior_alpha <- function(a,b,resSy2,Sn,pi){
  
  if (missing(pi)) pi <- rep(1/length(Sn),length(Sn))
  #Compute posterior alpha (vector of inclusion probability)
  #Compute vector of marginal log lik (and include the prior)
  mloglik <- vec_log_marginal_prob(a,b,resSy2,Sn)#+log(pi)
  #Just do a few numerical tricks: removing the maximum and doing logSumExp
  r.mloglik <- mloglik #-max(mloglik)
  nc <- matrixStats::logSumExp(r.mloglik)
  alpha <- exp(r.mloglik -nc)
  return(alpha)
}


posterior_precision <- function(a,b,S_ysq,n.plus){
  #Compute posterior parameters
  a.post <- a +n.plus/2
  b.post <- b + revcumsum(S_ysq)/2
  return(prec.post=a.post/b.post)
}







scale_effect <- function(prec.post.matrix,incl.prob.matrix,S_ysq,L){

  if(L!=0){
  #Allows to scale the effect of multiple loops. (it needs to be always one less that the one I am considering)
  scale_vec <- c()
  scaled.S_ysq <- S_ysq
  for (l in 1:L){
    scale_single <- cumsum(prec.post.matrix[l,]*incl.prob.matrix[l,]) #here it is the scale contribution of each single point
    alternative_model <-1-cumsum(incl.prob.matrix[l,])
    scale <- scale_single+alternative_model
    scaled.S_ysq <- scaled.S_ysq*scale
  }
  return(scaled.S_ysq=scaled.S_ysq)
  } else {
    return(scaled.S_ysq=S_ysq)
  }
}













build_cs <- function(alpha.vec,alpha,delta,lower=0){
  
  #Compute the cs for a given vector alpha of posterior weights
  #It is a recursion, allowing for some sort of binary segmentation in order to build the sets
  #alpha.vec - vector of posterior alpha_t
  #alpha - level of the credible set
  #delta - minimum spacing condition
  #lower - it is just an index using to build the recursion
  
  alpha.sorted <- sort(alpha.vec,decreasing=TRUE)
  id.sorted <- order(alpha.vec,decreasing=TRUE)
  id.set <- which.min(cumsum(alpha.sorted)<=alpha)
  
  cs <- sort(id.sorted[1:id.set])
  check.cs <- diff(cs)>delta
  if (sum(check.cs)==0){
    return(list(point=c(cs+lower)[which.max(alpha.vec[cs])],
                cs=(cs+lower),prob=sum(alpha.vec[cs])))
  } else { 
    #If you are here, you have not managed to construct the set at first place. 
    #a. Identify where the points that have distance more than delta are
    id.low <- cs[which(check.cs==1)]
    #b. Compute the sum of weights within each segment (you need to understand if you can still have cs). NOTE: This part is arbitrary, you need to decide where to allocate points in the middle
    par.sums <- cumsum(alpha.vec)[c(id.low,length(alpha.vec))]-c(0,cumsum(alpha.vec)[id.low])
    #c. If you can attempt to construct another cs (there can be maximum 1 because we assume alpha>0.5), go ahead, otherwise empty
    if (any(par.sums>alpha)){
      id <- which(par.sums>alpha)
      #Identify the exact segment
      if (id==1){
        n.alpha.vec=alpha.vec[1:id.low[1]]
        lower=lower
      } else if (id==length(par.sums)){
        n.alpha.vec=alpha.vec[(id.low[id-1]+1):length(alpha.vec)]
        lower=lower+id.low[id-1]
      } else{
        n.alpha.vec=alpha.vec[(id.low[id-1]+1):id.low[id]]
        lower=lower+id.low[id-1]
      }
      #Rerun the function
      return(build_cs(n.alpha.vec,alpha,delta,lower))
    } else{ #It is impossible to construct a cs
      return(list(point="Not available",cs="Not available",prob=NA))
    }
  }
}


nodelta_cs <- function(alpha.vec,alpha){
  
  "
  Compute the cs for a given vector alpha of posterior weights
  Inputs are
  -alpha.vec - vector of posterior alpha_t
  -alpha - level of the credible set
  Difference with build_cs is that we are building sets as is, not using the delta
  "
  
  alpha.sorted <- sort(alpha.vec,decreasing=TRUE)
  id.sorted <- order(alpha.vec,decreasing=TRUE)
  id.set <- which.min(cumsum(alpha.sorted)<=alpha)
  cs <- sort(id.sorted[1:id.set])
  
  return(list(point=cs[which.max(alpha.vec[cs])],cs=cs,prob=sum(alpha.vec[cs])))
}




overlapping_sets <- function(cs_list){
  
    "Function to remove the overlapping credible sets
    Arbitrary criteria: if sets overlap more than 50% length of the shortest, 
    pick the shortest
    "
    K <- length(cs_list)
    i <- 1
    remove <- c()
    if(K>1){
    for (i in 1:(K-1)){
      for (j in (i+1):K){
      minlength <- min(length(cs_list[[i]]$cs),length(cs_list[[j]]$cs))
      length.inter <- length(intersect(cs_list[[i]]$cs,cs_list[[j]]$cs))
      if (length.inter>=0.5*minlength){
        id <- which.min(c(length(cs_list[[i]]$cs),length(cs_list[[j]]$cs)))
        if (id==1) remove <- c(remove,j) #if 1, it measn that i is the shortest
        if (id==2) remove <- c(remove,i) #if 2, it means that j is the shortest
      }
      }
    }
    idOverlap <- unique(remove)
    if(length(idOverlap)>0) cs_list <- cs_list[-idOverlap]
    }
    return(cs_list)
}




trend_filter_bic <- function(y,k=0,prec=1){

  "It is a weighted least square trend filtering. The inputs are 
  - data sequence
  - order of the tf
  - a vector with the precision
  Selection of lambda through BIC"
  
  path = glmgen::trendfilter(y,k=0,weights = prec)
  n <- length(y)
  sigma_hat2 = sum(diff(y*sqrt(prec))^2)/( (n-1)*2 )
  bic = rep(0,length(path$lambda))
  for(j in 1:length(path$lambda))
  {
    k = length(which(abs(diff( path$beta[,j] ))>10^-3))
    bic[j] = sum(((y - path$beta[,j])^2)*sqrt(prec)) +  k* log(n)* sigma_hat2
  }
  best_ind = which.min(bic)
  return(mean=path$beta[,best_ind])
}  
  
