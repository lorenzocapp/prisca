#' solocp_single compute the marginal inclusion probabilities for all time instance in the case n_t=1
#'
#' @param y vector of observations
#' @param sigma2 sample variance
#' @param q sparsity inducing parameters (if not included set to 0.1)
#' @param tau2 shrinkage parameter (if not included set to 2/sqrt(n))
#' @param tau2.spike spike variance parameters (if not included set to 1/n)
#' @param tau2.slab slab variance parameters (if not included set to n)

#'
#' @return Weight of the spike compenent  \code{w.spike}
#'  Weight of the slab compenent  \code{w.slab}
#'   inclusion probability \code{ratio}
#' @export


solocp_single<-function(y,sigma,q=0.1,tau2=NULL,tau2.spike=NULL,tau2.slab=NULL){

  sigma2 <- sigma^2
  n<-length(y)
  grid <- seq(1,n)
  n.grid <- length(grid)-1
  n1<-rep(1,n.grid+1)

  if (is.null(tau2)){ tau2=2/sqrt(n)}
  if (is.null(tau2.spike)){tau2.spike=1/n}
  if (is.null(tau2.slab)){tau2.slab=n}

  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)

  #Fill the weight matrix
  lower<-tau2/(tau2+sigma2) #this is going down
  lower.local.means<-y[n.grid+1] #this is going down
  parsumLOW <- lower
  for (j1 in (n-1):1){
    new.lower<-tau2*((n-j1+1)-parsumLOW)^2/(tau2*((n-j1+1)-(parsumLOW))+sigma2)
    lower<-c(new.lower,lower)
    parsumLOW <- parsumLOW + new.lower
    # if (abs(new.lower-1)<.000001){
    #   lower <- c(rep(1,n-length(lower)),lower)
    #   break}
  }
  revlower <- revcumsum(lower)
  revy <- revcumsum(y)
  lower.local.means<-y[n.grid+1] #this is going down
  parsum <- lower.local.means*lower[n]
  for (j1 in (n-1):1){
    new.lower.mean <- (revy[j1]-parsum)/((n-j1+1)-revlower[j1+1])
    lower.local.means<-c(new.lower.mean,lower.local.means)
    parsum <- parsum + lower[j1]*new.lower.mean
  }

  sum.lower<-revcumsum(lower)
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-revcumsum(mean.disc)
  marg.mean <- c(sum.mean.disc[2:length(sum.mean.disc)],0)


  marg.weight <- c(sum.lower[2:length(sum.lower)],0)
  #Define M and GAMMA and Y matrices
  sum.inv<-seq(n,1)

  GAM<-matrix(NA,nrow=n,ncol=n)
  M<-matrix(NA,nrow=n,ncol=n)
  Y<-matrix(NA,nrow=n,ncol=n)

  sum.y<-rev(cumsum(rev(y)))

  #param
  GAM[,1] <- 1
  Y[,1] <- sum.y[1]-marg.mean
  M[,1]<-1/((sum.inv[1] - marg.weight)*GAM[,1]+sigma2*tau2^(-1))
  GAM[2:n,2]<-1-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]
  Y[2:n,2]<-sum.y[2]-marg.mean[2:n]-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]*GAM[2:n,1]*Y[2:n,1]
  #weight
  w.spike<-sqrt(tau2.spike^(-1)/(lower[1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.spike^(-1)))#maybe wrong but it does not matter
  w.slab<-sqrt(tau2.slab^(-1)/(lower[1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.slab^(-1)))#maybe wrong but it does not matter
  parsumMGAM2 <- M[1:n,1]*GAM[1:n,1]^2
  parsumMGAMY <- M[1:n,1]*GAM[1:n,1]*Y[1:n,1]
  for (i in 3:(n-1)){
    #param
    M[(i-1):n,(i-1)]<-1/((sum.inv[(i-1)] - marg.weight[(i-1):n])*GAM[(i-1):n,(i-1)]+sigma2*tau2^(-1))
    parsumMGAM2 <- parsumMGAM2 +  M[1:n,(i-1)]*GAM[1:n,(i-1)]^2
    GAM[i:n,i]<-1-(sum.inv[i] - marg.weight[i:n])*parsumMGAM2[i:n]
    parsumMGAMY <- parsumMGAMY + M[1:n,(i-1)]*GAM[1:n,(i-1)]*Y[1:n,(i-1)]
    Y[i:n,i]<-sum.y[i]-marg.mean[i:n]-(sum.inv[i] - marg.weight[i:n])*parsumMGAMY[i:n]
    #weight
  }

  #param
  M[(n-1):n,(n-1)]<-1/((sum.inv[(n-1)] - marg.weight[(n-1):n])*GAM[(n-1):n,(n-1)]+sigma2*tau2^(-1))
  GAM[n,n]<-1-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]^2)
  Y[n,n]<-sum.y[i]-marg.mean[n]-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]*Y[n,1:(n-1)])

  idGAM <- matrix(c(seq(2,n),seq(1,n-1)),ncol=2)
  w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.spike^(-1)))*exp(1/(2*sigma2)*diag(Y)[-1]^2/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.spike^(-1))))
  w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.slab^(-1)))*exp(1/(2*sigma2)*diag(Y)[-1]^2/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.slab^(-1))))


  ratio<-q*w.slab/(q*w.slab+(1-q)*w.spike)
  id.inf <- which(w.slab==Inf)
  ratio[id.inf] <- 1


  return(list(w.slab=w.slab,
              w.spike=w.spike,
              ratio=ratio))


}




revcumsum <- function(x){
  return(rev(cumsum(rev(x))))
}


