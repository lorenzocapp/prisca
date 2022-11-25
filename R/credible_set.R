#' Compute all the credible sets for a given matrix alpha of posterior weights (L x T matrix)
#'
#' @param alpha.mat  Output of prisca or autoprisca. It is a T x L matrix. Each row is a posterior of an effect \gamma_l (\alpha_l)
#' @param p Probability level to build credible sets. Default is 0.9 if not specified
#' @param rm_overlap Default is true. Check if there duplicated sets and remove them
#' @param usedelta Default is FALSE. If TRUE, it tries to construct set such that the minimum distance between two consecutive times is delta
#' @param delta Minimum distance to be used as input if usedelta = TRUE

#'
#' @return  \code{cs_list}:A list of credible sets. For each set, we provide: 1) point estimated 2) the set 3) promised coverage of the set
#' @export


credible_set <- function(alpha.mat,p=0.9,rm_overlap=TRUE,usedelta=FALSE,delta=3){
  
  #Compute all the credible sets for a given matrix alpha of posterior weights (L x T matrix)
  #It invokes the function build_cs for every vector.
  #alpha.vec - vector of posterior alpha_t
  #alpha - level of the credible set
  #delta - minimum spacing condition
  #lower - it is just an index using to build the recursion
  
  L <- dim(alpha.mat)[1]
  T <- dim(alpha.mat)[2]
  cs_list <- list()
  
  for (l in 1:L){
    alpha.vec <- alpha.mat[l,]
    if (usedelta==TRUE){
      cs_list[[l]]<- build_cs(alpha.mat[l,],p,delta,lower=0)
      cs_list[[l]]$index <- l
    } else {
      cs_list[[l]]<- nodelta_cs(alpha.mat[l,],p)
      cs_list[[l]]$index <- l     
    }
    
  }
  
  #Remove the "not available"
  idNA <- which(lapply(cs_list,function(y) (length(y$cs)==1 && y$cs=="Not available"))==TRUE)
  if(length(idNA)>0) cs_list <- cs_list[-idNA]
  #Remove sets that are too long. Arbitrary criteria: long means >1/4T
  idLong <- which(lapply(cs_list,function(y) length(y$cs)>0.25*T)==TRUE)
  if(length(idLong)>0) cs_list <- cs_list[-idLong]
  #Remove overlapping sets: 
  if (rm_overlap==TRUE) cs_list <- overlapping_sets(cs_list)
  
  return(cs_list=cs_list)
}