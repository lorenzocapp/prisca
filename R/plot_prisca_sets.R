#' Plot the posterior probability along with the credible sets. It is very much inspired by the corresponding function of the Susie package (Wang et al. 20)
#'
#' @param out output of prisca or autprisca
#' @param p  probability level to construct the set used for plotting
#' @param cex as in R plot (there is a default)
#' @param cex.lab as in R plot (there is a default)
#' @param cex.axis as in R plot (there is a default)
#' @param lwd as in R plot (there is a default)
#' @param usedelta Default is FALSE. If TRUE, it tries to construct set such that the minimum distance between two consecutive times is delta
#' @param delta Minimum distance to be used as input if usedelta = TRUE

#'
#' @return  \code{plot}
#' @export


plot_prisca_sets <- function(out,p=0.9,cex=1.4,
                             cex.lab=1.3,cex.axis=1.3, lwd=3,usedelta=FALSE,delta=3){
  
  "Plot function of the sets. It is very much inspired by the Susie function"
  
  
  cs <- credible_set(out$alpha.mat,p)
  
  color=c("chartreuse2","aquamarine2","chocolate2","bisque2","cornflowerblue",
          "brown2","darkgoldenrod2","blue2")
  
  
  alpha.mat <- out$alpha.mat
  id <- sapply(cs, function(y) y$index)
  ymax <- max(alpha.mat)*1.1
  plot(alpha.mat[id[1],],pch=1,cex=cex,lwd=lwd,ylim=c(0,ymax),
       ylab="Posterior",xlab="Time",cex.lab=cex.lab,cex.axis=cex.axis)
  for (i in id[-1]){
    points(alpha.mat[i,],pch=1,cex=cex,lwd=lwd)
    
  }
  for (i in 1:length(id)){
    points(cs[[i]]$cs,alpha.mat[cs[[i]]$index,cs[[i]]$cs],col=color[1],cex=cex,pch=16)
    color = c(color[-1],color[1])
  }
}