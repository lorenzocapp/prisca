#' Plot the data along with the credible sets. It is very much inspired by the corresponding function of the Susie package (Wang et al. 20)
#'
#' @param out output of prisca or autprisca
#' @param p  probability level to construct the set used for plotting
#' @param cex as in R plot (there is a default)
#' @param cex.lab as in R plot (there is a default)
#' @param cex.axis as in R plot (there is a default)
#' @param lwd as in R plot (there is a default)
#' @param usedelta Default is FALSE. If TRUE, it tries to construct set such that the minimum distance between two consecutive times is delta
#' @param delta Minimum distance to be used as input if usedelta = TRUE
#' @param option It defines where to color the credible sets. Default is "data", meaning that credible sets are coloured directly in the data points. If anything else, colors of the credible sets in the bottom

#'
#' @return  \code{plot}
#' @export




plot_prisca_data_sets <- function(out,p=0.9,cex=1.4,
                                  cex.lab=1.3,cex.axis=1.3, lwd=3,
                                  usedelta=FALSE,delta=3,
                                  option="data"){
  
  "Plot data together with credible sets. If option data, the color is inside
    the circle, if option is different, the color is at the bottom.
    "
  
  y <- out$y
  cs <- credible_set(out$alpha.mat,p)
  id <- sapply(cs, function(y) y$index)
  color=c("chartreuse2","aquamarine2","chocolate2","bisque2","cornflowerblue",
          "brown2","darkgoldenrod2","blue2")
  ymax <- max(y)
  ymin <- min(y)
  
  plot(y,pch=1,cex=cex,lwd=lwd,ylim=c(ymin-0.11*ymax,ymax+0.11*ymax),
       ylab="y",xlab="Time",cex.lab=cex.lab,cex.axis=cex.axis)
  
  if (option=="data"){
    for (i in 1:length(id)){
      points(cs[[i]]$cs,y[cs[[i]]$cs],col=color[1],cex=cex,pch=16)
      color = c(color[-1],color[1])
    }
  } else {
    loc <- ymin-0.11*ymax
    for (i in 1:length(id)){
      points(cs[[i]]$cs,rep(loc,length(cs[[i]]$cs)),bg=color[1],cex=cex,pch=22,col="black")
      color = c(color[-1],color[1])
    }
  }
}

