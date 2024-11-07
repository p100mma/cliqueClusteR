
#' Zero out weights in similarity matrix less than specified threshold and return thresholded matrix
#'
#' Convenience function for thresholding of similarity matrix in simplified syntax
#'
#' @param W Weight matrix to be thresholded
#' @param t A value below which weights will be zeroed out
#' @return A weight matrix `W` but with weights below `t` set to zero
#' @examples
#' M<-matrix(runif(100),nrow=10)
#' print(min(M))
#' M %thr% 0.5 -> M
#' print(min(M))
#' @export

`%thr%` <- function(W, t) {W[W<t]=0; W}


#' Extract values of matrix below the diagonal
#'
#' @param W A matrix
#' @return A vector of values under the diagonal sorted as [lower.tri()] function does
#' @examples
#' data(leukemia)
#' M<- matrix(1:9, nrow=3)
#' ltr_(M)
#' @export

ltr_<-function(M) M[lower.tri(M)]
