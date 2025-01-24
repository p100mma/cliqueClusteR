
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


#' Convert symmetric matrix to undirected weighted graph from igraph 
#'
#' Note: this discards loops (diagonal of the matrix)
#'
#' @param W Weight matrix of the graph to construct
#' @return Undirected weighted graph object from igraph with weights from `W`
#' @examples
#' data(leukemia)
#' hist(igraph::transitivity(igr_(W), type="weighted"), main="weighted clustering coefficient")
#' @export

igr_<-function(W) igraph::graph_from_adjacency_matrix( W, mode="undirected", weighted=TRUE, diag=FALSE)


#' Convert simmilarities to dissimilarities and vice versa
#'
#' Each entry of matrix `M[i,j]` in the output gets replaced by `max(M) - M[i,j]`
#'
#' @param W Numeric matrix
#' @return A numeric matrix "flipped" over its maximum
#' @examples
#' library(magrittr)
#' data(leukemia)
#' leukemia[1:4,1:4]
#' leukemia[1:4,1:4] %>% flip_()
#' @export

flip_<-function(W) max(W) - W

#' Index both rows and columns of a matrix
#'
#' @param W Numeric matrix
#' @param idx vector of indexes or logical mask
#' @return `W[idx,idx]`
#' @export

`%[%`<-function(W,idx) W[idx,idx, drop=FALSE]

#' Calculate outer degree of the cluster in network based ow weight matrix
#' 
#' @param W Numeric matrix
#' @param mask logical mask
#' @return `sum(W[mask, !mask])` 

`%odeg%`<- function(W,mask) sum(W[mask, !mask])

#' Calculate mean outer degree of the cluster in network based ow weight matrix
#' 
#' @param W Numeric matrix
#' @param mask logical mask
#' @return `sum(W[mask, !mask])/sum(mask)` 

`%modeg%`<- function(W,mask) sum(W[mask, !mask])/sum(mask)

#' Calculate inner degree of the cluster in network based ow weight matrix

#' 
#' @param W Numeric matrix
#' @param mask logical mask
#' @return `sum(W[mask, mask])` 

`%ideg%`<- function(W,mask) sum(W[mask, mask])

#' Calculate mean inner degree of the cluster in network based ow weight matrix

#' @param W Numeric matrix
#' @param mask logical mask
#' @return `sum(W[mask, mask])/sum(mask)` 

`%mideg%`<- function(W,mask) sum(W[mask, mask])/sum(mask)


`%modularity%` <- function(S, partition) {s_gr<-S %>% igr_(); igraph::modularity(s_gr, membership=partition, weights= igraph::E(s_gr)$weight) }
