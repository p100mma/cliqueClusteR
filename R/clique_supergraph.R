#' Compute percentages of all possible connections between cliques
#' 
#' Given a weight/adjacency matrix `A` and a membership vector the function counts the present edges between
#' each pair of cliques and divides it by maximum possible number of such edges. Result of the calculation
#' is stored in a symmetric matrix: `cS[i,j]` is the percentage of possible connections between clique number `i`
#' and `j`.
#'
#' @details
#' Function assumes groups are cliques, so diagonal of the resulting matrix is set to 1 (all possible edges inside
#' each group are present).
#'
#' `WorA` is always converted into binary adjacency (1 if entry `>0`).
#' 
#' If there are `N` unique non zero labels in `cl_mem`, the resulting matrix is of size `N x N`.
#' Labels in `cl_mem` might get converted to consecutive integers by [tidyUpLabels()] and based on the result of
#' such conversion the mapping from clique label to column/row number is defined.
#' If sorted unique non zero labels in `cl_mem` were not exactly equal to `1:max(cl_mem)`, the 
#' mapping from input to output labels is saved in a form of a 2 data.frame in an attribute
#' "labelMap". 1st column contains new consecutive labels and second column the previous, original labels.
#'
#' Note that if cliques of `cl_mem` were found in a thresholded version of the similarity matrix, input `WorA` should be the thresholded version.
#' @param cl_mem an integer vector of nonnegative clique labels where 0 encodes sssumed group of nodes outside of all cliques. Length of `cl_mem` must be equal to number of rows/columns of `worA`.
#' @param WorA An adjacency or weight matrix of a graph in which the cliques encoded by `cl_mem` were found.
#' @return A symmetric matrix giving clique to clique similarity, optionally with an additional attribute "labelMap", see details.
#' @seealso [cliquePartitioneR::tidyUpLabels()], [cliquePartitioneR::greedyCliquePartitioner()], [%thr%()]
#' @examples
#' data(leukemia)
#' opt_<- thr_optimizer(leukemia)
#' cliqueSimilarity(opt_$maximizer_partition, leukemia %thr% opt_$thr)-> cS
#' dim(cS)
#' table(opt_$maximizer_partition) 
#' @export

cliqueSimilarity<- function(cl_mem, WorA) {
	stopifnot( length(dim(WorA))==2)
	stopifnot (nrow(WorA)==ncol(WorA))
	stopifnot (length(cl_mem)==ncol(WorA)) 
        	
	non0<- cl_mem[ cl_mem!=0 ]
	labelsOK=all(sort(unique(non0))==c(1:max(cl_mem)) )
	if (!labelsOK) {
		old_mem<-cl_mem
		cl_mem<-tidyUpLabels(cl_mem)
		labelMap=data.frame(new_label= cl_mem, old_label= old_mem)
		labelMap= labelMap[!duplicated(LabelMap),]
		labelMap= labelMap[order(LabelMap$new_label),]
	}	
	diag(WorA)=0
	A= (WorA>0)*1.
	C_S= cs_matrix(cl_mem=cl_mem, A=A)
	if (!labelsOK) attr(C_S, "labelMap")=labelMap
	C_S
}


cs_matrix<- function(cl_mem, A){
  stopifnot(length(unique(cl_mem))<=ncol(A))
  stopifnot(all(A %in% c(0,1)))
  #all the tests below should pass
  #if they do not, use tidyUpLabels
  cl_mem_uq<-1:max(cl_mem)
  #to optimize, cbind is slow
  componentIndicator<-do.call(cbind,
                              lapply(cl_mem_uq, function(k) cl_mem==k)) #nxk
  n_con<-A %*% componentIndicator  # n_con[i,k] number of connections of node i to cluster k
  numer<- t(componentIndicator) %*% n_con   # number of connections from cluster ki to cluster kj
  compSizes<- colSums(componentIndicator)
  denom<- compSizes %o% compSizes #for 100% correctness we would have to put..
  #... diag(denom)[[i]]<- compSizes[[i]]*(compSizes[[i]]-1)/2
  # but we can put 1s on diagonals instead.
  diag(numer)<-diag(denom)<- 1
  numer/denom
}
