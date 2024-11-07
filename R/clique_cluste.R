
unlapply<- function(...) unlist(lapply(...))

#internal clique clustering functions

# at this point input is assumed correct
# non 0 labels of clq_mem are consecutive positive integers from 1 to max(clq_mem)
# cS is of size n_non0labels x n_non0labels 
# max(clq_mem)== n_non0labels
# clq_importances is of length max(clq_mem) and clq_importances[[i]] is importance of clique of label "i"
# output: 
#	clust_mem$node- cluster labels of each node
#	clust_mem$cliq- cluster labels of each clique
#	core$cliq- indicator if clique is a core clique
#	core$node- indicator if node belongs to a core clique
greedyCliqueJoin<- function(cS, clq_mem, clq_importances) {
	clq_untouched=rep(TRUE, nrow(cS))
	clq_isCore=rep(FALSE, nrow(cS))
	clq_clust_mem=rep(0, nrow(cS))
	n_clust=0
	while(any(clq_untouched)){
	#if clique was processed set its importance very low
	clq_importances[ !clq_untouched ] = -999
	#pick most important available as new core
	which.max( clq_importances )-> new_core_label
	clq_isCore[new_core_label]=TRUE
	#join all yet unprocessed cliques that are connected to new core
	clq_new_clust= (cS[ new_core_label , ] > 0) & (clq_untouched)	
	#(redundant if diag(cS) is >0)
	clq_new_clust[new_core_label] = TRUE
	#update cliq cluster membership, increment cluster label
	n_clust=n_clust +1
	clq_clust_mem[new_core_label]= n_clust
	#finally, mark new cluster as processed
	clq_untouched[clq_new_clust] = FALSE	
	}
	clust_mem=list()
	core=list()
	# clique lvl membership and core indicator
	clust_mem$cliq=clq_clust_mem
	core$cliq=clq_isCore
	#now set up node cluster labels and core indicator
	clust_mem$node=rep(0, length(clq_mem))
	core$node=rep(FALSE, length(clq_mem))
	for (i_c in 1:max(clq_mem)) {
		clust_mem$node[ clq_mem==i_c ] = clust_mem$cliq[[i_c]]
		if (core$cliq[[i_c]]) core$node[ clq_mem==i_c ] = TRUE
	}	
	attr(clust_mem, "core")= core
	return(clust_mem)
}

# at this point input is assumed correct
# non 0 labels of clq_mem are consecutive positive integers from 1 to max(clq_mem)
# cS is of size n_non0labels x n_non0labels 
# max(clq_mem)== n_non0labels
cS_components<- function(cS, clq_mem) {
	clust_mem=list()
	igraph::components(
	igraph::graph_from_adjacency_matrix(
	cS, mode="undirected", weighted=TRUE, diag=FALSE
	))$membership -> clust_mem$cliq
	node_clust_mem<- rep(0, length(clq_mem))
	for (i_c in 1:max(clq_mem)) 
		node_clust_mem[ clq_mem == i_c ] = clust_mem$cliq[[i_c]]
	clust_mem$node=node_clust_mem
	return(clust_mem)
	}	


#' Relax cliques by joining cliques connected by more edges than specified
#'
#' Given node similarity graph and group membership vector, function
#' will join groups into clusters if they are connected by more than
#' `frac` of total possible number edges between them.
#' Groups are assumed to be cliques.
#'
#' @details
#' If `CS_given` is `FALSE`, then `SorCS` is assumed to be similarity matrix of nodes,
#' and group similarity matrix is constructed by [cliqueSimilarity()].
#' This includes relabeling non zero clique labels to consecutive integers from 1 to number of cliques.
#' Note that `cliqueSimilarity` does not take into the account the weights so `SorCS` can be adjacency matrix for the same effect.
#'
#' Result of relabeling is stored in the attribute `labelMap` (if it was needed) inherited from output of `cliqueSimilarity`.
#' If `CS_given` is `TRUE`, then `SorCS` is assumed to be precomputed output of [cliqueSimilarity()] and `labelMap` is inherited from it.
#'
#' `relax_method` can be currently set to one of `c('greedyCliqueJoin', 'components')`: 
#' \itemize{
#' \item `greedyCliqueJoin` - sequentially initialize each cluster by choosing most important not-yet clustered clique as an initial `core` of the cluster, then join every clique to the core if it has at least `frac` of total possible connections to the `core`. Continue until all cliques are examined. `clq_importance` by default (when it is set to `NULL`) is given by clique sizes.
#' \item `components` - clusters are connected components of the weighted clique supergraph, where nodes are cliques and weights are clique similarities computed by `[n observed edges between cliques i,j]/[n total possible edges between cliques i,j]`. Weights less than `frac` are zeroed out.
#' }
#'
#' @param partition A clique membership indicator vector (of nodes)
#' @param SorCS A similarity matrix of nodes or of cliques of nodes
#' @param frac Minimum fraction of connections allowed between cliques to be joined, number in  `[0,1]`.
#' @param CS_given If `TRUE`, then `SorCS` is assumed to be clique similarity matrix and clique similarity calculation using `partition` is not performed. 
#' @param relax_method A clique relaxation method, character string giving the name of the variant. see details.
#' @param clq_importance If not `NULL` then this argument is used in clique relaxation. Vector of numerical importance score of each clique (of length equal to number of non zero labels in `partition`). See details (only used if `relax_method=="greedyCliqueJoin"`).
#' @return A list with two named components:
#'
#' \itemize{
#' \item `cliq` - a clique membership cluster indicator vector (cliques are in the same cluster if relaxation joined them together) 
#' \item `node` - a node membership cluster indicator vector (node is in a cluster `i` if clique to which it belongs is inside cluster `i`)
#' }
#'
#' Additionally, attribute `core` is set if `relax_method=="greedyCliqueJoin"`, a list of similar structure with two elements: 
#' \itemize{
#' \item `cliq` - logical vector, entry is `TRUE` if given clique is a `core` 
#' \item `node` - logical vector, entry is `TRUE` if given node belongs to `core` clique
#' }
#'
#' (see details for meaning of `core`)
#'
#' Attribute `labelMap` is set if it was set if non zero labels of cliques were not consecutive integers, see details.
#' @examples
#' data(leukemia)
#' data(leukemia_clusters)
#' #example 1: threshold selection by optimizing "clique complexity score"
#' opt_<- thr_optimizer(leukemia)
#' cliqueSimilarity(opt_$maximizer_partition, leukemia %thr% opt_$thr)-> cS
#' #initial threshold on leukemia matrix so good we do not really need thresholding of cS for reasonable separation
#' cluster_membership<- relax_cliques(opt_$maximizer_partition, cS,
#'                                    CS_given=TRUE,        #default FALSE, means 2nd arg is cS not orig. similarity
#' 				     relax_method="components",
#'					frac=0.0) #no thresholding of cS here
#' #node lvl memberships: clique and cluster
#' table(opt_$maximizer_partition)
#' table(cluster_membership$node)
#' #clique lvl membership:
#' table(cluster_membership$cliq)
#' igraph::compare(cluster_membership$node, leukemia_clusters, 'adjusted.rand')
#' heatmap(leukemia, Rowv=NA, Colv=NA, scale="none", main="original similarity matrix")
#' heatmap(leukemia %thr% opt_thr, Rowv=NA, Colv=NA, scale="none", main="similarity thresholded by optimizing complexity")
#' 
#' #example 2: threshold selection by considering MST
#' thr_1<- critical_mst_thr(leukemia)
#' p1<- greedyCliquePartitioner(leukemia %thr% thr_1)$membership
#' cl1<- relax_cliques(p1, SorCS=leukemia %thr% thr_1, #one can give similarity instead of cS matrix
#'                       relax_method="components", frac=0.5 )
#' #first application: not enough to find the natural clusters fully
#' heatmap(leukemia %thr% thr_1, Rowv=NA, Colv=NA, scale="none", main="similarity thresholded by lvl1 critical MST")
#' igraph::compare(cl1$node, leukemia_clusters, "adjusted.rand")
#' #increasing frac improves result
#' cl1<- relax_cliques(p1, SorCS=leukemia %thr% thr_1, #one can give similarity instead of cS matrix
#'                       relax_method="components", frac=0.95 )
#' igraph::compare(cl1$node, leukemia_clusters, "adjusted.rand")
#' #another strategy: apply initial thresholding more than once
#'  #iterate 3 times
#' thr_3<- critical_mst_thr(leukemia %thr% critical_mst_thr(leukemia %thr% thr_1))
#' #much more separation
#' heatmap(leukemia %thr% thr_3, Rowv=NA, Colv=NA, scale="none", main="similarity thresholded by lvl3 critical MST")
#' p3<- greedyCliquePartitioner(leukemia %thr% thr_3)$membership
#' #cliques alone can separate clusters pretty accurately
#' igraph::compare(p3, leukemia_clusters, "adjusted.rand")
#' #relaxation by connecting cliques similar above 0.1 solves problem almost perfectly 
#' cl3<- relax_cliques(p3, leukemia %thr% thr_3, relax_method="components", frac=0.95)
#' igraph::compare(cl3$node, leukemia_clusters, "adjusted.rand") 
#' @export

relax_cliques<- function(partition,
			 SorCS,
			 frac=0.5,
			 CS_given=FALSE,
			 relax_method="greedyCliqueJoin",
			 clq_importance=NULL) {
	stopifnot(relax_method %in% c('greedyCliqueJoin', 'components'))
	stopifnot( length(frac)==1)
	stopifnot( (frac >= 0) && (frac <=1))
	#get clique sim matrix
	if (!CS_given) {
	CS= cliqueSimilarity(partition, SorCS)
	} else {CS=SorCS}
	if (!is.null(clq_importance)) stopifnot( length(clq_importance)==nrow(CS))
	#relabel if required
	if (!is.null(attr(CS,"labelMap")))
	{
	np=rep(0,length(partition))
	lM= attr(CS, "labelMap")
	for (c_i in 1:nrow(lM))
		np[ partition== lM[c_i, 2] ] = lM[c_i, 1]
	partition<-np
	}
	CS[CS<frac]=0	
	if (relax_method=="greedyCliqueJoin"){
	 if (is.null(clq_importance) ) clq_importance=  unlapply( 1:nrow(CS), function(c_i) sum(partition==c_i))
	 clust_mem<- greedyCliqueJoin( cS=CS, clq_mem=partition, clq_importances= clq_importance) 
	}
	if (relax_method=='components') {
	 clust_mem<-cS_components(cS=CS, clq_mem=partition)
	}
	
	if (!is.null(attr(CS,"labelMap"))) attr(clust_mem, "labelMap")<- attr(CS, "labelMap")
	clust_mem
}



#objective functions for clique clustering threshold choice

#' Calculate Dunn Index goodness of clustering variant based on kNN 
#'
#' For a given dissimilarity graph `dS` and partition of `dS` into 
#' clusters, function will compute Dunn Index considering 
#' average of weights of all edges between all k nearest neighbors of a pair of farthest/closest
#' nodes located inside one cluster/in different clusters, see details.
#'
#' @details
#' In general, Dunn Index is defined as `min_between_cluster_distance/max_cluster_diameter` - the bigger
#' the value, the better the clustering.
#'
#' Note that it is formulated in terms of distances, not in similarities, so to keep interpretation consistent
#' with literature this function works with dissimilarities. Below by distance we mean dissimilarity.
#' Dissimilarity can be derived from similarity matrix for example by `dS= max(S) - S`.
#'
#' For finding non-spherical clusters, this function computes `between_cluster_distance` 
#' between clusters `i,j` as a minimum distance between any nodes inside clusters `i,j`.
#' And `cluster_diameter` is defined as maximal nearest-neighbor distance in the cluster. 
#'
#' To reduce sensivity to the outliers, a paramater `k` is used to compute averages of distances in
#' deriving diameters and between cluster distances.
#'
#' Mean kNN distance between nodes `g,h` in clusters `i,j` (`i` can be equal to `j`) is defined as:
#` \itemize{
#' \item `mean( weights of all edges between N_k(g) and N_k(h) )`
#' \item where `N_k(g)` is set of `k` nearest neighbors of `g` in cluster `i`
#' \item `N_k(h)` is set of `k` nearest neighbors of `h` in cluster `j`
#'}
#' Value of the index is undefined if there is only one cluster or all clusters are singletons, in that case 
#' `NA` is returned. 
#' For the case when maximum cluster diameter is 0 value of the index is infinity 
#' and if also minimum between cluster distance is 0 then the function returns 0 (there are no differences
#' between inter and intra cluster distances so a given clustering is not correct).
#'
#' If `outlier_mask` is not `NULL`, then it is assumed that it is a logcial vector of length equal to number of nodes in the network.
#' If `outlier_mask[[i]]` is `TRUE` then node `i` gets removed from the network before calculating Dunn Index.
#' This is an ad hoc method of masking some of the particularly problematic nodes that connect clusters
#' by very few but strong connections (lowering the `between_cluster_distance`). 
#' Criteria for being an outlier here can be based on low clustering coefficient (transitivity), see [igraph::transitivity()].
#'
#' Note that even though package offers various thresholding schemes for similarity matrices, for assesing
#' goodness of the clustering one should calculate this based on dissimlarity derived from the unthresholded version
#' of the graph.
#'
#' @param dS A dissimilarity matrix of nodes in the graph. The higher the `dS[i,j]` the more dissimilar nodes `i` and `j` are.
#' @param partition A vector of membership labels of clusters of each node, 0 labels get converted to singleton labels by [cliquePartitioneR::uniqueSingletonLabels()]
#' @param outlier_mask (experimental) if not null, if `outlier_mask[i]` is `TRUE` then node `i` is not included in the calculation, see details.
#' @return A scalar value indicating the goodness of clustering, the bigger the better. Can be `Inf` if max cluster diameter is 0. If both `min_between_cluster_distance` and `max_cluster_diameter` are 0, function returns 0 by convention.
#' @examples
#' data(leukemia)
#' data(leukemia_clusters)
#' thr_1<- critical_mst_thr(leukemia)
#' p1<- greedyCliquePartitioner(leukemia %thr% thr_1)$membership
#' thr_3<- critical_mst_thr(leukemia %thr% critical_mst_thr(leukemia %thr% thr_1))
#' p3<- greedyCliquePartitioner(leukemia %thr% thr_3)$membership
#' cl3<- relax_cliques(p3, leukemia %thr% thr_3, relax_method="components", frac=0.95)
#' print("ground truth discovery scores:")
#' igraph::compare(p1, leukemia_clusters, "adjusted.rand") 
#' igraph::compare(p3, leukemia_clusters, "adjusted.rand") 
#' igraph::compare(cl3$node, leukemia_clusters, "adjusted.rand") 
#' print("corresponding values of NN_DunnIndex")
#' NN_DunnIndex(p1, max(leukemia) - leukemia)
#' NN_DunnIndex(p3, max(leukemia) - leukemia)
#' NN_DunnIndex(cl3$node, max(leukemia) - leukemia)
#' @seealso [%thr%], [igraph::modularity()], [relax_cliques()], [critical_mst_thr()]
#' @export


NN_DunnIndex<- function(partition, dS, n_NN_inside=1, 
				       n_NN_btw=2,
					 outlier_mask=NULL) {
	
	stopifnot(length(partition)==ncol(dS))	
	partition<- uniqueSingletonLabels(partition)
	if (!is.null(outlier_mask)) {
		partition<- partition[!outlier_mask]; 
		fastTable(partition)-> tabulation
		if (all(tabulation$count==1)){ warning(" outlier_mask given makes all clusters singletons");
						return(NA) }
		if (length(tabulation$count)==1){ warning(" outlier_mask given results in one cluster remaining");
						return(NA) }
		dS<- dS[!outlier_mask, !outlier_mask] } else {
		fastTable(partition)-> tabulation
		if (all(tabulation$count==1)) { warning(" All clusters are singletons."); return(NA) }
		if (length(tabulation$count)==1) {warning(" There is only one cluster in given partition.");
						  return(NA) }
		}
	# distance to NN inside each group per node
	labels<- unique(partition)
	groupPlacements<- unlapply(labels, function(x) partition==x)
	groupSizes<- unlapply(groupPlacements, sum)
	d_local_NNs<- unlapply(1:ncol(dS), function(j) 
					 if ( groupSizes[[ which(labels==partition[[j]]) ]]==1 ) {
						NA
					} else { 
			mean(	sort(dS[,j][ groupPlacements[[which(labels==partition[[j]])]] ] )[2: min(groupSizes[[j]],n_NN_inside) ] )
				         }
			       )
	diameters<- unlapply( labels, 
	#if cl is singleton then diam is 0 else its max d_local_NNs[group]
	function(label) if(sum(partition==label)==1) return(0) else  max(d_local_NNs [ partition==label ])
			    )
	max_diam=max(diameters)
	denom_is_0= (max_diam==0)
	cluster_distances=vector()
	for (i in seq_along(labels))
		for (j in seq_along(labels))
			if (i<j) {
				 dS[ partition== labels[[i]],  
				cluster_distances[[length(cluster_distances) +1 ]]=min(
				as.vector(dS[ partition== labels[[i]], partition==labels[[j]] ])
										      )
				}
	min_c_dist= min(cluster_distances)
	numer_is_0=(min_c_dist==0)
	print(max_diam)
	print(min_c_dist)
	if (denom_is_0 && numer_is_0 ) return(0) else return (min_c_dist/max_diam)
}


exhaustive_relax_thr_search<- function( clq_mem, CS, n_divisions, 
					clScoreFun, S,
					relax_method,
					relaxator_otherArgs=list(), 
					...) {
	diag(CS)<-0
	CS_ltr<- CS[lower.tri(CS)]
	ord_CS_ltr<- sort(unique(CS_ltr))
	sr=range(ord_CS_ltr)
	seq(from=sr[[1]],to=sr[[length(sr)]], length.out=n_divisions)-> seqq
	X_t=seqq
	Y_t=vector()
	partitions=list()

	for (x_t in X_t){
	do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=x_t,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		)-> membership
	  Y_t[[length(Y_t)+1]]<- clScoreFun(membership$node, S, ...)
	  partitions[[ length(partitions)+1]]= membership
	} 

	max_idx<-which.max(Y_t)
	if (max_idx==1) {l_t= X_t[[1]]; r_t=X_t[[2]] 
	} else if (max_idx== length(Y_t)) {
	  l_t= X_t[[length(Y_t)-1]]
	  r_t= X_t[[length(Y_t)]]
	} else {
	  l_t= X_t[[ max_idx -1]]
	  r_t= X_t[[ max_idx +1]]
	}
	max_neighborhood= ord_CS_ltr[ (ord_CS_ltr >= l_t) & (ord_CS_ltr <= r_t) ]
	list( X_t= X_t, Y_t= Y_t, max_idx=max_idx, max_neighborhood= max_neighborhood, partitions=partitions)

}


thr_rel_f<- function(frac, clq_mem,  CS, ordCS_ltr, clScoreFun, S,
			relax_method,
			relaxator_otherArgs=list(), ...	
		     ){
  t=frac
  stopifnot((t >= ordCS_ltr[[1]]) &&(t <= ordCS_ltr[[length(ordCS_ltr)]]) )
  diag(CS)<-0
  if (t %in% ordCS_ltr) t_exact=TRUE else t_exact=FALSE
  if (t_exact) {
	do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=t,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		)$node-> clu
	 clScoreFun(clu, ...) -> y_t
		} else {
       r.s_IDX<-min( which(ordCS_ltr > t) )
       r.s<- ordCS_ltr[[r.s_IDX]]
       l.s<- ordCS_ltr[[r.s_IDX -1 ]]
	#calculate for l.s
	do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=l.s,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		)$node-> clu
	 clScoreFun(clu, ...) -> y_l
	#then for r.s
	do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=r.s,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		)$node-> clu
	 clScoreFun(clu, ...) -> y_r
	#then interpolate :)
       A= (y_r - y_l)/(r.s-l.s)
       b= y_l - A*l.s
       y_t= A*t + b
       print(c(t, y_t))
	
	}
 return(y_t)
}

precision_relax_thr_search<- function(max_neighborhood, clq_mem,
				       CS, clScoreFun, S, 
					relax_method,
			relaxator_otherArgs=list(), 
			tol=0.1, ...) {	

	rs= range(max_neighborhood)
	l_t= rs[[1]]
	r_t= rs[[2]]
	opt_res= optimize( f= thr_rel_f, interval=c(l_t,r_t),maximum=TRUE,
	clq_mem=clq_mem, CS=CS, ordCS_ltr=max_neighborhood, 
	clScoreFun=clScoreFun, S=S, relax_method=relax_method,
	relaxator_otherArgs=relaxator_otherArgs, tol=tol, ...)
	list( opt_res=opt_res,
	      partition= do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=opt_res$maximum,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		)
	  )
}

#' Optimal clique relaxation by clique supergraph weight threshold optimization
#'
#' Given a graph partitioned into cliques this function will join them into clusters by optimizing clique similarity threshold below which clique
#' similarities will be ignored.
#'
#' @details
#'
#' Function will consider a clique similarity supergraph with weight matrix `CS` where 
#' `CS_{ij} = [n observed edges between cliques i,j]/[n total possible edges between cliques i,j]`,
#' It will try to find optimal weight threshold `frac` of `CS` such that when similarity values below `frac` are set to zero, 
#' clusters (of original nodes) resulting from clique relaxation will maximize given `cs_thr_objective` objective function.
#'
#' Clique relaxation methods are specified by argument `relax_method`, see [relax_cliques()] for description of each available method.
#' This function basically performs [relax_cliques()] step for different values of `frac` until it find a maximum of `cs_thr_objective` function.
#' 
#' Optimization is done by searching ranges of weights in `CS` in the same manner as it is done in [thr_optimizer()] (but in that function, search and thresholding is done on original node simil. matrix `S`).
#' However, here, objective function is evaluated based on original similarity matrix `S`, not `CS`. 
#' So the search is performed in two steps like in [thr_optimizer()] and arguments: 
#'   `n_init_steps`, `keep_init_eval_history`, `keep_init_partitions`, `precis`, `tol`
#' are acting exactly the same as in [thr_optimizer()], but range of weights which is searched is that of `CS` matrix and scores are obtained by 
#' evaluating `cs_thr_objective` on resulting partition and original node similarity `S`.
#'
#' `cs_thr_objective` can be a character string, then it can either be `modularity` (optimizing [igraph::modularity()]) or `nn` (optimizing [NN_DunnIndex()], converting similarities to distances).
#' This argument can also be a function. Valid function here must accept node cluster membership vector as first argument and original node to node similarity matrix `S` as the second. Additional arguments can be passed by `...`. Return value should be scalar (the bigger the better grade).
#' 
#'  `relax_method`, `clq_importance` are passed to [relax_cliques()], see that function for description of these arguments.
#' 
#' If `CS` is `NULL`, then group similarity matrix of `partition` is constructed by [cliqueSimilarity()].
#' This includes relabeling non zero clique labels to consecutive integers from 1 to number of cliques.
#' Note that `cliqueSimilarity` does not take into the account the weights so `SorCS` can be adjacency matrix for the same effect.
#'
#' Result of relabeling is stored in the attribute `labelMap` (if it was needed) inherited from output of `cliqueSimilarity`.
#' If `CS` is not `NULL`, then it is assumed to be precomputed output of [cliqueSimilarity()] and `labelMap` is inherited from it.
#'
#' @param partition A clique membership indicator vector (of nodes), 0 denotes free nodes.
#' @param S A node to node similarity matrix that was used to get `partition`
#' @param CS A clique to clique similarity matrix (assumed to be computed by [cliqueSimilarity()]
#' @param relax_method,clq_importance Clique relaxation specification, see [relax_cliques()]
#' @param cs_thr_objective A character string specifying objective function to optimize when searching for threshold or an objective function, see details.
#' @param n_init_steps,keep_init_eval_history,keep_init_partitions,precis,tol,Arguments specifying threshold searching, see details and [thr_optimizer()] for in depth explaination.
#' @return A list with following named components:
#' \itemize{
#' \item `frac` - final chosen threshold of weights of `CS`
#' \item `objective` - A value of `cs_thr_objective` corresponding to partition obtained using `frac`
#' \item `membership` - a sublist of 2: 
#' \itemize{
#' \item `cliq` - a clique membership cluster indicator vector of final partition (cliques are in the same cluster if relaxation joined them together) 
#' \item `node` - a node membership cluster indicator vector of final (node is in a cluster `i` if clique to which it belongs is inside cluster `i`)
#' }
#' \item `init_search_points`, `init_scores` Vectors of length `n_init_steps` containing search points and corresponding scores computed in initial search (set if `keep_init_eval_history` is `TRUE` and `n_init_steps` is not `NULL`) 
#' \item `init_partitions` A list of length `n_init_steps` containing clique partitions corresponding to scores and search points mentioned above (set if `keep_init_partitions` is `TRUE` and `n_init_steps` is not `NULL`). Note that here each partition is a two el. list structured like `membership` element described above. 
#' }
#' 
#' Additionally, attribute `core` is set to element `membership` if `relax_method=="greedyCliqueJoin"`, a list with two elements: 
#' \itemize{
#' \item `cliq` - logical vector, entry is `TRUE` if given clique is a `core` 
#' \item `node` - logical vector, entry is `TRUE` if given node belongs to `core` clique
#' }
#' (see details of [relax_cliques()] for meaning of `core`)
#' 
#' Attribute `labelMap` is set if non zero labels of cliques were not consecutive integers, see details. 
#' @examples
#' print("A")
#' @export


relax_cliques_optimally<- function(partition,
			  S,
			  CS=NULL,
			  relax_method="greedyCliqueJoin",
			  clq_importance=NULL,
			  cs_thr_objective="modularity",
			  n_init_steps=10,
			  keep_init_eval_history=TRUE,
			  keep_init_partitions=FALSE,
			  precis=TRUE,
			  tol=0.1, ...) {
	stopifnot(length(dim(S))==2)	
	stopifnot(ncol(S)==length(partition))
	
	exh= !(is.null(n_init_steps))
	if ((!exh)&&(!precis)) 
		stop(" n_init_steps cannot be NULL when precis is FALSE")
	if (!(class(cs_thr_objective) %in% c('character', 'function')))
		stop("cs_thr_objective must be character string or a function object")
	
	#prepare objective function
	if (class(cs_thr_objective)=="character")
	{
		if (!(cs_thr_objective %in% c("nn", "modularity")))
			stop("if cs_thr_objective is a character string it has be either 'modularity' or 'nn'")
	if (cs_thr_objective=="mst") cst_of= function(partition, S, ...) NN_DunnIndex(partition, max(S) - S, ...)
	if (cs_thr_objective=="modularity") cst_of= function(partition, S, ...) igraph::modularity( 
       igraph::graph_from_adjacency_matrix(S, mode="undirected", 
					   weighted=TRUE, diag=FALSE)	
	, partition, directed=FALSE, ...) 
	}
	if (class(cs_thr_objective)=='function') cst_of=cs_thr_objective

	#prepare CS 

	if (is.null(CS))  CS= cliqueSimilarity(partition, S)
	#relabel if required
	if (!is.null(attr(CS,"labelMap")))
	{
	np=rep(0,length(partition))
	lM= attr(CS, "labelMap")
	for (c_i in 1:nrow(lM))
		np[ partition== lM[c_i, 2] ] = lM[c_i, 1]
	partition<-np
	}
	
	if (class(relax_method)=="character")  {
		if (relax_method=="greedyCliqueJoin")
			relaxator_otherArgs=list( clq_importance= clq_importance)   } else relaxator_otherArgs=NULL

	#initial exhaustive search	
	if (exh) {
	stopifnot(n_init_steps>1)
exhaustive_relax_thr_search( clq_mem=partition, CS=CS, 
			    n_divisions=n_init_steps, 
			    clScoreFun=cst_of, S=S,
			     relax_method=relax_method,
			     relaxator_otherArgs = relaxator_otherArgs,
					...)->init_search 
	}
	if (precis) {
	#precise search
	if (exh) maxNeigh= init_search$max_neighborhood else maxNeigh=range(S)
	precision_relax_thr_search(max_neighborhood=maxNeigh, 
				  clq_mem=clq_mem, CS=CS, 
				clScoreFun=cst_of, S=S, 
					relax_method=relax_method,
			relaxator_otherArgs=relaxator_otherArgs, 
			tol=tol, ...)-> prec_search 	
	}
	
	optimization_result<- list()
	if (precis) x_f= prec_search$opt_res$maximum else x_f=init_search$X_t[[ init_search$max_idx ]]
	if (precis) y_f= prec_search$opt_res$objective else y_f=init_search$Y_t[[ init_search$max_idx ]]
		
	optimization_result$frac = x_f
	optimization_result$objective= y_f

	if (precis) mp= prec_search$partition else mp= init_search$partitions[[ init_search$max_idx ]]

	optimization_result$membership= mp
	if (exh)
	{
		if (keep_init_eval_history) {
						optimization_result$init_search_points=init_search$X_t
						optimization_result$init_scores=init_search$Y_t
					    }
		if (keep_init_partitions)       optimization_result$init_partitions=init_search$partitions

	}
	if (!is.null(attr(CS, 'labelMap'))) 
		attr(optimization_result, 'labelMap') = attr(CS, 'labelMap')
	optimization_result
}


