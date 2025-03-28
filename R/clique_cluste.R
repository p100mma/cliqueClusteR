
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



covers_sNN<- function(cS) {
cover_indicators<- cS > 0 
diag(cover_indicators)<- 0
t(cover_indicators) %*% cover_indicators
}


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
	clq_clust_mem[clq_new_clust]= n_clust
	#finally, mark new cluster as processed
	clq_untouched[clq_new_clust] = FALSE	
	}
	
	#find the strongest pull of each clique to each cluster
clust_indicators<- matrix(nrow=ncol(cS), ncol=n_clust)
for ( i in 1:n_clust) clust_indicators[,i] = (clq_clust_mem==i)*1.
# c x c times c x cl -> c x cl
cS %*% clust_indicators %>% apply(1, which.max) -> attractor_label
#check which clique have strongest pull elsewhere and are not themselves cores
clustSizes<- colSums(clust_indicators)
sizeOfClustPerClq<- lapply(clq_clust_mem, function(k) clustSizes[[k]] ) %>% unlist()
pulled_on_cliques <- which( (clq_clust_mem!=attractor_label) &
			     (sizeOfClustPerClq>1)  &(!clq_isCore)
					)
n_swaps=0
while(length(pulled_on_cliques))
{

#calculate local gain on a swap 
lapply(pulled_on_cliques, function(i)
	{ old_label=clq_clust_mem[i]
	  hypothetical = clq_clust_mem	
	  new_label= attractor_label[i]
	  hypothetical[i] = new_label
	   cl1<-(clq_clust_mem==new_label)
	   cl2<-(clq_clust_mem==old_label)
	   hcl1<-(hypothetical==new_label)
            hcl2<-(hypothetical==old_label)
	  old_cl1_coef<- (cS %ideg% cl1)/(cS %odeg% cl2)
	  old_cl2_coef<- (cS %ideg% cl2)/(cS %odeg% cl2)
	  new_cl1_coef<- (cS %ideg% hcl1)/(cS %odeg% hcl1)
	  new_cl2_coef<- (cS %ideg% hcl2)/(cS %odeg% hcl2)
	  might_swap = (old_cl1_coef < new_cl1_coef) && (old_cl2_coef < new_cl2_coef)
	  ifelse(might_swap,
		 (new_cl1_coef - old_cl1_coef) + (-old_cl2_coef + new_cl2_coef),
		-9999)
	}
	) %>% unlist() -> swap_gains
if (all(swap_gains <=0 )) pulled_on_cliques<- c() else {
	 n_swaps=n_swaps+1
	 swap_to_do<- pulled_on_cliques [which.max(swap_gains)]
	 new_label <- attractor_label[ swap_to_do ]
	 clq_clust_mem[ swap_to_do ] = new_label  
	 clust_indicators<- matrix(nrow=ncol(cS), ncol=n_clust)
	 for ( i in 1:n_clust) clust_indicators[,i] = (clq_clust_mem==i)*1.
	 # c x c times c x cl -> c x cl
	 cS %*% clust_indicators %>% apply(1, which.max) -> attractor_label
	 sizeOfClustPerClq<- lapply(clq_clust_mem, function(k) clustSizes[[k]] ) %>% unlist()
	 #check which clique have strongest pull elsewhere and are not themselves cores
	 pulled_on_cliques <- which( (clq_clust_mem!=attractor_label) &
			     (sizeOfClustPerClq>1) & (!clq_isCore)
						 )
	}
#print(length(pulled_on_cliques))
}
print("n_swaps made:")
print(n_swaps)
# examine joins of clusters of size >1
#	 is [mean_in_degree/mean_out_degree ](joined) > [mean_in_degree/mean_out_degree] (both separate?)
	#cluster 2 cluster pull
	t(clust_indicators) %*% cS %*% clust_indicators -> cl2cl
	#pull inside clusters
	diag(cl2cl)-> cl2in
	cl2out<- cl2cl; diag(cl2out)<-0; 
	#cluster with highest aggregate connection strength to the cluster [[k]]
	attractor_label<- apply(cl2out,2, function(neigh) which.max(neigh) )
	#strength of that pull
	cl2out<- lapply(1:ncol(cl2out), function(k) cl2cl[k, attractor_label[[k]]] ) %>% unlist()
	clSizes<- colSums(clust_indicators)
	pulled_clusters<- which( (cl2in < cl2out) & (clSizes>1) )
	n_merges=0
	while (any(pulled_clusters)){
	lapply(pulled_clusters, 
		function(i)
		{ 
		  old_label=i
		  new_label= attractor_label[i]
		  hypothetical = clq_clust_mem	
		  hypothetical[clq_clust_mem==i] = new_label
		  #old clusters cl1 and cl2
		  cl1<- clq_clust_mem==old_label
		  cl2<- clq_clust_mem==new_label
		  #hypothetical cluster after merging cl1 and cl2
		  hcl<- hypothetical==new_label 
		  old_cl1_coef<- (cS %ideg%  cl1)/(cS %odeg% cl1) 
		  old_cl2_coef<- (cS %ideg%  cl2)/(cS %odeg% cl2) 
		  hcl_coef<- (cS %ideg% hcl)/(cS %odeg% hcl) 
	  	  might_merge = (old_cl1_coef < hcl_coef) && (old_cl2_coef < hcl_coef)
		  ifelse(might_merge,
			hcl_coef - mean(c(old_cl1_coef, old_cl2_coef)),
			-9999	
			)
			}
	      ) %>% unlist() -> merge_gains
	if (all(merge_gains <= 0)) pulled_clusters= c() else  {
	n_merges=n_merges+1
	 merge_to_do<- pulled_clusters [which.max(merge_gains)]
	 new_label <- attractor_label[ merge_to_do ]
	 clq_clust_mem[ clq_clust_mem==merge_to_do ] = new_label  
	 clust_indicators<- matrix(nrow=ncol(cS), ncol=n_clust)
	 for ( i in 1:n_clust) clust_indicators[,i] = (clq_clust_mem==i)*1.
	 # c x c times c x cl -> c x cl
	 t(clust_indicators) %*% cS %*% clust_indicators -> cl2cl
	 #pull inside clusters
	 diag(cl2cl)-> cl2in
	 cl2out<- cl2cl; diag(cl2out)<-0; 
	 #cluster with highest aggregate connection strength to the cluster [[k]]
	 attractor_label<- apply(cl2out,2, function(neigh) which.max(neigh) )
	 #strength of that pull
	 cl2out<- lapply(1:ncol(cl2out), function(k) cl2cl[k, attractor_label[[k]]] ) %>% unlist()
	 clSizes<- colSums(clust_indicators)
	 pulled_clusters<- which( (cl2in < cl2out) & (clSizes>1) )

		}

	}
	print("n merges made:")
	print(n_merges)
			
#	# check if single core merges improve the out deg min deg coefficients
#	clq_size_of_cluster<- lapply(1:nrow(cS), function(i) sum(clq_clust_mem==i)) %>% unlist() 	
#	isolated_cores<- which(clq_isCore & (clq_size_of_cluster==1))
#	clust_indicators<- matrix(nrow=nrow(cS), ncol=n_clust)
#	for ( i in 1:n_clust) clust_indicators[,i] = (clq_clust_mem==i)*1.
#	#we dont want to count contribution of j to itself
#	for ( j in isolated_cores) clust_indicators[j, clq_clust_mem[j] ]= 0
#	# i x c times c x cl -> i x cl
#	cS[isolated_cores,,drop=FALSE] %*% clust_indicators %>% apply(1, function(x) ifelse(max(x)>0, which.max(x),-999 )) -> attractor_label
#	pulled_cores<- isolated_cores[ attractor_label>0 ]
#	isolated_cores<- isolated_cores[ attractor_label>0 ]
#	attractor_label<- attractor_label [ attractor_label >0 ]
#	while (length(pulled_cores)) {
#	#calculate local gain on a merge on an isolated core
#	lapply(pulled_cores, function(i)
#		{ old_label=clq_clust_mem[i]
#		  hypothetical = clq_clust_mem	
#		  new_label= attractor_label[i]
#		  hypothetical[i] = new_label
#		  old_cl1_coef<- cS %ideg% (clq_clust_mem==new_label)
#		  old_cl1_coef<- old_cl1_coef/( cS %odeg% (clq_clust_mem==new_label))
#		  new_cl1_coef<- cS %ideg% (hypothetical==new_label)
#		  new_cl1_coef<-new_cl1_coef/( cS %odeg% (hypothetical==new_label))
#		  might_swap = (old_cl1_coef < new_cl1_coef) 
#		  ifelse(might_swap,
#			 (new_cl1_coef - old_cl1_coef) 
#			-9999)
#		}
#		) %>% unlist() -> swap_gains
#	print(swap_gains)
#	if (all(swap_gains <=0 )) pulled_cores<- c() else {
#		 swap_to_do<- pulled_cores [which.max(swap_gains)]
#		 new_label <- attractor_label[ swap_to_do ]
#		 clq_clust_mem[ swap_to_do ] = new_label  
#		 isolated_cores<- isolated_cores[-c(swap_to_do)]
#		 n_clust<- n_clust - 1
#	if (!length(isolated_cores)) { pulled_cores=c() } else {
#	clust_indicators<- matrix(nrow=nrow(cS), ncol=n_clust)
#	for ( i in 1:n_clust) clust_indicators[,i] = (clq_clust_mem==i)*1.
#	#we dont want to count contribution of j to itself
#	for ( j in isolated_cores) clust_indicators[j, clq_clust_mem[j] ]= 0
#	# i x c times c x cl -> i x cl
#	cS[isolated_cores,,drop=FALSE] %*% clust_indicators %>% apply(1, function(x) ifelse(max(x)>0, which.max(x),-999 )) -> attractor_label
#	pulled_cores<- isolated_cores[ attractor_label>0 ]
#	isolated_cores<- isolated_cores[ attractor_label>0 ]
#	attractor_label<- attractor_label [ attractor_label >0 ]
#		}
#	}
#	}	

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
#' heatmap(leukemia %thr% opt_$thr, Rowv=NA, Colv=NA, scale="none", main="similarity thresholded by optimizing complexity")
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
			 clq_importance=NULL,
			 do_signif=FALSE,
			 lvl=0.05
			 ) {
	stopifnot(relax_method %in% c('greedyCliqueJoin', 'components'))
	stopifnot( length(frac)==1)
	stopifnot( (frac >= 0) && (frac <=1))
	#get clique sim matrix
	if (!CS_given) {
	CS= cliqueSimilarity(partition, SorCS, do_signif, lvl)
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

#' Calculate Dunn Index goodness of clustering variant based on internal MST of clusters 
#'
#' For a given dissimilarity graph `dS` and partition of `dS` into 
#' clusters, function will compute Dunn Index considering 
#' minimum spanning tree (MST) of each cluster to get diameter
#' and optionally a Hausdorff distance for between cluster distances instead of 
#' original single linkage distance, see details.
#'
#' @details
#' In general, Dunn Index is defined as `min_between_cluster_distance/max_cluster_diameter` - the bigger
#' the value, the better the clustering.
#'
#' Note that it is formulated in terms of distances, not in similarities, so to keep interpretation consistent
#' with literature this function works with dissimilarities. Below by distance we mean dissimilarity.
#' Dissimilarity can be derived from similarity matrix for example by `dS= max(S) - S`.
#'
#' As first proposed by Dunn [ADD REFERENCE], `between_cluster_distance` 
#' between clusters `i,j` as a minimum distance between any nodes inside clusters `i,j`, like in
#' single linkage method of hierarchical clustering.
#' If argument `dX.Y` is set to `"hausdorff"`, then `between_cluster_distance` changes to
#' Hausdorff distance: 
#' \itemize{
#' \item `d_H(X,Y):= max ( max_x d(x,Y) , max_y d(y,X) )`
#' \item where `X,Y` are clusters, `x` points in `X` and `y` points in `Y`,
#' \item `d(x,Y):= min_y d(x,y)`.
#' }
#' This makes the index less influenced by few points of contact between clusters if in general they are 
#' well separated.
#'
#' And `cluster_diameter` is defined here as maximum distance on MST of a cluster. 
#' Such choices make the index suited for finding elongated clusters in which distances only between
#' adjacent neighbors are small but distances between indirect neighbors are allowed ot be big.
#'
#' Value of the index is undefined if there is only one cluster or all clusters are singletons, in that case 
#'  0 is returned by convention. 
#' For the case when maximum cluster diameter is 0 value of the index is infinity. 
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
#' @param partition A vector of membership labels of clusters of each node, 0 labels get converted to singleton labels by [cliquePartitioneR::uniqueSingletonLabels()]
#' @param dS A dissimilarity matrix of nodes in the graph. The higher the `dS[i,j]` the more dissimilar nodes `i` and `j` are.
#' @param dX.Y A character string specifying type of between cluster distance, see details.
#' @param outlier_mask (experimental) if not null, if `outlier_mask[i]` is `TRUE` then node `i` is not included in the calculation, see details.
#' @return A scalar value indicating the goodness of clustering, the bigger the better. Can be `Inf` if max cluster diameter is 0. If both `min_between_cluster_distance` and `max_cluster_diameter` are 0, function returns 0 by convention.
#' @examples
#' library(magrittr)
#' points<- matrix(rnorm(800*3), ncol=3)
#' points<-points/sqrt(rowSums(points^2)) #unit sphere
#' points[1:400,]= points[1:400,]*20 # outer
#' points[401:800,]= points[401:800,]*0.5 #inner
#' ref=rep(1,800); ref[401:800]=2 # reference labels for points
#' plot(points, col=ifelse(ref>1,"red","blue"))
#' dS<- as.matrix(dist(points))
#' flip_<- function(M) max(M)- M #convert similarities to distances and vice versa
#' dS %>% flip_() -> S
#' thr<- critical_mst_thr(S)
#' # first get cliques
#' pa<- greedyCliquePartitioner(S %thr% thr)$membership
#' #not perfect yet
#' igraph::compare(pa, ref, "adjusted.rand")
#' #corresponding dunn indices:
#' MinST_DunnIndex(ref,dS) #reference
#' MinST_DunnIndex(pa,dS) #cliques
#' print("hausdorff:")
#' MinST_DunnIndex(ref,dS, "hausdorff") #reference
#' MinST_DunnIndex(pa,dS, "hausdorff") #cliques
#' # clusters from clique relaxation
#' relax_cliques(pa %>% uniqueSingletonLabels(), S %thr% thr, relax_method = "components")-> cl
#' igraph::compare(ref, cl$node, "adjusted.rand")
#' MinST_DunnIndex(cl$node,dS)
#' MinST_DunnIndex(cl$node,dS,"hausdorff")
#' #fine tune frac of kept edges between cliques, better score:
#' relax_cliques(pa %>% uniqueSingletonLabels(), S %thr% thr, relax_method = "components", frac=0.1)-> cl
#' igraph::compare(ref, cl$node, "adjusted.rand")
#' #indexes are better
#' MinST_DunnIndex(cl$node,dS)
#' MinST_DunnIndex(cl$node,dS,"hausdorff")
#' @seealso [%thr%], [igraph::modularity()], [relax_cliques()], [critical_mst_thr()]
#' @importFrom magrittr %>%
#' @export

MinST_DunnIndex<-function(partition, dS,dX.Y="single_linkage", 
			    hausdorff_q=1,
			 outlier_mask=NULL){
  stopifnot(length(partition)==ncol(dS))	
  stopifnot(class(dX.Y)=="character")
  stopifnot(dX.Y %in% c("single_linkage", "hausdorff"))
  partition<- uniqueSingletonLabels(partition)
  if (!is.null(outlier_mask)) dS<- dS[!outlier_mask, !outlier_mask] 
  if (!is.null(outlier_mask)) partition<- partition[!outlier_mask ]
  labels=unique(partition)
  if (length(labels)==1) return(0)
  placements = lapply(labels, function(x) which(partition==x))
  lapply(placements, length) %>% unlist() -> sizes
  if (all(sizes==1)) return(0)
  lapply(seq_along(placements), function(i) 
  { if(sizes[[i]]==1) return(0) else{
        dS[placements[[i]],
           placements[[i]]] %>% igr_() %>% igraph::mst() %>%
            igraph::E() -> internal_skeleton
            quantile((internal_skeleton)$weight,hausdorff_q)
                                     }
    }
    ) %>% unlist() -> diameters
   max_diam= max(diameters)
   ij=0
   distances<-vector()
   #for efficiency, since if there are singletons and threshold is high
   #it might happen that almost all clusters are singletons
   #then following loop is very slow
   #so first handle clusters of size>1
   non_singles<- which(sizes!=1)
   for (i in non_singles)
     for (j in non_singles)
       if (i<j)
     {ij=ij+1; 
      if (dX.Y=="single_linkage")
     dS[ placements[[i]], 
         placements[[j]] ] %>% min() -> distances[[ij]] 
      if (dX.Y=="hausdorff")
      {
        dS[ placements[[i]], 
            placements[[j]], drop=FALSE ] %>% apply(2,min) -> dx.Y
        dS[ placements[[i]], 
            placements[[j]], drop=FALSE ] %>% apply(1, min)-> dy.X
        c(quantile(dx.Y, hausdorff_q),
	 quantile(dy.X, hausdorff_q)) %>% max() -> distances[[ij]]
      }
       }
   # singletons to singletons
   singles= which(sizes==1)
   placements[singles] %>% unlist () -> singleton_placement
   sdS<- dS[ singleton_placement, singleton_placement ]	
   distances= c(distances, sdS[ lower.tri(sdS) ]) 
   rm(sdS)
   ij<- length(distances)
   # singletons to non singletons
   for (i in singles)
     for (j in non_singles)
       if (i<j)
     {ij=ij+1; 
     dS[ placements[[i]], 
         placements[[j]] ]-> ij_distances #if i is singleton then this is just 1d vector
      if (dX.Y=="single_linkage")
	distances[[ij]] <- min(ij_distances)
      if (dX.Y=="hausdorff")
	distances[[ij]]<- max( min(ij_distances), #distance of singleton to big cluster is just the minimal distance
						  #quantile over 1 element set is just that element.
			        quantile(ij_distances, hausdorff_q)  #distance of each el. of big cluster to 
								     #singleton cluster is just distance of that
								     #element. there are as many of them as elements
								     #of big cluster so we take quantile over them
			      )
      }
       
   min_dist<-min(distances)
   if ((max_diam==0) && (min_dist==0)) return(0)
   min_dist/max_diam
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
	if (max_idx==1) {l_i= 1; r_i=2
	} else if (max_idx== length(Y_t)) {
	  l_i= length(Y_t)-1
	  r_i= length(Y_t)
	} else {
	  l_i=  max_idx -1
	  r_i=  max_idx +1
	}
	#print(c(X_t[[l_i]], X_t[[r_i]]))
	max_neighborhood= c( X_t[[l_i]], ord_CS_ltr[ (ord_CS_ltr >= X_t[[l_i]]) & (ord_CS_ltr <= X_t[[r_i]]) ],
			     X_t[[r_i]])
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
	 clScoreFun(clu, S, ...) -> y_t
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
	 clScoreFun(clu, S, ...) -> y_l
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
	 clScoreFun(clu, S, ...) -> y_r
	#then interpolate :)
       A= (y_r - y_l)/(r.s-l.s)
       b= y_l - A*l.s
       y_t= A*t + b
       #print(c(t, y_t))
	
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
	if (l_t==r_t) {
		opt_res=list()
		opt_res$maximum= l_t
			} else {
		 opt_res= optimize( f= thr_rel_f, interval=c(l_t,r_t),maximum=TRUE,
		 clq_mem=clq_mem, CS=CS, ordCS_ltr=max_neighborhood, 
		 clScoreFun=clScoreFun, S=S, relax_method=relax_method,
		 relaxator_otherArgs=relaxator_otherArgs, tol=tol, ...)
				}
	opt_pa= do.call(relax_cliques, 
	c(	list( partition= clq_mem,
		      SorCS=CS,
		      frac=opt_res$maximum,
		      CS_given=TRUE,
		      relax_method=relax_method),
		relaxator_otherArgs 
	  )
		       )  
	opt_res$objective= clScoreFun(opt_pa$node, S, ...) 
	list( opt_res=opt_res,
	      partition=opt_pa 	  )
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
# It will try to find optimal weight threshold `frac` of `CS` such that when similarity values below `frac` are set to zero, 
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
#' `cs_thr_objective` can be a character string, then it can either be `modularity` (optimizing [igraph::modularity()]) or `mst` (optimizing [MinST_DunnIndex()], converting similarities to distances).
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
#' @param S An original, non thresholded node similarity matrix
#' @param t_S A threshold that was used on `S` to find `partition`
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
#'
#' @examples
#' library(magrittr)
#' data(leukemia)
#' #example of not so good similarity threshold 
#' #optimal clique relaxation tries to get most of it anyway
#' opt_= list( thr=critical_mst_thr(leukemia) )
#' opt_$maximizer_partition=greedyCliquePartitioner(leukemia %tht% opt_$thr)
#' relax_cliques_optimally(opt_$maximizer_partition %>% uniqueSingletonLabels(), leukemia, opt_$thr,
#'                         n_init_steps=2)-> relax_result #toy example,init steps to 2 to make prec. search useful
#' data(leukemia_clusters) #to compare with discovered clusters
#' igraph::compare(leukemia_clusters, relax_result$membership$node, "adjusted.rand")
#' #plot the initial search history
#' plot(relax_result$init_search_points, relax_result$init_scores, xlab="thresholds checked in initial search",
#' ylab="objective function value", type="p", pch=16)
#' # mark the result of precision search also
#' abline(v=relax_result$frac)
#' abline(h=relax_result$objective)
#' @export


relax_cliques_optimally<- function(partition,
			  S,
			  t_S,
			  CS=NULL,
			  relax_method="greedyCliqueJoin",
			  clq_importance=NULL,
			  cs_thr_objective="modularity",
			  n_init_steps=10,
			  keep_init_eval_history=TRUE,
			  keep_init_partitions=FALSE,
			  precis=TRUE,
			  tol=0.1,
			  do_signif=TRUE,
			  lvl=0.05, ...) {
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
		if (!(cs_thr_objective %in% c("mst", "modularity")))
			stop("if cs_thr_objective is a character string it has be either 'modularity' or 'mst'")
	if (cs_thr_objective=="mst") cst_of= function(partition, S, ...) MinST_DunnIndex(partition, max(S) - S, ...)
	if (cs_thr_objective=="modularity") cst_of= function(partition, S, ...) {
	igr_(S)-> gS	
	igraph::modularity( gS 
	, partition, directed=FALSE, weight=igraph::E(gS)$weight,, ...) 
										}
	}
	if (class(cs_thr_objective)=='function') cst_of=cs_thr_objective

	#prepare CS 

	if (is.null(CS))  CS= cliqueSimilarity(partition, S %thr% t_S, do_signif, lvl)
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
		if (relax_method=="greedyCliqueJoin") {
			relaxator_otherArgs=list( clq_importance= clq_importance)   } else relaxator_otherArgs=NULL
						}

	#initial exhaustive search	
	if (exh) {
	#print("CS range:")
	#rs<-range(CS)
	#print(seq(rs[[1]], rs[[2]], length.out=n_init_steps))
	stopifnot(n_init_steps>1)
exhaustive_relax_thr_search( clq_mem=partition, CS=CS, 
			    n_divisions=n_init_steps, 
			    clScoreFun=cst_of, S=S,
			     relax_method=relax_method,
			     relaxator_otherArgs = relaxator_otherArgs,
					...)->init_search 
	#plot(init_search$X_t, init_search$Y_t)
	}
	if (precis) {
	#precise search
	if (exh) maxNeigh= init_search$max_neighborhood else maxNeigh=range(CS)
	#print("maxNeigh:")
	#print(maxNeigh)
	precision_relax_thr_search(max_neighborhood=maxNeigh, 
				  clq_mem=partition, CS=CS, 
				clScoreFun=cst_of, S=S, 
					relax_method=relax_method,
			relaxator_otherArgs=relaxator_otherArgs, 
			tol=tol, ...)-> prec_search 	
	print("init vs prec")
	print(init_search$Y_t[[ init_search$max_idx ]] )
	print(prec_search$opt_res$objective)
	if (exh)  #sometimes prec search fails
		if (init_search$Y_t[[ init_search$max_idx ]] > prec_search$opt_res$objective)
			{ 
			 print("prec failed")
			 prec_search$opt_res$maximum= init_search$X_t[[ init_search$max_idx ]]
			 prec_search$opt_res$objective=init_search$Y_t[[ init_search$max_idx ]]	
			 prec_search$partition= init_search$partitions[[ init_search$max_idx ]]
			}
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


