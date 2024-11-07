

n_clqScore<- function(clique_labels){
	length(unique(clique_labels[clique_labels!=0]))
}

complexity_clqScore<- function(clique_labels){
	clique_labels<- uniqueSingletonLabels(clique_labels)
	fastTable(clique_labels)-> tabulation
	# - [log n_cliques! + \sum_cliques log(clique size!)]
	-(lfactorial(length(tabulation$value)) + sum(lfactorial(tabulation$count)))
}


thr_cl_f<- function(t, S, ordS_ltr, clScoreFun,
			partitioner,
			partitioner_otherArgs=list(), ...	
		     ){
  stopifnot((t >= ordS_ltr[[1]]) &&(t <= ordS_ltr[[length(ordS_ltr)]]) )
  diag(S)<-0
  if (t %in% ordS_ltr) t_exact=TRUE else t_exact=FALSE
  if (t_exact) {
    S[S< t]=0
    do.call(partitioner, c(list(S), list(partitioner_otherArgs)))$membership -> clqs
    clScoreFun(clqs, ...) -> y_t
  #  prev_y=-9999
  } else {
       r.s_IDX<-min( which(ordS_ltr > t) )
       r.s<- ordS_ltr[[r.s_IDX]]
       l.s<- ordS_ltr[[r.s_IDX -1 ]]
       S[S<l.s]=0
       do.call(partitioner, c(list(S), partitioner_otherArgs))$membership -> clqs
       clScoreFun(clqs, ...) -> y_l
       S[S<r.s]=0
       do.call(partitioner, c(list(S), partitioner_otherArgs))$membership -> clqs
       #greedyCliquePartitioner(S, expansion_mode = expansion_mode, unique_singletonLabels = FALSE)$membership -> clqs
       clScoreFun(clqs, ...) -> y_r
       A= (y_r - y_l)/(r.s-l.s)
       b= y_l - A*l.s
       y_t= A*t + b
       print(c(t, y_t))
       #if (prev_y==y_t) n_eq_tries= n_eq_tries +1 
  }
  return(y_t)
    }

grid_search_exhaustive<- function( S, n_divisions, clScoreFun, 
			partitioner,
			partitioner_otherArgs=list(), ...	
			){
	  
	S_ltr<- S[lower.tri(S)]
	ord_S_ltr<- sort(unique(S_ltr))

	sr=range(ord_S_ltr)
	seq(from=sr[[1]],to=sr[[length(sr)]], length.out=n_divisions)-> seqq
	X_t=seqq
	Y_t=vector()
	partitions=list()

	Sx=S
	for (x_t in X_t){
	  Sx[Sx < x_t]=0
       	  do.call(partitioner, c(list(Sx), partitioner_otherArgs))$membership -> clqs
	  Y_t[[length(Y_t)+1]]<- clScoreFun(clqs, ...)
	  partitions[[ length(partitions)+1]]= clqs
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
	max_neighborhood= ord_S_ltr[ (ord_S_ltr >= l_t) & (ord_S_ltr <= r_t) ]
	list( X_t= X_t, Y_t= Y_t, max_idx=max_idx, max_neighborhood= max_neighborhood, partitions=partitions)
	}

precision_search_piecewiseLinear<- function(max_neighborhood, S, clScoreFun, 
			partitioner,
			partitioner_otherArgs=list(), 
			tol=0.1, ...) {	
	rs= range(max_neighborhood)
	l_t= rs[[1]]
	r_t= rs[[2]]
	opt_res=optimize(f = thr_cl_f, interval=c(l_t, r_t), maximum = TRUE,
         S=S, ordS_ltr= max_neighborhood, partitioner= partitioner, partitioner_otherArgs=partitioner_otherArgs, clScoreFun=clScoreFun, tol=tol, ...)
	S[ S< opt_res$maximum ] =0
	list(opt_res=opt_res,
       	     partition=do.call(partitioner, c(list(S), partitioner_otherArgs))$membership
	     )
}


#' Find a separating weight threshold of a dense undirected weighted graph by considering its maximum spanning tree
#' 
#' This is a heuristic that looks at the most extreme change in weight of the edges in a maximum spanning tree of the graph to deduce a threshold of weights under which edges are optimal to zero uut. See details for full description.
#' 
#' @details
#' Note that in the following we use MST to refer to maximum spanning tree (as opposed to minimum spanning tree).
#' Objective of the method is to find an optimal `thr` in a sense that in a graph described by
#' thresholded weight matrix `S_t` (where similarities lower than `thr` are zeroed out) connections inside
#' main natural clusters are perserved.
#' 
#' To deduce a separating threshold based on maximum spanning tree (MST), following steps are taken. 
#' First, a  maximum spanning tree of the graph described by input `S` is computed.
#' Then, weights of all edges of MST are sorted in the increasing order.
#' Differences between pairs of subsequent sorted weights are computed.
#' The threshold is chosen as an arithmetic mean of pairs of weights for which the computed difference was maximal.
#' 
#' The rationale for the method is the following. If there are natural clusters in the input graph, 
#' then the similarities between nodes inside those clusters are substantially bigger than similarities
#' between pairs of nodes coming from different clusters. 
#' If the natural clusters are truly separated then also the maximum between cluster similarity is substantially 
#' smaller than the minimum  of maximum similarities (computed per each node pair) in the cluster.
#' If we consider weights of MST starting from the smallest one, if one weight is suddenly very big,
#' hypothesis is that the previous weight is of the edge joining different natural clusters.
#' Moreover, if it is so, then it is also the tightest possible connection between distinct natural clusters
#' (since it is picked from MST).
#' Arithmetic mean of that weight and the next one in the sorted stack 
#' is surely bigger than this critical weight, so zeroing out connections with weights below that mean 
#' will delete only connections between natural clusters and keep tight connections inside clusters intact.
#' 
#' Note that this method always finds a threshold (we can awlays find maximum difference). 
#' Wether it is truly remarkably big difference, that is another question.
#' 
#' This method should be considered only as a initial filtering step, not an end all solution, because it might
#' happen that there are few strong connections between natural clusters (while most of the other ones
#' are very weak). This method by design will not be able to split such clusters by itself.
#' 
#' @param S A symmetric weight matrix of a graph of which threshold will be chosen based on its maximum spanning tree
#' @return A scalar value giving the found threshold
#' @seealso [igraph::mst()], [thr_optimizer()]
#' @examples
#' data(leukemia)
#' opt_thr<- critical_mst_thr(leukemia)
#' @export


critical_mst_thr<- function(S) {
#igraph gives minimum spanning tree so we switch to dissimilarity graph
dS= max(S) - S
G_ds= igraph::graph_from_adjacency_matrix(adjmatrix=dS, mode="undirected", weighted=TRUE, diag=FALSE)

igraph::mst(G_ds)-> ds_mst
#convert back to similarities so we sort weights of MAXIMIMUM spanning tree as in description
sort( max(S) - igraph::E(ds_mst)$weight, decreasing=TRUE)-> wMSTsorted
wMSTsorted<- wMSTsorted[ wMSTsorted!=0 ]
wdiffs= -1*(diff(wMSTsorted))
wmeans= (wMSTsorted[ 2: length(wMSTsorted) ] + wMSTsorted[ 1: (length(wMSTsorted)-1) ])/2
wmeans[[ which.max(wdiffs) ]]-> w_t
attr(w_t, "mst_weights")<- wMSTsorted
attr(w_t, "wdiffs")<- wdiffs
attr(w_t, "wmeans")<- wmeans
w_t

}

#' Find optimal threshold of weights of a graph based on partition scoring
#'
#' Given a predefined scoring strategy of a partition of a set of nodes into groups, function will pick 
#' such a threshold of weights `thr` such that for a thresholded graph 
#' (with connections with weight `< thr` zeroed out) the value of the score is maximized.
#' By default, the partitioning is done by splitting a graph into cliques using 
#' [cliquePartitioneR::greedyCliquePartitioner()].
#' However, by supplying partitioning function by `custom_partitioner` argument,
#' another strategy can be used, see details.
#' 
#' @details 
#' The search for the optimal threshold `thr` is performed in two stages. 
#' Either of them can be run on its own based on input arguments.
#'
#' If `n_init_steps` is not `NULL`, then the range of weights in `S` (excluding diagonal) gets sampled at 
#' `n_init_steps` uniformly spaced points and `cl_score` gets computed at each of those points.
#' 
#' If `precis` is `TRUE`, then the precise search is performed. 
#' If `n_init_steps` was not `NULL`, the precise search is done on interval containing maximum from first
#' exhaustive search: 
#' \itemize{
#' \item If maximum was found on any of the points in the interior of the range, say at step `i`, then precise search is performed on an interval between points `i-1` and `i+1`
#' \item On the other hamd, if maximum was found on the first (or last) point of the range, 
#' then the precise search is done between points `1` and `2` (or between `n_init_steps-1` and `n_init_steps` in the second case).
#' }
#' If `n_init_steps` was `NULL`, then precise search will be done on the whole range of weights of `S`.
#' 
#' Precise search works by a combination of golden section search and parabolic interpolation 
#' (see [stats::optimize()]. A function computing `cl_score` is assumed to be continuous and this is implemented
#' by linearly interpolating its values between values of weights actually present in `S`.
#' 
#' If `cl_score` is character it can be set to any of strings in `c('n_cliques', 'complexity')`: 
#' \itemize{
#' \item 'n_cliques' : the value of the score of a clique parition is number of cliques of size at least 3.
#' \item 'complexity' :the score of the partition is defined as `log` of the factorial 
#' describing clique sizes and their number: `K! * size_1! * size_2! ... size_K!` where 
#' `K` is the total number of cliques in the partition (including singletons and pairs) 
#' and `size_i` gives size of clique number `i`.
#' }
#' If `cl_score` is a function, then none of the above scores are used and the supplied custom function 
#' is used instead. A valid `cl_score` function should accept vector of membership labels of nodes to each clique
#' as first argument
#' (where 0 labels encodes free nodes outside of true cliques) and return a scalar value (the higher the better
#' the input clique partition is evaluated).
#' Other arguments for `cl_score` are passed by `...` (constant per each call in the optimization process).
#'
#' If `custom_partitioner` is not `NULL` it is assumed that it is a function.
#' It will be used for splitting nodes into groups instead of clique based method and resulting partitions
#' will be scored by `cl_score`. Valid function for `custom_partitioner` should accept similarity matrix
#' as its first argument.
#' Additional (constant throughout the optimization) arguments passed to `custom_partitioner` should be
#' supplied by `custom_partitioner_otherArgs` named list.
#' 
#' @param S A square matrix of weights of a undirected similarity graph. Assumed to be symmetric.
#' @param expansion_mode Passed to `cliquePartiotioneR::greedyCliquePartitioneR()`, a strategy for clique building.
#' @param cl_score A character string indicating a type of clique parition scoring to maximize OR a custom clique scoring function. See details for requirements for such a function.
#' @param n_init_steps If not `NULL`, number of initial steps of exhaustive search on uniformly spaced points on range of weights of `S` to compute.
#' @param keep_init_eval_history If `TRUE`, the search points and values of `cl_score` evaluated at those points are returned (applicable if `n_init_steps` is not `NULL`)
#' @param keep_init_partitions If `TRUE`, list of partitions of length `n_init_steps` is returned. Each element corresponds to partition computed at each initial search point (applicable if `n_init_steps` is not `NULL`).
#' @param precis If `TRUE`, then the function will perform precise search for a maximum of `cl_score` based on `stats::optimize()` function (possibly after initial exhaustive search if `n_init_steps` is not `NULL`).
#' @param tol Passed to `stats::optimize`, only matters if `precis` is `TRUE`. The smaller the value, the longer the precise search will take.
#' @param keep_final_partition If `TRUE`, the final partition maximizing the `cl_score` is returned along with the final score and threshold found.
#' @param custom_partitioner A custom node grouping function used instead of clique partitioning, see details.
#' @param custom_partitioner_otherArgs List of other named arguments passed to `custom_partitioner` if is used.
#' @param ... Additional arguments passed to partition scoring function (`cl_score` if it is a function).
#' @return A list with the following fields:
#' \itemize{
#' \item `thr` Threshold of weights in `S` that maximizes `cl_score` on a thresholded `S`
#' \item `objective` A value of `cl_score` corresponding to partition obtained using `thr`.
#' \item  `maximizer_partition` A partition which resulted at `cl_score` value in `objective`(set if `keep_final_partition` is `TRUE`) 
#' \item `init_search_points`, `init_scores` Vectors of length `n_init_steps` containing search points and corresponding scores computed in initial search (set if `keep_init_eval_history` is `TRUE` and `n_init_steps` is not `NULL`) 
#' \item `init_partitions` A list of length `n_init_steps` containing clique partitions corresponding to scores and search points mentioned above (set if `keep_init_partitions` is `TRUE` and `n_init_steps` is not `NULL`) 
#' }
#' @seealso [critical_mst_thr()], [%thr%()]
#' @examples
#' data(leukemia)
#' opt_res<-thr_optimizer(leukemia) 
#' print(opt_res$thr)
#' plot(opt_res$init_search_points, opt_res$init_scores)
#' abline(v= opt_res$thr)
#' abline(h= opt_res$objective)
#' table(opt_res$maximizer_partition)
#' heatmap(leukemia, scale="none", Rowv=NA, Colv=NA, main="original similarity matrix")
#' heatmap(leukemia %thr% opt_res$thr, scale="none", Rowv=NA, Colv=NA, main="thresholed similarity")
#' @export

thr_optimizer<- function( S, expansion_mode="basic", cl_score="complexity", 
			  n_init_steps=10,
			  keep_init_eval_history=TRUE,
			  keep_init_partitions=FALSE,
			  precis=TRUE,
			  tol=0.1, 
			  keep_final_partition=TRUE,
			  custom_partitioner=NULL,
			  custom_partitioner_otherArgs=NULL,
			 ...){
	exh= !(is.null(n_init_steps))
	if ((!exh)&&(!precis)) 
		stop(" n_init_steps cannot be NULL when precis is FALSE")
	if (!(class(cl_score) %in% c("character", "function")))
		stop(" cl_score must be either a character string or a function object")
	#default cl_score setups
	if (class(cl_score)=="character") 
	{
		valid_strs=c('complexity', 'n_cliques')
		if (!(cl_score %in% valid_strs))
		{
			info_str<- c(" if cl_score is a charcter string then it must be set to one of:\n",
			do.call(sprintf, c(list( paste(rep('%s', length(valid_strs)), collapse=" ")),
					 as.list(valid_strs))
			       )
				     )
			stop(info_str)
		}
	 if (cl_score=='n_cliques') clScoreFun= n_clqScore
	 if (cl_score=='complexity') clScoreFun=  complexity_clqScore
	}
	# custom cl_score
	if (class(cl_score)=='function') clScoreFun=cl_score


	# default clique partitioning setup
	if (is.null(custom_partitioner)) {
	partitioner=greedyCliquePartitioner
		partitioner_otherArgs= list( expansion_mode=expansion_mode, unique_singletonLabels=FALSE)
	} else {
	# custom partitioner given
	partitioner=custom_partitioner
		partitioner_otherArgs=custom_partitioner_otherArgs
	}		
	#initial exhaustive search	
	if (exh) {
	stopifnot(n_init_steps>1)
	grid_search_exhaustive( S=S, n_divisions=n_init_steps, 
				clScoreFun=clScoreFun,
				partitioner=partitioner,
				partitioner_otherArgs=partitioner_otherArgs,
				... 
				) -> init_search
	} 
	# precise piecewiseLinear optimization
	if (precis) {
	if (exh) maxNeigh= init_search$max_neighborhood else maxNeigh=range(S)
	precision_search_piecewiseLinear(max_neighborhood=maxNeigh, S=S, 
					clScoreFun=clScoreFun, 
					partitioner=partitioner,
					partitioner_otherArgs=partitioner_otherArgs,
					tol=tol,
					...
					)-> prec_search 

	}
	
	optimization_result<- list()
	if (precis) x_f= prec_search$opt_res$maximum else x_f=init_search$X_t[[ init_search$max_idx ]]
	if (precis) y_f= prec_search$opt_res$objective else y_f=init_search$Y_t[[ init_search$max_idx ]]
		
	optimization_result$thr = x_f
	optimization_result$objective= y_f
	
	if (keep_final_partition) 
	{
		if (precis) mp= prec_search$partition else mp= init_search$partitions[[ init_search$max_idx ]]
		optimization_result$maximizer_partition= mp
	} 		
	if (exh)
	{
		if (keep_init_eval_history) {
						optimization_result$init_search_points=init_search$X_t
						optimization_result$init_scores=init_search$Y_t
					    }
		if (keep_init_partitions)       optimization_result$init_partitions=init_search$partitions

	}
	optimization_result
}

