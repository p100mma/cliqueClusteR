#' Run Simple Clique Paritioning based Clustering of Similiarty Spaces
#'
#' Given an input matrix of pairwise node similarities or coordinates of objects in some metric space `XorS`,
#' function will cluster the nodes/points based on their similarity matrix by
#' partitioning similarity graph into cliques  and relaxing the cliques into clusters, both (optionally) in an automatic,
#' optmial way.
#'
#' @param XorS Either a matrix where each row gives coordinates of points between which to compute similarities or node similarity matrix, see details.
#' @param t_S A threshold choice for similarity matrix. Either character string determining strategy to choose this value optimally, a numeric value giving exact thr. to use or a function computing this thr. based on similarity matrix `S`, see details.
#' @param t_CS A threshold choice of clique similarity matrix used in clique relaxation step. Either a character string specifying strategy to use for finding this value or an exact numeric value to use, see details.
#' @param X_simil_fun If not `NULL`, then it is assumed `XorS` is a matrix of coordinates of points (one row per one point). If a character string, then it specifies distance metric to use for computing similarities. If a function, then it will be used  on `XorS` to produce similarity matrix for clique partitioning.
#' @param return_S,return_CS Logical values, if `TRUE`, then return similarity matrix `S` and clique similarity matrix `CS` utilized, respectively.
#' @param expansion_mode A clique partitioning strategy to use, see [cliquePartitioneR::greedyCliquePartitioner()].
#' @param t_S_init_steps,t_S_keep_init_eval_history,t_S_keep_init_partitions,t_S_precis,t_S_tol Passed to [thr_optimizer()] if t_S is set to one of `"complexity","n_cliques"`. These parameters control threshold optimization process for similarity matrix `S`, see [thr_optimizer()] for in depth description.
#' @param t_CS_relax_method Character string specifying clique relaxation method to use, see [relax_cliques()].
#' @param t_CS_n_init_steps,t_CS_keep_init_eval_history,t_CS_keep_init_partitions,t_CS_precis,t_CS_tol,t_CS_dX.Y  Passed to [relax_cliques_optimally()] if `t_CS` is set to one of `"mst","modularity"`. These parameters control threshold optimization process for clique similarity matrix `CS`, see [relax_cliques_optimally()] for in depth description.
#' @return A named list with the following components: 
#' \itemize{
#' \item `clique_membership` A clique membership vector of nodes. Singletons are given unique labels.
#' \item `t_S` A numerical value of threshold used on similarity matrix `S` before clique partitioning.
#' \item `t_CS` A numerical value of threshold used on clique similarity matrix in clique relaxation process.
#' \item `cluster_membership` A sublist of 2 elements, each giving a kind of membership vector:
#' \itemize{
#'  \item `cliq` gives cluster labes to which cliques belong
#'  \item `node` gives cluster labels to which nodes belong (cluster of node `i` is the cluster of clique in which `i` is situated)
#' }
#' \item `S` A symmetic matrix of node similarities used (non thresholded, set if `return_S` was `TRUE`) 
#' \item `t_S_objective` A value of objective function obtained for final used `t_S` (set if threshold optimization of matrix `S` was done)
#' \item `t_S_init_partitions` Partitions corresponing to inital search for `t_S` (set if threshold optimization of matrix `S` was done and `t_S_keep_init_partitions` was `TRUE`)
#' \item `t_S_init_search_points` Initial search points for `t_S`  (set if inital threshold optimization of matrix `S` was done and `t_S_keep_init_score` was `TRUE`)
#' \item `t_S_init_scores` Initial scores of each initial threshold in `t_S_init_search_points`  (set if inital threshold optimization of matrix `S` was done and `t_S_keep_init_score` was `TRUE`)
#' \item `CS` A symmetic matrix of clique similarities used (non thresholded, set if `return_CS` was `TRUE`) 
#' \item `t_CS_objective` A value of objective function obtained for final used `t_CS` (set if threshold optimization of matrix `CS` was done)
#' \item `t_CS_init_partitions` Partitions corresponing to inital search for `t_CS` (set if initial threshold optimization of matrix `CS` was done and `t_CS_keep_init_partitions` was `TRUE`)
#' \item `t_CS_init_scores` Initial scores of each initial threshold in `t_CS_init_search_points`  (set if inital threshold optimization of matrix `CS` was done and `t_CS_keep_init_score` was `TRUE`)
#' }
#' @export
run_scpcss<- function(XorS,                           
		     t_S="complexity",
		     t_CS="modularity",
		     X_simil_fun=NULL,
		     return_S=FALSE,
		     return_CS=FALSE,
		     expansion_mode="basic",
		     t_S_init_steps=10,
		     t_S_keep_init_eval_history=TRUE,
		     t_S_keep_init_partitions=FALSE,
		     t_S_precis=TRUE,
		     t_S_tol=0.1, 
		     t_CS_relax_method="greedyCliqueJoin",
		     t_CS_n_init_steps=10,
		     t_CS_keep_init_eval_history=TRUE,
		     t_CS_keep_init_partitions=FALSE,
		     t_CS_precis=TRUE,
		     t_CS_tol=0.1,
		     t_CS_dX.Y="hausdorff",
		     do_signif=FALSE,
		     lvl=0.05) {

	stopifnot(length(dim(XorS))==2)
	stopifnot(!is.null(ncol(XorS)))
	stopifnot(!is.null(nrow(XorS)))
	stopifnot( (ncol(XorS)>0) && (nrow(XorS)>0))
		
	#case I: X is given
	if ((nrow(XorS)!= ncol(XorS)) && (is.null(X_simil_fun)))
		stop(paste0("if XorS is not square matrix then\n",
		      "it is assumed its coordinate matrix of objects,",
		      "and X_simil_fun must be given"))

	stopifnot(class(t_S) %in% c("character", "function","numeric"))	
	stopifnot(class(t_CS) %in% c("character","numeric"))	
	
	if (class(t_S)=="character")
		if(!(t_S %in% c("mst", "complexity", "n_cliques")))
			stop(paste0("If t_S is a character string it has to be one of:\n",
				"'mst','complexity','n_cliques'"))
	
	if (class(t_CS)=="character")
		if (!(t_CS %in% c("mst", "modularity")))
			stop("if t_CS is a character string it has be either 'modularity' or 'mst'")

	if (class(t_S)=="numeric")
		stopifnot(!is.na(t_S))	
	
	if (class(t_CS)=="numeric")
		stopifnot(!is.na(t_CS))	

	if (!is.null(X_simil_fun)) {
	stopifnot(class(X_simil_fun) %in% c("character", "function"))
	

	if (class(X_simil_fun)=="character") {
		if(!(X_simil_fun %in% c("euclidean", "maximum", 
					"manhattan", "canberra", 
					 "binary") ) )
			stop("If X_simil_fun is a character vector,
			      then it is assumed it is to be computed \n
			      based on distances between rows in XorS. \n
			      It has to be one of: \n
			      'euclidean', 'maximum', 'manhattan', 
			      'canberra', 'binary' or 'minkowski'")
		if (nrow(XorS)==ncol(XorS)) 
		warning("XorS is square but X_simil_fun was given \n
			 assuming XorS is coordinate matrix of points.
			 Computing similarities between rows of XorS.")
	   dist(XorS, method=X_simil_fun) %>% as.matrix() %>% flip_() ->S
	} else { #X_simil_fun is a function
	   X_simil_fun(XorS)-> S
	  }

	} else {
	S=XorS; 
	}
	 rm(XorS)
	thr_optimizer_ran=FALSE	
	if (class(t_S)== "character"){
		if (t_S=="mst") { 
		t_S = critical_mst_thr(S)
		S %thr% t_S -> S_t
		} else {
		# in that case t_S is "complexity" or "n_cliques"
		# so it goes to thr_optimizer()
		 thr_optimizer( S=S, expansion_mode=expansion_mode,
				cl_score=t_S,
				n_init_steps=t_S_init_steps,
				keep_init_eval_history= t_S_keep_init_eval_history,
				keep_init_partitions= t_S_keep_init_partitions,
				precis=t_S_precis,
				tol=t_S_tol,
				keep_final_partition=TRUE,
				) -> t_Sopt
		  t_S= t_Sopt$thr
		  S %thr% t_S -> S_t
		  thr_optimizer_ran=TRUE 
		}	
	} else if (class(t_S)=="function") { 
	    t_S=t_S(S)	
	    S %thr% t_S -> S_t	
	} else { # must be numeric then
	   S %thr% t_S -> S_t
	}
      	 	
	# partition thresholded S
	Pa= greedyCliquePartitioner(S_t, expansion_mode=expansion_mode)$membership
	# unique labels per each singleton clique (not 0s)
	# + relabel to consecutive integers just in case
	Pa= tidyUpLabels(uniqueSingletonLabels(Pa))
	# clique similarity (CS) matrix
	CS<-cliqueSimilarity(cl_mem= Pa, WorA= S_t, do_signif=do_signif, lvl=lvl)
	relax_optimalization_done=FALSE

	if (class(t_CS)=="numeric") {
		C_i<-relax_cliques(partition=Pa, 
			      SorCS=CS,
			      frac=t_CS,
			      CS_given=TRUE,
			      relax_method=t_CS_relax_method)
	 } else if (class(t_CS)=="character") {
		if (t_CS=="mst")
		relax_cliques_optimally(partition=Pa,
					S=S,
					t_S=t_S,
					CS=CS,
					relax_method=t_CS_relax_method,
					#clq_importance=t_CS_clq_importance,
					cs_thr_objective=t_CS,
					n_init_steps=t_CS_n_init_steps,
					keep_init_eval_history=t_CS_keep_init_eval_history,
					keep_init_partitions=t_CS_keep_init_partitions,
					precis=t_CS_precis,
					tol=t_CS_tol,
					dX.Y=t_CS_dX.Y
					) -> t_CSopt
		else
		relax_cliques_optimally(partition=Pa,
					S=S,
					t_S=t_S,
					CS=CS,
					relax_method=t_CS_relax_method,
					#clq_importance=t_CS_clq_importance,
					cs_thr_objective=t_CS,
					n_init_steps=t_CS_n_init_steps,
					keep_init_eval_history=t_CS_keep_init_eval_history,
					keep_init_partitions=t_CS_keep_init_partitions,
					precis=t_CS_precis,
					tol=t_CS_tol
					) -> t_CSopt
		 C_i <- t_CSopt$membership
		 t_CS= t_CSopt$frac	
		 relax_optimalization_done=TRUE
	}




	result=list(clique_membership=Pa,
		    t_S= t_S,
		    t_CS=t_CS,
		    cluster_membership=C_i
		   )
	if (return_S) result$S=S
	if (thr_optimizer_ran) {
		result$t_S_objective=t_Sopt$objective
		if( (!is.null(t_S_init_steps)) && t_S_keep_init_partitions )
			result$t_S_init_partitions= t_Sopt$init_partitions
		if( (!is.null(t_S_init_steps)) && t_S_keep_init_eval_history) {
			result$t_S_init_search_points= t_Sopt$init_search_points
			result$t_S_init_scores = t_Sopt$init_scores
	}

	}

	if (return_CS) result$CS=CS

	if (relax_optimalization_done) {
		 result$t_CS_objective= t_CSopt$objective
		if( (!is.null(t_CS_n_init_steps)) && t_CS_keep_init_partitions )
			result$t_CS_init_partitions= t_CSopt$init_partitions
		if( (!is.null(t_CS_n_init_steps)) && t_CS_keep_init_eval_history) {
			result$t_CS_init_search_points= t_CSopt$init_search_points
			result$t_CS_init_scores = t_CSopt$init_scores
	}
		
	}
		


	return(result)
}
		     
		     
		     
		     
		      
		      
