% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clique_cluste.R
\name{MinST_DunnIndex}
\alias{MinST_DunnIndex}
\title{Calculate Dunn Index goodness of clustering variant based on internal MST of clusters}
\usage{
MinST_DunnIndex(
  partition,
  dS,
  dX.Y = "single_linkage",
  hausdorff_q = 1,
  outlier_mask = NULL
)
}
\arguments{
\item{partition}{A vector of membership labels of clusters of each node, 0 labels get converted to singleton labels by \code{\link[cliquePartitioneR:uniqueSingletonLabels]{cliquePartitioneR::uniqueSingletonLabels()}}}

\item{dS}{A dissimilarity matrix of nodes in the graph. The higher the \code{dS[i,j]} the more dissimilar nodes \code{i} and \code{j} are.}

\item{dX.Y}{A character string specifying type of between cluster distance, see details.}

\item{outlier_mask}{(experimental) if not null, if \code{outlier_mask[i]} is \code{TRUE} then node \code{i} is not included in the calculation, see details.}
}
\value{
A scalar value indicating the goodness of clustering, the bigger the better. Can be \code{Inf} if max cluster diameter is 0. If both \code{min_between_cluster_distance} and \code{max_cluster_diameter} are 0, function returns 0 by convention.
}
\description{
For a given dissimilarity graph \code{dS} and partition of \code{dS} into
clusters, function will compute Dunn Index considering
minimum spanning tree (MST) of each cluster to get diameter
and optionally a Hausdorff distance for between cluster distances instead of
original single linkage distance, see details.
}
\details{
In general, Dunn Index is defined as \code{min_between_cluster_distance/max_cluster_diameter} - the bigger
the value, the better the clustering.

Note that it is formulated in terms of distances, not in similarities, so to keep interpretation consistent
with literature this function works with dissimilarities. Below by distance we mean dissimilarity.
Dissimilarity can be derived from similarity matrix for example by \code{dS= max(S) - S}.

As first proposed by Dunn \link{ADD REFERENCE}, \code{between_cluster_distance}
between clusters \verb{i,j} as a minimum distance between any nodes inside clusters \verb{i,j}, like in
single linkage method of hierarchical clustering.
If argument \code{dX.Y} is set to \code{"hausdorff"}, then \code{between_cluster_distance} changes to
Hausdorff distance:
\itemize{
\item \verb{d_H(X,Y):= max ( max_x d(x,Y) , max_y d(y,X) )}
\item where \verb{X,Y} are clusters, \code{x} points in \code{X} and \code{y} points in \code{Y},
\item \verb{d(x,Y):= min_y d(x,y)}.
}
This makes the index less influenced by few points of contact between clusters if in general they are
well separated.

And \code{cluster_diameter} is defined here as maximum distance on MST of a cluster.
Such choices make the index suited for finding elongated clusters in which distances only between
adjacent neighbors are small but distances between indirect neighbors are allowed ot be big.

Value of the index is undefined if there is only one cluster or all clusters are singletons, in that case
0 is returned by convention.
For the case when maximum cluster diameter is 0 value of the index is infinity.
and if also minimum between cluster distance is 0 then the function returns 0 (there are no differences
between inter and intra cluster distances so a given clustering is not correct).

If \code{outlier_mask} is not \code{NULL}, then it is assumed that it is a logcial vector of length equal to number of nodes in the network.
If \code{outlier_mask[[i]]} is \code{TRUE} then node \code{i} gets removed from the network before calculating Dunn Index.
This is an ad hoc method of masking some of the particularly problematic nodes that connect clusters
by very few but strong connections (lowering the \code{between_cluster_distance}).
Criteria for being an outlier here can be based on low clustering coefficient (transitivity), see \code{\link[igraph:transitivity]{igraph::transitivity()}}.

Note that even though package offers various thresholding schemes for similarity matrices, for assesing
goodness of the clustering one should calculate this based on dissimlarity derived from the unthresholded version
of the graph.
}
\examples{
library(magrittr)
points<- matrix(rnorm(800*3), ncol=3)
points<-points/sqrt(rowSums(points^2)) #unit sphere
points[1:400,]= points[1:400,]*20 # outer
points[401:800,]= points[401:800,]*0.5 #inner
ref=rep(1,800); ref[401:800]=2 # reference labels for points
plot(points, col=ifelse(ref>1,"red","blue"))
dS<- as.matrix(dist(points))
flip_<- function(M) max(M)- M #convert similarities to distances and vice versa
dS \%>\% flip_() -> S
thr<- critical_mst_thr(S)
# first get cliques
pa<- greedyCliquePartitioner(S \%thr\% thr)$membership
#not perfect yet
igraph::compare(pa, ref, "adjusted.rand")
#corresponding dunn indices:
MinST_DunnIndex(ref,dS) #reference
MinST_DunnIndex(pa,dS) #cliques
print("hausdorff:")
MinST_DunnIndex(ref,dS, "hausdorff") #reference
MinST_DunnIndex(pa,dS, "hausdorff") #cliques
# clusters from clique relaxation
relax_cliques(pa \%>\% uniqueSingletonLabels(), S \%thr\% thr, relax_method = "components")-> cl
igraph::compare(ref, cl$node, "adjusted.rand")
MinST_DunnIndex(cl$node,dS)
MinST_DunnIndex(cl$node,dS,"hausdorff")
#fine tune frac of kept edges between cliques, better score:
relax_cliques(pa \%>\% uniqueSingletonLabels(), S \%thr\% thr, relax_method = "components", frac=0.1)-> cl
igraph::compare(ref, cl$node, "adjusted.rand")
#indexes are better
MinST_DunnIndex(cl$node,dS)
MinST_DunnIndex(cl$node,dS,"hausdorff")
}
\seealso{
\link{\%thr\%}, \code{\link[igraph:modularity.igraph]{igraph::modularity()}}, \code{\link[=relax_cliques]{relax_cliques()}}, \code{\link[=critical_mst_thr]{critical_mst_thr()}}
}
