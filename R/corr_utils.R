

#' Minimal optimized computation of Pearson correlation by matrix product
#'
#' This function will generally perform faster than base R \code{cor}, provided one uses OpenBLAS matrix multiplication backend. 
#'
#' @param D Must be an object of a subclass of \code{matrix}. Rows store samples and columns - variables. No \code{NA}s allowed.
#' @return A numeric square matrix of Pearson correlation coefficient compute between each pair of varables in \code{D}.
#' @examples
#' data(brca)
#' C<- corfast(brca)
#' @importFrom matrixStats colSds
#' @seealso [corpval()] and [fastPearsonData()] for p-values, \href{https://rmflight.github.io/posts/2023-11-01-installing-openblas-for-selfcompiled-r}{this link} [1] for how to set up OpenBLAS on Linux, 
#'	and \href{https://github.com/david-cortes/R-openblas-in-windows}{that one} for instructions for Windows.
#' @references [1]  M Flight, Robert. 2023. “Installing OpenBLAS for Self-Compiled R.” November 1, 2023. <https://rmflight.github.io/posts/2023-11-01-installing-openblas-for-selfcompiled-r>
#' @export

corfast<- function(D){    # n x p
stopifnot(any(class(D) %in% c('matrix')))   #to be safe, input has to be matrix
D_means<- colMeans(D)
# matrixStats required for efficient SD calculation
D_sds<-colSds(D)
# In R, if we subtract MATRIX - VECTOR, VECTOR gets recycled column wise to match MATRIX,
#  we can exploit this, hence t(D)
tZ<-  (t(D) - D_means)/( D_sds )   # p  x n
 (1/(nrow(D)-1))*(tZ %*% t(tZ))   # p x p
}


#' Minimal optimized computation of a p-value of a Pearson correlation coefficient
#'
#' This function will generally perform faster than base R \code{cor.test} but assumes much about the cleaniness of the data input.
#'
#' @param r A vector of correlation coefficients for which to compute p-values.
#' @param n A number of samples used to compute each statistic in \code{r}. Assumed to be the same for all entries of \code{r}.
#' @return A vector of corresponding 2-tailed p-values of statistics in \code{r}, where under the null hypothesis \code{r} assumes \code{T} distribution with \code{n-2} degrees of freedom.
#' @examples
#' data(brca)
#' C<- corfast(brca)
#' pv<- corpval(C[lower.tri(C)], nrow(brca))
#' @seealso [corfast()] and [fastPearsonData()] for fast correlation coefficient calculation
#' @export

corpval<- function(r,n){
# two tailed
t_stat= sqrt(n-2) * r/sqrt(1 - r^2)
2*pmin( pt( t_stat, n-2),
	pt(t_stat, n-2, lower.tail=FALSE ) ) 
}

#' Generate list storing all pairwise Pearson correlation coefficients between all variables in the dataset, and their (adjusted) 2-tailed p-values using \code{corfast}
#'
#' Function will run very fast if efficient matrix multiplication backed is used (like OpenBLAS). However, no input correctness checks are performed (apart from warning on presence of \code{NA} in input).
#'
#' @param D Must be an object of a subclass of \code{matrix}. Rows store samples and columns - variables. No checks on 0 variance are performed.
#' @param NA.warn Logical, if true (default), then a warning will be displayed if \code{D} contains any \code{NA} values.
#' @param p.adjust.method The method of p-value adjustment used. For \code{D} with \code{N} columns, number of tests is \code{N*(N-1)*0.5}. Gets passed to argument \code{p.adjust.method} in a call of \code{p.adjust} function.
#' @return A list storing all the important data about correlations in \code{D} with the following fields:
#' \itemize{
#' \item \code{N} Number of variables in \code{D}.
#' \item \code{r} A vector of all correlation coefficients between all pairs of variables in \code{D}.
#' \item \code{P} A vector of unadjusted p-values of \code{r}.
#' \item \code{Pa} A vector of adjusted p-values of \code{r}.
#' }
#' @examples
#' data(brca)
#' fastPearsonData(brca) -> brca_corData
#' C<- matrix(nrow= brca_corData$N, ncol= brca_corData$N) # restore matrices from vectors in corData like this:
#' C[lower.tri(C)]= brca_corData$r
#' C<-t(C)
#' C[lower.tri(C)]= brca_corData$r
#' diag(C)=1
#' @seealso [corfast()], [corpval()] for individual component functions, [rcorrData()] for alternative using \code{rcorr} package,
#'	[similarity_matrix()] for creating similarity matrix out of output of this function
#' @export

fastPearsonData<- function(D, NA.warn=TRUE, p.adjust.method="holm") {

if (NA.warn) if (any(is.na(D))) warn("warning: D contains NAs. Consider using rcorrData")
r= corfast(D)
ltr<- lower.tri(r)
corData<- list()
corData$N= ncol(D)
corData$r = r[ltr]
corData$P = corpval( corData$r , nrow(D))
corData$Pa= p.adjust( corData$P ,method=p.adjust.method  )
corData
}

#' Generate list storing all pairwise Pearson correlation coefficients between all variables in the dataset, and their (adjusted) 2-tailed p-values using \code{Hmisc::rcorr}
#'
#' Function utilizes \code{rcorr} from \code{Hmisc} package. It is less speed and memory efficient than \code{fastPearsonData} but is safer to use if input data is not properly cleaned.
#'
#' @param D Must be an object of a subclass of \code{matrix}. Rows store samples and columns - variables. No checks on 0 variance are performed.
#' @param NA.warn Logical, if true (default), then a warning will be displayed if \code{D} contains any \code{NA} values.
#' @param p.adjust.method The method of p-value adjustment used. For \code{D} with \code{N} columns, number of tests is \code{N*(N-1)*0.5}. Gets passed to argument \code{p.adjust.method} in a call of \code{p.adjust} function.
#' @return A list storing all the important data about correlations in \code{D} with the following fields:
#' \itemize{
#' \item \code{N} Number of variables in \code{D}.
#' \item \code{r} A vector of all correlation coefficients between all pairs of variables in \code{D}.
#' \item \code{P} A vector of unadjusted p-values of \code{r}.
#' \item \code{Pa} A vector of adjuste p-values of \code{r}.
#'}
#' @examples
#' data(brca)
#' rcorrData(brca) -> brca_corData
#' C<- matrix(nrow= brca_corData$N, ncol= brca_corData$N) # restore matrices from vectors in corData like this:
#' C[lower.tri(C)]= brca_corData$r
#' C<-t(C)
#' C[lower.tri(C)]= brca_corData$r
#' diag(C)=1
#' @importFrom Hmisc rcorr
#' @seealso [fastPearsonData()] for efficient version of correlation coefficient calculation for clean data, [Hmisc::rcorr()] for correlation function used here,
#'	[similarity_matrix()] for creating similarity matrix out of output of this function
#' @export

rcorrData<- function(D, rcorr.method='pearson',
			p.adjust.method='holm'){
rcorr_object<- rcorr(D, type= rcorr.method)
r<-rcorr_object$r
ltr<-lower.tri(r)
corData<- list()
corData$N= ncol(D)
corData$r = rcorr_object$r[ltr]
corData$P = rcorr_object$P[ltr]
corData$Pa= p.adjust( corData$P, method=p.adjust.method  )
corData$n=nrow(D)
corData
}

#' Compute similarity matrix from input correlation coefficients
#'
#' Function takes compressed data about symmetric correlation matrix and builds a \code{WGCNA} style similarity matrix out of it, raising elements to \code{power} and zeroing out unsignificant correlations.
#' 
#' @param corData A list storing lower triangular halves of correlation and p-value matrices, see output of \code{fastPearsonData} or \code{rcorrData} for how it should be structured.
#' @param addLoops Logical, if \code{TRUE}, diagonal of the similarity matrix will assume values of \code{1} nad \code{0} otherwise.
#' @param power A power to which raise the absolute value of the correlation coefficients.
#' @param level A significance level of adjusted p-values over which similarities are set to zero.
#' @param return.lower.tri  A logical, if \code{TRUE} then function returns lower triangular half of the similarity matrix. Otherwise, full (symmetric) matrix of similarities is returned.
#' @return Depending on the value of \code{return.lower.tri}, either a square, symmetric matrix of nonnegative similarities between objects ( \code{S[i,j]= abs(cor(X_i,X_j))^power} if \code{cor(X_i,X_j)} is significant at \code{level} after p-value adjustment), or a lower triangular half of that matrix (\code{S[lower.tri(S)]})
#' @examples
#' data(brca)
#' brca_corData= fastPearsonData(brca)
#' S= similarity_matrix( brca_corData)
#' ncol(S)==ncol(brca)
#' @seealso  [fastPearsonData()] or [Hmisc::rcorr()] for functions creating suitable input \code{corData} argument.
#' @export

similarity_matrix<- function( corData, addLoops=TRUE, 
				power=2, level=0.05,
			     return.lower.tri=FALSE) {
if (!return.lower.tri){ S<- matrix(nrow=corData$N,
	   			ncol=corData$N)
if (addLoops) diag(S)=1 else diag(S)=0 
ltr= lower.tri(S) }
S_ltr= corData$r 
S_ltr[ corData$Pa >= level ] = 0
thr_used= if (!length(left_weights)) max(corData$r[ S_ltr==0 ] ) + min(1/sqrt(corData$n),1e-5) else  min(S_ltr[S_ltr!=0])  - min(1/sqrt(corData$n),1e-5)
attr(S_ltr,"thr")=thr_used 
S_ltr= abs(S_ltr)^power
if (!return.lower.tri) {
 S[ ltr ] = S_ltr
 S<- t(S)
 S[ ltr] = S_ltr
 return(S)
			} else {
 return(S_ltr)
	}
}


