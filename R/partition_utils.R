#' Check if partition Q subdivides partition P
#'
#' @param P membership vector
#' @param Q membership vector
#' @return `TRUE` if `Q` subdivides `P`, `FALSE` if there are overlaps
#' @export
`%subdivides%`<- function(Q,P) {
stopifnot(length(P)==length(Q))
cbind(P,Q)-> pairs
n_unique_pairs<- sum(!duplicated(pairs))
n_unique_Q<- length(unique(Q))
n_unique_pairs == n_unique_Q
}
