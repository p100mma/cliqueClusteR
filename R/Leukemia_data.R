#' Leukemia dataset similarity matrix'
#'
#' Dataset was gathered from `FCPS` package [1]. Contains similarity matrix of 554 subjects calculated based on
#' gene expression profiles.
#'
#' Citing [2]: 
#' "The anonymized leukemia dataset consists of 12,692 gene expressions from 554 subjects (...)
#'  Of the subjects, 109 were healthy, 15 were diagnosed with acute promyelocytic leukemia (APL), 
#' 266 had chronic lymphocytic leukemia (CLL), and 164 had acute myeloid leukemia (AML).
#' The leukemia dataset was preprocessed, resulting in a high-dimensional dataset with 7747 variables 
#' and 554 data points separated into natural clusters as determined by the illness status and defined by the 
#' patterns of change in distance and density."
#'
#' @docType data
#'
#' @usage data(leukemia)
#'
#' @format A vector of integers encoding cluster labels of each gene. Zeroes denote clusterless genes.
#'
#' @keywords cluster_labels
#' 
#' @references [1] Michael Christoph, Stier Q (2021). “Fundamental clustering algorithms suite.” SoftwareX, 13, 100642. ISSN 2352-7110, \url{https://doi.org/10.3390/pr9101697} 
#' @references [2] Thrun, M. C., & Ultsch, A. (2020). Clustering benchmark datasets exploiting the fundamental clustering problems. Data in brief, 30, 105501.
#'
#' @seealso [leukemia_clusters] for the reference cluster labels
#' @examples 
#' data(leukemia)
#' dim(leukemia)
"leukemia"
