#' Example dataset from ecology to illustrate fitting partially additive GAMs
#'
#' Original dataset was copied from 
#' https://stackoverflow.com/questions/37037445/using-mob-trees-partykit-package-with-nls-model 
#' For illustrating the funtionality of package gamtree, two noise variables 
#' were added: \code{noise} (a continuous variable) and `cluster_id` (a 
#' categorical indicator for a random intercept). 
#'
#' \itemize{
#'   \item Species: a categorical covariate, indicator for specie (signal).
#'   \item PAR: a continuous covariate (signal).
#'   \item Pn: continuous response variable (signal).
#'   \item cluster_id: a categorical covariate; indicator for random intercept
#'   term (noise).
#'   \item noise: a continuous covariate (noise).
#' }
#'
#' @name eco
#' @usage data(eco)
#' @docType data
#' @format A data frame with 628 observations and 5 variables
#' @keywords datasets
#' @examples data("eco")
#' summary(eco)
#' 
NULL