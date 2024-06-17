## TODO: Assign additional class to spline-based (g)lmertrees? splinetree?
## Note, there is already a package on rprart-based spline trees and RF called splinetree
## But it partitions based on spline coefficients and does not allow for a mixed-effects
## specification. It has function (and perhaps class) splineTree, so splinetree can still
## be used.


#' Set up splines bases for use with function (g)lmertree.
#' 
#' \code{setup.spline} takes a dataset and spline specification as input, and returns the 
#' dataset with spline bases added.
#' 
#' @param spline character vector of length 1, describing the spline basis to be created.
#' Currently, functions \code{ns} and \code{bs} are supported. See Examples.
#' @param data a \code{data.frame} containing the variable referred to in \code{spline}.
#' @param ... additional arguments to be passed to the function specified in \code{spline}.
#' @return a \code{data.frame} with a many rows and columns as \code{data}. The spline 
#' consists of \code{df} basis functions, but is contained in a single column named 
#' \code{spline}, followed by the name of the predictor variable specified.
#' @examples 
#' data <- setup.spline("ns(PAR, df = 3)", data = eco)
#' head(data)
#' matplot(x = data$PAR[order(data$PAR)], 
#'         y = data$spline.PAR[order(data$PAR),], type = "l")
#' data <- setup.spline("bs(PAR, degree = 2, df = 4)", data = eco)
#' head(data)
#' matplot(x = data$PAR[order(data$PAR)], 
#'         y = data$spline.PAR[order(data$PAR),], type = "l")
#' @seealso \code{\link{plot.splinetree}} \code{\link{predict.splinetree}} 
#' \code{\link[glmertree]{lmertree}} \code{\link[glmertree]{glmertree}}
#' @import splines glmertree
#' @export
setup.spline <- function(spline, data, ...) {
  
  data <- as.data.frame(data)
  data[paste0("spline.", names(which(sapply(names(data), grepl, spline))))] <- 
    with(data, eval(parse(text = as.character(spline))))
  return(data)
  
}
## Can create basis using mgcv as follows:
#sc <- smoothCon(s(PAR, bs = "tp", k = 7), data = eco)
#sc[[1L]]$X ## 2nd-to-last column seems to be the intercept always


## TODO: Assign additional class to spline-based (g)lmertrees? splinetree?
## Note, there is already a package on rprart-based spline trees and RF called splinetree
## But it partitions based on spline coefficients and does not allow for a mixed-effects
## specification.


#' Fit a (g)lmertree using spline-based partitioning.
#' 
#' \code{splinetree} is a wrapper for functions \code{(g)lmertree} to simplify
#' fitting, visualizing and predicting spline-based trees.
#' 
#' @param formula A four-part function See Examples below, and
#' \code{\link[glmertree]{lmertree}} or \code{\link[glmertree]{glmertree}}.
#' @param data a \code{data.frame} containing the variable referred to in \code{spline}.
#' @param family family specification for \code{\link[glmertree]{glmertree}}. See 
#' \code{\link[stats]{glm}} documentation for families.
#' @param ... additional arguments to be passed to function \code{\link[glmertree]{lmertree}} 
#' (default, i.e., \code{family = "gaussian"})
#' or \code{\link[glmertree]{glmertree}} (\code{family} other than \code{gaussian}).
#' @return A object of class\code{"splinetree" }and \code{"lmertree"} or \code{"glmertree"}.
#' @examples  sp <- splinetree(Pn ~ ns(PAR, df = 5) | Specimen | Species, data = eco)
#' sp
#' @seealso \code{\link{plot.splinetree}} \code{\link{predict.splinetree}} 
#' \code{\link[glmertree]{lmertree}} \code{\link[glmertree]{glmertree}}
#' @export
splinetree <- function(formula, data, family = "gaussian", ...) {

  ## prepare formula and data
  ff <- as.Formula(formula)
  spline <- as.character(formula(ff, lhs = 0, rhs = 1))[2L]
  data <- setup.spline(spline, data = data)
  upd <- paste0(". ~ spline.", all.vars(ff)[2L], "| . | .")
  formula <- update(ff, as.Formula(upd))
  
  ## fit tree
  tree <- if (family == "gaussian") {
    lmertree(formula = formula, data = data)
  } else {
    glmertree(formula = formula, data = data, family = family, ...)
  }
    
  ## return result
  class(tree) <- c("splinetree", class(tree))
  return(tree)
}



#' Plotting function for visualization of spline-based (g)lmertrees.
#' 
#' \code{plot.splinetree} takes a fitted \code{(g)lmertree} with splines and plots it. 
#' It is a wrapper for \code{\link[glmertree]{plot.glmertree}} and  
#' \code{\link[glmertree]{plot.lmertree}}, with critical adjustments for
#' better visualization of spline models in the terminal nodes.
#' 
#' @param x fitted object of class \code{(g)lmertree} containing splines specified
#' by \code{\link{setup.spline}} in the terminal node model.
#' @param ... additional arguments to be passed to \code{\link[glmertree]{plot.lmertree}}
#' or \code{\link[glmertree]{plot.glmertree}}.
#' @examples
#' st <- splinetree(Pn ~ ns(PAR, df = 5) | Specimen | Species, data = eco, 
#'                cluster = Specimen)
#' plot(st)
#' @seealso \code{\link{setup.spline}} \code{\link{predict.splinetree}} 
#' \code{\link[glmertree]{plot.lmertree}} \code{\link[glmertree]{plot.glmertree}}
#' @export
plot.splinetree <- function(x, ...) {
  
  ## Replace spline basis functions with original x
  spline_name <- names(x$tree$data)[grep("spline", names(x$tree$data))]
  var_name <- substr(spline_name, 8L, stop = 10000)
  x$tree$data[ , spline_name] <- x$data[, var_name]
  
  ## Replace fitted with effects of node-specific splines 
  tree_node <- as.list(x$tree$node)
  for (i in unique(x$tree$fitted[["(fitted)"]])) {
    X <- cbind(1L, x$data[spline_name][x$tree$fitted["(fitted)"] == i, ])
    b <- tree_node[[i]]$info$coefficients#[grep("spline", names(tree_node[[i]]$info$coefficients))]
    tree_node[[i]]$info$object$fitted.values <- X%*%b
  }
  x$tree$node <- as.partynode(tree_node)
  class(x) <- ifelse(inherits(x, "lmertree"), "lmertree", "glmertree")
  plot(x$tree, ...)
  
}

#' Predict method for spline-based (g)lmertrees.
#' 
#' \code{predict.splinetree} computes predictions for a fitted 
#' \code{(g)lmertree} that is based on splines.
#' 
#' @param object fitted object of class \code{(g)lmertree} containing splines specified
#' by \code{\link{setup.spline}} in the terminal node model.
#' @param newdata \code{data.frame} with observations for which predictions should be computed.
#' @param ... additional arguments to be passed to \code{\link[glmertree]{predict.lmertree}}
#' or \code{\link[glmertree]{predict.glmertree}}.
#' @examples
#' st <- splinetree(Pn ~  ns(PAR, df = 5) | Specimen | Species, data = eco, 
#'                cluster = Specimen)
#' predict(st, newdata = eco[1L:5L, ])
#' @seealso \code{\link{setup.spline}} \code{\link{predict.splinetree}} 
#' \code{\link[glmertree]{lmertree}} \code{\link[glmertree]{glmertree}}
#' @export
predict.splinetree <- function(object, newdata, ...) {
  
  ## Set up spline basis functions for newdata  
  ff <- as.Formula(object)
  xname <- names(which(sapply(names(newdata), grepl, 
                              as.character(formula(ff, lhs = 0, rhs = 1))[2])))
  newdata[paste0("spline.", xname)] <- predict(object$data$spline, newx = newdata[ , xname])
  
  ## Compute predictions
  class(object) <- ifelse(inherits(object, "lmertree"), "lmertree", "glmertree")
  predict(object, newdata, ...)    
}
