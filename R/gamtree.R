utils::globalVariables(c(".tree", ".offset", ".global", ".weights", ".cluster", "current_offset"))

## TODO: Can enforce knot locations to apply to all nodes similarly? Yes, but only with bs = "cr"
## basis, and knots have to be passed to function gam / gamm4.
## C

## TODO: Can use a cubic spline basis in gamtree? Yes, just pass bs = "cr" to function s.

## TODO: Get the node-specific plots with CIs from mgcv and put them together in one plot



#' Recursively partition a dataset based on penalized GAMs.
#' 
#' \code{gamtree} recursively partitions a dataset into subgroups with 
#' penalized GAMs, characterized by differences in the parameter estimates. 
#' 
#' @param formula specifies the model formula, consisting of three 
#' parts: the response variable followed by a tilde ('~'); the terms for the 
#' node-specific GAMs, followed by a vertical bar ('|') and the potential 
#' partitioning variables (separated by a '+'). 
#' The 'by' argument of function
#' \code{\link[mgcv]{s}} may NOT be used in the node-specific GAM formulation. 
#' Refrain from
#' using the dot ('.') to specify all remaining variables in \code{data}, this may yield
#' unexpected results; make sure to specify each variable in the corresponding part
#' of the model formula. See Examples.
#' @param data \code{data.frame} containing the variables specified in \code{formula}.
#' @param weights numeric vector of length \code{nrow(data)}; optional case weights.
#' A weight of 2, for example, is equivalent to having made exactly the same 
#' observation twice.
#' @param offset numeric vector of length \code{nrow(data)}. Supplies model 
#' offset for use in fitting. Note that this offset will always be completely 
#' ignored when predicting.
#' @param cluster optional numeric or factor vector with a cluster ID to be 
#' employed for clustered covariances in the parameter stability tests.
#' This argument should be used when the partitioning variables are not measured
#' on the individual observation level, but on a higher level. E.g., when 
#' the response variables consists of repeated measurements on multiple 
#' @param verbose logical. Should progress be printed to the commande line in 
#' every iteration? If true, the iteration number, information on the 
#' splitting procedure, and the log-likelihood (with df) value of the fitted 
#' full mixed-effects gam model is printed.
#' @param alt_formula list with two elements, for specifying non-standard model formulae
#' for GAM. E.g., the formula list required for use of the \code{\link[mgcv]{multinom}}
#' family.
#' @param gam_ctrl a list of fit control parameters to replace defaults returned by 
#' \code{\link[mgcv]{gam.control}}. Values not set assume default values.
#' @param mob_ctrl a list with control parameters as returned by 
#' \code{\link[partykit]{mob_control}} to be passed to function 
#' \code{\link[partykit]{mob}}. Note: argument `xtype` is set to 
#' `data.frame`, by default, and cannot be adjusted.
#' @param ... additional arguments to be passed to function \code{\link[mgcv]{gam}}. 
#' The same arguments will be employed for the local
#' (node-specific) GAMs, and the global GAM (if specified in the model \code{formula}). 
#' @return Returns an object of class \code{"gamtree"}. This is a list, containing
#' (amongst others) the GAM-based recursive partition (in \code{$tree}), and the 
#' fitted full GAM model with local and/or global fitted effects (in \code{$gamm}. 
#' The following methods are available to extract information from the fitted object:
#' \code{\link{predict.gamtree}}, for obtaining predicted values for training and new
#' observations; \code{\link{plot.gamtree}} for plotting the tree and variables' effects; 
#' \code{\link{coef.gamtree}} for obtaining the estimated coefficients; 
#' \code{\link{summary.gamtree}} for a summary of the fitted model.
#' 
#' @examples
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$Specimen)
#' summary(gt)
#' 
#' @seealso \code{\link{predict.gamtree}} \code{\link{plot.gamtree}} 
#' \code{\link{coef.gamtree}} \code{\link{summary.gamtree}}
#' @import mgcv partykit gamm4 merDeriv
#' @importFrom stats as.formula formula logLik predict update terms complete.cases
#' @importFrom lme4 fixef VarCorr
#' @importFrom Formula as.Formula
#' @export
gamtree <- function(formula, data, weights = NULL, cluster = NULL, 
                    offset = NULL, verbose = FALSE, 
                    gam_ctrl = list(), 
                    mob_ctrl = mob_control(parm = c(1,2,4)), 
                    alt_formula = NULL, ...) {

  if (!inherits(data, "data.frame")) warning("Argument data should specify a data.frame")
  
  if (length(gam_ctrl) == 0L) gam_ctrl <- NULL
  ## Construct formulas for tree (tf), local gam (lgf) and global gam (ggf):
  ##  - tf is for passing response, predictors and part vars to mob()
  ##  - lgf is for fitting node-specific GAM in the nodes
  ##  - ggf is an intermediate formula, used for creating full GAM formula later
  ##  - fgf is for fitting full GAM
  ff <- as.Formula(formula)
  lgf <- formula(ff, lhs = NULL, rhs = 1)
  # TODO: Test if by-statements problematic for partitioning  
  local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
  part_vars <- all.vars(formula(ff, lhs = 0, rhs = 2))
  if (length(local_vars) == 0) local_vars <- 1L
  response <- ff[[2]] 
  if (length(response) > 1L && is.call(response[1])) {
    response <- as.character(response)
    response <- paste0(response[1L], "(", response[2L], ", ", response[3L], ")")
  }
  tf <- formula(paste(response, "~", 
                paste(local_vars, collapse = " + "), "|", 
                paste(part_vars, collapse = " + ")))
  
  q_cluster <- substitute(cluster)
  if (!is.null(q_cluster)) {
    data$.cluster <- eval(q_cluster, data)
    if (length(eval(q_cluster, data)) != nrow(data))
      warning("Variable lengths differ for 'cluster' and 'data'.", immediate. = TRUE)
  }
  if (!is.null(q_cluster) && !inherits(data$.cluster, c("numeric", "character", "factor", "integer"))) {
    warning("Argument 'cluster' should specify an object of class numeric, factor or character, or should be NULL.", immediate. = TRUE)
  }
  
  ## Prepare data
  data <- if (is.null(q_cluster)) {
    data[ , c(as.character(response), local_vars, part_vars)]
  } else {
    data[ , c(as.character(response), local_vars, part_vars, ".cluster")]
  }
  N <- nrow(data)
  data <- data[complete.cases(data), ]
  if (nrow(data) != N) {
    warning(paste(N - nrow(data), "observations were removed due to missing values."))
    N <- nrow(data)
  }
  data$.weights <- if (is.null(weights)) 1 else weights

  if (is.null(mob_ctrl)) {
    mob_ctrl <- mob_control(verbose = verbose, xtype = "data.frame",
                            ytype = "data.frame")
  } else {
    mob_ctrl <- do.call("mob_control", args = mob_ctrl)
    mob_ctrl$verbose <- verbose
    mob_ctrl$ytype <- mob_ctrl$xtype <- "data.frame"
  }
  
  ## remember call
  cl <- match.call()
  
  ## initialization
  data$.offset <- if (is.null(offset)) 0L else offset
  
  ## grow tree
  if (is.null(q_cluster)) {
    tree <- mob(tf, data = data, local_gam_form = lgf, fit = gamfit, 
                weights = .weights, offset = .offset, 
                control = mob_ctrl, gam_ctrl = gam_ctrl, ...)   
  } else {
    tree <- mob(tf, data = data, local_gam_form = lgf, fit = gamfit, 
                weights = .weights, offset = .offset, 
                control = mob_ctrl, cluster = .cluster, gam_ctrl = gam_ctrl, ...)       
  }

  ## collect results
  result <- list(
    formula = formula,
    data = data,
    call = cl,
    tree = tree
  )
  class(result) <- "gamtree"
  return(result)
}



##############################################################################
##
## gamm4-based fitting function
##
gamfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, 
                   estfun = NULL, object = NULL, ..., local_gam_form,
                   gam_ctrl = list()) {
  
  ## Prepare arguments to be passed to fitting function
  args <- list(...)
  if (is.null(x)) { 
    x <- matrix(1, nrow = NROW(y), ncol = 1L,
                dimnames = list( NULL, "(Intercept)") )
  }
  args <- c(list(formula = local_gam_form, data = cbind(x, y), weights = weights, 
                 offset = offset, control = gam_ctrl), args)
  
  ## Fit model
  object <- do.call("gamm4", args) 
  class(object) <- "gamm4"
  
  ## Return results
  list(object = object,
    coefficients = list(fixef = fixef(object$mer), ranef = VarCorr(object$mer)),
    objfun = logLik(object$mer), # minimized objective function
    estfun = estfun.lmerMod(object$mer, level = 1) # empir. estimating functions / scores
  )
}




##############################################################################
##
## mgcv-based fitting function for gam-based recursive partitioning with MOB.
##
# gamfit <- function (y, x, start = NULL, weights = NULL, offset = NULL, 
#                     ..., local_gam_form, gam_ctrl = list())
# {
#   
#   args <- list(...)
#   
#   if (is.null(x)) 
#     x <- matrix(1, nrow = NROW(y), ncol = 1L, 
#                 dimnames = list(NULL, "(Intercept)"))
#   
#   args <- c(list(formula = local_gam_form, data = cbind(x, y), 
#                  weights = weights, offset = offset, control = gam_ctrl), args)
#   
#   do.call("gam", args)
#   
#   ## TODO: in ?gamObject we read that deviance is "model deviance (not penalized deviance)".
#   ##
#   ## TODO: Estimation of par stab tests. From mob() documentation:
#   ## All tests also require an estimate of the corresponding variance-covariance
#   ## matrix of the estimating functions. The default is to compute this using an
#   ## outer-product-of-gradients (OPG) estimator. Alternatively, the corrsponding
#   ## information matrix or sandwich matrix can be used if: (a) the estimating
#   ## functions are actually maximum likelihood scores, and (b) a vcov() method
#   ## (based on an estimate of the information) is provided for the fitted model
#   ## objects. The corresponding option in mob_control() is to set vcov ="info"
#   ## or vcov="sandwich"rather than vcov="opg" (the default).
#   ## mgcv supplies a vcov.gam method, so condition b) is met, a) also met?
# }




#' Print method for a fitted GAM tree
#' 
#' Prints the local and/or global terms in a fitted GAM tree.  
#'
#' @param x object of class \code{gamtree}.
#' @param ... further arguments to be passed to \code{\link[partykit]{print.modelparty}}.
#' @export
#' @examples
#' ## GAM tree without global terms:
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$Specimen)
#' gt ## or: print(gt)
print.gamtree <- function(x, ...) {
  
  tree_node <- as.list(x$tree$node)
  for (i in unique(x$tree$fitted[["(fitted)"]])) {
    ## To get the correct (or at least one) model in the terminal nodes
    tree_node[[i]]$info$coefficients <- print_gam(tree_node[[i]]$info$object$gam)
  }
  x$tree$node <- as.partynode(tree_node)
  print(x$tree, ...)
  
}

print_gam <- function (x, ...) {
  cat("\nFormula:\n")
  if (is.list(x$formula)) 
    for (i in 1:length(x$formula)) print(x$formula[[i]]) else print(x$formula)
  n.smooth <- length(x$smooth)
  if (n.smooth == 0) 
    cat("Total model degrees of freedom", sum(x$edf), "\n") else {
      edf <- 0
      cat("\nEstimated degrees of freedom:\n")gt
      for (i in 1:n.smooth) edf[i] <- sum(x$edf[x$smooth[[i]]$first.para:x$smooth[[i]]$last.para])
      edf.str <- format(round(edf, digits = 4), digits = 3, scientific = FALSE)
      for (i in 1:n.smooth) {
        cat(edf.str[i], " ", sep = "")
        if (i%%7 == 0) cat("\n")
      }
      cat(" total =", round(sum(x$edf), digits = 2), "\n")
    }
  if (!is.null(x$method) && !(x$method %in% c("PQL", "lme.ML", "lme.REML"))) 
    cat("\n", x$method, " score: ", x$gcv.ubre, "     ", sep = "")
  if (!is.null(x$rank) && x$rank < length(x$coefficients)) 
    cat("rank: ", x$rank, "/", length(x$coefficients), sep = "")
  cat("\n")
  invisible(x)
}


#' Summary method for a fitted GAM tree
#' 
#' Prints a summary of the local and/or global terms in a fitted GAM tree.  
#'
#' @param object object of class \code{gamtree}.
#' @param ... further arguments to be passed to \code{summary.gam}.
#' @section Warning:
#' The printed results by default also provide standard error and significance
#' tests. These should be taken with a big grain of salt, because they do NOT 
#' account for the searching of the tree structure; they assume the tree structure 
#' was known in advance. They thus should be interpreted as overly optimistic and
#' with caution. 
#' @export
#' @examples
#' ## GAM tree without global terms:
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$Specimen)
#' summary(gt)
summary.gamtree <- function(object, ...) {

  tree_node <- as.list(object$tree$node)
  for (i in unique(object$tree$fitted[["(fitted)"]])) {
    ## To get the correct (or at least one) model in the terminal nodes
    tree_node[[i]]$info$coefficients <- paste0("node ", i, print_gam(tree_node[[i]]$info$object$gam))
  }
  object <- as.partynode(tree_node)
  summary(object$tree, ...)
  
}

#' Plotting method for GAM trees
#' 
#' Takes a fitted GAM tree and plots the smooth functions fitted in each of the
#' terminal nodes of the tree.
#' 
#' @param x object of class \code{gamtree}.
#' @param which character. The default (\code{"both"}) plots the tree structure, 
#' followed by the model fitted in the terminal nodes. Alternatively, \code{"tree"}
#' will plot the tree structure, and \code{"terms"} will plot the smooth (and 
#' parametric) terms from
#' the terminal-node-specific and global model. Note that the fitted curves in 
#' the tree do not convey a conditional function of the predictor on the $x$-axis
#' (as plotted when "terms" is specified). They are a function of the predictor 
#' on the $x$-axis, as well as all other predictors in the model and could thus be 
#' referred to as 'marginal' fitted curves.
#  @param which_terms character; \code{"local"}, \code{"global"} or \code{"both"}.
#  Only used when argument \code{which} equal \code{"global"} or \code{"both"}.
#  Specifies whether the local and/or global (smooth) terms should be plotted.
#' @param dim numeric vector of length two. Specifies how many rows and columns 
#' of plots should be fit on a single page.
#' \code{NULL} (the default) has the routine leave all settings as they are.
#' Using \code{par(mfrow = c( , ))} before plotting then provides control over 
#' the plot dimensions (and number of pages).
#' @param ylim \code{"firstplot"} (default), \code{NULL}, or a numeric vector of 
#' length 2. Only used for plotting the terminal-node models (not the tree). 
#' Specifies how the limits of the y-axes of the terminal node plots 
#' should be chosen. The default (\code{"firstnode"}) uses the observations 
#' in the first node to determine the limits of the y-axes for all plots. 
#' Alternatively, \code{NULL} will determine the limits of the y-axes 
#' separately for each plot. Alternatively, a numeric vector of length 
#' two may be specified, specifying the lower and upper limits of the 
#' y-axes.
#' @param treeplot_ctrl list of (named) arguments to be passed to 
#' \code{\link[partykit]{plot.party}}.
#' @param gamplot_ctrl list of (named) arguments to be passed to 
#' \code{\link[mgcv]{plot.gam}}. Note that not all arguments 
#' of \code{plot.gam} are supported. . 
#' @param ... further arguments, currently not used. 
#' 
#' @section Warning:
#' The plotted terms by default also represent confidence bands. These should
#' be taken with a big grain of salt, because they do NOT account for the 
#' searching of the tree structure; they assume the tree structure was
#' known in advance. They should be interpreted as overly optimistic and with
#' caution. 
#' @examples
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, 
#'                cluster = eco$Specimen) 
#' plot(gt, which = "tree") # default is which = 'both'
#' plot(gt, which = "terms")
#' @importFrom graphics plot par
#' @export
plot.gamtree <- function(x, which = "both", 
                         dim = NULL, ylim = "firstnode", treeplot_ctrl = list(), 
                         gamplot_ctrl = list(), ...) {
  
  ## Argument checking
  if (!which %in% c("tree", "terms", "both")) 
    warning("Argument which should specify one of 'tree', 'terms' or 'both'")
  
  ## Replace node info to allow for plotting
  tree_node <- as.list(x$tree$node)
  for (i in unique(x$tree$fitted[["(fitted)"]])) {
    ## To get the correct (or at least one) model in the terminal nodes
    tree_node[[i]]$info$object <- tree_node[[i]]$info$object$gam
  }
  x$tree$node <- as.partynode(tree_node)
  
  ## Plot tree if requested
  if (which != "terms") {
    treeplot_ctrl[["x"]] <- x$tree
    treeplot_ctrl[["terminal_panel"]] <- node_bivplot
    do.call(plot, treeplot_ctrl)
  }
  
  ## Plot GAM terms if requested
  if (which != "tree") {
    
    ## Check if there is a global part to the GAM
    #global_gam <- all(length(as.Formula(eval(x$call$formula))) == c(1, 3))
    #if (which_terms %in% c("global", "both") && !global_gam) {
    #  warning("Fitted GAM tree does not contain global terms. Plots of local smooths will be returned only.")
    #  which_terms <- "local"
    #}
    
    if (!is.null(dim)) par(mfrow = c(dim[1], dim[2]))
    
    # if (x$joint) {
    #   smooth_term_labels <- dimnames(summary(x$gamm)$s.table)[[1]]
    #   label_ids <- 1:length(smooth_term_labels)
    #   local_GAM_term_ids <- grep(".tree", smooth_term_labels)
    #   global_GAM_term_ids <- label_ids[!label_ids %in% local_GAM_term_ids]
    #   gamplot_ctrl[["x"]] <- x$gamm
    #   if (is.numeric(ylim) && length(ylim) == 2L) gamplot_ctrl[["ylim"]] <- ylim
    #   
    #   if (which_terms != "global") {
    #     ## Plot local smooths (and parametric terms)
    #     ## TODO: Allow for having an intercept-only model in the terminal nodes
    #     terminal_nodes <- levels(x$data$.tree)
    #     x_name <- as.character(x$tree[[1]]$node$info$object$pred.formula)[2]
    #     for (i in 1:length(local_GAM_term_ids)) {
    #       gamplot_ctrl[["select"]] <-  local_GAM_term_ids[i]
    #       gamplot_ctrl[["main"]] <- paste("node", terminal_nodes[i])
    #       if (i == 1L && !is.null(ylim)) {
    #         tmp <- do.call(plot, gamplot_ctrl)
    #         if (ylim[1L] == "firstnode") {
    #           c(min(tmp[[1L]]$fit - .5*tmp[[1L]]$se.mult*tmp[[1L]]$se), 
    #             max(tmp[[1L]]$fit + .5*tmp[[1L]]$se.mult*tmp[[1L]]$se))
    #         }
    #       } else {
    #         do.call(plot, gamplot_ctrl)        
    #       }
    #       ## add node-specific residuals to plot
    #       obs_ids <- predict(x$tree) == terminal_nodes[i]
    #       if (residuals) points(x = x$data[obs_ids, x_name], y = residuals(x$gamm)[obs_ids], pch = ".")
    #     }
    #   }
    #   
    #   if (which_terms != "local" && length(global_GAM_term_ids) > 0L) {
    #     ## Plot global smooths (and parametric terms)
    #     for (i in 1:length(global_GAM_term_ids)) {
    #       gamplot_ctrl[["select"]] <-  global_GAM_term_ids[i]
    #       gamplot_ctrl[["main"]] <- paste("global term", smooth_term_labels[global_GAM_term_ids[i]])
    #       do.call(plot, gamplot_ctrl)        
    #     }
    #   }
    #   
    # } else { ## joint is FALSE
    
    ## Plot global terms if requested
    #if (which_terms != "local") {
    #  gamplot_ctrl[["x"]] <- x$gamm
    #  do.call(plot, gamplot_ctrl)     
    #}
    
    ## Plot local terms if requested
    #if (which_terms != "global") {
    terminal_nodes <- c()
    for (i in 1:length(x$tree)) {
      if (is.null(x$tree[[i]]$node$kids)) terminal_nodes <- c(terminal_nodes, i)
    }
    for (i in terminal_nodes) {
      gamplot_ctrl[["x"]] <- x$tree[[i]]$node$info$object
      gamplot_ctrl[["main"]] <- paste("node", i)
      do.call(plot, gamplot_ctrl)   
    }
    #}
    #}
    
  }
}

#' @importFrom stats fitted
fitted.gamm4 <- function(object, ...) {
  fitted(object$gam, ...)
}

#' Extract estimated coefficients of a GAM tree
#' 
#' Returns the estimates global or local coefficients of a GAM tree.
#' 
#' @param object an object of class \code{gamtree}.
#' @param which character. Specifies whether local (\code{"local"}; 
#' node-specific estimates from the tree) or globally (\code{"global"})
#' estimated coefficients should be returned. 
#' @param ... currently not used.
#' 
#' @return Returns a matrix (if \code{which = "local"}) or a vector (if 
#' \code{which = "global"}) with estimated coefficients. 
#' ## GAM tree without global terms:
#' gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$Specimen)
#' coef(gt1)
#' 
#' ## GAM tree with global terms:
#' gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$Specimen)
#' coef(gt2)
#' coef(gt2, which = "global")
#' @export
#' @importFrom stats coef
coef.gamtree <- function(object, which = "local", ...) {

  if (which == "local") {
    if (object$joint) {
      ## Get coefs from GAM
      coefs <- coef(object$tree) # to obtain shape     
      coefs_gam <- coef(object$gamm)
      coefs_gam <- coefs_gam[grep(".tree", names(coefs_gam))]
      for (i in rownames(coefs)) {
        tmp <- coefs_gam[grep(paste0(".tree", i), names(coefs_gam))]
        intercept_id <- which(colnames(coefs) == "(Intercept)")
        coefs[i, intercept_id] <- tmp[paste0(".tree", i)]
        coefs[i, -intercept_id] <- tmp[grepl(paste0(":.tree"), names(tmp)) | 
                                         grepl(paste0(".tree", i, ":"), names(tmp)) ]
      }
    } else {
      ## Get coefs from tree
      coefs <- coef(object$tree)       
    }
  } else { # which == "global"
    ## Check if there is even a global components
    if (!all(length(as.Formula(object$call$formula)) == c(1, 3))) {
      warning("global coefficients were requested, but the specified GAM tree does not have a global model part")
      return(NULL)
    }
    if (object$joint) {
      ## Get coefs from GAM
      coefs <- coef(object$gamm)
      coefs <- coefs[-grep(paste0(".tree"), names(coefs))]
    } else {
      ## Get coefs from GAM
      coefs <- coef(object$gamm)
    }
  }
  return(coefs)
}


#' Predictions from fitted GAM tree
#'  
#' @description Takes a fitted GAM tree (of class \code{"gamtree"}) and produces 
#' predictions given a new set of values for the model covariates,
#' or for the original covariate values used for fitting the GAM tree.
#' 
#' @param object an object of class \code{gamtree}.
#' @param newdata a \code{data.frame} containing the values of the model
#' covariates for which predictions should be returned. The default 
#' (\code{NULL}) returns predictions for the original training data. 
#' @param type character vector of length 1, specifying the type of prediction
#' to be returned. \code{"link"} (the default) returns values on the scale
#' of the linear predictor. Alternatively, \code{"response"} returns predictions
#' on the scale of the response variable; \code{"node"} returns an integer vector
#' of node identifiers.
#' @param ... further arguments to be passed to \code{\link{mgcv}predict.gam}.
#' 
#' @return Returns a vector of predicted values.
#'  
#' @importFrom graphics par plot
#' @export
predict.gamtree <- function(object, newdata = NULL, type = "link", ...) {

  ## TODO: Implement possibility of excluding local and/or global effects
  ## TODO: Implement support for mgcv's gam.predict option type = "terms" or "iterms", 
  ## and standard errors can be requested with predictions through use of 
  ## argument "se.fit" 
  
  if (!inherits(object, "gamtree")) warning("predict.gamtree only works for objects of class gamtree")
  
  if (!type %in% c("link", "response", "node")) warning("Argument type should be one of 'link', 'response' or 'node'")
  
  ## TODO: foresee in 4 cases:
  ## joint=TRUE and global_gam=TRUE: reconstruct tree factor and predict with object$gamm
  ## joint=TRUE and global_gam=FALSE: reconstruct tree factor and predict with object$gamm
  ## joint=FALSE and global_gam=TRUE: get predictions from tree, add predictions from gamm
  ## joint=FALSE and global_gam=FALSE: get predictions from object$tree

  ## Check whether there is a global GAM
  global_gam <- all(length(as.Formula(object$call$formula)) == c(1, 3))
  
  ## Check whether newdata is supplied
  if (is.null(newdata)) {
    newdata_supplied <- FALSE
    newdata <- object$data
  } else {
    newdata_supplied <- TRUE    
  }
  
  ## Generate predictions
  if (type == "node") {
    preds <- factor(predict(object$tree, newdata = newdata, 
                   type = "node"), levels = levels(object$data$.tree), ...)
  } else if (object$joint) {
    if (!newdata_supplied) {
      newdata$.tree <- factor(predict(object$tree, newdata = newdata, 
                              type = "node"), levels = levels(object$data$.tree), ...)
    } 
    preds <- predict(object$gamm, newdata = newdata, type = type, ...)
  } else {
    if (global_gam) {
      preds <- predict(object$tree, newdata = newdata, type = "link", ...) + 
        predict(object$gamm, newdata = newdata, type = "link", ...)
      if (type == "response") {
        warning("For jointly estimated GAM trees with global terms, other prediction types than 'node' and 'link' are not yet supported")
      }
    } else {
      preds <- predict(object$tree, type = type, newdata = newdata, ...)
    }
  }

  ## Return predictions
  return(preds)

}




#' Sum of the observation-wise gradient contributions to the fitted model.
#' 
#' Computes the sum of the column-wise sums of the observation-wise gradient 
#' contributions to the fitted GAM in each of the tree nodes. This yields a 
#' sum for each coefficients in each node. This sum should be reasonably
#' close to zero. The more the column-wise sum of the observation-wise
#' gradient contribution differ
#' s from zero, the more the parameter 
#' stability tests w.r.t. to that coefficient will be overpowered.
#' 
#' @param object an object of class \code{gamtree}.
#' @param var logical. Should the variance of the observation-wise gradient
#' contributions be returned?
#' @param return_fits whether fitted GAMs should be returned.
#' 
#' @return Returns a matrix with the sum of the observation-wise gradient
#' contributions, with a row for each node of the tree and a column
#' for each coefficient. If \code{return_fits = TRUE}, a list with two
#' components is returned: The first component \code{$grad} contains the
#' a matrix with the column-wise sums of the gradient contributions. The
#' second component \code{$fits} contains a list of the node-specific GAM 
#' fits.
#' 
#' @importFrom sandwich estfun
#' @export
check_grad <- function(object, var = FALSE, return_fits = FALSE) {
  gamfits <- refit.modelparty(object$tree)
  if (inherits(gamfits, "list")) {
    ## Then at least one split was implemented and gamfits is a list of GAM fits
    est_fun <- estfun(gamfits[[1]]$mer)
    gradients <- colSums(est_fun)
    gradients <- as.data.frame(t(gradients))
    vars <- apply(est_fun, 2, var)
    vars <- as.data.frame(t(vars))
    for (i in 2L:length(object$tree)) {
      est_fun <- estfun(gamfits[[i]]$mer)
      gradients[i, ] <- colSums(est_fun)
      vars[i, ] <- apply(est_fun, 2, var)
    }
  } else {
    ## Then no splits were implemented and gamfits is a single GAM fit
    gradients <- colSums(estfun(gamfits))
    gradients <- as.data.frame(t(gradients))
    if (var) {
      gradients <- apply(estfun(gamfits), 2, var)
      gradients <- as.data.frame(t(vars))
    }
  }
  if (!return_fits && !var) {
    res <- gradients
  } else if (return_fits) {
    res <- list(grad = gradients, fits = gamfits)
    if (var) res$var <- vars
  } else if (var) {
    res <- list(grad = gradients, var = vars)
  }
  return(res)
}
