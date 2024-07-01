utils::globalVariables(c(".tree", ".offset", ".global", ".weights", ".cluster", "current_offset"))

## TODO: Get the node-specific plots with CIs from mgcv and put them together in one plot.
## TODO: Create cv.splinetree and cv.gamtree functions to obtain hypothesis tests.
## TODO: Allow splinetree to use fixed splines with mgcv.
## TODO: Build unit tests.
## TODO: Reduce computational load of merDeriv with 
## TODO: Speed up computations. Bottleneck of glmtree seems to be split selection. Check with AZ if 
## we can do other approach than exhaustive search. Summing scores within clusters in ytrafo may speed up
## ctree. Check with profiling whether merDeriv or ctree is the bottleneck, as ctree tests may involve lots
## of matrix operations as well. Finally, may be possible to improve speed of merDeriv.


## Impossible: Can enforce knot locations to apply to all nodes similarly? In mgcv::gam with smoothCon, not in gamm4.

## Done: Allow use of function s in the model formula for gamtree with method = "ctree" and "mob"
## Done: Can use a cubic spline basis in gamtree? Yes, just pass bs = "cr" to function s.
## Done: Use ctree for partitioning penalized GAMs.





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
#' @param REML logical, defaults to \code{TRUE}. Passed on to `gamm4` and in turn `lmer` (but 
#' not `glmer`) fitting routines to control whether REML or ML estimation is used.
#' @param method character, one of \code{"ctree"} or \code{"mob"}, indicates which
#' partitioning algorithm should be used. See details below. 
#' @param offset numeric vector of length \code{nrow(data)}. Supplies model 
#' offset for use in fitting. Note that this offset will always be completely 
#' ignored when predicting.
#' @param cluster optional, a name refering to a colum of \code{data}, or a 
#' numeric or factor vector with a cluster ID to be 
#' employed for clustered covariances in the parameter stability tests. 
#' Most useful if \code{method = "mob"}, for \code{method = "ctree"} probably
#' less so as it may yield overly conservative splitting. 
#' This argument should be used when the partitioning variables are not measured
#' on the individual observation level, but on a higher level. E.g., when 
#' the response variables consists of repeated measurements of the same
#' respondents. 
#' @param verbose logical. Should progress be printed to the commande line in 
#' every iteration? If true, the iteration number, information on the 
#' splitting procedure, and the log-likelihood (with df) value of the fitted 
#' full mixed-effects gam model is printed.
#' @param alt_formula list with two elements, for specifying non-standard model formulae
#' for GAM. E.g., the formula list required for use of the \code{\link[mgcv]{multinom}}
#' family.
#' @param gam_ctrl a list of fit control parameters to replace defaults returned by 
#' \code{\link[mgcv]{gam.control}}.
#' @param tree_ctrl a \code{list} of one or more control parameters as accepted by 
#' \code{\link[partykit]{mob_control}} (to be passed to function 
#' \code{\link[partykit]{mob}} if \code{method = "mob"}), or 
#' \code{\link[partykit]{ctree_control}} (to be passed to function
#' \code{\link[partykit]{ctree}} is \code{method = "ctree"}). 
#' Note: arguments \code{xtype} and \code{ytype} of \code{mob_control} are set to \code{"data.frame"}, 
#' by default, 
#' this cannot be changed. Argument \code{parm} of \code{mob_control} will be overruled
#' by the argument of the same name of the current function.
#' @param parm vector of one or more integers, indicating which parameters should be
#' included in the parameter stability tests. The default \code{c(1, 2, 4)} includes
#' the intercept, linear slope and error variance of the smoothing spline. The 3rd parameter is the variance
#' of smooth term. It is excluded by default, because its inclusion yields too high
#' power in many situations.
#' @param ... additional arguments to be passed to function \code{\link[gamm4]{gamm4}}. 
#' @return Returns an object of class \code{"gamtree"}. This is a list, containing
#' (amongst others) the GAM-based recursive partition (in \code{$tree}). 
#' The following methods are available to extract information from the fitted object:
#' \code{\link{predict.gamtree}}, for obtaining predicted values for training and new
#' observations; \code{\link{plot.gamtree}} for plotting the tree and variables' effects; 
#' \code{\link{coef.gamtree}}, \code{\link{fixef.gamtree}} and \code{\link{ranef.gamtree}}  
#' for extracting estimated coefficients. \code{\link{VarCorr.gamtree}} for extracting
#' random-effects (co)variances, \code{\link{summary.gamtree}} for a summary of the 
#' fitted models.
#' 
#' @examples
#' gt_m <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = Specimen)
#' summary(gt_m)
#' gt_c <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, method = "ctree")
#' summary(gt_c)
#' 
#' @details MOB is short for model-based recursive
#' partitioning, ctree is short for conditional inference tree. MOB is 
#' based more strongly on parametric theory, thereby allowing for easy inclusion 
#' of clustering structures into the estimation procedure (see also argument
#' \code{cluster}), yielding similar to a GEE-type approach for estimation of 
#' multilevel and longitudinal data structures. Yet, computation time for MOB is much
#' larger than for ctree, which is mostly due to how it searches for
#' the optimal splitting value, after the variable for splitting
#' has been selected. ctree uses tests based on permutation theory, 
#' and thereby offers a less parametrically oriented approach. It is much
#' faster than MOB, but does not provide a natural way of accounting
#' for multilevel or longitudinal data structures.
#' @seealso \code{\link{predict.gamtree}} \code{\link{plot.gamtree}} 
#' \code{\link{coef.gamtree}} \code{\link{summary.gamtree}}
#' @import mgcv partykit gamm4 merDeriv
#' @importFrom stats as.formula formula logLik predict update terms complete.cases
#' @importFrom lme4 fixef VarCorr
#' @importFrom Formula as.Formula Formula
#' @importFrom graphics par plot
#' @export
gamtree <- function(formula, data, weights = NULL, REML = TRUE, 
                    method = "mob",
                    cluster = NULL, 
                    offset = NULL, verbose = FALSE,
                    parm = c(1,2,4),
                    gam_ctrl = list(), 
                    tree_ctrl = list(), 
                    alt_formula = NULL, ...) {

  ## remember call
  cl <- match.call()
  
  if (!inherits(data, "data.frame")) warning("Argument data should specify a data.frame")
  
  tree_ctrl <- if (method == "mob") {
    tree_ctrl$parm <- parm
    tree_ctrl$xtype <- tree_ctrl$ytype <- "data.frame"
    do.call(partykit::mob_control, tree_ctrl)
  } else if (method == "ctree") {
    do.call(partykit::ctree_control, tree_ctrl)    
  }
  
  if (length(gam_ctrl) == 0L) gam_ctrl <- NULL
  ## Construct formulas for tree (tf), local gam (lgf) and global gam (ggf):
  ##  - tf is for passing response, predictors and part vars to mob()
  ##  - lgf is for fitting node-specific GAM in the nodes
  ##  - ggf is an intermediate formula, used for creating full GAM formula later
  ##  - fgf is for fitting full GAM
  ff <- Formula::as.Formula(formula)
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
    
  if (method == "mob") {

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
    
    ## initialization
    data$.offset <- if (is.null(offset)) 0L else offset
    
    ## grow tree
    if (is.null(q_cluster)) {
      tree <- mob(tf, data = data, local_gam_form = lgf, fit = gamfit, 
                  weights = .weights, offset = .offset, 
                  control = tree_ctrl, gam_ctrl = gam_ctrl, 
                  REML = REML, ...)   
    } else {
      tree <- mob(tf, data = data, local_gam_form = lgf, fit = gamfit, 
                  weights = .weights, offset = .offset, 
                  control = tree_ctrl, cluster = .cluster, gam_ctrl = gam_ctrl, 
                  REML = REML, ...)       
    }
    
  } else if (method == "ctree") {
    
    ## TODO: check if cluster argument is correctly passed, if not, implement

    ytrafo_gamm4 <- function(.lgf, .REML = TRUE, .response, .local_vars, 
                             .parm = c(1,2,4), ...) {
      
      gamm4_form <- gsub(as.character(.response), "y", gsub(.local_vars, "x", .lgf))
      gamm4_form <- formula(paste0(gamm4_form[2L], gamm4_form[1L], gamm4_form[3L])) 
      .REML <- eval(.REML)
      .parm <- eval(.parm)
      
      function(y, x, start = NULL, weights, offset, cluster = NULL, 
                            estfun = TRUE, object = TRUE, ...) {
        environment(gamm4_form) <- environment()
        mod <- gamm4::gamm4(gamm4_form, REML = .REML, ...)
        class(mod) <- "gamm4"
        list(object = mod, estfun = merDeriv::estfun.lmerMod(mod$mer, level = 1L)[, .parm])
      }
    }
    
    tree <- ctree(tf, data = data, 
                  ytrafo = ytrafo_gamm4(.lgf = lgf, .REML = REML, .response = response, 
                                        .local_vars = local_vars, .parm = parm), 
                  control = tree_ctrl, ...)
    
    ## Fill in terminal nodes if they have nobs < minsplit (ctree defaults: minsplit=20L, minbucket=7L)
    node_ids <- predict(tree, type = "node")
    tree_node <- as.list(tree$node)
    for (i in 1L:length(tree)) {
      if (is.null(tree[[i]]$node$info)) {
        
        ## TODO: Allow for passing further arguments to gamm4 
        node_data <- tree$data[node_ids == i, ]
        gamm4_args <- list(formula = lgf, REML = REML, #gam_ctrl = gam_ctrl, 
                           data = node_data, ...)
        mod <- do.call(gamm4::gamm4, gamm4_args) 
        class(mod) <- "gamm4"
        tree_node[[i]]$info <- list(
          criterion = NULL,
          p.value = NULL,
          object = mod,
          converged = TRUE,
          nobs = nrow(node_data))
      }
    }
    tree$node <- as.partynode(tree_node)
    tree$info$Formula <- Formula::Formula(tf)
  }

  ## collect results
  result <- list(
    tree = tree,
    formula = formula,
    data = data,
    call = cl,
    REML = REML,
    method = method
  )
  class(result) <- "gamtree"
  return(result)
}





#' Internal function for extracting fitted values from MOB-based GAM trees.
#' 
#' \code{fitted.gamm4} extract fitted values from objects of class \code{gamm4}.
#' 
#' @param object an object of class \code{gamm4}.
#' @param ... currently not used.
#' @importFrom stats fitted
#' @export
fitted.gamm4 <- function(object, ...) fitted(object$gam, ...)



#' Internal function for extracting predictions from MOB-based GAM trees.
#' 
#' \code{predict.gamm4} extract predictions from objects of class \code{gamm4}.
#' 
#' @param object an object of class \code{gamm4}.
#' @param newdata an optional \code{data.frame} in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... currently not used.
#' @export
predict.gamm4 <- function(object, newdata, ...) predict(object$gam, newdata = newdata, ...)



# #' @importFrom lme4 VarCorr
# #' @export
# VarCorr.gamm4 <- function(x, sigma=1, ...) VarCorr(x$mer, sigma=sigma, ...)

# #' @importFrom lme4 ranef
# #' @export
# ranef.gamm4 <- function(object, ...) ranef(object$mer, ...)

# #' @export
# summary.gamm4 <- function(object, ...) summary(object$gam, ...)





#' Extract random-effects covariance matrices from a GAM tree.
#' 
#' \code{VarCorr.gamtree} extracts fixed-effects random-effects covariance 
#' matrices from the nodes of a GAM tree. 
#' 
#' @param x an object of class \code{"gamtree"}.
#' @param sigma an optional numeric value used as a multiplier for the standard deviations.
#' @param which character. \code{"terminal"} (default) returns (co)variances for 
#' all terminal nodes, \code{"inner"} returns the (co)variances for all inner (splitting) 
#' nodes, \code{"all"} returns covariances for all nodes.
#' @param ... additional arguments to be passed to \code{\link[lme4]{VarCorr.merMod}}.
#' @importFrom lme4 VarCorr
#' @export
VarCorr.gamtree <- function(x, sigma=1, which = "terminal", ...) {
  
  nodes <- if (which == "all") {
    1L:length(x$tree)
    } else if (which == "terminal") {
      sort(unique(x$tree$fitted[["(fitted)"]])) 
    } else if (which == "inner") {
        (1L:length(x$tree))[-sort(unique(x$tree$fitted[["(fitted)"]]))]
    }

  vc <- list()
  counter <- 0L
  for (i in nodes) {
    counter <- counter + 1L
    vc[[counter]] <- as.data.frame(VarCorr(x$tree[[i]]$node$info$object$mer, ...))
    vc[[counter]] <- vc[[counter]][ , -which(colnames(vc[[counter]]) == "grp")]
    if (nrow(vc[[counter]]) == 2L) rownames(vc[[counter]]) <- c("Smooth", "Residual")
  }
  names(vc) <- paste("node", nodes)
  return(vc)
}




#' Extract fixed-effects coefficients from a GAM tree.
#' 
#' \code{fixef.gamtree} extracts fixed-effects coefficients from a GAM tree. 
#' 
#' @param object an object of class \code{"gamtree"}.
#' @param ... further arguments to be passed to \code{\link[lme4]{fixef.merMod}}.
#' @importFrom nlme fixef
#' @export
fixef.gamtree <- function(object, ...) {
  if (object$method == "mob") {
    fixeff <- do.call(rbind, coef(object$tree)[, "fixef"])
    colnames(fixeff) <- c("(Intercept)", 
                         all.vars(formula(Formula::as.Formula(object$formula), rhs = 1L, lhs = 0L)))
  } else if (object$method == "ctree") {
    rows <- sort(unique(object$tree$fitted[["(fitted)"]]))
    fixeff <- matrix(NA, nrow = length(rows), ncol = 2L, 
                     dimnames = list(as.character(rows), 
                                     c("(Intercept)", all.vars(formula(Formula::as.Formula(object$formula), rhs = 1L, lhs = 0L)))))
    for (i in rows) {
      fixeff[as.character(i), ] <- fixef(object$tree[[i]]$node$info$object$mer)  
    }
  }
  return(fixeff)
}
  



#' Extract coefficients from a GAM tree.
#' 
#' \code{coef.gamtree} extracts fixed- or random-effects coefficients from 
#' a GAM tree. 
#' 
#' @param object an object of class \code{"gamtree"}.
#' @param which character. Either \code{"fixed"} (default) or \code{"random"}, 
#' indicating that fixed- or random-effects coefficients should be returned,
#' respectively.
#' @param ... further arguments to be passed to \code{\link[lme4]{fixef.merMod}} or 
#' \code{\link[lme4]{ranef.merMod}}.
#' @importFrom stats coef
#' @export
coef.gamtree <- function(object, which = "fixed", ...) {
  if (which == "fixed") {
    fixef.gamtree(object)
  } else if (which == "random") {
    ranef.gamtree(object)
  }
}


#' Extract random-effects coefficients from a GAM tree.
#' 
#' \code{ranef.gamtree} extracts random-effects coefficients from a GAM tree. 
#' 
#' @param object an object of class \code{"gamtree"}.
#' @param ... further arguments to be passed to \code{\link[lme4]{ranef.merMod}}.
#' @importFrom nlme ranef
#' @export
ranef.gamtree <- function(object, ...) {
  rows <- sort(unique(object$tree$fitted[["(fitted)"]]))
  raneff <- ranef(object$tree[[1L]]$node$info$object$mer)$Xr
  x_name <- as.character(formula(object$tree$info$Formula, lhs = 0, rhs = 1L))[2L]
  raneff <- matrix(NA, nrow = length(rows), ncol = nrow(raneff), 
                   dimnames = list(as.character(rows), 
                                   paste0("s(", x_name, ").b", row.names(raneff))))
  for (i in rows) {
    raneff[as.character(i), ] <- t(ranef(object$tree[[i]]$node$info$object$mer)$Xr)
  }
  return(raneff)
}



##############################################################################
##
## gamm4-based fitting function, used only if method = "mob"
##
gamfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, 
                   estfun = NULL, object = NULL, ..., local_gam_form,
                   gam_ctrl = list(), REML = TRUE) {
  
  ## Prepare arguments to be passed to fitting function
  args <- list(...)
  args$REML <- REML
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







#' Print method for a fitted GAM tree
#' 
#' Prints the local and/or global terms in a fitted GAM tree.  
#'
#' @param x object of class \code{gamtree}.
#' @param ... further arguments to be passed to \code{\link[partykit]{print.modelparty}} or
#. \code{\link[partykit]{print.modelparty}}
#' @export
#' @examples
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = Specimen)
#' gt ## or: print(gt)
print.gamtree <- function(x, ...) {

  tree_node <- as.list(x$tree$node)  
  for (i in unique(predict(x$tree, type = "node")))
  if (x$method == "mob") {
    for (i in unique(x$tree$fitted[["(fitted)"]])) {
      ## To get the correct (or at least one) model in the terminal nodes
      tree_node[[i]]$info$coefficients <- print_gam(tree_node[[i]]$info$object$gam)
    }
    x$tree$node <- as.partynode(tree_node)
  } else if (x$method == "ctree") {
    ## TODO: Adjust ctree so it is printed in a more informative manner.
    ## Not sure if possible.
  }
  
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
      cat("\nEstimated degrees of freedom:\n")
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
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = Specimen)
#' summary(gt)
summary.gamtree <- function(object, ...) {

  #tree_node <- as.list(object$tree$node)
  #if (object$method == "mob") {
  #  for (i in unique(object$tree$fitted[["(fitted)"]])) {
  #    ## To get the correct (or at least one) model in the terminal nodes
  #    tree_node[[i]]$info$coefficients <- paste0("node ", i, print_gam(tree_node[[i]]$info$object$gam))
  #  }
  #  object <- as.partynode(tree_node)
  #} else {
  #  ## TODO: Adjust so that summary will be printed nicely
  #}
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
#'                cluster = Specimen) 
#' plot(gt, which = "tree") # default is which = 'both'
#' plot(gt, which = "terms")
#' @importFrom graphics plot par
#' @importFrom grid gpar
#' @export
plot.gamtree <- function(x, which = "both", 
                         ylim = "firstnode", treeplot_ctrl = list(), 
                         gamplot_ctrl = list(), ...) {
  
  ## TODO: Implement for method = "ctree"
  
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
    gamplot_data <- list()
    for (i in terminal_nodes) {
      gamplot_ctrl[["x"]] <- x$tree[[i]]$node$info$object
      gamplot_ctrl[["main"]] <- paste("node", i)
      gamplot_data[[i]] <- do.call(plot, gamplot_ctrl)   
    }
    #}
    #}
    invisible(gamplot_data)
    ## TODO: To do something with the returned result, see "Datasets and application.Rmd" rats example.
  }
}







#' Get predictions from fitted GAM tree
#'  
#' @description Takes a fitted GAM tree (of class \code{"gamtree"}) and returns 
#' predictions given a new set of values for the model covariates,
#' or for the original covariate values used for fitting the GAM tree.
#' 
#' @param object an object of class \code{gamtree}.
#' @param newdata a \code{data.frame} containing the values of the model
#' covariates for which predictions should be returned. The default 
#' (\code{NULL}) returns predictions for the original training data. 
#' @param type character vector of length 1, specifying the type of prediction
#' to be returned. \code{"response"} (the default) returns values on the scale
#' of the response variable. Alternatively, \code{"link"} (only available if 
#' \code{method = "mob"})returns values on the scale of the linear predictor;
#' \code{"node"} returns an integer vector of node identifiers.
#' @param ... further arguments to be passed to \code{\link[partykit]{predict.party}}.
#' 
#' @return Returns a vector of predicted values.
#'  
#' @export
predict.gamtree <- function(object, newdata = NULL, type = "link", ...) {

  ## TODO: Implement support for mgcv's gam.predict option type = "terms" or "iterms", 
  ## and standard errors can be requested with predictions through use of 
  ## argument "se.fit" 
  
  if (!inherits(object, "gamtree")) warning("predict.gamtree only works for objects of class gamtree")
  
  if (!type %in% c("link", "response", "node")) warning("Argument type should be one of 'link', 'response' or 'node'")

  ## Check whether there is a global GAM
  global_gam <- all(length(Formula::as.Formula(object$formula)) == c(1, 3))
  
  ## Check whether newdata is supplied
  if (is.null(newdata)) {
    newdata_supplied <- FALSE
    newdata <- object$data
  } else {
    newdata_supplied <- TRUE    
  }
  
  ## Generate predictions
  if (type == "node") {
    preds <- predict(object$tree, newdata = newdata, type = "node", ...)
  } else if (type %in% c("response", "link")) {
    if (object$method == "mob") {
      preds <- predict(object$tree, type = type, newdata = newdata, ...)
    } else if (object$method == "ctree") {
      
      ## Extract node memberships, initialize vector with predictions, recover 
      ## name of node-specific predictor variable and compute predictions node-by-node
      node_ids <- predict(object$tree, newdata = newdata, type = "node")
      preds <- rep(NA, times = length(node_ids))
      x_name <- try(all.vars(formula(Formula::as.Formula(object$formula), lhs = 0, rhs = 1L)))
      if (inherits(x_name, "try-error")) stop("Name of predictor variable could not be recovered.")
      if (length(x_name) != 1L) warning("There seem to be multiple predictor variables for the node-specific model, results may not be correct.")
      newdata$x <- newdata[ , x_name]
      for (i in unique(node_ids)) {
        preds[node_ids == i] <- predict(
          object$tree[[i]]$node$info$object$gam, newdata = newdata[node_ids == i, ],
          type = type, ...)
      }
    }
  }

  ## Return predictions
  return(preds)

}
