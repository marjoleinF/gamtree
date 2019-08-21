utils::globalVariables(c(".tree", ".offset", ".global", ".weights", ".cluster"))

## From mob() documentation:
##
## All tests also require an estimate of the corresponding variance-covariance
## matrix of the estimating functions. The default is to compute this using an  
## outer-product-of-gradients (OPG) estimator. Alternatively, the corrsponding  
## information matrix or sandwich matrix can be used if: (a) the estimating  
## functions are actually maximum likelihood scores, and (b) a vcov() method
## (based on an estimate of the information) is provided for the fitted model 
## objects. The corresponding option in mob_control() is to set vcov ="info" 
## or vcov="sandwich"rather than vcov="opg" (the default).
##
## mgcv supplies a vcov.gam method, so condition b) is met, but unclear if a) 
## is met.



#' Recursively partition a dataset based on a (local) GAM, while accounting for 
#' globally specified GAM terms.
#'
#' \code{gamtree} recursively partition a dataset based on a (local) GAM, in
#' which a smooth function (thin plate regression spline) of a single predictor
#' variable is estimated, while accounting for globally specified smooth 
#' functions and random effects.
#' 
#' @param formula specifies the model formula, consisting of three or four 
#' parts: the response variable followed by a tilde, the terms for the 
#' node-specific GAM followed by a vertical bar, the partitioning variables 
#' followed by a vertical bar, and optionally the terms for the global GAM. 
#' See examples below.
#' @param data specifies the dataset (must be a data.frame).
#' @param weights currently ignored! optional numeric vector of case weights.
#' @param cluster optional numeric or factor vector with a cluster ID to be 
#' passed employed for clustered covariances in the parameter stability tests.
#' @param offset currently ignored! numeric vector with an a priori known 
#' component to be included
#' in the model \code{y ~ x1 + x2 + ....}. Will be applied to both the local
#' and global model and will not be updated in consecutive iterations. 
#' @param abstol specifies the convergence criterion. If the log-likelihood 
#' values of the full mixed-effects gam, from two consecutive iterations differ 
#' less than abstol, estimation is halted.
#' @param verbose logical. Should progress be printed to the commande line in 
#' every iteration? If true, the iteration number, information on the 
#' splitting procedure, and the log-likelihood (with df) value of the fitted 
#' full mixed-effects gam model is printed.
#' @param method The smoothing parameter estimation method used by \code{gam()}.
#' See documentation of function \code{gam()} for details.
#' @param mob_ctrl a list with control parameters as returned by \code{mob_control}
#' to be passed to function \code{mob()}. (Note that argument `xtype` is set to 
#' `data.frame`, by default, and cannot be adjusted.)
#' @param ... additional arguments to be passed the locally fitted \code{gam}s
#' , through \code{mob()}'s fitting function (currently, the only option is the 
#' default \code{gam_fit()}).
#' 
#' @return Returns an object of class \code{"gamtree"}. This is a list, containing
#' (amongst others) the GAM-based recursive partition (in \code{$tree}), and the 
#' fitted full GAM with both local and/or global fitted effects (in \code{$gamm}). 
#' 
#' @examples
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$specimen, verbose = FALSE)
#' summary(gt)
#' 
#' @import mgcv partykit Formula
#' @importFrom stats as.formula formula logLik predict update terms
#' @export
# formula <- score ~ s(months) | s(CHILDID, bs = "re") | GENDER + RACE + WKSESL + 
#   C1GMOTOR + C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + 
#   AGEBASELINE 
# formula <- score ~ s(months) |  GENDER + RACE + WKSESL + 
#   C1GMOTOR + C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + 
#   AGEBASELINE 
# data <- traindata
# weights <- cluster <- offset <- mob_ctrl <- NULL
# abstol <- 0.001
# verbose <- TRUE
# method <- "REML"
gamtree <- function(formula, data, weights = NULL, cluster = NULL, 
                    offset = NULL, abstol = 0.001, 
                    verbose = TRUE, method = "REML", mob_ctrl = NULL, 
                    ...) {
  
  data$.weights <- if (is.null(weights)) 1 else weights
  if (!is.null(cluster)) data$.cluster <- cluster
  
  if (is.null(mob_ctrl)) {
    mob_ctrl <- mob_control(verbose = verbose, xtype = "data.frame",
                            ytype = "data.frame")
  } else {
    mob_ctrl$verbose <- verbose
    mob_ctrl$ytype <- mob_ctrl$xtype <- "data.frame"
  }
  
  ## Construct formulas for tree (tf), local gam (lgf) and global gam (ggf):
  ff <- as.Formula(formula)
  lgf <- formula(ff, lhs = NULL, rhs = 1)
  global_gam <- all(length(Formula(ff)) == c(1,3))
  if (global_gam) {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 3))
    ggf <- formula(ff, lhs = NULL, rhs = 2)
  } else {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 2))
    ggf <- formula(ff, lhs = NULL, rhs = 0)
  }
  response <- ff[[2]] 
  if (length(response) > 1L && is.call(response[1])) {
    response <- as.character(response)
    response <- paste0(response[1], "(", response[2], ", ", response[3], ")")
  } 
  tf <- formula(paste(response, "~", 
                paste(local_vars, collapse = " + "), "|", 
                paste(part_vars, collapse = " + ")))

  ## Construct formulas for full gam (fgf):
  new_alt <- new <- attr(terms(lgf), "term.labels")
  new <- gsub(")", ", by = .tree)", new)
  new[!grepl(")", new)] <- paste0(new[!grepl(")", new)], ":.tree") 
  new <- paste(new, collapse = " + ")
  new_alt <- paste(new_alt, collapse = " + ")
  new <- paste0("~ 0 + .tree + ", new, " + .")
  new_alt <- paste0("~ ", new_alt, "+ .")
  fgf <- update(old = ggf, new = formula(new))
  ## TODO: This should not omit intercept, but does in fgf_alt for 
  ## formula <- cbind(succ, fail) ~ s(PAR) | Species
  ## Only for cbind response? Only for local-only gamtree?
  fgf_alt <- update(old = ggf, new = formula(new_alt))
  
  ## remember call
  cl <- match.call()
  
  ## initialization
  data$.offset <- if (is.null(offset)) 0 else offset
  data$.global <- 0
  newloglik <- -Inf
  oldloglik <- c(-Inf, -Inf) # last element is the oldest
  iteration <- 0L
  continue <- TRUE
  
  ## iterate between tree and global model estimation
  while (continue) {
    
    iteration <- iteration + 1L
    if (verbose) print(paste("iteration", iteration))
    
    if (!is.null(offset)) data$.global <- data$.global + offset
    
    ## Grow tree and get node memberships:
    if (is.null(cluster)) {
      tree <- mob(tf, data = data, local_gam_form = lgf, 
                  fit = gam_fit, method = method, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, ...)   
    } else {
      tree <- mob(tf, data = data, local_gam_form = lgf, 
                  fit = gam_fit, method = method, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, cluster = .cluster, ...)       
    }
    data$.tree <- factor(predict(tree, newdata = data, type = "node"))
    
    ## Fit gam with global and local models:
    if (length(tree) > 1L) {
      gamm <- gam(fgf, data = data, method = method, weights = .weights,
                  offset = .offset, ...)
    } else {
      gamm <- gam(fgf_alt, data = data, method = method, weights = .weights,
                  offset = .offset, ...)
    }
    
    ## Obtain predictions of global model only:
    if (iteration == 1L) {
      ## TODO: account for possibility that the 1st tree has size 1 and there
      ##       is no .tree variable
      keep_terms <- c()
      ## Get global parametric terms:      
      if (sum(gamm$nsdf) > 0L) {
        if (length(gamm$nsdf) > 1L) {
          pstart <- attr(gamm$nsdf, "pstart")
          ind <- rep(0, 0)
          for (i in 1L:length(gamm$nsdf)) if (gamm$nsdf[i] > 0L) {
            ind <- c(ind, pstart[i]:(pstart[i] + gamm$nsdf[i] - 1L))
          }
        } else {
          pstart <- 1L
          ind <- 1L:gamm$nsdf
        }
        par_terms <- names(gamm$coefficients[ind])
        keep_terms <- par_terms[!grepl(".tree", par_terms)]
      }
      ## Get global smooth terms:
      for (i in 1:length(gamm$smooth)) {
        if (!grepl(":.tree", gamm$smooth[[i]]$label)) {
          keep_terms <-  c(keep_terms, gamm$smooth[[i]]$label)    
        }
      }
    }
    if (length(inds <- which(keep_terms == "(Intercept)")) > 0L) {
      keep_terms <- keep_terms[-inds]
    }
    data$.global <- rowSums(predict(gamm, newdata = data, type = "terms",
                                 terms = keep_terms))
    
    ## iteration information
    newloglik <- logLik(gamm)    
    continue <- abs(newloglik - oldloglik[1]) > abstol
    if (continue & (abs(newloglik - oldloglik[2]) < abstol)) {
      if (newloglik > oldloglik[1]) continue <- FALSE
    }
    oldloglik[2] <- oldloglik[1]
    oldloglik[1] <- newloglik
    if (verbose) print(newloglik)
    
  }
  
  ## collect results
  result <- list(
    formula = formula,
    gamm_form = fgf,
    call = cl,
    tree = tree,
    gamm = gamm,
    data = data,
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    iterations = iteration, 
    abstol = abstol
  )
  class(result) <- "gamtree"
  return(result)
}



#############################################################
##
## Fitting function for gam-based recursive partitioning with MOB.
##
gam_fit <- function(y, x, start = NULL, weights = NULL, 
                    offset = NULL, local_gam_form, ...) {
  
  ## TODO: Allow for taking weights and offset arguments
  gam(local_gam_form, data = cbind(x, y), ...)
  
}



#' Summary for a fitted GAM tree
#' 
#' Prints a summary of the local and/or global terms in a fitted GAM tree.  
#'
#' @param object object of class \code{gamtree}.
#' @param ... further arguments to be passed to \code{summary.gam}.
#' @export
summary.gamtree <- function(object, ...) {
  summary(object$gamm, ...)
}



#' Plotting of GAM trees
#' 
#' Takes a fitted GAM tree and plots the smooth functions fitted in each of the
#' terminal nodes of the tree.
#' 
#' @param x object of class \code{gamtree}.
#' @param which character. The default (\code{"both"}) plots the tree structure, 
#' followed by the model fitted in the terminal nodes. Alternatively, \code{"tree"}
#' will plot the tree structure only, and \code{"nodes"} will plot the terminal-node
#' models only.
#' @param treeplot_ctrl list of arguments to be passed to \code{plot.party()}.
#' @param gamplot_ctrl list of arguments to be passed to \code{plot,gam()}. 
#' @param ... further arguments, currently not used. 
#' 
#' @examples 
#' gt <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, verbose = FALSE, cluster = eco$specimen) 
#' plot(gt, which = "tree") # fdefault is which = 'both'
#' plot(gt, which = "nodes", residuals = TRUE)
#' 
#' @importFrom graphics plot par
#' @export
plot.gamtree <- function(x, which = "both", treeplot_ctrl = list(), 
                         gamplot_ctrl = list(), ...) {
  ## TODO: allow for plotting local (current default) and global smooths 
  if (which != "nodes") {
    ## Plot observed data in terminal nodes:
    treeplot_ctrl[["x"]] <- x$tree
    treeplot_ctrl[["terminal_panel"]] <- node_bivplot
    do.call(plot, treeplot_ctrl)
  }
  if (which != "tree") {
    ## Plot models in terminal nodes:
    terminal_nodes <- c()  
    for (i in 1:length(x$tree)) {
      if (is.null(x$tree[[i]]$node$kids)) terminal_nodes <- c(terminal_nodes, i)
    }
    nrow <- ncol <- ceiling(sqrt(length(terminal_nodes)))
    if (length(terminal_nodes) <= (nrow-1L)*ncol) nrow <- nrow - 1L
    par(mfrow = c(nrow, ncol))
    for (i in terminal_nodes) {
      ## TODO: base these plots on fitted gam model. 
      ## This requires specifying the 'select' argument, which appears to
      ## count smooth terms first, then parametric terms, but only of order 1,
      ## that is, only intercept terms are plotted.
      gamplot_ctrl[["x"]] <- x$tree[[i]]$node$info$object
      gamplot_ctrl[["main"]] <- paste("node", i)
      do.call(plot, gamplot_ctrl)
    }
  }
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
#' 
#' @export
#' @importFrom stats coef
coef.gamtree <- function(object, which = "local", ...) {
  coefs <- object$gamm$coefficients
  if (all(!grepl(".tree", names(coefs)))) {
    warning("No splits were made in the tree. All estimated coefficients are global and all will be returned.")
  } else {
    local_coef_ids <- grep(".tree", names(coefs))
    if (which == "local") {
      local_coefs <- coefs[local_coef_ids]
      no_nodes <- width(object$tree)
      coefs <- coef(object$tree)
      coefs[, 1] <- local_coefs[1:no_nodes]
      coefs[, -1] <- matrix(local_coefs[-(1:no_nodes)], byrow = TRUE, nrow = no_nodes)
    } else if (which == "global") {
      coefs <- coefs[-local_coef_ids]
    }
  }
  return(coefs)
}




#' Predictions from fitted GAM tree
#'  
#' Takes a fitted GAM tree (produced by \code{gamtree()} and produces 
#' predictions given a new set of values for the model covariates,
#' of for the original covariate values used for fitting the GAM tree.
#' 
#' @param object an object of class \code{gamtree}.
#' @param newdata a \code{data.frame} containing the values of the model
#' covariates for which predictions should be returned. The default 
#' (\code{NULL}) returns predictions for the original training data. 
#' @param ... further arguments to be passed to \code{predict.gam}.
#' 
#' @return Returns a vector of predicted values.
#'  
#' @importFrom graphics par plot
#' @export
predict.gamtree <- function(object, newdata = NULL, ...) {
  
  if (is.null(newdata)) {
    newdata <- object$data
  } else {
    newdata$.tree <- factor(predict(object$tree, newdata = newdata, 
                                    type = "node"), 
                            levels = levels(object$data$.tree))
  }
  ## TODO: Implement possibility of excluding local and/or global effects
  predict(object$gamm, newdata = newdata, ...)
}




#' Sum of the observation-wise gradient contributions to the fitted model.
#' 
#' Computes the sum of the observation-wise gradient contributions to the
#' fitted model in each of the terminal nodes. The sums should be reasonably
#' close to zero.
#' 
#' @param object an object of class \code{gamtree}.
#' 
#' @return Returns a matrix with the sum of the observation-wise gradient
#' contributions, with a row for each terminal node of the tree and a column
#' for each coefficient.
#' 
#' @importFrom sandwich estfun
#' @export
check_grad <- function(object) {
  gamfits <- refit.modelparty(object$tree)
  gradients <- colSums(estfun(gamfits[[1]]))
  gradients <- as.data.frame(t(gradients))
  for (i in 1L:length(object$tree)) {
    gradients[i, ] <- colSums(estfun(gamfits[[i]]))
  }
  return(gradients)
}
