utils::globalVariables(c(".tree", ".offset", ".global", ".weights", ".cluster"))

#' Recursively partition a dataset based on a (local) GAM, while accounting for 
#' global GAM terms.
#'
#' \code{gamtree} recursively partition a dataset based on a (local) GAM, in
#' which a smooth function (spline) of a predictor variable is estimated, 
#' while accounting for globally specified smooth functions and random effects.
#' 
#' @param formula specifies the model formula, consisting of three or four 
#' parts: the response variable followed by a tilde, the terms for the 
#' node-specific GAM followed by a vertical bar, optionally the terms for 
#' the global GAM followed by a vertical bar, and then the partitioning 
#' variables. See examples below.
#' @param data specifies the dataset (must be a data.frame).
#' @param weights currently ignored! optional numeric vector of case weights.
#' @param cluster optional numeric or factor vector with a cluster ID to be 
#' employed for clustered covariances in the parameter stability tests.
#' This argument should be used when the partitioning variables are not measured
#' on the individual observation level, but on a higher level. E.g., when 
#' the response variables consists of repeated measurements on multiple 
#' observations and the partitioning variables are time-invariant covariates.
#' @param joint Should local models be re-estimated jointly along with the
#' global model? Jointly re-estimating may yield a more accurate final model,
#' and may need less iterations to converge, but likely yield longer computation
#' times.
#' @param globalstart Should estimation initialize with the global model? Defaults 
#' \code{FALSE}, resulting in the tree (partition or subgroup structure) being 
#' estimated, first. If set to \code{TRUE}, (ignored if no global part has been
#' specified in the \code{formula}), the global part of the model will be estimated
#' first.  
#' @param offset currently ignored!
#' @param abstol numeric. Specifies the convergence criterion. If the log-likelihood 
#' values of the full mixed-effects gam, from two consecutive iterations differ 
#' less than abstol, estimation is halted.
#' @param verbose logical. Should progress be printed to the commande line in 
#' every iteration? If true, the iteration number, information on the 
#' splitting procedure, and the log-likelihood (with df) value of the fitted 
#' full mixed-effects gam model is printed.
#' @param method character. Smoothing parameter estimation method used by \code{gam()}.
#' See \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}} for details.
#' @param bigdata logical. Should function \code{\link[mgcv]{bam}} be used to fit
#' the GAMs? Defaults to \code{FALSE}, resulting in function \code{\link[mgcv]{gam}} 
#' being used for fitting the GAMs. \code{TRUE} will result in function 
#' \code{\link[mgcv]{bam}} being used for fitting the GAMs.
#' @param mob_ctrl a list with control parameters as returned by 
#' \code{\link[partykit]{mob_control}} to be passed to function 
#' \code{\link[partykit]{mob}}. Note: argument `xtype` is set to 
#' `data.frame`, by default, and cannot be adjusted.
#' @param ... additional arguments to be passed to the locally fitted 
#' \code{gam}s, through \code{mob()}'s fitting function (currently, the only 
#' option is the default \code{gam_fit()}).
#' 
#' @return Returns an object of class \code{"gamtree"}. This is a list, containing
#' (amongst others) the GAM-based recursive partition (in \code{$tree}), and the 
#' fitted full GAM with both local and/or global fitted effects (in \code{$gamm},
#' if terms for the global GAM were specified in \code{formula}). 
#' 
#' @examples
#' ## GAM tree without global terms:
#' gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$specimen)
#' summary(gt1)
#' 
#' ## GAM tree with global terms:
#' gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$specimen)
#' summary(gt2)
#' 
#' @seealso \code{\link{predict.gamtree}} \code{\link{plot.gamtree}} 
#' \code{\link{coef.gamtree}} \code{\link{summary.gamtree}}
#' @import mgcv partykit Formula
#' @importFrom stats as.formula formula logLik predict update terms complete.cases
#' @export
gamtree <- function(formula, data, weights = NULL, cluster = NULL, 
                    joint = TRUE, globalstart = FALSE, offset = NULL, 
                    abstol = 0.001, verbose = FALSE, method = "REML", 
                    bigdata = FALSE, mob_ctrl = mob_control(), ...) {

  ## Construct formulas.
  ## For tree (tf), local gam (lgf) and global gam (ggf).
  ##  - tf is for passing response, predictoors and part vars to mob()
  ##  - lgf is for fitting node-specific GAM in the nodes
  ##  - fgf is for fitting full GAM
  ff <- as.Formula(formula)
  lgf <- formula(ff, lhs = NULL, rhs = 1)
  global_gam <- all(length(Formula(ff)) == c(1,3))
  if (global_gam) {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 3))
    ggf <- formula(ff, lhs = NULL, rhs = 2)
    global_vars <- all.vars(formula(as.Formula(as.Formula(ggf), lhs = 0)))
  } else {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 2))
    #ggf <- formula(ff, lhs = NULL, rhs = 0)
  }
  if (length(local_vars) == 0) local_vars <- 1
  response <- ff[[2]] 
  if (length(response) > 1L && is.call(response[1])) {
    response <- as.character(response)
    response <- paste0(response[1], "(", response[2], ", ", response[3], ")")
  }
  tf <- formula(paste(response, "~", 
                paste(local_vars, collapse = " + "), "|", 
                paste(part_vars, collapse = " + ")))
  
  ## Construct formula for full gam (fgf).
  ##  - fgf_alt is for fitting the full GAM if mod() did not split
  ##  - If joint = FALSE, then fgf = fgf_alt = ggf + .offset term
  if (global_gam) {
    if (joint) {
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
    } else {
      fgf <- fgf_alt <- update(ggf, new = formula(" ~ . + offset(.offset)"))
    }
  }
  
  ## Prepare data
  if (global_gam) {
    data <- data[ , c(as.character(response), local_vars, part_vars, global_vars)]
  } else {
    data <- data[ , c(as.character(response), local_vars, part_vars)]    
  }
  N <- nrow(data)
  data <- data[complete.cases(data), ]
  if (nrow(data) != N) {
    warning(paste(N - nrow(data), "observations were removed due to missing values."))
    N <- nrow(data)
  }
  if (!is.null(offset)) {
    warning("Argument offset was specified but will be ignored.")
    offset <- NULL
  }
  data$.weights <- if (is.null(weights)) 1 else weights
  if (!is.null(cluster)) data$.cluster <- cluster
  
  if (is.null(mob_ctrl)) {
    mob_ctrl <- mob_control(verbose = verbose, xtype = "data.frame",
                            ytype = "data.frame")
  } else {
    mob_ctrl$verbose <- verbose
    mob_ctrl$ytype <- mob_ctrl$xtype <- "data.frame"
  }
  
  ## remember call
  cl <- match.call()
  
  ## initialization
  data$.offset <- if (is.null(offset)) 0 else offset
  ## TODO: Allow for starting with estimation of global model
  if (global_gam && globalstart) {
    if (bigdata) {
      gamm <- bam(fgf_alt, data = data, method = method, weights = .weights, 
                  ...)
    } else {
      gamm <- gam(fgf_alt, data = data, method = method, weights = .weights, 
                  ...)
    }
    data$.global <- predict(gamm, newdata = data)
  } else {
    data$.global <- 0    
  }
  newloglik <- -Inf
  oldloglik <- c(-Inf, -Inf) # 2nd element is the oldest
  iteration <- 0L
  continue <- TRUE
  
  ## iterate between tree and global model estimation
  while (continue) {
    
    iteration <- iteration + 1L
    if (verbose) print(paste("iteration", iteration))
    if (!is.null(offset)) data$.global <- data$.global + offset
    if (bigdata) fit <- bam_fit else fit <- gam_fit 
    ## grow tree and get node memberships
    if (is.null(cluster)) {
      tree <- mob(tf, data = data, local_gam_form = lgf, 
                  fit = fit, method = method, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, ...)   
    } else {
      tree <- mob(tf, data = data, local_gam_form = lgf, 
                  fit = fit, method = method, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, cluster = .cluster, ...)       
    }
    if (global_gam) {
      if (joint) {
        data$.tree <- factor(predict(tree, newdata = data, type = "node"))
      } else {
        data$.offset <- if (is.null(offset)) {
          predict(tree, newdata = data, type = "link")
        } else {
          offset + predict(tree, newdata = data, type = "link")
        }
      }
    
      ## Fit gam with global and local models:
      if (length(tree) > 1L) {
        if (bigdata) {
          gamm <- bam(fgf, data = data, method = method, weights = .weights, 
                      ...)          
        } else {
          gamm <- gam(fgf, data = data, method = method, weights = .weights, 
                      ...)
        }
      } else {
        if (bigdata) {
          gamm <- bam(fgf_alt, data = data, method = method, weights = .weights, 
                      ...)
        } else {
          gamm <- gam(fgf_alt, data = data, method = method, weights = .weights, 
                      ...)
        }
      }
    
      ## Obtain predictions of global parts of full model
      if (iteration == 1L) {
        keep_terms <- colnames(predict(gamm, newdata = data[1,], type = "terms"))
        keep_terms <- keep_terms[!grepl(".tree", keep_terms)]
      }
      if (length(tree) == 1L) {
        ## keep only global terms:
        keep_terms1 <- keep_terms
        for (i in local_vars) {
          keep_terms1 <- keep_terms[-grep(i, keep_terms)]
        }
        data$.global <- rowSums(predict(gamm, newdata = data, type = "terms",
                                        terms = keep_terms1))
      } else {
        data$.global <- rowSums(predict(gamm, newdata = data, type = "terms",
                                        terms = keep_terms))
      }
      
      ## iteration information
      newloglik <- logLik(gamm)    
      continue <- abs(newloglik - oldloglik[1]) > abstol
      if (continue & (abs(newloglik - oldloglik[2]) < abstol)) {
        if (newloglik > oldloglik[1]) continue <- FALSE
      }
      oldloglik[2] <- oldloglik[1]
      oldloglik[1] <- newloglik
      if (verbose) print(newloglik)
    } else continue <- FALSE 
  }
  
  ## collect results
  result <- list(
    formula = formula,
    gamm_form = if (global_gam) fgf else NULL,
    call = cl,
    tree = tree,
    gamm = if (global_gam) gamm else NULL,
    data = data,
    loglik = if (global_gam) as.numeric(newloglik) else NULL,
    df = if (global_gam) attr(newloglik, "df") else NULL,
    iterations = if (global_gam) iteration, NULL, 
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
  ## See partykit:::glmfit
  if (is.null(x)) x <- rep(1, length = length(y)) 
  gam(local_gam_form, data = cbind(x, y), ...)
  
  ## From mob() documentation:
  ## All tests also require an estimate of the corresponding variance-covariance
  ## matrix of the estimating functions. The default is to compute this using an  
  ## outer-product-of-gradients (OPG) estimator. Alternatively, the corrsponding  
  ## information matrix or sandwich matrix can be used if: (a) the estimating  
  ## functions are actually maximum likelihood scores, and (b) a vcov() method
  ## (based on an estimate of the information) is provided for the fitted model 
  ## objects. The corresponding option in mob_control() is to set vcov ="info" 
  ## or vcov="sandwich"rather than vcov="opg" (the default).
  ## mgcv supplies a vcov.gam method, so condition b) is met, a) also met?
}


#############################################################
##
## Fitting function for bam-based recursive partitioning with MOB.
##
bam_fit <- function(y, x, start = NULL, weights = NULL, 
                    offset = NULL, local_gam_form, ...) {
  
  ## TODO: Allow for taking weights and offset arguments
  ## See partykit:::glmfit
  if (is.null(x)) x <- rep(1, length = length(y)) 
  bam(local_gam_form, data = cbind(x, y), ...)

}




#' Summary method for a fitted GAM tree
#' 
#' Prints a summary of the local and/or global terms in a fitted GAM tree.  
#'
#' @param object object of class \code{gamtree}.
#' @param ... further arguments to be passed to \code{summary.gam}.
#' @export
#' @examples
#' ## GAM tree without global terms:
#' gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$specimen)
#' summary(gt1)
#' 
#' ## GAM tree with global terms:
#' gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$specimen)
#' summary(gt2)
summary.gamtree <- function(object, ...) {
  if (is.null(object$gamm)) {
    summary(object$tree, ...)
  } else {
    summary(object$gamm, ...)
  }
}



#' Plotting method for GAM trees
#' 
#' Takes a fitted GAM tree and plots the smooth functions fitted in each of the
#' terminal nodes of the tree.
#' 
#' @param x object of class \code{gamtree}.
#' @param which character. The default (\code{"both"}) plots the tree structure, 
#' followed by the model fitted in the terminal nodes. Alternatively, \code{"tree"}
#' will plot the tree structure only, and \code{"nodes"} will plot the smooths from
#' the terminal-node-specific and global model.
#' @param which_terms character; \code{"local"}, \code{"global"} or \code{"both"}.
#' Specifies whether the local and/or global smooth terms should be plotted.
#' @param ylim \code{"firstplot"} (default), \code{NULL}, or a numeric vector of 
#' length 2. Only used for plotting the terminal-node models (not the tree). 
#' Specifies how the limits of the y-axes of the terminal node plots 
#' should be chosen. The default (\code{"firstnode"}) uses the observations 
#' in the first node to determine the limits of the y-axes for all plots. 
#' Alternatively, \\code{NULL} will determine the limits of the y-axes 
#' separately for each plot. Alternatively, a numeric vector of length 
#' two may be specified, specifying the lower and upper limits of the 
#' y-axes.
#' @param treeplot_ctrl list of (named) arguments to be passed to 
#' \code{plot.party()}.
#' @param gamplot_ctrl list of (named) arguments to be passed to 
#' \code{plot.gam()}. 
#' @param ... further arguments, currently not used. 
#' 
#' @examples
#' gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, 
#'                cluster = eco$specimen) 
#' plot(gt1, which = "tree") # default is which = 'both'
#' plot(gt1, which = "nodes", gamplot_ctrl = list(residuals = TRUE)) 
#' gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$specimen) 
#' plot(gt2, which = "tree") # default is which = 'both'
#' plot(gt2, which = "nodes", gamplot_ctrl = list(residuals = TRUE))
#' 
#' 
#' @importFrom graphics plot par
#' @export
plot.gamtree <- function(x, which = "both", which_terms = "both", 
                         ylim = "firstnode", treeplot_ctrl = list(), 
                         gamplot_ctrl = list(), ...) {

  if (which != "nodes") {
    ## Plot observed data in terminal nodes:
    treeplot_ctrl[["x"]] <- x$tree
    treeplot_ctrl[["terminal_panel"]] <- node_bivplot
    do.call(plot, treeplot_ctrl)
  }
  if (which != "tree") {
    
    ##todo: if there is no global gam:
    if (is.null(x$gamm)) {
      terminal_nodes <- unique(predict(x$tree))
      for (i in terminal_nodes) {
        plot(x$tree[[i]]$node$info$object)
      }
    } else {
      joint <- is.null(x$call$joint)
      if (!joint) joint <- x$call$joint
      if (joint) {
        smooth_term_labels <- dimnames(summary(x$gamm)$s.table)[[1]]
        label_ids <- 1:length(smooth_term_labels)
        local_GAM_term_ids <- grep(".tree", smooth_term_labels)
        global_GAM_term_ids <- label_ids[!label_ids %in% local_GAM_term_ids]
        if (which_terms == "local") {
          number_of_plots <- length(terminal_nodes)
        } else if (which_terms == "both") {
          number_of_plots <- length(smooth_term_labels)
        } else if (which_terms == "global") {
          number_of_plots <- length(global_GAM_term_ids)
        }
        if (number_of_plots > 25) number_of_plots <- 25
        nrow <- ncol <- ceiling(sqrt(number_of_plots))
        if (number_of_plots <- (nrow-1L)*ncol) nrow <- nrow - 1L
        par(mfrow = c(nrow, ncol))
        gamplot_ctrl[["x"]] <- x$gamm
        if (is.numeric(ylim) && length(ylim) == 2L) gamplot_ctrl[["ylim"]] <- ylim
        if (which_terms != "global") {
          terminal_nodes <- c()  
          for (i in 1:length(x$tree)) {
            if (is.null(x$tree[[i]]$node$kids)) terminal_nodes <- c(terminal_nodes, i)
          }
          for (i in 1:length(local_GAM_term_ids)) {
            gamplot_ctrl[["select"]] <-  local_GAM_term_ids[i]
            gamplot_ctrl[["main"]] <- paste("node", terminal_nodes[i])
            if (i == 1L && !is.null(ylim)) {
              tmp <- do.call(plot, gamplot_ctrl)
              if (ylim[1L] == "firstnode") {
                c(min(tmp[[1L]]$fit - .5*tmp[[1L]]$se.mult*tmp[[1L]]$se), 
                  max(tmp[[1L]]$fit + .5*tmp[[1L]]$se.mult*tmp[[1L]]$se))
              }
            } else {
              do.call(plot, gamplot_ctrl)        
            }
          }
        }
        if (which_terms != "local" && length(global_GAM_term_ids) > 0L) {
          for (i in 1:length(global_GAM_term_ids)) {
            gamplot_ctrl[["select"]] <-  global_GAM_term_ids[i]
            gamplot_ctrl[["main"]] <- paste("global term:", smooth_term_labels[global_GAM_term_ids[i]])
            do.call(plot, gamplot_ctrl)        
          }
        }
      } else { ## joint is FALSE
        ## TODO: Prepare par(mfrow = c(.., ..))
        if (which_terms != "local") {
          gamplot_ctrl[["x"]] <- x$gamm
          do.call(plot, gamplot_ctrl)     
        }
        if (which_terms != "global") {
          terminal_nodes <- c()
          for (i in 1:length(x$tree)) {
            if (is.null(x$tree[[i]]$node$kids)) terminal_nodes <- c(terminal_nodes, i)
          }
          for (i in terminal_nodes) {
            gamplot_ctrl[["x"]] <- x$tree[[i]]$node$info$object
            gamplot_ctrl[["main"]] <- paste("node", i)
            do.call(plot, gamplot_ctrl)   
          }
        }
      }
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
#' ## GAM tree without global terms:
#' gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$specimen)
#' coef(gt1)
#' 
#' ## GAM tree with global terms:
#' gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
#'               data = eco, cluster = eco$specimen)
#' coef(gt2)
#' coef(gt2, which = "global")
#' @export
#' @importFrom stats coef
coef.gamtree <- function(object, which = "local", ...) {
  if (is.null(object$gamm)) {
    coefs <- coef(object$tree)
  } else {
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
  }
  return(coefs)
}



#' Predictions from fitted GAM tree
#'  
#' Takes a fitted GAM tree (produced by \code{gamtree()}) and produces 
#' predictions given a new set of values for the model covariates,
#' of for the original covariate values used for fitting the GAM tree.
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
#' @param ... further arguments to be passed to \code{predict.gam}.
#' 
#' @return Returns a vector of predicted values.
#'  
#' @importFrom graphics par plot
#' @export
predict.gamtree <- function(object, newdata = NULL, type = "link", ...) {

  if (is.null(newdata)) {
    newdata <- object$data
  } else {
    if (type != "node") {
      ## Reconstruct .tree (and .offset) variables
      newdata$.tree <- factor(predict(object$tree, newdata = newdata,
                                      type = "node"), 
                              levels = levels(object$data$.tree))
      joint <- is.null(object$call$joint)
      if (!joint) joint <- object$call$joint
      if (!joint) {
        newdata$.offset <- predict(object$tree, newdata = newdata, 
                                   type = "link")
      }
    }
  }
  if (type == "node") {
    factor(predict(object$tree, newdata = newdata, 
                   type = "node"), levels = levels(object$data$.tree))
  } else if (is.null(object$gamm)) {
    predict(object$tree, newdata = newdata, type = type, ...)
  } else {
    predict(object$gamm, newdata = newdata, type = type, ...)
  }
  ## TODO: Implement possibility of excluding local and/or global effects
  ## TODO: Implement support for mgcv's gam.predict option type = "terms" or "iterms", 
  ## and standard errors can be requested with predictions through use of 
  ## argument "se.fit" 
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
