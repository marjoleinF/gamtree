utils::globalVariables(c(".tree", ".offset", ".global", ".weights", ".cluster", "current_offset"))

## TODO: Check summary method for all joint and global combinations
## TODO: Check predict method for all joint and global combinations
## TODO: Check plot method for all joint and global combinations
## TODO: Check print method for all joint and global combinations


#' Recursively partition a dataset based on a (local) GAM, while accounting for 
#' global GAM terms.
#'
#' \code{gamtree} recursively partition a dataset based on a (local) GAM, in
#' which a smooth function (spline) of a predictor variable is estimated, 
#' while accounting for globally specified smooth functions and random effects.
#' 
#' @param formula specifies the model formula, consisting of three or four 
#' parts: the response variable followed by a tilde ('~'); the terms for the 
#' node-specific (local) GAMs, followed by a vertical bar ('|'); the global GAM terms followed 
#' by a vertical bar (this global part is optional!), and then the potential 
#' partitioning variables (separated by a '+'). The 'by' argument of function
#' \code{\link[mgcv]{s}} may NOT be used in the node-specific GAM formulation. Refrain from
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
#' observations and the partitioning variables are time-invariant covariates.
#' @param joint Should local models be re-estimated jointly along with the
#' global model? Jointly re-estimating may yield a more accurate final model,
#' and may need less iterations to converge. A possible advantage of joint 
#' re-estimation of the local, node-specific models is that the same smoothing and 
#' scale parameters will be employed within each smooth term. If the local spline 
#' models are estimated separately within each node, they will obtain a different
#' scale and smoothing parameter in each node. It however depends on the data
#' problem at hand if smoothing and scale parameters should be estimated 
#' locally (i.e., separately for each node), or globally.
#' @param globalstart Should estimation initialize with the global model? Defaults 
#' \code{FALSE}, resulting in the tree (partition or subgroup structure) being 
#' estimated, first. If set to \code{TRUE}, (ignored if no global part has been
#' specified in the \code{formula}), the global part of the model will be estimated
#' first.  
#' @param maxit numeric. The maximum number of iterations to be performed in 
#' estimation of the model tree.
#' @param abstol numeric. Specifies the convergence criterion. If the log-likelihood 
#' values of the full GAM (with both local and global components), differ less than 
#' abstol between two consecutive iterations, estimation is halted.
#' @param verbose logical. Should progress be printed to the commande line in 
#' every iteration? If true, the iteration number, information on the 
#' splitting procedure, and the log-likelihood (with df) value of the fitted 
#' full mixed-effects gam model is printed.
#' @param bam logical. Should function \code{\link[mgcv]{bam}} be used to fit
#' the GAMs? Defaults to \code{FALSE}, resulting in function \code{\link[mgcv]{gam}} 
#' being used for fitting the GAMs. \code{TRUE} will result in function 
#' \code{\link[mgcv]{bam}} being used for fitting the GAMs.
#' @param alt_formula list with two elements, for specifying non=standard model formulea
#' for GAM. E.g., the formula list required for use of the \code{\link[mgcv]{multinom}}
#' family.
#' @param gam_ctrl a list of fit control parameters to replace defaults returned by 
#' \code{\link[mgcv]{gam.control}}. Values not set assume default values.
#' @param mob_ctrl a list with control parameters as returned by 
#' \code{\link[partykit]{mob_control}} to be passed to function 
#' \code{\link[partykit]{mob}}. Note: argument `xtype` is set to 
#' `data.frame`, by default, and cannot be adjusted.
#' @param ... additional arguments to be passed to function \code{\link[mgcv]{gam}} 
#' (or \code{\link[mgcv]{bam}}). The same arguments will be employed for the local
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
                    abstol = 0.001, verbose = FALSE, maxit = 100, 
                    bam = FALSE, gam_ctrl = list(), mob_ctrl = mob_control(), 
                    alt_formula = NULL, ...) {

  ## TODO: Argument checking
  
  
  ## specify use of gam() versus bam()
  if (bam) {
    mob_gamfit <- bamfit
    gam_fit <- mgcv::bam
  } else {
    mob_gamfit <- gamfit # TODO: Replace bam_fit function
    gam_fit <- mgcv::gam
  }
  ## Construct formulas
  ##
  ## For tree (tf), local gam (lgf) and global gam (ggf):
  ##  - tf is for passing response, predictors and part vars to mob()
  ##  - lgf is for fitting node-specific GAM in the nodes
  ##  - ggf is an intermediate formula, used for creating full GAM formula later
  ##  - fgf is for fitting full GAM
  ff <- as.Formula(formula)
  lgf <- formula(ff, lhs = NULL, rhs = 1)
  if (any(grepl("by=", gsub(" ", "", as.character(lgf))))) warning(
    "It appears that the 'by' argument was used in the node-specific GAM specification. This may yield unexpected results if joint=TRUE."
    )    
  global_gam <- all(length(Formula(ff)) == c(1, 3))
  if (global_gam) {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 3))
    ggf <- formula(ff, lhs = NULL, rhs = 2)
    global_vars <- all.vars(ggf[[3]])
  } else {
    local_vars <- all.vars(formula(ff, lhs = 0, rhs = 1))
    part_vars <- all.vars(formula(ff, lhs = 0, rhs = 2))
    if (joint) {
      ggf <- formula(ff, lhs = NULL, rhs = 0)
    }
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
  ##  - fgf_alt is for fitting the full GAM if mob() did not split
  ##  - If joint = FALSE, then fgf = lgf but with by = .tree added, and fgf_alt = lgf
  new_alt <- new <- attr(terms(lgf), "term.labels")
  if (length(new) > 1L) warning("It appears more than 1 term has been specified for the node-specific (local) GAM. This may lead to unexpected results.")

  if (global_gam) {
    new <- gsub(")", ", by = .tree)", new)
    new[!grepl(")", new)] <- paste0(new[!grepl(")", new)], ":.tree") 
    new <- paste(new, collapse = " + ")
    new <- paste0("~ 0 + .tree + ", new, " + .")
    new_alt <- paste(new_alt, collapse = " + ")
    new_alt <- paste0("~ ", new_alt, "+ .")
    if (joint) {
      fgf <- update(old = ggf, new = formula(new))
      ## TODO: This should not omit intercept, but does in fgf_alt for 
      ## formula <- cbind(succ, fail) ~ s(PAR) | Species
      ## Only for cbind response? Only for local-only gamtree?
      fgf_alt <- update(old = ggf, new = formula(new_alt))
    } else {
      fgf <- fgf_alt <- ggf 
    }
  } else if (joint) {
    ## Global formula equal local formula, plus by = .tree for the smooth and 
    ##  separate factor .tree added
    new <- gsub(")", ", by = .tree)", new)
    new <- paste0("~ 0 + .tree + ", new)
    fgf <- formula(paste0(response, new))
    ## Global formula should equal local formula:
    fgf_alt <- lgf   
  }
  
  if (!is.null(alt_formula)) {
    lgf <- fgf_alt <- alt_formula[[1]]
    fgf <- alt_formula[[2]]
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
  data$.weights <- if (is.null(weights)) 1 else weights
  if (!is.null(cluster)) data$.cluster <- cluster
  ## TODO: Allow for cluster argument to refer to one of the columns of data
  
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
  data$.offset <- if (is.null(offset)) 0L else offset
  ## TODO: Allow for starting with estimation of global model
  if (global_gam && globalstart) {
    gamm <- gam_fit(fgf_alt, data = data, weights = .weights, 
                    offset = data$.offset, gam_ctrl = gam_ctrl, ...)
    data$.global <- predict(gamm, newdata = data) + data$.offset
  } else {
    data$.global <- data$.offset    
  }
  newloglik <- -Inf
  oldloglik <- c(-Inf, -Inf) # 2nd element is the oldest
  iteration <- 0L
  continue <- TRUE
  
  ## iterate between tree and global model estimation
  while (continue) {
    
    iteration <- iteration + 1L
    if (verbose) print(paste("iteration", iteration))
    
    ## grow tree and get node memberships
    if (is.null(cluster)) {
      tree <- mob(tf, data = data, local_gam_form = lgf, fit = mob_gamfit, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, gam_ctrl = gam_ctrl, ...)   
    } else {
      tree <- mob(tf, data = data, local_gam_form = lgf, fit = mob_gamfit, 
                  weights = .weights, offset = .global, 
                  control = mob_ctrl, cluster = .cluster, gam_ctrl = gam_ctrl, ...)       
    }
    if (joint || global_gam) {
      data$.tree <- if (joint) {
        factor(predict(tree, newdata = data, type = "node"))
      } else {
        predict(tree, newdata = data, type = "link")
      }
    
      ## if joint=FALSE and global_gam=TRUE: offset should be data$.offset + data$.tree
      ## if joint=FALSE and global_gam=FALsE: we did not enter this loop anyway
      ## if joint=TRUE and global_gam=TRUE: offset should be data$.offset
      ## if joint=TRUE and global_gam=FALSE: offset should be data$.offset
      ##
      ## Fit full GAM with global and local parts:
      if (joint) {
        if (length(tree) > 1L) {
          gamm <- gam_fit(fgf, data = data, weights = .weights, offset = .offset, 
                          gam_ctrl = gam_ctrl, ...)          
        } else {
          gamm <- gam_fit(fgf_alt, data = data, weights = .weights, offset = .offset, 
                          gam_ctrl = gam_ctrl, ...)
        }
      } else {
        data$current_offset <- data$.offset + data$.tree
        ## Note: gam() only accepts unquoted colnames from data for the offset argument  
        if (length(tree) > 1L) {
          gamm <- gam_fit(fgf, data = data, weights = .weights, offset = current_offset,
                          gam_ctrl = gam_ctrl, ...)          
        } else {
          gamm <- gam_fit(fgf_alt, data = data, weights = .weights, offset = current_offset, 
                          gam_ctrl = gam_ctrl, ...)
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
                                        terms = keep_terms1)) + data$.offset
      } else {
        data$.global <- rowSums(predict(gamm, newdata = data, type = "terms",
                                        terms = keep_terms)) + data$.offset
      }
      
      ## iteration information
      newloglik <- logLik(gamm)    
      continue <- (abs(newloglik - oldloglik[1]) > abstol) &
        (iteration < maxit)
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
    gamm_form = if (global_gam || joint) fgf else NULL,
    call = cl,
    tree = tree,
    gamm = if (global_gam || joint) gamm else NULL,
    data = data,
    loglik = if (global_gam || joint) as.numeric(newloglik) else NULL, # TODO: return a vector if no global model
    df = if (global_gam || joint) attr(newloglik, "df") else NULL, # TODO: return a vector if no global model
    iterations = if (global_gam) iteration else NULL, 
    abstol = abstol,
    joint = joint
  )
  class(result) <- "gamtree"
  return(result)
}



#############################################################
##
## Fitting function for gam-based recursive partitioning with MOB.
##
gamfit <- function (y, x, start = NULL, weights = NULL, offset = NULL, 
                    ..., local_gam_form, gam_ctrl = list())
{

  args <- list(...)
  
  if (is.null(x)) 
    x <- matrix(1, nrow = NROW(y), ncol = 1L, 
                dimnames = list(NULL, "(Intercept)"))
  
  args <- c(list(formula = local_gam_form, data = cbind(x, y), 
                 weights = weights, offset = offset, control = gam_ctrl), args)
  
  do.call("gam", args)
  
  ## TODO: in ?gamObject we read that deviance is "model deviance (not penalized deviance)".
  ##
  ## TODO: Estimation of par stab tests. From mob() documentation:
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
bamfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, 
                    ..., local_gam_form, gam_ctrl = list()) {
  
  args <- list(...)
  
  if (is.null(x)) 
    x <- matrix(1, nrow = NROW(y), ncol = 1L, 
                dimnames = list(NULL, "(Intercept)"))
  
  args <- c(list(formula = local_gam_form, data = cbind(x, y), 
                 weights = weights, offset = offset, control = gam_ctrl), args)
  
  do.call("bam", args)

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

  if (object$joint) {
    summary(object$gamm, ...)
  } else {
    ## joint=FALSE and global_gam=TRUE
    summary(object$tree, ...)
    summary(object$gamm, ...)
    ## joint=FALSE and global_gam=FALSE
    summary(object$tree, ...)    
    ## TODO: Make nicer summary from the tree and the GAM
    # sum <- summary(object$tree, ...)
    # num_nodes <- length(sum)
    # fam <- data.frame(matrix(NA, nrow = num_nodes, ncol = 3, 
    #                          dimnames = list(names(sum), c("Family", "Link", "Formula"))))
    # 
    # par_terms <- data.frame(matrix(NA, nrow = num_nodes, ncol = 2L, 
    #                         dimnames = list(names(sum), c("Estimate", "Std.Error"))))
    # nonpar_terms <- data.frame(matrix(NA, nrow = num_nodes, ncol = 2L, 
    #                                   dimnames = list(names(sum), c("edf", "Ref.df"))))
    # for (i in 1:num_nodes) {
    #   fam[i, "Family"] <- sum[[i]]$family$family
    #   fam[i, "Link"] <- sum[[i]]$family$link
    #   fam[i, "Formula"] <- deparse(sum[[i]]$formula)
    #   fam[i, "Explained_deviance"] <- sum[[i]]$dev.expl
    #   fam[i, "N"]<- sum[[i]]$n
    #   fam[i, "REML criterion"] <- sum[[i]]$sp.criterion
    #   
    #   ## Get significance of smooth terms. Extract edf. and Ref.df
    #   nonpar_terms[i, "edf"] <- sum[[i]]$edf
    #   nonpar_terms[i, "Ref.df"] <- sum[[i]]$s.table[ , "Ref.df"]
    #   
    #   ## Get parametric coefficients: estimated value and standard error
    #   par_terms[i, "Estimate"] <- sum[[i]]$p.coeff
    #   par_terms[i, "Std.Error"] <- sum[[i]]$se["(Intercept)"]
    # }
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
        coefs[i, -intercept_id] <- tmp[grep(paste0(":.tree"), names(tmp))]
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
