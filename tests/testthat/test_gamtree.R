## TODO: Implement tests for weight passing, to tree and full GAM
## TODO: check if we can use variety of gam families


## Fit GAM trees w/o global terms
gt1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco, cluster = eco$specimen)
#saveRDS(coef(gt1), "prev_res/coef_gt1.RDS")
saveRDS(gt1$tree[[1]]$node$info$test, "prev_res/gt1_teststat.RDS")
saveRDS(gt1$tree[[1]]$node$info$coefficients, "prev_res/tree_coefs_gt1.RDS")
saveRDS(gt1$gamm$coefficients, "prev_res/gamm_coefs_gt1.RDS")
saveRDS(predict(gt1), "prev_res/preds_gt1.RDS")

gt1.1 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, data = eco)

gt2 <- gamtree(Pn ~ s(PAR, k = 5L) | Species, method = "REML", data = eco, 
               cluster = eco$specimen, joint = FALSE)
#saveRDS(coef(gt2), "prev_res/coef_gt2.RDS")
#saveRDS(predict(gt2), "prev_res/preds_gt2.RDS")
#saveRDS(predict(gt2), "prev_res/preds_gt2.RDS")

## Fit GAM trees with global effects
gt3 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
               method = "ML", data = eco, cluster = eco$specimen)
#saveRDS(coef(gt3), "prev_res/coef_gt3.RDS")
#saveRDS(coef(gt3, "global"), "prev_res/coef_gt3_global.RDS")
#saveRDS(predict(gt3), "prev_res/preds_gt3.RDS")
#saveRDS(gt3$gamm$coefficients, "prev_res/gamm_coefs_gt3.RDS")
#saveRDS(predict(gt3), "prev_res/preds_gt3.RDS")

gt4 <- gamtree(Pn ~ s(PAR, k = 5L) | s(cluster_id, bs = "re") + noise | Species, 
               method = "REML", data = eco, cluster = eco$specimen, joint = FALSE,
               family = scat, offset = rep(100, times = nrow(eco)))
#saveRDS(coef(gt4), "prev_res/coef_gt4.RDS")
#saveRDS(gt4$tree[[1]]$node$info$test, "prev_res/gt4_teststat.RDS")
#saveRDS(coef(gt4, "global"), "prev_res/coef_gt4_global.RDS")
#saveRDS(gt4$tree[[1]]$node$info$coefficients, "prev_res/tree_coefs_gt4.RDS")
#saveRDS(predict(gt4), "prev_res/preds_gt4.RDS")

context("Test gamtree and associated methods")


test_that("are tree- and GAM formulas well constructed by gamtree?", {
  
  ## Full GAM formulas
  expect_equal(gt1$gamm_form, formula(Pn ~ 0 + .tree + s(PAR, k = 5, by = .tree))) 
  expect_equal(gt2$gamm_form, NULL)
  expect_equal(gt3$gamm_form, formula(Pn ~ .tree + s(PAR, k = 5, by = .tree) + 
                                        s(cluster_id, bs = "re") + noise - 1))
  expect_equal(gt4$gamm_form, formula(Pn ~ s(cluster_id, bs = "re") + noise))
  
  ## Node-specific GAM formulas
  expect_equal(gt1$tree[[1]]$node$info$object$call$formula, 
               formula(Pn ~ s(PAR, k = 5L)))
  expect_equal(gt2$tree[[1]]$node$info$object$call$formula, 
               formula(Pn ~ s(PAR, k = 5L)))
  expect_equal(gt3$tree[[1]]$node$info$object$call$formula, 
               formula(Pn ~ s(PAR, k = 5L)))
  expect_equal(gt4$tree[[1]]$node$info$object$call$formula, 
               formula(Pn ~ s(PAR, k = 5L)))
})


test_that("does gamtree pass arguments correctly to gam()?", {
  
  ## Check method
  expect_equal(gt2$tree[[1]]$node$info$object$method, "REML")
  expect_equal(gt3$tree[[1]]$node$info$object$method, "ML")

  ## Check family
  expect_equal(gt3$gamm$family$family, "gaussian")
  expect_equal(gt3$tree[[1]]$node$info$object$family$family, "gaussian")
  expect_equal(gt4$gamm$family$family, "Scaled t(Inf,1.199)")
  expect_equal(gt4$tree[[1]]$node$info$object$family$family, "Scaled t(47.389,1.365)")

  ## Check offset
  expect_equal(sum(gt3$gamm$offset), 0L, tolerance = 1.49e-08)
  expect_equal(unname(rowSums(predict(gt3$gamm, type = "terms")[ , c("noise", "s(cluster_id)")])),
               gt3$tree[[1]]$node$info$object$offset, tolerance = 1.49e-08)
  expect_equal(gt4$gamm$offset, 
               unname(100 + predict(gt4$tree, newdata = gt4$data, type = "response")),
               tolerance = 1.49e-08)
})


test_that("does gamtree yield the right tree structure and pass arguments correctly to mob()?", {
  
  ## Check tree size
  expect_equal(length(gt1$tree), 3L, tolerance = 1.49e-08)
  expect_equal(length(gt1.1$tree), 9L, tolerance = 1.49e-08)
  
  ## Check param stab tests
  expect_equal(gt1$tree[[1]]$node$info$test, readRDS("prev_res/gt1_teststat.RDS"), 
               tolerance = 1.49e-08)
  expect_equal(gt4$tree[[1]]$node$info$test, readRDS("prev_res/gt4_teststat.RDS"), 
               tolerance = 1.49e-08)
  
  ## Check tree coefficients
  expect_equal(gt1$tree[[1]]$node$info$coefficients, readRDS("prev_res/tree_Coefs_gt1.RDS"),
               tolerance = 1.49e-08)
  expect_equal(gt4$tree[[1]]$node$info$coefficients, readRDS("prev_res/tree_Coefs_gt4.RDS"),
               tolerance = 1.49e-08)
  
  ## Check global coefficients
  expect_equal(gt1$gamm$coefficients, readRDS("prev_res/gamm_coefs_gt1.RDS"),
               tolerance = 1.49e-08)
  expect_equal(gt3$gamm$coefficients, readRDS("prev_res/gamm_coefs_gt3.RDS"),
               tolerance = 1.49e-08)
  
})


test_that("test coef.gamtree", {
  
  ## Check local coefficients
  expect_equal(coef(gt1), readRDS("prev_res/coef_gt1.RDS"), tolerance = 1.49e-08)
  expect_equal(coef(gt2), readRDS("prev_res/coef_gt2.RDS"), tolerance = 1.49e-08)
  expect_equal(coef(gt3), readRDS("prev_res/coef_gt3.RDS"), tolerance = 1.49e-08)
  expect_equal(coef(gt4), readRDS("prev_res/coef_gt4.RDS"), tolerance = 1.49e-08)
  
  ## Check global coefficients
  expect_warning(expect_null(coef.gamtree(gt1, "global")))
  expect_warning(expect_null(coef.gamtree(gt2, "global")))
  expect_equal(coef(gt3, which = "global"), readRDS("prev_res/coef_gt3_global.RDS"), tolerance = 1.49e-08)
  expect_equal(coef(gt4, which = "global"), readRDS("prev_res/coef_gt4_global.RDS"), tolerance = 1.49e-08)
})



test_that("test predict.gamtree", {
  
  expect_equal(predict(gt1), readRDS("prev_res/preds_gt1.RDS"), tolerance = 1.49e-08)
  expect_equal(predict(gt2), readRDS("prev_res/preds_gt2.RDS"), tolerance = 1.49e-08)
  expect_equal(predict(gt3), readRDS("prev_res/preds_gt3.RDS"), tolerance = 1.49e-08)
  expect_equal(predict(gt4), readRDS("prev_res/preds_gt4.RDS"), tolerance = 1.49e-08)
  
})