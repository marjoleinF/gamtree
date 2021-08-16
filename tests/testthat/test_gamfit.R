## Full model to estimate is: ff <- response ~ s(time, k = 4) | treatment
tf <- response ~ time  | noise1 + noise2 + treatment # tree formula
lgf <- response ~ s(time, k = 4)  # local GAM formula
method <- "REML"
offset <- 100
tree <- mob(tf, data = rats, local_gam_form = lgf, alpha = .20,
            fit = gamfit, control = mob_control(xtype = "data.frame", 
                                                ytype = "data.frame"))  
offset_tree <- mob(tf, data = rats, local_gam_form = lgf, alpha = .20, 
                   offset = rep(offset, times = nrow(rats)), method = method,
                   family = "scat",
                   fit = gamfit, control = mob_control(xtype = "data.frame", 
                                                       ytype = "data.frame"))  

#saveRDS(tree, file = "prev_res/tree.RDS")
#saveRDS(offset_tree, file = "prev_res/offset_tree.RDS")

context("Test gamfit function")

test_that("does gamfit return the right tree structure and coefficients?", {
  
  expect_equal(length(tree), 3L, tolerance = 1.49e-08)
  
  expect_equal(coef(tree), coef(readRDS("prev_res/tree.RDS")), tolerance = 1.49e-08) 
  
  expect_equal(tree[[1]]$node$split, readRDS("prev_res/tree.RDS")[[1]]$node$split, 
               tolerance = 1.49e-08)
  
  expect_equal(tree[[1]]$node$info$coefficients, 
               readRDS("prev_res/tree.RDS")[[1]]$node$info$coefficients, 
               tolerance = 1.49e-08)
  
  expect_equal(tree[[2]]$node$info$coefficients, 
               readRDS("prev_res/tree.RDS")[[2]]$node$info$coefficients, 
               tolerance = 1.49e-08)
  
  expect_equal(tree[[3]]$node$info$coefficients, 
               readRDS("prev_res/tree.RDS")[[3]]$node$info$coefficients, 
               tolerance = 1.49e-08)
})

test_that("does gamfit pass gam() arguments correctly?", {
  
  ## Test offset
  expect_equal(offset_tree[[1]]$node$info$object$offset, 
               rep(100L, times = nrow(rats)), tolerance = 1.49e-08)
  
  ## Test family
  expect_equal(offset_tree[[1]]$node$info$object$family, 
               readRDS("prev_res/offset_tree.RDS")[[1]]$node$info$object$family, 
               tolerance = 1.49e-08)
  
  ## Test method
  expect_equal(offset_tree[[1]]$node$info$object$method, "REML")
  expect_equal(tree[[1]]$node$info$object$method, "GCV")
  expect_equal(tree[[1]]$node$info$object$method, "GCV")
  
})