test_that("MDI is 0 for a unique leaf", {
  leaf <- make_leaf()
  expect_equal(MDI_tree(leaf, 'X1', 10), 0)
})

test_that("Correct MDI for a 1-depth tree", {
  left <- make_leaf()
  right <- make_leaf()
  root <- make_node(n=10, pval=0.2, var='X1', left, right)
  
  # Correct covariate
  # Contribution = n/n_root*(1-p) = 10/10*(1-0.2)= 0.8
  expect_equal(MDI_tree(root, 'X1', n_root=10), 0.8)
  
  # Different covariate
  expect_equal(MDI_tree(root, 'X2', n_root=10), 0)
})

test_that("Correct MDI for multiple-depth tree", {
  left_leaf <- make_leaf()
  right_leaf <- make_leaf()
  left_node <- make_node(n=6, pval=0.4, var='X2', left_leaf, right_leaf)
  right_leaf_root <- make_leaf()
  root <- make_node(n=10, pval= 0.1, var='X1', left_node, right_leaf_root)
  
  expected_X1 <- 10/10*(1-0.1)
  expected_X2 <- 6/10*(1-0.4)
  expected_X3 <- 0
  
  expect_equal(MDI_tree(root, 'X1', 10), expected_X1)
  expect_equal(MDI_tree(root, 'X2', 10), expected_X2)
  expect_equal(MDI_tree(root, 'X3', 10), expected_X3)
})

test_that("Correct MDI for a forest", {
  tree1 <- make_node(n=10, pval=0.2, var='X1', 
                     left=make_leaf(), right=make_leaf())
  tree2 <- make_node(n=10, pval=0.5, var='X1', 
                     left=make_leaf(), right=make_leaf())
  forest <- list(tree1, tree2)
  expected <- mean(c(0.8, 0.5))
  expect_equal(MDI(forest, 'X1'), expected)
})

test_that("Correct MDI for a forest w/ different n_roots", {
  tree1 <- make_node(n=10, pval=0.2, var = 'X1',
                     left = make_leaf(), right = make_leaf())   
  right_subtree <- make_node(n=12, pval=0.2, var = 'X1',
                             left = make_leaf(), right = make_leaf())
  tree2 <- make_node(n=20, pval=0.1, var = 'X2', 
                     left = make_leaf(), right = right_subtree)
  forest <- list(tree1, tree2)
  
  expected <- mean(c(0.8, 12/20*0.8))
  expect_equal(MDI(forest, "X1"), expected)
})

test_that("Correct sorting of MDI", {
  tree1 <- make_node(n=10, pval=0.1, var='X1', 
                     left=make_leaf(), right=make_leaf())
  forest <- list(tree1)
  
  covariates <- data.frame(X1 = rnorm(10), X2 = rnorm(10))

  result <- MDI_all(forest, covariates)
  
  expect_named(result, c('X1', 'X2'))
  expect_true(result[['X1']] > result[['X2']])
  expect_equal(unname(result[['X2']]), 0)
  expect_equal(result, sort(result, decreasing = TRUE))
})