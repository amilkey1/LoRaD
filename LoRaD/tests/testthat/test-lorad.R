setwd("R")
devtools::load_all("../")
setwd("..")
library(testthat)
library(lorad)

test_that("colspec is set up correctly", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for no error
  expect_no_error(lorad_estimate(params, colspec, 0.5, "left", 0.1))
})


test_that("colspec is missing a param", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive")
  
  # test for any error
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1))
})


test_that("colspec is empty", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c()
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "colspec has 0 length")
})


test_that("training frac is not between 0 and 1", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 5, "left", 0.1), "training fraction must be between 0 and 1 but is 5")
})

test_that("colspec is missing a param", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "colspec does not match column names - colspec is too short")
})

test_that("colspec is too long", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive", "scooby"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "colspec does not match column names - colspec is too long")
})

test_that("colspec is not specified", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "object 'colspec' not found")
})

test_that("colspec has typos", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="itiration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "unknown colspec found: itiration")
})

test_that("colspec has typos in param names", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "logkernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "colspec does not match column names")
})

test_that("coverage frac is not between 0 and 1", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, .5, "left", 3), "coverage fraction must be between 0 and 1 but is 3")
})

test_that("coverage frac is below 0", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, .5, "left", -3), "coverage fraction must be between 0 and 1 but is -3")
})

test_that("coverage frac is a string", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, .5, "left", ".5"), "coverage fraction must be numeric")
})

test_that("training frac is a string", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, ".5", "left", .5), "training fraction must be numeric")
})

test_that("training mode is unknown", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, .5, "scooby", .5), "Unknown training mode scooby")
})

test_that("simplexfinal specified without simplex columns", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/GTR-I-G-samples.log', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("Iteration"="iteration", "Posterior"="posterior", "Likelihood"="ignore", "Prior"="ignore", "alpha"="positive", "edge_length_proportions.1."="simplexfinal", "edge_length_proportions.2."="simplex", "edge_length_proportions.3."="simplex", "edge_length_proportions.4."="simplex", "edge_length_proportions.5."="simplex", "edge_length_proportions.6."="simplex", "edge_length_proportions.7."="simplexfinal", "edgelens.1."="ignore", "edgelens.2."="ignore", "edgelens.3."="ignore", "edgelens.4."="ignore", "edgelens.5."="ignore", "edgelens.6."="ignore", "edgelens.7."="ignore", "er.1."="simplex", "er.2."="simplex", "er.3."="simplex", "er.4."="simplex", "er.5."="simplex", "er.6."="simplexfinal", "pi.1."="simplex", "pi.2."="simplex", "pi.3."="simplex", "pi.4."="simplexfinal", "pinvar"="proportion", "site_rates.1."="ignore", "site_rates.2."="ignore", "site_rates.3."="ignore", "site_rates.4."="ignore", "tree_length"="positive")
  
  # test for error message
  expect_error(lorad_estimate(params, colspec, .5, "left", .5), "simplexfinal column found with no preceding simplex columns")
})

test_that("no variation in params", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/jc-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for no error
  expect_error(lorad_estimate(params, colspec, 0.5, "left", 0.1), "there is not enough parameter variation to estimate the marginal likelihood")
})
