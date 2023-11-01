setwd("R")
devtools::load_all("../")
setwd("..")
library(testthat)
library(lorad)

#source("R/lorad-estimate.R")
#source("R/lorad-calc-log-sum.R")
#source("R/lorad-standardize-estimation-sample.R")
#source("R/lorad-standardize.R")
#source("R/lorad-transform.R")

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
  
  # test for any error
  expect_(lorad_estimate(params, colspec, 0.5, "left", 0.1), "colspec has 0 length")
  
  
})


test_that("training frac is not between 0 and 1", {
  # Create a data frame holding the parameter sample
  params <- read.table('test/k80-samples.txt', header=TRUE)
  
  # Create a data frame holding the column specifications
  colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")
  
  # test for any error
  expect_error(lorad_estimate(params, colspec, 5, "left", 0.1))
})

#params, colspec, training_frac, training_mode, coverage