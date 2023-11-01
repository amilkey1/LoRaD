devtools::load_all("../")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

source("R/lorad-estimate.R")
source("R/lorad-calc-log-sum.R")
source("R/lorad-standardize-estimation-sample.R")
source("R/lorad-standardize.R")
source("R/lorad-transform.R")

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