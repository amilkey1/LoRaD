# print out the current working directory
library(devtools)
print(getwd())

# load loRad (specify .. because that is where the R directory is located
devtools::load_all("../")

# create a data frame holding the parameter sample
params <- read.table('k80-samples.txt', header=TRUE)

# create a data frame holding the column specifications
colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")

# run the transform function in the lorad package, providing the params and colspec data frames
z <- lorad(params, colspec, .5, "left", .5, 'k80-samples.txt')