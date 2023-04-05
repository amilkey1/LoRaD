# print out the current working directory
library(devtools)
library(expint)	
library(expm)
print(getwd())

# Load loRad (specify .. because that is where the R directory is located
devtools::load_all("../")

# Create a data frame holding the parameter sample
params <- read.table('jc-samples.txt', header=TRUE)

# Create a data frame holding the column specifications
colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="ignore")

# Run the lorad function in the loRad package, providing 
# the params and colspec data frames as well as specifying
# the training fraction, the training sample selection method,
# and the coverage fraction
lorad(params, colspec, 0.5, "left", 0.5, "jc-samples.txt")
