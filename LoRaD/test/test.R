# print out the current working directory
library(devtools)
library(expm)
print(getwd())

# load loRad (specify .. because that is where the R directory is located
devtools::load_all("../")

# create a data frame holding the parameter sample
params <- read.table('mcmc-samples.txt', header=TRUE)

# create a data frame holding the column specifications
colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")

# run the transform function in the lorad package, providing the params and colspec data frames
#df <- transform(params, colspec)
#X <- standardize(df)
#print(X)
#print(cov(X))
#print(colMeans(X))
#print(df_list$means)
#print(df_list$inv_matrix%*%df_list$matrix)
#print(df_list$sqrt_matrix%*%df_list$sqrt_matrix)
#print(df_list$inv_matrix)

z <- lorad(params, colspec, .05, "random")
#print(z)
# print out the resulting transformed samples
#print(df)
#print(df_list)

# Here is the last line in mcmc-samples.txt
# 1000000	-459.882711187	0.207259314	3.043891290
# 1000000 is really the 10000th sample because we saved every 100

# Here is the last line of the output
# 10000  -460.3434  -1.5737845 1.1131367
#  -460.3434 = -459.882711187 + ln(0.207259314) + ln(3.043891290)
# -1.5737845 = ln(0.207259314)
#  1.1131367 = ln(3.043891290)
