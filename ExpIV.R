## Need to install Rcpp, RcppCNPy, fields packages
library(RcppCNPy)
library(fields)
library(caret)
library(dplyr)
library(tidyverse)
library(GPfit)
library(ggplot2)
library(Rcpp)
library(testthat)
library(nloptr)
library(lattice)
library(kergp)
library(DiceDesign)
library(DiceKriging)
library(gridExtra)

#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)

#Load data in
setwd("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification")

datafile <- read.csv("gpd1.csv", header = T)

# save as tibble
datafile <- base::data.frame(datafile)
#datafile <- tibble::tibble(datafile)

# Replace "rstor" with "vstor" in the paths of Column2
datafile$path <- sub("/mnt/rstor/", "/mnt/vstor/", datafile$path)


#filenames=datafile[,2]
filenames <- datafile[, 2]

y <- datafile[,3]
y <- data.frame(y)

indexall <- datafile[,1]
indexall <- data.frame(indexall)

n=length(filenames) # number of data points
resize_pxl=512 ## for resizing the images to 512*512 pixels for initial speed-up otherwise plotting is slow in my computer
data=matrix(0, n, resize_pxl^2) #initialization

resizeImage = function(im, w.out, h.out) {
  # function to resize an image 
  # im = input image, w.out = target width, h.out = target height
  # Bonus: this works with non-square image scaling.
  
  # initial width/height
  w.in = nrow(im)
  h.in = ncol(im)
  
  # Create empty matrix
  im.out = matrix(rep(0,w.out*h.out), nrow =w.out, ncol=h.out )
  
  # Compute ratios -- final number of indices is n.out, spaced over range of 1:n.in
  w_ratio = w.in/w.out
  h_ratio = h.in/h.out
  
  # Do resizing -- select appropriate indices
  im.out <- im[ floor(w_ratio* 1:w.out), floor(h_ratio* 1:h.out)]
  
  return(im.out)
}

for (i in 1:n){
  temp1=npyLoad(filenames[i])
  temp2=resizeImage(temp1,resize_pxl,resize_pxl) # can delete this if need to keep the original images
  data[i,]=as.vector(temp2) ## converting matrices to vectors and storing as rows
}

## Principal component analysis
pr=prcomp(data, center=F)

# Calculate the proportion of variance explained by each PC
variance_explained <- pr$sdev^2 / sum(pr$sdev^2)

# Calculate the cumulative proportion of variance explained
cumulative_variance <- cumsum(variance_explained)

# Determine how many PCs to keep to explain 95% of the variance
num_pcs_to_keep <- which(cumulative_variance >= 0.95)[1]
num_pcs_to_keep2 <- which(cumulative_variance >= 0.99)[1]
num_pcs_to_keep3 <- which(cumulative_variance >= 0.995)[1]
# Print the results
cat("Number of PCs to keep for 95% variance explained:", num_pcs_to_keep, "\n")
cat("Number of PCs to keep for 99% variance explained:", num_pcs_to_keep2, "\n")
cat("Number of PCs to keep for 99.5% variance explained:", num_pcs_to_keep3, "\n")

set.seed(26011993)
#Define dimensions
d <- 22; 
inputs <- paste("x", 1L:d, sep = "")
#Define 6 Kernels
#1. Matern 5/2 kernel
CovM52 <- covTS(inputs = inputs, d = d, kernel = "k1Matern5_2",
                dep = c(range = "input"))
CovM52

#2. Matern 3/2 kernel
CovM32 <- covTS(inputs = inputs, d = d, kernel = "k1Matern3_2",
                dep = c(range = "input"))
CovM32

#3. Exponential kernel = Matern 1/2
CovExp <- covTS(inputs = inputs, d = d, kernel = "k1Exp",
                dep = c(range = "input"))
CovExp

#4. Squared Exponential kernel = Gaussian
CovGauss <- covTS(inputs = inputs, d = d, kernel = "k1Gauss",
                  dep = c(range = "input"))
CovGauss 

#5. Power Exponential kernel
CovPowExp <- covTS(inputs = inputs, d = d, kernel = "k1PowExp",
                   dep = c(range = "input"))
CovPowExp

#6. Additive Exponential kernel
# Define the exponential kernel for matrices
expKernel <- function(x1, x2, par) {
  # Initialize the covariance matrix
  covMatrix <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))
  
  # Loop over each dimension to calculate the covariance additively
  for (i in 1:ncol(x1)) {
    # Extracting the parameters for each dimension: range and variance
    range = par[2 * i - 1]
    variance = par[2 * i]
    
    # Compute the distance for the current dimension
    distanceMatrix = outer(x1[, i], x2[, i], FUN = function(u, v) abs(u - v))
    
    # Compute the covariance for the current dimension
    covMatrix <- covMatrix + variance * exp(-distanceMatrix / range)
  }
  
  return(covMatrix)
}

# Setup the covariance manager using the exponential kernel
myCovMan <- covMan(
  kernel = expKernel,
  hasGrad = FALSE,  # Set to TRUE if gradient is needed for optimization
  d = 22,  # Updated number of dimensions
  acceptMatrix = TRUE,  # This must be TRUE for handling matrices
  parNames = rep(c("range", "variance"), each = d),
  par = rep(1, 2*d),  # Initial values for each parameter
  parLower = rep(c(range = 0.01, variance = 0.01), d),
  parUpper = rep(c(range = 10, variance = 10), d),
  label = "AdditiveExp"
)

#7. Additive Power Exponential kernel
# Define the powered exponential kernel for matrices
powExpKernel <- function(x1, x2, par) {
  # Initialize the covariance matrix
  covMatrix <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))
  
  # Loop over each dimension to calculate the covariance additively
  for (i in 1:ncol(x1)) {
    # Extracting the parameters for each dimension: range, variance, and power
    range = par[3 * i - 2]
    variance = par[3 * i - 1]
    power = par[3 * i]
    
    # Compute the distance for the current dimension
    distanceMatrix = outer(x1[, i], x2[, i], FUN = function(u, v) abs(u - v))
    
    # Compute the covariance for the current dimension
    covMatrix <- covMatrix + variance * exp(- (distanceMatrix / range) ^ power)
  }
  
  return(covMatrix)
}

# Setup the covariance manager using the powered exponential kernel
myCovManPower <- covMan(
  kernel = powExpKernel,
  hasGrad = FALSE,  # Set to TRUE if gradient is needed for optimization
  d = 22,  # Updated number of dimensions
  acceptMatrix = TRUE,  # This must be TRUE for handling matrices
  parNames = rep(c("range", "variance", "power"), each = d),
  par = rep(1, 3*d),  # Initial values for each parameter
  parLower = rep(c(range = 0.01, variance = 0.01, power = 0.1), d),
  parUpper = rep(c(range = 10, variance = 10, power = 2), d),
  label = "AdditivePowExp"
)

#8. Additive Power Exponential kernel
# Define the powered exponential kernel for matrices with fixed power 1.95
powExpKernelFixed <- function(x1, x2, par) {
  # Initialize the covariance matrix
  covMatrix <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))
  
  # Fixed power parameter
  power <- 1.95
  
  # Loop over each dimension to calculate the covariance additively
  for (i in 1:ncol(x1)) {
    # Extracting the parameters for each dimension: range and variance
    range = par[2 * i - 1]
    variance = par[2 * i]
    
    # Compute the distance for the current dimension
    distanceMatrix = outer(x1[, i], x2[, i], FUN = function(u, v) abs(u - v))
    
    # Compute the covariance for the current dimension
    covMatrix <- covMatrix + variance * exp(- (distanceMatrix / range) ^ power)
  }
  
  return(covMatrix)
}

# Setup the covariance manager using the powered exponential kernel with fixed power
myCovManPowerFixed <- covMan(
  kernel = powExpKernelFixed,
  hasGrad = FALSE,  # Set to TRUE if gradient is needed for optimization
  d = 22,  # Updated number of dimensions
  acceptMatrix = TRUE,  # This must be TRUE for handling matrices
  parNames = rep(c("range", "variance"), each = d),
  par = rep(1, 2*d),  # Initial values for each parameter
  parLower = rep(c(range = 0.01, variance = 0.01), d),
  parUpper = rep(c(range = 10, variance = 10), d),
  label = "AdditivePowExpFixed"
)


# Set seed for reproducibility
set.seed(26011993)
## Inputs for Guassian Process using PCA
X=pr$x[,1:d] ##Reducing dimension from 512*512 to pr_n
colnames(X) <- inputs  # Assuming you want column names as X1, X2, ..., Xd
# Rescale the values to the interval (0, 1)
X <- (X - min(X)) / (max(X) - min(X))
data <- data.frame(X, y = y, indexall = indexall)




# Set a random seed for reproducibility
set.seed(2503)

# Split the Xy into training, testing, and validation sets
train_prop <- 0.65  # Proportion of Xy for training
validation_prop <- 0.15  # Proportion of Xy for validation
test_prop <- 0.20   # Proportion of Xy for testing

# Calculate the number of rows for each set
n <- nrow(data)
n_train <- round(train_prop * n)
n_test <- round(test_prop * n)
n_validation <- n - n_train - n_test

# Create indices for each set
train_indices <- sample(1:n, n_train)
sorttrain <- sort(train_indices)
test_indices <- sample(setdiff(1:n, train_indices), n_test)
sorttest <- sort(test_indices)
validation_indices <- sample(setdiff(1:n, c(train_indices, test_indices)), n_validation)
sortvalidation <- sort(validation_indices)
# Create the subsets
train <- data[sorttrain, ]
test <- data[sorttest, ]
validation <- data[sortvalidation, ]

#Save the index for train, test and validation set in excel to export
library("writexl")
TrainindexAll <- data.frame(sorttrain)
write_xlsx(TrainindexAll, "~/TrainindexAll.xlsx")
TestindexAll <- data.frame(sorttest)
write_xlsx(TestindexAll, "~/TestindexAll.xlsx")
ValidationindexAll <- data.frame(sortvalidation)
write_xlsx(ValidationindexAll, "~/ValidationindexAll.xlsx")





#Fit Matern 5/2 and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.CovM52 <- gp(formula = y ~ 1, data = train,
                   cov = CovM52, noise = TRUE,
                   parCovIni = rep(1, 2*d),
                   parCovLower = c(rep(1e-4, 2*d)),
                   parCovUpper = c(rep(5, d), rep(10,d)))

# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.CovM52

#Fit Matern 3/2 and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.CovM32 <- gp(formula = y ~ 1, data = train,
                   cov = CovM32, noise = TRUE,
                   parCovIni = rep(1, 2*d),
                   parCovLower = c(rep(1e-4, 2*d)),
                   parCovUpper = c(rep(5, d), rep(10,d)))

# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.CovM32

#Fit Exponential and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.CovExp <- gp(formula = y ~ 1, data = train,
                   cov = CovExp, noise = TRUE,
                   parCovIni = rep(1, 2*d),
                   parCovLower = c(rep(1e-4, 2*d)),
                   parCovUpper = c(rep(5, d), rep(10,d)))


# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.CovExp

#Fit Gaussian and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.CovGauss <- gp(formula = y ~ 1, data = train,
                     cov = CovGauss, noise = TRUE,
                     parCovIni = rep(1, 2*d),
                     parCovLower = c(rep(1e-4, 2*d)),
                     parCovUpper = c(rep(5, d), rep(10,d)))


# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.CovGauss

#Fit Power Exponential and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.CovPowExp <- gp(formula = y ~ 1, data = train,
                      cov = CovPowExp, noise = TRUE,
                      parCovIni = rep(1, 2*d+1),
                      parCovLower = c(rep(1e-4, 2*d+1)),
                      parCovUpper = c(rep(5, d), rep(10,d+1)))

# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.CovPowExp

#Fit Additive Exponential and record time taken
# Record the start time
start_time <- Sys.time()

## Creation of a gp object
fitgp.myCovMan <- gp(formula = y ~ 1, data = train,
                     cov = myCovMan, noise = TRUE,
                     parCovIni = rep(1, 2*d),  # Initial parameter values updated for 10 dimensions
                     parCovLower = rep(1e-4, 2*d),  # Lower bounds for parameters updated
                     parCovUpper = rep(c(rep(5, d), rep(10, d)))  # Upper bounds for parameters updated
)
# Record the end time
end_time <- Sys.time()

# Calculate the time difference in seconds
time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours <- floor(time_taken / 3600)
minutes <- floor((time_taken %% 3600) / 60)
seconds <- time_taken %% 60

# Print the result
cat("Time taken: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds\n")
fitgp.myCovMan

# Fit Additive Powered Exponential Kernel with Variable Power and record time taken
# Record the start time
start_time_variable <- Sys.time()

# Creation of a gp object
fitgp.myCovManPower <- gp(
  formula = y ~ 1, data = train,
  cov = myCovManPower, noise = TRUE,
  parCovIni = rep(1, 3*d),  # Initial parameter values updated for 10 dimensions
  parCovLower = rep(1e-4, 3*d),  # Lower bounds for parameters updated
  parCovUpper = rep(c(rep(5, d), rep(10, d), rep(2, d)))  # Upper bounds for parameters updated
)

# Record the end time
end_time_variable <- Sys.time()

# Calculate the time difference in seconds
time_taken_variable <- as.numeric(difftime(end_time_variable, start_time_variable, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours_variable <- floor(time_taken_variable / 3600)
minutes_variable <- floor((time_taken_variable %% 3600) / 60)
seconds_variable <- time_taken_variable %% 60

# Print the result
cat("Time taken for Variable Power: ", hours_variable, " hours, ", minutes_variable, " minutes, ", seconds_variable, " seconds\n")
fitgp.myCovManPower


# Fit Additive Powered Exponential Kernel with Fixed Power and record time taken
# Record the start time
start_time_fixed <- Sys.time()

# Creation of a gp object
fitgp.myCovManPowerFixed <- gp(
  formula = y ~ 1, data = train,
  cov = myCovManPowerFixed, noise = TRUE,
  parCovIni = rep(1, 2*d),  # Initial parameter values updated for 10 dimensions
  parCovLower = rep(1e-4, 2*d),  # Lower bounds for parameters updated
  parCovUpper = rep(c(rep(5, d), rep(10, d)))  # Upper bounds for parameters updated
)

# Record the end time
end_time_fixed <- Sys.time()

# Calculate the time difference in seconds
time_taken_fixed <- as.numeric(difftime(end_time_fixed, start_time_fixed, units = "secs"))

# Convert seconds to hours, minutes, and seconds
hours_fixed <- floor(time_taken_fixed / 3600)
minutes_fixed <- floor((time_taken_fixed %% 3600) / 60)
seconds_fixed <- time_taken_fixed %% 60

# Print the result
cat("Time taken for Fixed Power: ", hours_fixed, " hours, ", minutes_fixed, " minutes, ", seconds_fixed, " seconds\n")
fitgp.myCovManPowerFixed

# RSME for training and testing data
#CovM52 Train
predCovM52.train <- predict(fitgp.CovM52, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMECovM52.train <- sqrt(mean((train$y - predCovM52.train$mean)^2))
RSMECovM52.train
#CovM52 Test
predCovM52.test <- predict(fitgp.CovM52, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMECovM52.test <- sqrt(mean((test$y - predCovM52.test$mean)^2))
RSMECovM52.test

#CovM32 Train
predCovM32.train <- predict(fitgp.CovM32, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMECovM32.train <- sqrt(mean((train$y - predCovM32.train$mean)^2))
RSMECovM32.train
#CovM32 Test
predCovM32.test <- predict(fitgp.CovM32, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMECovM32.test <- sqrt(mean((test$y - predCovM32.test$mean)^2))
RSMECovM32.test

# RSME for training and testing data
#CovExp Train
predCovExp.train <- predict(fitgp.CovExp, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMECovExp.train <- sqrt(mean((train$y - predCovExp.train$mean)^2))
RSMECovExp.train
#CovM52 Test
predCovExp.test <- predict(fitgp.CovExp, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMECovExp.test <- sqrt(mean((test$y - predCovExp.test$mean)^2))
RSMECovExp.test

#CovGauss Train
predCovGauss.train <- predict(fitgp.CovGauss, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMECovGauss.train <- sqrt(mean((train$y - predCovGauss.train$mean)^2))
RSMECovGauss.train
#CovGauss Test
predCovGauss.test <- predict(fitgp.CovGauss, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMECovGauss.test <- sqrt(mean((test$y - predCovGauss.test$mean)^2))
RSMECovGauss.test

#CovPowExp Train
predCovPowExp.train <- predict(fitgp.CovPowExp, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMECovPowExp.train <- sqrt(mean((train$y - predCovPowExp.train$mean)^2))
RSMECovPowExp.train
#CovPowExp Test
predCovPowExp.test <- predict(fitgp.CovPowExp, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMECovPowExp.test <- sqrt(mean((test$y - predCovPowExp.test$mean)^2))
RSMECovPowExp.test

#myCovMan Train
predmyCovMan.train <- predict(fitgp.myCovMan, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMEmyCovMan.train <- sqrt(mean((train$y - predmyCovMan.train$mean)^2))
RSMEmyCovMan.train
#myCovMan Test
predmyCovMan.test <- predict(fitgp.myCovMan, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMEmyCovMan.test <- sqrt(mean((test$y - predmyCovMan.test$mean)^2))
RSMEmyCovMan.test

#myCovManPower Train
predmyCovManPower.train <- predict(fitgp.myCovManPower, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMEmyCovManPower.train <- sqrt(mean((train$y - predmyCovManPower.train$mean)^2))
RSMEmyCovManPower.train
#myCovManPower Test
predmyCovManPower.test <- predict(fitgp.myCovManPower, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMEmyCovManPower.test <- sqrt(mean((test$y - predmyCovManPower.test$mean)^2))
RSMEmyCovManPower.test

#myCovManPowerFixed Train
predmyCovManPowerFixed.train <- predict(fitgp.myCovManPowerFixed, newdata = as.matrix(train[ , inputs]), type = "UK")
RSMEmyCovManPowerFixed.train <- sqrt(mean((train$y - predmyCovManPowerFixed.train$mean)^2))
RSMEmyCovManPowerFixed.train
#myCovManPowerFixed Test
predmyCovManPowerFixed.test <- predict(fitgp.myCovManPowerFixed, newdata = as.matrix(test[ , inputs]), type = "UK")
RSMEmyCovManPowerFixed.test <- sqrt(mean((test$y - predmyCovManPowerFixed.test$mean)^2))
RSMEmyCovManPowerFixed.test

#Plots CovM52
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predCovM52.train$mean[predCovM52.train$mean < 0] <- 0
predCovM52.train$lower95[predCovM52.train$lower95 < 0] <- 0
#Prepare Dataframe
predCovM52.plottrain <- cbind(predCovM52.train$mean, predCovM52.train$lower95, predCovM52.train$upper95, train$y, train$indexall)
predCovM52.plottrain <- data.frame(predCovM52.plottrain)
# Sort the dataframe by column X5
predCovM52.plottrain <- predCovM52.plottrain %>%
  arrange(X5)
predCovM52.plottrain

# Update the ggplot code
plottrainCovM52 <- ggplot(predCovM52.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Mattern 5/2"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.7, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovM52.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predCovM52.test$mean[predCovM52.test$mean < 0] <- 0
predCovM52.test$lower95[predCovM52.test$lower95 < 0] <- 0
#Prepare Dataframe
predCovM52.plottest <- cbind(predCovM52.test$mean, predCovM52.test$lower95, predCovM52.test$upper95, test$y, test$indexall)
predCovM52.plottest <- data.frame(predCovM52.plottest)
# Sort the dataframe by column X5
predCovM52.plottest <- predCovM52.plottest %>%
  arrange(X5)
predCovM52.plottest

# Update the ggplot code
plottestCovM52 <- ggplot(predCovM52.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Prediction Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Mattern 5/2"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern52)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.7, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovM52.plottest)

#Plots CovM32
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predCovM32.train$mean[predCovM32.train$mean < 0] <- 0
predCovM32.train$lower95[predCovM32.train$lower95 < 0] <- 0
#Prepare Dataframe
predCovM32.plottrain <- cbind(predCovM32.train$mean, predCovM32.train$lower95, predCovM32.train$upper95, train$y, train$indexall)
predCovM32.plottrain <- data.frame(predCovM32.plottrain)
# Sort the dataframe by column X5
predCovM32.plottrain <- predCovM32.plottrain %>%
  arrange(X5)
predCovM32.plottrain

# Update the ggplot code
plottrainCovM32 <- ggplot(predCovM32.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Mattern 3/2"
    #title = "Original Values and Predictions with Confidence Bounds for train set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovM32.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predCovM32.test$mean[predCovM32.test$mean < 0] <- 0
predCovM32.test$lower95[predCovM32.test$lower95 < 0] <- 0
#Prepare Dataframe
predCovM32.plottest <- cbind(predCovM32.test$mean, predCovM32.test$lower95, predCovM32.test$upper95, test$y, test$indexall)
predCovM32.plottest <- data.frame(predCovM32.plottest)
# Sort the dataframe by column X5
predCovM32.plottest <- predCovM32.plottest %>%
  arrange(X5)
predCovM32.plottest

# Update the ggplot code
plottestCovM32 <- ggplot(predCovM32.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Mattern 3/2"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovM32.plottest)


#Plots CovExp
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predCovExp.train$mean[predCovExp.train$mean < 0] <- 0
predCovExp.train$lower95[predCovExp.train$lower95 < 0] <- 0
#Prepare Dataframe
predCovExp.plottrain <- cbind(predCovExp.train$mean, predCovExp.train$lower95, predCovExp.train$upper95, train$y, train$indexall)
predCovExp.plottrain <- data.frame(predCovExp.plottrain)
# Sort the dataframe by column X5
predCovExp.plottrain <- predCovExp.plottrain %>%
  arrange(X5)
predCovExp.plottrain

# Update the ggplot code
plottrainCovExp <- ggplot(predCovExp.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Exponential"
    #title = "Original Values and Predictions with Confidence Bounds for train set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovExp.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predCovExp.test$mean[predCovExp.test$mean < 0] <- 0
predCovExp.test$lower95[predCovExp.test$lower95 < 0] <- 0
#Prepare Dataframe
predCovExp.plottest <- cbind(predCovExp.test$mean, predCovExp.test$lower95, predCovExp.test$upper95, test$y, test$indexall)
predCovExp.plottest <- data.frame(predCovExp.plottest)
# Sort the dataframe by column X5
predCovExp.plottest <- predCovExp.plottest %>%
  arrange(X5)
predCovExp.plottest

# Update the ggplot code
plottestCovExp <- ggplot(predCovExp.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Exponential"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovExp.plottest)


#Plots CovGauss
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predCovGauss.train$mean[predCovGauss.train$mean < 0] <- 0
predCovGauss.train$lower95[predCovGauss.train$lower95 < 0] <- 0
#Prepare Dataframe
predCovGauss.plottrain <- cbind(predCovGauss.train$mean, predCovGauss.train$lower95, predCovGauss.train$upper95, train$y, train$indexall)
predCovGauss.plottrain <- data.frame(predCovGauss.plottrain)
# Sort the dataframe by column X5
predCovGauss.plottrain <- predCovGauss.plottrain %>%
  arrange(X5)
predCovGauss.plottrain

# Update the ggplot code
plottrainCovGauss <- ggplot(predCovGauss.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Gaussian"
    #title = "Original Values and Predictions with Confidence Bounds for train set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovGauss.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predCovGauss.test$mean[predCovGauss.test$mean < 0] <- 0
predCovGauss.test$lower95[predCovGauss.test$lower95 < 0] <- 0
#Prepare Dataframe
predCovGauss.plottest <- cbind(predCovGauss.test$mean, predCovGauss.test$lower95, predCovGauss.test$upper95, test$y, test$indexall)
predCovGauss.plottest <- data.frame(predCovGauss.plottest)
# Sort the dataframe by column X5
predCovGauss.plottest <- predCovGauss.plottest %>%
  arrange(X5)
predCovGauss.plottest

# Update the ggplot code
plottestCovGauss <- ggplot(predCovGauss.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Gaussian"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovGauss.plottest)


#Plots CovPowExp
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predCovPowExp.train$mean[predCovPowExp.train$mean < 0] <- 0
predCovPowExp.train$lower95[predCovPowExp.train$lower95 < 0] <- 0
#Prepare Dataframe
predCovPowExp.plottrain <- cbind(predCovPowExp.train$mean, predCovPowExp.train$lower95, predCovPowExp.train$upper95, train$y, train$indexall)
predCovPowExp.plottrain <- data.frame(predCovPowExp.plottrain)
# Sort the dataframe by column X5
predCovPowExp.plottrain <- predCovPowExp.plottrain %>%
  arrange(X5)
predCovPowExp.plottrain

# Update the ggplot code
plottrainCovPowExp <- ggplot(predCovPowExp.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Power Exponential"
    #title = "Original Values and Predictions with Confidence Bounds for train set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovPowExp.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predCovPowExp.test$mean[predCovPowExp.test$mean < 0] <- 0
predCovPowExp.test$lower95[predCovPowExp.test$lower95 < 0] <- 0
#Prepare Dataframe
predCovPowExp.plottest <- cbind(predCovPowExp.test$mean, predCovPowExp.test$lower95, predCovPowExp.test$upper95, test$y, test$indexall)
predCovPowExp.plottest <- data.frame(predCovPowExp.plottest)
# Sort the dataframe by column X5
predCovPowExp.plottest <- predCovPowExp.plottest %>%
  arrange(X5)
predCovPowExp.plottest

# Update the ggplot code
plottestCovPowExp <- ggplot(predCovPowExp.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Power Exponential"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predCovPowExp.plottest)


#Plots myCovManPowerFixed
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predmyCovManPowerFixed.train$mean[predmyCovManPowerFixed.train$mean < 0] <- 0
predmyCovManPowerFixed.train$lower95[predmyCovManPowerFixed.train$lower95 < 0] <- 0
#Prepare Dataframe
predmyCovManPowerFixed.plottrain <- cbind(predmyCovManPowerFixed.train$mean, predmyCovManPowerFixed.train$lower95, predmyCovManPowerFixed.train$upper95, train$y, train$indexall)
predmyCovManPowerFixed.plottrain <- data.frame(predmyCovManPowerFixed.plottrain)
# Sort the dataframe by column X5
predmyCovManPowerFixed.plottrain <- predmyCovManPowerFixed.plottrain %>%
  arrange(X5)
predmyCovManPowerFixed.plottrain

# Update the ggplot code
plottrainmyCovManPowerFixed <- ggplot(predmyCovManPowerFixed.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Additive Power Exponential"
    #title = "Original Values and Predictions with Confidence Bounds for train set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predmyCovManPowerFixed.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predmyCovManPowerFixed.test$mean[predmyCovManPowerFixed.test$mean < 0] <- 0
predmyCovManPowerFixed.test$lower95[predmyCovManPowerFixed.test$lower95 < 0] <- 0
#Prepare Dataframe
predmyCovManPowerFixed.plottest <- cbind(predmyCovManPowerFixed.test$mean, predmyCovManPowerFixed.test$lower95, predmyCovManPowerFixed.test$upper95, test$y, test$indexall)
predmyCovManPowerFixed.plottest <- data.frame(predmyCovManPowerFixed.plottest)
# Sort the dataframe by column X5
predmyCovManPowerFixed.plottest <- predmyCovManPowerFixed.plottest %>%
  arrange(X5)
predmyCovManPowerFixed.plottest

# Update the ggplot code
plottestmyCovManPowerFixed <- ggplot(predmyCovManPowerFixed.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Additive Power Exponential"
    #title = "Original Values and Predictions with Prediction Bounds for test set (PB2 Mattern32)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "none",  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
dim(predmyCovManPowerFixed.plottest)

# Arrange the plots in a 2x3 grid
grid.arrange(plottrainCovM52, plottrainCovM32, plottrainCovExp, plottrainCovGauss, plottrainCovPowExp, plottestmyCovManPowerFixed, nrow = 2, ncol = 3)

# Arrange the plots in a 2x3 grid
grid.arrange(plottestCovM52, plottestCovM32, plottestCovExp, plottestCovGauss, plottestCovPowExp, plottestmyCovManPowerFixed, nrow = 2, ncol = 3)


################################################################################
#Plots myCovMan
#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)
# Create a plot to display the results for train
predmyCovMan.train$mean[predmyCovMan.train$mean < 0] <- 0
predmyCovMan.train$lower95[predmyCovMan.train$lower95 < 0] <- 0
#Prepare Dataframe
predmyCovMan.plottrain <- cbind(predmyCovMan.train$mean, predmyCovMan.train$lower95, predmyCovMan.train$upper95, train$y, train$indexall)
predmyCovMan.plottrain <- data.frame(predmyCovMan.plottrain)
# Sort the dataframe by column X5
predmyCovMan.plottrain <- predmyCovMan.plottrain %>%
  arrange(X5)
predmyCovMan.plottrain

# Update the ggplot code
ggplot(predmyCovMan.plottrain, aes(x = X5)) +
  geom_point(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Original Values and Predictions with Confidence Bounds for train set (PB2 AdditiveExp)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.6, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25),  # Increase legend title size
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25)  # Increase axis title size
  )

dim(predmyCovMan.plottrain)

# Create a plot to display the results for test
# Create a plot to display the results for test
predmyCovMan.test$mean[predmyCovMan.test$mean < 0] <- 0
predmyCovMan.test$lower95[predmyCovMan.test$lower95 < 0] <- 0
#Prepare Dataframe
predmyCovMan.plottest <- cbind(predmyCovMan.test$mean, predmyCovMan.test$lower95, predmyCovMan.test$upper95, test$y, test$indexall)
predmyCovMan.plottest <- data.frame(predmyCovMan.plottest)
# Sort the dataframe by column X5
predmyCovMan.plottest <- predmyCovMan.plottest %>%
  arrange(X5)
predmyCovMan.plottest

# Update the ggplot code
ggplot(predmyCovMan.plottest, aes(x = X5)) +
  geom_line(aes(y = X1, color = "Predicted Value"), size = 1) +  # Line for mean
  geom_ribbon(aes(ymin = X2, ymax = X3, fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = X4, color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "Original Values and Predictions with Confidence Bounds for test set (PB2 AdditiveExp)"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.6, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25),  # Increase legend title size
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25)  # Increase axis title size
  )

dim(predmyCovMan.plottest)


# This code can be repeated for other experiments (for the high dimensional case) by making a few changes such as the number of components (d) and as well as for Experiments I - IV