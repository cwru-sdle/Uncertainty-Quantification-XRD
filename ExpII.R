## Need to install Rcpp, RcppCNPy, fields packages
library(RcppCNPy)
library(fields)
library(caret)
library(dplyr)
library(tidyverse)
library(GPfit)
library(ggplot2)
library(gridExtra)

#SDLE PLOT STYLE
source("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification/SDLE_Style.R", echo=TRUE)

#Load data in
setwd("~/CSE_MSE_RXF131/cradle-members/mds3/aeo49/Git/Uncertainty-quantification")

datafile=read.csv("gpd_pb2.csv", header=T)
filenames=datafile[,2]
# Replace "rstor" with "vstor" in the paths of Column2
filenames <- sub("/mnt/rstor/", "/mnt/vstor/", filenames)
y=datafile[,3]
y <- data.frame(y)

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
## number of components to be kept
pr_n <- 5

# Plot eigenvalues
par(mfrow=c(1,2))
plot(variance_explained, type="b", xlab = "Number of PCs", xlim = c(1, 20), 
     ylab = "Proportion of Variance Explained", main = "Variance Explained by Each PC for Exp II")

# Add legend
legend("topright", legend=c("Variance Explained"), col="black", pch=1)

# Plot cumulative variance
plot(cumulative_variance, type="b", xlab = "Number of PCs", xlim = c(1, 20),  
     ylab = "Cumulative Proportion of Variance Explained", 
     main = "Cumulative Variance Explained for Exp II")

# Add horizontal line at y = 0.995, y = 0.995 and y = 0.95
abline(h = c(0.995, 0.99, 0.95), col = c("red", "blue","orange"), lty = c(2, 2, 2))

# Add legend for the horizontal lines
legend("bottomright", legend=c("99.5% Cumulative Variance Explained, n = 17", "99% Cumulative Variance Explained, n = 12", "95% Cumulative Variance Explained, n = 5"), col=c("red", "blue", "orange"), lty=2)

# Axis adjustments
axis(2, at = seq(0.5, 1, by = 0.1))

## plot the eigen vectors to show the dominant ones
for (i in 1:pr_n) { 
  par(mfrow=c(1,1))
  image.plot(matrix(pr$rotation[,i], nrow=resize_pxl), col = gray.colors(100))
}

# Calculate the proportion of variance explained by each PC
variance_explained <- pr$sdev^2 / sum(pr$sdev^2)

# Calculate the cumulative proportion of variance explained
cumulative_variance <- cumsum(variance_explained)

# Determine how many PCs to keep to explain 95% of the variance
num_pcs_to_keep <- which(cumulative_variance >= 0.95)[1]

# Print the results
cat("Number of PCs to keep for 95% variance explained:", num_pcs_to_keep, "\n")

## number of components to be kept
pr_n=5

## plot eigen values
par(mfrow=c(1,1))
plot(pr$sdev^2/sum(pr$sdev^2), type="b",xlab = "Number of PCs", xlim = c(1, 20), ylab = "Proportion of Variance Explained", main = "Variance Explained by Each PC")
plot(cumsum(pr$sdev^2)/sum(pr$sdev^2), type="b", xlab = "Number of PCs", xlim = c(1, 20), ylab = "Cumulative Proportion of Variance Explained", main = "Cumulative Variance Explained")
axis(2, at = seq(0.5, 1, by = 0.1))

## plot the eignen vectors to show the dominant ones
for (i in 1:5) {
  par(mfrow=c(1,1))
  image.plot(matrix(pr$rotation[,i], nrow=resize_pxl), col = gray.colors(100))
}

## data reconstruction from truncated PCA
data_c=pr$x[,1:pr_n]%*%t(pr$rotation[,1:pr_n])

##plot the first recontructed data and the original data
set.seed(1125031)

## Inputs for Guassian Process using PCA
X=pr$x[,1:pr_n] ##Reducing dimension from 512*512 to pr_n

# Rescale the values to the interval (0, 1)
X <- (X - min(X)) / (max(X) - min(X))
Xy <- cbind(X,y)
colnames(Xy) <- c("PC1", "PC2",  "PC3", "PC4", "PC5", "BetaFraction")
View(Xy)
summary(Xy)

# Calculate the 90th percentile of the 6th column
p96 = quantile(Xy[, 6], 0.96)

# Subset the data where values in the 6th column are above the 90th percentile
subset_data = subset(Xy, Xy[, 6] >= p96)

# Print to view index subset_data
print(subset_data)

###############
random_numbers <- sample(97:145, 49, replace = FALSE)

# Set up the plotting area outside the loop
par(mfrow = c(1, 2))

for (i in 1:49) {
  original_data <- matrix(data[random_numbers[i], ], nrow = resize_pxl)
  reconstructed_data <- matrix(data_c[random_numbers[i], ], nrow = resize_pxl)
  
  original_zlim <- range(original_data, na.rm = TRUE)
  reconstructed_zlim <- range(reconstructed_data, na.rm = TRUE)
  combined_zlim <- range(original_zlim, reconstructed_zlim)
  
  # Plot Images
  image.plot(original_data, main = "Original", col = gray.colors(12), zlim = combined_zlim)
  image.plot(reconstructed_data, main = "Reconstructed using 5 PCs", col = gray.colors(12), zlim = combined_zlim)
}

# Set a random seed for reproducibility
set.seed(2503)

# Split the Xy into training, testing, and validation sets
train_prop <- 0.65  # Proportion of Xy for training
validation_prop <- 0.15  # Proportion of Xy for validation
test_prop <- 0.20   # Proportion of Xy for testing

# Calculate the number of rows for each set
n <- nrow(Xy)
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
train <- Xy[sorttrain, ]
test <- Xy[sorttest, ]
validation <- Xy[sortvalidation, ]


# Print the number of rows in each set
cat("Number of rows in training set: ", nrow(train), "\n")
cat("Number of rows in testing set: ", nrow(test), "\n")
cat("Number of rows in validation set: ", nrow(validation), "\n")

# Record start time
start_time <- Sys.time()

#Fit GP_model
gp_model <- GPfit::GP_fit(train[, 1:5], train[, 6])

# Record end time
end_time <- Sys.time()

# Calculate total time
total_time <- end_time - start_time

# Print start, end, and total time
cat("Start time:", start_time, "\n")
cat("End time:", end_time, "\n")
cat("Total time:", total_time, "\n")

# Predictions for training data
predictions_train <- predict(gp_model, train[, 1:5])

# Calculate the standard error
train_standard_error <- sqrt(predictions_train$MSE)

# Calculate the degrees of freedom
train_df <- n_train - 1  # For a mean squared error (MSE), df is n - 1

# Calculate the critical value for a 95% confidence interval (two-tailed)
alpha <- 0.05
train_critical_value <- qt(1 - alpha / 2, train_df)

# Calculate the margin of error
train_margin_of_error <- train_critical_value * train_standard_error

# Calculate the lower and upper bounds of the confidence interval
train_lower_bound <- predictions_train$Y_hat - train_margin_of_error
train_upper_bound <- predictions_train$Y_hat + train_margin_of_error

# Create a data frame with the results
train_original <- train[, 6]
predictions_train <- cbind(predictions_train$complete_data, train_lower_bound, train_upper_bound, train_original)

# Print or view the confidence interval data
print(predictions_train)


# Predictions for testing data 
predictions_test <- predict(gp_model, test[, 1:5])

# Calculate the standard error
test_standard_error <- sqrt(predictions_test$MSE)

# Calculate the degrees of freedom
test_df <- n_train - 1  # For a mean squared error (MSE), df is n - 1

# Calculate the critical value for a 95% confidence interval (two-tailed)
alpha <- 0.05
test_critical_value <- qt(1 - alpha / 2, test_df)

# Calculate the margin of error
test_margin_of_error <- test_critical_value * test_standard_error

# Calculate the lower and upper bounds of the confidence interval
test_lower_bound <- predictions_test$Y_hat - test_margin_of_error
test_upper_bound <- predictions_test$Y_hat + test_margin_of_error

#Turn any value less than Zero to Zero
test_lower_bound[test_lower_bound < 0] <- 0

# Create a data frame with the results
test_original <- test[, 6]
predictions_test <- cbind(predictions_test$complete_data, test_lower_bound, test_upper_bound, test_original)

#Turn any value less than Zero to Zero
predictions_test[, 6][predictions_test[, 6]<0]<-0

# Print or view the confidence interval data
print(predictions_test)

# Create a plot to display the results for train
predictions_train <- data.frame(predictions_train)
trainplot <- ggplot(predictions_train, aes(x = sorttrain)) +
  geom_line(aes(y = predictions_train[, 6], color = "Predicted Value"), size = 1) +  # Line for y_hat
  geom_ribbon(aes(ymin = predictions_train[, 8], ymax = predictions_train[, 9], fill = "Confidence Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = predictions_train[, 10], color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "(a) Training Set"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(1, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )

# Calculate Mean Squared Prediction Error (MSPE) for train set
mspe_train <- mean((predictions_train[, 10] - predictions_train[, 6])^2)

# Print or use the value of MSPE for test set
print(mspe_train)

# Create a plot to display the results for test
predictions_test <- data.frame(predictions_test)
testplot <- ggplot(predictions_test, aes(x = sorttest)) +
  geom_line(aes(y = predictions_test[, 6], color = "Predicted Value"), size = 1) +  # Line for y_hat
  geom_ribbon(aes(ymin = predictions_test[, 8], ymax = predictions_test[, 9], fill = "Prediction Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = predictions_test[, 10], color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "Observation",
    y = "Beta Fraction",
    title = "(c) Test Set"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(1, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )
# Calculate Mean Squared Prediction Error (MSPE) for test set
mspe_test <- mean((predictions_test[, 10] - predictions_test[, 6])^2)

# Print or use the value of MSPE for test set
print(mspe_test)


trainintplot <- ggplot(predictions_train, aes(x = sorttrain)) +
  geom_line(aes(y = test_critical_value*sqrt(predictions_train[, 7]))) +  # Line for y_hat
  labs(
    x = "Observation",
    y = "Length of Confidence Interval",
    title = "(b) Training Set Error Length"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(1, 0.85),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )



#plot(sorttrain, test_critical_value*sqrt(predictions_train[, 7]), type = "l", ylab = "Lengh of Confidence Interval", xlab = "Observation", main = "Lengh of Confidence Interval for Training Set")
# lines(sorttrain, predictions_train[, 8], type = "l", col = "black")
# plot(predictions_train[, 6] - predictions_train[, 8], type = "l")
# plot(predictions_test[, 6] - predictions_t[, 8], type = "l")

# Create a plot to display the results for test
predictions_test <- data.frame(predictions_test)
test1pcplot <- ggplot(predictions_test, aes(x = predictions_test[, 1])) +
  geom_line(aes(y = predictions_test[, 6], color = "Predicted Value"), size = 1) +  # Line for y_hat
  geom_ribbon(aes(ymin = predictions_test[, 8], ymax = predictions_test[, 9], fill = "Prediction Interval"), alpha = 0.5) +  # Grey color bands for confidence interval
  geom_point(aes(y = predictions_test[, 10], color = "Original Value"), size = 2) +  # Points for y
  labs(
    x = "PC1",
    y = "Beta Fraction",
    title = "(d) Test Set [1 PC]"
  ) +
  SDLE_theme() +
  scale_color_manual(values = c("Predicted Value" = "deeppink3", "Original Value" = "steelblue3")) +
  scale_fill_manual(values = "grey") +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  theme(
    legend.position = c(1, 0.8),  # Position legend at the top
    legend.justification = "right",  # Justify legend to the right
    legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12),  # Increase legend title size
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),  # Adjust right margin to accommodate legend
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 14)  # Increase axis title size
  )

# Arrange the plots in a 2x3 grid
grid.arrange(trainplot, trainintplot, testplot, test1pcplot, nrow = 2, ncol = 2)
