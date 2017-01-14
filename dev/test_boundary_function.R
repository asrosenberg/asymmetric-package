library(asymmetric)
library(npbr)
library(quantreg)
library(frontier)
library(Rglpk)
source("dev/estimate_boundary.R")
source("R/functions.R")
x <- runif(75, min = 0, max = 1)
y <- 1/2 * x^(1/2) * runif(75, min = 0, max = 2)
plot(x, y)
test_data <- data.frame(input_black = x, output_blue = y)
estimate_boundary(test_data, input = test_data$input_black,
  output = test_data$output_blue, method = "QR", AOC = FALSE)
estimate_boundary(test_data, input = test_data$input_black,
  output = test_data$output_blue, method = "SFA", AOC = FALSE)
estimate_boundary(test_data, input = test_data$input_black,
  output = test_data$output_blue, method = "Kernel", AOC = FALSE)
estimate_boundary(test_data, input = test_data$input_black,
  output = test_data$output_blue, method = "Polynomial", AOC = FALSE)
# try to get it to work with sufficient data
# All the techniques require changing with sufficient data
# So, the best thing to do is add a call for sufficient condition
x <- (runif(75, min = 0, max = 1))
y <- (1/2 * x^(1/2) * runif(75, min = 0, max = 2))
sufx <- 1 - x
sufy <- 1 - y
suf_Dat <- data.frame(sufx, sufy)
plot(sufx, sufy)
test <- estimate_boundary(suf_Dat, input = 1 - sufx,
  output = 1 - sufy, method = "QR", AOC = FALSE)
points(x = sort(sufx), y = sort(1 - test), type = "l")
# need to flip both x and y
# then return 1 - boundary

plot(sufx, sufy)
test <- estimate_boundary(suf_Dat, input = sufx,
  output = sufy, method = "Polynomial", AOC = FALSE, sufficient = TRUE)
points(x = sort(sufx), y = sort(test), type = "l")




