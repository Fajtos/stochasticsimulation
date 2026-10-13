# Load the Adaptive Rejection Sampling package
library(ars)

# Defining the parameters of the GHD
mu <- 0
delta <- 1
alpha <- 2
beta <- 0.5

# Define the log density (h) and its derivative (h_prime)
h_func <- function(x) {
  -alpha * sqrt(delta^2 + (x - mu)^2) + beta * (x - mu)
}

h_prime_func <- function(x) {
  -alpha * (x - mu) / sqrt(delta^2 + (x - mu)^2) + beta
}

# Initial points for the envelope
init_points <- c(-2, 0, 2)

# Running the Adaptive Rejection Sampler
set.seed(123)
n <- 1000
sample_gh <- ars(n = n, f = h_func, fprima = h_prime_func, x = init_points)

# Printing the first few values
head(sample_gh)

# Visualization of Results
library(ghyp)

hist(sample_gh, breaks = 40, probability = TRUE, 
     main = "Adaptive Rejection Method: Generalized Hyperbolic Distribution", 
     xlab = "x", ylab = "Density", col = "lightblue")

# Create GHD object using the correct alpha/delta constructor
ghd_obj <- ghyp.ad(lambda = 1, alpha = 2, beta = 0.5, delta = 1, mu = 0)

# Overlay theoretical density
curve(dghyp(x, object = ghd_obj), add = TRUE, col = "red", lwd = 2)

legend("topright", legend = c("Simulated Sample", "Theoretical Density"), 
       col = c("lightblue", "red"), lwd = c(1, 2), fill = c("lightblue", NA))

# Visualizing the Tangent Envelope
x_seq <- seq(-4, 6, length.out = 500)
h_vals <- -2 * sqrt(1 + (x_seq - 0)^2) + 0.5 * (x_seq - 0)

# Ploting the log density
plot(x_seq, h_vals, type = "l", lwd = 2, col = "black",
     main = "Log Density and Tangent Envelope",
     xlab = "x", ylab = "h(x)", ylim = c(-15, 5))

# Drawing tangent lines at specific points (x = -2, 0, 2)
tangent_points <- c(-2, 0, 2)
for (tp in tangent_points) {
  h_tp <- -2 * sqrt(1 + (tp - 0)^2) + 0.5 * (tp - 0)
  h_prime_tp <- -2 * (tp - 0) / sqrt(1 + (tp - 0)^2) + 0.5
  abline(a = h_tp - h_prime_tp * tp, b = h_prime_tp, col = "red", lty = 2)
  points(tp, h_tp, pch = 19, col = "blue")
}
legend("topright", legend = c("Log Density h(x)", "Tangent Envelopes", "Tangent Points"), 
       col = c("black", "red", "blue"), lty = c(1, 2, NA), pch = c(NA, NA, 19))


# Monte Carlo Estimates from the simulated sample
mc_mean <- mean(sample_gh)
mc_var  <- var(sample_gh)

# Exact Theoretical values using the ghyp object created in Part 1
true_mean <- mean(ghd_obj)
true_var  <- as.numeric(vcov(ghd_obj))

# Display the comparison table
results <- data.frame(
  Metric = c("Mean", "Variance"),
  Monte_Carlo = c(round(mc_mean, 4), round(mc_var, 4)),
  Theoretical = c(round(true_mean, 4), round(true_var, 4)),
  Absolute_Difference = c(round(abs(mc_mean - true_mean), 4), round(abs(mc_var - true_var), 4))
)

print(results)
