# Parameters
num_simulations <- 20000  # Number of simulations
sample_size <- 30        # Size of each sample
population_mean <- 0     # Mean under the null hypothesis

# Initialize a vector to store p-values
p_values <- numeric(num_simulations)

# Perform the simulations
for (i in 1:num_simulations) {
  # Generate a random sample from a normal distribution with mean 0
  sample <- rnorm(sample_size, mean = population_mean, sd = 1)
  
  # Perform a one-sample t-test
  test_result <- t.test(sample, mu = population_mean)
  
  # Store the p-value
  p_values[i] <- test_result$p.value
}

# Plot the distribution of p-values
hist(p_values, breaks = 30, main = "Distribution of P-values Under the Null Hypothesis",
     xlab = "P-value", col = "lightblue", border = "black")
