
# Function to simulate LCA process with asymmetric inhibition

func_LCA <- function(S1, S2, bias, gamma, beta_1, beta_2, sigma, dt, threshold) {
  # Initialize accumulator states for both options
  x1 <- 0
  x2 <- 0
  
  for (t in 1:n_steps) {
    # Calculate inputs (evidence) for each option
    input1 <- S1
    input2 <- S2 + bias  # Bias is added to S2
    
    # Update accumulators with leakage, asymmetric competition, and noise
    dx1 <- (input1 - gamma * x1 - beta_1 * x2) * dt + sigma * rnorm(1, 0, sqrt(dt))
    dx2 <- (input2 - gamma * x2 - beta_2 * x1) * dt + sigma * rnorm(1, 0, sqrt(dt))
    
    x1 <- x1 + dx1
    x2 <- x2 + dx2
    
    # Check if any accumulator reaches the threshold
    if (x1 >= threshold) {
      return(1)  # Choose option 1 (sated odor)
    } else if (x2 >= threshold) {
      return(0)  # Choose option 2 (non-sated odor)
    }
  }
  # Default to choosing the option with higher strength if no threshold is reached
  return(ifelse(S1 > S2, 1, 0))
}
