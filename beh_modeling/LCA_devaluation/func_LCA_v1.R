
# Function to simulate LCA process with asymmetric inhibition

# S1, S2: input values
# bias1: bias term to option 1
# bias2: bias term to option 2
# gamma: leakage rate
# beta1: inhibition from option 2 to option 1
# beta2: inhibition from option 1 to option 2
# sigma: noise std
# dt: time step size
# n_steps: time steps for accumulation process
# threshold: decision threshold
# eps: identity precision, affecting how well cues are assigned to two accumulators

func_LCA <- function(S1, S2, bias1, bias2, gamma, beta_1, beta_2, 
                     sigma, dt, threshold, eps) {
  # Initialize accumulator states for both options
  x1 <- 0 # sated odor option accumulator
  x2 <- 0 # non-sated odor option accumulator
  out = NA # initiate out as NA
  
  # S1 and S2 will be assigned to the correct accumulators with a probability of eps
  flip = TRUE
  u <- runif(1)
  if (u <= eps){
    flip = FALSE
  }

  if (flip) {
    S3 = S1
    S1 = S2
    S2 = S3
  }
  
  for (t in 1:n_steps) {
    input1 <- S1 + bias1  # bias towards option 1, informed by odor pleasantness
    input2 <- S2 + bias2  # bias towards option 2, informed by odor pleasantness
    
    # Update accumulators with leakage, asymmetric competition, and noise
    dx1 <- (input1 - gamma * x1 - beta_1 * x2) * dt + sigma * rnorm(1, 0, sqrt(dt))
    dx2 <- (input2 - gamma * x2 - beta_2 * x1) * dt + sigma * rnorm(1, 0, sqrt(dt))
    
    x1 <- x1 + dx1
    x2 <- x2 + dx2
    
    # Check if any accumulator reaches the threshold
    if (x1 >= threshold) {
      out = 1  # Choose option 1 (sated odor)
    } else if (x2 >= threshold) {
      out = 0  # Choose option 2 (non-sated odor)
    }
  }
  # Default to choosing the option with higher strength if no threshold is reached
  if(is.na(out)){
  out = ifelse(S1 > S2, 1, 0)
  }
  
  if (flip){ # but eventually we want to report beh associated with the input
    out = 1 - out
  }
  
  return(out)
  
}
