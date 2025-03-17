
# this simulation wants to try to interpret selective satiation
# using different mechanisms of LCA
# why do people prefer non-sated odor after eating the meal?
# attempt 1: asymmetric inhibition (higher from non-sated to sated)
# attempt 2: bias (positive bias towards non-sated odor)

rm(list=ls()) 
library(openxlsx)
library(plyr)
library(tidyverse)
source('../Scripts_R/Setup.R')
source('func_LCA_v0.R')

set.seed(123)  # For reproducibility

# set up
# option1: sated odor
# option2: non-sated odor

# Define parameters for the LCA model
n_steps <- 1000        # Time steps for accumulation process
gamma <- 0.1           # Leakage rate
threshold <- 1.0       # Decision threshold
sigma <- 0.4          # Noise standard deviation
dt <- 0.01             # Time step size
S_grid <- seq(0.1, 0.9, length.out = 5)  # Grid of values for S1 and S2
n_simulations <- 50  # Number of repetitions per S1 and S2 combination


#####################################
# simulate with lateral inhibition
#####################################

betas = expand.grid(beta_1 = c(0.1,0.5,0.9), beta2 = c(0.1,0.5,0.9))
choice_outcomes_df_with_betas = NULL

for(b in 1:nrow(betas)){

beta_1 <- betas[b,1] # Inhibition from option 2 to option 1 (asymmetric)
beta_2 <- betas[b,2] # Inhibition from option 1 to option 2 (asymmetric)
bias <- 0            # Bias parameter (set to zero)
eps = 0.8            # precision of identity knowledge

# Initialize matrices to store choice outcomes and accuracy
choice_outcomes <- matrix(0, nrow = length(S_grid), ncol = length(S_grid))

# Simulate choices for varying levels of S1 and S2, with multiple repetitions
for (i in 1:length(S_grid)) {
  for (j in 1:length(S_grid)) {
    print(paste('i=',i,'j=',j))
    S1 <- S_grid[i]
    S2 <- S_grid[j]
    total_choices <- 0
    for (sim in 1:n_simulations) {
      n_trials <- 100  # Number of trials per simulation
      choices <- numeric(n_trials)
      for (trial in 1:n_trials) {
        choice <- func_LCA(S1, S2, bias, gamma, beta_1, beta_2, sigma, dt, 
                           threshold)
        choices[trial] <- choice
      }
      # Accumulate choices and correct choices
      total_choices <- total_choices + mean(choices)  # Average of choices (proportion choosing option 1)
    }
    # Calculate the average across simulations
    choice_outcomes[i, j] <- total_choices / n_simulations
    
  }
}

# Convert matrices to data frames for plotting
choice_outcomes_df <- expand.grid(S1 = S_grid, S2 = S_grid)
choice_outcomes_df$ChoiceOutcome <- as.vector(choice_outcomes)

choice_outcomes_df_with_betas[[b]] = choice_outcomes_df

}



xlab = 'Learned option 1 (SA) value'
ylab = 'Learned option 2 (NS) value'

p = NULL
for(b in 1:nrow(betas)){
  p[[b]] = ggplot(choice_outcomes_df_with_betas[[b]], 
                  aes(x = S1, y = S2, fill = ChoiceOutcome)) +
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
    scale_fill_gradient2(low = "red", high = "green", 
                         mid = 'white', midpoint = 0.5,
                         name = "P (sated odor)") +
    labs(title = paste('beta1 =',betas[b,1],'beta2 =',betas[b,2]), 
         x = xlab, y = ylab)
  
}

pdf('Simu_LCA_with_lateral_inhibition.pdf',12,12)
ggarrange(p[[1]],p[[2]],p[[3]],
          p[[4]],p[[5]],p[[6]],
          p[[7]],p[[8]],p[[9]],
          common.legend = T,legend = 'right',
          nrow = 3,ncol = 3)
dev.off()

saveRDS(choice_outcomes_df_with_betas,file = 'Simu_LCA_strength_inhibition.rds')


#####################################
# simulate with bias
#####################################

beta_1 <- 0.1          # Inhibition from option 2 to option 1 (asymmetric)
beta_2 <- 0.1          # Inhibition from option 1 to option 2 (asymmetric)
bias <- c(0.1,0.5,0.9)              # Bias parameter (set to zero)

choice_outcomes_df_with_bias = NULL

for(b in 1:length(bias)){
# Initialize matrices to store choice outcomes and accuracy
choice_outcomes <- matrix(0, nrow = length(S_grid), ncol = length(S_grid))

# Simulate choices for varying levels of S1 and S2, with multiple repetitions
for (i in 1:length(S_grid)) {
  for (j in 1:length(S_grid)) {
    print(paste('i=',i,'j=',j))
    S1 <- S_grid[i]
    S2 <- S_grid[j]
    
    # Initialize variables to accumulate results across simulations
    total_choices <- 0
    
    # Repeat simulations for more stable results
    for (sim in 1:n_simulations) {
      n_trials <- 100  # Number of trials per simulation
      choices <- numeric(n_trials)
      
      for (trial in 1:n_trials) {
        choice <- func_LCA(S1, S2, bias[b], gamma, beta_1, beta_2, sigma, dt, threshold)
        choices[trial] <- choice
      }
      
      # Accumulate choices and correct choices
      total_choices <- total_choices + mean(choices)  # Average of choices (proportion choosing option 1)
        }
    
    # Calculate the average across simulations
    choice_outcomes[i, j] <- total_choices / n_simulations
    
  }
}

# Convert matrices to data frames for plotting
choice_outcomes_df <- expand.grid(S1 = S_grid, S2 = S_grid)
choice_outcomes_df$ChoiceOutcome <- as.vector(choice_outcomes)

choice_outcomes_df_with_bias[[b]] = choice_outcomes_df

}


p = NULL
for(b in 1:length(bias)){
  p[[b]] = ggplot(choice_outcomes_df_with_bias[[b]], 
                  aes(x = S1, y = S2, fill = ChoiceOutcome)) +
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + 
    scale_fill_gradient2(low = "red", high = "green", 
                         mid = 'white', midpoint = 0.5,
                         name = "P (sated odor)") +
    labs(title = paste('bias =',bias[b]), 
         x = xlab, y = ylab)
  
}

saveRDS(choice_outcomes_df_with_bias,file = 'Simu_LCA_strength_bias.rds')

pdf('Simu_LCA_with_bias.pdf',12,4)
ggarrange(p[[1]],p[[2]],p[[3]],
            common.legend = T,nrow = 1,ncol = 3,
            legend = 'right')
dev.off()

# a quick comment on this simulation
# not surprisingly, both lateral inhibition and bias terms can cause a shift
# of choice towards non-sated odors
# but interestingly, for choices on the diagonal line, lateral inhibition
# influences options with higher values more, while bias influences options 
# with lower values more






