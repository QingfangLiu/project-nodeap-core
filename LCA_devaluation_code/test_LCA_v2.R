


rm(list=ls()) 
source('func_LCA_v2.R')

# set up
# option1: sated odor
# option2: non-sated odor

# Define parameters for the LCA model
n_steps <- 2000        # Time steps for accumulation process
gamma <- 0.1           # Leakage rate
threshold <- 4       # Decision threshold
sigma <- 0.4           # Noise standard deviation
dt <- 0.01             # Time step size
t0 = 2                 # non-decision time
n_simulations <- 1e2  # Number of repetitions per S1 and S2 combination

S1 <- 0.8
S2 <- 0.8

eps <- 1      # precision of identity knowledge
beta_1 <- 0.3  # Inhibition from option 2 to option 1 
beta_2 <- 0.1             # Inhibition from option 1 to option 2
bias1 <- 0                # Bias to option 1 (set to zero)
bias2 <- 0                # Bias to option 2 (set to zero)

simu1 = func_LCA(S1, S2, bias1, bias2, gamma, beta_1, beta_2, 
                sigma, dt, threshold, eps, t0)

beta_1 <- 0.1  # Inhibition from option 2 to option 1 
beta_2 <- 0.1             # Inhibition from option 1 to option 2
bias1 <- 0                # Bias to option 1 (set to zero)
bias2 <- 0.3                # Bias to option 2 (set to zero)

simu2 = func_LCA(S1, S2, bias1, bias2, gamma, beta_1, beta_2, 
                sigma, dt, threshold, eps, t0)

# keep some good simulation examples for better visualization
save.image(file = 'test_LCA_v2.RData')

load('test_LCA_v2.RData')
par(mfrow=c(1,2))
matplot(simu1$state,xlim = c(0,1000),type = 'l')
matplot(simu2$state,xlim = c(0,1000),type = 'l')

source('../Scripts_R/Setup.R')

colnames(simu1$state) <- c("Sated", "Non-sated")
time <- 1:nrow(simu1$state)  
mat_df <- data.frame(Time = time, simu1$state)
mat_long <- pivot_longer(mat_df, cols = -Time, 
                         names_to = "Option", values_to = "Value")
mat_long$Option = factor(mat_long$Option)
levels(mat_long$Option) = c('Non-sated','Sated')
p1=ggplot(mat_long, aes(x = Time, y = Value, color = Option)) +
  geom_line(aes(lty=Option),size=1.2) +
  labs(title = "Lateral inhibition",
       x = "Time step",
       y = "Evidence",
       color = "Option") +
  coord_cartesian(xlim = c(0,700)) +
  scale_linetype_manual(values = c('Non-sated'="solid", 'Sated'="dotted")) +
  scale_color_manual(values = c('Sated'='gray','Non-sated'='#832a78')) +
  annotate("segment",x = 500, xend = 500, y = 3.2, yend = 1.1,  
    arrow = arrow(type = "closed", length = unit(0.2, "cm")),
    color = "black", size = 1) +
  common + theme(legend.position = 'none')

colnames(simu2$state) <- c("Sated", "Non-sated")
mat_df <- data.frame(Time = time, simu2$state)
mat_long <- pivot_longer(mat_df, cols = -Time, 
                           names_to = "Option", values_to = "Value")
mat_long$Option = factor(mat_long$Option)
levels(mat_long$Option) = c('Non-sated','Sated')
p2=ggplot(mat_long, aes(x = Time, y = Value, color = Option)) +
    geom_line(aes(lty=Option),size=1.3) +
    labs(title = "Drift rate",
         x = "Time step",
         y = "Evidence",
         color = "Option") +
    coord_cartesian(xlim = c(0,700)) +
    scale_linetype_manual(values = c('Non-sated'="solid", 'Sated'="dotted")) +
    scale_color_manual(values = c('Sated'='gray','Non-sated'='#832a78')) +
    annotate("segment",x = 0, xend = 550, y = 0, yend = 4,  
           arrow = arrow(type = "closed", length = unit(0.2, "cm")),
           color = "#832a78", size = 1) +
    annotate("segment",x = 0, xend = 550, y = 0, yend = 3,  
           arrow = arrow(type = "closed", length = unit(0.2, "cm")),
           color = "gray", size = 1) +
    common +
    theme(legend.position = c(0.75, 0.2))

pdf('LCA_illustrations.pdf',8,4)
ggarrange(p2,p1,nrow = 1)
dev.off()



