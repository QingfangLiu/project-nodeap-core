
########################################################
# Purpose
# - Compute posterior means for:
#   (1) Learning rates α[j, c, p]  → saved as a 3D array (.rds)
#   (2) End-of-learning latent values w[j, c, p] at i = ntrials+1
#       → saved as a long CSV for downstream analyses/plots
# - Produce a PDF figure summarizing w by session (faceted by subject).
#
# Inputs
# - Setup.R: defines paths (e.g., processed_dir, beh_model_dir), colors (use.col.sess).
# - beh_model_dir/results/posterior_learning_model.rds  (R2jags “bugs” object)
#   • sims.list$alpha  [draw, sub, sess, cuepair]
#   • sims.list$w      [draw, sub, sess, trial_index, cuepair]
# - processed_dir/Conditioning.RData (for SubID/Sess indexing)
#
# Outputs
# - beh_model_dir/results/alpha_mean_array.rds                  (posterior mean α, 3D array)
# - beh_model_dir/results/df_w_per_sub_sess_cuepair.csv         (posterior mean w at i=121, long)
# - beh_model_dir/results/w_end_learning_by_sub_sess.pdf        (box+jitter by session, faceted by subject)
#
# Notes
# - “i = 121” corresponds to ntrials = 120 (model stores w at 1..ntrials+1).
# - This script uses posterior MEANS (not MAP).
########################################################


rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# Input posterior file
post_rds_name <- "posterior_learning_model.rds"
post_rds_path <- file.path(beh_model_dir, "results", post_rds_name)
fig_path <- file.path(beh_model_dir, "results", "w_end_learning_by_sub_sess.pdf")

# Output files
alpha_rds_path <- file.path(beh_model_dir, "results", "alpha_mean_array.rds")
w_csv_path     <- file.path(beh_model_dir, "results", "df_w_per_sub_sess_cuepair.csv")

# 1) Load posterior samples
samples <- readRDS(file = post_rds_path)

# Expected components in R2jags output
alpha <- samples$BUGSoutput$sims.list$alpha  # [draw, sub, sess, cuepair]
w     <- samples$BUGSoutput$sims.list$w      # [draw, sub, sess, trial_index, cuepair]

# 2) Summarize alpha (posterior mean) and save

# alpha: mean over MCMC draws → [sub, sess, cuepair]

alpha_mean <- apply(alpha, c(2, 3, 4), mean)
saveRDS(alpha_mean, alpha_rds_path)

# 3) Summarize end-of-learning w (posterior mean at i = 121) and save

# End-of-learning index = 121 = (ntrials + 1) with ntrials = 120
use_w      <- w[ , , , 121, ]                 # [draw, sub, sess, cuepair]
use_w_mean <- apply(use_w, c(2, 3, 4), mean)  # → [sub, sess, cuepair]

# Bring in indexing (SubID, Sess, ncuepair) from conditioning data
load(file = file.path(processed_dir, "Conditioning.RData")) 

use_dat = conditioning_dat  
sub = use_dat$SubID
nsubs = length(unique(sub))
sess = use_dat$Sess
nsess = max(sess)
ncuepair = 12

df_w = structure(use_w_mean,
                      .Dim = c(nsubs, nsess, ncuepair), 
                      .Dimnames = structure(list(SubID = unique(as.vector(sub)),
                                                 Sess = 1:nsess,
                                                 CuePair = 1:ncuepair
                                                 ), 
                                            .Names = c("SubID","Sess","CuePair")))
df_w = adply(df_w, c(1,2,3))
names(df_w) = c(names(df_w)[1:3],'w')
df_w = df_w %>% arrange(SubID,Sess)
write.csv(df_w, row.names = FALSE, file = w_csv_path)

# 4) Quick visualization 

pdf(fig_path, width = 8, height = 6)  # open PDF device

ggplot(df_w,aes(x=Sess,y=w,color=Sess,fill=Sess)) +
  geom_boxplot(alpha=0.4) +
  geom_jitter(alpha=0.5) +
    facet_wrap(~SubID) +
  scale_color_manual(values = use.col.sess) +
  scale_fill_manual(values = use.col.sess) +
  labs(y='Learned stimulus value',x='Session')

dev.off()  # close PDF device
