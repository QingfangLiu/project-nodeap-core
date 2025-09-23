
########################################################
# Estimate Hierarchical Learning Model (JAGS)
#
# Purpose
# Fit the project’s JAGS learning model and save posterior
# samples (including latent trajectories w) as an .rds file.
#
# Inputs (from Setup.R / data files)
# - Setup.R: defines `project_folder`, paths (e.g., `processed_dir`,
#   `beh_model_dir`), common libs/options.
# - processed_dir/Conditioning.RData:
#     * `conditioning_dat` with columns:
#       SubID, Sess, Cond_day1, CuePair (1..12), OdorChosen (0/1)
#
# Model file
# - beh_model_dir/model_code/learning_model.txt
#
# Key data passed to JAGS
# - Integers: nsubs, nsess, ntrials (=120), ncuepair (=12)
# - Arrays: CuePair[nsubs, nsess, ntrials], Resp[...], w[...,(ntrials+1),...]
#
# Outputs
# - results/posterior_learning_model.rds  (R2jags 'bugs' object)
# - results/posterior_design_arrays.rds         (CuePair, run_idx, meta_ss, etc.)
#
# Notes
# - w is tracked across trials; to limit object size, MCMC thinning
#   is increased (n.thin = 10).
# - Learning rate varies by cue pair; collapsing would force identical
#   learned values across pairs.
# - This run may take a while.
########################################################

rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

library(rjags)
library(R2jags)

# 1) Load behavioral data
load(file = file.path(processed_dir,'Conditioning.RData')) # this contains some NA trials
use_dat = conditioning_dat  

# 2) Prepare indices and core arrays for JAGS
sub = use_dat$SubID
sub = mapvalues(sub,from = unique(sub),to = 1:length(unique(sub))) # recode sub into consective numbers
nsubs = length(unique(sub))
sess = use_dat$Sess
nsess = max(sess)
cond = 3 - as.numeric(use_dat$Cond_day1) # sham-1; cTBS-2
ncond = length(unique(cond))

ntrials = 120  # 24 trials/run × 5 runs per session; treat as continuous learning
ncuepair = 12
  
# Allocate trial-wise arrays
CuePair = TMS = Resp = array(NA,c(nsubs,nsess,ntrials))
for(j in 1:nsubs){
  for(c in 1:nsess){
      CuePair[j,c,] = use_dat$CuePair[sub==j & sess==c] # vals: 1,2,...,12
      Resp[j,c,] = use_dat$OdorChosen[sub==j & sess==c] # vals: 0,1
      TMS[j,c,] = cond[sub==j & sess==c] # vals: 1,2
  }
}

# --- NEW: build per-session run index and a small metadata table ---
ntrials_per_run <- 24L
nruns           <- 5L
stopifnot(ntrials == ntrials_per_run * nruns)

# run_idx[sub, sess, trial] = 1..5 for trials 1..120 in correct order
run_idx <- array(NA_integer_, c(nsubs, nsess, ntrials))
run_templ <- rep(seq_len(nruns), each = ntrials_per_run)
for (j in seq_len(nsubs)) for (c in seq_len(nsess)) run_idx[j, c, ] <- run_templ

# (optional) per-(SubID,Sess) metadata you’ll want for plotting
meta_ss <- unique(
  data.frame(
    SubID     = use_dat$SubID,
    Sess      = as.integer(use_dat$Sess),
    Cond_day1 = use_dat$Cond_day1,
    StimLoc   = use_dat$StimLoc,
    stringsAsFactors = FALSE
  )
)

# 3) Initialize latent state w (trial 1 prior) and build data list

w = array(NA,c(nsubs,nsess,(ntrials+1),ncuepair))
w[,,1,] = 1/2 # initial value at trial 1
dlist <- list("nsubs","nsess","ntrials","ncuepair",
              "CuePair","Resp","w")

# 4) MCMC configuration

n.chains = 3
n.iter = 5e3
n.burnin = 2e3
n.thin = 10

# 5) Model file path

model.file = file.path(beh_model_dir,'model_code','learning_model.txt')

# 6) Initial values

myinits = rep(list(list(mu = rep(0.5,nsess), 
                        k = rep(1,nsubs),
                        a = array(runif(nsubs * nsess, 
                                        min = 0.1, max = 0.9), dim = c(nsubs,nsess)),
                        kappa = 1,
                        alpha = array(runif(nsubs * nsess * ncuepair, 
                        min = 0.1, max = 0.9), dim = c(nsubs,nsess,ncuepair)))),
               n.chains)

# Parameters to save (includes latent w)
parameters <- c("mu","k","a","kappa","alpha","w")

# 7) Run JAGS (this may take a while)

start_time <- proc.time()
samples <- jags(data = dlist,
                inits = myinits,
                parameters.to.save = parameters,
                model.file = model.file,
                n.chains = n.chains, 
                n.iter = n.iter,
                n.burnin = n.burnin, 
                n.thin = n.thin)
end_time <- proc.time()
print(end_time - start_time)

# 8) Save posterior object

out_filepath = file.path(beh_model_dir,'results','posterior_learning_model.rds')
saveRDS(samples, file = out_filepath)

# --- save design sidecar for plotting w without guessing structure ---
design_path <- file.path(beh_model_dir, "results", "posterior_design_arrays.rds")
saveRDS(
  list(
    CuePair          = CuePair,               # [sub × sess × trial]
    run_idx          = run_idx,               # [sub × sess × trial] (1..nruns)
    SubID_levels     = unique(use_dat$SubID),
    ntrials_per_run  = ntrials_per_run,
    nruns            = nruns,
    meta_ss          = meta_ss               # SubID, Sess, Cond_day1, StimLoc
  ),
  file = design_path
)
cat("Saved design arrays to:\n", design_path, "\n")


