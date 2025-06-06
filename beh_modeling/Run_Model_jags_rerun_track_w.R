
# modified from 'Run_Model_jags.R'
# but want to keep track of the learned strength w in the end
# but w has too large size so increase 'thin' parameter from 1 to 10

# also realize that learning rate has to vary with cue pair
# otherwise the modeled learned value of each cue pair will certainly be
# the same

rm(list=ls()) 
library(openxlsx)
library(rjags)
library(R2jags)
library(plyr)
library(tidyverse)
source('../Scripts_R/Setup.R')

# load data to use
load(file = '../ProcessedData/Conditioning.RData') # this contains some NA trials
use_dat = conditioning_dat  

# prepare model inputs as vectors
sub = use_dat$SubID
sub = mapvalues(sub,from = unique(sub),to = 1:length(unique(sub))) # recode sub into consective numbers
nsubs = length(unique(sub))
sess = use_dat$Sess
nsess = max(sess)
cond = 3 - as.numeric(use_dat$Cond_day1) # sham-1; cTBS-2
ncond = length(unique(cond))

ntrials = 120 # 24 trials/run * 5runs for each sess
# here we simply collapse across 5 runs as learning should be continous
  
CuePair = TMS = Resp = array(NA,c(nsubs,nsess,ntrials))
for(j in 1:nsubs){
  for(c in 1:nsess){
      CuePair[j,c,] = use_dat$CuePair[sub==j & sess==c] # vals: 1,2,...,12
      Resp[j,c,] = use_dat$OdorChosen[sub==j & sess==c] # vals: 0,1
      TMS[j,c,] = cond[sub==j & sess==c] # vals: 1,2
  }
}

# w track the change of strength across trials
ncuepair = 12
w = array(NA,c(nsubs,nsess,(ntrials+1),ncuepair))
w[,,1,] = 1/2

dlist <- list("nsubs","nsess","ntrials","ncuepair",
              "CuePair","Resp","w")

################ set up MCMC ##############

n.chains = 3
n.iter = 5e3
n.burnin = 2e3
n.thin = 10

##########################

# 3 models to test here
# 'bySess': learning rates are varied by session order
# 'same': constant learning rates for each subj
# 'byTMS': learning rates are varied by TMS condition: sham, cTBS

##########################

opt = 'bySess'
model.file = paste0('LearningModel_',opt,'_rerun_track_w.txt')
myinits = rep(list(list(mu = rep(0.5,nsess), 
                        kappa = 1,
                        alpha = array(runif(nsubs * nsess * ncuepair, 
                        min = 0.1, max = 0.9), dim = c(nsubs,nsess,ncuepair)))),
               n.chains)
parameters <- c("mu","alpha","kappa","w")
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
end_time - start_time
saveRDS(samples,file = paste0('./JagsPosteriors/',opt,'rerun_track_w.rds'))



