
# assume a Q-learning process
# learn the value of each pair of cue
# b/c cues are always learned in fixed pairs
# savory & sweet contain both value 1
# non-odor contains value 0
# each trial, updating the value of chosen cue
# choice rule: 
  # a simple soft-max rule to get started?
  # can also incorporate RT?

load(file = '../ProcessedData/Conditioning.RData') # load processed conditioning data
use.data = subset(conditioning_dat,SubID=='NODEAP_06' & Sess==1)
ntrials = nrow(use.data)
choice = numeric(ntrials)
ncues = 24

part.data = subset(use.data,Run==1)
w = matrix(NA,nrow = ncues,ncol = 3+ntrials)
# col1: cue ID
# col2: expectation of odor of choosing this cue

# pre-set cue ID
# note: we dont assume the pre-existence of the table, just easy to code this way
w[,1] = part.data[,c('CueLeft')]
w[,2] = 2 - part.data$OdorLR
w[,3] = 0.5 # initialize all expectation to be 0.5 (or random?)

alpha = 0.8 # learning rate
# choice rule
exp_choice = function(x,y){
  exp(x)/sum(exp(x),exp(y))
}

ratio_choice = function(x,y){
  x/sum(x,y)
}

for(i in 1:ntrials){
CueLeft = use.data$CueLeft[i]
CueRight = use.data$CueRight[i]
# choice generation
wloc_left = which(w[,1] == CueLeft)
wloc_right = which(w[,1] == CueRight)
p_left = exp_choice(w[wloc_left,i+2],w[wloc_right,i+2]) # p of choosing left
choice[i] = 2 - rbinom(1,1,p_left)           # 1 for left, 2 for right

wloc_chosen = ifelse(choice[i]==1,wloc_left,wloc_right)
# update w of the chosen cue
w[,i+3] = w[,i+2]
reward = as.numeric(choice[i]== use.data$OdorLR[i])
w[wloc_chosen,i+3] = w[wloc_chosen,i+3] + 
  alpha * (reward - w[wloc_chosen,i+3])
}
acc = choice == use.data$OdorLR
w[,c(1,2,ntrials+3)]

mean(acc[1:24])
mean(acc[25:48])
mean(acc[49:72])
mean(acc[73:96])
mean(acc[97:120])
