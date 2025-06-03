
# this code wants to simulate the learning process of discrimination task
# learning is on each pair of two cues (inseparable) 

Simu_learning = function(alpha){

nt = 10 # number of time points
S = numeric(nt) # strength (knowledge of which cue is correct)
S[1] = 0.5 # initialize as 0.5
delta = numeric(nt)

for(t in 1:nt){
  delta[t] = 1 - S[t]
  S[t+1] = S[t] + alpha * delta[t]
}

return(S)
}

# since S is 0-1, simply use S as the prop of selecting correctly
# shape mimicking the observed probability of selecting odors

alpha = 0.5
S=Simu_learning(alpha)
plot(S,type = 'l')


alpha_vec = seq(0.1,0.9,0.1)

out=sapply(alpha_vec,Simu_learning)
matplot(out,type = 'l',ylab = 'Prob of correct',xlab = 'Trials')


