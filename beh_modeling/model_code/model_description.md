# Hierarchical Learning Model (by Session, with Cue-Pair Tracking)

This document describes the JAGS model implemented in  
`model_code/learning_model.txt`.  
It explains the purpose, inputs, update rules, observation model,  
hierarchical priors, and overall generative story.

---

## Purpose
The model estimates how participants learn cue–outcome associations  
in a behavioral task. It tracks the evolution of associative strength `w`  
for each **subject × session × cue pair** across trials,  
and links these latent values to observed trial-by-trial responses.  
The model allows learning rates to vary hierarchically across cue pairs,  
sessions, and subjects.

---

## Data Structure & Inputs
- **Indices**
  - `j = 1…nsubs` : subjects  
  - `c = 1…nsess` : sessions  
  - `i = 1…ntrials` : trials  
  - `p = 1…ncuepair` : cue pairs  

- **Observed**
  - `Resp[j,c,i]` : binary trial response (0/1)  
  - `CuePair[j,c,i]` : which cue pair was presented at trial *i*  

- **Latent**
  - `w[j,c,i,p]` : associative strength of cue pair *p* at trial *i*  

---

## State Update Equation
On each trial, only the presented cue pair updates:

\[
w_{j,c,i+1,p} \;=\; w_{j,c,i,p} \;+\; \mathbf{1}\{ \text{CuePair}_{j,c,i}=p \}\;
\alpha_{j,c,p}\,[\,1 - w_{j,c,i,p}\,].
\]

- If pair *p* is shown: `w` moves toward 1 with rate `α[j,c,p]`.  
- If not: `w` simply carries over unchanged.  
- Initial `w` values must be provided as input.

---

## Observation Model
Each response is drawn from a Bernoulli distribution  
with success probability equal to the latent `w` of the shown pair:

\[
\text{Resp}_{j,c,i} \sim \text{Bernoulli}\big(w_{j,c,i,\;\text{CuePair}_{j,c,i}}\big).
\]

---

## Hierarchical Priors

- **Session means**  
  \[
  \mu_c \sim \text{Beta}(8,2), \quad \text{truncated to } (0.001,0.999)
  \]

- **Global concentration**  
  \[
  \kappa \sim \text{Gamma}(1,0.1)
  \]

- **Subject-level concentration**  
  \[
  k_j \sim \text{Gamma}(1,0.1)
  \]

- **Subject–session means**  
  \[
  a_{j,c} \sim \text{Beta}(\mu_c \kappa,\,(1-\mu_c)\kappa),
  \quad \text{truncated to } (0.001,0.999)
  \]

- **Cue-pair–specific learning rates**  
  \[
  \alpha_{j,c,p} \sim \text{Beta}(a_{j,c} k_j,\,(1-a_{j,c})k_j),
  \quad \text{truncated to } (0.001,0.999)
  \]

---

## Generative Story
1. Draw session-level mean learning rates (`μ_c`) and hyperparameters (`κ`, `k_j`).  
2. For each subject × session, draw a mean rate `a[j,c]`.  
3. For each cue pair, draw its learning rate `α[j,c,p]`.  
4. Initialize latent values `w`.  
5. For each trial: update `w` of the presented cue pair, and generate a binary response.  

---

## Outputs
- Posterior estimates of:  
  - Cue-pair–specific learning rates `α[j,c,p]`  
  - Subject–session means `a[j,c]`  
  - Session means `μ[c]`  
  - Concentrations `k_j` and `κ`  
- Posterior predictive distributions of responses.  

