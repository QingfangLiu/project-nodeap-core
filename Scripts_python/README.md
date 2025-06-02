# 🧠 TMS Effect Analysis via Variational Autoencoders (VAEs)

This repository contains code and analyses for uncovering the effects of Transcranial Magnetic Stimulation (TMS) using resting-state fMRI data and Variational Autoencoders (VAEs).

## 📌 Project Goal

To identify and characterize the effect of TMS on brain connectivity patterns by applying VAEs and Conditional VAEs (cVAEs) to resting-state fMRI-derived functional connectivity (FC) data.

## 📂 Data Overview

- **Functional Connectivity Types:**
  - AAL-based FC
  - AAL-sphere FC

- **Versions:**
  - Raw FC
  - Weighted FC (with/without resampling to handle class imbalance)

## 🧪 Model Variants

1. `VAE` on AAL-based FC  
2. `VAE` on AAL-sphere FC  
3. `Conditional VAE (cVAE)` on AAL-sphere FC, with subject ID added to both encoder and decoder

> ✅ Currently focusing on **Model 3**, which shows a clearer TMS effect and reduced subject clustering in latent space.

## 📊 Analyses

- Session-wise latent space distances to the **null (control)** condition
- Planned correlation with behavioral outcomes (behavioral summary data in progress)
- Exploratory use of `XGBoost` with one-hot encoded subject ID (not effective so far)

## 🛠️ Features

- Support for both raw and weighted input versions  
- Optional resampling to address TMS condition imbalance  
- Flexible integration of subject-level conditions  
- Modularized code for embedding, plotting, and model evaluation  

## ⏭️ Next Steps

- Generate behavioral effect summary for correlation analyses  
- Test alternative ROIs for latent space outputs  
- Improve robustness of conditional modeling  

## 📁 Repository Structure

```
├── data/                 # Preprocessed FC data (not public)
├── models/               # VAE and cVAE model definitions
├── analysis/             # Scripts for distance calculation and visualization
├── utils/                # Helper functions
├── results/              # Saved outputs (latent embeddings, plots, etc.)
├── README.md             # This file
```

## 🔒 Data Access

Note: Due to privacy and IRB restrictions, raw MRI data and some intermediate files are not publicly shared.

---

For questions or collaboration requests, feel free to contact the maintainer.
