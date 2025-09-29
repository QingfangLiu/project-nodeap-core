# Project: NODEAP Core - Distinct contributions of anterior and posterior OFC to outcome-guided behavior

This NODEAP study tests differential roles of the anterior (aOFC) and posterior (pOFC) lateral orbitofrontal cortex in outcome-guided (goal-directed) behavior. Participants completed odor-reward learning and outcome devaluation across multiple sessions. Network-targeted TMS was applied to aOFC or pOFC in separate sessions to probe effects on learning vs. inference.

This repository contains all code, behavioral data, MRI derivatives, and figures for the NODEAP project.  
Below is the folder structure with brief descriptions.

```
project-nodeap-core/
├─ atlas_masks/                     # Atlas files and ROI/tissue masks for fMRI analysis
│  ├─ atlas_aal116/                 # AAL116 atlas in MNI space
│  ├─ ofc_connectivity_masks/       # aOFC/pOFC seeds and targets for connectivity analysis
│  └─ tissue_masks/                 # Tissue probability / binary masks for preprocessing and denoising
│
├─ beh_data_found/             		# Raw behavioral data collected after organization
│
├─ beh_data_processed/         		# Cleaned & processed behavioral datasets
│
├─ beh_modeling/               		# Model fitting scripts and outputs 
│
├─ Figs_final/                		# Final publication-ready figures
│
├─ Figs_generated/             		# Figures saved directly by R/Python code
│
├─ scripts/                         # All analysis and preprocessing code
│  ├─ beh_analysis/              	# Behavioral statistical analysis
│  ├─ beh_format/              		# Format behavioral data for downstream analysis
│  ├─ fc_AAL_spheres/           	# Functional connectivity analyses using AAL ROI spheres
│  ├─ mri_postproc_VAE/             # MRI analysis from VAE outputs 
│  ├─ get_coordinates/           	# Extract/convert coordinates for analysis
│  ├─ fmri_preprocessing/        	# Preprocessing pipeline
│  ├─ utils/                        # Helper functions, configs
└─ README.md                   		# Project overview (this file)
```

- For details on the hierarchical Bayesian behavioral modeling used in this project, see the dedicated [README in `beh_modeling/`](./beh_modeling/README.md).
- The resting-state fMRI data have been publicly deposited on [OpenNeuro (ds006693)](https://openneuro.org/datasets/ds006693).  
- The cVAE analysis of resting-state fMRI data is available in a [standalone repository](https://github.com/QingfangLiu/project-nodeap-fmri-cvae) containing session-wise processed functional connectivity data.  


## Citation
Please cite our paper if you use this dataset or code:
Liu et al., *Distinct contributions of anterior and posterior orbitofrontal cortex to outcome-guided behavior*,  
*Current Biology* (in press). DOI will be added when available.

## Data Notes
Some sessions contain irregularities (missing data, experiment errors, or special conditions for specific subjects).  
See [DATA_NOTES.md](notes/DATA_NOTES.md) for the full log.

## Contact
Qingfang Liu  
Research Fellow, National Institute on Drug Abuse (NIDA) 
Email: psychliuqf@gmail.com


