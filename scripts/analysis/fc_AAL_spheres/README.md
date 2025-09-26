# FuncConn_AAL_spheres

This folder contains MATLAB scripts for computing functional connectivity (FC) matrices in the NODEAP study. The pipeline uses resliced AAL116 ROIs, tissue masks, and OFC seed/target spheres to compute ROI–ROI, sphere–sphere, and ROI–sphere connectivity across subjects and sessions.

---

## Dependencies
- MATLAB (R2020+ recommended)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) on MATLAB path  
- Statistics Toolbox (`zscore`, `corr`)  
- Preprocessed fMRI in **2×2×2 mm MNI space**  

---

## Paths used in scripts
Update these in each script if necessary:

```matlab
study_folder   = '/Volumes/X9Pro/NODEAP';
dat_folder     = '/Volumes/X9Pro/NODEAP/MRI';
project_folder = '/Users/liuq13/project-nodeap-core';
count_xlsx     = '/Volumes/X9Pro/NODEAP/experiment_metadata/MRI_func_count.xlsx';
```

---

## Order of scripts

### Step 0: Setup masks
1. **`step0_reslice_aal116.m`**  
   - Reslices `atlas_aal116/aal.nii` to the subject data grid (2 mm).  
   - Output: `atlas_masks/atlas_aal116/raal.nii`.  

2. **`step0_make_tissue_masks.m`**  
   - Reslices SPM TPMs to 2 mm, thresholds GM/WM/CSF.  
   - Outputs in `atlas_masks/tissue_masks/`:  
     - `gm_0.1_2mm.nii`  
     - `wm_0.9_2mm.nii`  
     - `csf_0.9_2mm.nii`  
     - `r2mmTPM.nii`  

---

### Step 1: Main AAL FC pipeline
3. **`step1_build_AAL_FC.m`**  
   - Loads preprocessed fMRI (`s6w2*.nii`).  
   - Extracts GM voxels, regresses nuisance (motion + tissue signals + drift).  
   - Saves filtered GM voxel time series.  
   - Averages signals into AAL ROIs → ROI×ROI FC matrix.   

---

### Step 2: Spheres analyses
4. **`step2_fc_sphere2sphere.m`**  
   - Uses filtered GM data.  
   - Extracts time series from 4 OFC spheres (aOFC seed/target, pOFC seed/target).  
   - Computes 4×4 correlation matrix.  

5. **`step2_fc_roi2sphere.m`**  
   - Correlates each AAL ROI signal with each of the 4 OFC spheres.  
   - Produces ROI×4 correlation matrix.  

---

## Notes
- Session availability is controlled via `count_xlsx`. Sessions with value ≠ 1 are skipped.  
- Run scripts in order: **step0 → step1 → step2**.  
