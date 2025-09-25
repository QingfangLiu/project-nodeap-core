# fMRI Preprocessing – README

This folder contains the MATLAB/SPM pipeline used to preprocess the NODEAP resting-state fMRI (multi-echo) data, build subject-level means, normalize/smooth, and generate motion/slice nuisance regressors.

> For study-wide irregularities (missing sessions, special handling), **see `notes/DATA_NOTES.md`**.

---

## Requirements

* **MATLAB** (Statistics & Machine Learning Toolbox recommended)
* **SPM12** on your MATLAB path
* (Only for Step 1) **dcm2niix** installed and reachable from MATLAB via `system()`
  If you start from the OpenNeuro echo-wise NIfTIs, you can **skip Step 1**.

Directory assumptions (editable in scripts):

* `dat_folder` points at the dataset root (e.g., `/Volumes/X9Pro/NODEAP/MRI` or your OpenNeuro download).
* Per subject: `<dat_folder>/<SubID>/nifti/`
* Echo-wise sessions live in `<dat_folder>/<SubID>/nifti/<Sess>_me/` with files like
  `S1D1_rest_..._e1.nii`, `..._e2.nii`, `..._e3.nii`.

Session labels (columns throughout): `D0, S1D1, S1D2, S2D1, S2D2, S3D1, S3D2`.

---

## Run order (happy path)

1. **Convert DICOM → NIfTI (4D)**
   `Do_1_convert_dicoms_to_nifti_rest4D.m`

   * Converts T1 (`anat/`) and each rest session (`<Sess>_rest`) to NIfTI.
   * **Skip if using OpenNeuro NIfTIs.**

2. **Realign echo 1 across *all* sessions (estimate only)**
   `Do_2_realign_rest_e1_all_sessions.m`

   * Runs SPM **Realign (Estimate)** on echo-1 from all sessions together.
   * **No `r*` files are written;** SPM updates the header transform and writes `rp_*.txt`.

3. **Apply echo-1 transform to other echoes; reslice; PAID combine**
   `Do_3_apply_e1_realign_reslice_combine_PAID.m`

   * Copies the voxel-to-world mapping from e1 → e2/e3 (header update).
   * **Reslices** all echoes (writes `r*.nii`).
   * Computes **tSNR** per echo, tSNR×TE weights, and **PAID-weighted** combined volumes:
     `functional/<Sess>/fvol_###.nii` and a merged `fvol_4d.nii`.

4. **Count functional sessions (availability)**
   `Do_4_count_functional_sessions.m`

   * Writes `experiment_metadata/MRI_func_count.xlsx` with per-session availability
     (1 = 310 vols, 0.5 = 205 vols, 0 = missing).

5. **Build subject mean across sessions (native space)**

   * **Most subjects:** `Do_5_build_subject_mean_fMRI.m` → `functional/mean_fvol.nii`
   * **Special: `NODEAP_41`** (S3D2 has 104×104×78):
     `Do_5_build_subject_mean_NODEAP_41_split.m` →
     `functional/mean_fvol_1.nii` (D0..S3D1) and `functional/mean_fvol_2.nii` (S3D2)

6. **Coregister to T1, Normalize to MNI, Smooth**

   * **Most subjects:** `Do_6_post_PAID_coreg_norm_smooth.m`

     * Coreg (estimate) `mean_fvol.nii` → T1 (bring along session 4D).
     * Normalize **3 mm** (`w3fvol_4d.nii`) and **2 mm** (`w2fvol_4d.nii`) → Smooth **6 mm**
       (`s6w3fvol_4d.nii`, `s6w2fvol_4d.nii`).
   * **Special: `NODEAP_41`**: `Do_6_post_PAID_coreg_norm_smooth_NODEAP_41.m`

     * Coreg **both** means (`mean_fvol_1.nii`, `mean_fvol_2.nii`) to T1 and normalize/smooth sessions accordingly.
     * This resolves the S3D2 slice-count mismatch once in MNI space.

7. **Nuisance regressors (motion + slice metrics)**
   `make_nuisance_regressors_motion_slice.m`

   * From `rp_*.txt` (e1), builds 24 motion regressors (6, Δ6, 6², Δ6²).
   * Computes odd–even **slice diff** and across-slice **variance** from e1, flags outlier volumes, and writes one-hot regressors.
   * Outputs per session: `NRegressor/<SubID>/nuisance_regressors_<Sess>.txt` and diagnostics.

8. **Motion post-hoc (cTBS vs. sham)**
   `motion_posthoc_ctbs_vs_sham.m`

   * Recreates the manuscript motion analysis (mixed-effects; means ± SD).

---

## Special case: `NODEAP_41`

* **Issue:** Session **S3D2** was acquired at **104×104×78** (others are 104×104×60).
* **Handling in this pipeline:**

  * Step **5** uses the *split* script to create two means (`mean_fvol_1`, `mean_fvol_2`).
  * Step **6** uses the `NODEAP_41` variant to coregister/normalize each block separately.
  * After normalization, all sessions are on the same MNI grid; downstream analyses are standard.
* **If you insist on native-space cross-session work**, you must reslice S3D2 to the 60-slice grid before concatenation (not recommended—use the normalized outputs).

---

## Reproducing from **OpenNeuro** echo-wise NIfTIs

Because the **echo-wise NIfTI files** are publicly available, you can reproduce the processed fMRI from **Step 2 onward**:

1. Download the dataset and set `dat_folder` to your local OpenNeuro path.
2. Run **Step 2 → Step 6** (use the `NODEAP_41` variants where applicable).
3. Output products will match this pipeline’s naming (e.g., `functional/fvol_4d.nii`, `s6w3fvol_4d.nii`, `s6w2fvol_4d.nii`).
4. `dcm2niix` is **not required** in this route; **SPM12** is.

> Minor path edits: all scripts define `dat_folder` (and a few have `NR_dir` or metadata paths). Change those to your environment; no other code changes should be necessary.

---

## Outputs (key files)

* Per session (native): `functional/<Sess>/fvol_###.nii`, `fvol_4d.nii`
* Per subject (native): `functional/mean_fvol.nii` (or `_1` / `_2` for `NODEAP_41`)
* Normalized & smoothed: `w3fvol_4d.nii`, `s6w3fvol_4d.nii`, `w2fvol_4d.nii`, `s6w2fvol_4d.nii`
* Weights/maps: `nifti/tSNPmaps/tSNR_echo_<e>_<Sess>.nii`, `nifti/w_tSNR/w_tSNR_TE_echo_<e>_<Sess>.nii`
* Motion/nuisance: `NRegressor/<SubID>/nuisance_regressors_<Sess>.txt`, `SliceDiffVar.bmp`, `BadVolumes.mat`
* Session availability: `experiment_metadata/MRI_func_count.xlsx`

---

## Notes & provenance

* Detailed exceptions (missing/partial sessions, D0≡S1D1 substitutions, the `NODEAP_41` anomaly) are documented in **`notes/DATA_NOTES.md`**.
* Step 2 realignment is **estimate-only**; it updates header transforms but **does not** write `r*` images.
* Step 3 writes `r*` files (reslicing) and produces PAID-weighted combined data.

If anything looks off (e.g., missing `rp_*.txt`, session not found), check `DATA_NOTES.md` and confirm `dat_folder` points to the correct tree.
