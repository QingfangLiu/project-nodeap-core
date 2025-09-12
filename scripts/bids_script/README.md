# NODEAP BIDS Conversion & Upload Pipeline

This folder contains all scripts needed to prepare and upload the NODEAP fMRI dataset to [OpenNeuro](https://openneuro.org/).

---

## 1. Generate Participant Files (`make_participants.py`)

- Reads participant demographic and experimental information from:
  - `Data Organization_QL.xlsx`
  - `NODEAP_DataCollectionSheet_QL.xlsx`
- Produces:
  - `participants.tsv` — required BIDS participant table
  - `participants.json` — metadata descriptions for each column
  - `dataset_description.json` — required BIDS dataset metadata

Run this script **first** before generating any imaging files.

---

## 2. Prepare Anatomical Data (`make_anat_bids.py`)

- Locates T1-weighted structural scans (`T1`, `MPRAGE`) under each subject’s `nifti/anat` folder
- Copies the best candidate `.nii` to the BIDS `anat/` folder
- Copies sidecar JSON if available
- **Defaces** the anatomical image using `pydeface` on the uncompressed `.nii`
- **Compresses** with system `gzip -n` (avoids timestamp differences for reproducibility)
- Skips subjects already processed (`*_T1w.nii.gz` present)

This script ensures consistent header handling and minimizes OpenNeuro validation errors.

---

## 3. Prepare Functional Multi-Echo Data (`make_func_bids_me.py`)

- Finds all `*_me` session folders under each subject’s raw `nifti` directory
- Groups runs by stem (excluding echo number), preserving echo numbers
- Copies each echo to the BIDS `func/` folder
- **Compresses** using the same `gzip -n` approach as anatomical data
- Copies sidecar JSONs if available

This step mirrors the anatomical workflow, ensuring that headers remain valid and avoiding “cannot parse NIfTI header” validator errors.

---

## 4. Final Cleanup Before Upload (`preupload_cleanup.py`)

- Removes `"AcquisitionDuration"` field from all `*_bold.json` files
- Deletes all `.DS_Store` files recursively
- Prints a summary count of fixed JSONs and removed `.DS_Store` files

Run this script **immediately before uploading** to guarantee a clean BIDS tree.

---

## 5. Upload to OpenNeuro

Use the OpenNeuro CLI (via Deno):

```bash
deno run -A jsr:@openneuro/cli upload \
  --dataset ds006671 \
  --affirmDefaced \
  /Users/liuq13/nodeap-bids/
```

---

## Notes & Best Practices

- **Defacing:** `pydeface` is safer and more BIDS-friendly than `fsl_deface`.
- **Compression:** always use system `gzip -n` for reproducibility and validator compliance.
- **No double compression:** never write `.nii.gz` then gzip again — this corrupts headers.

