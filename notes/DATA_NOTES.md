## Summary of Behavioral Data Issues
 
### Subject: `NODEAP_17`

- **Condition Assignment**: SC–CS–SS  
- **Behavioral Data Issue**:
  - Only completed **two runs of conditioning** in **Session 2** (cTBS–sham)
  - Behavioral data from Session 2 is **partial** and saved separately
  - For analysis, **Sessions 1 and 3** are included
  - Not used for analyzing behavioral data during conditioning task.
  - Could slightly enhance power for **sham–cTBS vs sham–sham** comparison
- **Neural Data Handling**:
  - **Neural data from Session 2 is included**, *except when correlating MRI with behavior*
- **Summary Action**:
  - Entire Session 2 treated as **missing for behavior**
  - Kept in MRI-level analyses where behavioral correlation is not required


### Subject: `NODEAP_44`

#### Session Schedule & Conditions

| Session | Date       | Stimulation | Devaluation       |
|---------|------------|-------------|-------------------|
| S1D1    | 2.13.23    | cTBS        | —                 |
| S1D2    | 2.14.23    | Sham        | Sweet             |
| S2D1    | 2.28.23    | Sham        | —                 |
| S2D2    | 3.1.23     | cTBS        | Sweet             |
| S3D1    | 3.13.23    | Sham        | —                 |
| S3D2    | 3.14.23    | Sham        | Savory            |

- **Note**:
  - Follows a **1–2–3 stimulation order**, which is consistent with counterbalancing.
  - However, the **devaluation order is unusual**:
    - Sweet was devalued twice (Sessions 1 and 2),
    - Savory only in Session 3.
- **Status**:
  - Have taken into this account in data analysis.

### Conditioning Data: General Notes

- One subject **missed one conditioning run**, but the data was **recovered from intermediate variables**.
- For **all other subjects and runs**, conditioning data in the **final dataset match** the intermediate data after checking.

---

# Summary of MRI Data Issues

## Missing Data

* `NODEAP_30`: Missing **S3D2** (Sham)
* `NODEAP_83`: Missing **S3D1** (cTBS)
* `NODEAP_44`: **S1D1** has partial acquisition (only 205 volumes) during cTBS scan
* `NODEAP_87`, `NODEAP_88`: Identical data for **D0** and **S1D1**

  * Resting-state data before S1D1 was not collected.
  * Both subjects had Sham on S1D1, so S1D1 was reused as D0 to calculate stimulation coordinates.
  * These scans should only be used **once** when analyzing across the 7-session timeline.
  * Only the compressed `D0_rest` folder is retained; D0 data should be treated as missing.
* `NODEAP_41`:

  * **What happened:** Session **S3D2** was acquired at **104×104×78**, while all other sessions for this subject are **104×104×60**.
  * **How we addressed it in the current pipeline:**

    1. Split means into two blocks for `NODEAP_41`:

       * Block 1: **D0..S3D1** → `functional/mean_fvol_1.nii`
       * Block 2: **S3D2** → `functional/mean_fvol_2.nii`
    2. **Coregister (estimate only)** each block’s mean to T1 and apply the transform to the corresponding 4D session files.

       * No `r*` files are created; SPM writes transform sidecars; original NIfTI headers remain intact.
    3. **Normalize (write)** all sessions to MNI at **3 mm** and **2 mm**, then **smooth (6 mm FWHM)**.

       * Per-session outputs: `w3fvol_4d.nii` → `s6w3fvol_4d.nii`, and `w2fvol_4d.nii` → `s6w2fvol_4d.nii`.
  * **Status:** `NODEAP_41` anomaly is resolved under the normalize-then-smooth pipeline.

## Handling in Analysis

* **Missing sessions** are treated as missing data when standardizing motion-related metrics (e.g., diff and var) across sessions.
* For `NODEAP_44` **S1D1**:

  * The 205 available volumes were included in the standardization.



