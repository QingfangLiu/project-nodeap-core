# 🧠 NODEAP Project data issues

## 📌 Subject: `NODEAP_17`
- **Issue**: Only completed **two runs of conditioning** in **Session 2**.
- **Action Taken**:
  - Entire **Session 2** treated as **missing data**.
  - **Sessions 1 and 3** are still included in the analysis.
  - This subject was assigned to the **SC–CS–SS** condition, so the **CS (Session 2)** data is missing.

---

## 📌 Subject: `NODEAP_44`

### 🗓️ Session Schedule & Conditions

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

---

## 🧪 Conditioning Data: General Notes

- One subject **missed one conditioning run**, but the data was **recovered from intermediate variables**.
- For **all other subjects and runs**, conditioning data in the **final dataset match** the intermediate data after checking.

---

## 🧲 Summary of MRI Data Issues

### 📉 Missing Data

- `NODEAP_30`: Missing **S3D2** (Sham)
- `NODEAP_83`: Missing **S3D1** (cTBS)
- `NODEAP_44`: **S1D1** has partial acquisition (only 205 volumes) during cTBS scan
- `NODEAP_87`, `NODEAP_88`: Identical data for **D0** and **S1D1**
  - Resting-state data before S1D1 was not collected
  - Both subjects had Sham on S1D1, so S1D1 was reused as D0 to calculate stimulation coordinates
  - These scans should only be used **once** when analyzing across the 7-session timeline
  - Only the compressed `D0_rest` folder is retained; D0 data should be treated as missing
- `NODEAP_41`: **S3D2** has different image dimensions (**104×104×78**)
  - Should be fine if a **normalization** step is applied
  - For **realignment**, handle this session separately
  - If using **native-space processing**, add an additional **reslicing step** (i.e., apply `+'r'` prefix)

### ⚙️ Handling in Analysis

- **Missing sessions** are treated as missing data when standardizing motion-related metrics (e.g., diff and var) across sessions
- For `NODEAP_44` **S1D1**:
  - The 205 available volumes were included in the standardization
  - Global connectivity was computed based on the available data

