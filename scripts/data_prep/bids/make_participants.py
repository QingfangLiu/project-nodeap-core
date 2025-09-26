# pip install pandas openpyxl
import pandas as pd, re, json
from pathlib import Path

# ====== INPUT FILE PATHS ======
f_org  = Path("/Volumes/X9Pro/NODEAP/experiment_metadata/Data Organization_QL.xlsx")
f_coll = Path("/Volumes/X9Pro/NODEAP/experiment_metadata/NODEAP_DataCollectionSheet_QL.xlsx")

# ====== OUTPUT FOLDER ======
out_dir = Path("/Users/liuq13/nodeap-bids")
out_dir.mkdir(parents=True, exist_ok=True)

out_tsv  = out_dir / "participants.tsv"
out_json = out_dir / "participants.json"
out_desc = out_dir / "dataset_description.json"

# ---- Read Excel ----
org  = pd.read_excel(f_org)
coll = pd.read_excel(f_coll)

# ---- Strip header spaces ----
org.columns  = org.columns.str.strip()
coll.columns = coll.columns.str.strip()

# ---- Select relevant columns ----
org  = org[["Subject", "Sex", "Age"]].copy()
coll = coll[["Subject", "Anterior/Posterior", "Sweet/Savory", "Order"]].copy()

# Normalize Subject cells
org["Subject"]  = org["Subject"].astype(str).str.strip()
coll["Subject"] = coll["Subject"].astype(str).str.strip()

# ---- Merge ----
df = pd.merge(org, coll, on="Subject", how="inner")

# ---- Rename to BIDS fields ----
df = df.rename(columns={
    "Subject": "original_id",
    "Sex": "sex",
    "Age": "age",
    "Anterior/Posterior": "group",
    "Sweet/Savory": "meal_order",     # will become 3-session sequence
    "Order": "tms_order"
})

# ---- Normalize values ----
df["sex"] = df["sex"].astype(str).str.upper().str[0]      # M/F
df["age"] = df["age"].round(0).astype(int)

def map_group(v):
    s = str(v).lower()
    if "anterior" in s:  return "aOFC"
    if "posterior" in s: return "pOFC"
    return v
df["group"] = df["group"].map(map_group)

def map_start_meal(v):
    s = str(v).strip().lower()
    if s.startswith("sweet"):  return "Sweet"
    if s.startswith("sav"):    return "Savory"
    return pd.NA
start_meal = df["meal_order"].map(map_start_meal)

# ---- Convert to 3-session sequence (alternate) ----
def alt3(start):
    if start == "Sweet":  return "Sweet-Savory-Sweet"
    if start == "Savory": return "Savory-Sweet-Savory"
    return pd.NA
df["meal_order"] = start_meal.map(alt3)

# ---- Special case: NODEAP_44 (devaluation Sweet, Sweet, then Savory) ----
mask_44 = df["original_id"].str.contains(r"NODE?E?AP[_-]?0?44$", case=False, regex=True)
df.loc[mask_44, "meal_order"] = "Sweet-Sweet-Savory"


# ---- Map numeric TMS order to letter codes ----
def map_tms_order(code):
    stim_map = {
        123: [['C','S'], ['S','C'], ['S','S']],  # session 1,2,3 (D1,D2)
        132: [['C','S'], ['S','S'], ['S','C']],
        213: [['S','C'], ['C','S'], ['S','S']],
        231: [['S','C'], ['S','S'], ['C','S']],
        312: [['S','S'], ['C','S'], ['S','C']],
        321: [['S','S'], ['S','C'], ['C','S']],
    }
    try:
        seq = stim_map[int(code)]
        # join pairs as two-letter strings, then join sessions with '-'
        return "-".join("".join(p) for p in seq)
    except Exception:
        return "NA"

df["tms_order"] = df["tms_order"].map(map_tms_order)

# ---- Assign participant_id ----
df["num"] = df["original_id"].str.extract(r"(\d+)$").astype(int)
df = df.sort_values("num").reset_index(drop=True).drop(columns="num")
df.insert(0, "participant_id", [f"sub-{i:02d}" for i in range(1, len(df)+1)])

# ---- Save participants.tsv ----
df.to_csv(out_tsv, sep="\t", index=False)

# ---- Save participants.json ----
participants_json = {
  "participant_id": {"Description": "BIDS-compliant subject label"},
  "original_id":    {"Description": "Original subject code (e.g., NODEAP_06)"},
  "sex":            {"Description": "Biological sex of participant", "Levels": {"M": "Male","F": "Female","O": "Other"}},
  "age":            {"Description": "Age of participant at time of experiment","Units": "years"},
  "group":          {"Description": "Stimulation target group","Levels": {"aOFC": "Anterior OFC","pOFC": "Posterior OFC"}},
  "meal_order":     {
    "Description": "Devaluation meal across the three sessions (Session1-Session2-Session3). Alternates from the initial meal unless overridden by subject-specific notes.",
    "Levels": {
      "Sweet-Savory-Sweet": "Sweet devalued in Sessions 1 & 3, Savory in Session 2",
      "Savory-Sweet-Savory": "Savory devalued in Sessions 1 & 3, Sweet in Session 2",
      "Sweet-Sweet-Savory":  "Special case (e.g., NODEAP_44): Sweet in Sessions 1 & 2, Savory in Session 3"
    }
  },
  "tms_order": {
    "Description": "TMS condition per session pair (D1+D2) across Session1-Session2-Session3, using C=cTBS and S=sham (e.g., 'CS-SC-SS').",
    "Levels": {
      "CS-SC-SS": "S1: cTBS+sham, S2: sham+cTBS, S3: sham+sham",
      "CS-SS-SC": "S1: cTBS+sham, S2: sham+sham, S3: sham+cTBS",
      "SC-CS-SS": "S1: sham+cTBS, S2: cTBS+sham, S3: sham+sham",
      "SC-SS-CS": "S1: sham+cTBS, S2: sham+sham, S3: cTBS+sham",
      "SS-CS-SC": "S1: sham+sham, S2: cTBS+sham, S3: sham+cTBS",
      "SS-SC-CS": "S1: sham+sham, S2: sham+cTBS, S3: cTBS+sham",
      "NA":        "Unspecified or invalid code"
    }
  }
}

with open(out_json, "w") as f:
    json.dump(participants_json, f, indent=2)

print(f"✅ Saved {out_tsv}")
print(f"✅ Saved {out_json}")

# ---- dataset_description.json (unchanged from your draft, shown for completeness) ----
desc = {
  "Name": "NODEAP fMRI dataset",
  "BIDSVersion": "1.9.0",
  "License": "CC0",
  "Authors": [
    "Qingfang Liu","Daria Porter","Hadeel Damra","Yao Zhao",
    "Joel L. Voss","Geoffrey Schoenbaum","Thorsten Kahnt"
  ],
  "Acknowledgements": "We thank all study participants and lab members for their contributions. This work was supported by National Institute on Deafness and Other Communication Disorders grant R01DC015426 (to T.K.) and the Intramural Research Program at the National Institute on Drug Abuse (ZIA DA000642 to T.K. and DA000587 to G.S.). This research was supported in part by the Intramural Research Program of the National Institutes of Health (NIH). The contributions of the NIH authors are considered Works of the United States Government. The findings and conclusions presented in this dataset are those of the authors and do not necessarily reflect the views of the NIH or the U.S. Department of Health and Human Services.",
  "HowToAcknowledge": "Please cite our paper: Liu et al., 'Distinct contributions of anterior and posterior orbitofrontal cortex to outcome-guided behavior', *Current Biology* (in press). Add DOI when available.",
  "Funding": [
    "NIH NIDCD R01DC015426 (to T.K.)",
    "NIH NIDA ZIA DA000642 (to T.K.)",
    "NIH NIDA DA000587 (to G.S.)",
    "Intramural Research Program of the NIH"
  ],
  "ReferencesAndLinks": ["https://openneuro.org","https://bids.neuroimaging.io/"]
}
with open(out_desc, "w") as f:
    json.dump(desc, f, indent=2)

print(f"✅ Saved {out_desc}")
