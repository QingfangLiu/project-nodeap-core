# pip install pandas openpyxl
import pandas as pd, re, json
from pathlib import Path

# ====== INPUT FILE PATHS ======
f_org  = Path("/Users/liuq13/NODEAP/Data Organization_QL.xlsx")
f_coll = Path("/Users/liuq13/NODEAP/NODEAP_DataCollectionSheet_QL.xlsx")

# ====== OUTPUT FOLDER ======
out_dir = Path("/Users/liuq13/nodeap-bids")
out_dir.mkdir(parents=True, exist_ok=True)

out_tsv  = out_dir / "participants.tsv"
out_json = out_dir / "participants.json"

# ---- Read Excel ----
org  = pd.read_excel(f_org)
coll = pd.read_excel(f_coll)

# ---- Strip header spaces ----
org.columns = org.columns.str.strip()
coll.columns = coll.columns.str.strip()

# ---- Select relevant columns ----
org  = org[["Subject", "Sex", "Age"]].copy()
coll = coll[["Subject", "Anterior/Posterior", "Sweet/Savory", "Order"]].copy()

# normalize Subject *cell values* to prevent merge mismatches (e.g., trailing spaces)
org["Subject"]  = org["Subject"].astype(str).str.strip()
coll["Subject"] = coll["Subject"].astype(str).str.strip()

# ---- Merge ----
df = pd.merge(org, coll, on="Subject", how="inner")

# ---- Rename to BIDS fields ----
df = df.rename(columns={
    "Subject": "original_id",
    "Sex": "sex",
    "Age": "age",
    "Race": "race",
    "Ethnicity": "ethnicity",
    "Anterior/Posterior": "group",
    "Sweet/Savory": "meal_order",
    "Order": "tms_order"
})

# ---- Normalize values ----
df["sex"] = df["sex"].astype(str).str.upper().str[0]   # M/F

# Round age to nearest integer
df["age"] = df["age"].round(0).astype(int)

def map_group(v):
    s = str(v).lower()
    if "anterior" in s: return "aOFC"
    if "posterior" in s: return "pOFC"
    return v
df["group"] = df["group"].map(map_group)

def map_meal(v):
    s = str(v).lower()
    if s.startswith("sweet"): return "Sweet"
    if s.startswith("sav"): return "Savory"
    return v
df["meal_order"] = df["meal_order"].map(map_meal)

# ---- Assign participant_id ----
df["num"] = df["original_id"].str.extract(r"(\d+)$").astype(int)
df = df.sort_values("num").reset_index(drop=True)
df.insert(0, "participant_id", [f"sub-{i:02d}" for i in range(1, len(df)+1)])
df = df.drop(columns="num")

# ---- Save participants.tsv ----
df.to_csv(out_tsv, sep="\t", index=False)

# ---- Save participants.json ----
participants_json = {
  "participant_id": {"Description": "BIDS-compliant subject label"},
  "original_id": {"Description": "Original subject code (e.g., NODEAP_06)"},
  "sex": {"Description": "Biological sex of participant", "Levels": {"M": "Male","F": "Female","O": "Other"}},
  "age": {"Description": "Age of participant at time of experiment","Units": "years"},
  "group": {"Description": "Stimulation target group","Levels": {"aOFC": "Anterior OFC","pOFC": "Posterior OFC"}},
  "meal_order": {"Description": "First meal condition","Levels": {"Sweet": "Sweet meal","Savory": "Savory meal"}},
  "tms_order": {"Description": "Counterbalanced TMS order code"}
}
with open(out_json, "w") as f:
    json.dump(participants_json, f, indent=2)

print(f"✅ Saved {out_tsv}")
print(f"✅ Saved {out_json}")


desc = {
  "Name": "NODEAP fMRI dataset",
  "BIDSVersion": "1.9.0",
  "License": "CC0",
  "Authors": [
    "Qingfang Liu",
    "Daria Porter",
    "Hadeel Damra",
    "Yao Zhao",
    "Joel L. Voss",
    "Geoffrey Schoenbaum",
    "Thorsten Kahnt"
  ],
  "Acknowledgements": "We thank all study participants and lab members for their contributions. This work was supported by National Institute on Deafness and Other Communication Disorders grant R01DC015426 (to T.K.) and the Intramural Research Program at the National Institute on Drug Abuse (ZIA DA000642 to T.K. and DA000587 to G.S.). This research was supported in part by the Intramural Research Program of the National Institutes of Health (NIH). The contributions of the NIH authors are considered Works of the United States Government. The findings and conclusions presented in this dataset are those of the authors and do not necessarily reflect the views of the NIH or the U.S. Department of Health and Human Services.",
  "HowToAcknowledge": "Please cite our paper: Liu et al., 'Distinct contributions of anterior and posterior orbitofrontal cortex to outcome-guided behavior', *Current Biology* (in press). Add DOI when available.",
  "Funding": [
    "NIH NIDCD R01DC015426 (to T.K.)",
    "NIH NIDA ZIA DA000642 (to T.K.)",
    "NIH NIDA DA000587 (to G.S.)",
    "Intramural Research Program of the NIH"
  ],
  "ReferencesAndLinks": [
    "https://openneuro.org",
    "https://bids.neuroimaging.io/"
  ]
}


out_desc = out_dir / "dataset_description.json"
with open(out_desc, "w") as f:
    json.dump(desc, f, indent=2)

print(f"✅ Saved {out_desc}")

