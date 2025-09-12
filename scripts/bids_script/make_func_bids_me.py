# make_func_bids_me.py
# Copy multi-echo RAW fMRI from /Volumes/X9Pro/NODEAP/MRI/<ID>/nifti/*_me
# → /Users/liuq13/nodeap-bids/sub-XX/ses-<SESSION>/func with BIDS names.

import re, gzip, shutil
from pathlib import Path
import pandas as pd

# ---- paths ----
RAW_ROOT  = Path("/Volumes/X9Pro/NODEAP/MRI")
BIDS_ROOT = Path("/Users/liuq13/nodeap-bids")
PARTICIPANTS_TSV = BIDS_ROOT / "participants.tsv"
TASK = "rest"  # change if you want another task label

BIDS_ROOT.mkdir(parents=True, exist_ok=True)

# subject mapping: original_id -> participant_id
pt = pd.read_csv(PARTICIPANTS_TSV, sep="\t")
mapping = dict(zip(pt["original_id"].astype(str), pt["participant_id"].astype(str)))

# regex to find echo number (_e1.nii / -e1.nii)
ECHO_RE = re.compile(r"(?:^|[_\-])e(\d+)(?:[_\-]|\.|$)", re.IGNORECASE)

def is_raw_nii(p: Path) -> bool:
    n = p.name
    if not n.lower().endswith(".nii"):
        return False
    # skip true SPM-derived prefixes ONLY (lowercase)
    if n.startswith(("r", "w", "y", "s")):
        return False
    # optional: also skip obvious SPM aux files
    if n.lower().startswith(("rp_", "spm_")):
        return False
    return True

def stem_without_echo(p: Path) -> str:
    """remove the `_e<n>` token from the stem for run grouping"""
    return ECHO_RE.sub("", p.stem)

for orig_id, sub_id in mapping.items():
    subj_root = RAW_ROOT / orig_id / "nifti"
    if not subj_root.exists():
        print(f"⚠️  {orig_id}: {subj_root} not found, skipping")
        continue

    # find session folders like D0_me, S1D1_me, ...
    session_dirs = [d for d in subj_root.iterdir() if d.is_dir() and d.name.endswith("_me")]
    if not session_dirs:
        print(f"⚠️  {orig_id}: no *_me session folders under {subj_root}")
        continue

    for sdir in sorted(session_dirs, key=lambda p: p.name):
        session = sdir.name[:-3]  # strip "_me" → D0, S1D1, ...
        dest = BIDS_ROOT / sub_id / f"ses-{session}" / "func"
        dest.mkdir(parents=True, exist_ok=True)

        # candidate raw NIfTIs in this session
        cands = [p for p in sdir.glob("*.nii") if is_raw_nii(p) and ECHO_RE.search(p.name)]
        if not cands:
            print(f"ℹ️  {orig_id} ses-{session}: no raw multi-echo files found")
            continue

        # group by run using stem-without-echo
        groups = {}
        for p in cands:
            key = stem_without_echo(p)
            groups.setdefault(key, []).append(p)

        # write each run
        run_idx = 1
        for key, files in sorted(groups.items()):
            # map echo -> file
            echo_map = {}
            for f in files:
                m = ECHO_RE.search(f.name)
                if not m: 
                    continue
                echo = int(m.group(1))
                echo_map[echo] = f

            # copy each echo
            for echo in sorted(echo_map.keys()):
                src = echo_map[echo]
                out_bold = dest / f"{sub_id}_ses-{session}_task-{TASK}_run-{run_idx:02d}_echo-{echo}_bold.nii.gz"
                out_json = dest / f"{sub_id}_ses-{session}_task-{TASK}_run-{run_idx:02d}_echo-{echo}_bold.json"

                # skip if already present
                if out_bold.exists() and out_json.exists():
                    continue

                # write .nii.gz
                with open(src, "rb") as fin, gzip.open(out_bold, "wb") as fout:
                    shutil.copyfileobj(fin, fout)

                # copy sidecar json if available
                sidecar = src.with_suffix(".json")
                if sidecar.exists():
                    shutil.copy2(sidecar, out_json)
                else:
                    print(f"ℹ️  {orig_id} ses-{session} run-{run_idx:02d} echo-{echo}: missing sidecar JSON")

            print(f"✓ {orig_id} → {sub_id}: ses-{session} run-{run_idx:02d} with {len(echo_map)} echoes")
            run_idx += 1

print("Done (functional copy).")
