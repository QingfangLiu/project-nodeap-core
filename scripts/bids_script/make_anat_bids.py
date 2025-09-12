# make_anat_bids.py
# Updated: copy → deface (on uncompressed .nii) → system gzip -n
# Requires: FSL on PATH (fsl_deface), pandas

import pandas as pd
from pathlib import Path
import subprocess, shutil, sys, os

# -------- PATHS (edit if needed) --------
raw_root  = Path("/Volumes/X9Pro/NODEAP/MRI")          # where NODEAP_* live
bids_root = Path("/Users/liuq13/nodeap-bids")          # BIDS dataset root
participants_tsv = bids_root / "participants.tsv"      # created earlier

bids_root.mkdir(parents=True, exist_ok=True)

# -------- Load subject mapping --------
pt = pd.read_csv(participants_tsv, sep="\t")          # needs participant_id, original_id
mapping = dict(zip(pt["original_id"], pt["participant_id"]))

# -------- helper: choose a T1w from anat/ --------
bad_tokens = ("seg", "echo", "tsnr", "inverse", "deformation")
bad_prefix = ("w", "y", "m")  # warped, inverse-warp, bias-corrected, etc.

def looks_like_t1(path: Path) -> bool:
    n = path.name.lower()
    if not n.endswith(".nii"):
        return False
    if any(tok in n for tok in bad_tokens):
        return False
    if path.name[0].lower() in bad_prefix:
        return False
    return ("t1" in n) or ("mprage" in n)

def pick_best_t1(choices):
    def score(p: Path): return len(p.name) + p.name.count("_")
    return sorted(choices, key=score)[0] if choices else None

def deface_inplace(nifti_path: Path):
    # write UNCOMPRESSED temp, then replace original .nii
    tmp = nifti_path.with_name(nifti_path.stem + "_defaced.nii")   # <- .nii, not .nii.gz
    subprocess.run(
        [sys.executable, "-m", "pydeface", str(nifti_path),
         "--outfile", str(tmp), "--force"],
        check=True
    )
    shutil.move(str(tmp), str(nifti_path))

# -------- iterate subjects --------
for orig_id, sub_id in mapping.items():
    dest_anat = bids_root / sub_id / "anat"
    dest_anat.mkdir(parents=True, exist_ok=True)

    bids_nii   = dest_anat / f"{sub_id}_T1w.nii"       # uncompressed
    bids_nii_gz = dest_anat / f"{sub_id}_T1w.nii.gz"   # final compressed
    bids_json  = dest_anat / f"{sub_id}_T1w.json"

    if bids_nii_gz.exists():
        print(f"• {orig_id} → {sub_id}: T1w already exists (gz); skipping")
        continue

    subj_dir = raw_root / orig_id / "nifti" / "anat"
    if not subj_dir.exists():
        print(f"⚠️  Skipping {orig_id}: no {subj_dir}")
        continue

    nii_candidates = [p for p in subj_dir.glob("*.nii") if looks_like_t1(p)]
    if not nii_candidates:
        print(f"⚠️  Skipping {orig_id}: no T1/MPRAGE NIfTI in {subj_dir}")
        continue

    t1 = pick_best_t1(nii_candidates)
    print(f"✓ {orig_id} -> {sub_id}: using {t1.name}")

    # Copy original .nii
    shutil.copy2(t1, bids_nii)

    # copy sidecar if present
    sidecar = t1.with_suffix(".json")
    if sidecar.exists():
        shutil.copy2(sidecar, bids_json)
    else:
        print(f"   ℹ️  no sidecar JSON for {t1.name} (OK; can add later)")

    # deface uncompressed .nii
    try:
        deface_inplace(bids_nii)
        print(f"   ✓ defaced {bids_nii.name}")
    except Exception as e:
        print(f"   ⚠️ deface failed for {sub_id}: {e}")
        continue

    # gzip with system tool (-n avoids timestamps)
    try:
        subprocess.run(["gzip", "-n", str(bids_nii)], check=True)
        print(f"   ✓ gzipped {bids_nii_gz.name}")
    except Exception as e:
        print(f"   ⚠️ gzip failed for {sub_id}: {e}")

print("Done.")
