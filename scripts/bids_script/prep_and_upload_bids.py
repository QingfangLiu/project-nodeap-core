#!/usr/bin/env python3
"""
Run BIDS preparation and upload pipeline sequentially.

Steps:
1. Run make_anat_bids.py to prepare anatomical data.
2. Run make_func_bids_me.py to prepare functional data.
3. Run preupload_cleanup.py to remove redundant files (e.g., .DS_Store).
4. Run OpenNeuro CLI (via Deno) to upload dataset.

Modify dataset_id and dataset_dir as needed before running.
"""

import subprocess
from pathlib import Path

# --- USER CONFIGURATION ---
dataset_id = "ds006671"
dataset_dir = Path("/Users/liuq13/nodeap-bids/")

# --- STEP 1: MAKE ANAT BIDS ---
print(">>> Running make_anat_bids.py...")
subprocess.run(["python3", "make_anat_bids.py"], check=True)

# --- STEP 2: MAKE FUNC BIDS ---
print(">>> Running make_func_bids_me.py...")
subprocess.run(["python3", "make_func_bids_me.py"], check=True)

# --- STEP 3: CLEANUP ---
print(">>> Running preupload_cleanup.py...")
subprocess.run(["python3", "preupload_cleanup.py"], check=True)

# --- STEP 4: UPLOAD TO OPENNEURO ---
print(">>> Uploading dataset to OpenNeuro...")
upload_cmd = [
    "deno", "run", "-A", "jsr:@openneuro/cli", "upload",
    "--dataset", dataset_id,
    "--affirmDefaced",
    str(dataset_dir)
]
subprocess.run(upload_cmd, check=True)

print(">>> All steps completed successfully!")
