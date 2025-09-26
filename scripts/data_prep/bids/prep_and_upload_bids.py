#!/usr/bin/env python3
"""
Run BIDS preparation and upload pipeline sequentially with retry logic.

Steps:
1. Run make_anat_bids.py to prepare anatomical data.
2. Run make_func_bids_me.py to prepare functional data.
3. Run preupload_cleanup.py to remove redundant files (e.g., .DS_Store).
4. Run OpenNeuro CLI (via Deno) to upload dataset, with up to 3 retries if it fails.

Modify dataset_id and dataset_dir as needed before running.
"""

import subprocess
import sys
import time
from pathlib import Path

# --- USER CONFIGURATION ---
dataset_id = "ds006693"
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

# --- STEP 4: UPLOAD TO OPENNEURO (with retry) ---
upload_cmd = [
    "deno", "run", "-A", "jsr:@openneuro/cli", "upload",
    "--dataset", dataset_id,
    "--affirmDefaced",
    str(dataset_dir)
]

max_retries = 3
attempt = 1
success = False

while attempt <= max_retries and not success:
    print(f">>> Attempt {attempt}/{max_retries}: Uploading dataset to OpenNeuro...")
    try:
        subprocess.run(upload_cmd, check=True)
        success = True
        print(">>> Upload completed successfully!")
    except subprocess.CalledProcessError:
        print(f"!!! Upload failed (attempt {attempt}).")
        if attempt < max_retries:
            wait_time = 5  # seconds before retry
            print(f"Retrying in {wait_time} seconds...")
            time.sleep(wait_time)
        attempt += 1

if not success:
    sys.exit("!!! Upload failed after maximum retries. Please check your Deno setup or network connection.")
