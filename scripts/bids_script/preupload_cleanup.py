#!/usr/bin/env python3
"""
preupload_cleanup.py

Run this script before uploading to OpenNeuro:
- Removes "AcquisitionDuration" from all *_bold.json files.
- Deletes all .DS_Store files recursively.
"""

import json
import os
from pathlib import Path

bids_root = Path("/Users/liuq13/nodeap-bids")

fixed_json_count = 0
removed_dsstore_count = 0

# --- 1) Fix *_bold.json files ---
for jf in bids_root.rglob("*_bold.json"):
    with open(jf, "r") as f:
        data = json.load(f)
    if "AcquisitionDuration" in data:
        del data["AcquisitionDuration"]
        with open(jf, "w") as f:
            json.dump(data, f, indent=2)
        fixed_json_count += 1
        print("Fixed:", jf)

# --- 2) Remove .DS_Store files ---
for dsf in bids_root.rglob(".DS_Store"):
    try:
        os.remove(dsf)
        removed_dsstore_count += 1
        print("Removed:", dsf)
    except OSError as e:
        print("Could not remove:", dsf, e)

print("\nSummary:")
print(f"  JSON files fixed: {fixed_json_count}")
print(f"  .DS_Store files removed: {removed_dsstore_count}")
