import json
from pathlib import Path

bids_root = Path("/Users/liuq13/nodeap-bids")

for jf in bids_root.rglob("*_bold.json"):
    with open(jf, "r") as f:
        data = json.load(f)
    if "AcquisitionDuration" in data:
        del data["AcquisitionDuration"]
        with open(jf, "w") as f:
            json.dump(data, f, indent=2)
        print("Fixed:", jf)
