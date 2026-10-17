#!/usr/bin/env python3
"""Optional shallow clones; preserve licenses and record actual commit IDs."""
from pathlib import Path
import argparse
import datetime
import json
import subprocess

course = Path(__file__).resolve().parents[1]
catalog = json.loads((course / "resources/tutorials.json").read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("action", nargs="?", choices=["list", "fetch"], default="list")
parser.add_argument("name", nargs="?", choices=list(catalog))
args = parser.parse_args()
if args.action == "list":
    for key, value in catalog.items():
        print(key, value)
else:
    if not args.name:
        parser.error("fetch requires one tutorial name")
    dest = course / "resources/upstream" / args.name
    if dest.exists():
        raise SystemExit(f"Keeping existing directory without changes: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "clone", "--depth", "1", catalog[args.name], str(dest)], check=True)
    commit = subprocess.check_output(["git", "-C", str(dest), "rev-parse", "HEAD"], text=True).strip()
    log = course / "resources/upstream_receipts.json"
    data = json.loads(log.read_text()) if log.exists() else []
    data.append({"name": args.name, "url": catalog[args.name], "commit": commit,
                 "time": datetime.datetime.now(datetime.timezone.utc).isoformat()})
    log.write_text(json.dumps(data, indent=2) + "\n")
    print("Cloned with original licenses retained:", dest, commit)
