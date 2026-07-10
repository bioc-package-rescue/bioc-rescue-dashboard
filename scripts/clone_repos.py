#!/usr/bin/env python3
"""
clone_repos.py

Clones each package in packages.csv that has build status ERROR/TIMEOUT from
the bioc-package-rescue GitHub organisation into the local workspace.

Run sync_packages.py first to keep packages.csv up to date, and
rescue_repos.py to ensure the GitHub repos exist before cloning.
"""
import csv
import os
import re
import subprocess
import urllib.request

BUILD_DATABASES = {
    "Software":       "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt",
}

PACKAGES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "packages.csv")
WORKSPACE_ROOT = "/Users/Levi/git/bioc-package-rescue"


def load_packages():
    with open(PACKAGES_CSV, newline="") as f:
        return [r["package"] for r in csv.DictReader(f)]


def load_build_status():
    status_db = {}
    for pkg_type, url in BUILD_DATABASES.items():
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req) as resp:
                content = resp.read().decode("utf-8")
            for line in content.splitlines():
                m = re.match(r"^([^#]+)#([^#]+)#([^:]+):\s*(.*)$", line.strip())
                if m:
                    pkg, _, _, st = [x.strip() for x in m.groups()]
                    if st in ("ERROR", "TIMEOUT"):
                        status_db[pkg] = "error"
                    elif pkg not in status_db:
                        status_db[pkg] = "OK"
        except Exception as e:
            print(f"  Warning: could not fetch {pkg_type} build DB: {e}")
    return status_db


def main():
    packages = load_packages()
    print(f"Loaded {len(packages)} packages from packages.csv.")

    print("Loading build status databases...")
    status = load_build_status()

    targets = [pkg for pkg in packages if status.get(pkg) == "error"]
    print(f"Found {len(targets)} packages with ERROR status.\n")

    success = skipped = failed = 0
    for i, pkg in enumerate(targets, 1):
        local = os.path.join(WORKSPACE_ROOT, pkg)
        print(f"[{i}/{len(targets)}] {pkg}...")
        if os.path.exists(local):
            print(f"  [Skip] {local} already exists.")
            skipped += 1
            continue

        ssh = f"git@github.com:bioc-package-rescue/{pkg}.git"
        print(f"  Cloning {ssh}...")
        r = subprocess.run(["git", "clone", ssh, local], capture_output=True, text=True)
        if r.returncode == 0:
            print(f"  Cloned successfully.")
            success += 1
        else:
            print(f"  Clone failed: {r.stderr.strip()}")
            failed += 1

    print(f"\nDone!  Cloned: {success}  Skipped: {skipped}  Failed: {failed}")


if __name__ == "__main__":
    main()
