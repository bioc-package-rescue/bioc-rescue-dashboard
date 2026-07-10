#!/usr/bin/env python3
"""
update_reusable_workflows.py

For each package in packages.csv that has build status ERROR/TIMEOUT, writes
the minimal caller stub for the centralized check-bioc GHA workflow to the
local checkout and pushes to the remote.

Run sync_packages.py first to keep packages.csv up to date, and
clone_repos.py to ensure local checkouts exist.
"""
import csv
import os
import re
import subprocess
import time
import urllib.request

BUILD_DATABASES = {
    "Software":       "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt",
}

PACKAGES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "packages.csv")
WORKSPACE_ROOT = "/Users/Levi/git/bioc-package-rescue"

WORKFLOW_CONTENT = """\
name: R-CMD-check-bioc

on:
  push:
    branches:
      - main
      - master
      - devel
  pull_request:

jobs:
  run-check:
    uses: bioc-package-rescue/workflows/.github/workflows/check-bioc.yml@main
"""


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


def run(args, cwd=None):
    env = {k: v for k, v in os.environ.items() if k != "GITHUB_TOKEN"}
    return subprocess.run(args, capture_output=True, text=True, cwd=cwd, env=env)


def main():
    packages = load_packages()
    print(f"Loaded {len(packages)} packages from packages.csv.")

    print("Loading build status databases...")
    status = load_build_status()

    targets = [pkg for pkg in packages if status.get(pkg) == "error"]
    print(f"Found {len(targets)} packages with ERROR status.\n")

    success = failed = 0
    for i, pkg in enumerate(targets, 1):
        local = os.path.join(WORKSPACE_ROOT, pkg)
        print(f"[{i}/{len(targets)}] {pkg}...")

        if not os.path.exists(local):
            print(f"  [Skip] {local} does not exist — run clone_repos.py first.")
            failed += 1
            continue

        branch = run(["git", "branch", "--show-current"], cwd=local).stdout.strip()
        if not branch:
            print(f"  [Error] Could not determine active branch.")
            failed += 1
            continue
        print(f"  Branch: {branch}")

        wf_dir = os.path.join(local, ".github", "workflows")
        os.makedirs(wf_dir, exist_ok=True)
        wf_path = os.path.join(wf_dir, "check-bioc.yml")
        with open(wf_path, "w") as f:
            f.write(WORKFLOW_CONTENT)

        run(["git", "add", ".github/workflows/check-bioc.yml"], cwd=local)
        changed = ".github/workflows/check-bioc.yml" in run(["git", "status", "--porcelain"], cwd=local).stdout

        if changed:
            msg = "Use centralized check-bioc GHA workflow\n\nCo-authored-by: Antigravity <gemini@google.com>"
            r = run(["git", "commit", "-m", msg], cwd=local)
            if r.returncode != 0:
                print(f"  [Error] Commit failed: {r.stderr.strip()}")
                failed += 1
                continue
            print(f"  Committed.")

        r = run(["git", "push", "origin", branch], cwd=local)
        if r.returncode != 0:
            print(f"  [Error] Push failed: {r.stderr.strip()}")
            failed += 1
            continue

        print(f"  Done.")
        success += 1
        time.sleep(0.5)

    print(f"\nDone!  Updated: {success}  Failed: {failed}")


if __name__ == "__main__":
    main()
