#!/usr/bin/env python3
"""
rescue_repos.py

For each package in packages.csv with build status ERROR/TIMEOUT, creates a
rescue repository in the bioc-package-rescue GitHub organisation:
  - If the package has a GitHub upstream, it is forked.
  - Otherwise it is cloned from git.bioconductor.org and pushed as a new repo.

Run sync_packages.py first to keep packages.csv up to date.
"""
import csv
import os
import re
import shutil
import subprocess
import time
import urllib.request

BUILD_DATABASES = {
    "Software":       "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt",
}

PACKAGES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "packages.csv")
TEMP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp_clones")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_packages():
    """Returns list of (package, type) from packages.csv."""
    with open(PACKAGES_CSV, newline="") as f:
        return [(r["package"], r["type"]) for r in csv.DictReader(f)]


def load_build_status():
    """Returns {pkg: "error"|"OK"|"NA"} across all build databases."""
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


def get_github_url(package_name, pkg_type):
    if pkg_type == "Software":
        segment = "bioc"
    elif pkg_type == "ExperimentData":
        segment = "data/experiment"
    elif pkg_type == "AnnotationData":
        segment = "data/annotation"
    else:
        segment = "bioc"
    landing = f"https://bioconductor.org/packages/release/{segment}/html/{package_name}.html"
    try:
        req = urllib.request.Request(landing, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp:
            html = resp.read().decode("utf-8")
        links = re.findall(r'href="(https?://github\.com/[^"]+)"', html)
        if links:
            m = re.match(r"(https?://(?:www\.)?github\.com/[^/]+/[^/#?]+)", links[0])
            if m:
                url = m.group(1).rstrip("/")
                return url[:-4] if url.endswith(".git") else url
    except Exception:
        pass
    return f"https://git.bioconductor.org/packages/{package_name}"


def run(args, cwd=None):
    env = {k: v for k, v in os.environ.items() if k != "GITHUB_TOKEN"}
    return subprocess.run(args, capture_output=True, text=True, cwd=cwd, env=env)


def repo_exists(pkg):
    return run(["gh", "repo", "view", f"bioc-package-rescue/{pkg}"]).returncode == 0


def fork_repo(repo_url, pkg):
    print(f"  Forking {repo_url}...")
    r = run(["gh", "repo", "fork", repo_url, "--org", "bioc-package-rescue", "--clone=false"])
    if r.returncode == 0:
        print(f"  Forked successfully.")
        return True
    print(f"  Fork failed: {r.stderr.strip()}")
    return False


def clone_and_push(pkg):
    print(f"  Cloning from git.bioconductor.org...")
    os.makedirs(TEMP_DIR, exist_ok=True)
    tmp = os.path.join(TEMP_DIR, pkg)
    if os.path.exists(tmp):
        shutil.rmtree(tmp)

    r = run(["git", "clone", f"https://git.bioconductor.org/packages/{pkg}", tmp])
    if r.returncode != 0:
        print(f"  Clone failed: {r.stderr.strip()}")
        return False

    r = run(["gh", "repo", "create", f"bioc-package-rescue/{pkg}", "--public"])
    if r.returncode != 0:
        if repo_exists(pkg):
            print(f"  Repository already exists.")
            shutil.rmtree(tmp)
            return True
        print(f"  Create failed: {r.stderr.strip()}")
        shutil.rmtree(tmp)
        return False

    ssh = f"git@github.com:bioc-package-rescue/{pkg}.git"
    run(["git", "remote", "rename", "origin", "upstream"], cwd=tmp)
    run(["git", "remote", "add", "origin", ssh], cwd=tmp)
    run(["git", "push", "-u", "origin", "--all"], cwd=tmp)
    run(["git", "push", "-u", "origin", "--tags"], cwd=tmp)
    shutil.rmtree(tmp)
    print(f"  Cloned, created, and pushed successfully.")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    packages = load_packages()
    print(f"Loaded {len(packages)} packages from packages.csv.")

    print("Loading build status databases...")
    status = load_build_status()

    targets = [(pkg, ptype) for pkg, ptype in packages if status.get(pkg) == "error"]
    print(f"Found {len(targets)} packages with ERROR status.\n")

    success = skipped = failed = 0
    for i, (pkg, ptype) in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] {pkg}...")
        if repo_exists(pkg):
            print(f"  [Skip] bioc-package-rescue/{pkg} already exists.")
            skipped += 1
            continue

        repo_url = get_github_url(pkg, ptype)
        ok = fork_repo(repo_url, pkg) if "github.com" in repo_url else clone_and_push(pkg)
        (success if ok else failed).__class__  # just for reference
        if ok:
            success += 1
        else:
            failed += 1
        time.sleep(2)

    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR, ignore_errors=True)

    print(f"\nDone!  Forked/Created: {success}  Skipped: {skipped}  Failed: {failed}")


if __name__ == "__main__":
    main()
