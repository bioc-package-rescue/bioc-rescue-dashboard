#!/usr/bin/env python3
"""
sync_repos.py

Syncs each rescue repository with Bioconductor before any rescue changes are
made. The `devel` branch on git.bioconductor.org is the reference: it is never
rewritten, and the rescue default branch is rebuilt on top of it.

For each local checkout (see clone_repos.py):

  1. Fetch `devel` from https://git.bioconductor.org/packages/<pkg> into the
     `bioc` remote, and fetch `origin` (the rescue repo).
  2. Collect the commits on origin/<default> that are not on bioc/devel,
     skipping merge commits and commits whose patch is already on bioc/devel
     (e.g. the same fix merged twice through different PRs).
  3. Rebuild <default> as bioc/devel plus those commits, replayed oldest first
     with `git cherry-pick`. The result is linear: no merge commits above
     bioc/devel.
  4. Without --push, leave the result on a local `bioc-sync` branch for
     review and restore the checkout. With --push: save the old tip as `pre-sync/<YYYYMMDD>` on origin,
     force-push <default> with --force-with-lease, and turn off merge commits
     for the repo's PRs (only rebase and squash merging stay on) so the
     history stays linear.

Without --push nothing leaves the local checkout and <default> is not
touched. If a replayed commit
conflicts, the checkout is restored and the package is reported for a manual
sync. A dropped merge commit that resolved conflicts by hand may carry changes
that are not replayed, so review the printed commit list before pushing.

Usage:
  python scripts/sync_repos.py [--push] [--workspace DIR] [PKG ...]

With no PKG, all packages in packages.csv with build status ERROR/TIMEOUT are
synced, as in the other scripts.
"""
import argparse
import csv
import datetime
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

BIOC_REMOTE = "bioc"
BIOC_BRANCH = "devel"
ORG = "bioc-package-rescue"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def git(cwd, *args):
    return run(["git", *args], cwd=cwd)


def git_out(cwd, *args):
    r = git(cwd, *args)
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {r.stderr.strip()}")
    return r.stdout.strip()


def default_branch(local):
    r = git(local, "symbolic-ref", "--short", "refs/remotes/origin/HEAD")
    if r.returncode != 0:
        git(local, "remote", "set-head", "origin", "--auto")
        r = git(local, "symbolic-ref", "--short", "refs/remotes/origin/HEAD")
    if r.returncode != 0:
        return None
    return r.stdout.strip().split("/", 1)[1]


def ensure_bioc_remote(local, pkg):
    url = f"https://git.bioconductor.org/packages/{pkg}"
    remotes = git_out(local, "remote").splitlines()
    if BIOC_REMOTE not in remotes:
        git_out(local, "remote", "add", BIOC_REMOTE, url)
    else:
        current = git_out(local, "remote", "get-url", BIOC_REMOTE)
        if current.rstrip("/").removesuffix(".git") != url:
            print(f"  Note: remote '{BIOC_REMOTE}' points to {current}, not {url}.")


def is_ancestor(local, a, b):
    return git(local, "merge-base", "--is-ancestor", a, b).returncode == 0


def is_clean(local):
    return git_out(local, "status", "--porcelain", "--untracked-files=no") == ""


def replay(local, base, commits):
    """Cherry-pick commits onto base on a detached HEAD.

    Returns (ok, picked, dropped, failed_sha).
    """
    git_out(local, "checkout", "-q", "--detach", base)
    picked, dropped = [], []
    for sha in commits:
        r = git(local, "cherry-pick", "--allow-empty-message", sha)
        if r.returncode == 0:
            picked.append(sha)
            continue
        # A commit that becomes empty duplicates a change already applied.
        no_changes = (git(local, "diff", "--quiet").returncode == 0 and
                      git(local, "diff", "--cached", "--quiet").returncode == 0)
        if no_changes:
            git(local, "cherry-pick", "--skip")
            dropped.append(sha)
            continue
        git(local, "cherry-pick", "--abort")
        return False, picked, dropped, sha
    return True, picked, dropped, None


def disable_merge_commits(pkg):
    r = run(["gh", "api", "-X", "PATCH", f"repos/{ORG}/{pkg}",
             "-F", "allow_merge_commit=false",
             "-F", "allow_rebase_merge=true",
             "-F", "allow_squash_merge=true"])
    if r.returncode != 0:
        print(f"  Warning: could not turn off merge commits: {r.stderr.strip()}")


# ---------------------------------------------------------------------------
# Per-package sync
# ---------------------------------------------------------------------------

def sync_package(pkg, workspace, push):
    """Returns one of: 'synced', 'updated', 'skipped', 'failed'."""
    local = os.path.join(workspace, pkg)
    if not os.path.isdir(local):
        print(f"  [Skip] {local} does not exist; run clone_repos.py first.")
        return "skipped"
    if not is_clean(local):
        print(f"  [Skip] {local} has uncommitted changes.")
        return "skipped"

    ensure_bioc_remote(local, pkg)
    r = git(local, "fetch", "-q", BIOC_REMOTE, BIOC_BRANCH)
    if r.returncode != 0:
        print(f"  [Error] fetch {BIOC_REMOTE} {BIOC_BRANCH}: {r.stderr.strip()}")
        return "failed"
    git_out(local, "fetch", "-q", "origin")

    branch = default_branch(local)
    if not branch:
        print("  [Error] Could not determine the default branch of origin.")
        return "failed"

    bioc = f"{BIOC_REMOTE}/{BIOC_BRANCH}"
    remote_ref = f"origin/{branch}"
    old_tip = git_out(local, "rev-parse", remote_ref)
    bioc_tip = git_out(local, "rev-parse", bioc)
    print(f"  {remote_ref} {old_tip[:8]}  {bioc} {bioc_tip[:8]}")

    r = git(local, "rev-parse", "--verify", "-q", branch)
    if r.returncode == 0 and r.stdout.strip() != old_tip:
        print(f"  [Skip] local '{branch}' differs from {remote_ref}; push or reset it first.")
        return "skipped"

    merges = int(git_out(local, "rev-list", "--count", "--merges", f"{bioc}..{remote_ref}"))
    if is_ancestor(local, bioc, remote_ref) and merges == 0:
        ahead = git_out(local, "rev-list", "--count", f"{bioc}..{remote_ref}")
        print(f"  In sync: {ahead} commit(s) on top of {bioc}, no merges.")
        return "synced"

    original = git_out(local, "rev-parse", "--abbrev-ref", "HEAD")
    target = branch if push else "bioc-sync"

    if is_ancestor(local, remote_ref, bioc):
        # Behind only: fast-forward, nothing to replay.
        new_tip = bioc_tip
        print(f"  Behind {bioc}: fast-forward.")
        git_out(local, "checkout", "-q", "-B", target, new_tip)
        force = False
    else:
        commits = git_out(local, "rev-list", "--reverse", "--topo-order", "--no-merges",
                          "--cherry-pick", "--right-only",
                          f"{bioc}...{remote_ref}").split()
        print(f"  Replaying {len(commits)} commit(s) onto {bioc}; "
              f"dropping {merges} merge commit(s).")
        ok, picked, dropped, failed_sha = replay(local, bioc, commits)
        if not ok:
            git(local, "checkout", "-q", original if original != "HEAD" else old_tip)
            subject = git_out(local, "log", "-1", "--format=%h %s", failed_sha)
            print(f"  [Conflict] {subject}")
            print("  Restored the checkout. Sync this package by hand.")
            return "failed"
        new_tip = git_out(local, "rev-parse", "HEAD")
        git_out(local, "checkout", "-q", "-B", target, new_tip)
        for sha in dropped:
            print("    dropped (already applied): " + git_out(local, "log", "-1", "--format=%h %s", sha))
        force = True

    print("  New history on top of " + bioc + ":")
    log = git_out(local, "log", "--format=%h %s", f"{bioc}..{target}")
    for line in (log.splitlines() or ["(none)"]):
        print("    " + line)

    if not push:
        git(local, "checkout", "-q", original if original != "HEAD" else old_tip)
        print(f"  Dry run: result is on local branch '{target}'; '{branch}' and origin "
              "are unchanged. Re-run with --push to apply.")
        return "updated"

    if force:
        backup = f"pre-sync/{datetime.date.today():%Y%m%d}"
        r = git(local, "push", "-q", "origin", f"{old_tip}:refs/heads/{backup}")
        if r.returncode != 0:
            print(f"  [Error] backup push failed: {r.stderr.strip()}")
            return "failed"
        print(f"  Saved old tip as origin/{backup}.")
        r = git(local, "push", "-q", f"--force-with-lease={branch}:{old_tip}", "origin", branch)
    else:
        r = git(local, "push", "-q", "origin", branch)
    if r.returncode != 0:
        print(f"  [Error] push failed: {r.stderr.strip()}")
        return "failed"
    disable_merge_commits(pkg)
    print(f"  Pushed {branch}. Rebase open PR branches onto it.")
    return "updated"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("packages", nargs="*", help="packages to sync (default: ERROR/TIMEOUT packages)")
    ap.add_argument("--push", action="store_true", help="push the synced default branch to origin")
    ap.add_argument("--workspace", default=WORKSPACE_ROOT, help=f"local checkouts (default: {WORKSPACE_ROOT})")
    args = ap.parse_args()

    targets = args.packages
    if not targets:
        packages = load_packages()
        print(f"Loaded {len(packages)} packages from packages.csv.")
        print("Loading build status databases...")
        status = load_build_status()
        targets = [pkg for pkg in packages if status.get(pkg) == "error"]
        print(f"Found {len(targets)} packages with ERROR status.\n")

    counts = {"synced": 0, "updated": 0, "skipped": 0, "failed": 0}
    for i, pkg in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] {pkg}...")
        try:
            counts[sync_package(pkg, args.workspace, args.push)] += 1
        except RuntimeError as e:
            print(f"  [Error] {e}")
            counts["failed"] += 1

    print(f"\nDone!  Already in sync: {counts['synced']}  Updated: {counts['updated']}  "
          f"Skipped: {counts['skipped']}  Failed: {counts['failed']}")


if __name__ == "__main__":
    main()
