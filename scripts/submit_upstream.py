#!/usr/bin/env python3
"""
submit_upstream.py

Automates preparing and opening a Pull Request from a rescued fork back to the
original upstream repository (parent), ensuring that rescue-specific files
(like .github/workflows/check-bioc.yml) are excluded from the PR diff.

Rescue repos are synced with git.bioconductor.org devel before any fixes
(see sync_repos.py), so the rescue changes are the diff between
upstream/devel and the fix branch. Only that diff is applied to the GitHub
author's default branch; commits that came from git.bioconductor.org are not
proposed as fixes.

Remotes in the local checkout:
  origin    the rescue fork (bioc-package-rescue/<pkg>)
  upstream  git.bioconductor.org: fetches over HTTPS
            (https://git.bioconductor.org/packages/<pkg>), pushes over SSH
            (git@git.bioconductor.org:packages/<pkg>)
  author    the package author's GitHub repository, which the rescue repo was
            forked from (its GitHub "parent")

Usage:
  python scripts/submit_upstream.py <package_name> <fix_branch_name>
"""

import sys
import os
import json
import subprocess
from datetime import datetime

WORKSPACE_ROOT = "/Users/Levi/git/bioc-package-rescue"

def run(cmd, cwd=None, capture=False):
    """Helper to run shell commands."""
    if capture:
        result = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return result.stdout.strip()
    else:
        subprocess.run(cmd, cwd=cwd, check=True)

def main():
    if len(sys.argv) < 3:
        print("Usage: python scripts/submit_upstream.py <package_name> <fix_branch_name>")
        sys.exit(1)
        
    pkg = sys.argv[1]
    fix_branch = sys.argv[2]
    local_path = os.path.join(WORKSPACE_ROOT, pkg)
    
    if not os.path.isdir(local_path):
        print(f"Error: Local directory not found for package {pkg} at {local_path}")
        sys.exit(1)
        
    print(f"Resolving parent repository for bioc-package-rescue/{pkg}...")
    try:
        repo_info_json = run(["gh", "repo", "view", f"bioc-package-rescue/{pkg}", "--json", "parent,defaultBranchRef"], capture=True)
        repo_info = json.loads(repo_info_json)
    except Exception as e:
        print(f"Error calling gh CLI: {e}")
        sys.exit(1)
        
    parent = repo_info.get("parent")
    if not parent:
        print(f"Package {pkg} is not a fork on GitHub (no parent repository). Cannot submit upstream PR.")
        sys.exit(1)
        
    parent_owner = parent["owner"]["login"]
    parent_name = parent["name"]
    parent_slug = f"{parent_owner}/{parent_name}"
    print(f"Found parent upstream repository: {parent_slug}")
    
    # Get upstream default branch
    try:
        parent_info_json = run(["gh", "repo", "view", parent_slug, "--json", "defaultBranchRef"], capture=True)
        parent_info = json.loads(parent_info_json)
        upstream_default_branch = parent_info["defaultBranchRef"]["name"]
    except Exception:
        upstream_default_branch = "master" # standard fallback
        
    print(f"Upstream default branch: {upstream_default_branch}")
    
    # Ensure the author's GitHub repo (the fork's parent) is added locally
    remotes = run(["git", "remote"], cwd=local_path, capture=True).splitlines()
    if "author" not in remotes:
        author_url = f"https://github.com/{parent_slug}.git"
        print(f"Adding remote 'author' -> {author_url}")
        run(["git", "remote", "add", "author", author_url], cwd=local_path)
    else:
        print("Remote 'author' is already configured.")
        
    print("Fetching author...")
    run(["git", "fetch", "author"], cwd=local_path)

    # The rescue changes are measured from the Bioconductor devel branch that
    # the rescue repo was synced to, not from the author's GitHub repo, which may
    # lag behind or differ from git.bioconductor.org.
    bioc_url = f"https://git.bioconductor.org/packages/{pkg}"
    bioc_push_url = f"git@git.bioconductor.org:packages/{pkg}"
    if "upstream" not in remotes:
        print(f"Adding remote 'upstream' -> {bioc_url} (push: {bioc_push_url})")
        run(["git", "remote", "add", "upstream", bioc_url], cwd=local_path)
    else:
        current = run(["git", "remote", "get-url", "upstream"], cwd=local_path, capture=True)
        if current.rstrip("/").removesuffix(".git") != bioc_url:
            print(f"Error: remote 'upstream' points to {current}, not {bioc_url}. "
                  f"Fix it with: git remote set-url upstream {bioc_url}")
            sys.exit(1)
    current_push = run(["git", "remote", "get-url", "--push", "upstream"], cwd=local_path, capture=True)
    if current_push != bioc_push_url:
        run(["git", "remote", "set-url", "--push", "upstream", bioc_push_url], cwd=local_path)
    print("Fetching upstream devel...")
    run(["git", "fetch", "upstream", "devel"], cwd=local_path)
    base_ref = "upstream/devel"
    try:
        run(["git", "merge-base", "--is-ancestor", base_ref, fix_branch], cwd=local_path)
    except subprocess.CalledProcessError:
        print(f"Error: '{fix_branch}' does not contain {base_ref}. "
              f"Run scripts/sync_repos.py {pkg} and rebase the fix branch first.")
        sys.exit(1)
    
    # Create a temporary submission branch off the author's default branch
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    submit_branch = f"submit-upstream-{timestamp}"
    print(f"Creating clean local branch '{submit_branch}' off 'author/{upstream_default_branch}'...")
    run(["git", "checkout", "-b", submit_branch, f"author/{upstream_default_branch}"], cwd=local_path)
    
    success = False
    try:
        # Find files modified in fix_branch compared to upstream/devel
        print(f"Identifying package files changed between {base_ref} and '{fix_branch}'...")
        diff_files = run(["git", "diff", "--name-only", f"{base_ref}..{fix_branch}"], cwd=local_path, capture=True).splitlines()
        
        # Filter out rescue-specific configuration files
        files_to_checkout = []
        for f in diff_files:
            if f.startswith(".github/") or f == ".gitignore" or f == ".gitattributes":
                print(f"  Skipping rescue-specific file: {f}")
            else:
                files_to_checkout.append(f)
                
        if not files_to_checkout:
            print("No package-specific code changes found to submit.")
            return
            
        print("Applying rescue changes to:")
        for f in files_to_checkout:
            print(f"  + {f}")
        patch = subprocess.run(
            ["git", "diff", "--binary", f"{base_ref}..{fix_branch}", "--", *files_to_checkout],
            cwd=local_path, stdout=subprocess.PIPE, check=True).stdout
        applied = subprocess.run(["git", "apply", "--3way", "--index"],
                                 cwd=local_path, input=patch)
        if applied.returncode != 0:
            print(f"Error: the rescue changes do not apply cleanly to author/{upstream_default_branch}. "
                  "Resolve the conflicts by hand.")
            run(["git", "reset", "-q", "--hard"], cwd=local_path)
            return
        if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=local_path).returncode == 0:
            print(f"author/{upstream_default_branch} already has these changes; nothing to submit.")
            success = True
            return
            
        # Commit the changes
        commit_msg = f"Fix R CMD check errors for {pkg}\n\nCo-authored-by: Antigravity <gemini@google.com>"
        print("Committing changes...")
        run(["git", "commit", "-m", commit_msg], cwd=local_path)
        
        # Push to origin (fork)
        print(f"Pushing submit branch to origin/{submit_branch}...")
        run(["git", "push", "origin", f"{submit_branch}:{submit_branch}"], cwd=local_path)
        
        # Create GitHub Pull Request back to upstream parent
        pr_title = f"Fix R CMD check build/check errors for {pkg}"
        pr_body = (
            f"This pull request contains AI-generated fixes for `{pkg}` that have been verified to resolve current Bioconductor build/check errors.\n\n"
            f"Please review these changes, and if they look OK, merge them back into your upstream git.bioconductor.org repository.\n\n"
        )
        print(f"Creating pull request to {parent_slug}...")
        pr_url = run([
            "gh", "pr", "create",
            "--repo", parent_slug,
            "--head", f"bioc-package-rescue:{submit_branch}",
            "--base", upstream_default_branch,
            "--title", pr_title,
            "--body", pr_body
        ], cwd=local_path, capture=True)
        
        print("\n" + "="*60)
        print(f"SUCCESS! Pull request successfully created at:\n{pr_url}")
        print("="*60 + "\n")
        success = True
        
    finally:
        # Cleanup: switch back to the original fix branch and delete local submission branch
        print(f"Cleaning up: switching back to '{fix_branch}'...")
        run(["git", "checkout", fix_branch], cwd=local_path)
        run(["git", "branch", "-D", submit_branch], cwd=local_path)
        if not success:
            print("Failed to complete upstream PR submission.")

if __name__ == "__main__":
    main()
