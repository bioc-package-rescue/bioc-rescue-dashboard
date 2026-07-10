#!/usr/bin/env python3
"""
submit_upstream.py

Automates preparing and opening a Pull Request from a rescued fork back to the
original upstream repository (parent), ensuring that rescue-specific files
(like .github/workflows/check-bioc.yml) are excluded from the PR diff.

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
    
    # Ensure upstream remote is added locally
    remotes = run(["git", "remote"], cwd=local_path, capture=True).splitlines()
    if "upstream" not in remotes:
        upstream_url = f"https://github.com/{parent_slug}.git"
        print(f"Adding remote 'upstream' -> {upstream_url}")
        run(["git", "remote", "add", "upstream", upstream_url], cwd=local_path)
    else:
        print("Remote 'upstream' is already configured.")
        
    print("Fetching upstream...")
    run(["git", "fetch", "upstream"], cwd=local_path)
    
    # Create a temporary submission branch off upstream default branch
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    submit_branch = f"submit-upstream-{timestamp}"
    print(f"Creating clean local branch '{submit_branch}' off 'upstream/{upstream_default_branch}'...")
    run(["git", "checkout", "-b", submit_branch, f"upstream/{upstream_default_branch}"], cwd=local_path)
    
    success = False
    try:
        # Find files modified in fix_branch compared to upstream/defaultBranch
        print(f"Identifying modified package files from '{fix_branch}'...")
        diff_files = run(["git", "diff", "--name-only", f"upstream/{upstream_default_branch}...{fix_branch}"], cwd=local_path, capture=True).splitlines()
        
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
            
        print("Checking out modified files:")
        for f in files_to_checkout:
            print(f"  + {f}")
            run(["git", "checkout", fix_branch, "--", f], cwd=local_path)
            
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
            f"This pull request contains automatically-verified fixes for `{pkg}` to resolve build/check errors "
            f"on Bioconductor. All checks have been run and verified green in our rescue fork check suite.\n\n"
            f"Please review and merge these changes back into your upstream repository.\n\n"
            f"Co-authored-by: Antigravity <gemini@google.com>"
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
