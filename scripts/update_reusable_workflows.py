import os
import subprocess
import re
import urllib.request
import time

HELP_WANTED_URL = "https://bioconductor.org/developers/help_wanted/"

BUILD_DATABASES = {
    "Software": "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt"
}

parsed_status_db = {}
workspace_root = "/Users/Levi/git/bioc-package-rescue"

workflow_content = """name: R-CMD-check-bioc

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

def fetch_help_wanted_packages():
    print(f"Fetching Help Wanted page from {HELP_WANTED_URL}...")
    try:
        req = urllib.request.Request(HELP_WANTED_URL, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            html = resp.read().decode('utf-8')
            
        release_start = html.find('<h3>RELEASE</h3>')
        if release_start == -1:
            raise ValueError("Could not find RELEASE section on Help Wanted page")
            
        devel_start = html.find('<h3>DEVEL</h3>', release_start)
        if devel_start == -1:
            devel_start = len(html)
            
        release_html = html[release_start:devel_start]
        
        sec_mappings = {
            "Software Packages": "Software",
            "Experiment Packages": "ExperimentData",
            "Annotation Packages": "AnnotationData"
        }
        
        packages = []
        
        # Parse Deprecated section
        for sec_header, pkg_type in sec_mappings.items():
            header_pattern = rf"<h6>{sec_header}</h6>"
            header_match = re.search(header_pattern, release_html)
            if not header_match:
                continue
                
            start_idx = header_match.end()
            next_headers = [m.start() for m in re.finditer(r"<h6>|<h3>", release_html[start_idx:])]
            end_idx = start_idx + (next_headers[0] if next_headers else len(release_html) - start_idx)
            
            sub_html = release_html[start_idx:end_idx]
            links = re.findall(r'href="/packages/([^/]+)/"', sub_html)
            for pkg in sorted(list(set(links))):
                packages.append((pkg, pkg_type))
                
        # Parse Voluntarily Listed section
        start_pos = html.find('id="voluntarily-listed"')
        if start_pos != -1:
            end_pos = html.find('id="deprecated-packages"', start_pos)
            if end_pos == -1:
                end_pos = len(html)
            sub_html = html[start_pos:end_pos]
            
            v_links = re.findall(r'href="/packages/([^"/]+)(?:/")?"', sub_html)
            for pkg in sorted(list(set(v_links))):
                packages.append((pkg, None))
                
        return packages
    except Exception as e:
        print(f"Error fetching/parsing Help Wanted page: {e}")
        return []

def load_build_databases():
    for pkg_type, url in BUILD_DATABASES.items():
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as resp:
                content = resp.read().decode('utf-8')
            for line in content.split('\n'):
                line = line.strip()
                if not line or line == '---' or line.startswith('Title:') or line.startswith('Description:') or line.startswith('Source:'):
                    continue
                match = re.match(r"^([^#]+)#([^#]+)#([^:]+):\s*(.*)$", line)
                if match:
                    pkg, builder, stage, status = match.groups()
                    pkg = pkg.strip()
                    builder = builder.strip()
                    stage = stage.strip()
                    status = status.strip()
                    
                    if pkg not in parsed_status_db:
                        parsed_status_db[pkg] = {}
                    if builder not in parsed_status_db[pkg]:
                        parsed_status_db[pkg][builder] = {}
                    parsed_status_db[pkg][builder][stage] = status
        except Exception as e:
            print(f"Error loading {pkg_type} DB: {e}")

def get_package_build_status(package_name):
    if package_name not in parsed_status_db:
        return "NA"
    builders = parsed_status_db[package_name]
    has_error = False
    for builder, stages in builders.items():
        for stage, status in stages.items():
            if status in ["ERROR", "TIMEOUT"]:
                has_error = True
    return "error" if has_error else "OK"

def run_command(args, cwd=None):
    env = os.environ.copy()
    if 'GITHUB_TOKEN' in env:
        del env['GITHUB_TOKEN']
    result = subprocess.run(args, capture_output=True, text=True, cwd=cwd, env=env)
    return result

def main():
    print("Loading build databases...")
    load_build_databases()
    
    # Fetch all help wanted packages
    packages = fetch_help_wanted_packages()
    
    # Filter and find all packages with status == error
    error_packages = []
    for pkg, pkg_type in packages:
        status = get_package_build_status(pkg)
        if status == "error":
            error_packages.append(pkg)
                
    total = len(error_packages)
    print(f"Found {total} packages with ERROR status.")
    
    success_count = 0
    failed_count = 0
    
    for idx, pkg in enumerate(error_packages, 1):
        local_path = os.path.join(workspace_root, pkg)
        print(f"\n[{idx}/{total}] Processing {pkg}...")
        
        if not os.path.exists(local_path):
            print(f"  [Error] Local repository path {local_path} does not exist.")
            failed_count += 1
            continue
            
        # Get active branch name
        branch_res = run_command(["git", "branch", "--show-current"], cwd=local_path)
        if branch_res.returncode != 0:
            print(f"  [Error] Could not get active branch for {pkg}: {branch_res.stderr.strip()}")
            failed_count += 1
            continue
        active_branch = branch_res.stdout.strip()
        print(f"  Active branch: {active_branch}")
        
        # Write workflow file
        workflows_dir = os.path.join(local_path, ".github", "workflows")
        os.makedirs(workflows_dir, exist_ok=True)
        workflow_path = os.path.join(workflows_dir, "check-bioc.yml")
        
        with open(workflow_path, "w") as f:
            f.write(workflow_content)
            
        # Stage, commit, and push
        run_command(["git", "add", ".github/workflows/check-bioc.yml"], cwd=local_path)
        
        # Check if anything changed
        status_res = run_command(["git", "status", "--porcelain"], cwd=local_path)
        if ".github/workflows/check-bioc.yml" in status_res.stdout:
            commit_res = run_command(["git", "commit", "-m", "Use centralized check-bioc GHA workflow", "--author=Antigravity <gemini@google.com>"], cwd=local_path)
            if commit_res.returncode != 0:
                print(f"  [Error] Failed to commit for {pkg}: {commit_res.stderr.strip()}")
                failed_count += 1
                continue
            print(f"  Committed workflow.")
        else:
            print(f"  Workflow already committed or no changes.")
            
        # Push to remote
        push_res = run_command(["git", "push", "origin", active_branch], cwd=local_path)
        if push_res.returncode != 0:
            print(f"  [Error] Failed to push {pkg} branch {active_branch}: {push_res.stderr.strip()}")
            failed_count += 1
            continue
            
        print(f"  Successfully updated and pushed workflow for {pkg}.")
        success_count += 1
        time.sleep(0.5) # simple rate limit politeness
        
    print(f"\nBatch processing complete!")
    print(f"  Total error packages: {total}")
    print(f"  Successfully updated and pushed: {success_count}")
    print(f"  Failed: {failed_count}")

if __name__ == "__main__":
    main()
