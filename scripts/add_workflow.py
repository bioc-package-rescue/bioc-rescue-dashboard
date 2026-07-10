import os
import subprocess
import re
import urllib.request
import time

PACKAGES = {
    "Software": [
        "BubbleTree", "CuratedAtlasQueryR", "IONiseR", "hmdbQuery", "rols", "scviR", "SGCP",
        "APAlyzer", "ballgown", "bamsignals", "barcodetrackR", "basecallQC", "biobroom",
        "biodbChebi", "BiRewire", "BPRMeth", "CelliD", "ChIPQC", "ccrepe", "cummeRbund",
        "debCAM", "DeconRNASeq", "geneXtendeR", "GEOexplorer", "GNET2", "granulator",
        "hca", "IMAS", "Melissa", "MetaNeighbor", "MethReg", "mfa", "microSTASIS",
        "MineICA", "motifcounter", "MSPrep", "nearBynding", "netprioR", "normr",
        "Organism.dplyr", "partCNV", "RcisTarget", "receptLoss", "RgnTX", "RiboProfiling",
        "RTCGA", "shiny.gosling", "soGGi", "SigFuge", "spatzie", "SQLDataFrame",
        "supersigs", "tLOH"
    ],
    "ExperimentData": [
        "curatedBreastData", "Fletcher2013b", "rRDPData"
    ],
    "AnnotationData": [
        "hpAnnot"
    ]
}

BUILD_DATABASES = {
    "Software": "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt"
}

parsed_status_db = {}
workspace_root = "/Users/Levi/git/bioc-package-rescue"

workflow_content = """on:
  push:
    branches:
      - main
      - devel
  pull_request:

name: R-CMD-check-bioc

jobs:
  check:
    runs-on: ubuntu-latest
    container: ${{ matrix.config.image }}

    name: ${{ matrix.config.bioc }}

    strategy:
      fail-fast: false
      matrix:
        config:
          - { bioc: 'release', image: 'bioconductor/bioconductor_docker:RELEASE_3_21' }
          - { bioc: 'devel',   image: 'bioconductor/bioconductor_docker:devel'         }

    env:
      R_REMOTES_NO_ERRORS_FROM_WARNINGS: true
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      NOT_CRAN: true

    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::BiocCheck, any::rcmdcheck
          needs: check
          cache: true

      - name: R CMD CHECK
        uses: r-lib/actions/check-r-package@v2
        with:
          args: 'c("--no-manual", "--no-vignettes", "--timings")'
          error-on: '"error"'

      - name: BiocCheck
        shell: Rscript {0}
        run: |
          BiocCheck::BiocCheck(
            package        = ".",
            `new-package`  = FALSE,
            `no-check-bioc-help` = TRUE
          )
"""

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
    
    # Find all packages with status == error
    error_packages = []
    for pkg_type, package_list in PACKAGES.items():
        for pkg in sorted(package_list):
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
        
        # Check if anything changed (in case it already exists/committed)
        status_res = run_command(["git", "status", "--porcelain"], cwd=local_path)
        if ".github/workflows/check-bioc.yml" in status_res.stdout:
            commit_res = run_command(["git", "commit", "-m", "Add GitHub Actions workflow for Bioconductor checks", "--author=Antigravity <gemini@google.com>"], cwd=local_path)
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
            
        print(f"  Successfully added and pushed workflow for {pkg}.")
        success_count += 1
        time.sleep(0.5) # simple rate limit politeness
        
    print(f"\nBatch processing complete!")
    print(f"  Total error packages: {total}")
    print(f"  Successfully updated and pushed: {success_count}")
    print(f"  Failed: {failed_count}")

if __name__ == "__main__":
    main()
