import os
import subprocess
import re
import urllib.request

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
    skipped_count = 0
    failed_count = 0
    
    for idx, pkg in enumerate(error_packages, 1):
        local_path = os.path.join(workspace_root, pkg)
        print(f"\n[{idx}/{total}] Processing {pkg}...")
        
        # Check if local directory exists
        if os.path.exists(local_path):
            print(f"  [Skip] Local directory {local_path} already exists.")
            skipped_count += 1
            continue
            
        github_ssh_url = f"git@github.com:bioc-package-rescue/{pkg}.git"
        print(f"  Cloning {github_ssh_url} to {local_path}...")
        
        res = subprocess.run(["git", "clone", github_ssh_url, local_path], capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  Successfully cloned {pkg}.")
            success_count += 1
        else:
            print(f"  Error cloning {pkg}: {res.stderr.strip()}")
            failed_count += 1
            
    print(f"\nLocal cloning complete!")
    print(f"  Total packages: {total}")
    print(f"  Successfully cloned: {success_count}")
    print(f"  Skipped (already exist): {skipped_count}")
    print(f"  Failed: {failed_count}")

if __name__ == "__main__":
    main()
