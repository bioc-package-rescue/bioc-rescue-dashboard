import os
import subprocess
import re
import urllib.request
import time
import shutil

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
temp_dir_base = "/Users/Levi/git/bioc-package-rescue/bioc-rescue-dashboard/scripts/temp_clones"

def get_landing_page_url(package_name, pkg_type):
    if pkg_type == "Software":
        segment = "bioc"
    elif pkg_type == "ExperimentData":
        segment = "data/experiment"
    elif pkg_type == "AnnotationData":
        segment = "data/annotation"
    else:
        segment = "bioc"
    return f"https://bioconductor.org/packages/release/{segment}/html/{package_name}.html"

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

def get_repo_url(package_name, pkg_type):
    landing_url = get_landing_page_url(package_name, pkg_type)
    fallback_url = f"https://git.bioconductor.org/packages/{package_name}"
    try:
        req = urllib.request.Request(landing_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            html = resp.read().decode('utf-8')
        github_links = re.findall(r'href="(https?://github\.com/[^"]+)"', html)
        if github_links:
            github_url = github_links[0]
            match = re.match(r'(https?://(?:www\.)?github\.com/[^/]+/[^/#?]+)', github_url)
            if match:
                clean_url = match.group(1).rstrip('/')
                if clean_url.endswith('.git'):
                    clean_url = clean_url[:-4]
                return clean_url
        return fallback_url
    except Exception:
        return fallback_url

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
    # Copy env and delete GITHUB_TOKEN to fallback to keyring
    env = os.environ.copy()
    if 'GITHUB_TOKEN' in env:
        del env['GITHUB_TOKEN']
    
    result = subprocess.run(args, capture_output=True, text=True, cwd=cwd, env=env)
    return result

def check_repo_exists(pkg_name):
    res = run_command(["gh", "repo", "view", f"bioc-package-rescue/{pkg_name}"])
    return res.returncode == 0

def fork_github_repo(repo_url, pkg_name):
    print(f"  Forking {repo_url} to bioc-package-rescue...")
    res = run_command(["gh", "repo", "fork", repo_url, "--org", "bioc-package-rescue", "--clone=false"])
    if res.returncode == 0:
        print(f"  Successfully forked {pkg_name}.")
        return True
    else:
        print(f"  Error forking {pkg_name}: {res.stderr.strip()}")
        return False

def clone_and_push_bioc_repo(pkg_name):
    print(f"  Cloning from git.bioconductor.org...")
    # Create temp dir
    os.makedirs(temp_dir_base, exist_ok=True)
    temp_pkg_dir = os.path.join(temp_dir_base, pkg_name)
    if os.path.exists(temp_pkg_dir):
        shutil.rmtree(temp_pkg_dir)
        
    bioc_git_url = f"https://git.bioconductor.org/packages/{pkg_name}"
    
    # Git clone
    res_clone = run_command(["git", "clone", bioc_git_url, temp_pkg_dir])
    if res_clone.returncode != 0:
        print(f"  Error cloning {pkg_name} from Bioconductor: {res_clone.stderr.strip()}")
        return False
        
    print(f"  Creating public repository bioc-package-rescue/{pkg_name}...")
    res_create = run_command(["gh", "repo", "create", f"bioc-package-rescue/{pkg_name}", "--public"])
    if res_create.returncode != 0:
        # Check if it failed because it exists (race condition)
        if check_repo_exists(pkg_name):
            print(f"  Repository bioc-package-rescue/{pkg_name} already exists.")
            shutil.rmtree(temp_pkg_dir)
            return True
        print(f"  Error creating repository for {pkg_name}: {res_create.stderr.strip()}")
        shutil.rmtree(temp_pkg_dir)
        return False
        
    print(f"  Configuring remote and pushing to GitHub...")
    # Add remote
    github_ssh_url = f"git@github.com:bioc-package-rescue/{pkg_name}.git"
    run_command(["git", "remote", "rename", "origin", "upstream"], cwd=temp_pkg_dir)
    run_command(["git", "remote", "add", "origin", github_ssh_url], cwd=temp_pkg_dir)
    
    # Push all branches and tags
    res_push = run_command(["git", "push", "-u", "origin", "--all"], cwd=temp_pkg_dir)
    if res_push.returncode != 0:
        print(f"  Warning during git push --all: {res_push.stderr.strip()}")
        
    res_push_tags = run_command(["git", "push", "-u", "origin", "--tags"], cwd=temp_pkg_dir)
    if res_push_tags.returncode != 0:
        print(f"  Warning during git push --tags: {res_push_tags.stderr.strip()}")
        
    # Clean up
    shutil.rmtree(temp_pkg_dir)
    print(f"  Successfully cloned, created, and pushed {pkg_name}.")
    return True

def main():
    print("Loading build databases...")
    load_build_databases()
    
    # Find all packages with status == error
    error_packages = []
    for pkg_type, package_list in PACKAGES.items():
        for pkg in sorted(package_list):
            status = get_package_build_status(pkg)
            if status == "error":
                repo = get_repo_url(pkg, pkg_type)
                error_packages.append((pkg, repo))
                
    total = len(error_packages)
    print(f"Found {total} packages with ERROR status.")
    
    success_count = 0
    skipped_count = 0
    failed_count = 0
    
    for idx, (pkg, repo) in enumerate(error_packages, 1):
        print(f"\n[{idx}/{total}] Processing {pkg}...")
        
        # Check if already exists in bioc-package-rescue
        if check_repo_exists(pkg):
            print(f"  [Skip] bioc-package-rescue/{pkg} already exists on GitHub.")
            skipped_count += 1
            continue
            
        is_github = "github.com" in repo
        if is_github:
            success = fork_github_repo(repo, pkg)
        else:
            success = clone_and_push_bioc_repo(pkg)
            
        if success:
            success_count += 1
        else:
            failed_count += 1
            
        # Polite sleep to avoid hitting secondary rate limits on GitHub
        time.sleep(2)
        
    print(f"\nProcessing Complete!")
    print(f"  Total: {total}")
    print(f"  Forked/Created: {success_count}")
    print(f"  Skipped (already exists): {skipped_count}")
    print(f"  Failed: {failed_count}")
    
    # Clean up temp clones directory if exists
    if os.path.exists(temp_dir_base):
        try:
            shutil.rmtree(temp_dir_base)
        except Exception:
            pass

if __name__ == "__main__":
    main()
