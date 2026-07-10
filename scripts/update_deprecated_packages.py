import csv
import urllib.request
import re
import datetime
import os

PACKAGES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "packages.csv")

BUILD_DATABASES = {
    "Software": "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt"
}

STATS_DATABASES = {
    "Software": "https://bioconductor.org/packages/stats/bioc/bioc_pkg_stats.tab",
    "ExperimentData": "https://bioconductor.org/packages/stats/data-experiment/bioc_pkg_stats.tab",
    "AnnotationData": "https://bioconductor.org/packages/stats/data-annotation/bioc_pkg_stats.tab"
}

# Parsed status database in memory:
# { pkg_name: { builder: { stage: status } } }
parsed_status_db = {}

# Parsed download stats in memory:
# { pkg_name: distinct_ips }
parsed_stats_db = {}

# Package name -> package type mapping derived from build databases
package_types = {}

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

def get_build_page_url(package_name, pkg_type):
    if pkg_type == "Software":
        segment = "bioc-LATEST"
    elif pkg_type == "ExperimentData":
        segment = "data-experiment-LATEST"
    elif pkg_type == "AnnotationData":
        segment = "data-annotation-LATEST"
    else:
        segment = "bioc-LATEST"
    return f"https://bioconductor.org/checkResults/release/{segment}/{package_name}/"

def load_build_databases():
    print("Loading build status databases...")
    for pkg_type, url in BUILD_DATABASES.items():
        try:
            print(f"  Fetching {pkg_type} DB from {url}...")
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
                    
                    if pkg not in package_types:
                        package_types[pkg] = pkg_type
            print(f"  Successfully parsed {pkg_type} DB.")
        except Exception as e:
            print(f"  Warning: Failed to load build status DB for {pkg_type} from {url}: {e}")

def load_download_stats(target_year="2025"):
    print(f"Loading download statistics for {target_year}...")
    for pkg_type, url in STATS_DATABASES.items():
        try:
            print(f"  Fetching stats from {url}...")
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as resp:
                content = resp.read().decode('utf-8')
            for line in content.split('\n'):
                line = line.strip()
                if not line or line.startswith("Package\t"):
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    pkg, year, month, distinct_ips, downloads = parts[0], parts[1], parts[2], parts[3], parts[4]
                    if year == target_year and month == "all":
                        parsed_stats_db[pkg] = int(distinct_ips)
            print(f"  Successfully parsed stats for {pkg_type}.")
        except Exception as e:
            print(f"  Warning: Failed to load download stats for {pkg_type} from {url}: {e}")

def get_repo_url(package_name, pkg_type):
    landing_url = get_landing_page_url(package_name, pkg_type)
    fallback_url = f"https://git.bioconductor.org/packages/{package_name}"
    try:
        req = urllib.request.Request(landing_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            html = resp.read().decode('utf-8')
        
        # Search for github.com links in HTML
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
    except Exception as e:
        return fallback_url

def get_package_build_status(package_name, pkg_type):
    if package_name not in parsed_status_db:
        return "NA"
    
    builders = parsed_status_db[package_name]
    
    has_error = False
    has_warning = False
    
    for builder, stages in builders.items():
        for stage, status in stages.items():
            if status in ["ERROR", "TIMEOUT"]:
                has_error = True
            elif status == "WARNINGS":
                has_warning = True
                
    if has_error:
        return "error"
    if has_warning:
        return "warning"
        
    # If all OK/skipped, check detailed checksrc logs for NOTES
    has_note = False
    for builder, stages in builders.items():
        check_status = stages.get("checksrc")
        if check_status == "OK":
            if pkg_type == "Software":
                segment = "bioc-LATEST"
            elif pkg_type == "ExperimentData":
                segment = "data-experiment-LATEST"
            elif pkg_type == "AnnotationData":
                segment = "data-annotation-LATEST"
            else:
                segment = "bioc-LATEST"
                
            log_url = f"https://bioconductor.org/checkResults/release/{segment}/{package_name}/{builder}-checksrc.html"
            try:
                req = urllib.request.Request(log_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req) as resp:
                    log_html = resp.read().decode('utf-8')
                status_matches = re.findall(r'Status:\s*([^\n<]+)', log_html)
                for status_text in status_matches:
                    if 'NOTE' in status_text:
                        has_note = True
                        break
            except Exception:
                pass
        if has_note:
            break
            
    if has_note:
        return "note"
    return "OK"

def load_packages_csv():
    """Read packages.csv and return two grouped dicts:
    deprecated_packages and voluntarily_packages, both {type: [pkg, ...]}.
    Also includes 'manual' entries in the voluntarily bucket."""
    deprecated = {"Software": [], "ExperimentData": [], "AnnotationData": []}
    voluntary  = {"Software": [], "ExperimentData": [], "AnnotationData": []}
    with open(PACKAGES_CSV, newline="") as f:
        for row in csv.DictReader(f):
            pkg  = row["package"]
            ptype = row["type"] if row["type"] in deprecated else "Software"
            src   = row["source"]
            if src == "deprecated":
                deprecated[ptype].append(pkg)
            else:  # voluntarily-listed or manual
                voluntary[ptype].append(pkg)
    return deprecated, voluntary


def generate_markdown_table(packages_dict):
    table_lines = []
    table_lines.append("| Package | Type | Build Page | Downloads (2025 IPs) | Repository | Bioconductor Build Status | Rescue Status | Rescue release | Rescue devel |")
    table_lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    
    total_packages = sum(len(pkgs) for pkgs in packages_dict.values())
    processed_count = 0
    
    for pkg_type, package_list in packages_dict.items():
        print(f"  Processing {pkg_type} list...")
        for pkg in sorted(package_list):
            processed_count += 1
            
            landing_url = get_landing_page_url(pkg, pkg_type)
            build_page_url = get_build_page_url(pkg, pkg_type)
            repo_url = get_repo_url(pkg, pkg_type)
            build_status = get_package_build_status(pkg, pkg_type)
            
            # Get 2025 download statistics
            distinct_ips = parsed_stats_db.get(pkg, "NA")
            
            # Construct build status badge link or use shields.io
            if pkg_type == "Software":
                build_badge_url = f"https://www.bioconductor.org/shields/build/release/bioc/{pkg}.svg"
                build_badge_markdown = f"![Build]({build_badge_url})"
            elif pkg_type == "ExperimentData":
                build_badge_url = f"https://www.bioconductor.org/shields/build/release/data-experiment/{pkg}.svg"
                build_badge_markdown = f"![Build]({build_badge_url})"
            else:
                # AnnotationData has no build status badge, use shields.io
                color_map = {
                    "OK": "brightgreen",
                    "note": "blue",
                    "warning": "orange",
                    "error": "red",
                    "NA": "lightgrey"
                }
                color = color_map.get(build_status, "lightgrey")
                label = build_status.upper() if build_status != "note" else "NOTE"
                build_badge_url = f"https://img.shields.io/badge/build-{label}-{color}"
                build_badge_markdown = f"![Build]({build_badge_url})"
                
            # Build Row
            package_link = f"[{pkg}]({landing_url})"
            build_page_link = f"[Build Page]({build_page_url})"
            repo_link = f"[Repo]({repo_url})"
            
            # Rescue Status is set to the organization repo if build status is ERROR
            if build_status == "error":
                rescue_status = f"[Rescue Repo](https://github.com/bioc-package-rescue/{pkg})"
                rescue_release = f"[![Release](https://github.com/bioc-package-rescue/{pkg}/actions/workflows/check-bioc.yml/badge.svg?branch=RELEASE_3_23)](https://github.com/bioc-package-rescue/{pkg}/actions/workflows/check-bioc.yml?query=branch%3ARELEASE_3_23)"
                rescue_devel = f"[![Devel](https://github.com/bioc-package-rescue/{pkg}/actions/workflows/check-bioc.yml/badge.svg?branch=devel)](https://github.com/bioc-package-rescue/{pkg}/actions/workflows/check-bioc.yml?query=branch%3Adevel)"
            else:
                rescue_status = "NA"
                rescue_release = "NA"
                rescue_devel = "NA"
            
            row = f"| {package_link} | {pkg_type} | {build_page_link} | {distinct_ips} | {repo_link} | {build_badge_markdown} | {rescue_status} | {rescue_release} | {rescue_devel} |"
            table_lines.append(row)
            
    return "\n".join(table_lines)

def main():
    print("Starting deprecated and voluntarily listed packages table generation...")
    print(f"Reading package list from {PACKAGES_CSV}...")
    deprecated_packages, voluntarily_packages = load_packages_csv()

    load_build_databases()
    load_download_stats(target_year="2025")

    markdown_lines = []
    markdown_lines.append("# Bioconductor Package Rescue Dashboard")

    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    markdown_lines.append(f"*Last updated: {now_str}*")
    markdown_lines.append("")

    markdown_lines.append("## Deprecated Packages")
    markdown_lines.append("These packages are scheduled for deprecation in Bioconductor.")
    markdown_lines.append("")
    print("\nGenerating Deprecated Packages Table...")
    markdown_lines.append(generate_markdown_table(deprecated_packages))
    markdown_lines.append("")

    markdown_lines.append("## Voluntarily Listed Packages")
    markdown_lines.append("These packages are voluntarily listed by their maintainers as needing transfer or assistance.")
    markdown_lines.append("")
    print("\nGenerating Voluntarily Listed Packages Table...")
    markdown_lines.append(generate_markdown_table(voluntarily_packages))

    output_filepath = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "README.md")
    )
    with open(output_filepath, "w") as f:
        f.write("\n".join(markdown_lines) + "\n")

    print(f"\nDone! Dashboard successfully written to {output_filepath}")

if __name__ == "__main__":
    main()
