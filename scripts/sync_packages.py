#!/usr/bin/env python3
"""
sync_packages.py

Fetches the current list of deprecated and voluntarily-listed packages from
the Bioconductor Help Wanted page and merges them into packages.csv.

Existing rows are preserved (including source="manual" entries).
New packages discovered on the Help Wanted page are added automatically.
Packages no longer on the Help Wanted page are NOT removed — they stay with
their original source value so manual additions and past entries are retained.

Usage:
    python sync_packages.py
"""
import csv
import re
import urllib.request
import os

HELP_WANTED_URL = "https://bioconductor.org/developers/help_wanted/"
BUILD_DATABASES = {
    "Software":       "https://bioconductor.org/checkResults/release/bioc-LATEST/BUILD_STATUS_DB.txt",
    "ExperimentData": "https://bioconductor.org/checkResults/release/data-experiment-LATEST/BUILD_STATUS_DB.txt",
    "AnnotationData": "https://bioconductor.org/checkResults/release/data-annotation-LATEST/BUILD_STATUS_DB.txt",
}

PACKAGES_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "packages.csv")
CSV_FIELDS = ["package", "type", "source"]


def fetch_html(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        return resp.read().decode("utf-8")


def fetch_help_wanted():
    """Returns {pkg: (type, source)} from the Help Wanted page."""
    html = fetch_html(HELP_WANTED_URL)
    result = {}

    # --- Deprecated (RELEASE section) ---
    release_start = html.find("<h3>RELEASE</h3>")
    if release_start == -1:
        raise ValueError("Could not find RELEASE section")
    devel_start = html.find("<h3>DEVEL</h3>", release_start)
    release_html = html[release_start: devel_start if devel_start != -1 else len(html)]

    sec_mappings = {
        "Software Packages":   "Software",
        "Experiment Packages": "ExperimentData",
        "Annotation Packages": "AnnotationData",
    }
    for sec_header, pkg_type in sec_mappings.items():
        m = re.search(rf"<h6>{sec_header}</h6>", release_html)
        if not m:
            continue
        start = m.end()
        nexts = [x.start() for x in re.finditer(r"<h6>|<h3>", release_html[start:])]
        end = start + (nexts[0] if nexts else len(release_html) - start)
        for pkg in re.findall(r'href="/packages/([^/]+)/"', release_html[start:end]):
            result[pkg] = (pkg_type, "deprecated")

    # --- Voluntarily listed ---
    v_start = html.find('id="voluntarily-listed"')
    if v_start != -1:
        v_end = html.find('id="deprecated-packages"', v_start)
        sub = html[v_start: v_end if v_end != -1 else len(html)]
        for pkg in re.findall(r'href="/packages/([^"/]+)(?:/")?"', sub):
            if pkg not in result:
                result[pkg] = (None, "voluntarily-listed")  # type resolved later

    return result


def resolve_types_from_build_db():
    """Returns {pkg: type} by scanning all three build databases."""
    pkg_type_map = {}
    for pkg_type, url in BUILD_DATABASES.items():
        try:
            content = fetch_html(url)
            for line in content.splitlines():
                m = re.match(r"^([^#]+)#", line)
                if m:
                    pkg_type_map[m.group(1).strip()] = pkg_type
        except Exception as e:
            print(f"  Warning: could not fetch {pkg_type} build DB: {e}")
    return pkg_type_map


def load_csv():
    """Returns an ordered dict of {package: {type, source}}."""
    if not os.path.exists(PACKAGES_CSV):
        return {}
    with open(PACKAGES_CSV, newline="") as f:
        reader = csv.DictReader(f)
        return {row["package"]: {"type": row["type"], "source": row["source"]} for row in reader}


def save_csv(packages):
    """Writes the packages dict to packages.csv, sorted by source then name."""
    order = {"deprecated": 0, "voluntarily-listed": 1, "manual": 2}
    rows = sorted(packages.items(), key=lambda kv: (order.get(kv[1]["source"], 9), kv[0]))
    with open(PACKAGES_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for pkg, info in rows:
            writer.writerow({"package": pkg, "type": info["type"], "source": info["source"]})


def main():
    print("Fetching Help Wanted packages...")
    help_wanted = fetch_help_wanted()
    print(f"  Found {len(help_wanted)} packages on Help Wanted page.")

    print("Resolving package types from build databases...")
    pkg_type_map = resolve_types_from_build_db()

    # Resolve types for voluntarily-listed packages whose type is unknown
    for pkg, (pkg_type, source) in help_wanted.items():
        if pkg_type is None:
            help_wanted[pkg] = (pkg_type_map.get(pkg, "Software"), source)

    print("Loading existing packages.csv...")
    existing = load_csv()

    added = []
    for pkg, (pkg_type, source) in help_wanted.items():
        if pkg not in existing:
            existing[pkg] = {"type": pkg_type, "source": source}
            added.append(pkg)
        # Never overwrite manual entries or change source of already-known packages

    save_csv(existing)

    print(f"\nDone! packages.csv updated.")
    print(f"  Total packages: {len(existing)}")
    if added:
        print(f"  Newly added: {', '.join(sorted(added))}")
    else:
        print("  No new packages added.")


if __name__ == "__main__":
    main()
