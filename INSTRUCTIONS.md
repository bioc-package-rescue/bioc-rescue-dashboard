# Bioconductor Package Rescue Instructions

This document provides instructions for agents and developers on how to maintain the `bioc-rescue-dashboard`, manage rescue repositories, and update the centralized GitHub Actions check workflow.

---

## Package List: `packages.csv`

**`packages.csv`** (in the root of this repo) is the single source of truth for all tracked packages. Every script reads from it — nothing is hard-coded elsewhere.

### Format

```csv
package,type,source
BPRMeth,Software,deprecated
VariantAnnotation,Software,voluntarily-listed
myCustomPkg,Software,manual
```

| Column    | Values                                          |
|-----------|-------------------------------------------------|
| `package` | Bioconductor package name                       |
| `type`    | `Software`, `ExperimentData`, or `AnnotationData` |
| `source`  | `deprecated`, `voluntarily-listed`, or `manual` |

### Adding a Package Manually

Edit `packages.csv` directly and add a row with `source=manual`. No other files need to be changed — the next dashboard run will pick it up automatically.

---

## Typical Workflow

Run in order when setting up fresh or after Bioconductor publishes updates:

```bash
# 1. Sync package list from Bioconductor Help Wanted page
python scripts/sync_packages.py

# 2. Create rescue repos on GitHub (fork or clone from git.bioconductor.org)
python scripts/rescue_repos.py

# 3. Clone repos locally under /Users/Levi/git/bioc-package-rescue/
python scripts/clone_repos.py

# 4. Apply/update centralized GHA workflow to all local checkouts
python scripts/update_reusable_workflows.py

# 5. Regenerate the README.md dashboard
python scripts/update_deprecated_packages.py
```

Steps 2–5 are idempotent — already-existing repos, clones, and workflows are skipped automatically.

---

## 1. Syncing the Package List (`sync_packages.py`)

Fetches the Bioconductor [Help Wanted page](https://bioconductor.org/developers/help_wanted/) and updates `packages.csv`:
- New packages from the **Deprecated** and **Voluntarily Listed** sections are added automatically.
- Existing rows (including `source=manual` entries) are **never removed or overwritten**.

```bash
python scripts/sync_packages.py
```

After running, commit the updated `packages.csv` to the repo.

---

## 2. Creating Rescue Repos (`rescue_repos.py`)

For each package in `packages.csv` with Bioconductor build status `ERROR` or `TIMEOUT`, creates a public repository in the `bioc-package-rescue` GitHub organization:

- If the package has an existing **GitHub repository** (scraped from its Bioconductor landing page), it is **forked** into the organization.
- Otherwise, it is **cloned from `git.bioconductor.org`** and pushed as a new repo.
- Already-existing repos are skipped.

Requires the `gh` CLI authenticated to an account with org-level repo creation permissions.

```bash
python scripts/rescue_repos.py
```

---

## 3. Cloning Repos Locally (`clone_repos.py`)

Clones each rescue repository from GitHub into the local workspace under `/Users/Levi/git/bioc-package-rescue/<PackageName>`. Only packages with build status ERROR are cloned. Already-existing local directories are skipped.

```bash
python scripts/clone_repos.py
```

---

## 4. Centralized GHA Workflow

### Central repository
- **URL**: `https://github.com/bioc-package-rescue/workflows`
- **File**: `.github/workflows/check-bioc.yml`

This repository holds the full workflow definition. To update the Bioconductor release version (e.g., `RELEASE_3_23` → `RELEASE_3_24`):
1. Edit `bioconductor/bioconductor_docker:RELEASE_3_23` in the matrix config.
2. Commit and push to `main`. All 49+ package repositories immediately use the new version.

### Package-level caller stub

Each package repository contains a minimal `.github/workflows/check-bioc.yml`:

```yaml
name: R-CMD-check-bioc

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
```

### Batch applying / refreshing stubs (`update_reusable_workflows.py`)

Writes the caller stub to every local checkout that has build status ERROR, commits, and pushes:

```bash
python scripts/update_reusable_workflows.py
```

---

## 5. Regenerating the Dashboard (`update_deprecated_packages.py`)

Reads `packages.csv`, fetches current Bioconductor build status and 2025 download statistics, and rewrites `README.md` with two tables (Deprecated Packages and Voluntarily Listed / Manual packages).

```bash
python scripts/update_deprecated_packages.py
```
