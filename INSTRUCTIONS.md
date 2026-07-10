# Bioconductor Package Rescue Instructions

This document describes how to maintain the `bioc-rescue-dashboard`, manage rescue repositories, and update the centralized GitHub Actions check workflow.

---

## Package List: `packages.csv`

**`packages.csv`** (in the root of this repo) is the single source of truth for all tracked packages. Every script reads from it — nothing is hard-coded in the scripts themselves.

### Format

```csv
package,type,source
BPRMeth,Software,deprecated
VariantAnnotation,Software,voluntarily-listed
myCustomPkg,Software,manual
```

| Column    | Values                                            |
|-----------|---------------------------------------------------|
| `package` | Bioconductor package name                         |
| `type`    | `Software`, `ExperimentData`, or `AnnotationData` |
| `source`  | `deprecated`, `voluntarily-listed`, or `manual`   |

### Adding a Package Manually

Edit `packages.csv` directly, add a row with `source=manual`, and commit. The next dashboard run will pick it up automatically.

---

## Configuration

Two scripts hard-code the local workspace root path. If running on a different machine, update `WORKSPACE_ROOT` near the top of:

- `scripts/clone_repos.py`
- `scripts/update_reusable_workflows.py`

Default: `/Users/Levi/git/bioc-package-rescue`

---

## Typical Workflow

Run in order when setting up fresh or after Bioconductor publishes updates:

```bash
# 1. Sync package list from Bioconductor Help Wanted page → packages.csv
python scripts/sync_packages.py
git add packages.csv && git commit -m "Sync packages.csv"

# 2. Create rescue repos on GitHub (fork or clone from git.bioconductor.org)
python scripts/rescue_repos.py

# 3. Clone rescue repos locally
python scripts/clone_repos.py

# 4. Apply centralized GHA workflow to all local checkouts
python scripts/update_reusable_workflows.py

# 5. Regenerate the README.md dashboard
python scripts/update_deprecated_packages.py
git add README.md && git commit -m "Update dashboard"
```

Steps 2–5 are idempotent — already-existing repos, clones, and workflows are skipped automatically.

---

## Script Reference

### `sync_packages.py` — Update the package list

Fetches the Bioconductor [Help Wanted page](https://bioconductor.org/developers/help_wanted/) and merges it into `packages.csv`:
- New packages from the **Deprecated** and **Voluntarily Listed** sections are added.
- Existing rows (including `source=manual` entries) are **never removed or overwritten**.

### `rescue_repos.py` — Create GitHub repos

For each package in `packages.csv` with build status `ERROR` or `TIMEOUT`, creates a public repository in the `bioc-package-rescue` organization:
- If the package has an upstream **GitHub repository**, it is **forked**.
- Otherwise it is **cloned from `git.bioconductor.org`** and pushed as a new repo.
- Packages whose org repo already exists are skipped.

Requires the `gh` CLI authenticated to an account with organization-level repository creation permissions.

### `clone_repos.py` — Clone repos locally

Clones each rescue repo (ERROR/TIMEOUT packages only) from `git@github.com:bioc-package-rescue/<pkg>` into the local workspace. Already-existing local directories are skipped.

### `update_reusable_workflows.py` — Push GHA workflow stubs

Writes the minimal caller stub (see below) to `.github/workflows/check-bioc.yml` in every local checkout with ERROR/TIMEOUT status, commits with `Antigravity <gemini@google.com>` as author, and pushes. Already-up-to-date repos are skipped.

### `update_deprecated_packages.py` — Regenerate dashboard README

Reads `packages.csv`, fetches the current Bioconductor build status and most recent complete-year download statistics, and rewrites `README.md` with two tables: **Deprecated Packages** and **Voluntarily Listed / Manual packages**.

---

## Centralized GHA Workflow

### Central repository

- **URL**: `https://github.com/bioc-package-rescue/workflows`
- **File**: `.github/workflows/check-bioc.yml`

This holds the full workflow definition. To update the Bioconductor release version:
1. Edit the `image:` lines in the matrix config:
   ```yaml
   - { bioc: 'release', image: 'bioconductor/bioconductor_docker:RELEASE_3_23' }
   ```
2. Commit and push to `main`. All package repositories pick up the change on their next run — no per-repo edits needed.

### Package-level caller stub

Each rescued package repository contains a minimal `.github/workflows/check-bioc.yml`:

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

Run `update_reusable_workflows.py` to batch-apply or refresh this stub across all local checkouts.
