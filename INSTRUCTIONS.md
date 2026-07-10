# Bioconductor Package Rescue Instructions

This document provides instructions for agents and developers on how to maintain the `bioc-package-rescue` dashboard, rescue new packages, manage local checkouts, and update the centralized GitHub Actions check workflows.

---

## 1. Updating the Dashboard
The dashboard `README.md` displays two main tables of packages (Deprecated and Voluntarily Listed) with details on their type, Bioconductor build status, 2025 download statistics (distinct IPs), and links to their rescue repository and build status.

To regenerate and update the dashboard:
1. Ensure your current working directory is `bioc-rescue-dashboard/`.
2. Run the update script:
   ```bash
   python scripts/update_deprecated_packages.py
   ```
3. Commit and push the updated `README.md` to GitHub.

---

## 2. Rescuing New Packages (Creating Repos)
For packages on the Bioconductor "Help Wanted" list that have a build status of `ERROR` or `TIMEOUT`, we create a public rescue repository in the `bioc-package-rescue` GitHub organization.

The script `scripts/rescue_repos.py` automates this process:
1. **GitHub Auth**: Ensure the `gh` CLI is installed and authenticated to an account with permissions to create repositories in the `bioc-package-rescue` organization.
2. Run the rescue script:
   ```bash
   python scripts/rescue_repos.py
   ```
   *Behavior*:
   - If the package has an existing GitHub repository (scraped from the Bioconductor package landing page), the script will **fork** it into the `bioc-package-rescue` organization.
   - If no GitHub repository is found, the script will clone from `git.bioconductor.org` and push it to a **new public repository** in the organization.

---

## 3. Cloning Repositories Locally
To manage local checkouts of the rescued repositories under the workspace root (`/Users/Levi/git/bioc-package-rescue/`):
1. Run the clone script:
   ```bash
   python scripts/clone_repos.py
   ```
   *Behavior*:
   - Scrapes the build databases to identify packages with errors.
   - Clones their rescue repositories into `/Users/Levi/git/bioc-package-rescue/{package_name}` if they do not already exist locally.

---

## 4. Centralized Bioconductor Check Workflow (`check-bioc.yml`)
To maintain consistency and make future updates easy, the GitHub Actions check workflow is centralized.

### Central Workflow Repository
- **Repo URL**: `https://github.com/bioc-package-rescue/workflows`
- **File Path**: `.github/workflows/check-bioc.yml`

This workflow defines the actual execution environment (runs inside the `bioconductor/bioconductor_docker` container) and checks R CMD check and BiocCheck against the release and devel versions in parallel.

### Updating the central check version (e.g. updating RELEASE_3_23)
1. Clone the `workflows` repository:
   ```bash
   git clone git@github.com:bioc-package-rescue/workflows.git
   ```
2. Edit `.github/workflows/check-bioc.yml`. Locate the strategy matrix config:
   ```yaml
   matrix:
     config:
       - { bioc: 'release', image: 'bioconductor/bioconductor_docker:RELEASE_3_23' }
       - { bioc: 'devel',   image: 'bioconductor/bioconductor_docker:devel'         }
   ```
3. Update the `image` string for release or devel as needed.
4. Commit and push the changes to `main` branch.
   *Note: All package checkouts referencing `@main` will instantly use the updated image configuration on their next run.*

### Package-level Workflow Setup
Each individual package repository has a minimal stub at `.github/workflows/check-bioc.yml` that delegates execution to the central workflow:
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

### Batch Updating Packages
If you need to batch-apply or reset the package-level workflow files across all rescued repositories:
1. Ensure all local checkouts exist under the workspace root.
2. Run the update script:
   ```bash
   python scripts/update_reusable_workflows.py
   ```
   *Behavior*:
   - Iterates through all local package directories.
   - Overwrites `.github/workflows/check-bioc.yml` in each repository with the caller stub.
   - Commits and pushes the change to the active remote branch using `Antigravity <gemini@google.com>` as author.
