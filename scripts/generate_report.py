import os
import subprocess
import re
import csv

packages = [
    "BPRMeth", "BiRewire", "BubbleTree", "CSAR", "CelliD", "DeconRNASeq", 
    "GEOexplorer", "GNET2", "IMAS", "IONiseR", "MSPrep", "MetaNeighbor", 
    "MethReg", "MineICA", "Organism.dplyr", "RTCGA", "RcisTarget", "RgnTX", 
    "RiboProfiling", "SGCP", "SQLDataFrame", "SigFuge", "ballgown", 
    "barcodetrackR", "basecallQC", "biobroom", "biodbChebi", "ccrepe", 
    "debCAM", "netZooR", "netprioR", "partCNV", "rols", "scviR", "soGGi", 
    "tLOH", "tidytof"
]

def run_git_cmd(cwd, cmd):
    try:
        result = subprocess.run(
            cmd, cwd=cwd, shell=True, text=True, capture_output=True, check=True, errors='replace'
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        return ""

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dashboard_dir = os.path.dirname(script_dir)
    workspace = os.path.dirname(dashboard_dir)
    output_file = os.path.join(dashboard_dir, "package_fix_details.md")
    csv_file = os.path.join(dashboard_dir, "package_fix_stats.csv")
    
    with open(output_file, "w") as f, open(csv_file, "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Package", "Files Changed", "Insertions", "Deletions"])
        
        f.write("# Package Fix Details\n\n")
        f.write("This document contains the complete diffs, commit summaries, and line change statistics for the 37 Bioconductor packages rescued in this workspace.\n\n")
        
        for pkg in sorted(packages):
            repo_dir = os.path.join(workspace, pkg)
            if not os.path.isdir(repo_dir):
                print(f"Warning: {pkg} not found.")
                continue
            
            # Find our fix commits. We look for the co-authored trailer or we can just look at the last few commits
            # since we branched from devel. 
            # A robust way: get all commits, find the ones with our trailer or that we authored.
            log_cmd = "git log --pretty=format:'%H%x09%s%x09%b'"
            log_output = run_git_cmd(repo_dir, log_cmd)
            
            our_commits = []
            oldest_hash = None
            
            for line in log_output.split('\n'):
                if not line: continue
                parts = line.split('\t', 2)
                if len(parts) >= 2:
                    c_hash = parts[0]
                    c_subj = parts[1]
                    c_body = parts[2] if len(parts) > 2 else ""
                    
                    c_lower = (c_subj + "\n" + c_body).lower()
                    # Filter out GHA workflow commits
                    if "github action" in c_lower or "gha " in c_lower or "check-bioc.yml" in c_lower or "workflow" in c_lower:
                        continue
                        
                    # Check if it's our commit
                    if "antigravity" in c_lower or "gemini@google.com" in c_lower:
                        our_commits.append((c_hash, c_subj))
                        oldest_hash = c_hash
            
            if not our_commits:
                print(f"No substantive fix commits found for {pkg}. Trying fallback (last commit).")
                # Fallback: just diff HEAD^ to HEAD if we didn't use the trailer
                fallback_cmd = "git log -1 --pretty=format:'%H%x09%s'"
                fallback_out = run_git_cmd(repo_dir, fallback_cmd).split('\t')
                if len(fallback_out) >= 2:
                    oldest_hash = fallback_out[0]
                    our_commits.append((fallback_out[0], fallback_out[1]))
            
            if not oldest_hash:
                continue
                
            # To get the diff of all our changes, we diff the parent of the oldest fix commit vs HEAD.
            diff_base = f"{oldest_hash}^"
            
            # Substantive commit subjects
            # Reverse to chronological order
            our_commits.reverse()
            subjects = [subj for (h, subj) in our_commits]
            
            # Diff and stats
            diff_cmd = f"git diff {diff_base} HEAD -- . ':(exclude).github'"
            diff_content = run_git_cmd(repo_dir, diff_cmd)
            
            stat_cmd = f"git diff --shortstat {diff_base} HEAD -- . ':(exclude).github'"
            stat_content = run_git_cmd(repo_dir, stat_cmd).strip()
            
            if not stat_content:
                stat_content = "0 files changed, 0 insertions(+), 0 deletions(-)"
                
            files = insertions = deletions = 0
            match_files = re.search(r'(\d+) file', stat_content)
            if match_files: files = int(match_files.group(1))
            match_ins = re.search(r'(\d+) insertion', stat_content)
            if match_ins: insertions = int(match_ins.group(1))
            match_del = re.search(r'(\d+) deletion', stat_content)
            if match_del: deletions = int(match_del.group(1))
            
            csv_writer.writerow([pkg, files, insertions, deletions])
                
            # Write to markdown
            f.write(f"## {pkg}\n\n")
            f.write("**Substantive Commits:**\n")
            for subj in subjects:
                f.write(f"- {subj}\n")
            f.write("\n")
            
            f.write("**Line Changes:**\n")
            f.write(f"`STAT_LINES_CHANGED: {pkg} | {stat_content}`\n\n")
            
            f.write("**Complete Diff:**\n")
            f.write("```diff\n")
            f.write(diff_content + "\n")
            f.write("```\n\n")
            
            print(f"Processed {pkg}")

if __name__ == "__main__":
    main()
