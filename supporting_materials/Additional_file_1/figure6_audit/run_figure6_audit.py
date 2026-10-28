"""Run the Figure 6E audit from an extracted Additional file 4 in isolation."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

parser=argparse.ArgumentParser()
parser.add_argument("--attachment",type=Path,required=True,help="Extracted Additional_file_4 directory")
parser.add_argument("--workspace",type=Path,required=True,help="New empty writable audit workspace")
parser.add_argument("--rscript",default="Rscript")
args=parser.parse_args()
attachment=args.attachment.resolve()
workspace=args.workspace.resolve()
if workspace.exists() and any(workspace.iterdir()):
    raise SystemExit("Use a new empty workspace; existing evidence will not be overwritten.")
workspace.mkdir(parents=True,exist_ok=True)
inputs={
 "audit_input/spots_galnt_annot.csv.gz":"02_intermediate/fig6_visium/spots_galnt_annot.csv.gz",
 "tables/H1_panin_enzymes_robustness.csv":"04_results/2026-10-12_Figure6_GALNT_run01/H1_panin_enzymes_robustness.csv",
}
for source,target in inputs.items():
    p=workspace/target;p.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(attachment/source,p)
env=os.environ.copy();env["PROJECT_ROOT"]=str(workspace)
subprocess.run([args.rscript,str(attachment/"provenance/verify_figure6.R")],env=env,check=True)
print("Audit complete. Outputs are in the isolated workspace's dated report folder.")
