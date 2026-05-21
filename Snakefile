import os

configfile: "config/config.yaml"

TARGET_FILES = []

PROJECTS = config.get("projects", [])

include: "workflow/rules/3_Align_Filter.smk"
include: "workflow/rules/4_DEseq.smk"

for proj in PROJECTS:
    species = proj["species"]
    modules = proj.get("modules", {})
    
    if modules.get("3_Align_Filter", False):
        TARGET_FILES.append(os.path.join("result", species, "3_Align_Filter", "Merge_finished.txt"))
    if modules.get("4_DEseq", False):
        TARGET_FILES.append(os.path.join("result", species, "4_DEseq", "DEseq_finished.txt"))

rule all:
    input:
        TARGET_FILES
