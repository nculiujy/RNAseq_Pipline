import os

def get_rawdata_dir(wildcards):
    for p in config["projects"]:
        if p["species"] == wildcards.species:
            return p["rawdata_dir"]
    return ""

def get_hisat2_index(wildcards):
    for p in config["projects"]:
        if p["species"] == wildcards.species:
            return p.get("hisat2_index", "")
    return ""

rule Align:
    input:
        script="workflow/scripts/3_1_Align.py"
    output:
        marker="result/{species}/3_Align_Filter/Align_finished.txt"
    params:
        inputdir=get_rawdata_dir,
        outputdir="result/{species}/3_Align_Filter",
        threads=config.get("align_threads", 4),
        max_parallel=config.get("align_parallel", 2),
        hisat2_index=get_hisat2_index
    log:
        "logs/{species}_align.log"
    threads: 8
    shell:
        """
        python {input.script} \
            --inputdir {params.inputdir} \
            --outputdir {params.outputdir} \
            --threads {params.threads} \
            --maxparallel {params.max_parallel} \
            --hisat2_index {params.hisat2_index} > {log} 2>&1
        touch {output.marker}
        """

def get_gtf_file(wildcards):
    for p in config["projects"]:
        if p["species"] == wildcards.species:
            return p.get("gtf_file", "")
    return ""

rule Quant:
    input:
        script="workflow/scripts/3_2_Quant.py",
        align_flag="result/{species}/3_Align_Filter/Align_finished.txt"
    output:
        marker="result/{species}/3_Align_Filter/Quant_finished.txt"
    params:
        inputdir="result/{species}/3_Align_Filter",
        outputdir="result/{species}/3_Align_Filter",
        threads=config.get("quant_threads", 4),
        max_parallel=config.get("quant_parallel", 2),
        gtf_base=config.get("gtf_base", "/home/public_software_annotation"),
        gtf_file=get_gtf_file
    log:
        "logs/{species}_quant.log"
    threads: 8
    shell:
        """
        python {input.script} \
            --inputdir {params.inputdir} \
            --outputdir {params.outputdir} \
            --threads {params.threads} \
            --maxparallel {params.max_parallel} \
            --gtf_base {params.gtf_base} \
            --gtf_file {params.gtf_file} > {log} 2>&1
        touch {output.marker}
        """

rule Filter:
    input:
        script="workflow/scripts/3_3_Filter.py",
        quant_flag="result/{species}/3_Align_Filter/Quant_finished.txt"
    output:
        marker="result/{species}/3_Align_Filter/Filter_finished.txt",
        csv="result/{species}/3_Align_Filter/alignment_quality.csv"
    params:
        inputdir="result/{species}/3_Align_Filter",
        outputdir="result/{species}/3_Align_Filter"
    log:
        "logs/{species}_filter.log"
    threads: 1
    shell:
        """
        python {input.script} \
            --inputdir {params.inputdir} \
            --outputdir {params.outputdir} > {log} 2>&1
        """

rule Merge:
    input:
        script="workflow/scripts/3_4_merge.py",
        filter_flag="result/{species}/3_Align_Filter/Filter_finished.txt",
        qc_csv="result/{species}/3_Align_Filter/alignment_quality.csv"
    output:
        marker="result/{species}/3_Align_Filter/Merge_finished.txt"
    params:
        inputdir="result/{species}/3_Align_Filter",
        outputdir="result/{species}/3_Align_Filter/Matrices",
        gtf_base=config.get("gtf_base", "/home/public_software_annotation")
    log:
        "logs/{species}_merge.log"
    threads: 1
    shell:
        """
        python {input.script} \
            --inputdir {params.inputdir} \
            --outputdir {params.outputdir} \
            --filter_csv {input.qc_csv} \
            --gtf_base {params.gtf_base} > {log} 2>&1
        touch {output.marker}
        """
