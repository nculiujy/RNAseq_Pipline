rule QC:
    input:
        script="workflow/scripts/2_QC_fqfile.pl",
        download_flag=os.path.join(config["result_dir"], "1_download_geo", "finished.txt")
    output:
        marker=os.path.join(config["result_dir"], "2_QC", "QC_finished.txt")
    params:
        inputdir=os.path.join(config["result_dir"], "1_download_geo"),
        outputdir=os.path.join(config["result_dir"], "2_QC"),
        threads=config.get("qc_threads", 4),
        max_parallel=config.get("qc_parallel", 4)
    log:
        os.path.join(config["log_dir"], "2_QC", "qc.log")
    threads: 4
    shell:
        """
        perl {input.script} \
            --inputdir {params.inputdir} \
            --outputdir {params.outputdir} \
            --threads {params.threads} \
            --maxparallel {params.max_parallel} > {log} 2>&1
        """
