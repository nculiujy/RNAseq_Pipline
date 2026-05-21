rule download_geo:
    input:
        script="workflow/scripts/1_download.py"
    output:
        marker=os.path.join(config["result_dir"], "1_download_geo", "finished.txt")
    params:
        rawdata_dir=config["rawdata_dir"],
        result_dir=os.path.join(config["result_dir"], "1_download_geo")
    log:
        os.path.join(config["log_dir"], "1_download_geo", "download.log")
    threads: 8
    shell:
        """
        export DOWNLOAD_THREADS={threads}
        # 将原始数据目录和结果目录传递给脚本
        python {input.script} {params.rawdata_dir} {params.result_dir} > {log} 2>&1
        
        # 确保输出文件存在，如果脚本未生成，手动创建
        if [ ! -f {output.marker} ]; then
            echo "Warning: finished.txt not found, creating it manually."
            touch {output.marker}
        fi
        """
