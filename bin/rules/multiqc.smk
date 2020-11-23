#############################################################################
##### 			Data quality control and cleaning  		#####
#############################################################################

rule multiqc:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "R1 R2".split()),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "pR1 pR2 uR1 uR2".split()),
        expand(str(OUT / "log/trimmomatic/trimmomatic_{sample}.log"), sample = SAMPLES)
    output:
        str(OUT / "MultiQC/multiqc.html"),
    conda:
        "../../envs/multiqc.yaml"
    benchmark:
        str(OUT / "log/benchmark/MultiQC.txt")
    threads: config["threads"]["fastqc"]
    resources:
        mem_mb=config["mem_mb"]["fastqc"]
    params:
        output_dir=str(OUT / "MultiQC"),
        config_file="files/multiqc_config.yaml"
    log:
        str(OUT / "log/multiqc/MultiQC_report.log")
    shell:
        """
multiqc --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
    """
