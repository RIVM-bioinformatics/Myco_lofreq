#############################################################################
##### 			Low-frequency variants analysis  		#####
#############################################################################

rule vcf_annotator:
    input:
        str(OUT / "lofreq/{sample}.vcf")
    output:
        var=temp(str(OUT / "annotated_vcf/{sample}.vcf")),
        annot_var=str(OUT / "annotated_vcf/{sample}_annotated.vcf")
    conda:
        "../../envs/vcf-annotator.yaml"
    benchmark:
        str(OUT / "log/benchmark/annotated_vcf_{sample}.txt")
    threads: config["threads"]["annotated_vcf"]
    resources:
        mem_mb=config["mem_mb"]["annotated_vcf"]
    log:
        str(OUT / "log/annotated_vcf/annotated_vcf_{sample}.log")
    params:
        genebank_file=config["vcf_annotator"]["annotated_genome"]
    shell:
        """
sed 's/^NC_000962\.3\t/NC_000962\t/g' {input} > {output.var} 2> {log}
vcf-annotator --output {output.annot_var} {output.var} {params} 2> {log}
        """